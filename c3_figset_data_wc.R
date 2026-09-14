# ==============================================================================
# Chapter 3 — consolidate data for the figure set (genuine-53, v3.2.3):
#  (A) LAI-space: which LAI product predicts observed microclimate (cor).
#  (B/D) per-scenario MuSICA validation: R2/RMSE/bias/SD-recovery, ΔTmax & slope.
#  (C) magnitude axis: warm bias vs mean LAI fed.
# Writes CSVs into NC_Full/manuscripts/ch3/tables/ (figX_*.csv).
#   Rscript c3_figset_data.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(terra);library(randomForest);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); dopt<-CFG_C3$d_opt_m
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
# --- obs microclimate + macro ---
FMACRO<-"in_files/musica_in_Blois_pblh.nc"  # macro MUST match the forcing the windcorr runs used (pblh), NOT CFG_C3$forcing_file (FR-Blo) — see pipeline audit
ncf<-nc_open(FMACRO);MA<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
sdo_d<-sd(OBS$do); sdo_s<-sd(OBS$so,na.rm=TRUE)

# ============ (B/D) per-scenario MuSICA validation ============
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); d1<-"out_files/Chapter3/nc_genuine53_windcorr"
SC<-list(
 STATIC_ALS      =list(d=d3,lab="LiDAR full (static)", fam="LiDAR"),
 DYN_ALS         =list(d=d3,lab="LiDAR mag x S2 timing",fam="LiDAR"),
 CONST_ALS       =list(d=d3,lab="LiDAR const (year)",  fam="LiDAR"),
 STATIC_S2_ATBD  =list(d=d3,lab="S2 ATBD (static)",    fam="Sentinel-2"),
 DYN_S2_ATBD     =list(d=d3,lab="S2 ATBD (dynamic)",   fam="Sentinel-2"),
 NAIVE_S2_FORMSH =list(d=d3,lab="S2 + FORMS-H height", fam="Sentinel-2"),
 STATIC_S2_DOPT  =list(d=d3,lab="S2 at d_opt depth",   fam="Sentinel-2"),
 FUSION_H        =list(d=d1,lab="Layered fusion",      fam="Fusion"),
 STATIC_ALS_DOPT =list(d=d3,lab="LiDAR truncated d_opt",fam="d_opt truncated"))
getpp<-function(scn,dir){fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
 m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr,Tm)],by="hr")
  sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
 merge(m,OBS,by="id_plot")}
B<-rbindlist(lapply(names(SC),function(s){d<-getpp(s,SC[[s]]$d)
 data.table(scenario=SC[[s]]$lab, family=SC[[s]]$fam,
  dT_R2=cor(d$ds,d$do)^2, dT_RMSE=sqrt(mean((d$ds-d$do)^2)), dT_bias=mean(d$ds-d$do), dT_SDrec=sd(d$ds)/sdo_d,
  sl_R2=cor(d$ss,d$so,use="complete.obs")^2, sl_RMSE=sqrt(mean((d$ss-d$so)^2,na.rm=TRUE)), sl_bias=mean(d$ss-d$so,na.rm=TRUE), sl_SDrec=sd(d$ss,na.rm=TRUE)/sdo_s)}))
fwrite(B,file.path(TAB,"figB_scenario_metrics.csv")); cat("wrote figB_scenario_metrics.csv\n")

# ============ (C) magnitude axis: mean LAI fed vs bias ============
magsc<-list(STATIC_ALS="LiDAR full",STATIC_S2_ATBD="S2 raw",STATIC_S2_DOPT="S2 at d_opt",STATIC_ALS_DOPT="LiDAR trunc. d_opt")
magcol<-c(STATIC_ALS="LAI_ALS",STATIC_S2_ATBD="LAI_S2_ATBD",STATIC_S2_DOPT="LAI_S2_DOPT",STATIC_ALS_DOPT="LAI_ALS_DOPT")
C<-rbindlist(lapply(names(magsc),function(s){d<-getpp(s,d3)
 data.table(scenario=magsc[[s]], meanLAI=mean(df[[magcol[[s]]]]), dT_bias=mean(d$ds-d$do), dT_RMSE=sqrt(mean((d$ds-d$do)^2)))}))
fwrite(C,file.path(TAB,"figC_magnitude_axis.csv")); cat("wrote figC_magnitude_axis.csv\n")

# ============ (A) LAI-space: product predicts observed microclimate ============
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"; pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
dynsum<-function(pat){files<-list.files(NM,pattern=pat,full.names=TRUE);dts<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
 vals<-as.data.frame(terra::extract(rast(files),pts))[,-1,drop=FALSE];long<-rbindlist(lapply(seq_along(dts),function(i)data.table(plot_id=df$pid,doy=as.integer(format(dts[i],"%j")),lai=pmax(vals[[i]],0))))
 ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3); sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
df$S2opt<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$"); df$S2atbd<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$"); df$hfrac<-pmax(df$FORMS_H,0.1)/dopt
Dl<-as.data.frame(merge(OBS,df[,.(id_plot,LAI_ALS,S2opt,S2atbd,FORMS_H,hfrac)],by="id_plot"))
looRF<-function(fs){set.seed(42);f<-as.formula(paste("LAI_ALS~",paste(fs,collapse="+")));sapply(1:nrow(Dl),function(i)predict(randomForest(f,Dl[-i,],ntree=400),Dl[i,]))}
Dl$LAIcorr<-looRF(c("S2opt","S2atbd","FORMS_H","hfrac"))
prods<-list("LiDAR full"="LAI_ALS","Corrected S2 (operational)"="LAIcorr","FORMS-H height"="FORMS_H","S2 opt (dynamic)"="S2opt","S2 ATBD (dynamic)"="S2atbd")
A<-rbindlist(lapply(names(prods),function(nm){v<-Dl[[prods[[nm]]]];data.table(product=nm, cor_dTmax=cor(v,Dl$do), cor_slope=cor(v,Dl$so))}))
fwrite(A,file.path(TAB,"figA_laispace_cor.csv")); cat("wrote figA_laispace_cor.csv\n")
cat("\n=== A ===\n");print(A);cat("\n=== C ===\n");print(C)
