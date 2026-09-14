# ==============================================================================
# Chapter 3 — data for additional figures E/F/G/H (genuine-53, v3.2.3).
#  E premise: observed & MuSICA buffering vs LiDAR LAI.
#  F why-S2-fails: LiDAR vs S2 LAI, coloured by openness.
#  G two-tier ablation: LOO recovery + microclimate skill by feature set.
#  H validation scatter: observed vs simulated ΔTmax (LiDAR full vs dynamic S2).
#   Rscript c3_figset_data2.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(terra);library(randomForest);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); dopt<-CFG_C3$d_opt_m
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
ncf<-nc_open(CFG_C3$forcing_file);MA<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn,dir){fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr,Tm)],by="hr")
  sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
d3<-file.path(CFG_C3$out_dir,"nc_genuine53")
S_ALS<-getpp("STATIC_ALS",d3); S_S2<-getpp("DYN_S2_ATBD",d3)
# S2 dynamic summer means
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"; pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
dynsum<-function(pat){files<-list.files(NM,pattern=pat,full.names=TRUE);dts<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
 vals<-as.data.frame(terra::extract(rast(files),pts))[,-1,drop=FALSE];long<-rbindlist(lapply(seq_along(dts),function(i)data.table(plot_id=df$pid,doy=as.integer(format(dts[i],"%j")),lai=pmax(vals[[i]],0))))
 ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3); sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
df$S2opt<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$"); df$S2atbd<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$"); df$hfrac<-pmax(df$FORMS_H,0.1)/dopt

## E premise
E<-merge(merge(OBS,df[,.(id_plot,LAI_ALS)],by="id_plot"),S_ALS[,.(id_plot,sim_do=ds,sim_so=ss)],by="id_plot")
fwrite(E,file.path(TAB,"figE_premise.csv"))
## F why S2 fails
F<-df[,.(id_plot,LAI_ALS,S2atbd,S2opt,fCover,Hmax)]; fwrite(F,file.path(TAB,"figF_s2_vs_lidar.csv"))
## H validation scatter
H<-merge(merge(OBS[,.(id_plot,do)],S_ALS[,.(id_plot,LiDAR=ds)],by="id_plot"),S_S2[,.(id_plot,S2=ds)],by="id_plot"); fwrite(H,file.path(TAB,"figH_scatter.csv"))
## G ablation (LOO RF)
Dl<-as.data.frame(merge(OBS,df[,.(id_plot,LAI_ALS,S2opt,S2atbd,FORMS_H,hfrac,Hmax,fCover,VCI,LCV)],by="id_plot"))
looRF<-function(fs){set.seed(42);f<-as.formula(paste("LAI_ALS~",paste(fs,collapse="+")));sapply(1:nrow(Dl),function(i)predict(randomForest(f,Dl[-i,],ntree=400),Dl[i,]))}
sets<-list("S2 only"=c("S2opt","S2atbd"),"S2 + FORMS-H"=c("S2opt","S2atbd","FORMS_H","hfrac"),"S2 + LiDAR structure"=c("S2opt","S2atbd","Hmax","fCover","VCI","LCV"),"LiDAR structure"=c("Hmax","fCover","VCI","LCV"))
G<-rbindlist(lapply(names(sets),function(nm){p<-looRF(sets[[nm]]);data.table(feature_set=nm,operational=grepl("FORMS|only",nm)&!grepl("LiDAR",nm),LOO_R2_LAI=cor(p,Dl$LAI_ALS)^2,cor_dTmax=cor(p,Dl$do),cor_slope=cor(p,Dl$so))}))
G<-rbind(G,data.table(feature_set="LiDAR LAI (direct)",operational=FALSE,LOO_R2_LAI=1,cor_dTmax=cor(Dl$LAI_ALS,Dl$do),cor_slope=cor(Dl$LAI_ALS,Dl$so)))
fwrite(G,file.path(TAB,"figG_ablation.csv"))
cat("wrote figE/figF/figG/figH csv\n"); print(G)
EOF_MARK<-1
