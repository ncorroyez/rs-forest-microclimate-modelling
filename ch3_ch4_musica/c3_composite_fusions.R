# ==============================================================================
# Chapter 3 — explore composite fusion scenarios built by per-plot recombination
# of EXISTING MuSICA outputs (no new simulation). At the crossover LAI = 3.86:
#  * FUSION_H      (baseline): dense = LiDAR full, open = S2 raw (full LAD)
#  * FUSION_LAD    (S2-informed): dense = LiDAR full, open = S2 raw on d_opt LAD
#                  (LADOPT = leaf area concentrated in the depth S2 optically senses)
#  * FUSION_LAD_C  (fully S2-consistent open): dense = LiDAR full,
#                  open = S2 magnitude corrected to d_opt AND LAD truncated to d_opt
#  * ORACLE_2S     upper bound: per plot pick whichever of {LiDAR, S2 raw} is
#                  closer to observed dTmax (uses obs -> ceiling, not deployable)
# Also reports the per-regime de-biased RMSE floor for FUSION_H (offset removes
# the systematic warm bias; R2 unchanged).  Same windcorr nc / conventions.
#   -> chapter3_S2_LAI/tables/figK_composite_fusions.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"; CROSS<-3.86
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn){fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
R2<-function(a,b) suppressWarnings(cor(a,b,use="complete.obs")^2)
df[,st:=ifelse(LAI_ALS<CROSS,"open","dense")]
key<-df[,.(id_plot,st)]
# --- source scenarios (per plot ds=dTmax, ss=slope) --------------------------
SRC<-c("STATIC_ALS","STATIC_S2_ATBD","STATIC_S2_ATBD_LADOPT","STATIC_S2_DOPT_LADOPT")
S<-setNames(lapply(SRC,function(s)merge(getpp(s),key,by="id_plot")),SRC)
# --- composite builder: pick per-plot rows by stratum from two sources -------
compose<-function(dense_src,open_src){
  d<-S[[dense_src]][st=="dense"]; o<-S[[open_src]][st=="open"]
  merge(OBS,rbind(d,o)[,.(id_plot,ds,ss,st)],by="id_plot")}
# --- oracle: per plot pick source closer to observed dTmax -------------------
oracle<-function(a,b){
  A<-merge(OBS,S[[a]][,.(id_plot,dsA=ds,ssA=ss)],by="id_plot")
  B<-merge(A,S[[b]][,.(id_plot,dsB=ds,ssB=ss)],by="id_plot")
  B[,pick:=abs(dsA-do)<=abs(dsB-do)]
  B[,.(id_plot,do,so,ds=ifelse(pick,dsA,dsB),ss=ifelse(pick,ssA,ssB),st=key$st[match(id_plot,key$id_plot)])]}
metr<-function(M,lab){
  D<-M[st=="dense"];O<-M[st=="open"]
  data.table(scenario=lab,
    dT_dense=R2(D$ds,D$do),dT_open=R2(O$ds,O$do),dT_pool=R2(M$ds,M$do),
    dT_bias=mean(M$ds-M$do),dT_RMSE=sqrt(mean((M$ds-M$do)^2)),dT_SDrec=sd(M$ds)/sd(M$do),
    sl_pool=R2(M$ss,M$so))}
# per-regime de-biased RMSE floor (subtract stratum-mean bias; R2 unaffected)
rmse_debias<-function(M){Mc<-copy(M);Mc[,dsc:=ds-ave(ds-do,st),by=NULL];sqrt(mean((Mc$dsc-Mc$do)^2))}
FH  <-compose("STATIC_ALS","STATIC_S2_ATBD")
FL  <-compose("STATIC_ALS","STATIC_S2_ATBD_LADOPT")
FLC <-compose("STATIC_ALS","STATIC_S2_DOPT_LADOPT")
OR  <-oracle("STATIC_ALS","STATIC_S2_ATBD")
out<-rbind(metr(FH,"FUSION_H (dense LiDAR / open S2 full-LAD)"),
           metr(FL,"FUSION_LAD (open S2 on d_opt LAD, S2-informed)"),
           metr(FLC,"FUSION_LAD_C (open S2+LAD both d_opt)"),
           metr(OR,"ORACLE_2S (per-plot best of LiDAR/S2, ceiling)"))
out[,dT_RMSE_debias:=c(rmse_debias(FH),rmse_debias(FL),rmse_debias(FLC),rmse_debias(OR))]
fwrite(out,file.path(TAB,"figK_composite_fusions.csv"))
cat(sprintf("=== composite fusions, split LAI %.2f (25 dense / 28 open) ===\n",CROSS))
print(out[,.(scenario,dTd=round(dT_dense,2),dTo=round(dT_open,2),pool=round(dT_pool,2),bias=round(dT_bias,2),RMSE=round(dT_RMSE,2),RMSEdb=round(dT_RMSE_debias,2),SDrec=round(dT_SDrec,2),sl=round(sl_pool,2))],row.names=FALSE)
cat("\nwrote figK_composite_fusions.csv\n")
