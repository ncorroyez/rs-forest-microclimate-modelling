# ==============================================================================
# Chapter 3 §3.4 — stratified diagnostic (windcorr genuine-53). Split the 53 plots
# at the median LiDAR LAI (open vs dense) and compute, per stratum, the between-plot
# R2 / bias / RMSE for LiDAR full vs S2 ATBD, for both buffer metrics (ΔTmax, slope).
# This is the STRONG §3.4 argument: the pooled S2 R2 advantage is an open-canopy
# effect; in dense (saturated) canopy the LiDAR wins the ranking. Model-free split.
#   Rscript c3_stratified_wc.R  -> chapter3_S2_LAI/tables/figJ_stratified.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
FORC<-"in_files/musica_in_Blois_pblh.nc"
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
# macro
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
# obs
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn){fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
# scenarios to contrast
SC<-c(LiDAR="STATIC_ALS", "Sentinel-2"="STATIC_S2_ATBD")
med<-median(df$LAI_ALS)
strata<-function(pids){ df[match(pids,id_plot), ifelse(LAI_ALS<med,"Open (LAI<med)","Dense (LAI>=med)")] }
out<-rbindlist(lapply(names(SC),function(sensor){P<-getpp(SC[[sensor]]);M<-merge(OBS,P,by="id_plot")
  M[,stratum:=strata(id_plot)]
  M[,.(sensor=sensor, n=.N,
       dT_R2=cor(ds,do)^2, dT_bias=mean(ds-do), dT_RMSE=sqrt(mean((ds-do)^2)),
       sl_R2=cor(ss,so,use="complete.obs")^2, sl_bias=mean(ss-so,na.rm=TRUE)),by=stratum]}))
# also pooled
poolr<-rbindlist(lapply(names(SC),function(sensor){P<-getpp(SC[[sensor]]);M<-merge(OBS,P,by="id_plot")
  M[,.(sensor=sensor,stratum="Pooled (all 53)",n=.N,dT_R2=cor(ds,do)^2,dT_bias=mean(ds-do),dT_RMSE=sqrt(mean((ds-do)^2)),
       sl_R2=cor(ss,so,use="complete.obs")^2,sl_bias=mean(ss-so,na.rm=TRUE))]}))
out<-rbind(out,poolr)
fwrite(out,file.path(TAB,"figJ_stratified.csv"))
cat(sprintf("median LiDAR LAI split = %.2f\n",med)); print(out[order(sensor,stratum)],digits=3)
cat("\nwrote figJ_stratified.csv\n")
