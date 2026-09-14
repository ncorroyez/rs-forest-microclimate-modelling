# ==============================================================================
# Chapter 3 RQ4 — stratify the PRODUCT-LEVEL fusion (RF-corrected S2 magnitude,
# LAI_RF = random forest recovering LiDAR LAI from S2 + height features, forced
# on the full LiDAR PAD profile) at the sensor crossover LAI = 3.86 (25 dense /
# 28 open), exactly as Fig J. Tests whether calibrating the S2 magnitude toward
# the LiDAR LAI rescues the dense regime, or whether it inherits S2's open-canopy
# training bias and fails dense like raw S2 / NAIVE (a Simpson effect).
# Anchors: STATIC_ALS (reference), STATIC_S2_ATBD (raw S2), NAIVE_S2_FORMSH,
# FUSION_H (the switching-rule fusion already in the paper).
# Same windcorr nc, same conventions as c3_naive_stratified.R.
#   Rscript c3_rf_stratified.R  -> chapter3_S2_LAI/tables/figH_rf_stratified.csv
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
strata<-function(pids) df[match(pids,id_plot), ifelse(LAI_ALS<CROSS,"Open (LAI<3.86)","Dense (LAI>=3.86)")]
metr<-function(M) M[,.(n=.N, dT_R2=cor(ds,do)^2, dT_bias=mean(ds-do), dT_RMSE=sqrt(mean((ds-do)^2)),
                       dT_SDrec=sd(ds)/sd(do), sl_R2=cor(ss,so,use="complete.obs")^2, sl_bias=mean(ss-so,na.rm=TRUE))]
SCN<-c("STATIC_RF","STATIC_RF_DOPT","DYN_RF","STATIC_ALS","STATIC_S2_ATBD","NAIVE_S2_FORMSH","FUSION_H")
out<-rbindlist(lapply(SCN,function(scn){
  P<-getpp(scn); M<-merge(OBS,P,by="id_plot"); M[,stratum:=strata(id_plot)]
  rbind(cbind(scenario=scn,stratum="Dense (LAI>=3.86)",metr(M[stratum=="Dense (LAI>=3.86)"])),
        cbind(scenario=scn,stratum="Open (LAI<3.86)", metr(M[stratum=="Open (LAI<3.86)"])),
        cbind(scenario=scn,stratum="Pooled (all 53)", metr(M)))}))
fwrite(out,file.path(TAB,"figH_rf_stratified.csv"))
cat(sprintf("split at LAI %.2f\n",CROSS)); print(out,digits=3); cat("\nwrote figH_rf_stratified.csv\n")
