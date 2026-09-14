# ==============================================================================
# Chapter 3 — stratify EVERY available forcing scenario at the sensor crossover
# LAI = 3.86 (25 dense / 28 open) and rank them by the all-rounder criterion
# (maximin = min of the two within-stratum ΔTmax R2). Answers: which scenario
# wins where, which is the best overall (good in both regimes), and gives the
# full landscape to spot gaps. Same windcorr nc / conventions as
# c3_rf_stratified.R.  -> chapter3_S2_LAI/tables/figI_all_stratified.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config_CHS41.R")
TAB<-file.path(CFG_C3$out_dir,"tables"); ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc"); FORC<-CFG_C3$forcing_file; CROSS<-3.86
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
strata<-function(pids) df[match(pids,id_plot), ifelse(LAI_ALS<CROSS,"open","dense")]
R2<-function(a,b) suppressWarnings(cor(a,b,use="complete.obs")^2)
SCN<-list.dirs(WC,recursive=FALSE,full.names=FALSE)
SCN<-SCN[sapply(SCN,function(s)length(list.files(file.path(WC,s),pattern="\\.nc$"))==53)]
out<-rbindlist(lapply(SCN,function(scn){
  M<-merge(OBS,getpp(scn),by="id_plot"); M[,st:=strata(id_plot)]
  D<-M[st=="dense"]; O<-M[st=="open"]
  data.table(scenario=scn,
    dT_dense=R2(D$ds,D$do), dT_open=R2(O$ds,O$do), dT_pool=R2(M$ds,M$do),
    dT_bias=mean(M$ds-M$do), dT_SDrec=sd(M$ds)/sd(M$do),
    sl_pool=R2(M$ss,M$so))}))
out[,maximin:=pmin(dT_dense,dT_open)]
setorder(out,-maximin)
fwrite(out,file.path(TAB,"figI_all_stratified.csv"))
cat(sprintf("=== %d scenarios, split LAI %.2f, ranked by maximin(dense,open) ΔTmax R2 ===\n",nrow(out),CROSS))
print(out[,.(scenario,dT_dense=round(dT_dense,2),dT_open=round(dT_open,2),maximin=round(maximin,2),dT_pool=round(dT_pool,2),dT_bias=round(dT_bias,2),SDrec=round(dT_SDrec,2),sl=round(sl_pool,2))],row.names=FALSE)
cat("\nwrote figI_all_stratified.csv\n")
