# ==============================================================================
# Where should the dense/open split fall? (windcorr genuine-53)
#  (1) recompute stratified R2 (LiDAR STATIC_ALS vs S2 STATIC_S2_ATBD) at the Ch2
#      crossover LAI=3.86 (consistent with the fusion gate);
#  (2) robustness: sweep the LAI split threshold -> LiDAR-dense-R2 & S2-open-R2;
#  (3) locate the boundary PHYSICALLY from the buffering regime: 2-piece breakpoint
#      of observed slope~LAI and observed ΔTmax~LAI (base R grid search), plus the
#      slope distribution. NB: we split on a STRUCTURAL var (LAI), never on the
#      observed outcome, to avoid circularity; slope is used only to justify where.
#   Rscript c3_split_explore_wc.R
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
L<-getpp("STATIC_ALS");setnames(L,c("ds","ss"),c("dL","sL")); S<-getpp("STATIC_S2_ATBD");setnames(S,c("ds","ss"),c("dS","sS"))
M<-merge(merge(merge(OBS,L,by="id_plot"),S,by="id_plot"),df[,.(id_plot,LAI_ALS,FORMS_H)],by="id_plot")
r2<-function(a,b) if(sum(is.finite(a)&is.finite(b))>=4) cor(a,b,use="complete.obs")^2 else NA

## (1) stratified at the chosen regime split 3.86
strat<-function(thr){de<-M$LAI_ALS>=thr;rbind(
  data.table(thr=thr,stratum="dense",n=sum(de),LiDAR_dT=r2(M$dL[de],M$do[de]),S2_dT=r2(M$dS[de],M$do[de]),LiDAR_sl=r2(M$sL[de],M$so[de]),S2_sl=r2(M$sS[de],M$so[de])),
  data.table(thr=thr,stratum="open", n=sum(!de),LiDAR_dT=r2(M$dL[!de],M$do[!de]),S2_dT=r2(M$dS[!de],M$do[!de]),LiDAR_sl=r2(M$sL[!de],M$so[!de]),S2_sl=r2(M$sS[!de],M$so[!de])))}
cat("=== (1) stratified at regime split LAI=3.86 ===\n");print(strat(CROSS),digits=3)
fwrite(strat(CROSS),file.path(TAB,"figJ_stratified_cross386.csv"))

## (2) threshold-sweep robustness
sw<-rbindlist(lapply(seq(2.0,5.5,0.25),function(t){de<-M$LAI_ALS>=t
  data.table(thr=t,n_dense=sum(de),n_open=sum(!de),LiDAR_dense_dT=r2(M$dL[de],M$do[de]),S2_open_dT=r2(M$dS[!de],M$do[!de]),
             S2_dense_dT=r2(M$dS[de],M$do[de]),LiDAR_open_dT=r2(M$dL[!de],M$do[!de]))}))
cat("\n=== (2) split-threshold sweep (reversal robustness) ===\n");print(sw,digits=2)
fwrite(sw,file.path(TAB,"figJ_threshold_sweep.csv"))

## (3) physical boundary from observed buffering ~ LAI (base-R 2-piece breakpoint)
bp<-function(y){x<-M$LAI_ALS;ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  cand<-seq(quantile(x,.2),quantile(x,.8),length.out=40)
  rss<-sapply(cand,function(c){xl<-pmin(x,c);xr<-pmax(x-c,0);sum(resid(lm(y~xl+xr))^2)})
  cand[which.min(rss)]}
cat(sprintf("\n=== (3) buffering-regime breakpoint (observed) ===\n"))
cat(sprintf("breakpoint of observed slope ~ LAI : LAI = %.2f\n",bp(M$so)))
cat(sprintf("breakpoint of observed ΔTmax ~ LAI : LAI = %.2f\n",bp(M$do)))
cat(sprintf("observed slope: median=%.2f  IQR[%.2f,%.2f]  range[%.2f,%.2f]\n",median(M$so,na.rm=T),quantile(M$so,.25,na.rm=T),quantile(M$so,.75,na.rm=T),min(M$so,na.rm=T),max(M$so,na.rm=T)))
cat(sprintf("median LAI (old split)=%.2f | split used=%.2f\n",median(M$LAI_ALS),CROSS))
cat("\nwrote figJ_stratified_cross386.csv + figJ_threshold_sweep.csv\n")
