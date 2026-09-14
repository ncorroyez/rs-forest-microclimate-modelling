# ==============================================================================
# Regime-switching hybrid: LiDAR sim in dense canopy, S2 sim in open canopy.
# Tests whether picking the best sensor per stratum beats either alone (pooled).
# Two switch rules:
#   oracle  = switch on TRUE LiDAR-LAI median (upper bound; needs LAI to decide)
#   FORMS-H = switch on FORMS-H height median (operational, non-LiDAR switch)
# Also a LOO variant of the FORMS-H threshold to check overfitting.
#   Rscript c3_hybrid_switch_wc.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"
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
L<-getpp("STATIC_ALS"); setnames(L,c("ds","ss"),c("dL","sL"))
S<-getpp("STATIC_S2_ATBD"); setnames(S,c("ds","ss"),c("dS","sS"))
M<-merge(merge(merge(OBS,L,by="id_plot"),S,by="id_plot"),df[,.(id_plot,LAI_ALS,FORMS_H)],by="id_plot")
sdo_d<-sd(M$do); sdo_s<-sd(M$so,na.rm=TRUE)
mets<-function(dd,ss,tag){data.table(product=tag,
  dT_R2=cor(dd,M$do)^2,dT_bias=mean(dd-M$do),dT_RMSE=sqrt(mean((dd-M$do)^2)),dT_SDrec=sd(dd)/sdo_d,
  sl_R2=cor(ss,M$so,use="complete.obs")^2,sl_RMSE=sqrt(mean((ss-M$so)^2,na.rm=TRUE)),sl_SDrec=sd(ss,na.rm=TRUE)/sdo_s)}
# oracle switch on true LAI median
mL<-median(M$LAI_ALS); denseL<-M$LAI_ALS>=mL
hd_o<-ifelse(denseL,M$dL,M$dS); hs_o<-ifelse(denseL,M$sL,M$sS)
# operational switch on FORMS-H median
mH<-median(M$FORMS_H); denseH<-M$FORMS_H>=mH
hd_h<-ifelse(denseH,M$dL,M$dS); hs_h<-ifelse(denseH,M$sL,M$sS)
# LOO FORMS-H threshold (leave-one-out median of the rest) -> robustness to overfit
denseH_loo<-sapply(seq_len(nrow(M)),function(i){thr<-median(M$FORMS_H[-i]);M$FORMS_H[i]>=thr})
hd_l<-ifelse(denseH_loo,M$dL,M$dS); hs_l<-ifelse(denseH_loo,M$sL,M$sS)
out<-rbind(
  mets(M$dL,M$sL,"LiDAR full (pooled)"),
  mets(M$dS,M$sS,"S2 ATBD (pooled)"),
  mets(hd_o,hs_o,"Hybrid switch: LAI-oracle"),
  mets(hd_h,hs_h,"Hybrid switch: FORMS-H"),
  mets(hd_l,hs_l,"Hybrid switch: FORMS-H (LOO thr)"))
cat(sprintf("LAI median=%.2f | FORMS-H median=%.1f | agreement of the two splits: %.0f%%\n",
    mL,mH,100*mean(denseL==denseH)))
print(out,digits=3)
fwrite(out,file.path(TAB,"figK_hybrid_switch.csv")); cat("\nwrote figK_hybrid_switch.csv\n")
