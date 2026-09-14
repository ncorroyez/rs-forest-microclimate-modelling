# ==============================================================================
# Regime-switching hybrid with Ch2-corrected sensor products per regime.
#   Dense (saturated) -> LiDAR full (the reference; Ch2 has no LiDAR correction).
#   Open              -> S2, tried in 3 Ch2 flavours: raw ATBD, d_opt-corrected,
#                        rescaled-to-LiDAR-peak. Which S2 correction helps in open?
# Switch rules: LAI-absolute (split at 3.86, transferable oracle) and
# FORMS-H-absolute (operational; threshold = height that reproduces the LAI split).
#   Rscript c3_hybrid_ch2_wc.R  -> chapter3_S2_LAI/tables/figK_hybrid_ch2.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"
CROSS<-3.86   # regime split, chosen: near the median LiDAR LAI of the 53 plots (3.80),
              # inside the 3.0-4.5 plateau where the stratified reversal is insensitive
              # (Fig. J2). NOT a Ch2-derived crossover: no such value exists in the RSE
              # manuscript, and no derivation here reproduces it. Do not re-label.
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
grab<-function(scn,suff){P<-getpp(scn);setnames(P,c("ds","ss"),paste0(c("d","s"),suff));P}
M<-Reduce(function(a,b)merge(a,b,by="id_plot"),list(OBS,
  grab("STATIC_ALS","L"),grab("STATIC_S2_ATBD","Sa"),grab("STATIC_S2_DOPT","Sd"),grab("STATIC_S2_RESCALED","Sr")))
M<-merge(M,df[,.(id_plot,LAI_ALS,FORMS_H)],by="id_plot")
sdo_d<-sd(M$do); sdo_s<-sd(M$so,na.rm=TRUE)
mets<-function(dd,ss,tag){data.table(product=tag,
  dT_R2=cor(dd,M$do)^2,dT_bias=mean(dd-M$do),dT_RMSE=sqrt(mean((dd-M$do)^2)),dT_SDrec=sd(dd)/sdo_d,
  sl_R2=cor(ss,M$so,use="complete.obs")^2,sl_SDrec=sd(ss,na.rm=TRUE)/sdo_s)}
# absolute Ch2 switch on LAI; FORMS-H threshold that best reproduces it
denseL<-M$LAI_ALS>=CROSS
# pick FORMS-H abs threshold maximizing agreement with the LAI-3.86 split
cand<-sort(unique(M$FORMS_H)); agr<-sapply(cand,function(h)mean((M$FORMS_H>=h)==denseL)); hthr<-cand[which.max(agr)]
denseH<-M$FORMS_H>=hthr
hyb<-function(dsel,ssel,dense){data.table(dd=ifelse(dense,M$dL,dsel),ss=ifelse(dense,M$sL,ssel))}
mk<-function(open_d,open_s,dense,lbl){h<-hyb(open_d,open_s,dense);mets(h$dd,h$ss,lbl)}
out<-rbind(
  mets(M$dL,M$sL,"LiDAR full (pooled)"),
  mets(M$dSa,M$sSa,"S2 ATBD (pooled)"),
  mk(M$dSa,M$sSa,denseL,"Hybrid LAI-abs: LiDAR+S2 ATBD"),
  mk(M$dSd,M$sSd,denseL,"Hybrid LAI-abs: LiDAR+S2 d_opt"),
  mk(M$dSr,M$sSr,denseL,"Hybrid LAI-abs: LiDAR+S2 rescaled"),
  mk(M$dSa,M$sSa,denseH,"Hybrid FORMS-H: LiDAR+S2 ATBD"),
  mk(M$dSd,M$sSd,denseH,"Hybrid FORMS-H: LiDAR+S2 d_opt"),
  mk(M$dSr,M$sSr,denseH,"Hybrid FORMS-H: LiDAR+S2 rescaled"))
cat(sprintf("regime split LAI=%.2f -> %d dense / %d open | FORMS-H abs thr=%.1f m (agreement %.0f%%)\n",
    CROSS,sum(denseL),sum(!denseL),hthr,100*max(agr)))
print(out,digits=3)
fwrite(out,file.path(TAB,"figK_hybrid_ch2.csv")); cat("\nwrote figK_hybrid_ch2.csv\n")
