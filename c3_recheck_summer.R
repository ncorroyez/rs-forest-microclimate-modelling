# ==============================================================================
# Summer results re-check (windcorr genuine-53) before writing the Ch3 article.
# Re-derives the headline numbers straight from the nc (source of truth) and
# confirms them against the delivered CSVs. Also resolves the CONST_ALS vs
# DYN_ALS "identical rows" flag (expected in summer? or a bug?).
#   Rscript c3_recheck_summer.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
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
r2<-function(a,b) if(sum(is.finite(a)&is.finite(b))>=4) cor(a,b,use="complete.obs")^2 else NA
grab<-function(s,suf){P<-getpp(s);setnames(P,c("ds","ss"),paste0(c("d","s"),suf));P}
cat(sprintf("OBS: n=%d  SD(ΔTmax)=%.2f range[%.1f,%.1f]  SD(slope)=%.3f\n\n",nrow(OBS),sd(OBS$do),min(OBS$do),max(OBS$do),sd(OBS$so,na.rm=T)))
M<-Reduce(function(a,b)merge(a,b,by="id_plot"),list(OBS,
  grab("STATIC_ALS","L"),grab("STATIC_S2_ATBD","Sa"),grab("CONST_ALS","C"),grab("DYN_ALS","Dy"),
  grab("STATIC_S2_DOPT","Sd"),grab("STATIC_S2_RESCALED","Sr"),grab("STATIC_ALS_DOPT","Lt")))
M<-merge(M,df[,.(id_plot,LAI_ALS,FORMS_H)],by="id_plot")
# RQ1 premise
cat("== RQ1 premise (expect obs slope -0.92, MuSICA slope -0.91, dTmax -0.71) ==\n")
cat(sprintf("  obs slope~LAI r=%.2f | MuSICA(STATIC_ALS) slope r=%.2f | dTmax r=%.2f\n",
  cor(M$so,M$LAI_ALS),cor(M$sL,M$LAI_ALS),cor(M$dL,M$LAI_ALS)))
# RQ2 stratified 3.86
de<-M$LAI_ALS>=CROSS
cat(sprintf("\n== RQ2 stratified @3.86 (%d dense/%d open) ==\n",sum(de),sum(!de)))
cat(sprintf("  DENSE  dTmax R2: LiDAR %.2f  S2 %.2f | slope: LiDAR %.2f S2 %.2f\n",r2(M$dL[de],M$do[de]),r2(M$dSa[de],M$do[de]),r2(M$sL[de],M$so[de]),r2(M$sSa[de],M$so[de])))
cat(sprintf("  OPEN   dTmax R2: LiDAR %.2f  S2 %.2f | slope: LiDAR %.2f S2 %.2f\n",r2(M$dL[!de],M$do[!de]),r2(M$dSa[!de],M$do[!de]),r2(M$sL[!de],M$so[!de]),r2(M$sSa[!de],M$so[!de])))
cat(sprintf("  POOLED dTmax R2: LiDAR %.2f  S2 %.2f\n",r2(M$dL,M$do),r2(M$dSa,M$do)))
# RQ3 magnitude bias
cat("\n== RQ3 magnitude warm bias (expect trunc +1.14, full 0.73, S2 0.79) ==\n")
cat(sprintf("  LiDAR full %.2f | S2 ATBD %.2f | LiDAR trunc d_opt %.2f\n",mean(M$dL-M$do),mean(M$dSa-M$do),mean(M$dLt-M$do)))
# RQ4 fusion
hyb<-function(od){ifelse(de,M$dL,od)};hys<-function(os){ifelse(de,M$sL,os)}
sdo_d<-sd(M$do)
cat("\n== RQ4 regime fusion (expect ATBD 0.573/0.67/1.60/0.262/0.781; d_opt 0.48; resc 0.247) ==\n")
for(nm in c("ATBD","d_opt","rescaled")){od<-switch(nm,ATBD=M$dSa,d_opt=M$dSd,rescaled=M$dSr);os<-switch(nm,ATBD=M$sSa,d_opt=M$sSd,rescaled=M$sSr)
  h<-hyb(od);hs<-hys(os);cat(sprintf("  +S2 %-9s dTmax R2=%.3f bias=%.2f RMSE=%.2f SDrec=%.3f | slope R2=%.3f\n",nm,r2(h,M$do),mean(h-M$do),sqrt(mean((h-M$do)^2)),sd(h)/sdo_d,r2(hs,M$so)))}
cat(sprintf("  baselines: LiDAR dTmax R2=%.3f bias=%.2f | S2 dTmax R2=%.3f bias=%.2f\n",r2(M$dL,M$do),mean(M$dL-M$do),r2(M$dSa,M$do),mean(M$dSa-M$do)))
# ANOMALY: CONST vs DYN_ALS
cat("\n== CONST_ALS vs DYN_ALS (audit flagged identical summer rows) ==\n")
cat(sprintf("  max|ΔTmax diff| per plot = %.5f | max|slope diff| = %.5f\n",max(abs(M$dC-M$dDy)),max(abs(M$sC-M$sDy))))
cat(sprintf("  => %s\n", if(max(abs(M$dC-M$dDy))<1e-4)"IDENTICAL in summer (expected: DYN_ALS at summer plateau = CONST_ALS)" else "DIFFER (rows were coincidence)"))
cat("\nDONE recheck\n")
