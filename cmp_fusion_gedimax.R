# ==============================================================================
# cmp_fusion_gedimax.R — Chapter 4: does the GEDI-constrained autumn limb improve the
# microclimate? DYN_ALS_GEDIFALL vs DYN_S2_ANNUAL (identical except DOY >= 250)
# + CONST_ALS / STATIC_ALS / DYN_S2_ATBD context, on two windows:
#   summer JJAS 2021 (should be ~unchanged: control) and autumn Oct-Nov 2021
#   (where the correction lives). Paired bootstrap (B = 3000) on ΔRMSE(ΔTmax)
#   and plot-resampled ΔR² for GEDIFALL - ANNUAL. Contrast reported on the 47
#   plots that carry a real annual series (fallback plots excluded from the
#   paired test, kept in the context table).
# Writes NC_Full/manuscripts/ch4/tables/Table22_c4_fusion_gedimax.csv (+ _boot.csv).
#   Rscript cmp_fusion_gedimax.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
WINDOWS<-list(summer=seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"),
              autumn=seq(as.Date("2021-10-01"),as.Date("2021-11-30"),by="day"))
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
tb<-readRDS(file.path(CFG_C3$out_dir,"lai_prep","ts_by_plot.rds"))
with_series<-as.character(df[pid%in%names(tb$annual),id_plot])

FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MAall<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
hball<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hball[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hball<-hball[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
SC<-c("FUSION_GEDIMAX","FUSION_H","DYN_ALS_GEDIMAX","STATIC_ALS","STATIC_S2_ATBD")

run_window<-function(ds){
  MA<-MAall[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")]
  mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
  hb<-hball[as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
  hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
  OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],
             hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
  getpp<-function(scn){fs<-list.files(file.path(d3,scn),pattern="\\.nc$",full.names=TRUE)
   m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
    s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))]
    if(!nrow(s))return(NULL)
    sm<-merge(s,MA[,.(hr,Tm)],by="hr")
    sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date")
    data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
   merge(m,OBS,by="id_plot")}
  setNames(lapply(SC,getpp),SC)
}
met<-function(d) data.table(n=nrow(d),
  dT_R2=cor(d$ds,d$do)^2, dT_RMSE=sqrt(mean((d$ds-d$do)^2)), dT_bias=mean(d$ds-d$do),
  sl_R2=cor(d$ss,d$so,use="complete.obs")^2)

set.seed(42); B<-3000L
out<-list();boot<-list()
for(w in names(WINDOWS)){
  PP<-run_window(WINDOWS[[w]])
  out[[w]]<-rbindlist(lapply(names(PP),function(s)cbind(window=w,scenario=s,met(PP[[s]]))))
  # paired contrast on plots with a real annual series
  a<-PP$FUSION_GEDIMAX[id_plot%in%with_series];b<-PP[["FUSION_H"]][id_plot%in%with_series]
  m<-merge(a[,.(id_plot,ds_a=ds,ss_a=ss,do,so)],b[,.(id_plot,ds_b=ds,ss_b=ss)],by="id_plot")
  drmse<-replicate(B,{i<-sample(nrow(m),replace=TRUE)
    sqrt(mean((m$ds_a[i]-m$do[i])^2))-sqrt(mean((m$ds_b[i]-m$do[i])^2))})
  dr2  <-replicate(B,{i<-sample(nrow(m),replace=TRUE)
    cor(m$ds_a[i],m$do[i])^2-cor(m$ds_b[i],m$do[i])^2})
  boot[[w]]<-data.table(window=w,n_paired=nrow(m),
    dRMSE=mean(drmse),dRMSE_lo=quantile(drmse,.025),dRMSE_hi=quantile(drmse,.975),
    dR2=mean(dr2),dR2_lo=quantile(dr2,.025),dR2_hi=quantile(dr2,.975))
}
res<-rbindlist(out);num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
bt<-rbindlist(boot);num<-names(bt)[sapply(bt,is.numeric)];bt[,(num):=lapply(.SD,round,3),.SDcols=num]
for(w in names(WINDOWS)){cat("\n===",w,"===\n");print(res[window==w][order(-dT_R2)])}
cat("\n=== paired bootstrap FUSION_GEDIMAX - FUSION_H (47 series plots) ===\n");print(bt)
fwrite(res,file.path(TAB,"Table22_c4_fusion_gedimax.csv"))
fwrite(bt,file.path(TAB,"Table22b_c4_fusion_gedimax_boot.csv"))
cat("\nwrote Table15 / Table15b\n")
