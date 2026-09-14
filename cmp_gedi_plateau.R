# ==============================================================================
# cmp_gedi_plateau.R — Chapter 4, "et l'été ?": month-by-month summer
# decomposition of the GEDI-relevant temporal effects (no new sims).
# Two paired contrasts per month (Jun/Jul/Aug/Sep 2021 + JJAS), 47 series plots:
#   (a) RAMP  = DYN_S2_ANNUAL_FIX - CONST_ALS : is the real June 2021 leaf-up
#       ramp (S2-timed, GEDI-corroborated) worth anything vs a flat plateau?
#   (b) GEDIFALL = DYN_ALS_GEDIFALL - DYN_S2_ANNUAL_FIX : where in the summer
#       does the GEDI autumn-limb correction act (expected: September only)?
# Writes NC_Full/manuscripts/ch4/tables/Table18_c4_plateau.csv (+ _boot).
#   Rscript cmp_gedi_plateau.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
WINDOWS<-list(autumn=seq(as.Date("2021-10-01"),as.Date("2021-11-30"),by="day"),jun=seq(as.Date("2021-06-01"),as.Date("2021-06-30"),by="day"),
              jul=seq(as.Date("2021-07-01"),as.Date("2021-07-31"),by="day"),
              aug=seq(as.Date("2021-08-01"),as.Date("2021-08-31"),by="day"),
              sep=seq(as.Date("2021-09-01"),as.Date("2021-09-30"),by="day"),
              jjas=seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"))
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
SC<-c("DYN_ALS_GEDIPLATEAU","DYN_ALS_GEDIFALL","DYN_S2_ANNUAL_FIX","CONST_ALS","STATIC_ALS")

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
pboot<-function(PPa,PPb,B=3000L){
  m<-merge(PPa[id_plot%in%with_series,.(id_plot,ds_a=ds,do)],
           PPb[id_plot%in%with_series,.(id_plot,ds_b=ds)],by="id_plot")
  drmse<-replicate(B,{i<-sample(nrow(m),replace=TRUE)
    sqrt(mean((m$ds_a[i]-m$do[i])^2))-sqrt(mean((m$ds_b[i]-m$do[i])^2))})
  dr2<-replicate(B,{i<-sample(nrow(m),replace=TRUE)
    cor(m$ds_a[i],m$do[i])^2-cor(m$ds_b[i],m$do[i])^2})
  data.table(n_paired=nrow(m),
    dRMSE=mean(drmse),dRMSE_lo=quantile(drmse,.025),dRMSE_hi=quantile(drmse,.975),
    dR2=mean(dr2),dR2_lo=quantile(dr2,.025),dR2_hi=quantile(dr2,.975))
}
set.seed(42)
out<-list();boot<-list()
for(w in names(WINDOWS)){
  PP<-run_window(WINDOWS[[w]])
  out[[w]]<-rbindlist(lapply(names(PP),function(s)cbind(window=w,scenario=s,met(PP[[s]]))))
  boot[[w]]<-rbind(
    cbind(window=w,contrast="PLATEAU - GEDIFALL (flatten core)",
          pboot(PP[["DYN_ALS_GEDIPLATEAU"]],PP[["DYN_ALS_GEDIFALL"]])),
    cbind(window=w,contrast="PLATEAU - ANNUAL_FIX",
          pboot(PP[["DYN_ALS_GEDIPLATEAU"]],PP[["DYN_S2_ANNUAL_FIX"]])),
    cbind(window=w,contrast="PLATEAU - STATIC_ALS",
          pboot(PP[["DYN_ALS_GEDIPLATEAU"]],PP[["STATIC_ALS"]])))
}
res<-rbindlist(out);num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
bt<-rbindlist(boot);num<-names(bt)[sapply(bt,is.numeric)];bt[,(num):=lapply(.SD,round,3),.SDcols=num]
for(w in names(WINDOWS)){cat("\n===",w,"===\n");print(res[window==w][order(-dT_R2)])}
cat("\n=== paired bootstraps (47 series plots) ===\n");print(bt)
fwrite(res,file.path(TAB,"Table18_c4_plateau.csv"))
fwrite(bt,file.path(TAB,"Table18b_c4_plateau_boot.csv"))
cat("\nwrote Table16 / Table16b\n")
