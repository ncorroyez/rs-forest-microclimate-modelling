# ==============================================================================
# cmp_gedi_matrix.R — Chapter 4: the sensor-combination MATRIX, summer (JJAS,
# primary) and autumn (context). {ALS, S2, GEDI} x {static, temporal}, all on
# the same genuine-53 windcorr machinery, HOBO-validated (n = 53).
# Writes NC_Full/manuscripts/ch4/tables/Table20_c4_matrix.csv.
#   Rscript cmp_gedi_matrix.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
WINDOWS<-list(summer=seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"),
              autumn=seq(as.Date("2021-10-01"),as.Date("2021-11-30"),by="day"))

FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MAall<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
hball<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hball[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hball<-hball[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")

# scenario -> (sensor combo, static/temporal) mapping
SC<-rbindlist(list(
  list("CONST_ALS",          "ALS",        "constant"),
  list("STATIC_ALS",         "ALS",        "static"),
  list("STATIC_S2_ATBD",     "S2",         "static"),
  list("DYN_S2_ATBD",        "S2",         "temporal"),
  list("STATIC_GEDI_DIRECT", "GEDI",       "static"),
  list("DYN_GEDIONLY",       "GEDI",       "temporal"),
  list("STATIC_S2_RESCALED", "ALS+S2",     "static"),
  list("DYN_S2_ANNUAL_FIX",  "ALS+S2",     "temporal"),
  list("FUSION_H",           "ALS+S2",     "static (layered)"),
  list("DYN_ALS_GEDIONLY",   "ALS+GEDI",   "temporal"),
  list("STATIC_GEDI_RATIO",  "S2+GEDI",    "static"),
  list("STATIC_GEDI_RF",     "S2+GEDI",    "static (RF)"),
  list("DYN_S2GEDI",         "S2+GEDI",    "temporal"),
  list("NAIVE_S2_FORMSH",    "S2+height",  "static"),
  list("DYN_ALS_GEDIMAX",    "ALS+S2+GEDI","temporal")))
setnames(SC,c("scenario","combo","type"))

run_window<-function(ds){
  MA<-MAall[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")]
  mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
  hb<-hball[as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
  hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
  OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],
             hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
  sdo<-sd(OBS$do)
  getpp<-function(scn){fs<-list.files(file.path(d3,scn),pattern="\\.nc$",full.names=TRUE)
   m<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
    s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))]
    if(!nrow(s))return(NULL)
    sm<-merge(s,MA[,.(hr,Tm)],by="hr")
    sd<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date")
    data.table(id_plot=id,ds=mean(sd$Tmx_s-sd$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
   merge(m,OBS,by="id_plot")}
  rbindlist(lapply(SC$scenario,function(s){d<-getpp(s)
    data.table(scenario=s,n=nrow(d),
      dT_R2=cor(d$ds,d$do)^2, dT_RMSE=sqrt(mean((d$ds-d$do)^2)),
      dT_bias=mean(d$ds-d$do), dT_SDrec=sd(d$ds)/sdo,
      sl_R2=cor(d$ss,d$so,use="complete.obs")^2)}))
}
out<-rbindlist(lapply(names(WINDOWS),function(w) cbind(window=w,run_window(WINDOWS[[w]]))))
out<-merge(SC,out,by="scenario")
num<-names(out)[sapply(out,is.numeric)];out[,(num):=lapply(.SD,round,3),.SDcols=num]
for(w in names(WINDOWS)){cat("\n===",w,"===\n")
  print(out[window==w][order(-dT_R2),.(combo,type,scenario,dT_R2,dT_RMSE,dT_bias,dT_SDrec,sl_R2,n)],nrows=20)}
fwrite(out,file.path(TAB,"Table20_c4_matrix.csv"))
cat("\nwrote",file.path(TAB,"Table20_c4_matrix.csv"),"\n")
