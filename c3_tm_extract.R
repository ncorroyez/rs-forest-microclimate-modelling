# ==============================================================================
# Chapter 3 — HARMONIZATION with Chapter 1's canonical ΔTmax convention.
# Re-extract per-plot ΔTmax (and night ΔTmin) for OBS + every windcorr scenario
# under BOTH conventions, from existing NetCDFs (no re-simulation):
#   old : each series takes its OWN daily max (Ch3 legacy, c3_figset_data_wc.R)
#   new : TIME-MATCHED (Ch1, R/dtmax_convention.R): sub-canopy read at the hour
#         of the MACRO daily max. Clock: the pblh macro forcing is on the SAME
#         clock as the HOBO loggers (median daily-max hour 14 h for both; hourly
#         cross-correlation peak at 0/+1 h = physical understory lag), so the
#         Ch1 −1 h logger repair does NOT apply here: offset 0 for obs AND sims.
# Night analog: ΔTmin at the hour of the MACRO daily minimum (new) vs own-min (old).
# Out: out_files/Chapter3/tables/perplot_dtmax_conventions.csv
#   Rscript c3_tm_extract.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr)
  library(data.table);library(parallel);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source));source("Chapter3_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
FORC<-"in_files/musica_in_Blois_pblh.nc"
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
OUT<-file.path(CFG_C3$out_dir,"tables"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

# ---- macro references (max and min), forcing clock ---------------------------
MREF<-macro_ref(FORC,ds)                                   # date, Tmax_macro, t_max
macro_ref_min<-function(forcing_nc,dates,varname="Tair"){
  nc<-nc_open(forcing_nc);tu<-ncatt_get(nc,"time","units")$value
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC");th<-ncvar_get(nc,"time")
  Tm<-as.numeric(ncvar_get(nc,varname))-273.15;nc_close(nc)
  d<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=Tm)
  d<-d[as.Date(time)%in%dates]
  d[,.(Tmin_macro=min(Tmac,na.rm=TRUE),t_min=time[which.min(Tmac)]),by=.(date=as.Date(time))]}
MREFN<-macro_ref_min(FORC,ds)
# time-matched ΔTmin per series
delta_tmin_mean<-function(micro,mrefn,dates=NULL,min_days=30){
  m<-as.data.table(micro)[,.(time=floor_date(time,"hour"),Tmic)]
  m<-m[,.(Tmic=mean(Tmic,na.rm=TRUE)),by=time][,date:=as.Date(time)]
  j<-merge(m,mrefn,by="date")[time==t_min]
  if(!is.null(dates)) j<-j[date%in%dates]
  if(nrow(j)<min_days) return(NA_real_)
  mean(j$Tmic-j$Tmin_macro,na.rm=TRUE)}
# old-convention daily macro max / min (identical to the c3_* production recipe)
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds]
mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE),Tmn=min(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]

# ---- one series -> 4 metrics --------------------------------------------------
met4<-function(m){ # m: data.table(time,Tmic) hourly-able
  mm<-as.data.table(m)[,.(time=floor_date(time,"hour"),Tmic)]
  mm<-mm[,.(Tmic=mean(Tmic,na.rm=TRUE)),by=time][as.Date(time)%in%ds]
  dd<-mm[,.(Tmx_s=max(Tmic,na.rm=TRUE),Tmn_s=min(Tmic,na.rm=TRUE)),by=.(date=as.Date(time))]
  dd<-merge(dd,mday,by="date")
  list(d_old=mean(dd$Tmx_s-dd$Tmx,na.rm=TRUE),
       n_old=mean(dd$Tmn_s-dd$Tmn,na.rm=TRUE),
       d_new=delta_tmax_mean(mm,MREF,ds),
       n_new=delta_tmin_mean(mm,MREFN,ds))}

# ---- observations --------------------------------------------------------------
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(datetime)%in%ds]
hb[,time:=floor_date(datetime,"hour")]
HO<-hb[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot=as.character(id_plot),time)]
OBS<-HO[,{r<-met4(.SD[,.(time,Tmic)]);as.data.table(r)},by=id_plot][,scenario:="OBS"]
cat(sprintf("OBS done: n=%d | mean d_old=%+.3f d_new=%+.3f n_old=%+.3f n_new=%+.3f\n",
  nrow(OBS),mean(OBS$d_old),mean(OBS$d_new),mean(OBS$n_old),mean(OBS$n_new)))

# ---- scenarios ------------------------------------------------------------------
SCN<-list.dirs(WC,recursive=FALSE,full.names=FALSE)
SCN<-SCN[sapply(SCN,function(s)length(list.files(file.path(WC,s),pattern="\\.nc$"))==53)]
cat(sprintf("%d scenarios\n",length(SCN)))
one_scn<-function(scn){
  fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
  r<-rbindlist(lapply(fs,function(f){
    id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)")
    nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    tr<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc)
    if(is.null(tr)||!nrow(tr))return(NULL)
    m<-as.data.table(tr)[,.(time,Tmic=Tair_sim)]
    cbind(data.table(id_plot=id),as.data.table(met4(m)))}),fill=TRUE)
  r[,scenario:=scn]
  cat(sprintf("  %s done (n=%d)\n",scn,nrow(r))); r}
ALL<-rbindlist(mclapply(SCN,one_scn,mc.cores=4),fill=TRUE)
ALL<-rbind(ALL,OBS,fill=TRUE)
setcolorder(ALL,c("scenario","id_plot"))
fwrite(ALL,file.path(OUT,"perplot_dtmax_conventions.csv"))
cat(sprintf("wrote %s (%d rows)\n",file.path(OUT,"perplot_dtmax_conventions.csv"),nrow(ALL)))
