# ==============================================================================
# Hot-day re-extraction: DeltaTmax and micro-macro slope on the hottest 10% of days,
# for the observations and for all 16 forward-inclusion coalitions, on the aligned clock.
#
# Reads : out_files/musica_native20_forward/<bit>/musica_out_HOBO_<id>.nc  (16 coalitions x 53)
#         in_files/FR-Blo_2021_v2.nc / in_files/Blois_data_temperature.csv
# Writes: out_files/Chapter1/hot_extract.rds        (list(obs =, coal =))
#   Rscript c1_hot_extract.R
# ==============================================================================
# Extract HOT-day (top 10%) DeltaTmax + buffering slope: obs (HOBO) per logger,
# and sim per coalition (16 bits) per logger. Feeds the supplementary hot-day
# companions S3/S4/S5 (c1_hot_figures.R).
#
# REBUILT 2026-07-29. The previous version was stale on three counts at once:
#   (1) it read `musica_hobo_windcorr_forward`, not the native20 coalition set
#       that produces every main figure;
#   (2) it took each side's OWN daily maximum (the superseded convention A),
#       instead of the time-matched convention of R/dtmax_convention.R;
#   (3) it applied no clock offset, so observations (UTC) were compared against
#       a forcing on solar time (UTC+1).
# Simulations are now scored on the forcing clock (MREF) and observations one
# hour earlier (MOBS), exactly as the headline validation does.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(parallel)
  src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("pipeline/00_config.R")})
FORC<-"in_files/FR-Blo_2021_v2.nc"; COALD<-"out_files/musica_native20_forward"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1
MREF<-macro_ref(FORC,ds)          # forcing clock -> SIMULATIONS
MOBS<-macro_ref_obs(FORC,ds)      # logger clock  -> OBSERVATIONS
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
mac<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macH<-mac[as.Date(time)%in%ds][,.(Tm=mean(Tm,na.rm=TRUE)),by=time]
macD<-macH[,.(Tmx=max(Tm)),by=.(date=as.Date(time))]
hot<-macD[Tmx>=quantile(Tmx,0.90),date]
cat(sprintf("hot days: %d\n",length(hot)))

# Time-matched DeltaTmax over the hot days only, plus the micro-macro slope.
# `mref` selects the clock: MREF for simulations, MOBS for observations.
#' Hot-day DeltaTmax and micro-macro slope for one hourly series
#' @param mic data.table(time, Tmic), the sub-canopy or logger series
#' @param mref macro reference clock: MREF for simulations, MOBS for observations
#' @return list(dt, sl); sl is NA when fewer than 30 paired hours are available
metr<-function(mic,mref){
  dt<-delta_tmax_mean(mic,mref,hot,min_days=5)
  m<-merge(mic[as.Date(time)%in%hot],macH,by="time")
  list(dt=dt, sl=if(nrow(m)>30) as.numeric(coef(lm(Tmic~Tm,m))[2]) else NA_real_)}

# ---- observations --------------------------------------------------------------
hobo<-as.data.table(read.csv(CFG$hobo_temp_csv));hobo[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hobo<-hobo[position_sensor=="a"&as.Date(datetime)%in%ds&!id_plot%in%CFG$ids_to_remove]
ho<-hobo[,.(Tmic=mean(t_hobo,na.rm=TRUE)),by=.(id_plot,time=floor_date(datetime,"hour"))]
OBS<-ho[,{v<-metr(.SD[,.(time,Tmic)],MOBS);data.table(obs_dt_hot=v$dt,obs_sl_hot=v$sl)},by=id_plot]

# ---- simulations, per coalition -------------------------------------------------
BITS<-PIPE$BITS
ids<-OBS$id_plot                                   # the validated loggers, fixed
jobs<-CJ(bit=BITS,id=ids)
res<-mclapply(seq_len(nrow(jobs)),function(i){
  d<-micro_hourly_at(file.path(COALD,jobs$bit[i],sprintf("musica_out_HOBO_%s.nc",jobs$id[i])),Z)
  if(is.null(d))return(c(NA_real_,NA_real_));v<-metr(d,MREF);c(v$dt,v$sl)},mc.cores=6)
jobs[,`:=`(dt_hot=sapply(res,`[`,1),sl_hot=sapply(res,`[`,2))]
saveRDS(list(obs=OBS,coal=jobs),"out_files/Chapter1/hot_extract.rds")
cat(sprintf("DONE: obs %d loggers, coal %d rows (%d NA)\n",nrow(OBS),nrow(jobs),sum(is.na(jobs$dt_hot))))
