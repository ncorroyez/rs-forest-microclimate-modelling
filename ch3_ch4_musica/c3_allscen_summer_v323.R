# ==============================================================================
# Chapter 3 — ALL scenarios, summer per-plot validation under v3.2.3-iter
# (cache A = out_files/Chapter3/nc, FR-Blo_2021.nc yoyo iter, 2021). Same per-plot
# metric as Table4f (between-plot ΔTmax R² and slope-and-equilibrium R², n=53 or
# 47), 1 m, no shift, all summer days + hottest decile. Ranks the whole matrix.
# Run from z_Example root:  Rscript c3_allscen_summer_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
ncp <- file.path(CFG_C3$out_dir, "nc")                 # cache A (v3.2.3, FR-Blo_2021)
scen_dirs <- sort(basename(list.dirs(ncp, recursive=FALSE)))
ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")

ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, hr:=floor_date(time,"hour")]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]
thr <- quantile(mday$Tmax_macro, 0.90, na.rm=TRUE); hot <- mday[Tmax_macro>=thr]$date

hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo, date=as.Date(time))]
hb[, hr:=floor_date(time,"hour")]; hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
hb_day <- merge(hb[, .(Tmax_obs=max(Tobs,na.rm=TRUE)), by=.(id_plot,date)], mday, by="date")

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  as.data.table(res)[as.Date(time) %in% ds_all, .(time, Tsim=Tair_sim, hr=floor_date(time,"hour"), date=as.Date(time))] }
agg <- function(s,o){ d<-data.table(s=s,o=o)[is.finite(s)&is.finite(o)]
  if(!nrow(d)) return(data.table(n=0,r2=NA,rmse=NA,bias=NA))
  data.table(n=nrow(d), r2=cor(d$s,d$o)^2, rmse=sqrt(mean((d$s-d$o)^2)), bias=mean(d$s-d$o)) }

res <- list()
for(scn in scen_dirs){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  if(!length(fs)) next
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
  if(!nrow(sim)) next
  smh <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
  sday <- merge(sim[, .(Tmax_sim=max(Tsim,na.rm=TRUE)), by=.(id_plot,date)], mday, by="date")
  for(per in c("all","hot")){ dys <- if(per=="all") ds_all else hot
    dT <- merge(sday[date %in% dys, .(ds=mean(Tmax_sim-Tmax_macro,na.rm=TRUE)), by=id_plot],
                hb_day[date %in% dys, .(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)), by=id_plot], by="id_plot")
    sl <- merge(smh[date %in% dys, .(ss=slp(Tsim,Tmacro)), by=id_plot],
                hb[date %in% dys, .(so=slp(Tobs,Tmacro)), by=id_plot], by="id_plot")
    aT<-agg(dT$ds,dT$do); aS<-agg(sl$ss,sl$so)
    res[[length(res)+1]] <- data.table(scenario=scn, period=per, n=aT$n,
      dT_r2=aT$r2, dT_rmse=aT$rmse, dT_bias=aT$bias, sl_r2=aS$r2, sl_rmse=aS$rmse, sl_bias=aS$bias) }
}
R <- rbindlist(res); numc <- setdiff(names(R), c("scenario","period","n")); R[, (numc):=lapply(.SD,round,3),.SDcols=numc]
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4_allscen_summer_v323.csv")
cat(sprintf("\nbinary=v3.2.3-iter | cache=%s | forcing=%s | %d scenarios\n", ncp, CFG_C3$forcing_file, length(scen_dirs)))
cat("\n=== ALL scenarios, SUMMER 'all', ranked by between-plot ΔTmax R² ===\n")
print(R[period=="all"][order(-dT_r2), .(scenario,n,dT_r2,dT_rmse,dT_bias,sl_r2)])
cat("\nDONE -> Table4_allscen_summer_v323.csv (also has period='hot')\n")
