# ==============================================================================
# Chapter 3 — COMPLETE scenario ranking at n=53, summer, v3.2.3-iter.
# Consolidates the two v3.2.3 caches (all 53/53 now):
#   Chapter3/nc_v323iter  = make_all_scenarios_c3 matrix (26)
#   Chapter1/nc_v323iter  = the _NM / timing / FUSION_H / unifLAD set (8 unique)
# Per-plot between-plot metric (ΔTmax R² + slope-and-equilibrium R²), 1 m, no
# shift, all summer days + hottest decile. Every scenario now on the full 53 plots.
# Run from z_Example root:  Rscript c3_allscen_summer_v323_full.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(dplyr); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

d3 <- file.path(CFG_C3$out_dir, "nc_genuine53")                # FULLY GENUINE 53/53 (real df + un-masked series)
d1 <- "out_files/Chapter1/nc_v323iter"                          # Chapter1 (_NM set)
# scenario -> dir map: Chapter3 wins name collisions; add Chapter1-only scenarios
map <- setNames(file.path(d3, basename(list.dirs(d3, recursive=FALSE))),
                basename(list.dirs(d3, recursive=FALSE)))
for (s in basename(list.dirs(d1, recursive=FALSE)))
  if (!s %in% names(map)) map[s] <- file.path(d1, s)
cat(sprintf("%d scenarios (%d from Chapter3, %d unique from Chapter1)\n",
            length(map), length(list.dirs(d3, recursive=FALSE)), length(map)-length(list.dirs(d3, recursive=FALSE))))

ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, hr:=floor_date(time,"hour")]
mday <- macro_h[, .(Tmax_macro=max(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]
thr <- quantile(mday$Tmax_macro, 0.90, na.rm=TRUE); hot <- mday[Tmax_macro>=thr]$date

# GENUINE 53/53: every scenario now has real inputs for all 53 plots (real df +
# un-masked series), so no plot exclusion is needed.
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) &
         as.Date(time) %in% ds_all,
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
for(scn in names(map)){
  fs <- list.files(map[[scn]], pattern="\\.nc$", full.names=TRUE); if(!length(fs)) next
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
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4_allscen_summer_real53genuine.csv")
cat("\n=== ALL scenarios, SUMMER 'all', ranked by ΔTmax R² (GENUINE 53/53, real data everywhere) ===\n")
print(R[period=="all"][order(-dT_r2), .(scenario,n,dT_r2,dT_rmse,dT_bias,sl_r2,sl_rmse,sl_bias)])
cat("\nDONE -> Table4_allscen_summer_v323_real45.csv (period='hot' included)\n")
