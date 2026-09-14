# ==============================================================================
# Chapter 3 — nighttime buffering ΔTmin. Microclimatic buffering is bidirectional:
# the canopy also keeps nights warmer (ΔTmin = sub-canopy minus macroclimate daily
# MINIMUM, expected positive). We validate the same 6 scenarios on ΔTmin across the
# 53 HOBO over all summer days and the coldest nights (macro Tmin <= p10), to test
# whether the daytime LiDAR>S2 ranking also holds at night. 1 m, no shift.
# Run from z_Example root:  Rscript deltatmin.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt",
          STATIC_ALS="LiDAR full", STATIC_ALS_DOPT="LiDAR d_opt")
ds_all <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")

ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]
mday <- macro_h[, .(Tmin_macro=min(Tmacro,na.rm=TRUE)), by=.(date=as.Date(time))]
thr <- quantile(mday$Tmin_macro, 0.10, na.rm=TRUE); cold <- mday[Tmin_macro<=thr]$date
cat(sprintf("coldest-nights threshold (macro Tmin p10) = %.2f °C ; %d nights\n", thr, length(cold)))

hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, date=as.Date(time))]
hb_day <- merge(hb[, .(Tmin_obs=min(Tobs,na.rm=TRUE)), by=.(id_plot,date)], mday, by="date")
hb_day[, dmin_obs:=Tmin_obs-Tmin_macro]

read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res)[as.Date(time) %in% ds_all, .(date=as.Date(time), Tsim=Tair_sim)]; r }
agg <- function(d){ d<-d[is.finite(ds)&is.finite(do)]
  data.table(n=nrow(d), r2=cor(d$ds,d$do)^2, rmse=sqrt(mean((d$ds-d$do)^2)),
             mae=mean(abs(d$ds-d$do)), bias=mean(d$ds-d$do)) }
res <- list()
for(scn in names(scen)){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f); if(is.null(s))return(NULL)
    sd <- merge(s[, .(Tmin_sim=min(Tsim,na.rm=TRUE)), by=date], mday, by="date")
    sd[, .(id_plot=id, date, dmin_sim=Tmin_sim-Tmin_macro)] }), fill=TRUE)
  d <- merge(sim, hb_day[,.(id_plot,date,dmin_obs)], by=c("id_plot","date"))
  for(per in c("all","cold")){ dd <- if(per=="all") d else d[date %in% cold]
    pp <- dd[, .(ds=mean(dmin_sim,na.rm=TRUE), do=mean(dmin_obs,na.rm=TRUE)), by=id_plot]
    a <- agg(pp); res[[length(res)+1]] <- data.table(scenario=scen[[scn]], period=per, a) }
}
R <- rbindlist(res); numc<-c("r2","rmse","mae","bias"); R[, (numc):=lapply(.SD,round,3), .SDcols=numc]
R[, period:=factor(period, levels=c("all","cold"))]
fwrite(R, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4o_deltatmin.csv")
cat(sprintf("\nObserved mean ΔTmin across 53 plots: %.2f °C (all) / %.2f °C (coldest)\n",
            mean(hb_day$dmin_obs,na.rm=TRUE), mean(hb_day[date %in% cold]$dmin_obs,na.rm=TRUE)))
cat("\n=== ΔTmin (nighttime) sim-vs-obs across 53 plots ===\n")
print(R[order(period,-r2), .(scenario,period,n,r2,rmse,mae,bias)])
cat("\nDONE -> Table4o_deltatmin.csv\n")
