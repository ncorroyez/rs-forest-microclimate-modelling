# ==============================================================================
# Chapter 3 — SIMPLE summer validation on the 53 HOBO plots, TWO metrics ×
# TWO periods, per plot. Reuses the 7 scenarios' MuSICA nc (no re-sim).
#   ΔTmax : daily sub-canopy max offset (1 m, -2 h convention) -> mean per plot.
#   slope : the SLOPE-AND-EQUILIBRIUM slope (Gril et al. 2023) = lm(T_micro ~ T_macro)
#           over the FULL HOURLY series at 1 m with the SAME -2h shift as ΔTmax
#           (consistent everywhere), ALL temperatures. slope < 1 = buffering.
#   periods : all summer days ; hottest days (macro daily Tmax >= p75).
# Per plot we get sim & obs ΔTmax and slope; we validate sim-vs-obs ACROSS the
# 53 plots (R²/RMSE/bias) for each metric × period.
# Run from z_Example root:  Rscript extend_summer_metrics.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(purrr); library(stringr); library(data.table)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dmd  <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds))           # date, Tmax_macro
ncp  <- file.path(CFG_C3$out_dir, "nc")
scen <- c("DYN_S2_ATBD_NM","DYN_S2_OPT_NM","DYN_ATBD_rfull_NM","DYN_OPT_rdopt_NM","STATIC_ALS","STATIC_ALS_DOPT")

# hottest days (by macro daily Tmax)
thr     <- as.numeric(quantile(dmd$Tmax_macro, 0.90, na.rm=TRUE))
hotdays <- dmd[Tmax_macro >= thr, date]
cat(sprintf("hottest-day set: macro Tmax >= %.1f °C (p90, top 10%%) -> %d of %d days\n",
            thr, length(hotdays), nrow(dmd)))

# ---- hourly macro forcing temperature (for the slope) ----
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time = force_utc_nc(CFG_C3$forcing_file, "time"),
                      Tmacro = as.numeric(ncvar_get(ncf, "Tair")) - 273.15)
nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds]
macro_h[, date := as.Date(time)]

# ---- hourly HOBO observations (position "a") ----
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[, time := as.POSIXct(datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds,
         .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
# HOBO is on-the-hour (XX:00), macro/sim are on-the-half-hour (XX:30): join on floored hour
hb[, hr := floor_date(time, "hour")]
macro_hr <- macro_h[, .(hr=floor_date(time, "hour"), Tmacro, date)]
hb <- merge(hb, macro_hr, by="hr")

# daily ΔTmax (sim) per plot per scenario, and observed Δ
dm_daily <- dmd[, .(date, Tmax_macro)]
hd <- as.data.table(read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dmd, CFG_C3$ids_to_remove))
hd[, id_plot := as.character(id_plot)]

slope <- function(y, x) if (sum(is.finite(x)&is.finite(y)) >= 10) as.numeric(coef(lm(y~x))[2]) else NA_real_

# ---- observed per-plot ΔTmax & slope (all / hot) ----
obs_pp <- function() {
  dT <- hd[, .(dT_all=mean(Delta_obs,na.rm=TRUE),
               dT_hot=mean(Delta_obs[date %in% hotdays],na.rm=TRUE)), by=id_plot]
  sl <- hb[, .(sl_all=slope(Tobs,Tmacro),
               sl_hot=slope(Tobs[date %in% hotdays], Tmacro[date %in% hotdays])), by=id_plot]
  merge(dT, sl, by="id_plot")
}
OBS <- obs_pp()

read_sim_1m <- function(f){
  nc <- tryCatch(nc_open(f), error=function(e) NULL); if (is.null(nc)) return(NULL)
  res <- tryCatch(get_tair_at_z(nc, z_target=1.0), error=function(e) NULL); nc_close(nc)
  if (is.null(res) || !nrow(res)) return(NULL)
  r <- as.data.table(res); r[, time := time - lubridate::hours(2)]   # 1 m + same -2h shift as ΔTmax
  r[, .(time, Tsim=Tair_sim)]
}

# ---- simulated per-plot, per scenario ----
sim_pp <- function(scn){
  dir <- file.path(ncp, scn); fs <- list.files(dir, pattern="\\.nc$", full.names=TRUE)
  if (!length(fs)) return(NULL)
  # ΔTmax (daily max, 1m/-2h) via existing extractor
  sim_dT <- as.data.table(extract_deltatmax_hobo_scenario(dir, dmd, ds))
  dT <- sim_dT[, .(dT_all=mean(Delta_Tmax,na.rm=TRUE),
                   dT_hot=mean(Delta_Tmax[date %in% hotdays],na.rm=TRUE)), by=.(id_plot=as.character(id_plot))]
  # slope (hourly, nair==1, all temps)
  sl <- rbindlist(lapply(fs, function(f){
    id <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)")
    s <- read_sim_1m(f); if (is.null(s)||!nrow(s)) return(NULL)
    s <- merge(s, macro_h[, .(time, Tmacro, date)], by="time")
    data.table(id_plot=id, sl_all=slope(s$Tsim,s$Tmacro),
               sl_hot=slope(s$Tsim[s$date %in% hotdays], s$Tmacro[s$date %in% hotdays]))
  }), fill=TRUE)
  merge(dT, sl, by="id_plot")
}

vmetr <- function(sim, period){
  k <- if(period=="all") "all" else "hot"
  d <- merge(sim[, .(id_plot, ds=get(paste0("dT_",k)), ss=get(paste0("sl_",k)))],
             OBS[, .(id_plot, do=get(paste0("dT_",k)), so=get(paste0("sl_",k)))], by="id_plot")
  dd <- d[is.finite(ds)&is.finite(do)]; sd <- d[is.finite(ss)&is.finite(so)]
  data.table(period=period, n_plot=nrow(dd),
    dT_r2=cor(dd$ds,dd$do)^2, dT_rmse=sqrt(mean((dd$ds-dd$do)^2)), dT_bias=mean(dd$ds-dd$do),
    sl_r2=cor(sd$ss,sd$so)^2, sl_rmse=sqrt(mean((sd$ss-sd$so)^2)), sl_bias=mean(sd$ss-sd$so),
    sl_sim=mean(sd$ss), sl_obs=mean(sd$so))
}
res <- rbindlist(lapply(scen, function(scn){
  sim <- sim_pp(scn); if (is.null(sim)) { cat("  [skip]",scn,"\n"); return(NULL) }
  rbind(cbind(scenario=scn, vmetr(sim,"all")), cbind(scenario=scn, vmetr(sim,"hot")))
}), fill=TRUE)
num <- names(res)[sapply(res,is.numeric)]; res[, (num):=lapply(.SD,function(x) round(x,3)), .SDcols=num]
cat("\n=== SUMMER 53-HOBO: ΔTmax (daily) & micro/macro slope (hourly, all temps) × all/hot ===\n")
print(as.data.frame(res[, .(scenario,period,n_plot,dT_r2,dT_rmse,dT_bias,sl_r2,sl_rmse,sl_bias)]), row.names=FALSE)
cat(sprintf("\nObserved HOBO slope (all temps): mean across plots = %.2f (all) / %.2f (hot)\n",
            mean(OBS$sl_all,na.rm=TRUE), mean(OBS$sl_hot,na.rm=TRUE)))
fwrite(res, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4f_summer_2metrics_2periods.csv")
cat("\nDONE -> Table4f_summer_2metrics_2periods.csv\n")
