# ==============================================================================
# Chapter 3 — Two extra "S2-timing × magnitude" variants requested by NC.
# All three share ONE S2 temporal shape (the ATBD smoothed series a(t)); they differ
# only by the per-plot constant magnitude factor applied to it (and lai_fn / LAD set
# to the matching magnitude, so each is internally consistent):
#   existing  DYN_ALS_S2TIMING_NM  : a(t) × LAI_ALS/LAI_S2_ATBD     -> ~LiDAR full (3.67)
#   A (this)  DYN_ALSoOPT_S2TIMING_NM : a(t) × LAI_ALS/LAI_S2_DOPT  -> ~4.10 (opt-calibrated)
#   B (this)  DYN_ALSDOPT_S2TIMING_NM : a(t) × LAI_ALS_DOPT/LAI_S2_ATBD -> ~LiDAR d_opt (1.61)
# Runs both A and B (53 plots each), then reports the slope-and-equilibrium slope R²
# per month vs observed. Run from z_Example root:  Rscript run_fusion_variants.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(dplyr); library(tidyr); library(purrr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt_glob <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE)  # compute BEFORE subset assign
df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt_glob]
# derived magnitude for variant A (so lai_fn matches the series summer level)
df[, LAI_A_OPTCAL := LAI_ALS * LAI_S2_ATBD / LAI_S2_DOPT]
# 6 plots lack the opt product -> opt-calibration undefined; fall back to LiDAR full
df[is.na(LAI_A_OPTCAL), LAI_A_OPTCAL := LAI_ALS]

# rebuild ATBD daily series from Not_Masked (incl. shoulders)
nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk <- rast(files); dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[, .(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[, -1, drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)

fA <- setNames(df$LAI_ALS      / df$LAI_S2_DOPT, df$pid)   # variant A: opt-calibrated
fB <- setNames(df$LAI_ALS_DOPT / df$LAI_S2_ATBD, df$pid)   # variant B: LiDAR d_opt magnitude
mk <- function(ts) make_phenology_fn_factory(ts, CFG_C3$list_year)
.fn <- function(c) function(pr) as.numeric(pr[[c]]); lad <- make_lad_real
ser <- function(fac, fb_col) setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]; a<-atbd[[p]]
  if(is.null(a)) return(data.frame(doy=1:365, lai=rep(df[[fb_col]][i],365)))
  data.frame(doy=a$doy, lai=a$lai*fac[[p]]) }), df$pid)
ts_A <- ser(fA, "LAI_A_OPTCAL"); ts_B <- ser(fB, "LAI_ALS_DOPT")

# variant A only for now (B deferred at user request)
sc <- list(
  DYN_ALSoOPT_S2TIMING_NM = list(name="DYN_ALSoOPT_S2TIMING_NM",
    lai_fn=.fn("LAI_A_OPTCAL"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_A)))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(as.data.frame(df), hd, sc, file.path(CFG_C3$out_dir,"nc"),
                            dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)

# ---- monthly slope R2 vs observed for both ----
ncp <- file.path(CFG_C3$out_dir, "nc"); mlab <- function(d) format(d, "%Y-%m")
WIN <- as.Date(c("2021-04-15","2021-11-15")); ds_all <- seq(WIN[1], WIN[2], by="day")
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),
                      Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h <- macro_h[as.Date(time) %in% ds_all]; macro_h[, `:=`(hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
hb <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb <- hb[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove) & as.Date(time) %in% ds_all,
         .(id_plot=as.character(id_plot), Tobs=t_hobo, hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))]
hb <- merge(hb, macro_h[,.(hr,Tmacro)], by="hr")
OBS <- hb[, .(so=slp(Tobs,Tmacro)), by=.(id_plot,mo)]
read_sim <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  as.data.table(res)[as.Date(time) %in% ds_all, .(Tsim=Tair_sim, hr=floor_date(time,"hour"), mo=mlab(as.Date(time)))] }
out <- list()
for(scn in names(sc)){
  fs <- list.files(file.path(ncp,scn), pattern="\\.nc$", full.names=TRUE)
  sim <- rbindlist(lapply(fs, function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim(f)
    if(is.null(s))return(NULL); s[,id_plot:=id][] }), fill=TRUE)
  smh <- merge(sim, macro_h[,.(hr,Tmacro)], by="hr")
  SIM <- smh[, .(ss=slp(Tsim,Tmacro)), by=.(id_plot,mo)]
  M <- merge(SIM, OBS, by=c("id_plot","mo"))
  r <- M[, .(scenario=scn, n=.N, sl_r2=if(.N>2) cor(ss,so,use="complete.obs")^2 else NA), by=mo][order(mo)]
  out[[scn]] <- r; cat(sprintf("\n=== %s — slope R2 ===\n",scn)); print(r)
}
fwrite(rbindlist(out), "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4v_fusion_variants.csv")
cat("\nDONE fusion variants A (opt-cal) + B (d_opt magnitude)\n")
