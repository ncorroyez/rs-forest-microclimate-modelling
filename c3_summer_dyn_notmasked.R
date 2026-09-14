# ==============================================================================
# Chapter 3 — rebuild the per-plot DYNAMIC S2 LAI series from the NOT_MASKED daily
# rasters (15 dates, data over the whole site) instead of the Deciduous_Only mask
# (which dropped 8 open/edge HOBO plots). Whittaker/GAM smoothing, min_obs=3 so all
# 53 plots get a series. Re-runs the 4 dynamic S2 scenarios under "_NM" names.
# LiDAR (STATIC_ALS, STATIC_ALS_DOPT) is unchanged and reused.
# Run from z_Example root:  Rscript c3_summer_dyn_notmasked.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots); dopt <- CFG_C3$d_opt_m
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt <- mean(df$LAI_S2_DOPT,na.rm=TRUE)/mean(df$LAI_S2_ATBD,na.rm=TRUE); df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt]
rA <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD); bod <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT)
ropt_p <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)

# ---- rebuild ATBD daily series from NOT_MASKED rasters (all 53 plots) ----
nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk <- rast(files); dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[, .(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[, -1, drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)        # named by pid
cat(sprintf("rebuilt ATBD series: %d / 53 plots (Not_Masked, %d dates)\n", length(atbd), length(dates)))

mk <- function(ts) make_phenology_fn_factory(ts, CFG_C3$list_year); .fn <- function(c) function(pr) as.numeric(pr[[c]]); lad <- make_lad_real
mkser <- function(scale_fun, fb_col){ setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]; a<-atbd[[p]]
  if(is.null(a)) return(data.frame(doy=1:365, lai=rep(df[[fb_col]][i],365)))
  data.frame(doy=a$doy, lai=a$lai*scale_fun(i,p)) }), df$pid) }
ts_atbd  <- atbd
ts_opt   <- mkser(function(i,p) ropt_p[[p]],      "LAI_S2_DOPT")
ts_rfull <- mkser(function(i,p) rA,               "LAI_ALS")
ts_ordopt<- mkser(function(i,p) ropt_p[[p]]*bod,  "LAI_ALS_DOPT")

sc <- list(
  DYN_S2_ATBD_NM   = list(name="DYN_S2_ATBD_NM",   lai_fn=.fn("LAI_S2_ATBD"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_atbd)),
  DYN_S2_OPT_NM    = list(name="DYN_S2_OPT_NM",    lai_fn=.fn("LAI_S2_DOPT"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_opt)),
  DYN_ATBD_rfull_NM= list(name="DYN_ATBD_rfull_NM",lai_fn=.fn("LAI_ALS"),     hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_rfull)),
  DYN_OPT_rdopt_NM = list(name="DYN_OPT_rdopt_NM", lai_fn=.fn("LAI_S2_DOPT"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_ordopt)))
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds); hd <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dm, CFG_C3$ids_to_remove)
validate_scenarios_at_hobos(as.data.frame(df), hd, sc, file.path(CFG_C3$out_dir,"nc"), dm, ds, CFG_C3$forcing_file, CFG_C3$musica_cmd, FALSE)
cat("\nDONE Not_Masked dynamic runs\n")
