# ==============================================================================
# Chapter 3 robustness control — re-run the 8 monthly-slope scenarios under
# MuSICA v3.2.3 + ABL_flag='iter' (yoyo, real MERRA-2 PBLH), to test whether the
# v3.2.0 headline (LiDAR collapse at shoulders / fusion wins April; LiDAR wins
# the summer plateau) survives the ABL-coupled binary. The c1 benchmark showed
# v3.2.3 over-smooths the between-plot field (ΔTmax r 0.93->0.63); since the Ch3
# headline IS a between-plot ranking, expect absolute R² to drop for ALL
# scenarios — read the result on the summer-plateau survival + relative ordering.
#
# Scenarios reconstructed IDENTICALLY to the v3.2.0 production
# (c3_summer_dyn_notmasked.R + run_fusion_variants.R + R/scenarios_c3.R); the
# ONLY changes are: binary (3.2.3), forcing (PBLH), abl_flag='iter', and a
# 2020/366 phenology seed row (iter reads the day before 2021-01-01). Magnitudes
# are untouched: run_musica_one doubles the phenology Leaf_area and calc_phenology
# is linear in LAI_max_per_cohort, so a one-sided parametric curve doubled == the
# existing STATIC fallback (LAI_max=2*LAI, no doubling).
#   Rscript run_v323iter_ch3_full.R   (writes out_files/Chapter3/nc_v323iter/<sc>/)
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC))  # RELATIVE (staged as ../forcing.nc)
ABL  <- list("abl_flag" = '"iter"')
LY   <- 2020:2022
NCROOT <- file.path(CFG_C3$out_dir, "nc_v323iter"); dir.create(NCROOT, recursive=TRUE, showWarnings=FALSE)
NCORES <- 10L

# ---- seed-row injection (the only phenology change vs v3.2.0) -----------------
inject_seed <- function(ph) {
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]
  if (nrow(d) != 1L) stop("missing 2020/365 row to seed 2020/366")
  d$Julian_day <- 366L; rbind(ph, d)
}
mk_dyn  <- function(ts) function(pr) { pid <- plot_id_from_row(pr)
  if (!pid %in% names(ts)) return(NULL); inject_seed(make_phenology_from_s2(ts[[pid]], LY)) }
mk_stat <- function(col) function(pr) inject_seed(calc_phenology(list.year = LY, nleafage = 1,
  budburst_date = 115, leaf_age_max_in = 0.56, relative_age_firstmax = 0.10,
  relative_age_lastmax = 0.75, LAI_max_per_cohort = as.numeric(pr[[col]])))

# ---- data + magnitudes (identical to production) -----------------------------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE)/mean(df$LAI_S2_ATBD, na.rm=TRUE)
df[is.na(LAI_S2_DOPT), LAI_S2_DOPT := LAI_S2_ATBD*ropt]
df[, LAI_A_OPTCAL := LAI_ALS * LAI_S2_ATBD / LAI_S2_DOPT]
df[is.na(LAI_A_OPTCAL), LAI_A_OPTCAL := LAI_ALS]            # 6 near-bare plots fall back to LiDAR full
rA  <- mean(df$LAI_ALS)/mean(df$LAI_S2_ATBD)
bod <- mean(df$LAI_ALS_DOPT)/mean(df$LAI_S2_DOPT)
ropt_p <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
fp     <- setNames(df$LAI_ALS/df$LAI_S2_ATBD, df$pid)
fA     <- setNames(df$LAI_ALS/df$LAI_S2_DOPT, df$pid)

# ---- ATBD daily series from NOT_MASKED rasters (identical to production) ------
nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk <- rast(files); dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[, .(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[, -1, drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
cat(sprintf("ATBD series: %d/53 plots, %d dates\n", length(atbd), length(dates)))

mkser <- function(scale_fun, fb_col) setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]; a<-atbd[[p]]
  if (is.null(a)) return(data.frame(doy=1:365, lai=rep(df[[fb_col]][i],365)))
  data.frame(doy=a$doy, lai=a$lai*scale_fun(i,p)) }), df$pid)
ts_atbd   <- atbd
ts_opt    <- mkser(function(i,p) ropt_p[[p]],     "LAI_S2_DOPT")
ts_rfull  <- mkser(function(i,p) rA,              "LAI_ALS")
ts_ordopt <- mkser(function(i,p) ropt_p[[p]]*bod, "LAI_ALS_DOPT")
ts_fp     <- mkser(function(i,p) fp[[p]],         "LAI_ALS")
ts_fA     <- mkser(function(i,p) fA[[p]],         "LAI_A_OPTCAL")

# ---- the 8 scenarios ---------------------------------------------------------
.lf <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
SC <- list(
  STATIC_ALS              = list(name="STATIC_ALS",              lai_fn=.lf("LAI_ALS"),      phenology_fn=mk_stat("LAI_ALS")),
  STATIC_ALS_DOPT         = list(name="STATIC_ALS_DOPT",         lai_fn=.lf("LAI_ALS_DOPT"), phenology_fn=mk_stat("LAI_ALS_DOPT")),
  DYN_S2_ATBD_NM          = list(name="DYN_S2_ATBD_NM",          lai_fn=.lf("LAI_S2_ATBD"),  phenology_fn=mk_dyn(ts_atbd)),
  DYN_S2_OPT_NM           = list(name="DYN_S2_OPT_NM",           lai_fn=.lf("LAI_S2_DOPT"),  phenology_fn=mk_dyn(ts_opt)),
  DYN_ATBD_rfull_NM       = list(name="DYN_ATBD_rfull_NM",       lai_fn=.lf("LAI_ALS"),      phenology_fn=mk_dyn(ts_rfull)),
  DYN_OPT_rdopt_NM        = list(name="DYN_OPT_rdopt_NM",        lai_fn=.lf("LAI_S2_DOPT"),  phenology_fn=mk_dyn(ts_ordopt)),
  DYN_ALS_S2TIMING_NM     = list(name="DYN_ALS_S2TIMING_NM",     lai_fn=.lf("LAI_ALS"),      phenology_fn=mk_dyn(ts_fp)),
  DYN_ALSoOPT_S2TIMING_NM = list(name="DYN_ALSoOPT_S2TIMING_NM", lai_fn=.lf("LAI_A_OPTCAL"), phenology_fn=mk_dyn(ts_fA)))
for (nm in names(SC)) { SC[[nm]]$hmax_fn <- .lf("Hmax"); SC[[nm]]$fcover_fn <- .lf("fCover"); SC[[nm]]$lad_fn <- lad }

# ---- flattened task list (scenario x plot), resumable ------------------------
tasks <- rbindlist(lapply(names(SC), function(s) data.table(scn=s, i=seq_len(nrow(df)))))
for (s in names(SC)) dir.create(file.path(NCROOT, s), recursive=TRUE, showWarnings=FALSE)
cat(sprintf("launching %d runs (%d scenarios x %d plots) on %d cores\n", nrow(tasks), length(SC), nrow(df), NCORES))

run_one <- function(k) {
  s <- tasks$scn[k]; i <- tasks$i[k]; prow <- as.data.frame(df[i])
  nc <- file.path(NCROOT, s, sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return(sprintf("skip %s/%s", s, prow$id_plot))
  if (file.exists(nc)) file.remove(nc)   # remove a previous failed stub
  ok <- tryCatch({ run_musica_one(prow, SC[[s]], nc, FORC, BIN, extra_setup = ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error = function(e) FALSE)
  sprintf("%s %s/%s", if (ok) "OK" else "FAIL", s, prow$id_plot)
}
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(tasks)), run_one, mc.cores = NCORES, mc.preschedule = FALSE)
res <- unlist(res)
dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))
tab <- table(sub(" .*", "", res))
cat(sprintf("\nDONE in %.1f min | %s\n", dt, paste(names(tab), tab, sep="=", collapse=" ")))
fails <- grep("^FAIL", res, value=TRUE)
if (length(fails)) { cat("FAILURES:\n"); cat(fails, sep="\n"); cat("\n") }
saveRDS(res, file.path(NCROOT, "_run_status.rds"))
cat("status saved -> _run_status.rds\n")
