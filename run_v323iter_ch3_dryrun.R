# ==============================================================================
# DRY RUN — Chapter 3 robustness control: MuSICA v3.2.3 + ABL_flag='iter' (yoyo,
# real MERRA-2 PBLH) on the Ch3 monthly-slope scenarios. This script runs ONE
# plot for ONE dynamic (S2-driven) and ONE static (parametric) scenario, to
# confirm a NetCDF writes with no STOP 16 / STOP 6 and to time a single iter run
# (iter ≈ ×2, itermax=2) before any full fan-out.
#
# iter needs the phenology of the day before sim start (2021-01-01) -> 2020/366
# (2020 is a leap year). calc_phenology only emits 365 d/yr, so we extend
# list_year to 2020:2022 and inject a 2020/366 row (copy of 2020/365 = winter
# leaf-off, the physically correct seed). This mirrors c1's mk_phen but keeps the
# Ch3 magnitudes untouched: run_musica_one() doubles the phenology Leaf_area, and
# calc_phenology is linear in LAI_max_per_cohort, so a one-sided parametric curve
# doubled == the existing v3.2.0 STATIC fallback (LAI_max = 2*LAI, no doubling).
#   Rscript run_v323iter_ch3_dryrun.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

# ---- v3.2.3 iter setup -------------------------------------------------------
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)  # v3.2.3 official (command, not staged)
# forcing must be a path RELATIVE to project root: setup_musica stages it as
# `ln -s ../<forcing> tempdir/forcing.nc`, so an absolute path makes a dangling
# symlink -> STOP 6. (Same convention as the v3.2.0 production runs.)
FORC <- "in_files/musica_in_Blois_pblh.nc"
stopifnot(file.exists(FORC))
ABL  <- list("abl_flag" = '"iter"')
LY   <- 2020:2022   # 2020 included so the day-before-sim seed (2020/366) exists
NCDIR <- file.path(CFG_C3$out_dir, "nc_v323iter_dryrun"); dir.create(NCDIR, recursive=TRUE, showWarnings=FALSE)

# ---- seed-row injection (the only difference vs the v3.2.0 phenology) ---------
inject_seed <- function(ph) {
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]
  if (nrow(d) != 1L) stop("could not find 2020/365 row to seed 2020/366")
  d$Julian_day <- 366L
  rbind(ph, d)
}
mk_dyn_iter <- function(ts_by_plot) function(pr) {
  pid <- plot_id_from_row(pr); if (!pid %in% names(ts_by_plot)) return(NULL)
  inject_seed(make_phenology_from_s2(ts_by_plot[[pid]], LY))
}
mk_static_iter <- function(lai_col) function(pr) {
  L <- as.numeric(pr[[lai_col]])   # one-sided; run_musica_one doubles Leaf_area
  inject_seed(calc_phenology(list.year = LY, nleafage = 1, budburst_date = 115,
              leaf_age_max_in = 0.56, relative_age_firstmax = 0.10,
              relative_age_lastmax = 0.75, LAI_max_per_cohort = L))
}

# ---- data + one ATBD per-plot series (same construction as production) --------
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
df[, pid := sprintf("X%d_Y%d", round(x), round(y))]
df[is.na(LAI_ALS_DOPT), LAI_ALS_DOPT := LAI_ALS]
nmdir <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(nmdir, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
stk <- rast(files); dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[, .(x,y)]), geom=c("x","y"), crs=crs(stk))
vals <- as.data.frame(terra::extract(stk, pts))[, -1, drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
cat(sprintf("ATBD series rebuilt: %d/53 plots\n", length(atbd)))

# pick the first plot that HAS a series
i_dyn <- which(df$pid %in% names(atbd))[1]
prow  <- as.data.frame(df[i_dyn])
cat(sprintf("dry-run plot: id=%s pid=%s  LAI_ALS=%.2f LAI_S2_ATBD=%.2f Hmax=%.1f fCover=%.2f\n",
            prow$id_plot, prow$pid, prow$LAI_ALS, prow$LAI_S2_ATBD, prow$Hmax, prow$fCover))

# ---- sanity: phenology has the 2020/366 seed and 2021/2022 full years --------
ph_chk <- mk_dyn_iter(atbd)(prow)
cat(sprintf("pheno rows: %d | 2020/366 present: %s | yrs: %s\n",
            nrow(ph_chk), any(ph_chk$year==2020 & ph_chk$Julian_day==366),
            paste(sort(unique(ph_chk$year)), collapse=",")))

# ---- run the two scenarios, timed -------------------------------------------
sc_dyn <- list(name="DYN_S2_ATBD_NM",
  lai_fn=function(pr) as.numeric(pr$LAI_S2_ATBD), hmax_fn=function(pr) as.numeric(pr$Hmax),
  fcover_fn=function(pr) as.numeric(pr$fCover), lad_fn=make_lad_real, phenology_fn=mk_dyn_iter(atbd))
sc_static <- list(name="STATIC_ALS",
  lai_fn=function(pr) as.numeric(pr$LAI_ALS), hmax_fn=function(pr) as.numeric(pr$Hmax),
  fcover_fn=function(pr) as.numeric(pr$fCover), lad_fn=make_lad_real, phenology_fn=mk_static_iter("LAI_ALS"))

for (sc in list(sc_dyn, sc_static)) {
  nc <- file.path(NCDIR, sprintf("musica_out_HOBO_%s__%s.nc", prow$id_plot, sc$name))
  if (file.exists(nc)) file.remove(nc)
  t0 <- Sys.time()
  run_musica_one(prow, sc, nc, FORC, BIN, extra_setup = ABL)
  dt <- as.numeric(difftime(Sys.time(), t0, units="secs"))
  ok <- file.exists(nc) && file.size(nc) > 1000
  cat(sprintf("[%s] nc written: %s  size: %s  time: %.1f s\n",
              sc$name, ok, if (file.exists(nc)) format(structure(file.size(nc), class="object_size"), units="auto") else "NA", dt))
  if (ok) {
    ncc <- nc_open(nc); has <- all(c("Tair_z","relative_height") %in% names(ncc$var)); nc_close(ncc)
    cat(sprintf("       readable, has Tair_z/relative_height: %s\n", has))
  }
}
cat("\nDRY RUN DONE\n")
