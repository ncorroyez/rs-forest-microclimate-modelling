# ==============================================================================
# Is the simulated understory ΔTmax sensitive to the prescribed surface-boundary-layer
# height (h_sbl) that the iterative ABL ("yoyo") coupling takes as input?
#
# WHY THIS EXISTS. Appendix H and Section 2.4 both lean on the answer: the canonical
# forcing FR-Blo_2021_v2.nc carries an h_sbl series whose provenance is external to this
# work (hourly, median 48.8 m; NOT the MERRA-2 series of build_forcing_pblh.R, r = 0.54
# against it). The chapter's defence is that the exact level does not matter. That test
# was run once in 2026-07-23 and its outputs were deleted, leaving a load-bearing claim
# with no producer. This script is the producer.
#
# METHOD. For a sample of real cLHS plots, run the reference (full-canopy) configuration
# twice, identical in every respect except a CONSTANT h_sbl of 100 m and of 5000 m, two
# orders of magnitude apart, spanning the median of the archived series. Compare the summer-mean
# ΔTmax at 1 m. Everything else follows the native20 recipe: v3.2.3, abl_flag='iter',
# wind correction at the plot's base height, real LiDAR LAD, phenology carrying the
# day-before-simulation that iter requires.
#
# THIS SCRIPT RUNS MuSICA. It writes only under out_files/Chapter1/nc_blh_insensitivity
# and its own forcing cache; it touches no canonical output. 2 runs per plot.
#   Rscript scripts/c1_blh_insensitivity.R [n_plots]        (default 8)
# Reads :
#         in_files/FR-Blo_2021_v2.nc / in_files/model-3.2.3/musica
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
# Out: out_files/Chapter1/tables/tab_blh_insensitivity.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4); library(lubridate); library(data.table)
  library(parallel); library(rmusica); library(musica.tools)})
args <- commandArgs(trailingOnly = TRUE)
NP <- if (length(args)) as.integer(args[1]) else 8L
src <- list.files("R", "\\.R$", full.names = TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("Chapter3_config.R"); source("R/wind_correction.R"); source("R/cluster_relabel.R")
CFG_C3$musica_cmd <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC <- "in_files/FR-Blo_2021_v2.nc"; ABL <- list("abl_flag" = '"iter"')
# NOTE on the yoyo gradient diagnostic, tried 2026-07-31 and NOT retained. gradm_hsbl and
# gradc_hsbl would discriminate between "ΔTmax is insensitive to the boundary-layer height"
# and "the coupling is inert so h_sbl has no pathway", but they are not in the default
# history_variables list and come out as pure _FillValue (9.97e36). Adding them through
# extra_setup made MuSICA reject the namelist (STOP 6, 32-byte stub outputs), so the two
# are reported below only when actually written. The question stays open on the model side;
# it does not affect the ΔTmax result, which is measured directly.
HSBL <- c(low = 100, high = 5000)                      # two orders of magnitude apart
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
MREF <- macro_ref(FORC, ds)
NCDIR <- "out_files/Chapter1/nc_blh_insensitivity"; dir.create(NCDIR, recursive = TRUE, showWarnings = FALSE)

# Plot sample: real cLHS plots, spread across the four archetypes so the answer is not
# read off one canopy. Same seed discipline as c1_noise_floor.R.
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
samp[, P := relabel_cluster(Cluster)]
set.seed(7); pick <- samp[, .SD[sample(.N, max(1L, ceiling(NP / 4)))], by = P]
pick[, pid := sprintf("B%03d", .I)]

# A forcing identical to the wind-corrected one, with h_sbl replaced by a constant.
# Project-relative, like windcorr_forcing: MuSICA stages the forcing via `ln -s ../<forcing>`.
#' A copy of the wind-corrected forcing with h_sbl overwritten by a constant
#' @param hmax canopy height, selecting which wind-corrected forcing to start from
#' @param hsbl the constant surface-boundary-layer height to write (m)
#' @param cache_dir where the derived forcings are cached
#' @return path to the forcing NetCDF; stops if h_sbl did not write as a constant
blh_forcing <- function(hmax, hsbl, cache_dir = "out_files/blh_forc") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  base <- windcorr_forcing(hmax)
  ff <- file.path(cache_dir, sprintf("%s_h%04d.nc", sub("\\.nc$", "", basename(base)), hsbl))
  if (!file.exists(ff) || file.size(ff) < 1e5) {
    file.copy(base, ff, overwrite = TRUE)
    nc <- nc_open(ff, write = TRUE)
    v <- nc$var[["h_sbl"]]                      # write on the variable's declared shape
    ncvar_put(nc, "h_sbl", array(hsbl, dim = sapply(v$dim, function(d) d$len)))
    nc_close(nc)
    chk <- nc_open(ff); got <- as.numeric(ncvar_get(chk, "h_sbl")); nc_close(chk)
    if (!all(got == hsbl))                      # assert, do not assume, before spending runs
      stop("blh_forcing: h_sbl not written as a constant in ", ff,
           " (min ", min(got), ", max ", max(got), ")")
  }
  ff
}
# Build every forcing SERIALLY before the fan-out. windcorr_forcing caches by the wind
# factor rounded to 2 dp, so two plots can target the same file; letting mclapply workers
# copy and nc_open(write=TRUE) it concurrently would silently corrupt it and still produce
# output. c1_noise_floor.R has the same latent race.
invisible(lapply(seq_len(nrow(pick)), function(i)
  lapply(HSBL, function(h) blh_forcing(pick$Hmax[i], h))))
#' Phenology builder carrying the extra day-before-simulation that abl_flag='iter' needs
#' @param l one-sided LAI max per cohort
#' @return a FUNCTION(p) returning the phenology data.frame with 2020/366 appended
mk_phen <- function(l) function(p) {                    # iter needs the day before the sim
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
    leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
    LAI_max_per_cohort = l))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366; rbind(ph, d)
}
# The yoyo's own gradient terms, straight from the output. A near-zero ΔTmax difference is
# consistent with two different worlds: ΔTmax is genuinely insensitive to boundary-layer
# height, or the coupling is inert here so h_sbl has no pathway to act at all. These two
# variables discriminate, and they cost nothing.
#' Mean absolute yoyo gradient terms from one run, _FillValue dropped
#' @param nc path to a MuSICA output NetCDF
#' @return named numeric c(gm, gc); NA where the variable was not written
grads <- function(nc) {
  if (!file.exists(nc)) return(c(gm = NA_real_, gc = NA_real_))
  h <- nc_open(nc); on.exit(nc_close(h))
  gv <- function(v) {
    if (!v %in% names(h$var)) return(NA_real_)
    z <- as.numeric(ncvar_get(h, v))
    z <- z[is.finite(z) & abs(z) < 1e30]        # drop _FillValue (9.97e36): "not written"
    if (!length(z)) NA_real_ else mean(abs(z))
  }
  c(gm = gv("gradm_hsbl"), gc = gv("gradc_hsbl"))
}
#' Run (or re-read) both h_sbl levels for one plot and score them
#' @param i row index in the sampled plot table `pick`
#' @return one-row data.table: traits, dTmax at 100 m and 5000 m, and both gradient terms
one <- function(i) {
  p <- as.data.frame(pick[i])
  g <- function(tag, hsbl) {
    nc <- file.path(NCDIR, sprintf("%s_%s.nc", p$pid, tag))
    sc <- list(lai_fn = function(x) p$LAI, hmax_fn = function(x) p$Hmax,
               fcover_fn = function(x) p$fCover, lad_fn = make_lad_real,
               phenology_fn = mk_phen(p$LAI))
    if (!file.exists(nc) || file.size(nc) < 1000)
      run_musica_one(p, sc, nc, blh_forcing(p$Hmax, hsbl), CFG_C3$musica_cmd, extra_setup = ABL)
    m <- micro_hourly_at(nc, 1.0)
    c(dt = if (is.null(m)) NA_real_ else delta_tmax_mean(m, MREF, ds), grads(nc))
  }
  a <- g("h0100", HSBL[["low"]]); b <- g("h5000", HSBL[["high"]])
  data.table(pid = p$pid, P = p$P, LAI = p$LAI, Hmax = p$Hmax, fCover = p$fCover,
             dt_100 = a[["dt"]], dt_5000 = b[["dt"]],
             gradm_100 = a[["gm"]], gradm_5000 = b[["gm"]],
             gradc_100 = a[["gc"]], gradc_5000 = b[["gc"]])
}
R <- rbindlist(mclapply(seq_len(nrow(pick)),
       function(i) tryCatch(one(i), error = function(e) {cat("ERR", i, conditionMessage(e), "\n"); NULL}),
       mc.cores = 4), fill = TRUE)
R <- R[is.finite(dt_100) & is.finite(dt_5000)][, diff := dt_5000 - dt_100]
dir.create("out_files/Chapter1/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(R, "out_files/Chapter1/tables/tab_blh_insensitivity.csv")
cat(sprintf("\n=== BLH INSENSITIVITY (h_sbl = %g m vs %g m, n = %d plots) ===\n",
            HSBL[["low"]], HSBL[["high"]], nrow(R)))
print(R[, .(pid, P, dt_100 = round(dt_100, 5), dt_5000 = round(dt_5000, 5), diff = round(diff, 6))])
cat(sprintf("\n  |dTmax(5000 m) - dTmax(100 m)| : median %.6f | max %.6f degC\n",
            median(abs(R$diff)), max(abs(R$diff))))
cat(sprintf("  for scale, the LAI +0.5 lever in P4 is about 0.26 degC\n"))
cat("\n=== yoyo gradient terms (is the coupling active at all?) ===\n")
print(R[, .(pid, gradm_100 = signif(gradm_100, 3), gradm_5000 = signif(gradm_5000, 3),
            gradc_100 = signif(gradc_100, 3), gradc_5000 = signif(gradc_5000, 3))])
cat(sprintf("  mean|gradm|: %.4g (100 m) vs %.4g (5000 m) | mean|gradc|: %.4g vs %.4g\n",
            mean(R$gradm_100, na.rm = TRUE), mean(R$gradm_5000, na.rm = TRUE),
            mean(R$gradc_100, na.rm = TRUE), mean(R$gradc_5000, na.rm = TRUE)))
cat("  If these are ~0 or identical across the two h_sbl levels, the insensitivity is\n",
    "  STRUCTURAL (h_sbl has no pathway) rather than a physical near-cancellation.\n", sep = "")
cat("DONE\n")
