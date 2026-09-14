# ==============================================================================
# c1_radius_clhs_perturb.R
#
# Footprint robustness of the ATTRIBUTION (Fig 4/5 levers), radius by radius on
# the 400 cLHS plots. For one clipping radius, re-clip each plot, recompute the
# canopy variables and the real LAD profile at that radius, and run the native-unit
# perturbation design of c1_sensitivity_perplot_units.R around that radius-r
# baseline: 8 MuSICA sims per plot (base, uniform-LAD, LAI +/-, Hmax +/-, fCover +/-),
# levers reported in degrees per unit with the + and - side kept separate, plus the
# real-versus-uniform profile contrast dT_LAD = base - unif.
#
# NOT the published Fig 4/5. Three things differ from the shipped 3200 runs and the
# new levers will not reproduce the reported P4 numbers (0.26/0.27 degC) -- this is
# correct, not a bug:
#   (1) forcing  : CHS41 station-hsbl-Rmerge (aligned with the radius sweep), NOT FR-Blo_2021_v2;
#   (2) wind     : NO wind correction (the radius sweep convention), whereas Fig 4/5 used
#                  windcorr_forcing(Hmax), an Hmax-coupled forcing;
#   (3) baseline : variables re-clipped in a circle of this radius, NOT the stored
#                  native-20 m-grid values.
#
# dTmax convention = the delivered radius sweep (c1_radius_clhs_extract.R):
# macro_ref + delta_tmax_mean, time-matched to the macro daily max, NO clock shift.
# Plot key = clhs_%03d by row index, so the table joins radius_clhs_perplot.csv.
# nc cache is per-radius (out_files/radius_test_clhs_perturb/<rp>/), so no radius
# ever reads another radius's cached sim.
#
# Inputs : out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
#          <LiDAR catalogue>, out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc
# Output : out_files/Chapter1/tables/sensitivity_radius_clhs/r<rp>.csv
#
#   Rscript c1_radius_clhs_perturb.R 5           # one radius, all plots (~2.5 h)
#   Rscript c1_radius_clhs_perturb.R 5 2         # DRY RUN: radius 5 m, 2 plots
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(ncdf4); library(lubridate)
  library(data.table); library(parallel); library(rmusica); library(musica.tools) })
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("pipeline/00_config.R")
options(lidR.progress = FALSE)

args   <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript c1_radius_clhs_perturb.R <radius_m> [nmax]")
RADIUS <- as.numeric(args[1])
CHUNK  <- if (length(args) >= 2) as.integer(args[2]) else 1L
K      <- if (length(args) >= 3) as.integer(args[3]) else 1L
NMAX   <- { e <- Sys.getenv("DRY_NMAX"); if (nzchar(e)) as.integer(e) else Inf }  # cap plots for dry runs
rp     <- gsub("\\.", "p", as.character(RADIUS))

CTG  <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
DZ   <- 0.5; Z0 <- 0.5; KEXT <- 0.5                     # identical to c1_radius_clhs.R
MB   <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
ABL  <- list("abl_flag" = '"iter"')
FORC <- "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"   # Rmerge, no wind corr
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
Z_FIX <- 1.0
MREF <- macro_ref(FORC, ds)                             # canonical macro (daily max + its hour)

# fixed native steps and physical bounds -- identical to c1_sensitivity_perplot_units.R
D_LAI <- 0.5; D_HMAX <- 1.0; D_FCOV <- 0.10
LAI_FLOOR <- 0.1; HMAX_FLOOR <- 2; HMAX_CEIL <- 40; FCOV_MIN <- 0.5; FCOV_MAX <- 1

mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
        leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
        LAI_max_per_cohort = lai1))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366; rbind(ph, d) }

CL  <- readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds")
xy  <- as.matrix(CL[, c("x", "y")])
ids <- sprintf("clhs_%03d", seq_len(nrow(CL)))
ARCH <- as.character(relabel_cluster(CL$Cluster))
np  <- min(nrow(xy), NMAX)
ctg <- readLAScatalog(CTG); opt_select(ctg) <- "*"; opt_progress(ctg) <- FALSE

# Parallelism is at the SHELL level (one process per chunk), NOT mclapply: forking
# after lidR/terra open the LiDAR catalogue deadlocks the workers (0 nc, 0 % CPU,
# observed). This mirrors the working c1_sensitivity_perplot_units.R chunk pattern.
# Each invocation processes its modulo-K chunk sequentially; the orchestrator runs
# K invocations in parallel per radius.

NCDIR <- file.path("out_files/radius_test_clhs_perturb", rp); dir.create(NCDIR, recursive = TRUE, showWarnings = FALSE)
metric_of <- function(nc) { m <- micro_hourly_at(nc, Z_FIX); if (is.null(m)) return(NA_real_); delta_tmax_mean(m, MREF, ds) }
per_unit  <- function(sim, base, step) if (!is.finite(sim) || !is.finite(base) || abs(step) < 1e-9) NA_real_ else (sim - base) / step

# Re-clip one plot at RADIUS and build a prow carrying the radius-r variables and
# LAD_Layer_* shape columns that make_lad_real / make_lad_uniform consume.
clip_prow <- function(i) {
  las <- tryCatch(suppressMessages(clip_circle(ctg, xy[i, 1], xy[i, 2], RADIUS)), error = function(e) NULL)
  if (is.null(las) || is.empty(las) || nrow(las@data) < 50) return(NULL)
  las <- tryCatch(filter_duplicates(las), error = function(e) las)
  las <- tryCatch(suppressMessages(normalize_height(las, knnidw())), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  fin <- is.finite(las@data$Z) & las@data$Z >= 0; Z <- las@data$Z[fin]; if (length(Z) < 50) return(NULL)
  sa  <- if ("ScanAngle" %in% names(las@data)) las@data$ScanAngle[fin] else NULL
  sec <- if (!is.null(sa) && any(is.finite(sa))) mean(1 / cos(sa[is.finite(sa)] * pi / 180)) else 1
  d   <- tryCatch(lidR::LAD(Z, dz = DZ, k = KEXT, z0 = Z0), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || sum(d$lad, na.rm = TRUE) <= 0) return(NULL)
  lai  <- sum(d$lad, na.rm = TRUE) * DZ / sec
  hmax <- min(max(Z, na.rm = TRUE), 40)
  fcov <- max(mean(Z > 2, na.rm = TRUE), 0.5)
  prow <- data.frame(Cluster = CL$Cluster[i], x = xy[i, 1], y = xy[i, 2],
                     id_plot = ids[i], P = ARCH[i], Hmax = hmax, LAI = lai, fCover = fcov,
                     stringsAsFactors = FALSE)
  # LAD_Layer_* carry the profile SHAPE; make_lad_real renormalises to the target LAI.
  lad <- as.data.frame(t(d$lad)); names(lad) <- sprintf("LAD_Layer_%s", d$z)
  cbind(prow, lad)
}

# One perturbation arm: set the scalar levers and the LAD constructor, run, score.
run_get <- function(prow, tag, lai, hmax, fcov, ladf = make_lad_real) {
  nc <- file.path(NCDIR, sprintf("%s_%s.nc", prow$id_plot, tag))
  sc <- list(lai_fn = function(p) lai, hmax_fn = function(p) hmax, fcover_fn = function(p) fcov,
             lad_fn = ladf, phenology_fn = mk_phen(lai))
  if (!file.exists(nc) || file.size(nc) < 1000)
    tryCatch(run_musica_one(prow, sc, nc, FORC, MB, extra_setup = ABL), error = function(e) NULL)
  metric_of(nc) }

one_plot <- function(i) {
  prow <- tryCatch(clip_prow(i), error = function(e) NULL)
  if (is.null(prow)) return(NULL)
  laip_v <- prow$LAI + D_LAI;                 laim_v <- max(prow$LAI - D_LAI, LAI_FLOOR)
  hmp_v  <- min(prow$Hmax + D_HMAX, HMAX_CEIL); hmm_v <- max(prow$Hmax - D_HMAX, HMAX_FLOOR)
  fcp_v  <- min(prow$fCover + D_FCOV, FCOV_MAX); fcm_v <- max(prow$fCover - D_FCOV, FCOV_MIN)
  base <- run_get(prow, "base",    prow$LAI, prow$Hmax, prow$fCover)
  unif <- run_get(prow, "unifLAD", prow$LAI, prow$Hmax, prow$fCover, ladf = make_lad_uniform)
  laip <- run_get(prow, "LAIp",  laip_v,   prow$Hmax, prow$fCover)
  laim <- run_get(prow, "LAIm",  laim_v,   prow$Hmax, prow$fCover)
  hmp  <- run_get(prow, "Hmaxp", prow$LAI, hmp_v,     prow$fCover)
  hmm  <- run_get(prow, "Hmaxm", prow$LAI, hmm_v,     prow$fCover)
  fcp  <- run_get(prow, "fCovp", prow$LAI, prow$Hmax, fcp_v)
  fcm  <- run_get(prow, "fCovm", prow$LAI, prow$Hmax, fcm_v)
  data.table(radius = RADIUS, id_plot = prow$id_plot, P = prow$P, Cluster = prow$Cluster,
             LAI = prow$LAI, Hmax = prow$Hmax, fCover = prow$fCover, base = base,
             LAI_up  = per_unit(laip, base, laip_v - prow$LAI),  LAI_dn  = per_unit(laim, base, laim_v - prow$LAI),
             Hmax_up = per_unit(hmp,  base, hmp_v  - prow$Hmax), Hmax_dn = per_unit(hmm,  base, hmm_v  - prow$Hmax),
             fCov_up = per_unit(fcp,  base, fcp_v  - prow$fCover), fCov_dn = per_unit(fcm, base, fcm_v - prow$fCover),
             step_LAI = D_LAI, step_Hmax = D_HMAX, step_fCov = D_FCOV,
             dT_LAD = if (is.finite(base) && is.finite(unif)) base - unif else NA_real_)
}

my <- which(((seq_len(np) - 1L) %% K) == (CHUNK - 1L))
cat(sprintf("RADIUS %s m | chunk %d/%d | %d plots | forcing=Rmerge/no-wind | 8 sims/plot\n",
            RADIUS, CHUNK, K, length(my)))
t0 <- Sys.time()
R  <- rbindlist(lapply(my, function(i) {
        r <- tryCatch(one_plot(i), error = function(e) { cat(sprintf(" ERR %s: %s\n", ids[i], conditionMessage(e))); NULL })
        if (!is.null(r)) cat(sprintf(" r%s c%d: %s done\n", RADIUS, CHUNK, ids[i]))
        r }), fill = TRUE)
dt_h <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
cat(sprintf("RADIUS %s chunk %d DONE: %d plots in %.2f h\n", RADIUS, CHUNK, nrow(R), dt_h))

odir <- "out_files/Chapter1/tables/sensitivity_radius_clhs"; dir.create(odir, recursive = TRUE, showWarnings = FALSE)
of   <- if (K == 1L) file.path(odir, sprintf("r%s.csv", rp)) else file.path(odir, sprintf("r%s_part%02dof%02d.csv", rp, CHUNK, K))
fwrite(R, of); cat("wrote", of, "\n")

cat("\n=== median lever per archetype (degC per unit; +/- separate) ===\n")
print(R[, .(n = .N, LAI_up = round(median(LAI_up, na.rm = TRUE), 3),
            fCov_up = round(median(fCov_up, na.rm = TRUE), 3),
            Hmax_up = round(median(Hmax_up, na.rm = TRUE), 4),
            dT_LAD  = round(median(dT_LAD, na.rm = TRUE), 3)), by = P][order(P)])
cat("DONE\n")
