# ==============================================================================
# c1_radius_clhs_vars.R
#
# Structural-variable range of the cLHS radius sweep. For each cLHS plot and each
# clipping radius, re-clip the ALS cloud and recompute the four canopy variables
# that feed MuSICA, with EXACTLY the runner's parameters (c1_radius_clhs.R):
# scan-angle-corrected one-sided LAI, Hmax capped at 40 m, fCover floored at 0.5,
# and the LAD profile (summarised here by the relative height of its leaf-area
# centroid). No MuSICA, no scoring: this answers only "how far does the clipping
# radius move the extracted structure", read archetype by archetype.
#
# The 20 m radius is the reference reported throughout, so each variable is also
# expressed as a paired difference against its own 20 m value.
#
# Inputs : out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
#          <LiDAR catalogue on the external drive>
# Output : out_files/Chapter1/tables/radius_clhs_vars.csv
#
#   Rscript c1_radius_clhs_vars.R                 # full sweep, all radii
#   Rscript c1_radius_clhs_vars.R 5,20 20         # DRY RUN: radii 5 & 20, 20 plots
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(data.table); library(parallel) })
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
source("pipeline/00_config.R")
options(lidR.progress = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
RADII <- if (length(args) >= 1) as.numeric(strsplit(args[1], ",")[[1]]) else
         c(5, 10, 12.5, 15, 20, 25, 50)
NMAX  <- if (length(args) >= 2) as.integer(args[2]) else Inf
REF   <- 20

CTG <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
DZ  <- 0.5; Z0 <- 0.5; KEXT <- 0.5                 # identical to c1_radius_clhs.R
ctg <- readLAScatalog(CTG); opt_select(ctg) <- "*"; opt_progress(ctg) <- FALSE

CL   <- readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds")
xy   <- as.matrix(CL[, c("x", "y")])
ids  <- sprintf("clhs_%03d", seq_len(nrow(CL)))
ARCH <- as.character(relabel_cluster(CL$Cluster))
np   <- min(nrow(xy), NMAX)
cat(sprintf("radii: %s | plots: %d | cores: 10\n",
            paste0(RADII, "m", collapse = ","), np))

# One plot at one radius: the four canopy variables MuSICA is given, plus the
# relative centroid height of the LAD profile as a scalar summary of shape.
one <- function(i, radius) {
  las <- tryCatch(suppressMessages(clip_circle(ctg, xy[i, 1], xy[i, 2], radius)),
                  error = function(e) NULL)
  if (is.null(las) || is.empty(las) || nrow(las@data) < 50) return(NULL)
  las <- tryCatch(filter_duplicates(las), error = function(e) las)
  las <- tryCatch(suppressMessages(normalize_height(las, knnidw())), error = function(e) NULL)
  if (is.null(las)) return(NULL)
  fin <- is.finite(las@data$Z) & las@data$Z >= 0; Z <- las@data$Z[fin]
  if (length(Z) < 50) return(NULL)
  sa  <- if ("ScanAngle" %in% names(las@data)) las@data$ScanAngle[fin] else NULL
  sec <- if (!is.null(sa) && any(is.finite(sa))) mean(1 / cos(sa[is.finite(sa)] * pi / 180)) else 1
  d   <- tryCatch(lidR::LAD(Z, dz = DZ, k = KEXT, z0 = Z0), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  lai  <- sum(d$lad, na.rm = TRUE) * DZ / sec
  hmax <- min(max(Z, na.rm = TRUE), 40)
  fcov <- max(mean(Z > 2, na.rm = TRUE), 0.5)
  # relative height of the leaf-area centroid (0 = ground, 1 = top): a scalar the
  # profile-shape lever can be read against, independent of Hmax.
  cen  <- sum(d$z * d$lad, na.rm = TRUE) / sum(d$lad, na.rm = TRUE)
  data.table(id_plot = ids[i], radius = radius, P = ARCH[i],
             LAI = lai, Hmax = hmax, fCover = fcov,
             lad_centroid_rel = cen / hmax, n_points = length(Z))
}

t0 <- Sys.time()
V <- rbindlist(lapply(RADII, function(r) {
  cat(sprintf("  radius %5s m ...\n", r))
  rbindlist(mclapply(seq_len(np), one, radius = r, mc.cores = 10, mc.preschedule = FALSE))
}), fill = TRUE)
dt_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("extraction: %.1f min for %d plot-radii (%.2f s/plot-radius)\n",
            dt_min, nrow(V), dt_min * 60 / max(nrow(V), 1)))

out <- "out_files/Chapter1/tables/radius_clhs_vars.csv"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
fwrite(V, out); cat("wrote ", out, "\n")

# ---- range per radius, and paired difference against the 20 m reference -------
rng <- V[, .(LAI_med = median(LAI, na.rm = TRUE),
             LAI_lo = quantile(LAI, .1, na.rm = TRUE), LAI_hi = quantile(LAI, .9, na.rm = TRUE),
             fCov_med = median(fCover, na.rm = TRUE), Hmax_med = median(Hmax, na.rm = TRUE)),
         by = .(radius, P)][order(P, radius)]
cat("\n=== variable range per radius, by archetype (median [p10,p90] LAI) ===\n")
print(rng)

if (REF %in% RADII) {
  R0 <- V[radius == REF, .(id_plot, LAI0 = LAI, Hmax0 = Hmax, fCov0 = fCover)]
  D  <- merge(V[radius != REF], R0, by = "id_plot")
  D[, `:=`(dLAI = LAI - LAI0, dHmax = Hmax - Hmax0, dfCov = fCover - fCov0)]
  pd <- D[, .(n = .N,
              dLAI = round(median(dLAI, na.rm = TRUE), 3),
              dHmax = round(median(dHmax, na.rm = TRUE), 2),
              dfCov = round(median(dfCov, na.rm = TRUE), 3)),
          by = .(radius, P)][order(P, radius)]
  cat("\n=== paired median shift vs 20 m reference ===\n"); print(pd)
}
cat("DONE\n")
