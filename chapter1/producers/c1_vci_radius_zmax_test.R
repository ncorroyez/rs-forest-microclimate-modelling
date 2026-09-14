# ==============================================================================
# Why do our VCI values have no points near 0, unlike Fig. 4 of Gril et al. 2023?
# JO, 2026-08-20. Two candidate causes, tested together on the 53 loggers:
#   RADIUS  Gril uses a 5 m buffer (and reports the predictive power of VCI
#           dropping beyond 10 m); the chapter's VCI is a 10 m-pixel product.
#   ZMAX    lidR::VCI(z, zmax) normalizes the Shannon entropy by log(zmax/by).
#           We pass the PLOT's own max height, so a 3 m open stand has ~3 bins
#           and lands mid-range. A fixed site-wide ceiling gives it ~40 bins with
#           returns in 3, and the entropy collapses toward 0.
# 40.5 m is Gril's PAD slice ceiling for Blois, used here as the fixed ceiling.
# This does NOT assert what she passed to lidR::VCI; the paper does not say.
#
# No `zmax < 2` guard here, unlike c1_vci_tile.R: that guard drops exactly the
# plots most likely to sit near 0, so it would hide the effect being tested.
#
# Reads : /media/corroyez/MyPassport/.../2-las_utm  (leaf-on catalog, all returns)
#         CFG$hobo_geojson
# Writes: out_files/Chapter1/tables/tab_vci_radius_zmax.csv
#   Rscript scripts/c1_vci_radius_zmax_test.R
# ==============================================================================
suppressPackageStartupMessages({library(lidR); library(data.table); library(sf); library(parallel)})
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source))
CTG <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
ZFIX <- 40.5                              # Gril's PAD ceiling for the Blois forest
stopifnot(dir.exists(CTG))
ctg <- readLAScatalog(CTG, select = "xyzc"); opt_progress(ctg) <- FALSE
set_lidr_threads(1)

hp <- st_read(CFG$hobo_geojson, quiet = TRUE); hp <- hp[!hp$id_plot %in% CFG$ids_to_remove, ]
co <- st_coordinates(st_geometry(hp))
cat(sprintf("%d loggers\n", nrow(hp)))

vci_pair <- function(pt) {
  if (is.null(pt) || npoints(pt) < 50) return(c(own = NA, fixed = NA, zmax = NA, n = 0))
  ptn <- tryCatch(normalize_height(pt, tin()), error = function(e) NULL)
  if (is.null(ptn)) return(c(own = NA, fixed = NA, zmax = NA, n = 0))
  z <- ptn$Z; z <- z[is.finite(z) & z >= 0]
  if (!length(z)) return(c(own = NA, fixed = NA, zmax = NA, n = 0))
  zm <- max(z)
  c(own   = tryCatch(VCI(z, zmax = zm),   error = function(e) NA_real_),
    fixed = tryCatch(VCI(z, zmax = ZFIX), error = function(e) NA_real_),
    zmax = zm, n = length(z))
}
out <- rbindlist(lapply(c(5, 10, 25), function(rad) {
  cl <- clip_circle(ctg, co[, 1], co[, 2], rad)
  if (inherits(cl, "LAS")) cl <- list(cl)
  v <- do.call(rbind, mclapply(cl, vci_pair, mc.cores = 4))
  cat(sprintf("  radius %2d m done\n", rad))
  data.table(id_plot = hp$id_plot, radius = rad, vci_own = v[, "own"],
             vci_fixed = v[, "fixed"], zmax = v[, "zmax"], npts = v[, "n"]) }))
fwrite(out, "out_files/Chapter1/tables/tab_vci_radius_zmax.csv")

cat("\n=== VCI at the 53 loggers: own-max vs fixed 40.5 m ceiling ===\n")
print(out[, .(n = sum(is.finite(vci_own)),
              own_min = round(min(vci_own, na.rm = TRUE), 3),
              own_med = round(median(vci_own, na.rm = TRUE), 3),
              own_max = round(max(vci_own, na.rm = TRUE), 3),
              fix_min = round(min(vci_fixed, na.rm = TRUE), 3),
              fix_med = round(median(vci_fixed, na.rm = TRUE), 3),
              fix_max = round(max(vci_fixed, na.rm = TRUE), 3),
              own_below_0.2 = sum(vci_own < 0.2, na.rm = TRUE),
              fix_below_0.2 = sum(vci_fixed < 0.2, na.rm = TRUE)), by = radius][order(radius)])
cat("\nchapter's VCI (native 20 m raster, own-max): ")
C <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
cat(sprintf("%.3f - %.3f (median %.3f), %d below 0.2\n",
            min(C$VCI), max(C$VCI), median(C$VCI), sum(C$VCI < 0.2)))
cat("DONE\n")
