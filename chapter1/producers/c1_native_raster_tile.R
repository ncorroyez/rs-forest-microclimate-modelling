# ==============================================================================
# Native trait-raster worker: ONE tile, fresh R process (driver = c1_native_raster_site.R).
# Fresh process per tile (page-cache OOM guard, same rationale as c1_vci_tile.R).
# Computes at a given native resolution, from the leaf-on ALS point cloud, the SAME
# traits the pipeline uses: Hmax (max + p95), fCover, VCI, <sec theta>, and the
# MacArthur-Horn LAD profile (k=0.5, dz=1, z0=1, 39 bins 1.5..39.5 m) with LAI as its
# vertical integral. Method = foliage-profile inversion documented in §2.2 (lidR::LAD).
#   Rscript c1_native_raster_tile.R <infile.las> <outfile.tif> <res_m>
# ==============================================================================
suppressPackageStartupMessages({ library(lidR); library(terra) })
set_lidr_threads(1); options(lidR.progress = FALSE, lidR.verbose = FALSE)
a <- commandArgs(trailingOnly = TRUE); infile <- a[1]; outfile <- a[2]; RES <- as.numeric(a[3])
DZ <- 1; K <- 0.5; Z0 <- 1; ZB <- seq(1.5, 39.5, by = 1)          # fixed 39-bin grid

las <- tryCatch(readLAS(infile, select = "xyzca"), error = function(e) NULL)  # c = class (ground), a = scan angle
if (is.null(las) || is.empty(las)) quit(save = "no", status = 0)
las <- tryCatch(normalize_height(las, tin()), error = function(e) NULL)
if (is.null(las)) quit(save = "no", status = 2)
ang_name <- intersect(c("ScanAngle","ScanAngleRank"), names(las@data))
ang_name <- if (length(ang_name)) ang_name[1] else NA

trait_fun <- function(z, ang) {
  ok <- is.finite(z) & z >= 0; z <- z[ok]
  out <- c(hmax = NA_real_, hp95 = NA_real_, fcov = NA_real_, vci = NA_real_,
           sec = NA_real_, lai = NA_real_, setNames(rep(NA_real_, length(ZB)), sprintf("LAD_Layer_%.1f", ZB)))
  if (length(z) < 50) return(as.list(out))
  zmax <- max(z)
  out["hmax"] <- zmax; out["hp95"] <- as.numeric(stats::quantile(z, 0.95)); out["fcov"] <- mean(z > 2)
  if (zmax >= 2) out["vci"] <- lidR::VCI(z, zmax = zmax)
  if (!is.na(ang_name)) { av <- ang[ok]; av <- av[is.finite(av)]; if (length(av)) out["sec"] <- mean(1/cos(av*pi/180)) }
  d <- tryCatch(lidR::LAD(z, dz = DZ, k = K, z0 = Z0), error = function(e) NULL)
  if (!is.null(d) && nrow(d)) {
    idx <- match(round(d$z, 1), ZB); keep <- !is.na(idx)
    if (any(keep)) out[6 + idx[keep]] <- d$lad[keep]
    out["lai"] <- sum(d$lad, na.rm = TRUE) * DZ
  }
  as.list(out)
}
form <- if (!is.na(ang_name)) as.formula(sprintf("~trait_fun(Z, %s)", ang_name)) else ~trait_fun(Z, Z)
r <- pixel_metrics(las, form, res = RES)
terra::writeRaster(r, outfile, overwrite = TRUE)
quit(save = "no", status = 0)
