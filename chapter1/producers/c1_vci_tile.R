# ==============================================================================
# VCI worker: ONE tile, fresh R process (driver = c1_vci_raster_recompute_site.R).
# A fresh process per tile is REQUIRED: in-loop the page cache saturates with LiDAR
# data and a 2 GB tile's normalize_height(tin()) transient OOM-kills on a 30 GB box;
# exiting after each tile lets the OS reclaim all memory (proven: 2.0 GB tile, 52 s).
#   Rscript c1_vci_tile.R <infile.las> <outfile.tif>
# Reads 2-las_utm (absolute Z) -> normalize_height(tin()) -> VCI(Z,zmax) on ALL
# returns above ground (van Ewijk 2011, 1 m bins) at 10 m. Validated recipe.
# ==============================================================================
suppressPackageStartupMessages({ library(lidR); library(terra) })
set_lidr_threads(1); options(lidR.progress = FALSE)

a <- commandArgs(trailingOnly = TRUE)
infile <- a[1]; outfile <- a[2]

las <- tryCatch(readLAS(infile, select = "xyzc"), error = function(e) NULL)
if (is.null(las) || is.empty(las)) quit(save = "no", status = 0)
las <- tryCatch(normalize_height(las, tin()), error = function(e) NULL)
if (is.null(las)) quit(save = "no", status = 2)

vci_fun <- function(z) {
  z <- z[is.finite(z) & z >= 0]
  if (length(z) < 50) return(NA_real_)
  zmax <- max(z); if (zmax < 2) return(NA_real_)
  lidR::VCI(z, zmax = zmax)
}
r <- pixel_metrics(las, ~vci_fun(Z), res = 10)
terra::writeRaster(r, outfile, overwrite = TRUE)
quit(save = "no", status = 0)
