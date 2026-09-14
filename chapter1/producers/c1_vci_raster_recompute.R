# ==============================================================================
# VCI RASTER RECOMPUTE on NORMALIZED heights (fixes the absolute-elevation bug,
# confirmed by c1_vci_normalize_check.R: raw 0.58 vs normalized 0.86).
# Per 10 m pixel: normalize_height(tin()) within a buffered chunk, then
# lidR::VCI(Z, zmax=max(Z)) on height-above-ground (van Ewijk 2011).
# Writes a NEW file (keeps the buggy vci_res_10_m.tif for provenance).
#   Rscript c1_vci_raster_recompute.R
# Out: in_files/vci_res_10_m_norm.tif
# ==============================================================================
suppressPackageStartupMessages({ library(lidR); library(terra); library(future); library(sf) })
plan(multisession, workers = 2)        # 2 workers = bounded RAM on 2 GB tiles (4 OOM'd; memory not chunk-size)

CTG <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
ctg <- readLAScatalog(CTG, select = "xyzc")
# DEFAULT per-FILE chunking (each tile read ONCE = one pass; 150 m re-chunking caused redundant I/O)
opt_chunk_buffer(ctg) <- 20            # clean DTM at tile edges; catalog_map crops it
opt_progress(ctg)     <- TRUE

vci_fun <- function(z) {               # van Ewijk normalised entropy, 1 m bins, per-pixel zmax
  z <- z[is.finite(z) & z >= 0]
  if (length(z) < 50) return(NA_real_)
  zmax <- max(z); if (zmax < 2) return(NA_real_)
  lidR::VCI(z, zmax = zmax)
}
chunk_fun <- function(las) {
  las <- normalize_height(las, tin())
  pixel_metrics(las, ~vci_fun(Z), res = 10)
}

cat("Recomputing VCI raster on normalized heights (160 tiles, per-tile chunks, 2 workers)...\n")
t0 <- Sys.time()
r  <- catalog_map(ctg, chunk_fun)
out <- "in_files/vci_res_10_m_norm.tif"
writeRaster(r, out, overwrite = TRUE)
cat(sprintf("DONE in %.1f min -> %s\n", as.numeric(difftime(Sys.time(), t0, units="mins")), out))

# ---- sanity: global stats + extract at plots vs the per-plot check ----------
gs <- global(r, c("min","mean","max"), na.rm = TRUE)
cat(sprintf("raster norm VCI: min %.3f  mean %.3f  max %.3f  (buggy raster mean was 0.53)\n",
            gs[1,1], gs[1,2], gs[1,3]))
g  <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
pts <- vect(cbind(g$coord_x_utm31n, g$coord_y_utm31n), crs = crs(r))
ev  <- terra::extract(r, pts)[, 2]
cat(sprintf("extracted at %d plots: median %.3f (per-plot check gave 0.857)\n",
            sum(is.finite(ev)), median(ev, na.rm = TRUE)))
