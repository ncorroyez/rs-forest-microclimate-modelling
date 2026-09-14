# ==============================================================================
# Stage 2: mosaic the per-tile native-20 m rasters (c1_native_raster_tile.R output)
# into the 5 pipeline rasters, written to in_files_native20/ with the SAME filenames
# the pipeline expects, so load_lidar_rasters(in_files_native20, agg_factor=1) uses
# the native traits directly (aggregation becomes a no-op). Nothing is overwritten.
# ==============================================================================
suppressPackageStartupMessages({ library(terra) })
TDIR <- "out_files/native20/tiles"; OUT <- "in_files_native20"; dir.create(OUT, showWarnings = FALSE)
tifs <- list.files(TDIR, "\\.tif$", full.names = TRUE)
cat(sprintf("mosaicking %d tiles...\n", length(tifs)))
rl <- lapply(tifs, rast)
lyr <- names(rl[[1]])                                   # hmax,hp95,fcov,vci,sec,lai,LAD_Layer_1.5..39.5
moslyr <- function(nm) {                                # mosaic one layer across all tiles (mean on overlap)
  sc <- sprc(lapply(rl, function(r) r[[nm]]))
  terra::mosaic(sc, fun = "mean")
}
m <- rast(lapply(lyr, moslyr)); names(m) <- lyr
cat(sprintf("mosaic: %d x %d, res %g, %d layers\n", nrow(m), ncol(m), res(m)[1], nlyr(m)))
# write the 5 pipeline rasters (keep the *_res_10_m.tif names load_lidar_rasters expects)
writeRaster(m[["lai"]],  file.path(OUT, "lai_z1_res_10_m.tif"), overwrite = TRUE)
writeRaster(m[["vci"]],  file.path(OUT, "vci_res_10_m.tif"),   overwrite = TRUE)
writeRaster(m[["hmax"]], file.path(OUT, "max_res_10_m.tif"),   overwrite = TRUE)
writeRaster(m[["fcov"]], file.path(OUT, "fCover_res_10_m.tif"),overwrite = TRUE)
lad <- m[[grep("LAD_Layer", lyr, value = TRUE)]]
writeRaster(lad, file.path(OUT, "lad_profiles_z1_res_10_m.tif"), overwrite = TRUE)
# scan-angle secant raster (kept for the per-plot global correction step)
writeRaster(m[["sec"]], file.path(OUT, "sec_theta_res_10_m.tif"), overwrite = TRUE)
cat("wrote 5 pipeline rasters + sec_theta to", OUT, "\n")
cat(sprintf("LAI[%.1f,%.1f] Hmax[%.1f,%.1f] fCover[%.2f,%.2f] VCI[%.2f,%.2f]\n",
  global(m[["lai"]],"min",na.rm=T)[[1]], global(m[["lai"]],"max",na.rm=T)[[1]],
  global(m[["hmax"]],"min",na.rm=T)[[1]], global(m[["hmax"]],"max",na.rm=T)[[1]],
  global(m[["fcov"]],"min",na.rm=T)[[1]], global(m[["fcov"]],"max",na.rm=T)[[1]],
  global(m[["vci"]],"min",na.rm=T)[[1]], global(m[["vci"]],"max",na.rm=T)[[1]]))
