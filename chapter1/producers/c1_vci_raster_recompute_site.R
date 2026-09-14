# ==============================================================================
# VCI RASTER RECOMPUTE on NORMALIZED heights — DRIVER (one fresh R worker per tile).
# Fixes the absolute-elevation bug (VCI on raw Z ASL is compressed by ~ln(H)/ln(H+elev);
# c1_vci_normalize_check.R: raw 0.58 vs normalized 0.86). Wall-to-wall version of the
# validated point-wise recipe: per tile, 2-las_utm -> normalize_height(tin()) ->
# VCI(Z,zmax) on ALL returns above ground (van Ewijk 2011, 1 m bins) at 10 m.
#
# Each tile is processed in a SEPARATE Rscript (c1_vci_tile.R) so memory is fully
# reclaimed between tiles — in-loop processing OOM-kills on big interior tiles when
# the page cache is full (30 GB box). Per-tile .tif written + SKIPPED if present =>
# resumable: a reap/OOM loses at most one tile; just re-run. terra::vrt() mosaics.
#   Rscript c1_vci_raster_recompute_site.R <Blois|Aigoual|Mormal>
# Out: in_files/vci_tiles_<site>/<tile>_vci.tif ; in_files/vci_res_10_m_norm_<site>.tif
#      (Blois also overwrites the canonical in_files/vci_res_10_m_norm.tif)
# ==============================================================================
suppressPackageStartupMessages({ library(terra) })

site <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(site %in% c("Blois", "Aigoual", "Mormal"))

CFG <- list(
  Blois   = list(ctg = "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm",
                 geojson = "in_files/data_Blois_utm31n.geojson",
                 gx = "coord_x_utm31n", gy = "coord_y_utm31n"),
  Aigoual = list(ctg = "/media/corroyez/MyPassport/01_DATA/Aigoual/LiDAR/leaf_on/2-las_utm",
                 geojson = NA),
  Mormal  = list(ctg = "/media/corroyez/MyPassport/01_DATA/Mormal/LiDAR/leaf_on/2-las_utm",
                 geojson = NA)
)[[site]]

TILEDIR <- file.path("in_files", paste0("vci_tiles_", site))
OUT     <- file.path("in_files", paste0("vci_res_10_m_norm_", site, ".tif"))
WORKER  <- "c1_vci_tile.R"
dir.create(TILEDIR, recursive = TRUE, showWarnings = FALSE)

tiles <- list.files(CFG$ctg, pattern = "\\.las$", full.names = TRUE)
cat(sprintf("[%s] %d tiles ; one fresh R worker per tile ; res 10 m\n", site, length(tiles)))
t0 <- Sys.time(); ndone <- 0L; nskip <- 0L; nfail <- 0L
for (i in seq_along(tiles)) {
  f    <- tiles[i]
  outf <- file.path(TILEDIR, paste0(tools::file_path_sans_ext(basename(f)), "_vci.tif"))
  if (file.exists(outf)) { nskip <- nskip + 1L; next }            # resume: skip done tiles

  st <- system2("Rscript", c(WORKER, shQuote(f), shQuote(outf)),
                stdout = FALSE, stderr = FALSE)
  if (st == 0 && file.exists(outf)) {
    ndone <- ndone + 1L
    cat(sprintf("[%d/%d] %s  ok  (done %d, skip %d, fail %d, %.1f min)\n",
                i, length(tiles), basename(f), ndone, nskip, nfail,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  } else {                                                        # empty tile or worker error -> note, continue
    nfail <- nfail + 1L
    cat(sprintf("[%d/%d] %s  FAIL/empty (status %d)\n", i, length(tiles), basename(f), st))
  }
}
cat(sprintf("[%s] computed %d, skipped %d, failed/empty %d\n", site, ndone, nskip, nfail))

# ---- mosaic -----------------------------------------------------------------
done_files <- list.files(TILEDIR, pattern = "_vci\\.tif$", full.names = TRUE)
cat(sprintf("[%s] mosaicking %d tile rasters -> %s\n", site, length(done_files), OUT))
v   <- terra::vrt(done_files, overwrite = TRUE)
mos <- terra::writeRaster(v, OUT, overwrite = TRUE)
if (site == "Blois")
  terra::writeRaster(mos, "in_files/vci_res_10_m_norm.tif", overwrite = TRUE)

# ---- sanity: must land in the normalized band ~0.7-0.9, NOT compressed ~0.5 --
gs <- terra::global(mos, c("min", "mean", "max"), na.rm = TRUE)
cat(sprintf("[%s] norm VCI raster: min %.3f  mean %.3f  max %.3f  (buggy abs-Z mean ~0.53)\n",
            site, gs[1, 1], gs[1, 2], gs[1, 3]))

# ---- Blois only: validate against the point-wise ground truth (expect ~0.857) -
if (!is.na(CFG$geojson) && file.exists(CFG$geojson)) {
  suppressPackageStartupMessages(library(sf))
  g   <- sf::st_read(CFG$geojson, quiet = TRUE)
  pts <- terra::vect(cbind(g[[CFG$gx]], g[[CFG$gy]]), crs = terra::crs(mos))
  ev  <- terra::extract(mos, pts)[, 2]
  cat(sprintf("[%s] extracted at %d plots: median %.3f  (point-wise recompute gave 0.857)\n",
              site, sum(is.finite(ev)), median(ev, na.rm = TRUE)))
}
cat(sprintf("[%s] DONE in %.1f min\n", site, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
