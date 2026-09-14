# ==============================================================================
# ROBUST corrected VCI at the 400 cLHS pixels — tile-by-tile, SINGLE-THREAD,
# RESUMABLE. Avoids the OOM that killed the parallel runs (one tile in memory at
# a time) and the slowness of random clip_circle (each tile read at most once).
# Per owning tile: read once, normalize_height(tin()), clip 25 m around each
# cLHS pixel it owns, VCI on ALL returns above ground (van Ewijk 2011).
# Writes results after EVERY tile → killing it loses nothing; just re-run to resume.
#   Rscript c1_vci_points_tilewise.R
# Out (partial, appended): out_files/Chapter1/tables/vci_points_partial.csv
# ==============================================================================
suppressPackageStartupMessages({ library(lidR); library(data.table) })
set_lidr_threads(1)                                   # single-thread; no parallel OOM
CTG  <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
PART <- "out_files/Chapter1/tables/vci_points_partial.csv"

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, pid := .I]; pts <- samp[, .(pid, x, y)]

done <- if (file.exists(PART)) fread(PART) else data.table(pid=integer(), vci_norm=double(), tile=character())
tiles <- list.files(CTG, pattern="\\.las$", full.names=TRUE)
cat(sprintf("%d pixels, %d tiles ; %d pixels déjà faits\n", nrow(pts), length(tiles), nrow(done)))

vci_pt <- function(las, x0, y0, r=25) {              # VCI in r-m circle around (x0,y0)
  sub <- filter_poi(las, (X-x0)^2 + (Y-y0)^2 <= r^2)
  if (is.null(sub) || npoints(sub) < 200) return(NA_real_)
  z <- sub$Z; z <- z[is.finite(z) & z >= 0]
  if (length(z) < 200) return(NA_real_)
  zmax <- max(z); if (zmax < 2) return(NA_real_)
  tryCatch(VCI(z, zmax = zmax), error=function(e) NA_real_)
}

for (i in seq_along(tiles)) {
  h  <- readLASheader(tiles[i]); bb <- h@PHB
  own <- pts[x >= bb[["Min X"]] & x <= bb[["Max X"]] & y >= bb[["Min Y"]] & y <= bb[["Max Y"]]]
  own <- own[!pid %in% done$pid]
  if (!nrow(own)) next
  las <- tryCatch(readLAS(tiles[i], select="xyzc"), error=function(e) NULL)
  if (is.null(las) || is.empty(las)) next
  las <- tryCatch(normalize_height(las, tin()), error=function(e) NULL)
  if (is.null(las)) next
  res <- own[, .(pid, vci_norm = vapply(seq_len(.N), function(k) vci_pt(las, x[k], y[k]), 0.0),
                 tile = basename(tiles[i]))]
  fwrite(res, PART, append = file.exists(PART))
  done <- rbind(done, res)
  cat(sprintf("[%d/%d] %s : +%d pixels (total %d/%d)\n",
              i, length(tiles), basename(tiles[i]), nrow(res), nrow(done), nrow(pts)))
  rm(las); gc(FALSE)
}
cat(sprintf("FINI compute VCI : %d/%d pixels\n", nrow(done), nrow(pts)))
