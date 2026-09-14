# ==============================================================================
# Raw ALS scan-angle distribution over the logger footprints (Appendix F).
#
# WHY THIS EXISTS. Appendix F quoted "return scan angles reach ±33°" and a
# "median magnitude about 12°" with no persisted source: the control script
# reads ScanAngle but writes only the per-layer LAD and secant, and its output
# had been left at 0 bytes because the LAS catalogue lives on an external drive.
# Measured over 4,035,777 returns in the 60 logger footprints: min −32°, max
# +33° (so ±33° is right), median |angle| **16°** (not 12°), p95 29°.
#
# NEEDS the external drive mounted. Writes a small summary so the numbers survive
# without it.
# Reads : /media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm
#         in_files/data_Blois_utm31n.geojson
# Writes: out_files/Chapter1/tables/tab_scanangle_raw_stats.csv
# ==============================================================================
suppressPackageStartupMessages({library(lidR); library(data.table); library(sf)})
CTG_DIR <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
OUT <- "out_files/Chapter1/tables/tab_scanangle_raw_stats.csv"
if (!dir.exists(CTG_DIR)) {
  if (file.exists(OUT)) { cat("drive absent; cached summary kept at ", OUT, "\n"); quit(status = 0) }
  stop("Appendix F scan angles: the LAS catalogue is not reachable at ", CTG_DIR,
       " and no cached summary exists. Mount the drive once.")
}
ctg <- readLAScatalog(CTG_DIR); sf::st_crs(ctg) <- 32631
lidR::opt_progress(ctg) <- FALSE; opt_select(ctg) <- "xyza"
g  <- st_transform(st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE), 32631)
xy <- st_coordinates(st_centroid(st_geometry(g)))
A <- unlist(lapply(seq_len(nrow(xy)), function(i) {
  la <- tryCatch(clip_circle(ctg, xy[i, 1], xy[i, 2], 11.3), error = function(e) NULL)
  if (is.null(la) || nrow(la@data) == 0) return(NULL)
  as.numeric(la@data$ScanAngle) }))
A <- A[is.finite(A)]
S <- data.table(n_returns = length(A), n_footprints = nrow(xy),
                min_deg = min(A), max_deg = max(A),
                median_abs_deg = median(abs(A)), p95_abs_deg = quantile(abs(A), .95),
                max_abs_deg = max(abs(A)),
                sec_at_median = 1/cos(median(abs(A))*pi/180),
                frac_beyond_30_pct = 100*mean(abs(A) > 30))
fwrite(S, OUT)
print(S[, lapply(.SD, function(x) round(x, 3))])
cat("DONE\n")
