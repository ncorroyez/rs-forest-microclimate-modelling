# ==============================================================================
# Apply the scan-angle correction (Eq. 2, l = dz*<sec theta>) to the cLHS attribution
# sample: LAI_corr = LAI / <sec theta> per pixel, <sec theta> = mean(1/cos(ScanAngle)).
# Per-pixel <sec theta> from a 10 m mean-sec-theta raster (one catalog pass), extracted
# at the 400 cLHS (x,y). Global-mean fallback for pixels with no coverage. Cross-checked
# against the 60 field-plot sec_theta (Blois_lad_z05_sacorr_r25.csv). OVERWRITES the rds
# (original already backed up in _preScanAngle_backup). LAD_Layer_* untouched (shape is
# renormalised to LAI downstream in R/lad.R) — scale LAI ONLY (avoids double-correction).
# ==============================================================================
suppressMessages({ library(lidR); library(sf); library(data.table) })
options(lidR.progress = FALSE)
CTG_DIR <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
RDS <- "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"
RSEC <- 10                                       # sec-theta footprint (m) ~ cLHS pixel

# build spatial index (.lax) so clip_circle reads only the relevant tile (else it scans all 160)
las_files <- list.files(CTG_DIR, "\\.las$", full.names = TRUE)
todo <- las_files[!file.exists(sub("\\.las$", ".lax", las_files))]
if (length(todo)) { message("indexing ", length(todo), " tiles ..."); invisible(lapply(todo, function(f) tryCatch(rlas::writelax(f), error = function(e) NULL))) }

ctg <- readLAScatalog(CTG_DIR); opt_select(ctg) <- "xyza"; opt_progress(ctg) <- FALSE
# per-point <sec theta> by small clip (no normalize, no LAD — just the scan angle)
sec_at <- function(x, y) {
  las <- tryCatch(suppressMessages(clip_circle(ctg, x, y, RSEC)), error = function(e) NULL)
  if (is.null(las) || is.empty(las) || nrow(las@data) < 50) return(NA_real_)
  sa <- if ("ScanAngle" %in% names(las@data)) las@data$ScanAngle
        else if ("ScanAngleRank" %in% names(las@data)) las@data$ScanAngleRank else NULL
  if (is.null(sa)) return(NA_real_)
  sa <- sa[is.finite(sa)]; if (!length(sa)) return(NA_real_)
  mean(1/cos(sa * pi/180))
}

s <- as.data.table(readRDS(RDS))
message(sprintf("Computing <sec theta> at %d cLHS points (r=%dm) ...", nrow(s), RSEC))
s[, sec_theta := vapply(seq_len(.N), function(i) sec_at(x[i], y[i]), numeric(1))]
gmean <- mean(s$sec_theta, na.rm = TRUE)
nfill <- s[!is.finite(sec_theta), .N]
s[!is.finite(sec_theta), sec_theta := gmean]                 # global-mean fallback for empty pixels
cat(sprintf("cLHS pixels: %d ; sec_theta NA->global-fallback: %d ; global mean %.4f\n", nrow(s), nfill, gmean))
cat("cLHS sec_theta summary:\n"); print(round(summary(s$sec_theta), 4))
LAI0 <- s$LAI
s[, LAI := LAI / sec_theta]                                  # scan-angle correction (LAI only)
cat(sprintf("LAI change: mean %+.1f%% (sd %.1f%%) ; before mean %.2f -> after %.2f\n",
            100*(mean(s$LAI/LAI0)-1), 100*sd(s$LAI/LAI0), mean(LAI0), mean(s$LAI)))

# --- cross-check vs the 60 field-plot sec_theta (25 m clip, from the LAD recompute) ---
fp <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")[, .(id_plot, sec_theta_clip = sec_theta)]
g <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE); g <- st_transform(g, 32631)
gxy <- st_coordinates(g); gid <- if ("id_plot" %in% names(g)) g$id_plot else seq_len(nrow(g))
cc <- data.table(id_plot = gid, sec_r10 = vapply(seq_len(nrow(gxy)), function(i) sec_at(gxy[i,1], gxy[i,2]), numeric(1)))
cc <- merge(cc, fp, by = "id_plot")
cat(sprintf("\nCROSS-CHECK (field plots): r=10m point vs 25 m clip sec_theta, n=%d, r=%.3f, mean|diff|=%.4f\n",
            nrow(cc), cor(cc$sec_r10, cc$sec_theta_clip, use="complete.obs"),
            mean(abs(cc$sec_r10 - cc$sec_theta_clip), na.rm=TRUE)))

saveRDS(s, RDS)                                              # OVERWRITE (original backed up)
cat(sprintf("\nOVERWROTE %s with scan-angle-corrected LAI (baseline adoption)\n", RDS))
