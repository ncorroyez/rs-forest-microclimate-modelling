# ==============================================================================
# Recompute per-plot LAD/LAI profiles at FINE resolution for field plots.
# Same method as current (lidR::LAD, MacArthur-Horn, k=0.5) but:
#   z0 (h_min) = 0.5 m   (was 1)     → captures the 0.5–1 m understory stratum
#   dz (bins)  = 0.5 m   (was 1)
#   profiles padded/normalised to 40 m max  → grid z = 0.75, 1.25, … 39.75 (80 bins)
#   points < z0 = ground (0), excluded from LAD (lidR handles via z0)
# Footprint per plot: circular clip, radius R (default 25 m = build_hobo "buffer25").
# Source: leaf_on / 3-las_normalized_utm (height-normalised summer clouds).
# Sites: Blois, Mormal (Aigoual pending its field-plot file).
# Output: in_files/lad_z05/<site>_lad_z05_r<R>.csv  (id_plot,x,y,LAI,Hmax,LAD_Layer_*)
# ==============================================================================

suppressMessages({ library(lidR); library(sf); library(data.table); library(here) })
options(lidR.progress = FALSE)           # silence clip progress bars

R_FOOT <- 25                    # clip radius (m) — match current buffer25; change freely
# Parameterisable: default z0.5/dz0.5; set env LAD_Z0/LAD_DZ/LAD_SUF for z1-raw etc.
Z0  <- as.numeric(Sys.getenv("LAD_Z0", "0.5"))
DZ  <- as.numeric(Sys.getenv("LAD_DZ", "0.5"))
SUF <- Sys.getenv("LAD_SUF", "z05")
# Scan-angle path-length correction (Eq. 2 of the S2/LiDAR paper): OFF by default so the
# baseline output is unchanged (a diagnostic sec_theta column is added but LAD/LAI are not
# touched). Set SCAN_ANGLE_CORR=TRUE to divide LAD/LAI by the per-plot <sec theta> and write
# to a separate *_sacorr file, keeping the uncorrected baseline for comparison.
# UNTESTED: the raw LAS (MyPassport drive) were offline when this was written; validate the
# per-plot LAD against lidR before trusting the corrected profiles.
SEC_CORR <- as.logical(Sys.getenv("SCAN_ANGLE_CORR", "FALSE"))
if (isTRUE(SEC_CORR)) SUF <- paste0(SUF, "_sacorr")
KEXT <- 0.5; HMAX_CAP <- 40
z_grid <- seq(Z0 + DZ/2, HMAX_CAP - DZ/2, by = DZ)
lyr_names <- sprintf("LAD_Layer_%.2f", z_grid)
OUT <- here::here("in_files/lad_z05"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
message(sprintf("LAD config: z0=%.2f dz=%.2f suffix=%s (%d bins %.2f-%.2f)",
                Z0, DZ, SUF, length(z_grid), z_grid[1], tail(z_grid,1)))

DRIVE <- "/media/corroyez/MyPassport/01_DATA"
SITES <- list(
  # Use RAW 2-las_utm + on-the-fly normalisation (preserves sub-1.5 m returns that
  # the pre-normalised 3-las_normalized_utm cloud had been filtered out) → real
  # understory capture at z0=0.5. Normalisation: knnidw on classified ground.
  Blois  = list(ctg = file.path(DRIVE,"Blois/LiDAR/leaf_on/2-las_utm"),
                plots = here::here("in_files/data_Blois_utm31n.geojson"), normalize = TRUE),
  Aigoual = list(ctg = file.path(DRIVE,"Aigoual/LiDAR/leaf_on/2-las_utm"),
                plots = "/home/corroyez/Documents/NC_Full/01_DATA/Aigoual/Geo_Files/data_utm31n.geojson",
                normalize = TRUE),
  Mormal = list(ctg = file.path(DRIVE,"Mormal/LiDAR/leaf_on/2-las_utm"),
                plots = file.path(DRIVE,"Mormal/Geo_Files/data_utm31n.geojson"), normalize = TRUE)
)

lad_one_plot <- function(ctg, x, y, normalize = FALSE) {
  las <- tryCatch(suppressMessages(clip_circle(ctg, x, y, R_FOOT)), error = function(e) NULL)
  if (is.null(las) || is.empty(las)) return(NULL)
  las <- tryCatch(filter_duplicates(las), error = function(e) las)   # overlapping tiles
  if (nrow(las@data) < 50) return(NULL)
  if (normalize) {
    las <- tryCatch(suppressMessages(normalize_height(las, knnidw())), error = function(e) NULL)
    if (is.null(las)) return(NULL)
  }
  fin <- is.finite(las@data$Z) & las@data$Z >= 0
  Z <- las@data$Z[fin]
  if (length(Z) < 50) return(NULL)
  # <sec theta> over the same returns: l = dz * <sec theta> (Eq. 2), so LAD/LAI scale by
  # 1/<sec theta>. ScanAngle (LAS 1.4) or ScanAngleRank (older), both in degrees.
  sa_all <- if ("ScanAngle" %in% names(las@data)) las@data$ScanAngle
            else if ("ScanAngleRank" %in% names(las@data)) las@data$ScanAngleRank else NULL
  sa <- if (!is.null(sa_all)) sa_all[fin] else numeric(0)
  sec_bar <- if (length(sa) && any(is.finite(sa))) mean(1/cos(sa[is.finite(sa)] * pi/180)) else NA_real_
  scale <- if (isTRUE(SEC_CORR) && is.finite(sec_bar)) 1/sec_bar else 1
  d <- tryCatch(lidR::LAD(Z, dz = DZ, k = KEXT, z0 = Z0), error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  prof <- setNames(numeric(length(z_grid)), lyr_names)            # 0 above canopy
  idx <- match(round(d$z, 2), round(z_grid, 2))
  ok <- !is.na(idx) & is.finite(d$lad)
  prof[idx[ok]] <- d$lad[ok] * scale
  list(LAI = sum(prof, na.rm = TRUE) * DZ,                        # ∫LAD dz
       Hmax = min(max(Z, na.rm = TRUE), HMAX_CAP),
       prof = prof, n = length(Z), sec_bar = sec_bar)
}

for (site in names(SITES)) {
  cfg <- SITES[[site]]
  fout <- file.path(OUT, sprintf("%s_lad_%s_r%d.csv", site, SUF, R_FOOT))
  if (file.exists(fout)) { message(sprintf("[%s] output exists — skip", site)); next }
  n_las <- length(list.files(cfg$ctg, pattern = "\\.la[sz]$", ignore.case = TRUE))
  if (!dir.exists(cfg$ctg) || n_las == 0) { message(sprintf("[%s] no LAS in %s — skip", site, cfg$ctg)); next }
  message(sprintf("\n==== %s : reading catalog (%d tiles, normalize=%s) ====", site, n_las, cfg$normalize))
  ctg <- readLAScatalog(cfg$ctg); opt_progress(ctg) <- FALSE
  g <- st_read(cfg$plots, quiet = TRUE); g <- st_transform(g, 32631)
  xy <- st_coordinates(g); ids <- if ("id_plot" %in% names(g)) g$id_plot else seq_len(nrow(g))
  message(sprintf("[%s] %d plots ; radius %dm", site, nrow(g), R_FOOT))

  rows <- vector("list", nrow(g)); t0 <- Sys.time()
  for (i in seq_len(nrow(g))) {
    r <- lad_one_plot(ctg, xy[i,1], xy[i,2], normalize = cfg$normalize)
    if (is.null(r)) { message(sprintf("  [%s] plot %s : empty/insufficient", site, ids[i])); next }
    rows[[i]] <- data.table(id_plot = ids[i], x = xy[i,1], y = xy[i,2],
                            LAI = r$LAI, Hmax = r$Hmax, n_points = r$n, sec_theta = r$sec_bar,
                            t(setNames(r$prof, lyr_names)))
    if (i %% 10 == 0) message(sprintf("  ... %d/%d", i, nrow(g)))
  }
  DT <- rbindlist(rows, fill = TRUE)
  fout <- file.path(OUT, sprintf("%s_lad_%s_r%d.csv", site, SUF, R_FOOT))
  fwrite(DT, fout)
  message(sprintf("[%s] DONE %d/%d plots in %.1f min → %s",
                  site, nrow(DT), nrow(g), as.numeric(difftime(Sys.time(),t0,units="mins")), basename(fout)))
  message(sprintf("[%s] LAI range [%.2f, %.2f] ; Hmax range [%.1f, %.1f]",
                  site, min(DT$LAI,na.rm=T), max(DT$LAI,na.rm=T), min(DT$Hmax,na.rm=T), max(DT$Hmax,na.rm=T)))
}
