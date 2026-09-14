# ==============================================================================
# Radius sensitivity (SQUARE) — clip_rectangle, side = 2*radius around each HOBO,
# LAD + traits, run MuSICA. Same core as c1_radius_circle.R but square footprint
# (matches the existing radius_test_square/, so the 7 done radii are reused; this
# only needs the missing 12.5 and 40). Headless, resumable.
# Usage:  Rscript c1_radius_square.R <radius_m> [nmax_plots]
# Out:    out_files/radius_test_square/<r>m/musica_out_Blois_pt_<id>_radius_<r>.nc
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(dplyr); library(lubridate)
  library(musica.tools); library(rmusica)
})
args <- commandArgs(trailingOnly = TRUE)
radius <- as.numeric(args[1]); nmax <- if (length(args) >= 2) as.integer(args[2]) else Inf

ctg_path     <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
geojson_path <- "in_files/data_Blois_utm31n.geojson"
forcing_file <- "./in_files/musica_in_Blois.nc"
musica.cmd   <- "bash -i -c musica"
base_out_dir <- "./out_files/radius_test_square"

ctg <- readLAScatalog(ctg_path)
plots_sf <- st_read(geojson_path, quiet = TRUE)
df_plots <- plots_sf %>% st_coordinates() %>% as.data.frame() %>% rename(x = X, y = Y) %>%
  mutate(id_plot = plots_sf$id_plot)

radius_dir <- file.path(base_out_dir, paste0(radius, "m")); dir.create(radius_dir, recursive = TRUE, showWarnings = FALSE)
np <- min(nrow(df_plots), nmax)
cat(sprintf("RADIUS %s m (square, side=%s m) — %d plots -> %s\n", radius, 2*radius, np, radius_dir))

for (i in seq_len(np)) {
  id <- df_plots$id_plot[i]
  hf <- file.path(radius_dir, sprintf("musica_out_Blois_pt_%s_radius_%s.nc", id, radius))
  if (file.exists(hf)) next
  xc <- df_plots$x[i]; yc <- df_plots$y[i]
  pt <- tryCatch(clip_rectangle(ctg, xc - radius, yc - radius, xc + radius, yc + radius), error = function(e) NULL)
  if (is.null(pt) || npoints(pt) < 10) { message(sprintf(" skip %s (pts)", id)); next }
  dtm <- rasterize_terrain(pt, res = 1, algorithm = tin())
  pt_norm <- normalize_height(pt, dtm)
  lad_data <- LAD(pt_norm$Z, z0 = 1); allometry <- data.frame(height = lad_data$z, density = lad_data$lad)
  pai_val <- round(2 * sum(allometry$density, na.rm = TRUE), 2)
  h_max <- max(pt_norm$Z, na.rm = TRUE)
  chm <- rasterize_canopy(pt_norm, res = 0.5, pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1)))
  vp <- sum(!is.na(values(chm))); clumping_30h <- if (vp > 0) round(sum(values(chm) > 0.30 * h_max, na.rm = TRUE) / vp, 3) else 0
  pheno <- calc_phenology(list.year = c(2021, 2022), nleafage = 1, budburst_date = 115, leaf_age_max_in = 0.56,
                          relative_age_firstmax = 0.10, relative_age_lastmax = 0.75, LAI_max_per_cohort = pai_val)
  tryCatch(callmusica(
    musica.param = list("setupctl" = list("clumping_factor" = clumping_30h, "forcing_filename" = forcing_file,
                                          "history_filename" = hf, "forcing_height" = h_max + 2)),
    leaf.param = list("musica_veg1" = list("phenology" = pheno, "allometry" = allometry,
                                           "leafphenologyctl" = list("lai_max_per_cohort" = pai_val),
                                           "leafallometryctl" = list("canopy_height_top" = h_max, "canopy_height_bottom" = 2),
                                           "leafmusicactl" = list("canopy_height_top" = h_max))),
    musica.cmd = musica.cmd, keep.tmp = FALSE, out.netcdf = TRUE),
    error = function(e) message(sprintf(" MuSICA err %s: %s", id, e$message)))
  if (file.exists(hf)) message(sprintf(" done %s r%s", id, radius))
}
cat(sprintf("RADIUS %s DONE\n", radius))
