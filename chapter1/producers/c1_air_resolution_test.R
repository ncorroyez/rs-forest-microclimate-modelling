# ==============================================================================
# ONE test simulation at AIR_RESOLUTION_LEVEL = 1 (20 vegetation / 30 air layers
# instead of 10 / 15), on the tallest logger plot, to measure where the lowest
# air node actually lands. Everything else is the chapter's baseline
# (c1_hobo_native20_validation.R): v3.2.3, abl_flag = 'iter', wind-corrected
# forcing per plot, scan-angle-corrected LAI, native-20 m traits, fCover floored.
#
# WHY. JO (2026-08-20) asked for a 0.1 m reading. With the current grid the
# lowest node is 0.0285 x Hmax = 1.083 m on the tallest plot, so the deepest
# height available on all 53 loggers is 1.083 m. The grid is geometric (thickness
# ratio 1.12, node 1 = half the first thickness), which predicts 0.0069 x Hmax
# = 0.26 m at 20 vegetation layers. This run checks that prediction rather than
# trusting the extrapolation.
#
# Writes into a scratch directory, NOT into any chapter output.
#   CH1_RUN_MUSICA=TRUE Rscript scripts/c1_air_resolution_test.R
# ==============================================================================
if (!identical(Sys.getenv("CH1_RUN_MUSICA"), "TRUE"))
  stop("gated stage: re-run with CH1_RUN_MUSICA=TRUE")
suppressPackageStartupMessages({
  library(ncdf4); library(data.table); library(sf); library(terra)
  library(rmusica); library(musica.tools) })
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]; invisible(lapply(src, source)); source("pipeline/00_config.R")
FORC <- "in_files/FR-Blo_2021_v2.nc"
MB   <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
OUT  <- file.path(Sys.getenv("SCRATCH", "/tmp"), "airres_test"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
source("R/wind_correction.R")

rasters <- load_lidar_rasters("in_files_native20", 1)
hob <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
                                       hobo_buffer_mode = "buffer25"))
hob[fCover < 0.5, fCover := 0.5]
sac <- fread("in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv")[, .(id_plot, sec_theta)]
hob <- merge(hob, sac, by = "id_plot", all.x = TRUE)
hob[, sec_theta := mean(hob$sec_theta, na.rm = TRUE)]
pr <- as.data.frame(hob[which.max(Hmax)])                     # the tallest logger plot
cat(sprintf("test plot %s : Hmax = %.2f m, LAI = %.2f, fCover = %.2f\n",
            pr$id_plot, pr$Hmax, pr$LAI, pr$fCover))

mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
    leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
    LAI_max_per_cohort = lai1))
  d <- ph[ph$year == 2020 & ph$Julian_day == 365, , drop = FALSE]; d$Julian_day <- 366; rbind(ph, d) }
lai <- pr$LAI / pr$sec_theta
sc  <- list(lai_fn = function(p) lai, hmax_fn = function(p) pr$Hmax,
            fcover_fn = function(p) pr$fCover, lad_fn = make_lad_real,
            phenology_fn = mk_phen(lai))
prc <- pr; prc$LAI <- lai
fw  <- windcorr_forcing(pr$Hmax)

for (lvl in c(2L, 1L)) {                                      # 2 = current, 1 = doubled grid
  nc <- file.path(OUT, sprintf("airres%d_%s.nc", lvl, pr$id_plot))
  cat(sprintf("\n--- AIR_RESOLUTION_LEVEL = %d ---\n", lvl))
  run_musica_one(prc, sc, nc, fw, MB,
                 extra_setup = list("abl_flag" = '"iter"', "air_resolution_level" = lvl))
  if (!file.exists(nc)) { cat("  no output\n"); next }
  o <- nc_open(nc); rh <- ncvar_get(o, "relative_height")
  H <- median(ncvar_get(o, "veget_height_top"), na.rm = TRUE); nc_close(o)
  cat(sprintf("  nair = %d | node 1 relative = %.5f -> %.3f m (Hmax %.2f m)\n",
              length(rh), rh[1], rh[1] * H, H))
  cat(sprintf("  first five nodes (m): %s\n", paste(round(rh[1:5] * H, 3), collapse = "  ")))
}
cat("\nDONE\n")
