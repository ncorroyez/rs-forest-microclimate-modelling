# ==============================================================================
# PIPELINE STAGE 01 — Data processing: LiDAR rasters → forest pixels → FPCA.
# Builds the structural dataframe and the LAD-shape FPCA used for typology.
# Idempotent: skips if cache exists. Cache: OUT_DATA/forest_fpca.rds
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
suppressMessages({ library(terra); library(fda) })   # fda: create.bspline.basis for FPCA
cli_h1("STAGE 01 — data processing (rasters → forest df → FPCA)")

cache <- file.path(PIPE$OUT_DATA, "forest_fpca.rds")
if (file.exists(cache)) { cli_alert_success("forest_fpca.rds present — skip (delete to rebuild)."); }  else {
  rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
  fb  <- build_forest_dataframe(rasters$stack)        # list(df, mat_lad, z_breaks)
  fdf <- fb$df
  cli_alert("Forest pixels: {nrow(fdf)}")
  fpca <- tryCatch(compute_fpca(fb$mat_lad, hmax_vec = fdf$Hmax, z_breaks = fb$z_breaks, n_harm = 3),
                   error = function(e) { cli_alert_warning("FPCA skipped: {e$message}"); NULL })
  saveRDS(list(forest = fdf, fpca = fpca, z_break = fb$z_breaks), cache)
  cli_alert_success("Saved forest_fpca.rds ({nrow(fdf)} pixels, FPCA {if(is.null(fpca)) 'failed' else 'ok'}).")
}
