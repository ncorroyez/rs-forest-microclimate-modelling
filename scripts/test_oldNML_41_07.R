# ==============================================================================
# Test : re-run 41_07 REF with the OLD musica.nml (from .musica_2ce74bb751598/)
# and corrected Hmax (buffer-MAX 25 m), to isolate the contribution of the
# namelist drift to the residual ΔTmax bias.
#
# Safety : the current `musica.nml` is backed up to `musica.nml.pretest_<ts>`
# before being temporarily swapped. On exit (success OR error), we restore
# from the backup and verify md5 matches the original.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))

cli_h1("Test : 41_07 REF with old musica.nml + corrected Hmax")

# ---- 1. Backup current musica.nml --------------------------------------------
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
backup_path <- here::here(sprintf("musica.nml.pretest_%s", ts))
file.copy(here::here("musica.nml"), backup_path, overwrite = FALSE)
md5_original <- tools::md5sum(here::here("musica.nml"))
cli_alert("Backup: {.path {basename(backup_path)}}")
cli_alert("Original md5: {md5_original}")

# Setup restore on exit
restore_nml <- function() {
  if (file.exists(backup_path)) {
    file.copy(backup_path, here::here("musica.nml"), overwrite = TRUE)
    md5_now <- tools::md5sum(here::here("musica.nml"))
    if (md5_now == md5_original) {
      cli_alert_success("Restored musica.nml from backup (md5 matches original)")
    } else {
      cli_alert_danger("Restore failed! md5 mismatch — original is at {backup_path}")
    }
  }
}
on.exit(restore_nml(), add = TRUE)

# ---- 2. Swap to old musica.nml ----------------------------------------------
old_nml_src <- here::here(".musica_2ce74bb751598/musica.nml")
if (!file.exists(old_nml_src)) stop("Old musica.nml not found at ", old_nml_src)
file.copy(old_nml_src, here::here("musica.nml"), overwrite = TRUE)
cli_alert("Swapped musica.nml -> old snapshot version (.musica_2ce74bb751598)")

# ---- 3. Build HOBO inputs with corrected Hmax for 41_07 ---------------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
  hmax_raw_raster = file.path(CFG$in_dir, "max_res_10_m.tif"),
  hmax_buffer = 25))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5

plot_row <- df_hobo[df_hobo$id_plot == "41_07", , drop = FALSE]
cli_alert("41_07 inputs: LAI={round(plot_row$LAI,2)}, "
              ~ "Hmax={round(plot_row$Hmax,2)}m, fCover={round(plot_row$fCover,2)}")

# ---- 4. Build REF scenario ---------------------------------------------------
ref_scenario <- list(
  name      = "REF",
  bit_code  = "1111",
  lai_fn    = function(pr) as.numeric(pr$LAI),
  hmax_fn   = function(pr) as.numeric(pr$Hmax),
  fcover_fn = function(pr) as.numeric(pr$fCover),
  lad_fn    = make_lad_real)

# ---- 5. Run MuSICA ----------------------------------------------------------
out_dir <- here::here("out_files/musica_hobo_test_oldNML")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_nc <- file.path(out_dir, "musica_out_HOBO_41_07.nc")
if (file.exists(out_nc)) file.remove(out_nc)
cli_alert("Running MuSICA (1 sim, ~30 s)...")
t_start <- Sys.time()
res <- tryCatch(
  run_musica_one(plot_row, ref_scenario, out_nc,
                  CFG$forcing_file, CFG$musica_cmd),
  error = function(e) {
    cli_alert_danger("ERROR : {e$message}")
    NULL
  })
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
cli_alert_success("Done in {round(elapsed, 1)} s")

# ---- 6. Compute new ΔTmax and compare ----------------------------------------
if (file.exists(out_nc) && file.size(out_nc) > 1e6) {
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  res_d <- extract_deltatmax_one(out_nc, df_macro, CFG$date_seq,
                                    z_target = CFG$tair_target_height)
  cli_h2("Result")
  cat(sprintf("  Tmax_micro mean : %.3f °C\n", mean(res_d$Tmax_micro, na.rm = TRUE)))
  cat(sprintf("  Tmax_macro mean : %.3f °C\n", mean(res_d$Tmax_macro, na.rm = TRUE)))
  cat(sprintf("  Delta_Tmax mean : %+.3f °C\n", mean(res_d$Delta_Tmax, na.rm = TRUE)))
  cli_h2("Comparison table")
  cat("                            41_07  Tmax_micro  Delta_Tmax\n")
  cat("  Oct 2025 ERA5  Hmax 10.5   ->     25.34       +1.26\n")
  cat("  Apr 2026 H1f_4  Hmax 7.5   ->     27.52       +3.44 (current pipeline)\n")
  cat(sprintf("  Today new nml  Hmax 11.09  ->     %5.2f       %+5.2f (this test)\n",
              mean(res_d$Tmax_micro, na.rm = TRUE),
              mean(res_d$Delta_Tmax, na.rm = TRUE)))
}
