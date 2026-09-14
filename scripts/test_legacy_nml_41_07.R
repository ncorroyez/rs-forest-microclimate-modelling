# ==============================================================================
# Test : run 41_07 REF with the LEGACY musica.nml (`.musica_2ce74bb751598/`),
# but with HISTORY_VARIABLES replaced by current-compatible values, and
# the corrected Hmax (buffer-MAX 25 m).
#
# Backup-restore : the current `musica.nml` is moved to `musica.nml.before_test`
# before the run. At exit (success OR error) we restore it.
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

cli_h1("Test : 41_07 with LEGACY musica.nml + corrected Hmax")

# ---- 1. Backup current musica.nml -------------------------------------------
backup_path <- here::here("musica.nml.before_test")
file.copy(here::here("musica.nml"), backup_path, overwrite = TRUE)
md5_orig <- tools::md5sum(here::here("musica.nml"))
cli_alert("Backup at {.path {basename(backup_path)}}, md5={substr(md5_orig, 1, 8)}…")

# Restore on exit
on.exit({
  if (file.exists(backup_path)) {
    file.copy(backup_path, here::here("musica.nml"), overwrite = TRUE)
    if (tools::md5sum(here::here("musica.nml")) == md5_orig) {
      cli_alert_success("musica.nml restored (md5 OK)")
    } else {
      cli_alert_danger("musica.nml restore FAILED, backup is at {.path {basename(backup_path)}}")
    }
  }
}, add = TRUE)

# ---- 2. Build modified musica.nml ------------------------------------------
legacy_nml <- here::here(".musica_2ce74bb751598/musica.nml")
lines <- readLines(legacy_nml)
# Replace HISTORY_VARIABLES with current-compatible
old_hv <- grep('^[[:space:]]*HISTORY_VARIABLES[[:space:]]*=', lines)
if (length(old_hv) > 0) {
  lines[old_hv] <- '  HISTORY_VARIABLES = "layer_thickness", "z_soil", "dz_soil", "t_soil", "w_soil", "t_air", "w_air", "wind"'
}
# Remove HISTORY_DAILY_VARIABLES (replace with empty)
old_hdv <- grep('^[[:space:]]*HISTORY_DAILY_VARIABLES[[:space:]]*=', lines)
if (length(old_hdv) > 0) {
  lines <- lines[-old_hdv]
}
# Set FORCING_FILENAME to match current ERA5
old_ff <- grep('^[[:space:]]*FORCING_FILENAME[[:space:]]*=', lines)
if (length(old_ff) > 0) {
  lines[old_ff] <- '  FORCING_FILENAME = "in_files/musica_in_Blois.nc"'
}
# Write
writeLines(lines, here::here("musica.nml"))
cli_alert("Wrote legacy-derived musica.nml ({length(lines)} lines)")

# ---- 3. Prepare 41_07 input -------------------------------------------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
  hmax_raw_raster = file.path(CFG$in_dir, "max_res_10_m.tif"),
  hmax_buffer = 25))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
plot_row <- df_hobo[df_hobo$id_plot == "41_07", , drop = FALSE]
cli_alert("41_07 : LAI={round(plot_row$LAI,2)}, Hmax={round(plot_row$Hmax,2)}, "
            ~ "fCover={round(plot_row$fCover,2)}")

ref_scenario <- list(
  name      = "REF", bit_code = "1111",
  lai_fn    = function(pr) as.numeric(pr$LAI),
  hmax_fn   = function(pr) as.numeric(pr$Hmax),
  fcover_fn = function(pr) as.numeric(pr$fCover),
  lad_fn    = make_lad_real)

# ---- 4. Run ----------------------------------------------------------------
out_dir <- here::here("out_files/musica_hobo_test_legacynml")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_nc <- file.path(out_dir, "musica_out_HOBO_41_07.nc")
if (file.exists(out_nc)) file.remove(out_nc)

cli_alert("Running MuSICA (1 sim)...")
t_start <- Sys.time()
res <- tryCatch(
  run_musica_one(plot_row, ref_scenario, out_nc,
                  CFG$forcing_file, CFG$musica_cmd),
  error = function(e) { cli_alert_danger("ERROR: {e$message}"); NULL })
cli_alert("Done in {round(as.numeric(difftime(Sys.time(), t_start, units='secs')), 1)} s")

# ---- 5. Verify and report --------------------------------------------------
if (file.exists(out_nc) && file.size(out_nc) > 1e6) {
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  res_d <- extract_deltatmax_one(out_nc, df_macro, CFG$date_seq,
                                    z_target = CFG$tair_target_height)
  cli_h2("Result")
  cat(sprintf("  Tmax_micro mean  = %.3f °C\n", mean(res_d$Tmax_micro, na.rm = TRUE)))
  cat(sprintf("  Tmax_macro mean  = %.3f °C\n", mean(res_d$Tmax_macro, na.rm = TRUE)))
  cat(sprintf("  Delta_Tmax mean  = %+.3f °C\n", mean(res_d$Delta_Tmax, na.rm = TRUE)))
  cli_h2("Comparison")
  cat("  Apr 2026 H1f_4  current nml Hmax=7.5    : dT = +3.44 °C\n")
  cat("  Today hmaxfix  current nml Hmax=11.09   : dT = +0.28 °C\n")
  cat(sprintf("  Today legacy nml  Hmax=%.2f             : dT = %+.2f °C  (this run)\n",
              plot_row$Hmax, mean(res_d$Delta_Tmax, na.rm = TRUE)))
  cat("  TARGET (Oct 2025 Safran)                : dT = -0.55 °C  (buffer)\n")
} else {
  cli_alert_danger("NC file not produced — MuSICA crashed")
}
