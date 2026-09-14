# ==============================================================================
# Test which `musica.nml` parameter drives the residual ΔTmax bias on 41_07.
# Strategy : keep the current `musica.nml` on disk untouched. For each
# candidate parameter that changed between Oct 2025 and now, run MuSICA with
# extra_setup overriding just that one parameter and measure ΔTmax.
#
# Tested parameters :
#   - FORCING_TIMESTEP : 3600 (current) vs 1800 (old snapshot)
#   - fCA_soil         : 30 (current) vs absent (old)
#   - WATER_VOLUME_PERLAI : c(10,10) (current) vs absent (old)
#   - CO2_H2O_EQ_DEGREE_INLEAF : c(1.0,1.0) (current) vs absent (old)
#
# Note : "absent" in extra_setup means we still send the current value, since
# we can't easily unset a param via extra_setup. So this is a one-direction
# sensitivity test : how much does the parameter VALUE matter, not the
# presence of the line itself.
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

cli_h1("Test : isolate which .nml parameter drives the residual bias on 41_07")

# ---- 1. Build 41_07 input with corrected Hmax ------------------------------
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

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

# ---- 2. Helper to run + extract ΔTmax ---------------------------------------
out_dir <- here::here("out_files/musica_hobo_test_nmlparams")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_variant <- function(label, extra) {
  out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_41_07_%s.nc", label))
  if (file.exists(out_nc)) file.remove(out_nc)
  cli_alert("[{label}] running...")
  res <- tryCatch(
    run_musica_one(plot_row, ref_scenario, out_nc,
                    CFG$forcing_file, CFG$musica_cmd,
                    extra_setup = extra),
    error = function(e) { cli_alert_danger("  ERROR: {e$message}"); NULL })
  if (!file.exists(out_nc) || file.size(out_nc) < 1e6) {
    cli_alert_warning("  [{label}] failed")
    return(NULL)
  }
  d <- extract_deltatmax_one(out_nc, df_macro, CFG$date_seq,
                                z_target = CFG$tair_target_height)
  list(label = label,
        Tmax_micro = mean(d$Tmax_micro, na.rm = TRUE),
        dTmax = mean(d$Delta_Tmax, na.rm = TRUE))
}

# ---- 3. Run variants --------------------------------------------------------
variants <- list(
  "baseline"          = list(),  # current nml defaults
  "ftstep_1800"       = list(FORCING_TIMESTEP = 1800),
  "fCA_off"           = list(fCA_soil = 0),
  "fCA_3"             = list(fCA_soil = 3),
  "wvol_1"            = list(WATER_VOLUME_PERLAI = c(1, 1)),
  "wvol_100"          = list(WATER_VOLUME_PERLAI = c(100, 100)),
  "co2eq_off"         = list(CO2_H2O_EQ_DEGREE_INLEAF = c(0, 0)),
  "dalim_on"          = list(D_ALIM_FLAG = TRUE, D_ALIM = c(1, 1)))

results <- list()
for (nm in names(variants)) {
  results[[nm]] <- run_variant(nm, variants[[nm]])
}

# ---- 4. Report --------------------------------------------------------------
cli_h2("Results — 41_07 sensitivity to single-param changes")
cat(sprintf("\n%-22s  Tmax_micro  ΔTmax\n", "variant"))
cat("------------------------------------------------\n")
for (nm in names(results)) {
  r <- results[[nm]]
  if (is.null(r)) {
    cat(sprintf("%-22s  (failed)\n", nm))
  } else {
    cat(sprintf("%-22s  %6.2f °C   %+5.2f °C\n", nm, r$Tmax_micro, r$dTmax))
  }
}
cat("\nReferences :\n")
cat("  Oct 2025 ERA5 Hmax=10.5   : ΔTmax = +1.26 °C\n")
cat("  Apr 2026 H1f_4 Hmax=7.5   : ΔTmax = +3.44 °C\n")
cat("  Today HmaxFix Hmax=11.09  : ΔTmax = +0.28 °C  (target : reach < 0)\n")
