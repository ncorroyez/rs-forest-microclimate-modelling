# ==============================================================================
# Re-run MuSICA REF (1111 = all real) for the 53 HOBO sensors with the
# corrected Hmax (buffer-MAX 25 m on the raw 10 m raster).
#
# Goal : validate that the Hmax fix resolves the systematic ΔTmax > 0 bias.
# Compare new ΔTmax with the legacy values from
#   out_files/musica_hobo_validation/H1f_4_Full_real/
#
# Output :
#   out_files/musica_hobo_hmaxfix/REF/musica_out_HOBO_<id>.nc  (53 files)
#   outputs/hmax_fix_validation.csv  (per-sensor comparison)
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

cli_h1("Re-run MuSICA REF on 53 HOBO sensors with Hmax fix")

# ---- 1. Build HOBO inputs with corrected Hmax (buffer-MAX 25 m) -------------
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.table(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
  hmax_raw_raster = file.path(CFG$in_dir, "max_res_10_m.tif"),
  hmax_buffer = 25
))
# Floor fCover at 0.5 to match the floor05 convention used in the deck
df_hobo[fCover < 0.5, fCover := 0.5]
cli_alert("n HOBO inputs : {nrow(df_hobo)}")
cli_alert("Hmax (corrected) mean = {round(mean(df_hobo$Hmax), 2)} m, "
              ~ "range [{round(min(df_hobo$Hmax), 2)}, {round(max(df_hobo$Hmax), 2)}]")

# ---- 2. Build REF scenario (1111 : all real) ---------------------------------
ref_scenario <- list(
  name      = "REF",
  bit_code  = "1111",
  lai_fn    = function(pr) as.numeric(pr$LAI),
  hmax_fn   = function(pr) as.numeric(pr$Hmax),
  fcover_fn = function(pr) as.numeric(pr$fCover),
  lad_fn    = make_lad_real
)

# ---- 3. Run MuSICA for each HOBO sensor -------------------------------------
out_dir <- here::here("out_files/musica_hobo_hmaxfix/REF")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
t_start <- Sys.time()

df_hobo_df <- as.data.frame(df_hobo)
for (i in seq_len(nrow(df_hobo_df))) {
  plot_row <- df_hobo_df[i, , drop = FALSE]
  sim_id   <- sprintf("HOBO_%s", plot_row$id_plot)
  out_nc   <- file.path(out_dir, paste0("musica_out_", sim_id, ".nc"))
  if (file.exists(out_nc) && file.size(out_nc) > 1e6) {
    cat(sprintf("  [%d/%d] %s : already done (skip)\n",
                i, nrow(df_hobo_df), plot_row$id_plot))
    next
  }
  cat(sprintf("  [%d/%d] %s : Hmax=%.2fm, LAI=%.2f, fCover=%.2f\n",
              i, nrow(df_hobo_df), plot_row$id_plot,
              plot_row$Hmax, plot_row$LAI, plot_row$fCover))
  res <- tryCatch(
    run_musica_one(plot_row, ref_scenario, out_nc,
                    CFG$forcing_file, CFG$musica_cmd),
    error = function(e) {
      cli_alert_warning("    ERROR : {e$message}")
      NULL
    })
}
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
cli_alert_success("All sims done in {round(elapsed, 1)} min")

# ---- 4. Compute ΔTmax and compare with legacy --------------------------------
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

compute_delta <- function(nc_dir, ids) {
  out <- list()
  for (sid in ids) {
    f <- file.path(nc_dir, sprintf("musica_out_HOBO_%s.nc", sid))
    if (!file.exists(f)) next
    res <- tryCatch(
      extract_deltatmax_one(f, df_macro, CFG$date_seq,
                              z_target = CFG$tair_target_height),
      error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    out[[sid]] <- data.table(id_plot = sid,
                                 Tmax_micro_mean = mean(res$Tmax_micro, na.rm = TRUE),
                                 dTmax_mean = mean(res$Delta_Tmax, na.rm = TRUE))
  }
  rbindlist(out)
}

cli_alert("Computing Delta_Tmax (legacy + corrected)...")
df_legacy <- compute_delta(
  here::here("out_files/musica_hobo_validation/H1f_4_Full_real"),
  df_hobo$id_plot)
df_new <- compute_delta(out_dir, df_hobo$id_plot)

# Merge
m <- merge(df_legacy[, .(id_plot, Tmax_legacy = Tmax_micro_mean,
                            dT_legacy = dTmax_mean)],
            df_new[, .(id_plot, Tmax_new = Tmax_micro_mean,
                          dT_new = dTmax_mean)],
            by = "id_plot")
m[, Hmax := df_hobo[match(id_plot, df_hobo$id_plot), Hmax]]
m[, dT_change := dT_new - dT_legacy]
setorder(m, dT_change)

cli_h2("Per-sensor comparison (legacy vs corrected Hmax)")
print(m)

# Summary
cli_h2("Summary")
cat(sprintf("Tmax_micro mean : legacy = %.2f °C, corrected = %.2f °C, diff = %+.2f °C\n",
            mean(m$Tmax_legacy), mean(m$Tmax_new),
            mean(m$Tmax_new - m$Tmax_legacy)))
cat(sprintf("Delta_Tmax mean : legacy = %+.2f °C, corrected = %+.2f °C\n",
            mean(m$dT_legacy), mean(m$dT_new)))
cat(sprintf("Sensors with dT_new < 0 (buffering)  : %d / %d (was 0 / 53 in legacy)\n",
            sum(m$dT_new < 0), nrow(m)))
cat(sprintf("Sensors with dT_new > 0 (amplifying) : %d / %d (was 53 / 53 in legacy)\n",
            sum(m$dT_new > 0), nrow(m)))

# Save CSV
fwrite(m, here::here("outputs/hmax_fix_validation.csv"))
cli_alert_success("Saved comparison to outputs/hmax_fix_validation.csv")
