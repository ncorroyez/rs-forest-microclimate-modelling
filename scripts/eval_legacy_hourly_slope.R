# ==============================================================================
# Compute the LEGACY (main_compare_era5_ign.R) buf/amp metric on:
#   1. Oct 2025 original NCs   (out_files/musica/ERA5/)
#   2. V8 new sims              (out_files/musica_hobo_v8_truelegacy/)
#
# Metric : hourly Tair_z (nair == 1) ~ hourly ERA5 Tair → linear slope per plot.
#          buffer (slope < 1) vs amplification (slope > 1).
# Expected : 45 buffer / 8 amplification (user's memory).
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))

extract_nair1_hourly <- function(nc_file, date_seq) {
  if (!file.exists(nc_file)) return(NULL)
  nc <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc), add = TRUE)
  raw <- try(musica.tools::get_variable(nc, "Tair_z"), silent = TRUE)
  if (inherits(raw, "try-error") || is.null(raw) || nrow(raw) == 0) return(NULL)
  df <- raw %>%
    filter(as.Date(time) %in% date_seq, nair == 1) %>%
    mutate(Tair_z = Tair_z - 273.15, time = floor_date(time, "hour"))
  if (nrow(df) == 0) return(NULL)
  df
}

# ERA5 hourly
era5_nc <- nc_open(CFG$forcing_file)
era5_tair <- ncvar_get(era5_nc, "Tair") - 273.15
era5_time <- force_utc_nc(CFG$forcing_file, "time")
nc_close(era5_nc)
era5_df <- data.frame(time = floor_date(era5_time, "hour"), Tair_era5 = era5_tair) %>%
  filter(as.Date(time) %in% CFG$date_seq) %>%
  distinct(time, .keep_all = TRUE)
cli_alert("ERA5 hourly rows : {nrow(era5_df)}")

compute_regime <- function(nc_dir, pattern, label) {
  files <- list.files(nc_dir, pattern = pattern, full.names = TRUE)
  cli_h2("{label} ({length(files)} NCs)")
  out <- list()
  for (f in files) {
    id <- sub(pattern, "\\1", basename(f))
    if (id %in% CFG$ids_to_remove) next
    df <- extract_nair1_hourly(f, CFG$date_seq)
    if (is.null(df)) next
    m <- merge(df, era5_df, by = "time")
    if (nrow(m) < 24) next
    fit <- lm(Tair_z ~ Tair_era5, data = m)
    sl <- coef(fit)[2]
    out[[id]] <- data.frame(id_plot = id, slope = sl, log_slope = log(abs(sl)),
                             dT_mean = mean(m$Tair_z - m$Tair_era5, na.rm = TRUE))
  }
  DT <- as.data.table(do.call(rbind, out))
  cat(sprintf("  Buffer (slope < 1)        : %d / %d\n",
              sum(DT$slope < 1), nrow(DT)))
  cat(sprintf("  Amplification (slope > 1) : %d / %d\n",
              sum(DT$slope > 1), nrow(DT)))
  cat(sprintf("  slope mean=%.4f median=%.4f range=[%.4f, %.4f]\n",
              mean(DT$slope), median(DT$slope), min(DT$slope), max(DT$slope)))
  cat(sprintf("  dT_mean mean=%+.3f median=%+.3f\n",
              mean(DT$dT_mean), median(DT$dT_mean)))
  DT
}

cli_h1("LEGACY hourly-slope metric (nair == 1 ~ ERA5 hourly)")

DT_oct <- compute_regime("out_files/musica/ERA5",
                          "musica_out_ERA5_Blois_2021_pt_(.+)\\.nc$",
                          "Oct 2025 original NCs")

DT_v8  <- compute_regime("out_files/musica_hobo_v8_truelegacy",
                          "musica_out_HOBO_(.+)\\.nc$",
                          "V8 new sims (true legacy binary)")

both <- merge(DT_oct, DT_v8, by = "id_plot", suffixes = c("_oct","_v8"))
cli_h2("Cross-check (common plots: {nrow(both)})")
cat(sprintf("  slope correlation : %.4f\n", cor(both$slope_oct, both$slope_v8)))
cat(sprintf("  slope mean diff (V8 - Oct) : %+.4f\n",
            mean(both$slope_v8 - both$slope_oct)))
disc <- both[(slope_oct < 1) != (slope_v8 < 1)]
cat(sprintf("  Plots that flipped regime  : %d\n", nrow(disc)))
if (nrow(disc) > 0) print(disc[, .(id_plot, slope_oct, slope_v8)])

# Save
dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)
fwrite(DT_oct, "outputs/audit/regime_legacy_hourly_oct2025.csv")
fwrite(DT_v8,  "outputs/audit/regime_legacy_hourly_v8.csv")
cli_alert_success("Saved to outputs/audit/regime_legacy_hourly_{{oct2025,v8}}.csv")
