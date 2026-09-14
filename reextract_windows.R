# ==============================================================================
# Multi-window re-extraction of HOBO validation metrics (NO re-simulation).
#
# The MuSICA .nc outputs are full-year (2021-01 .. 2022-08). Only the metric
# *extraction* is windowed via date_seq. This script re-extracts metrics for all
# already-simulated scenarios over several validation windows to test whether
# dynamic (S2-driven) phenology beats static (parametric) phenology once the
# window includes shoulder seasons (leaf-out / early senescence), where the
# canopy actually changes — vs summer, where it is flat and dynamics add only
# noise.
#
# Run from the z_Example repo root:  Rscript reextract_windows.R
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})

src <- list.files("R", pattern = "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]      # SAFE: skip auto-executing scripts
invisible(lapply(src, source))
source("Chapter3_config.R")

# ---- Load prep + scenarios (same as Chapter3_main.R) -------------------------
prep         <- load_lai_prep(CFG_C3)
df_hobo      <- prep$df_plots
ts_list      <- prep$ts_by_plot
scenarios_c3 <- make_all_scenarios_c3(ts_list,
                                      list_year = CFG_C3$list_year,
                                      d_opt     = CFG_C3$d_opt_m,
                                      mode      = CFG_C3$scenarios_mode)

cat(sprintf("Loaded %d plots, %d scenarios\n",
            nrow(df_hobo), length(scenarios_c3)))

# SAFETY: only keep scenarios whose nc dir already exists and is non-empty, so
# re-extraction can never trigger a full (heavy) simulation of a missing dir.
nc_parent <- file.path(CFG_C3$out_dir, "nc")
scenarios_c3 <- Filter(function(sc) {
  d <- file.path(nc_parent, sc$name)
  dir.exists(d) && length(list.files(d, pattern = "\\.nc$")) > 0
}, scenarios_c3)
cat(sprintf("Re-extracting %d scenarios with existing nc: %s\n",
            length(scenarios_c3),
            paste(vapply(scenarios_c3, function(s) s$name, ""), collapse = ", ")))

# ---- Windows to test ---------------------------------------------------------
# S2 time series spans 2021-04-01 .. 2021-10-31 (ts_date_range), so DYN phenology
# is only well-defined inside that range. Nov senescence would need ts rebuild +
# re-sim (phase 2), so the autumn window is capped at Oct 31.
mkseq <- function(a, b) seq(as.Date(a), as.Date(b), by = "day")
windows <- list(
  summer    = mkseq("2021-06-01", "2021-09-30"),  # baseline (flat canopy)
  leafout   = mkseq("2021-04-01", "2021-05-31"),  # leaf-out: DYN should help most
  autumn    = mkseq("2021-10-01", "2021-10-31"),  # early senescence (S2 ts cap)
  shoulders = c(mkseq("2021-04-01", "2021-05-31"),
                mkseq("2021-10-01", "2021-10-31")),# leaf-out + early senescence
  extended  = mkseq("2021-04-01", "2021-10-31"),  # full S2-ts span
  # full year + leaf-off winter. NOTE: outside Apr-Oct the DYN (S2) phenology is
  # EXTRAPOLATED (no S2 obs) so winter LAI may be a spurious residual, whereas
  # STATIC parametric phenology models leaf-off (LAI->0). Interpret with care.
  winter    = c(mkseq("2021-01-01", "2021-03-31"),
                mkseq("2021-11-01", "2021-12-31")),
  fullyear  = mkseq("2021-01-01", "2021-12-31")
)

# ---- Re-extract per window ---------------------------------------------------
all_metrics <- list()
for (wn in names(windows)) {
  ds <- windows[[wn]]
  cat(sprintf("\n===== window %-10s (%s .. %s, %d days) =====\n",
              wn, min(ds), max(ds), length(ds)))
  df_macro      <- extract_macro_daily(CFG_C3$forcing_file, ds)
  df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds,
                                   df_macro, CFG_C3$ids_to_remove)
  val_out <- validate_scenarios_at_hobos(
    df_hobo_inputs = df_hobo,
    df_hobo_daily  = df_hobo_daily,
    scenarios      = scenarios_c3,
    parent_dir     = file.path(CFG_C3$out_dir, "nc"),
    df_macro       = df_macro,
    date_seq       = ds,
    forcing_file   = CFG_C3$forcing_file,
    musica_cmd     = CFG_C3$musica_cmd,
    force          = FALSE       # reuse existing nc; never re-simulate
  )
  m <- val_out$metrics
  m$window <- wn
  all_metrics[[wn]] <- m
}

res <- bind_rows(all_metrics) %>%
  dplyr::select(window, scenario, n, r2, rmse, mae, bias) %>%
  arrange(window, rmse)

out_csv <- file.path(CFG_C3$out_dir, "tables", "c3_metrics_multiwindow.csv")
write.csv(res, out_csv, row.names = FALSE)

cat("\n\n################ MULTI-WINDOW METRICS ################\n")
for (wn in names(windows)) {
  cat(sprintf("\n--- %s ---\n", wn))
  sub <- res[res$window == wn, c("scenario", "n", "r2", "rmse", "bias")]
  print(as.data.frame(sub), row.names = FALSE, digits = 3)
}
cat(sprintf("\nWritten: %s\n", out_csv))
