# ==============================================================================
# d20 sensitivity: re-simulate the d_opt-dependent scenarios with d_opt = 20 m
# (config must already be set to dopt_variant="d20", d_opt_m=20 and the d20 prep
# must be on disk). Sims are written to a SEPARATE parent dir (nc_d20/) so the
# d7 ("common") baseline nc are never clobbered. Then metrics are re-extracted
# over the same 5 windows as reextract_windows.R for a like-for-like comparison.
#
# env vars:
#   D20_NPLOTS : limit number of HOBO plots (timing dry-run); empty = all
# ==============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(rmusica); library(musica.tools)
})
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
source("Chapter3_config.R")

stopifnot(CFG_C3$dopt_variant == "d20", CFG_C3$d_opt_m == 20)

prep         <- load_lai_prep(CFG_C3)
df_hobo      <- prep$df_plots
ts_list      <- prep$ts_by_plot
scenarios_c3 <- make_all_scenarios_c3(ts_list,
                                      list_year = CFG_C3$list_year,
                                      d_opt     = CFG_C3$d_opt_m,
                                      mode      = CFG_C3$scenarios_mode)

# d_opt-dependent scenarios only (DOPT + LADOPT both contain "DOPT")
dep <- Filter(function(sc) grepl("DOPT", sc$name), scenarios_c3)
cat(sprintf("d20: anchor LAI_ALS_DOPT mean=%.2f (d7 was 1.56, full 4.04)\n",
            mean(df_hobo$LAI_ALS_DOPT, na.rm = TRUE)))
cat(sprintf("d_opt-dependent scenarios (%d): %s\n",
            length(dep), paste(vapply(dep, function(s) s$name, ""), collapse = ", ")))

# optional plot-limit for timing dry-run
npl <- Sys.getenv("D20_NPLOTS", "")
df_run <- if (nzchar(npl)) df_hobo[seq_len(min(as.integer(npl), nrow(df_hobo))), ] else df_hobo
cat(sprintf("Simulating %d plots into nc_d20/\n", nrow(df_run)))

parent_d20 <- file.path(CFG_C3$out_dir, "nc_d20")

windows <- local({
  mkseq <- function(a, b) seq(as.Date(a), as.Date(b), by = "day")
  list(
    summer    = mkseq("2021-06-01", "2021-09-30"),
    leafout   = mkseq("2021-04-01", "2021-05-31"),
    autumn    = mkseq("2021-10-01", "2021-10-31"),
    shoulders = c(mkseq("2021-04-01", "2021-05-31"), mkseq("2021-10-01", "2021-10-31")),
    extended  = mkseq("2021-04-01", "2021-10-31")
  )
})

t0 <- Sys.time()
all_metrics <- list()
for (wn in names(windows)) {
  ds <- windows[[wn]]
  df_macro      <- extract_macro_daily(CFG_C3$forcing_file, ds)
  df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, ds, df_macro, CFG_C3$ids_to_remove)
  val_out <- validate_scenarios_at_hobos(
    df_hobo_inputs = df_run,
    df_hobo_daily  = df_hobo_daily,
    scenarios      = dep,
    parent_dir     = parent_d20,
    df_macro       = df_macro,
    date_seq       = ds,
    forcing_file   = CFG_C3$forcing_file,
    musica_cmd     = CFG_C3$musica_cmd,
    force          = FALSE
  )
  m <- val_out$metrics; m$window <- wn; all_metrics[[wn]] <- m
  cat(sprintf("[window %s] done (%.1f min elapsed)\n", wn,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- bind_rows(all_metrics) %>%
  mutate(scenario = paste0(scenario, "@d20")) %>%
  dplyr::select(window, scenario, n, r2, rmse, mae, bias) %>%
  arrange(window, rmse)
out_csv <- file.path(CFG_C3$out_dir, "tables", "c3_metrics_d20.csv")
write.csv(res, out_csv, row.names = FALSE)

cat("\n\n################ d20 METRICS ################\n")
for (wn in names(windows)) {
  cat(sprintf("\n--- %s ---\n", wn))
  sub <- res[res$window == wn, c("scenario", "n", "r2", "rmse", "bias")]
  print(as.data.frame(sub), row.names = FALSE, digits = 3)
}
cat(sprintf("\nWritten: %s\n", out_csv))
