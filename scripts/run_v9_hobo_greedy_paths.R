# ==============================================================================
# Extend V9 HOBO sims to the bit coalitions needed for the GREEDY forward path
# (responding to JB Féret's request for full forward feature selection).
#
# Missing coalitions vs the existing slope-LOO-ordered path :
#   - 0010 (+fCover only)         : greedy-slope step 1
#   - 1010 (+LAI+fCover)          : greedy-slope/Tmax step 2-3
#   - 1110 (+LAI+fCover+Hmax)     : greedy step 3-4 (also Hmax_fCover_Lad bit reversed)
# Total : 3 coalitions × 53 HOBO = 159 sims, ~8 min with 4 workers.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(future); library(furrr)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

cli_h1("V9 HOBO greedy-path completion ({CFG$musica_cmd})")

bits_needed <- c("0010", "1010", "1110")
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
cli_alert("Coalitions to add : {paste(bits_needed, collapse=', ')}")
cli_alert("Mapped scenario names : {paste(bit2scn[bits_needed], collapse=', ')}")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs_hobo <- build_factorial_scenarios_archetypes(df_floor)

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
cli_alert("HOBO sensors : {nrow(df_hobo)}")

out_root_hobo <- here::here("out_files/musica_hobo_v9")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (bit in bits_needed) {
    sc_name <- bit2scn[bit]; if (is.na(sc_name)) next
    sc <- fac_scs_hobo[[sc_name]]
    out_dir <- file.path(out_root_hobo, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("HOBO jobs to run : {length(hobo_jobs)}")

t0 <- Sys.time()
future_walk(hobo_jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc,
                   CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("Greedy HOBO sims done in {round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1)} min")
