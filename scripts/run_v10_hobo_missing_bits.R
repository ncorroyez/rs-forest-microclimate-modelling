# ==============================================================================
# Run the 5 missing HOBO coalitions to enable TRUE greedy forward selection on
# HOBO data (answers JB Féret's question about feature selection).
#
# Missing bits : 0011, 0100, 0110, 1001, 1011
#   - bits where fCover=0 (0100, 1001) → V10 baseline matters
#   - bits where fCover=1 (0011, 0110, 1011) → baseline-independent
# All run with V10 fcov_b=mean for consistency.
# 5 × 53 = 265 sims, ~12 min with 4 workers.
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

cli_h1("V10 HOBO factorial completion — 5 missing bits")

bits_needed <- c("0011", "0100", "0110", "1001", "1011")
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
cli_alert("Bits to add : {paste(bits_needed, collapse=', ')}")
cli_alert("Scenarios   : {paste(bit2scn[bits_needed], collapse=', ')}")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs <- build_factorial_scenarios_archetypes(df_floor, fcov_b = "mean")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5

out_root <- here::here("out_files/musica_hobo_v10_fcovmean")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (bit in bits_needed) {
    sc_name <- bit2scn[bit]; if (is.na(sc_name)) next
    sc <- fac_scs[[sc_name]]
    out_dir <- file.path(out_root, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("HOBO jobs : {length(hobo_jobs)}")

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
cli_alert_success("V10 HOBO completion done in {round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1)} min")
