# ==============================================================================
# V10 — Sensitivity test : fCover baseline = mean(fCover) instead of 1
#
# Eva's question : "pourquoi avoir choisi une canopée complètement fermée comme
# baseline, et pas la moyenne comme les autres variables ?"
#
# Re-runs only the COALITIONS WHERE fCover=0 (i.e. the baseline value is used)
# from the V9 set, with fCover_b = mean(fCover) ≈ 0.869 instead of 1.0.
#
#   Archetypes : 4 × 8 affected coalitions = 32 sims
#   HOBO       : 4 affected bits in inverse path × 53 = 212 sims
#   TOTAL      : 244 sims, ~12 min with 4 workers
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

cli_h1("V10 — fCover_b = mean (sensitivity test)")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
FCOV_B_NEW <- mean(df_floor$fCover, na.rm = TRUE)
cli_alert("fCover baseline : OLD=1.000 → NEW={round(FCOV_B_NEW, 3)}")

fac_scs <- build_factorial_scenarios_archetypes(df_floor, fcov_b = FCOV_B_NEW)

# ----------------------------------------------------------------------------
# Phase 1 : Archetype coalitions where fCover=0 (8 of 16)
# ----------------------------------------------------------------------------
cli_h2("Phase 1 — Archetypes (8 affected coalitions × 4 = 32 sims)")

df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
LADMAP <- .ARCH_COAL_MAP
bits_fcov0 <- names(LADMAP)[as.integer(substr(LADMAP, 3, 3)) == 0]
cli_alert("Coalitions affected (fCover=0) : {length(bits_fcov0)}")

out_root_arch <- here::here("out_files/H1_archetypes_v10_fcovmean")
dir.create(out_root_arch, recursive = TRUE, showWarnings = FALSE)

jobs <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch <- df_archetypes[i, , drop = FALSE]
  arch_lbl <- sprintf("Arch_C%s", arch$Cluster)
  out_dir <- file.path(out_root_arch, arch_lbl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (sc_nm in bits_fcov0) {
    sc <- fac_scs[[sc_nm]]
    out_nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs) + 1L]] <- list(arch = arch, sc = sc, out_nc = out_nc,
                                      lbl = paste(arch_lbl, sc_nm))
  }
}
cli_alert("Archetype jobs : {length(jobs)}")

t0 <- Sys.time()
future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$arch, job$sc, job$out_nc,
                   CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("Phase 1 done in {round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1)} min")

# ----------------------------------------------------------------------------
# Phase 2 : HOBO sims for inverse path bits where fCover=0
#           (0000, 0001, 0101, 1101) × 53 = 212 sims
# ----------------------------------------------------------------------------
cli_h2("Phase 2 — HOBO (4 affected bits × 53 = 212 sims)")

bits_needed <- c("0000", "0001", "0101", "1101")
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
cli_alert("HOBO coalitions to rerun : {paste(bits_needed, collapse=', ')}")
cli_alert("Scenario names           : {paste(bit2scn[bits_needed], collapse=', ')}")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5

out_root_hobo <- here::here("out_files/musica_hobo_v10_fcovmean")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (bit in bits_needed) {
    sc_name <- bit2scn[bit]; if (is.na(sc_name)) next
    sc <- fac_scs[[sc_name]]
    out_dir <- file.path(out_root_hobo, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("HOBO jobs : {length(hobo_jobs)}")

t1 <- Sys.time()
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
cli_alert_success("Phase 2 done in {round(as.numeric(difftime(Sys.time(), t1, units='mins')), 1)} min")
cli_alert_success("V10 total in {round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1)} min")
