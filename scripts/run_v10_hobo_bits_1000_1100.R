# ==============================================================================
# Complete the mean-fCover-baseline (v10) HOBO factorial: run the ONLY two
# coalitions missing for an exact 16-coalition Shapley → 1000, 1100.
# (Both are fCover=0 i.e. baseline-dependent, so they MUST use fcov_b=mean.)
# 2 × 53 = 106 sims, legacy binary, mclapply (fork). Output: musica_hobo_v10_fcovmean/<bit>/
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(parallel)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/lad.R"))
source(here::here("R/validation.R")); source(here::here("R/h1_shapley_archetypes.R"))

stopifnot(file.exists(CFG$musica_cmd))
cli_h1("v10 fcovmean — complete bits 1000 & 1100 (legacy binary)")

bits_needed <- c("1000", "1100")
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs <- build_factorial_scenarios_archetypes(df_floor, fcov_b = "mean")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5

out_root <- here::here("out_files/musica_hobo_v10_fcovmean")
jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (bit in bits_needed) {
    sc <- fac_scs[[bit2scn[bit]]]
    out_dir <- file.path(out_root, bit); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs)+1L]] <- list(pr = pr, sc = sc, out_nc = out_nc, lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("Jobs: {length(jobs)}")
if (length(jobs) == 0) { cli_alert_success("Nothing to do."); quit(save="no") }

t0 <- Sys.time()
mclapply(jobs, function(job) tryCatch(
  run_musica_one(job$pr, job$sc, job$out_nc, CFG$forcing_file, CFG$musica_cmd),
  error = function(e) cat(sprintf("[%s] ERROR: %s\n", job$lbl, e$message))),
  mc.cores = 4, mc.preschedule = FALSE)
cli_alert_success("Done in {round(as.numeric(difftime(Sys.time(), t0, units='mins')),1)} min")
for (bit in bits_needed) cli_alert("{bit}: {length(list.files(file.path(out_root,bit),'\\\\.nc$'))} NC")
