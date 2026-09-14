# ==============================================================================
# PIPELINE STAGE 03 — MuSICA runs (HOBO × 16 coalitions, mean fCover baseline)
# Self-contained single-source set in PIPE$MUSICA_DIR_HOBO. Idempotent: only
# missing coalition×plot NCs are run (run_musica_one skips existing). Legacy binary.
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
suppressMessages({ library(terra); library(parallel); library(rmusica); library(musica.tools) })

cli_h1("STAGE 03 — MuSICA HOBO × 16 (fCover baseline = mean)")
stopifnot(file.exists(CFG$musica_cmd))

bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
df_floor <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))
fac_scs  <- build_factorial_scenarios_archetypes(df_floor, fcov_b = PIPE$FCOV_BASELINE)

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5            # canonical fCover floor

jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (bit in PIPE$BITS) {
    out_dir <- file.path(PIPE$MUSICA_DIR_HOBO, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs)+1L]] <- list(pr = pr, sc = fac_scs[[bit2scn[bit]]],
                                    out_nc = out_nc, lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("Missing sims to run : {length(jobs)} (of {length(PIPE$BITS)*nrow(df_hobo)})")

if (length(jobs) > 0) {
  t0 <- Sys.time()
  mclapply(jobs, function(job) tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc, CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR: %s\n", job$lbl, e$message))),
    mc.cores = 4, mc.preschedule = FALSE)
  cli_alert_success("Ran {length(jobs)} sims in {round(as.numeric(difftime(Sys.time(),t0,units='mins')),1)} min")
}

# completeness check
cnt <- sapply(PIPE$BITS, function(b) length(list.files(file.path(PIPE$MUSICA_DIR_HOBO, b), "\\.nc$")))
print(data.frame(bit = PIPE$BITS, n = cnt))
pipe_assert(all(cnt >= nrow(df_hobo)),
            sprintf("HOBO set complete: 16 coalitions × %d plots in %s", nrow(df_hobo), basename(PIPE$MUSICA_DIR_HOBO)))
