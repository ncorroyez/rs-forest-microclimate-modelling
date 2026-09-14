# ==============================================================================
# PIPELINE STAGE 03z (BRANCH z05) — MuSICA HOBO × 16 coalitions using the new
# UNDERSTORY-AWARE z0.5 LAD profiles (in_files/lad_z05/Blois_lad_z05_r25.csv).
# Same scenario design (fCover baseline = mean), legacy binary, idempotent.
# Per-plot LAD/LAI/Hmax come from z0.5; fCover from the raster (floored 0.5).
# Output: out_files/musica_hobo_z05/<bit>/
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(parallel)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/musica.R")); source(here::here("R/lad.R"))
source(here::here("R/validation.R")); source(here::here("R/h1_shapley_archetypes.R"))

MUSICA_BIN <- Sys.getenv("MUSICA_BIN", CFG$musica_cmd)   # default legacy v3.2.0; set to model-3.2.3 to compare
stopifnot(file.exists(MUSICA_BIN))
OUT_SUF <- Sys.getenv("OUT_SUF", Sys.getenv("LAD_SRC","z05"))   # output dir suffix (≠ input when comparing binaries)
cli_h1(sprintf("STAGE 03z — MuSICA HOBO × 16 | bin=%s | out=%s", basename(dirname(MUSICA_BIN)), OUT_SUF))
OUT_ROOT <- here::here(sprintf("out_files/musica_hobo_%s", OUT_SUF)); dir.create(OUT_ROOT, recursive=TRUE, showWarnings=FALSE)
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

# ---- build z0.5 HOBO inputs : LAD/LAI/Hmax from z05 csv, fCover from raster ---
SRC <- Sys.getenv("LAD_SRC","z05"); z05 <- fread(here::here(sprintf("in_files/lad_z05/Blois_lad_%s_r25.csv",SRC)))
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
hob_r <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
hp <- sf::st_read(CFG$hobo_geojson, quiet=TRUE) %>% filter(!id_plot %in% CFG$ids_to_remove)
hob_r$id_plot <- hp$id_plot
fcov <- hob_r[, .(id_plot, fCover, VCI = if ("VCI" %in% names(hob_r)) VCI else NA_real_)]

df <- merge(z05, fcov, by="id_plot")           # LAD_Layer_* (0.75..), LAI, Hmax, x, y, fCover
df[fCover < 0.5, fCover := 0.5]                 # canonical floor
# LAI is kept ONE-SIDED here. The single one-sided -> total two-sided doubling
# (× 2) happens in exactly ONE place: run_musica_one() (R/musica.R), applied
# identically to the real per-plot LAI AND to the cLHS baselines, so the
# real-vs-baseline toggle stays on the same scale. Do NOT double here (that was
# the historical double-doubling bug: real LAI reached MuSICA at ~17, baselines
# at ~13.5). LAD column scale is irrelevant — .lad_rescale() renormalizes the
# profile to the (doubled) LAI inside lad_fn, so only the LAD shape is used.
cli_alert("z0.5 LAI kept one-sided (median {round(median(df$LAI),2)}); ×2 to two-sided applied once in run_musica_one")
df <- as.data.frame(df)
cli_alert("z0.5 HOBO inputs: {nrow(df)} plots ; LAD cols {sum(grepl('LAD_Layer_',names(df)))} ; LAI med {round(median(df$LAI),2)}")

# scenarios: BASELINES from the SAME cLHS sample as legacy (lai_b~3.13) for a
# clean attribution comparison — only the per-plot REAL values come from z0.5.
df_clhs <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs <- build_factorial_scenarios_archetypes(as.data.frame(df_clhs), fcov_b = "mean")

jobs <- list()
for (i in seq_len(nrow(df))) {
  pr <- df[i, , drop=FALSE]
  for (bit in unname(.ARCH_COAL_MAP)) {
    out_dir <- file.path(OUT_ROOT, bit); dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs)+1L]] <- list(pr=pr, sc=fac_scs[[bit2scn[bit]]], out_nc=out_nc, lbl=paste(bit,pr$id_plot))
  }
}
cli_alert("Missing sims: {length(jobs)} / {16*nrow(df)}")
if (length(jobs) > 0) {
  t0 <- Sys.time()
  mclapply(jobs, function(job) tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc, CFG$forcing_file, MUSICA_BIN),
    error=function(e) cat(sprintf("[%s] ERR: %s\n", job$lbl, e$message))),
    mc.cores=4, mc.preschedule=FALSE)
  cli_alert_success("Ran {length(jobs)} sims in {round(as.numeric(difftime(Sys.time(),t0,units='mins')),1)} min")
}
cnt <- sapply(unname(.ARCH_COAL_MAP), function(b) length(list.files(file.path(OUT_ROOT,b),"\\.nc$")))
cli_alert("Coalition completeness: {min(cnt)}–{max(cnt)} / {nrow(df)}")
