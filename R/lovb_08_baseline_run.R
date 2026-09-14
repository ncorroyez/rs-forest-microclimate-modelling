# ==============================================================================
# LOVB analysis — Module 8 : Q10/Q90 baseline robustness (Option B)
#
# Tests whether the LOVB ranking is robust to the choice of baseline values.
# Defines TWO new baselines from cLHS sample :
#   - Q10 baseline : LAI_b = Q10(LAI), Hmax_b = Q10(Hmax), fCover_b = Q10(fCover),
#                    LAD_b = uniform   (canopée minimale)
#   - Q90 baseline : LAI_b = Q90(LAI), Hmax_b = Q90(Hmax), fCover_b = Q90(fCover),
#                    LAD_b = uniform   (canopée maximale)
# Reference mean baseline (already done) : LAI=3.13, Hmax=21, fCover=1, uniform.
#
# Subset : 50 cLHS plots stratified per cluster (~12-13 plots/cluster).
# Coalitions per baseline : NULL + 4 LOVB = 5 coalitions
# Total new MuSICA runs : 5 coalitions x 50 plots x 2 baselines = 500 (~2h).
#
# Output : out_files/H1_lovb_baselines/<baseline>_<bit>/musica_out_Sim_*.nc
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ncdf4)
  library(tidyverse)
  library(sf); library(terra)
  library(musica.tools); library(rmusica)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/forest.R"))
source(here::here("R/scenarios.R"))
source(here::here("R/musica.R"))

# Coalitions needed per baseline : NULL + 4 LOVB
.LOVB_BASELINE_BITS <- c("0000",
                          "0111",   # LOVB_LAI
                          "1011",   # LOVB_Hmax
                          "1101",   # LOVB_fCover
                          "1110")   # LOVB_LAD

# Build factorial scenarios with arbitrary (lai_b, hmax_b, fcov_b) baseline
build_lovb_baseline_scenarios <- function(label, lai_b, hmax_b, fcov_b) {
  cli_alert(sprintf("Baseline %s : LAI_b=%.3f  Hmax_b=%g  fCover_b=%.3f  LAD=uniform",
                     label, lai_b, hmax_b, fcov_b))
  mk_lai    <- function(real) if (real) function(pr) as.numeric(pr$LAI)    else function(pr) lai_b
  mk_hmax   <- function(real) if (real) function(pr) as.numeric(pr$Hmax)   else function(pr) hmax_b
  mk_fcover <- function(real) if (real) function(pr) as.numeric(pr$fCover) else function(pr) fcov_b
  mk_lad    <- function(real) if (real) make_lad_real                       else make_lad_uniform

  scs <- lapply(.LOVB_BASELINE_BITS, function(bit) {
    b <- strsplit(bit, "", fixed = TRUE)[[1]]
    list(
      name      = sprintf("BASE_%s_%s", label, bit),
      bit_code  = bit,
      baseline  = label,
      lai_fn    = mk_lai(b[1] == "1"),
      hmax_fn   = mk_hmax(b[2] == "1"),
      fcover_fn = mk_fcover(b[3] == "1"),
      lad_fn    = mk_lad(b[4] == "1")
    )
  })
  names(scs) <- vapply(scs, `[[`, "", "name")
  scs
}

# Subset stratified : keep first 12-13 plots per cluster from cLHS
select_baseline_subset <- function(df_sample, n_per_cluster = 13L) {
  setDT(df_sample)
  df_sample[, row_id := .I]
  set.seed(123)   # different from other seeds
  idx <- df_sample[, .SD[sample(.N, min(.N, n_per_cluster))], by = Cluster]
  df_sub <- as.data.frame(idx)
  cli_alert("Subset : {nrow(df_sub)} plots ({n_per_cluster}/cluster)")
  df_sub
}

# ---- Main runner -------------------------------------------------------------
lovb_run_baselines <- function(n_per_cluster = 13L) {
  cli_h1("Baseline robustness — launch MuSICA with Q10/Q90 baselines")

  df_sample <- readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample))
    df_sample <- df_sample %>% dplyr::rename(Cluster = Archetype)

  # Compute Q10/Q90 baselines
  lai_q10  <- quantile(df_sample$LAI,    0.10, type = 7, na.rm = TRUE)
  lai_q90  <- quantile(df_sample$LAI,    0.90, type = 7, na.rm = TRUE)
  hmax_q10 <- round(quantile(df_sample$Hmax,   0.10, type = 7, na.rm = TRUE))
  hmax_q90 <- round(quantile(df_sample$Hmax,   0.90, type = 7, na.rm = TRUE))
  fcov_q10 <- as.numeric(quantile(df_sample$fCover, 0.10, type = 7, na.rm = TRUE))
  fcov_q90 <- as.numeric(quantile(df_sample$fCover, 0.90, type = 7, na.rm = TRUE))

  cli_alert(sprintf("Q10 : LAI=%.2f Hmax=%g fCover=%.2f", lai_q10, hmax_q10, fcov_q10))
  cli_alert(sprintf("Q90 : LAI=%.2f Hmax=%g fCover=%.2f", lai_q90, hmax_q90, fcov_q90))

  scs_q10 <- build_lovb_baseline_scenarios("Q10", lai_q10, hmax_q10, fcov_q10)
  scs_q90 <- build_lovb_baseline_scenarios("Q90", lai_q90, hmax_q90, fcov_q90)
  scs_all <- c(scs_q10, scs_q90)

  df_sub <- select_baseline_subset(df_sample, n_per_cluster = n_per_cluster)

  # Run MuSICA for each scenario at each plot in subset
  out_root <- here::here("out_files/H1_lovb_baselines")
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(scs_all)) {
    sc <- scs_all[[nm]]
    out_dir <- file.path(out_root, sc$name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    cli_alert("[{sc$name}] {nrow(df_sub)} plots ->  {basename(out_dir)}")
    for (i in seq_len(nrow(df_sub))) {
      pr <- df_sub[i, , drop = FALSE]
      x_i <- as.integer(round(pr$x)); y_i <- as.integer(round(pr$y))
      out_nc <- file.path(out_dir, sprintf("musica_out_Sim_%03d_X%d_Y%d.nc", i, x_i, y_i))
      if (file.exists(out_nc)) next
      run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
    }
  }
  cli_alert_success("Baseline runs done. Folders : {.path {out_root}/BASE_*}")
  invisible(out_root)
}

if (sys.nframe() == 0) {
  lovb_run_baselines(n_per_cluster = 13L)
}
