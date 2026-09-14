# ==============================================================================
# Phase 2 — Re-run 468 MuSICA simulations with floored fCover.
#
# Outputs:
#   out_files/H1_factorial_floor05/H1F_*  (54 plots × 8 coalitions = 432 NC)
#   out_files/H1_archetypes_floor05/Arch_C1/* (1 archetype × 16 coalitions = 16 NC)
#   out_files/musica_hobo_validation_floor05/HOBO_*/* (4 sensors × 5 coalitions = 20 NC)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(ncdf4); library(tidyverse); library(sf); library(terra)
  library(musica.tools); library(rmusica)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/forest.R"))
source(here::here("R/scenarios.R"))
source(here::here("R/musica.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

cli_h1("Phase 2 — Re-run 468 MuSICA sims with floored fCover")
t_start <- Sys.time()

# ---- Load floored inputs -----------------------------------------------------
df_floor <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
df_orig  <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds")))
if ("Archetype" %in% names(df_orig) && !"Cluster" %in% names(df_orig))
  df_orig[, Cluster := Archetype]
affected_clhs <- readRDS(here::here("out_files/Sensitivity_Analysis/floor05_affected_clhs.rds"))

# Subset df_floor to the 54 affected plots (preserve original df_floor row order)
plot_keys <- affected_clhs$plots[, paste(x, y, sep = "_")]
df_floor_aff <- df_floor[paste(x, y, sep = "_") %in% plot_keys]
cli_alert("cLHS subset : {nrow(df_floor_aff)} affected plots")

# ---- Part A : cLHS factorial (54 plots × 8 coalitions) -----------------------
cli_h2("Part A — cLHS factorial (8 coalitions where fCover bit=1)")
scenarios_all <- scenarios_h1_factorial(df_floor)  # uses df_floor (with floored fCover) as data source
out_root_clhs <- here::here("out_files/H1_factorial_floor05")
dir.create(out_root_clhs, recursive = TRUE, showWarnings = FALSE)

# Convert to data.frame so make_lad_real can do plot_row[lad_cols] (data.table fails)
df_floor_aff_df <- as.data.frame(df_floor_aff)
for (sc_name in affected_clhs$scenarios) {
  sc <- scenarios_all[[sc_name]]
  if (is.null(sc)) { cli_alert_warning("Scenario not found : {sc_name}"); next }
  out_dir <- file.path(out_root_clhs, sc_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cli_alert("[{sc_name}] : {nrow(df_floor_aff_df)} plots")
  for (i in seq_len(nrow(df_floor_aff_df))) {
    pr <- df_floor_aff_df[i, , drop = FALSE]
    sim_id <- sprintf("Sim_%03d_X%d_Y%d", i, round(pr$x), round(pr$y))
    out_nc <- file.path(out_dir, paste0("musica_out_", sim_id, ".nc"))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
  }
}
cli_alert_success("cLHS factorial done at t = {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")

# ---- Part B : C1 archetype (16 coalitions, centroid recomputed from df_floor)
cli_h2("Part B — C1 archetype (16 coalitions, floored centroid)")
df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
# We only need to re-run C1 (C2/3/4 centroids unchanged because no plot in them was floored)
df_arch_c1 <- as.data.frame(df_archetypes[df_archetypes$Cluster == "1", , drop = FALSE])
fac_scs_arch <- build_factorial_scenarios_archetypes(df_floor)
out_root_arch <- here::here("out_files/H1_archetypes_floor05")
dir.create(out_root_arch, recursive = TRUE, showWarnings = FALSE)

for (sc in fac_scs_arch) {
  out_dir <- file.path(out_root_arch, "Arch_C1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_nc <- file.path(out_dir,
                        sprintf("musica_out_Arch_C1_%s.nc", sc$name))
  if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
  run_musica_one(df_arch_c1, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
}
cli_alert_success("C1 archetype done at t = {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")

# ---- Part C : HOBO sensors (4 sensors × 5 coalitions) ------------------------
cli_h2("Part C — HOBO sensors (4 with fCover < 0.5)")
rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo_inputs <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                                    CFG$ids_to_remove))
# Floor fCover in HOBO inputs
df_hobo_floor <- copy(df_hobo_inputs)
df_hobo_floor[fCover < 0.5, fCover := 0.5]
affected_hobo <- readRDS(here::here("outputs/lovb/data/floor05_affected_hobo.rds"))
df_hobo_aff <- df_hobo_floor[id_plot %in% affected_hobo$sensors$id_plot]
cli_alert("HOBO subset : {nrow(df_hobo_aff)} sensors")

# Build the same scenarios used in lovb_06_hobo_run.R
build_hobo_scenarios <- function(df_sample) {
  lai_b   <- mean(df_sample$LAI, na.rm = TRUE)
  hmax_b  <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  fcov_b  <- 1
  mk_lai    <- function(real) if (real) function(pr) as.numeric(pr$LAI)    else function(pr) lai_b
  mk_hmax   <- function(real) if (real) function(pr) as.numeric(pr$Hmax)   else function(pr) hmax_b
  mk_fcover <- function(real) if (real) function(pr) as.numeric(pr$fCover) else function(pr) fcov_b
  mk_lad    <- function(real) if (real) make_lad_real                       else make_lad_uniform
  needed <- c("0010", "0111", "1011", "1110", "1111")   # fCover bit=1 in HOBO 10-set
  lapply(needed, function(bit) {
    b <- strsplit(bit, "", fixed = TRUE)[[1]]
    list(name = paste0("HOBO_", bit), bit_code = bit,
          lai_fn = mk_lai(b[1]=="1"), hmax_fn = mk_hmax(b[2]=="1"),
          fcover_fn = mk_fcover(b[3]=="1"), lad_fn = mk_lad(b[4]=="1"))
  }) -> scs
  names(scs) <- vapply(scs, `[[`, "", "name")
  scs
}
hobo_scenarios <- build_hobo_scenarios(df_floor)
out_root_hobo <- here::here("out_files/musica_hobo_validation_floor05")
dir.create(out_root_hobo, recursive = TRUE, showWarnings = FALSE)
for (sc in hobo_scenarios) {
  out_dir <- file.path(out_root_hobo, sc$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  df_hobo_aff_df <- as.data.frame(df_hobo_aff)
  cli_alert("[{sc$name}] : {nrow(df_hobo_aff_df)} sensors")
  for (i in seq_len(nrow(df_hobo_aff_df))) {
    pr <- df_hobo_aff_df[i, , drop = FALSE]
    out_nc <- file.path(out_dir,
                          sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
  }
}
cli_alert_success("HOBO done at t = {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")
cli_alert_success("Phase 2 TOTAL : {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")
