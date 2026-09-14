# ==============================================================================
# test_forcing_height_30.R
#
# Quick test : 10 HOBO plots × 2 conditions
#   - fh_dyn  : forcing_height = hmax + 2  (comportement actuel)
#   - fh_30   : forcing_height = 30        (constant, hauteur réelle ERA5)
#
# Question : si le spread inter-plot de Tmax(nair=1) s'effondre en fh_30,
# alors la variance inter-plot que voit le Chapitre 1 est principalement
# un artefact d'ajustement de hauteur, pas un effet canopée.
#
# Output  : out_files/forcing_height_test/{fh_dyn,fh_30}/REF_all_real/*.nc
#           run_forcing_height_test.log
# ==============================================================================

setwd("/home/corroyez/Documents/z_Example_rmusica_31012025")

suppressPackageStartupMessages({
  library(musica.tools)
  library(rmusica)
  library(ncdf4)
  library(tidyverse)
  library(sf); library(terra)
})

source("R/config.R")
source("R/io.R")
source("R/forest.R")
source("R/lad.R")
source("R/scenarios.R")
source("R/musica.R")
source("R/validation.R")   # build_hobo_inputs

# ---- Setup partagé ----------------------------------------------------------
.BLOIS_SETUP <- list(
  "site_latitude"     = 47.57,
  "site_longitude"    = 1.26,
  "site_altitude"     = 39.18,
  "forcing_timestep"  = 3600,
  "n_leaf_age"        = 1,
  "param_files_path"  = '"./in_files/Blois/in_files/"',
  "species_names"     = '"musica_veg1"',
  "soil_file_name"    = '"./in_files/Blois/in_files/musica_soil.nml"',
  "history_variables" = c("layer_thickness", "z_soil", "dz_soil", "t_soil",
                          "w_soil", "t_air", "w_air", "wind"),
  "abl_flag"          = '"none"'
)

MUSICA_BIN <- normalizePath("in_files/model-3.2.3/musica")
MUSICA_NML <- normalizePath("in_files/model-3.2.3/musica.nml")
MUSICA_VAR <- normalizePath("in_files/model-3.2.3/variables.csv")

CONDS <- list(
  list(label = "fh_dyn", extra_setup = .BLOIS_SETUP),                          # legacy : hmax+2
  list(label = "fh_30",  extra_setup = c(.BLOIS_SETUP, list("forcing_height" = 30)))
)

OUT_ROOT <- "out_files/forcing_height_test"

# ---- Inputs ------------------------------------------------------------------
cat("[1/3] Loading inputs...\n")
rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                     CFG$ids_to_remove)
df_sample      <- readRDS("out_files/Sensitivity_Analysis/clhs_sample.rds")
ref_sc         <- scenario_reference(df_sample)

# 10 plots — pris à intervalles réguliers pour couvrir l'éventail de Hmax
N_TEST <- 10
idx    <- round(seq(1, nrow(df_hobo_inputs), length.out = N_TEST))
df_test <- df_hobo_inputs[idx, ]
cat(sprintf("  Testing %d plots: %s\n", N_TEST,
            paste(df_test$id_plot, collapse = ", ")))

# ---- Run ---------------------------------------------------------------------
cat("[2/3] Running 10 plots x 2 conditions...\n")
run_one_cond <- function(cond) {
  out_dir <- file.path(OUT_ROOT, cond$label, ref_sc$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(nrow(df_test))) {
    pr <- df_test[i, ]
    of <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(of)) next
    cat(sprintf("  [%s] HOBO %s (Hmax=%.1f, LAI=%.2f, fcover=%.2f)\n",
                cond$label, pr$id_plot, pr$Hmax, pr$LAI, pr$fcover))
    run_musica_one(pr, ref_sc, of, CFG$forcing_file, MUSICA_BIN,
                   musica_nml       = MUSICA_NML,
                   musica_variables = MUSICA_VAR,
                   extra_setup      = cond$extra_setup)
  }
  out_dir
}
out_dirs <- sapply(CONDS, run_one_cond)
names(out_dirs) <- sapply(CONDS, `[[`, "label")

# ---- Analyse -----------------------------------------------------------------
cat("[3/3] Extracting Tmax(nair=1) per plot, per condition...\n")
extract <- function(nc_file) {
  if (!file.exists(nc_file) || file.size(nc_file) < 1e6)
    return(c(hmax = NA, tmax_low = NA, tmax_high = NA))
  nc   <- nc_open(nc_file)
  tair <- ncvar_get(nc, "Tair_z")
  hmax <- ncvar_get(nc, "veget_height_top")[1]
  nc_close(nc)
  d <- dim(tair)
  c(hmax     = hmax,
    tmax_low = max(tair[1, ],   na.rm = TRUE) - 273.15,
    tmax_high= max(tair[d[1], ], na.rm = TRUE) - 273.15)
}

results <- map_dfr(CONDS, function(cond) {
  out_dir <- file.path(OUT_ROOT, cond$label, ref_sc$name)
  map_dfr(df_test$id_plot, function(id) {
    nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", id))
    r  <- extract(nc)
    tibble(condition = cond$label, id_plot = id,
           hmax = r["hmax"], tmax_low = r["tmax_low"], tmax_high = r["tmax_high"])
  })
})

cat("\n=== Resultats par plot ===\n")
print(results %>% pivot_wider(names_from = condition,
                              values_from = c(tmax_low, tmax_high),
                              id_cols = c(id_plot, hmax)),
      n = N_TEST)

cat("\n=== Spread inter-plot par condition ===\n")
summary <- results %>%
  group_by(condition) %>%
  summarise(
    n             = sum(!is.na(tmax_low)),
    mean_tmax_low = mean(tmax_low, na.rm = TRUE),
    sd_tmax_low   = sd(tmax_low,   na.rm = TRUE),
    range_low     = max(tmax_low, na.rm = TRUE) - min(tmax_low, na.rm = TRUE),
    sd_tmax_high  = sd(tmax_high,  na.rm = TRUE),
    .groups = "drop"
  )
print(summary)

cat("\n=== Interpretation ===\n")
ratio <- summary$sd_tmax_low[summary$condition == "fh_30"] /
         summary$sd_tmax_low[summary$condition == "fh_dyn"]
cat(sprintf("  Ratio SD(Tmax_low) fh_30 / fh_dyn = %.2f\n", ratio))
if (!is.na(ratio) && ratio < 0.3)
  cat("  -> La variance inter-plot s'effondre en fh=30 :\n",
      "     la variation observée en fh_dyn vient majoritairement de l'ajustement\n",
      "     de hauteur de reference (hmax+2), PAS de la canopée.\n")
if (!is.na(ratio) && ratio >= 0.7)
  cat("  -> La variance reste : la canopée joue vraiment.\n")
if (!is.na(ratio) && ratio >= 0.3 && ratio < 0.7)
  cat("  -> Contribution mixte : la canopée joue partiellement.\n")

cat("\n[done]\n")
