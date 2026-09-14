# Re-run the 2 corrupt sims for 41_39 (NULL and LVA_LAI) + the few sensors
# missing from new HOBO_0000 / HOBO_1000 / etc dirs (check). Output goes to
# H1f_0_Null_baseline and H1f_1_LAI_only (overwrites the 2 KB corrupt files).
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

cli_h1("Re-run corrupt HOBO sims for 41_39")
rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo_inputs <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                                    CFG$ids_to_remove))
df_hobo_inputs <- as.data.frame(df_hobo_inputs[id_plot == "41_39"])
cli_alert("Loaded 41_39 inputs : {nrow(df_hobo_inputs)} row")

df_sample_floor <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
lai_b   <- mean(df_sample_floor$LAI, na.rm = TRUE)
hmax_b  <- round(mean(floor(df_sample_floor$Hmax) + 1, na.rm = TRUE))
fcov_b  <- 1
sc_NULL <- list(name = "H1f_0_Null_baseline", bit_code = "0000",
                  lai_fn = function(pr) lai_b,
                  hmax_fn = function(pr) hmax_b,
                  fcover_fn = function(pr) fcov_b,
                  lad_fn = make_lad_uniform)
sc_LVA_LAI <- list(name = "H1f_1_LAI_only", bit_code = "1000",
                    lai_fn = function(pr) as.numeric(pr$LAI),
                    hmax_fn = function(pr) hmax_b,
                    fcover_fn = function(pr) fcov_b,
                    lad_fn = make_lad_uniform)
for (sc in list(sc_NULL, sc_LVA_LAI)) {
  out_nc <- here::here("out_files/musica_hobo_validation", sc$name,
                          sprintf("musica_out_HOBO_%s.nc", df_hobo_inputs$id_plot))
  cli_alert("Running {sc$name} for 41_39 (overwrite corrupt 2KB file)")
  run_musica_one(df_hobo_inputs, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
  cli_alert_success("Saved {out_nc}  size = {file.size(out_nc)}")
}
