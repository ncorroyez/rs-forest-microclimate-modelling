# Re-run 21 missing HOBO legacy sims (NULL + LVA_LAI for 10-12 sensors, skipping 41_39).
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

cli_h1("Re-run missing HOBO NULL + LVA_LAI sims")
t_start <- Sys.time()

# Sensors needing each coalition (excluding 41_39 which crashes)
need_null <- c("41_04","41_09","41_10","41_11","41_17","41_26","41_40","41_54","41_59")  # 9 (41_39 excluded)
need_lvalai <- c("41_04","41_09","41_10","41_11","41_17","41_26","41_36","41_40","41_42","41_54","41_59")  # 11 (41_39 excluded)

rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo_inputs <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                                    CFG$ids_to_remove))
df_sample_floor <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
lai_b   <- mean(df_sample_floor$LAI, na.rm = TRUE)
hmax_b  <- round(mean(floor(df_sample_floor$Hmax) + 1, na.rm = TRUE))
fcov_b  <- 1
sc_NULL <- list(name = "H1f_0_Null_baseline",
                  lai_fn = function(pr) lai_b,
                  hmax_fn = function(pr) hmax_b,
                  fcover_fn = function(pr) fcov_b,
                  lad_fn = make_lad_uniform)
sc_LVA_LAI <- list(name = "H1f_1_LAI_only",
                    lai_fn = function(pr) as.numeric(pr$LAI),
                    hmax_fn = function(pr) hmax_b,
                    fcover_fn = function(pr) fcov_b,
                    lad_fn = make_lad_uniform)

run_set <- function(sensors, sc) {
  cli_alert("[{sc$name}] : {length(sensors)} sensors")
  for (s in sensors) {
    pr <- as.data.frame(df_hobo_inputs[id_plot == s])
    out_nc <- here::here("out_files/musica_hobo_validation", sc$name,
                            sprintf("musica_out_HOBO_%s.nc", s))
    cli_alert("Run {s}")
    run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
    cli_alert("  -> size {file.size(out_nc)} bytes")
  }
}
run_set(need_null, sc_NULL)
run_set(need_lvalai, sc_LVA_LAI)
cli_alert_success("Done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
