# Force-delete the 22 corrupt legacy NCs and re-run them fresh.
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

cli_h1("Force-delete + rerun 22 missing HOBO sims")
t_start <- Sys.time()

# Sensors needing each coalition
need_null     <- c("41_04","41_09","41_10","41_11","41_17","41_26",
                     "41_39","41_40","41_54","41_59")  # 10
need_lvalai   <- c("41_04","41_09","41_10","41_11","41_17","41_26",
                     "41_36","41_39","41_40","41_42","41_54","41_59")  # 12

# Verify corrupt before delete : check time_len in NC
verify_corrupt <- function(nc_path) {
  if (!file.exists(nc_path)) return(TRUE)
  n <- tryCatch(nc_open(nc_path), error = function(e) NULL)
  if (is.null(n)) return(TRUE)
  ok <- n$dim$time$len > 0
  nc_close(n)
  !ok
}

# Delete corrupt files
for (s in need_null) {
  f <- here::here("out_files/musica_hobo_validation/H1f_0_Null_baseline",
                    sprintf("musica_out_HOBO_%s.nc", s))
  if (verify_corrupt(f)) {
    file.remove(f)
    cli_alert("Deleted corrupt {basename(f)} (NULL)")
  } else cli_alert("Skip {s} NULL : already valid")
}
for (s in need_lvalai) {
  f <- here::here("out_files/musica_hobo_validation/H1f_1_LAI_only",
                    sprintf("musica_out_HOBO_%s.nc", s))
  if (verify_corrupt(f)) {
    file.remove(f)
    cli_alert("Deleted corrupt {basename(f)} (LVA_LAI)")
  } else cli_alert("Skip {s} LVA_LAI : already valid")
}

# Rerun
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

for (sc in list(sc_NULL, sc_LVA_LAI)) {
  sensors <- if (sc$name == "H1f_0_Null_baseline") need_null else need_lvalai
  cli_alert("[{sc$name}] {length(sensors)} sims")
  for (s in sensors) {
    out_nc <- here::here("out_files/musica_hobo_validation", sc$name,
                            sprintf("musica_out_HOBO_%s.nc", s))
    if (file.exists(out_nc)) next   # Was either valid (kept) or already deleted+regen done
    pr <- as.data.frame(df_hobo_inputs[id_plot == s])
    cli_alert("  Run {s}")
    run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
    if (file.exists(out_nc)) {
      n <- nc_open(out_nc); tl <- n$dim$time$len; nc_close(n)
      cli_alert("    size={file.size(out_nc)}, time_len={tl}")
    }
  }
}
cli_alert_success("Total time : {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
