# Run the 6 missing HOBO coalitions × 53 sensors = 318 MuSICA sims.
# These are the 2-trait-real coalitions : 0011, 0101, 0110, 1001, 1010, 1100.
# Output : out_files/musica_hobo_validation_floor05/HOBO_<bit>/musica_out_HOBO_<id>.nc
# 4 floor-affected sensors use floored fCover ; 49 others use real fCover.
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

cli_h1("HOBO missing 6 coalitions × 53 sensors = 318 sims")
t_start <- Sys.time()

# Load HOBO inputs and apply fCover floor for the 4 affected sensors
rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo_inputs <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                                    CFG$ids_to_remove))
df_hobo_floor <- copy(df_hobo_inputs)
df_hobo_floor[fCover < 0.5, fCover := 0.5]
df_hobo_floor <- as.data.frame(df_hobo_floor)
cli_alert("HOBO sensors : {nrow(df_hobo_floor)}  (4 with floored fCover)")

# df_sample for baseline computation
df_sample_floor <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))

# Build the 6 missing coalitions (same logic as scenarios_h1_factorial / lovb_06_hobo_run)
MISSING_BITS <- c("0011","0101","0110","1001","1010","1100")
lai_b   <- mean(df_sample_floor$LAI, na.rm = TRUE)
hmax_b  <- round(mean(floor(df_sample_floor$Hmax) + 1, na.rm = TRUE))
fcov_b  <- 1
cli_alert("Baselines : LAI_b={round(lai_b,3)}  Hmax_b={hmax_b}  fCover_b={fcov_b}  LAD_b=uniform")

mk_lai    <- function(real) if (real) function(pr) as.numeric(pr$LAI)    else function(pr) lai_b
mk_hmax   <- function(real) if (real) function(pr) as.numeric(pr$Hmax)   else function(pr) hmax_b
mk_fcover <- function(real) if (real) function(pr) as.numeric(pr$fCover) else function(pr) fcov_b
mk_lad    <- function(real) if (real) make_lad_real                       else make_lad_uniform

scs <- lapply(MISSING_BITS, function(bit) {
  b <- strsplit(bit, "", fixed = TRUE)[[1]]
  list(name = paste0("HOBO_", bit), bit_code = bit,
        lai_fn = mk_lai(b[1]=="1"), hmax_fn = mk_hmax(b[2]=="1"),
        fcover_fn = mk_fcover(b[3]=="1"), lad_fn = mk_lad(b[4]=="1"))
})
names(scs) <- vapply(scs, `[[`, "", "name")
cli_alert("Scenarios : {.val {names(scs)}}")

out_root <- here::here("out_files/musica_hobo_validation_floor05")
for (sc in scs) {
  out_dir <- file.path(out_root, sc$name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cli_alert("[{sc$name}] : {nrow(df_hobo_floor)} sensors")
  for (i in seq_len(nrow(df_hobo_floor))) {
    pr <- df_hobo_floor[i, , drop = FALSE]
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    run_musica_one(pr, sc, out_nc, CFG$forcing_file, CFG$musica_cmd)
  }
  cli_alert_success("Done {sc$name} at t = {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
}
cli_alert_success("All 318 HOBO sims done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
