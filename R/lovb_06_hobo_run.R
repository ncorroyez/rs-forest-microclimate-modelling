# ==============================================================================
# LOVB analysis — Module 6 : run LOVB/LVA coalitions at 53 HOBO locations
#
# Reuses validate_scenarios_at_hobos() from R/validation.R to launch MuSICA at
# HOBO points for the 6 coalitions missing from the existing forward suite :
#   already simulated as forward H1 :
#     0000  ARCH_m_m_a_u   (H1f_0_Null_baseline)
#     1000  ARCH_r_m_a_u   (H1f_1_LAI_only)
#     1110  ARCH_r_r_r_u   (H1f_3_LAI_Hmax_fCover  =  LOVB_LAD)
#     1111  ARCH_r_r_r_r   (H1f_4_Full_real == REF_all_real)
#   needed for LOVB+LVA :
#     0001  LVA_LAD     |   0010  LVA_fCover  |   0100  LVA_Hmax
#     0111  LOVB_LAI    |   1011  LOVB_Hmax   |   1101  LOVB_fCover
#
# Total new MuSICA runs : 6 coalitions x 53 plots = 318 (~80 min single thread).
# Saved as out_files/musica_hobo_validation/HOBO_<bit_code>/musica_out_HOBO_*.nc
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
source(here::here("R/validation.R"))

# Build scenarios named "HOBO_<bit_code>" for the 6 missing coalitions.
# Baseline values follow the cLHS factorial convention exactly.
build_lovb_hobo_scenarios <- function(df_sample) {
  lai_b   <- mean(df_sample$LAI, na.rm = TRUE)
  hmax_b  <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  fcov_b  <- 1
  cli_alert(sprintf("Baselines  LAI=%.3f  Hmax=%g  fCover=%g  LAD=uniform",
                     lai_b, hmax_b, fcov_b))

  # Helper builders
  mk_lai    <- function(real) if (real) function(pr) as.numeric(pr$LAI)    else function(pr) lai_b
  mk_hmax   <- function(real) if (real) function(pr) as.numeric(pr$Hmax)   else function(pr) hmax_b
  mk_fcover <- function(real) if (real) function(pr) as.numeric(pr$fCover) else function(pr) fcov_b
  mk_lad    <- function(real) if (real) make_lad_real                       else make_lad_uniform

  # Coalitions needed (bit code "L H F D")
  needed <- c(
    "0001",  # LVA_LAD
    "0010",  # LVA_fCover
    "0100",  # LVA_Hmax
    "0111",  # LOVB_LAI
    "1011",  # LOVB_Hmax
    "1101"   # LOVB_fCover
  )

  scs <- lapply(needed, function(bit) {
    b <- strsplit(bit, "", fixed = TRUE)[[1]]
    list(
      name      = paste0("HOBO_", bit),
      bit_code  = bit,
      lai_fn    = mk_lai(b[1] == "1"),
      hmax_fn   = mk_hmax(b[2] == "1"),
      fcover_fn = mk_fcover(b[3] == "1"),
      lad_fn    = mk_lad(b[4] == "1")
    )
  })
  names(scs) <- vapply(scs, `[[`, "", "name")
  scs
}

# ---- Main runner -------------------------------------------------------------
lovb_run_hobo <- function(force = FALSE) {
  cli_h1("LOVB at HOBO sites — launch MuSICA on 6 missing coalitions")

  # Inputs
  rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
  df_macro       <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove)
  df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                     df_macro, CFG$ids_to_remove)
  df_sample      <- readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample.rds"))
  if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample)) {
    df_sample <- df_sample %>% dplyr::rename(Cluster = Archetype)
  }
  cli_alert("HOBO plots : {nrow(df_hobo_inputs)}  |  cLHS plots : {nrow(df_sample)}")

  scenarios <- build_lovb_hobo_scenarios(df_sample)
  cli_alert("Coalitions to run : {.val {names(scenarios)}}")

  out_root <- here::here("out_files/musica_hobo_validation")
  res <- validate_scenarios_at_hobos(
    df_hobo_inputs = df_hobo_inputs,
    df_hobo_daily  = df_hobo_daily,
    scenarios      = scenarios,
    parent_dir     = out_root,
    df_macro       = df_macro,
    date_seq       = CFG$date_seq,
    forcing_file   = CFG$forcing_file,
    musica_cmd     = CFG$musica_cmd,
    force          = force
  )

  cli_alert_success("HOBO LOVB runs done. Folders : {.path {out_root}/HOBO_*}")
  invisible(res)
}

# Run directly when sourced as a script :
if (sys.nframe() == 0) {
  lovb_run_hobo(force = FALSE)
}
