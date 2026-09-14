# ==============================================================================
# Chapter 3 — Shared configuration
# Sourced by both Chapter3_01_lai_corrections.R and Chapter3_main.R
# ==============================================================================

NC_FULL <- path.expand("~/Documents/NC_Full")
NC_BL   <- file.path(NC_FULL, "03_RESULTS/Blois/Metrics/Deciduous_Only")

CFG_C3 <- list(
  # Local paths
  in_dir              = "in_files",
  out_dir             = "out_files/Chapter3",
  forcing_file        = "in_files/FR-Blo_2021.nc",                 # FR-Blo (has h_sbl → iter)
  hobo_geojson        = "in_files/data_Blois_utm31n.geojson",
  hobo_temp_csv       = "in_files/Blois_data_temperature.csv",
  musica_cmd          = normalizePath("in_files/model-3.2.3/musica", mustWork = FALSE),  # v3.2.3
  abl_flag            = "iter",                                    # yoyo/iter (ABL coupling)

  # NC_Full data (read-only)
  nc_results_blois    = NC_BL,
  nc_lai_als_dopt     = file.path(NC_FULL, "output/intermediate/lai_als_dopt/Blois"),
  nc_s2_dopt          = file.path(NC_FULL, "output/intermediate/sm6/Blois"),

  # d_opt parameters (Blois, DSM_keepTrees, Pareto method)
  # "common" = all_sites d_opt = 7 m (aligned with Chapter 2 / sm6 pipeline)
  dopt_variant        = "common",     # "per_site" | "common" | "fixed_4" | "d20" (sensitivity)
  d_opt_m             = 7,            # effective depth in metres (all_sites Blois DSM)

  # S2 parameterisation tags
  # time series: Blois-specific parametrisation (per-date files)
  # summer snapshot LAI_S2_DOPT comes from dopt_variant = "common" (d_opt = 7 m)
  s2_dopt_distrib     = "optim_Blois",  # for DOPT time series

  # ATBD summer snapshot: prosail 3.0.0 ATBD_T (sm6a script 25), aligned with Ch2
  nc_s2_atbd          = file.path(NC_FULL, "output/intermediate/sm6/Blois",
                                   "s2lai_summer_atbd_T_res_10_m.tif"),

  # Temporal windows
  date_seq            = seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day"),
  ts_date_range       = c(as.Date("2021-04-01"), as.Date("2021-10-31")),
  list_year           = c(2021L, 2022L),

  # Processing
  agg_factor          = 2L,
  ids_to_remove       = c("41_13", "41_14", "41_20", "41_41",
                           "41_50", "41_51", "41_53"),

  # Scenario mode — controls which scenarios make_all_scenarios_c3() returns:
  #   "minimal" : 5 scenarios (STATIC_ALS reference + 2×2 ATBD/DOPT × static/dyn)
  #   "medium"  : adds DYN_ALS (H2 test) + DYN_S2_DOPT_LADOPT (LAD shape)
  #   "full"    : all 13-15 scenarios (original design)
  scenarios_mode      = "full",

  # Chapter 1 thermal sensitivity (°C per m²/m² LAI)
  # Fill in from the Chapter 1 GAMM marginal effect after running.
  # Negative = more LAI → more cooling.
  ch1_sensitivity     = NULL
)

# Convenience sub-paths
CFG_C3$lai_prep_dir <- file.path(CFG_C3$out_dir, "lai_prep")
CFG_C3$figs_dir     <- file.path(CFG_C3$out_dir, "figs")
CFG_C3$nc_dir       <- file.path(CFG_C3$out_dir, "nc")
CFG_C3$tables_dir   <- file.path(CFG_C3$out_dir, "tables")
