# ==============================================================================
# Chapter 1 — Configuration : paths and run flags
# ==============================================================================

# ---- Paths -------------------------------------------------------------------
CFG <- list(
  in_dir          = "in_files",
  out_dir         = "out_files/Sensitivity_Analysis",
  forcing_file    = "in_files/musica_in_Blois.nc",
  hobo_geojson    = "in_files/data_Blois_utm31n.geojson",
  hobo_temp_csv   = "in_files/Blois_data_temperature.csv",
  # MuSICA binary path. The /home/corroyez/Documents/musica/musica binary
  # (md5 53072278, Nov 2024) is the one used for the legacy Oct 2025 simulations
  # (45/8 buf/amp pattern matching HOBO observations). The in_files/model-3.2.3/musica
  # (md5 8b5139e0, May 2026) is a DIFFERENT model variant — substituting it shifts
  # ΔTmax by up to 1.4 °C and inverts the buf/amp regime ratio.
  musica_cmd      = "/home/corroyez/Documents/musica/musica",
  date_seq        = seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day"),
  ids_to_remove   = c("41_13", "41_14", "41_20", "41_41", "41_50", "41_51", "41_53"),
  agg_factor      = 2,        # 10 m → 20 m aggregation
  n_per_cluster   = 100,      # cLHS sample size per cluster
  k_clusters      = NULL,     # K-means cluster assignment per plot
  n_fpc_total     = 6,        # FPCs calculés par compute_fpca (scree plot complet)
  fpc_min_marginal = 0.10,    # seuil gain marginal : inclure FPC_k si varprop[k] >= seuil
  heatwave_thr    = 30,       # °C threshold for heatwave subset
  tair_target_height = 1.0,   # m above ground — hauteur fixe d'extraction
                              # de Tair (alignee sur HOBO). Utilisee par
                              # extract_deltatmax_one() + timing.R via
                              # get_tair_at_z() — interpolation plot-invariante.
  out_h2          = "out_files/H2_uniform_vs_real",
  out_h1_forward  = "out_files/H1_forward",
  out_h1_loo      = "out_files/H1_loo",
  out_h1_factorial= "out_files/H1_factorial",
  out_hobo        = "out_files/musica_hobo_validation",
  out_s2_annex    = "out_files/Annex_S2_FORMSH",
  out_figures     = "outputs/figures"
)

# ---- Run flags ---------------------------------------------------------------
FLAGS <- list(
  RUN_CLHS                  = FALSE,  # FALSE = réutilise clhs_sample.rds existant
  SKIP_EXISTING_SCENARIOS   = F,  # TRUE = skip tout répertoire scénario déjà peuplé
  RUN_H2_UNIFORM_VS_REAL    = F,  # H2 first pass
  RUN_H2_CLUSTER_TYPE       = F,  # H2 cluster-type LAD vs real LAD
  RUN_H1_FORWARD            = F,  # forward inclusion
  RUN_H1_LOO                = F,  # leave-one-out

  RUN_GAMM_EMULATOR         = F,  # Statistical emulator
  # ERA5 drivers to include in the GAMM. Remove any variable to exclude it.
  # Leave empty (character(0)) to run structure-only (no macroclimate).
  # Note: VPD_mean has high concurvity with Rad_mean (0.717) and Tmax_macro (0.881).
  GAMM_CLIM_VARS            = c("Tmax_macro", "Wind_mean", "Rad_mean"),
  GAMM_SHAPE_VAR            = "FPC",  # Shape predictor: "FPC", "CLUSTER" (cat.), "H_MEDIAN", "NONE"

  RUN_VERTICAL_PROFILES     = F,  # profils verticaux Tmax par cluster médian (H2)
  RUN_H2_STRUCTURE_ANALYSIS = TRUE,  # post-processing léger (< 30 s, pas de simulation)
  RUN_HOBO_VALIDATION       = TRUE,  # MuSICA at HOBO locations
  RUN_S2_ANNEX              = TRUE,  # Sentinel-2 annex (H3)
  RUN_H1_FACTORIAL          = TRUE,  # heavy, off by default
  RUN_SHAPLEY_BOOTSTRAP     = TRUE   # bootstrap Shapley stability (requires RUN_H1_FACTORIAL)
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
