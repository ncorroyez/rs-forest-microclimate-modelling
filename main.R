# ==============================================================================
# Chapter 1 — Entry point
# ==============================================================================

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

suppressPackageStartupMessages({
  library(musica.tools)
  library(rmusica)
  library(ncdf4)
  library(tidyverse)
  library(lubridate)
  library(broom.mixed)
  library(ggeffects)
  library(patchwork)
  library(viridis)
  library(mgcv)
  library(terra)
  library(fda)
  library(corrplot)
  library(sf)
  library(clhs)
  library(vip)
  library(Metrics)
  library(tidyterra)
})

theme_set(
  theme_bw() +
    theme(
      text             = element_text(size = 12),
      plot.title       = element_text(face = "bold", size = 14),
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95"),
      strip.text       = element_text(face = "bold")
    )
)

set.seed(42)

source("R/config.R")
source("R/io.R")
source("R/forest.R")
source("R/lad.R")
source("R/scenarios.R")
source("R/musica.R")
source("R/h2_analysis.R")
source("R/h1_analysis.R")
source("R/h1_shapley_conditional.R")
source("R/h1_shapley_archetypes.R")
source("R/fpca.R")
source("R/gamm.R")
source("R/validation.R")
source("R/sentinel.R")
source("R/profiles.R")
source("R/tests.R")

# ==============================================================================
# main() — orchestrator
# ==============================================================================

main <- function() {

  FIG     <- CFG$out_figures
  FIG_ANN <- file.path(CFG$out_figures, "annex")

  run_ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  setup_outputs_dirs()   # crée outputs/ et sous-répertoires silencieusement

  # ---- Logging : tee de tout le stdout + messages vers un fichier horodaté ----
  # sink(split=TRUE) : cat/print/writeLines → console ET fichier simultanément
  # sink(type="message") : message() → fichier (ne s'affiche plus en console)
  log_path <- file.path("outputs", sprintf("run_%s.log", run_ts))
  log_con  <- file(log_path, open = "wt")
  sink(log_con, split = TRUE)
  # NB : on ne redirige PAS les messages (type="message") pour que les erreurs
  # R restent visibles dans la console. Les message() sont rares (5 appels) et
  # apparaissent quand même dans le log via stderr si l'IDE le capture.
  on.exit({
    cat(sprintf("\n[log] Sauvegardé → %s\n", log_path))
    sink()
    close(log_con)
  }, add = TRUE)

  cat(sprintf("[0/7] Répertoires créés | log → %s\n", log_path))
  options(fig_review_dir = file.path("outputs", "figures_review", run_ts))
  cat(sprintf("  [review] Figures mirrorées → outputs/figures_review/%s/\n", run_ts))

  cat("[1/7] Loading LiDAR rasters...\n")
  rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)

  cat("[2/7] Building forest dataframe...\n")
  fd       <- build_forest_dataframe(rasters$stack)
  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

  # FPCA on the full forest (before clustering and cLHS) so that FPC1/FPC2
  # stratify both the clusters and the sample.
  cat("[2b] FPCA sur la forêt entière (avant clustering)...\n")
  fpca_res <- compute_fpca(fd$mat_lad, fd$df$Hmax, fd$z_breaks, n_harm = CFG$n_fpc_total)
  fpc_var   <- fpca_res$varprop
  n_fpc_sel <- select_n_fpc(fpc_var, min_marginal = CFG$fpc_min_marginal)
  cat(sprintf("    FPC variance (forêt entière) : %s\n",
              paste(round(100 * fpc_var, 1), collapse = " / ")))
  cat(sprintf("    FPCs retenus pour clustering/cLHS : %d (seuil gain >= %.0f%%)\n",
              n_fpc_sel, 100 * CFG$fpc_min_marginal))

  writeLines(
    c(sprintf("n_fpc_total: %d | n_fpc_selected: %d (min_marginal: %.0f%%)",
              CFG$n_fpc_total, n_fpc_sel, 100 * CFG$fpc_min_marginal),
      sprintf("FPC%d: %.1f%% (cumul %.1f%%) %s",
              seq_along(fpc_var),
              100 * fpc_var,
              100 * cumsum(fpc_var),
              ifelse(seq_along(fpc_var) <= n_fpc_sel, "<-- retenu", ""))),
    "outputs/fpca/fpca_variance.txt"
  )

  save_plot(fpca_res$gcv_plot, file.path(FIG_ANN, "A01_fpca_gcv.png"), width = 8, height = 5)
  print(fpca_res$gcv_plot)

  p_fpca_harm <- plot_fpc_harmonics(fpca_res, fd$mat_lad, fd$z_breaks, fd$df$Hmax)
  save_plot(p_fpca_harm, file.path(FIG, "09_fpca_harmonics.png"), width = 12, height = 6)
  print(p_fpca_harm)

  p_fpca_load <- plot_fpc_loadings(fpca_res)
  save_plot(p_fpca_load, file.path(FIG_ANN, "A03_fpca_loadings.png"), width = 10, height = 5)
  print(p_fpca_load)

  fpc_col_names <- paste0("FPC", seq_len(n_fpc_sel))
  fpc_mat       <- fpca_res$fpca$scores[, seq_len(n_fpc_sel), drop = FALSE]
  colnames(fpc_mat) <- fpc_col_names
  df_forest_raw <- fd$df %>% bind_cols(as.data.frame(fpc_mat))

  cluster_vars <- c("LAI", "Hmax", "fCover", fpc_col_names)
  cl_res    <- label_clusters(df_forest_raw, k = CFG$k_clusters, vars = cluster_vars)
  df_forest <- cl_res$df
  if (!is.null(cl_res$elbow_plot)) {
    save_plot(cl_res$elbow_plot, file.path(FIG_ANN, "A04_kmeans_elbow.png"), width = 8, height = 5)
    print(cl_res$elbow_plot)
  }

  p_fpca_scat <- plot_fpc_scatters(fpca_res, cluster_vec = df_forest$Cluster)
  save_plot(p_fpca_scat, file.path(FIG_ANN, "A02_fpca_scatters.png"), width = 12, height = 5)
  print(p_fpca_scat)

  p_fpca_recon <- plot_fpca_reconstruction(
    fpca_res,
    n_sites     = 6,
    cluster_vec = as.character(df_forest$Cluster),
    out_path    = file.path(FIG, "10_fpca_reconstruction.png")
  )
  print(p_fpca_recon)

  p_clust_prof <- plot_cluster_mean_profiles(df_forest, fpca_res$mat_rel, fpca_res$z_rel)
  save_plot(p_clust_prof, file.path(FIG, "08_clusters_lad_profiles.png"), width = 12, height = 6)
  print(p_clust_prof)

  # ---- 3. cLHS sampling -------------------------------------------------------
  if (FLAGS$RUN_CLHS) {
    cat("[3/7] cLHS sampling...\n")
    set.seed(42)
    df_sample <- sample_clhs_per_cluster(df_forest, CFG$n_per_cluster,
                                          vars = cluster_vars)
    saveRDS(df_sample, file.path(CFG$out_dir, "clhs_sample.rds"))
  } else {
    df_sample <- readRDS(file.path(CFG$out_dir, "clhs_sample.rds"))
    if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample)) {
      df_sample <- df_sample %>% rename(Cluster = Archetype)
    }
  }
  cat(sprintf("    sample size: %d plots\n", nrow(df_sample)))

  # Inspection : tous les profils LAD individuels par cluster (sample cLHS)
  p_lad_all <- plot_cluster_all_profiles(df_sample)
  save_plot(p_lad_all,
             file.path(FIG_ANN, "A_lad_profiles_per_cluster.png"),
             width = 14, height = 7)

  # Round df_sample coordinates to integers so they match NC filename integers
  # (musica.R creates filenames with round(plot$x), terra returns half-integer centroids).
  # Doing this once here fixes all downstream joins (join_with_sample, h2_analysis, gamm).
  df_sample <- df_sample %>% mutate(x = round(x), y = round(y))

  idx_sample <- match(paste(df_sample$x, df_sample$y),
                      paste(round(df_forest$x), round(df_forest$y)))
  n_na_idx <- sum(is.na(idx_sample))
  if (n_na_idx > 0)
    stop(sprintf("[idx_sample] %d/%d sample plots not found in df_forest — coordinate mismatch (rebuild cLHS with RUN_CLHS=TRUE)",
                 n_na_idx, nrow(df_sample)))
  fpc_scores_sample <- fpca_res$fpca$scores[idx_sample, seq_len(CFG$n_fpc_total), drop = FALSE]

  cat("[3b] Matrice de correlations structurelles (cLHS sample)...\n")
  corr_res <- plot_structural_correlations(
    df_sample,
    vars     = c("LAI", "Hmax", "fCover", "FPC1"),
    out_path = file.path(FIG_ANN, "A22_structural_correlations.png")
  )
  rho_lai_fcover <- if (!is.null(corr_res))
    corr_res$cor_matrix["LAI", "fCover"]
  else
    NA_real_

  ref_sc <- scenario_reference(df_sample)

  # ---- 4. H2 — uniform vs real LAD --------------------------------------------
  if (FLAGS$RUN_H2_UNIFORM_VS_REAL) {
    cat("[4/7] H2: uniform vs real LAD...\n")
    h2_scs    <- scenarios_h2_uniform_vs_real(df_sample)
    h2_paths  <- run_scenarios(df_sample, h2_scs, CFG$out_h2,
                               skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS)
    df_h2    <- extract_all_scenarios(h2_paths, df_macro, CFG$date_seq)
    df_h2_w  <- summarise_h2(df_h2)

    p_h2_dist <- plot_h2_distributions(df_h2)
    save_plot(p_h2_dist, file.path(FIG_ANN, "A06_h2_distributions.png"))
    print(p_h2_dist)

    p_h2_pair <- plot_h2_paired(df_h2_w)
    save_plot(p_h2_pair, file.path(FIG, "05_h2_paired_real_vs_uniform.png"))
    print(p_h2_pair)

    cat(sprintf("    mean per-plot diff (real - uniform): %.3f °C\n", mean(df_h2_w$diff, na.rm = TRUE)))
    p_h2_hw <- analyse_h2_distribution(df_h2_w, df_macro, CFG$heatwave_thr)
    save_plot(p_h2_hw, file.path(FIG_ANN, "A07_h2_heatwave.png"))

    cat("[H2 seasonal] Stabilité saisonnière de l'effet profil LAD...\n")
    plot_h2_seasonal(df_h2_w, df_macro, threshold = CFG$heatwave_thr,
                     out_path = file.path(FIG_ANN, "A17_h2_seasonal.png"))

    cat("[H2 cluster conditional] Résultat nul H2 conditionnel par cluster...\n")
    plot_h2_cluster_conditional(df_h2_w, df_sample,
                                 out_path = file.path(FIG_ANN, "A20_h2_cluster_conditional.png"))

    if (FLAGS$RUN_H2_STRUCTURE_ANALYSIS) {
      cat("[audit] Analyse structure de la diff H2...\n")
      fpc_sc_h2 <- if (exists("fpca_res") && !is.null(fpca_res))
                     fpca_res$fpca$scores[, 1:3] else NULL
      hmed_h2   <- if (exists("h_median_vec") && !is.null(h_median_vec))
                     h_median_vec else NULL
      analyse_h2_structure(
        df_h2_w, df_sample, df_macro,
        h_median_vec = hmed_h2,
        fpc_scores   = fpc_sc_h2,
        out_dir      = "outputs/audit"
      )
    }

    if (FLAGS$RUN_VERTICAL_PROFILES) {
      cat("[vertical] Profils verticaux Tmax par cluster médian...\n")
      plot_vertical_tmax_profiles(
        df_sample,
        real_lad_dir    = file.path(CFG$out_h2, "H2_real_LAD"),
        uniform_lad_dir = file.path(CFG$out_h2, "H2_uniform_LAD"),
        date_seq        = CFG$date_seq,
        out_dir         = "outputs/h2",
        summary_mode    = "mean_summer"
      )
      plot_vertical_tmax_profiles(
        df_sample,
        real_lad_dir    = file.path(CFG$out_h2, "H2_real_LAD"),
        uniform_lad_dir = file.path(CFG$out_h2, "H2_uniform_LAD"),
        date_seq        = CFG$date_seq,
        out_dir         = "outputs/h2",
        summary_mode    = "canicule_date"
      )
      plot_vertical_tmax_profiles(
        df_sample,
        real_lad_dir    = file.path(CFG$out_h2, "H2_real_LAD"),
        uniform_lad_dir = file.path(CFG$out_h2, "H2_uniform_LAD"),
        date_seq        = CFG$date_seq,
        out_dir         = "outputs/h2",
        summary_mode    = "median_date"
      )

      cat("[vertical] Profils LAD par cluster médian (tous clusters)...\n")
      medians_all <- select_cluster_median_plots(df_sample)
      if (!is.null(medians_all) && nrow(medians_all) > 0) {
        for (cl in unique(medians_all$Cluster)) {
          cl_median <- medians_all %>% filter(Cluster == cl)
          if (nrow(cl_median) == 1) {
            cl_full <- df_sample %>%
              filter(abs(x - cl_median$x) < 1 & abs(y - cl_median$y) < 1) %>%
              slice(1)
            if (nrow(cl_full) == 1) {
              plot_lad_profile_for_plot(
                cl_full,
                sprintf(file.path(FIG_ANN, "A08_c%s_median_lad_profile.png"), cl)
              )
            }
          }
        }
      }

      cat("[vertical] Diagnostic pic Cluster 1 (verify_cluster1_peak)...\n")
      verify_cluster1_peak(
        df_sample,
        real_lad_dir    = file.path(CFG$out_h2, "H2_real_LAD"),
        uniform_lad_dir = file.path(CFG$out_h2, "H2_uniform_LAD"),
        date_seq        = CFG$date_seq,
        out_dir         = "outputs/h2"
      )

      cat("[vertical] ΔTmax vertical par cluster (indépendant vs heure globale)...\n")
      plot_vertical_delta_by_cluster(
        df_sample,
        real_lad_dir = file.path(CFG$out_h2, "H2_real_LAD"),
        df_macro     = df_macro,
        date_seq     = CFG$date_seq,
        out_path     = file.path(FIG, "06_h2_vertical_delta_by_cluster.png")
      )

      cat("[vertical] Profils verticaux 3 scénarios H2 (canicule)...\n")
      plot_vertical_tmax_profiles_multi(
        scenarios_dict = list(
          "Real LAD"         = file.path(CFG$out_h2, "H2_real_LAD"),
          "Cluster-type LAD" = file.path(CFG$out_h2, "H2_cluster_type_LAD"),
          "Uniform LAD"      = file.path(CFG$out_h2, "H2_uniform_LAD")
        ),
        df_sample    = df_sample,
        date_seq     = CFG$date_seq,
        summary_mode = "canicule_date",
        diff_from    = "Uniform LAD",
        out_path     = file.path(FIG, "07_h2_vertical_3way_canicule.png")
      )
    }
  }

  # ---- 4b. H2 — cluster-type LAD vs real LAD ----------------------------------
  if (FLAGS$RUN_H2_CLUSTER_TYPE) {
    cat("[4b/7] H2: cluster-type LAD vs real LAD...\n")
    h2_ct_scs   <- scenarios_h2_cluster_type(df_sample)
    h2_ct_paths <- run_scenarios(df_sample, h2_ct_scs, CFG$out_h2,
                                 skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS)
    df_h2_ct    <- extract_all_scenarios(h2_ct_paths, df_macro, CFG$date_seq)
    df_h2_ct_w  <- df_h2_ct %>%
      pivot_wider(id_cols = c(x, y, date),
                  names_from  = scenario,
                  values_from = Delta_Tmax) %>%
      rename_with(~ "Real",        matches("real_LAD")) %>%
      rename_with(~ "ClusterType", matches("cluster_type_LAD")) %>%
      mutate(diff = Real - ClusterType)

    cat(sprintf("    mean per-plot diff (real - cluster-type): %.3f °C\n",
                mean(df_h2_ct_w$diff, na.rm = TRUE)))
    cat(sprintf("    sd diff: %.3f °C\n",
                sd(df_h2_ct_w$diff, na.rm = TRUE)))

    p_h2_ct_pair <- plot_h2_paired(df_h2_ct_w, col_x = "Real", col_y = "ClusterType",
                                    title    = "H2 — Per-plot, per-day ΔTmax: real vs cluster-type LAD",
                                    subtitle = "Paired comparison: impact of LAD shape typology on cooling")
    save_plot(p_h2_ct_pair, file.path(FIG_ANN, "A09_h2_ct_cluster_vs_real.png"))
    print(p_h2_ct_pair)

    summary_ct <- df_h2_ct_w %>%
      left_join(df_sample %>% select(x, y, Cluster), by = c("x", "y")) %>%
      group_by(Cluster) %>%
      summarise(
        mean_diff = mean(diff, na.rm = TRUE),
        sd_diff   = sd(diff,   na.rm = TRUE),
        n         = n(),
        .groups   = "drop"
      )
    cat("[4b] Diff (real - cluster-type) par cluster :\n")
    print(summary_ct)
    write.csv(summary_ct,
              "outputs/h2/h2_ct_summary_by_cluster.csv",
              row.names = FALSE)
  }

  # ---- 5. H1 — forward / LOO / factorial --------------------------------------
  scenario_paths <- list()
  ref_path <- file.path(CFG$out_h1_forward, ref_sc$name)
  run_musica_scenario(df_sample, ref_sc, ref_path,
                      skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS)
  scenario_paths[[ref_sc$name]] <- ref_path

  if (FLAGS$RUN_H1_FORWARD) {
    cat("[5a/7] H1 forward inclusion...\n")
    fw_scs <- scenarios_h1_forward(df_sample)

    sc0_params <- data.frame(
      variable = c("LAI_mean", "Hmax_mean", "fCover_absent", "LAD_absent"),
      valeur   = c(mean(df_sample$LAI, na.rm = TRUE),
                   round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE)),
                   1,
                   NA_real_),
      note     = c("moyenne échantillon cLHS",
                   "moyenne échantillon cLHS",
                   "canopée fermée homogène (pas de clumping)",
                   "profil uniforme 1 m → Hmax, densité = LAI_plot / depth")
    )
    write.csv(sc0_params, "outputs/h1/h1_sc0_baseline_params.csv", row.names = FALSE)
    cat(sprintf("  [SC0] LAI_mean=%.3f  Hmax_mean=%.2f m  fCover=1  LAD=uniforme\n",
                sc0_params$valeur[1], sc0_params$valeur[2]))
    scenario_paths <- c(scenario_paths,
      run_scenarios(df_sample, fw_scs, CFG$out_h1_forward,
                    skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS))
  }
  if (FLAGS$RUN_H1_LOO) {
    cat("[5b/7] H1 leave-one-out...\n")
    loo_scs <- scenarios_h1_loo(df_sample)
    scenario_paths <- c(scenario_paths,
      run_scenarios(df_sample, loo_scs, CFG$out_h1_loo,
                    skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS))
  }
  if (FLAGS$RUN_H1_FACTORIAL) {
    cat("[5c/7] H1 full factorial (HEAVY)...\n")
    fac_scs <- scenarios_h1_factorial(df_sample)

    factorial_preload_map <- c(
      H1F_m_m_a_u = file.path(CFG$out_h1_forward, "H1f_0_Null_baseline"),
      H1F_r_m_a_u = file.path(CFG$out_h1_forward, "H1f_1_LAI_only"),
      H1F_r_r_a_u = file.path(CFG$out_h1_forward, "H1f_2_LAI_Hmax"),
      H1F_r_r_r_u = file.path(CFG$out_h1_forward, "H1f_3_LAI_Hmax_fCover"),
      H1F_r_r_r_r = file.path(CFG$out_h1_forward, "H1f_4_Full_real"),
      H1F_m_r_r_r = file.path(CFG$out_h1_loo,     "H1l_dropLAI_mean"),
      H1F_r_m_r_r = file.path(CFG$out_h1_loo,     "H1l_dropHmax_mean"),
      H1F_r_r_a_r = file.path(CFG$out_h1_loo,     "H1l_dropfCover_mean")
    )
    preload_factorial_from_existing(factorial_preload_map, CFG$out_h1_factorial)

    scenario_paths <- c(scenario_paths,
      run_scenarios(df_sample, fac_scs, CFG$out_h1_factorial,
                    skip_existing = FLAGS$SKIP_EXISTING_SCENARIOS))
  }

  df_all_scenarios <- extract_all_scenarios(scenario_paths, df_macro, CFG$date_seq)
  df_scores        <- score_scenarios_vs_reference(df_all_scenarios, ref_sc$name)
  print(df_scores)
  write.csv(df_scores, "outputs/h1/h1_scores.csv", row.names = FALSE)

  # ---- Parsimony audit: conference-ready coalition stats (A24) ----------------
  # Encoding: bit order LAI | Hmax | fCover | LAD
  #   1000 = H1F_r_m_a_u = H1f_1_LAI_only
  #   1010 = H1F_r_m_r_u = LAI + fCover (S2+GEDI accessible)
  #   1110 = H1F_r_r_r_u = H1f_3_LAI_Hmax_fCover
  #   1111 = H1F_r_r_r_r = H1f_4_Full_real
  PARSIMONY_MAP <- data.frame(
    coalition = c("0000",                "1000",
                  "1010",                "1110",
                  "1111"),
    scenario  = c("H1f_0_Null_baseline", "H1f_1_LAI_only",
                  "H1F_r_m_r_u",         "H1f_3_LAI_Hmax_fCover",
                  "H1f_4_Full_real"),
    label     = c("Null baseline",       "LAI only (S2+GEDI)",
                  "LAI + fCover (S2+GEDI)", "LAI + Hmax + fCover",
                  "Full LiDAR"),
    stringsAsFactors = FALSE
  )

  parsimony_present <- PARSIMONY_MAP$scenario[PARSIMONY_MAP$scenario %in% df_scores$scenario]
  if (length(parsimony_present) >= 4) {
    df_parsimony <- df_scores %>%
      dplyr::filter(scenario %in% PARSIMONY_MAP$scenario) %>%
      dplyr::left_join(PARSIMONY_MAP, by = "scenario") %>%
      dplyr::select(scenario, coalition, label, n,
                    r2, rmse, mae, bias = mean_diff) %>%
      dplyr::arrange(match(scenario, PARSIMONY_MAP$scenario))

    cat("\n── Parsimony scenarios — RMSE & R² vs REF_all_real ──\n")
    print(as.data.frame(df_parsimony), digits = 4)

    # Key console numbers for conference Q&A
    get_rmse <- function(coal) {
      r <- df_parsimony$rmse[df_parsimony$coalition == coal]
      if (length(r) == 0) NA_real_ else r
    }
    rmse_null   <- get_rmse("0000")
    rmse_lai    <- get_rmse("1000")
    rmse_laifcv <- get_rmse("1010")
    rmse_full   <- get_rmse("1111")

    total_range <- rmse_null - rmse_full
    cat(sprintf("\nLAI alone vs null:         RMSE %+.3f°C (%+.1f%%)\n",
                rmse_lai - rmse_null,
                100 * (rmse_lai - rmse_null) / rmse_null))
    cat(sprintf("Adding fCover to LAI:      RMSE %+.3f°C (%+.1f%% of LAI-only)\n",
                rmse_laifcv - rmse_lai,
                100 * (rmse_laifcv - rmse_lai) / rmse_lai))
    if (!is.na(total_range) && total_range > 0) {
      pct_captured <- 100 * (rmse_null - rmse_laifcv) / total_range
      cat(sprintf("LAI+fCover captures:       %.1f%% of full LiDAR signal\n",
                  pct_captured))
    }

    dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)
    write.csv(df_parsimony, "outputs/audit/parsimony_scenarios.csv", row.names = FALSE)
    cat("  [parsimony] CSV: outputs/audit/parsimony_scenarios.csv\n")

    p_ladder <- plot_parsimony_ladder(df_scores)
    if (!is.null(p_ladder)) {
      save_plot(p_ladder, file.path(FIG_ANN, "A24_parsimony_ladder.png"),
                width = 9, height = 5)
      cat(sprintf("  [parsimony] PNG: %s\n",
                  file.path(FIG_ANN, "A24_parsimony_ladder.png")))
    }
  } else {
    message(sprintf(
      "[parsimony] Only %d/%d scenarios present — run with FLAGS$RUN_H1_FACTORIAL=TRUE to include H1F_r_m_r_u (1010).",
      length(parsimony_present), nrow(PARSIMONY_MAP)
    ))
  }

  # REF reconciliation: REF_all_real and H2_real_LAD share the same full-LiDAR
  # inputs — any divergence flags a reproducibility issue in the pipeline.
  # test_ref_h2_reconciliation() stops() the pipeline if RMSE > 1e-6.
  if (FLAGS$RUN_H2_UNIFORM_VS_REAL && exists("df_h2")) {
    cat("[audit] REF reconciliation: REF_all_real ↔ H2_real_LAD...\n")
    test_ref_h2_reconciliation(df_all_scenarios, df_h2, ref_name = ref_sc$name)

    df_ref_rec     <- df_all_scenarios %>%
      filter(scenario == ref_sc$name) %>%
      dplyr::select(x, y, date, Delta_REF = Delta_Tmax)
    df_h2_real_rec <- df_h2 %>%
      filter(scenario == "H2_real_LAD") %>%
      dplyr::select(x, y, date, Delta_H2 = Delta_Tmax)
    df_rec <- inner_join(df_ref_rec, df_h2_real_rec, by = c("x", "y", "date"))

    if (nrow(df_rec) > 0) {
      r2_rec   <- cor(df_rec$Delta_REF, df_rec$Delta_H2, use = "complete.obs")^2
      bias_rec <- mean(df_rec$Delta_H2 - df_rec$Delta_REF, na.rm = TRUE)
      rmse_rec <- sqrt(mean((df_rec$Delta_H2 - df_rec$Delta_REF)^2, na.rm = TRUE))

      lims_rec <- range(c(df_rec$Delta_REF, df_rec$Delta_H2), na.rm = TRUE)
      p_rec <- ggplot(df_rec, aes(x = Delta_REF, y = Delta_H2)) +
        geom_hex(bins = 80) +
        scale_fill_viridis_c(trans = "log10", name = "Count\n(log)") +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    colour = "black", linewidth = 0.8) +
        coord_fixed(xlim = lims_rec, ylim = lims_rec) +
        annotate("text", x = -Inf, y = Inf,
                 label = sprintf("R²=%.4f\nBias=%.4f°C\nRMSE=%.4f°C",
                                 r2_rec, bias_rec, rmse_rec),
                 hjust = -0.1, vjust = 1.2, fontface = "bold", size = 4) +
        labs(
          title    = "Contrôle cohérence : REF_all_real ↔ H2_real_LAD",
          subtitle = "Mêmes paramètres LiDAR — écart attendu ≈ 0",
          x        = expression("REF_all_real" ~ Delta * T[max] ~ (degree*C)),
          y        = expression("H2_real_LAD"  ~ Delta * T[max] ~ (degree*C))
        ) +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold"))
      save_plot(p_rec, file.path(FIG_ANN, "A16_ref_h2_reconciliation.png"),
                width = 7, height = 6)
      print(p_rec)
    }
  }

  p_h1_hier <- plot_scenario_hierarchy(df_scores, ref_sc$name)
  save_plot(p_h1_hier, file.path(FIG_ANN, "A05_h1_hierarchy.png"))
  print(p_h1_hier)

  # Variante dual : meme ordre Y (tri par RMSE), barres RMSE | barres R^2
  plot_scenario_hierarchy_dual(df_scores, ref_sc$name,
                                out_path = file.path(FIG_ANN, "A05b_h1_hierarchy_dual.png"))

  if (FLAGS$RUN_H1_FORWARD) {
    forward_order <- vapply(scenarios_h1_forward(df_sample), `[[`, "", "name")
    p_h1_fwd <- plot_forward_curve(df_scores, forward_order)
    save_plot(p_h1_fwd, file.path(FIG, "01_h1_forward_curve.png"))
    print(p_h1_fwd)

    forward_sc_names <- setdiff(forward_order, "H1f_0_Null_baseline")
    p_h1_sc <- plot_h1_scenario_scatters(
      df_all_scenarios,
      ref_scenario_name = ref_sc$name,
      forward_names     = forward_sc_names,
      out_path          = file.path(FIG, "04_h1_scatters_vs_ref.png")
    )
    print(p_h1_sc)

    p_h1_rank <- compare_h1_rankings(
      df_scores,
      forward_order = forward_order,
      out_path      = file.path(FIG, "02_h1_ranking.png")
    )
    print(p_h1_rank)

    # Wilcoxon test: is the marginal gain of +Profil LAD (H1f_3 -> H1f_4) significant?
    fwd_test_res <- test_forward_last_step(df_all_scenarios, ref_sc$name)
    cat(sprintf(
      "  [H1 fwd last step] Wilcoxon p=%.4f | Cohen's d=%.3f | median gain=%.4f°C | n=%d\n",
      fwd_test_res$p_value, fwd_test_res$cohen_d,
      fwd_test_res$median_gain_degC, fwd_test_res$n_pairs
    ))
    write.csv(
      data.frame(fwd_test_res),
      "outputs/h1/h1_fwd_last_step_test.csv",
      row.names = FALSE
    )
  }

  if (FLAGS$RUN_H1_LOO) {
    p_loo_2x2 <- plot_loo_scatters_2x2(
      df_all_scenarios,
      ref_scenario_name = ref_sc$name,
      out_path          = file.path(FIG, "03_h1_loo_scatters_2x2.png")
    )
    print(p_loo_2x2)
  }

  # ---- 5d. Shapley attribution (requires full factorial) ----------------------
  if (FLAGS$RUN_H1_FACTORIAL) {
    cat("[5d/7] Shapley exact attribution on 2^4 lattice...\n")
    fac_scs_named   <- scenarios_h1_factorial(df_sample)
    shap_coalitions <- extract_shapley_coalitions(df_scores, fac_scs_named)
    shap_res        <- compute_shapley_exact(shap_coalitions)
    plot_shapley_attribution(
      shap_res, df_scores,
      out_path = file.path(FIG, "02_shapley_attribution.png")
    )
    plot_shapley_attribution(
      shap_res, df_scores,
      out_path = file.path(FIG, "02_shapley_A.png"),
      panels   = "A"
    )
    plot_shapley_attribution(
      shap_res, df_scores,
      out_path = file.path(FIG, "02_shapley_C.png"),
      panels   = "C"
    )

    # ANOVA factorial decomposition — methodological cross-check (same 16 runs)
    cat("[5d+] ANOVA 2^4 factorial decomposition (triangulation Shapley)...\n")
    anova_res <- compute_anova_factorial(
      shap_coalitions,
      shap_res        = shap_res,
      rho_lai_fcover  = rho_lai_fcover,
      out_path        = file.path(FIG_ANN, "A21_anova_vs_shapley.png"),
      df_all          = df_all_scenarios,
      coalition_map   = setNames(
        vapply(fac_scs_named, function(sc) sc$name, character(1)),
        vapply(fac_scs_named, function(sc) {
          g <- sc$grid_row
          paste0(as.integer(g$LAI    == "real"),
                 as.integer(g$Hmax   == "real"),
                 as.integer(g$fCover == "real"),
                 as.integer(g$LAD    == "real"))
        }, character(1))
      ),
      cluster_vec     = setNames(df_sample$Cluster, paste(df_sample$x, df_sample$y)),
      n_boot          = 50
    )

    # Unit test: Shapley efficiency axiom
    run_unit_tests(df_sample, shap_res = shap_res,
                   df_all = if (exists("df_h2")) df_all_scenarios else NULL,
                   df_h2  = if (exists("df_h2")) df_h2 else NULL,
                   ref_name = ref_sc$name)

    if (FLAGS$RUN_SHAPLEY_BOOTSTRAP) {
      cat("[5e/7] Shapley bootstrap stability (n_boot=50, stratified per cluster)...\n")
      set.seed(42)
      df_boot_shap <- bootstrap_shapley(
        df_all_scenarios    = df_all_scenarios,
        df_sample           = df_sample,
        ref_scenario_name   = ref_sc$name,
        n_boot              = 50,
        factorial_scenarios = fac_scs_named,
        metric              = "rmse"
      )
      if (nrow(df_boot_shap) > 0) {
        write.csv(df_boot_shap,
                  "outputs/h1/shapley_bootstrap.csv", row.names = FALSE)
        plot_shapley_bootstrap(
          df_boot_shap, shap_exact = shap_res,
          out_path = file.path(FIG_ANN, "A18_shapley_bootstrap.png")
        )
      }
    }

    # ---- 5f. Shapley on R² (variance-explained attribution) -------------------
    # Bonus run : meme treillis 2^4, meme bootstrap, mais sur R^2 au lieu de RMSE.
    # v(N) = R^2_full_real - R^2_null. Chaque phi_i = gain de R^2 attribuable a i.
    cat("[5f/7] Shapley exact attribution on R^2 (metric: r2)...\n")
    shap_coalitions_r2 <- extract_shapley_coalitions(df_scores, fac_scs_named,
                                                       metric = "r2")
    shap_res_r2        <- compute_shapley_exact(shap_coalitions_r2)
    plot_shapley_attribution(
      shap_res_r2, df_scores,
      out_path = file.path(FIG, "02_shapley_attribution_r2.png")
    )

    if (FLAGS$RUN_SHAPLEY_BOOTSTRAP) {
      cat("[5f+] Shapley bootstrap R^2 (n_boot=50)...\n")
      set.seed(42)
      df_boot_shap_r2 <- bootstrap_shapley(
        df_all_scenarios    = df_all_scenarios,
        df_sample           = df_sample,
        ref_scenario_name   = ref_sc$name,
        n_boot              = 50,
        factorial_scenarios = fac_scs_named,
        metric              = "r2"
      )
      if (nrow(df_boot_shap_r2) > 0) {
        write.csv(df_boot_shap_r2,
                  "outputs/h1/shapley_bootstrap_r2.csv", row.names = FALSE)
        plot_shapley_bootstrap(
          df_boot_shap_r2, shap_exact = shap_res_r2,
          out_path = file.path(FIG_ANN, "A18b_shapley_bootstrap_r2.png")
        )
      }
    }

    # ---- 5g. Shapley 3-joueurs {LAI, fCover, Structure3D=(Hmax,LAD)} ----------
    # Sanity check : respecte le couplage structurel Hmax-LAD impose par
    # .lad_rescale(). Utilise 8 des 16 simulations factorielles existantes.
    cat("[5g/7] Shapley 3-joueurs {LAI, fCover, Structure3D=(Hmax,LAD)}...\n")
    shap_coal_3p     <- extract_shapley_coalitions_3p(df_scores, fac_scs_named,
                                                        metric = "rmse")
    shap_res_3p      <- compute_shapley_exact_3p(shap_coal_3p)
    plot_shapley_3p(shap_res_3p,
                     out_path = file.path(FIG, "02_shapley_3players.png"))

    shap_coal_3p_r2  <- extract_shapley_coalitions_3p(df_scores, fac_scs_named,
                                                        metric = "r2")
    shap_res_3p_r2   <- compute_shapley_exact_3p(shap_coal_3p_r2)
    plot_shapley_3p(shap_res_3p_r2,
                     out_path = file.path(FIG, "02_shapley_3players_r2.png"))

    # CSV combine 3p RMSE + R^2
    df_3p_combined <- data.frame(
      variable          = shap_res_3p$variables,
      shapley_rmse      = as.numeric(shap_res_3p$shapley),
      shapley_pct_rmse  = 100 * as.numeric(shap_res_3p$shapley) /
                                 shap_res_3p$total_effect,
      shapley_r2        = as.numeric(shap_res_3p_r2$shapley),
      shapley_pct_r2    = 100 * as.numeric(shap_res_3p_r2$shapley) /
                                 shap_res_3p_r2$total_effect
    )
    write.csv(df_3p_combined,
              "outputs/h1/shapley_3players.csv", row.names = FALSE)
    cat("  [Shapley 3p] CSV combine : outputs/h1/shapley_3players.csv\n")

    # ---- 5h. Conditional Shapley : phi_LAD stratified by Hmax level ----------
    # Decompose phi_LAD = (phi_A + phi_B) / 2 where Group A = Hmax coherent
    # (real) and Group B = Hmax incoherent (mean, .lad_rescale compresses).
    cat("[5h/7] Conditional Shapley for LAD stratified by Hmax...\n")
    shap_cond <- compute_shapley_conditional_LAD(shap_coalitions)
    cat(sprintf("  phi_LAD_global             = %+.4f C\n", shap_cond$phi_global))
    cat(sprintf("  phi_LAD | Hmax = real      = %+.4f C  (Group A : coherent domain)\n",
                shap_cond$phi_A))
    cat(sprintf("  phi_LAD | Hmax = baseline  = %+.4f C  (Group B : compressed by rescale)\n",
                shap_cond$phi_B))
    cat(sprintf("  Sanity check (A+B)/2 vs global : delta = %.2e\n\n",
                shap_cond$check_diff))

    if (FLAGS$RUN_SHAPLEY_BOOTSTRAP) {
      cat("[5h+] Conditional Shapley bootstrap (n=50, stratified by cluster)...\n")
      df_boot_cond <- bootstrap_shapley_conditional_LAD(
        df_all_scenarios    = df_all_scenarios,
        df_sample           = df_sample,
        ref_scenario_name   = ref_sc$name,
        n_boot              = 50,
        factorial_scenarios = fac_scs_named,
        metric              = "rmse"
      )
      if (nrow(df_boot_cond) > 0) {
        write.csv(df_boot_cond,
                  "outputs/h1/shapley_LAD_conditional_bootstrap.csv",
                  row.names = FALSE)
        plot_shapley_conditional_LAD(
          shap_cond, df_boot_cond,
          out_path = file.path(FIG_ANN, "A_shapley_LAD_conditional.png"),
          metric   = "rmse"
        )
        # Combined summary CSV (point estimate + CI95)
        df_summary_cond <- data.frame(
          quantity = c("phi_LAD_global",
                       "phi_LAD_given_Hmax_real",
                       "phi_LAD_given_Hmax_baseline"),
          phi      = c(shap_cond$phi_global, shap_cond$phi_A, shap_cond$phi_B),
          ci95_lo  = c(quantile(df_boot_cond$phi_global, 0.025, na.rm = TRUE),
                       quantile(df_boot_cond$phi_A,      0.025, na.rm = TRUE),
                       quantile(df_boot_cond$phi_B,      0.025, na.rm = TRUE)),
          ci95_hi  = c(quantile(df_boot_cond$phi_global, 0.975, na.rm = TRUE),
                       quantile(df_boot_cond$phi_A,      0.975, na.rm = TRUE),
                       quantile(df_boot_cond$phi_B,      0.975, na.rm = TRUE)),
          n_boot   = c(nrow(df_boot_cond), nrow(df_boot_cond), nrow(df_boot_cond))
        )
        write.csv(df_summary_cond,
                  "outputs/h1/shapley_LAD_conditional.csv", row.names = FALSE)
        cat("  [Shapley conditional LAD] CSV : outputs/h1/shapley_LAD_conditional.csv\n")
      }
    }

    # ---- 5i. Shapley archetypes (model-centric, 4 synthetic plots) ----------
    # 4 archetypes (1 par cluster) x 16 coalitions = 64 MuSICA runs.
    # Isole l'effet structural pur (pas de variance inter-plot intra-cluster).
    cat("[5i/7] Shapley on 4 synthetic archetypes...\n")
    df_archetypes <- make_synthetic_archetypes(df_sample)
    fac_scs_arch  <- build_factorial_scenarios_archetypes(df_sample)

    # ---- Choix du binaire MuSICA pour les archetypes ------------------------
    # Decoupe (binaire <-> NML+variables) :
    #   - Binaire = code Fortran (physique). Blois = pre-3.2.3, model-3.2.3 = recent.
    #   - NML + variables.csv = configuration. On utilise toujours Blois projet-root.
    #
    # POUR CHANGER LE BINAIRE : decommenter la ligne souhaitee ci-dessous.
    arch_musica_cmd <- normalizePath("in_files/Blois/musica",         mustWork = TRUE)
    # arch_musica_cmd <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
    cat(sprintf("  [archetypes] Binaire MuSICA : %s\n", arch_musica_cmd))

    arch_out_root <- "out_files/H1_archetypes"
    run_archetype_factorial(
      df_archetypes      = df_archetypes,
      fac_scs            = fac_scs_arch,
      out_root           = arch_out_root,
      forcing_file       = CFG$forcing_file,
      musica_cmd         = arch_musica_cmd,
      musica_nml         = "musica.nml",         # Blois NML (projet-root)
      musica_variables   = "./variables.csv"     # Blois variables (projet-root)
    )

    df_arch_delta <- extract_archetype_deltas(arch_out_root, df_macro, CFG$date_seq)
    cat(sprintf("  [archetype] %d (archetype, scenario, day) rows extracted\n",
                nrow(df_arch_delta)))

    shap_arch_point  <- compute_archetype_shapley(df_arch_delta)
    shap_pooled_point <- compute_pooled_archetype_shapley(df_arch_delta)
    shap_all_point   <- rbind(shap_arch_point, shap_pooled_point)
    cat("\n=== Shapley archetypes + pooled (point estimate) ===\n")
    print(shap_all_point)

    shap_arch_boot   <- bootstrap_archetype_shapley_temporal(df_arch_delta,
                                                              n_boot = 200)
    shap_pooled_boot <- bootstrap_pooled_archetype_shapley_temporal(df_arch_delta,
                                                                     n_boot = 200)
    shap_all_boot    <- rbind(shap_arch_boot, shap_pooled_boot)
    cat("\n[Note] Bootstrap on days, not on plots — reflects temporal variability\n",
        "       for a given archetype, not structural uncertainty within Blois.\n\n")

    # Summary CSV avec CI (inclut pooled)
    df_arch_summary <- merge(
      shap_all_point,
      shap_all_boot %>%
        dplyr::group_by(archetype, variable) %>%
        dplyr::summarise(ci95_lo = quantile(phi, 0.025, na.rm = TRUE),
                          ci95_hi = quantile(phi, 0.975, na.rm = TRUE),
                          n_boot  = dplyr::n(),
                          .groups = "drop"),
      by = c("archetype", "variable")
    )
    dir.create("outputs/h1", recursive = TRUE, showWarnings = FALSE)
    write.csv(df_arch_summary,
              "outputs/h1/shapley_archetypes_values.csv", row.names = FALSE)
    cat("  [archetypes] CSV : outputs/h1/shapley_archetypes_values.csv\n")

    # Comparaison cLHS vs archetypes (per-cluster + pooled + mean)
    df_clhs <- read.csv("outputs/figures/02_shapley_attribution.csv",
                          stringsAsFactors = FALSE)

    df_comp <- data.frame(
      variable          = df_clhs$variable,
      phi_cLHS_400plots = df_clhs$shapley,
      stringsAsFactors  = FALSE
    )
    # Per-cluster archetypes
    for (arch in sort(unique(shap_arch_point$archetype))) {
      col_nm <- paste0("phi_", arch)
      df_comp[[col_nm]] <- shap_arch_point$phi[
        match(df_comp$variable,
              shap_arch_point$variable[shap_arch_point$archetype == arch])]
    }
    # Mean across the 4 archetypes
    arch_cols <- grep("^phi_Arch_C", names(df_comp), value = TRUE)
    df_comp$phi_mean_archetypes <- rowMeans(df_comp[, arch_cols, drop = FALSE],
                                              na.rm = TRUE)
    # POOLED archetypes (4 archetypes treated as one sample)
    df_comp$phi_pooled_4archetypes <- shap_pooled_point$phi[
      match(df_comp$variable, shap_pooled_point$variable)]

    cat("\n=== Comparison cLHS vs per-archetype vs pooled vs mean ===\n")
    print(df_comp)
    write.csv(df_comp, "outputs/h1/shapley_archetypes_vs_cLHS.csv",
              row.names = FALSE)

    # Plots (now include pooled as 5th panel/column)
    plot_shapley_archetypes_barplot(
      shap_all_point, shap_all_boot,
      out_path = file.path(FIG_ANN, "A_shapley_archetypes_barplot.png")
    )
    plot_shapley_archetypes_comparison(
      shap_all_point, df_clhs,
      out_path = file.path(FIG_ANN, "A_shapley_archetypes_comparison.png")
    )
  }

  # ---- 6. GAMM emulator -------------------------------------------------------
  if (FLAGS$RUN_GAMM_EMULATOR) {
    cat("[6/7] GAMM emulator on reference scenario...\n")
    df_ref_daily <- df_all_scenarios %>% filter(scenario == ref_sc$name) %>% join_with_sample(df_sample)

    # Always compute both shape predictors so compare_gamm_variants() can run.
    mat_s        <- df_sample %>% dplyr::select(starts_with("LAD_Layer_")) %>% as.matrix()
    mat_s[is.na(mat_s)] <- 0
    h_median_vec_all <- compute_h_median(mat_s, as.numeric(gsub("LAD_Layer_", "", colnames(mat_s))))

    df_gamm <- prepare_gamm_data(df_ref_daily, df_sample,
                                  fpc_scores   = fpc_scores_sample,
                                  h_median_vec = h_median_vec_all)

    # Compare 3 shape-predictor variants: NONE / H_MEDIAN / FPC
    gamm_variants <- compare_gamm_variants(
      df_gamm,
      clim_vars = FLAGS$GAMM_CLIM_VARS,
      out_path  = file.path(FIG_ANN, "A10b_gamm_variants.png")
    )
    write.csv(gamm_variants$comparison,
              "outputs/gamm/gamm_variants_comparison.csv", row.names = FALSE)

    # Primary model for effects + diagnostics (driven by FLAGS$GAMM_SHAPE_VAR)
    gam_ref <- gamm_variants$models[[FLAGS$GAMM_SHAPE_VAR]]
    if (is.null(gam_ref)) {
      available <- names(gamm_variants$models)
      warning(sprintf("[GAMM] GAMM_SHAPE_VAR='%s' not available — falling back to '%s'",
                      FLAGS$GAMM_SHAPE_VAR, available[1]))
      gam_ref <- gamm_variants$models[[available[1]]]
    }
    capture.output(summary(gam_ref), file = "outputs/gamm/gamm_summary.txt")
    print(summary(gam_ref))

    p_gamm_eff <- plot_gamm_marginal_effects(gam_ref,
                                              shape_type = FLAGS$GAMM_SHAPE_VAR,
                                              clim_vars  = FLAGS$GAMM_CLIM_VARS)
    save_plot(p_gamm_eff, file.path(FIG_ANN, "A10_gamm_effects.png"), width = 12, height = 8)
    print(p_gamm_eff)

    png(file.path(FIG_ANN, "A11_gamm_residuals.png"), width = 1800, height = 900, res = 120)
    diag_res <- diagnose_residuals(gam_ref, df_gamm)
    dev.off()
    print(diag_res)
    write.csv(diag_res, "outputs/gamm/gamm_residuals_stats.csv", row.names = FALSE)

    conc_res <- run_concurvity_audit(gam_ref, out_dir = "outputs/audit")

    # ---- FPC1 robustness to Hmax collinearity (A23) ----------------------------
    # ρ(FPC1, Hmax) = 0.65 in our sample; GAMM concurvity s(FPC1_sc)~s(Hmax_sc) = 0.56.
    # Null H: the near-zero FPC1 effect is an artefact absorbed by Hmax.
    # Test: refit the FPC model without s(Hmax_sc) and compare the FPC1 smooth.
    gam_fpc_full <- gamm_variants$models[["FPC"]]
    if (!is.null(gam_fpc_full)) {
      cat("  [GAMM] FPC1 robustness test: fitting no-Hmax variant...\n")
      gam_fpc_nohmax <- fit_gamm_without_hmax(df_gamm, clim_vars = FLAGS$GAMM_CLIM_VARS)
      rob_res <- plot_fpc1_robustness_to_hmax(gam_fpc_full, gam_fpc_nohmax)

      save_plot(rob_res$plot, file.path(FIG_ANN, "A23_fpc1_robustness_to_hmax.png"),
                width = 12, height = 6)
      dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)
      write.csv(rob_res$df_smooth, "outputs/audit/fpc1_robustness.csv", row.names = FALSE)
      cat(sprintf(
        "FPC1 amplitude with Hmax = %.2f°C | without Hmax = %.2f°C | ratio = %.2f\n",
        rob_res$amp_with, rob_res$amp_without, rob_res$amp_with / rob_res$amp_without
      ))
    }

    # ---- LAI robustness to fCover collinearity (A23b) --------------------------
    # ρ(LAI, fCover) = 0.842 in our sample; GAMM concurvity worst pairwise = 0.842.
    # Null H: the LAI smooth is inflated because fCover absorbs overlapping variance.
    # Test: refit the NONE model without s(fCover_sc) and compare the LAI smooth.
    # Expected: LAI smooth shape and amplitude stable or slightly larger — qualitative
    # conclusion (LAI dominates) holds.
    gam_none_full <- gamm_variants$models[["NONE"]]
    if (!is.null(gam_none_full)) {
      cat("  [GAMM] LAI robustness test: fitting no-fCover variant...\n")
      gam_none_nofcover <- fit_gamm_without_fcover(df_gamm, clim_vars = FLAGS$GAMM_CLIM_VARS)
      rob_fcover <- plot_lai_robustness_to_fcover(gam_none_full, gam_none_nofcover)

      save_plot(rob_fcover$plot, file.path(FIG_ANN, "A23b_lai_robustness_to_fcover.png"),
                width = 12, height = 6)
      dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)
      write.csv(rob_fcover$df_smooth, "outputs/audit/lai_robustness.csv", row.names = FALSE)
      cat(sprintf(
        "LAI amplitude with fCover = %.2f°C | without fCover = %.2f°C | ratio = %.2f\n",
        rob_fcover$amp_with, rob_fcover$amp_without,
        rob_fcover$amp_with / rob_fcover$amp_without
      ))
      cat(sprintf(
        "Concurvity diagnostic: ρ(LAI,fCover) = 0.842 | amplitude ratio %.2f -> %s\n",
        rob_fcover$amp_with / rob_fcover$amp_without,
        ifelse(rob_fcover$amp_with / rob_fcover$amp_without > 0.7,
               "LAI smooth robust — qualitative conclusion holds",
               "WARNING: LAI amplitude drops >30% without fCover — interpret with caution")
      ))
    }
  }

  # ---- OPTIONAL: HOBO validation ----------------------------------------------
  if (FLAGS$RUN_HOBO_VALIDATION) {
    cat("[HOBO] Running validation of forward scenarios at HOBO locations...\n")
    df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq, df_macro, CFG$ids_to_remove)
    df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove)
    cat(sprintf("    %d HOBO plots retained\n", nrow(df_hobo_inputs)))

    hobo_scs <- if (FLAGS$RUN_H1_FORWARD) c(list(REF = ref_sc), scenarios_h1_forward(df_sample)) else list(REF = ref_sc)

    hobo_res <- validate_scenarios_at_hobos(df_hobo_inputs, df_hobo_daily, hobo_scs,
                                             parent_dir = CFG$out_hobo,
                                             df_macro   = df_macro,
                                             date_seq   = CFG$date_seq)
    print(hobo_res$metrics)
    write.csv(hobo_res$metrics, "outputs/hobo/hobo_metrics.csv", row.names = FALSE)

    p_hobo_val <- plot_hobo_validation(hobo_res$metrics)
    save_plot(p_hobo_val, file.path(FIG_ANN, "A12_hobo_validation.png"))
    print(p_hobo_val)

    fw_names  <- vapply(scenarios_h1_forward(df_sample), `[[`, "", "name")
    rep_hobos <- pick_representative_hobos(df_hobo_daily)

    df_real <- hobo_res$daily %>% filter(scenario == fw_names["Full_real"]) %>% mutate(Tmax_micro = Delta_sim + Tmax_macro)
    df_unif <- hobo_res$daily %>% filter(scenario == fw_names["LAI_Hmax_fCover"]) %>% mutate(Tmax_micro = Delta_sim + Tmax_macro)

    p_ts <- plot_timeseries_faceted(rep_hobos, df_hobo_daily, df_real, df_unif, df_macro)
    save_plot(p_ts, file.path(FIG_ANN, "A13_hobo_timeseries.png"), width = 12, height = 10)
    print(p_ts)

    cat("[HOBO] Scatter sim vs obs coloré par cluster structurel...\n")
    plot_hobo_sim_obs_by_cluster(
      hobo_res       = hobo_res,
      df_hobo_inputs = df_hobo_inputs,
      df_forest      = df_forest,
      scenario_name  = ref_sc$name,
      threshold      = CFG$heatwave_thr,
      out_path       = file.path(FIG_ANN, "A19_hobo_cluster_scatter.png")
    )

    cat("[HOBO] Scatter sim vs obs facette par cluster (A19b)...\n")
    plot_hobo_per_cluster(
      hobo_res       = hobo_res,
      df_hobo_inputs = df_hobo_inputs,
      df_forest      = df_forest,
      scenario_name  = ref_sc$name,
      threshold      = CFG$heatwave_thr,
      out_path       = file.path(FIG_ANN, "A19b_hobo_per_cluster.png")
    )
  }

  # ---- OPTIONAL: Sentinel-2 + FORMS-H annex -----------------------------------
  if (FLAGS$RUN_S2_ANNEX) {
    cat("[annex] Running S2 + FORMS-H naive scenario...\n")
    annex    <- make_s2_formsh_scenario(df_sample, CFG$in_dir, CFG$agg_factor)
    p_annex  <- run_musica_scenario(annex$df_sample, annex$scenario, CFG$out_s2_annex)
    df_annex <- extract_deltatmax_scenario(p_annex, df_macro, CFG$date_seq)

    cat("    [annex] Computing metrics S2 vs Full LiDAR...\n")
    df_ref_compare <- df_all_scenarios %>%
      filter(scenario == ref_sc$name) %>%
      dplyr::select(x, y, date, Delta_ref = Delta_Tmax)
    df_s2_matched  <- df_annex %>% inner_join(df_ref_compare, by = c("x", "y", "date"))

    s2_r2   <- cor(df_s2_matched$Delta_Tmax, df_s2_matched$Delta_ref, use = "complete.obs")^2
    s2_rmse <- sqrt(mean((df_s2_matched$Delta_Tmax - df_s2_matched$Delta_ref)^2, na.rm = TRUE))
    s2_bias <- mean(df_s2_matched$Delta_Tmax - df_s2_matched$Delta_ref, na.rm = TRUE)

    stats_label_s2 <- sprintf("R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", s2_r2, s2_rmse, s2_bias)

    df_s2_eval  <- bind_rows(df_annex, df_all_scenarios %>% filter(scenario == ref_sc$name))
    df_plot_s2  <- df_s2_eval %>%
      mutate(scenario_label = ifelse(scenario == ref_sc$name,
                                     "LiDAR (Full 3D Real)",
                                     "Sentinel-2 + FORMS-H (2D Proxy)"))

    p_s2_density <- ggplot(df_plot_s2,
                            aes(x = Delta_Tmax, fill = scenario_label, colour = scenario_label)) +
      geom_density(alpha = 0.3, linewidth = 0.8) +
      scale_fill_manual(values   = c("LiDAR (Full 3D Real)" = "#31688e",
                                     "Sentinel-2 + FORMS-H (2D Proxy)" = "#fca50a")) +
      scale_colour_manual(values = c("LiDAR (Full 3D Real)" = "#31688e",
                                     "Sentinel-2 + FORMS-H (2D Proxy)" = "#fca50a")) +
      labs(title    = "Annex H3 — Optical Proxy vs Structural Reality (Density)",
           subtitle = sprintf("R² = %.2f | RMSE = %.2f°C | Bias = %.2f°C", s2_r2, s2_rmse, s2_bias),
           x        = expression(Delta * T[max] ~ (degree*C)),
           y        = "Density",
           fill     = "Input Source",
           colour   = "Input Source") +
      theme(legend.position = "bottom")

    save_plot(p_s2_density, file.path(FIG_ANN, "A14_s2_density.png"))
    print(p_s2_density)

    axis_lims_s2 <- range(c(df_s2_matched$Delta_ref, df_s2_matched$Delta_Tmax), na.rm = TRUE)

    p_s2_paired <- ggplot(df_s2_matched, aes(x = Delta_ref, y = Delta_Tmax)) +
      geom_hex(bins = 100) +
      scale_fill_gradient(low = "grey80", high = "midnightblue",
                          trans = "log10", name = "Count\n(log)") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = "black", linewidth = 0.8) +
      geom_smooth(method = "lm", colour = "#fca50a", se = FALSE, linewidth = 1.2) +
      coord_fixed(xlim = axis_lims_s2, ylim = axis_lims_s2) +
      annotate("text", x = -Inf, y = Inf, label = stats_label_s2,
               hjust = -0.1, vjust = 1.2, fontface = "bold", size = 5, colour = "black") +
      labs(title    = "Annex H3 — Per-plot, per-day ΔTmax: S2 Proxy vs LiDAR",
           subtitle = "Paired comparison: impact of 2D optical proxy on spatial prediction",
           x        = expression("LiDAR (Full 3D Real)" ~ Delta * T[max] ~ (degree*C)),
           y        = expression("Sentinel-2 + FORMS-H" ~ Delta * T[max] ~ (degree*C))) +
      theme_bw(base_size = 14) +
      theme(legend.position = "right", plot.title = element_text(face = "bold"))

    save_plot(p_s2_paired, file.path(FIG_ANN, "A15_s2_paired.png"))
    print(p_s2_paired)

    writeLines(
      sprintf("R2: %.4f\nRMSE: %.4f deg C\nBias: %.4f deg C", s2_r2, s2_rmse, s2_bias),
      "outputs/s2_annex/s2_metrics.txt"
    )
  }

  cat("\n[done]\n")
  invisible(NULL)
}

if (interactive() || identical(Sys.getenv("RUN_MAIN"), "1")) {
  main()
}
