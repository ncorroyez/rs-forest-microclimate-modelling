# ==============================================================================
# Chapter 3 — Dynamic mapping of microclimatic attenuation
# Step 2 of 2: MuSICA simulations, validation, and spatial mapping
#
# PREREQUISITE: Chapter2bis_lai_corrections.R must have been run first.
#   It writes lai_prep_dir/{df_plots_lai.rds, ts_by_plot.rds, rf_loo_cv.rds}.
#
# Three FLAGS control execution:
#   DRY_RUN         = TRUE   → 3 plots only (verify MuSICA is callable)
#   RUN_SIMULATIONS = FALSE  → skip simulations (load existing nc/ outputs)
#   RUN_SPATIAL_MAP = FALSE  → skip domain-wide mapping
# ==============================================================================

# ---- 0. Libraries & source ---------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ncdf4)
  library(lubridate)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(mgcv)
  library(stringr)
  library(purrr)
  library(scales)
})

suppressPackageStartupMessages({
  library(rmusica)
  library(musica.tools)
})

# exclude h1_/lovb_ scripts: they carry top-level Chapter-1 re-run side effects
invisible(lapply(grep("/(h1_|lovb_)", list.files("R", pattern = "\\.R$", full.names = TRUE), value = TRUE, invert = TRUE), source))
source("Chapter3_config.R")

# ---- 1. Flags ----------------------------------------------------------------

FLAGS_C3 <- list(
  DRY_RUN          = FALSE,  # TRUE = 3 plots only
  RUN_SIMULATIONS  = TRUE,    # set TRUE when MuSICA is available
  RUN_TIMING_SHIFT = TRUE,    # heatwave temporal shift analysis (Bonus 2 / F6)
  RUN_SPATIAL_MAP  = FALSE,   # set TRUE after validation run
  SAVE_FIGURES     = TRUE
)

# ---- 2. Output directories ---------------------------------------------------

for (d in c(CFG_C3$out_dir,
            file.path(CFG_C3$out_dir, c("figs", "nc", "tables")))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 3. Load LAI prep outputs (from Chapter3_01_lai_corrections.R) -----------

cat("== Loading LAI preparation outputs ==\n")
prep       <- load_lai_prep(CFG_C3)
df_hobo    <- prep$df_plots
ts_list    <- prep$ts_by_plot

cat(sprintf("  %d HOBO plots loaded\n", nrow(df_hobo)))
cat(sprintf("  Time series variants: %s\n",
            paste(names(Filter(Negate(is.null), ts_list)), collapse = ", ")))
if (!is.null(prep$loo_cv)) {
  cat("  RF LOO-CV metrics:\n")
  print(attr(prep$loo_cv, "metrics"))
}

# ---- 4. Build scenario matrix ------------------------------------------------

cat("== Building Chapter 3 scenario matrix ==\n")
scenarios_c3 <- make_all_scenarios_c3(ts_list,
                                       list_year = CFG_C3$list_year,
                                       d_opt     = CFG_C3$d_opt_m,
                                       mode      = CFG_C3$scenarios_mode)

# Dry-run guard: 3 plots only
df_run <- if (FLAGS_C3$DRY_RUN) {
  cat("  [DRY_RUN] Limiting to 3 HOBO plots\n")
  df_hobo[seq_len(min(3L, nrow(df_hobo))), ]
} else {
  df_hobo
}

# ---- 5. MuSICA simulations ---------------------------------------------------

if (!FLAGS_C3$RUN_SIMULATIONS) {
  cat(paste0("\n[SKIP] FLAGS_C3$RUN_SIMULATIONS = FALSE\n",
             "  Set to TRUE to run MuSICA at all HOBO locations.\n\n"))
} else {
  cat("== Running MuSICA at HOBO locations ==\n")

  df_macro      <- extract_macro_daily(CFG_C3$forcing_file, CFG_C3$date_seq)
  df_hobo_daily <- read_hobo_daily(CFG_C3$hobo_temp_csv, CFG_C3$date_seq,
                                    df_macro, CFG_C3$ids_to_remove)

  val_out <- validate_scenarios_at_hobos(
    df_hobo_inputs = df_run,
    df_hobo_daily  = df_hobo_daily,
    scenarios      = scenarios_c3,
    parent_dir     = file.path(CFG_C3$out_dir, "nc"),
    df_macro       = df_macro,
    date_seq       = CFG_C3$date_seq,
    forcing_file   = CFG_C3$forcing_file,
    musica_cmd     = CFG_C3$musica_cmd,
    force          = FALSE,
    abl_flag       = CFG_C3$abl_flag        # v3.2.3 yoyo/iter
  )

  # ---- 6. Metrics table --------------------------------------------------------

  cat("\n== Validation metrics ==\n")
  print(val_out$metrics, n = Inf)
  write.csv(val_out$metrics,
            file.path(CFG_C3$out_dir, "tables", "c3_metrics_all.csv"),
            row.names = FALSE)

  # ---- 7. Validation figures ---------------------------------------------------

  cat("== Generating validation figures ==\n")

  sc_names <- unique(val_out$daily$scenario)
  is_dyn   <- grepl("^DYN", sc_names)
  n_stat   <- sum(!is_dyn); n_dyn <- sum(is_dyn)
  pal_sc   <- setNames(
    c(colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(n_stat + 2)[-(1:2)],
      colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n_dyn + 2)[-(1:2)]),
    c(sc_names[!is_dyn], sc_names[is_dyn])
  )
  figs_dir <- file.path(CFG_C3$out_dir, "figs")

  # RMSE bar chart
  p_rmse <- plot_hobo_validation(
    val_out$metrics,
    title = sprintf("Chapter 3 — RMSE per scenario (%d total)", nrow(val_out$metrics))
  )

  # Scatter matrix: obs vs sim ΔTmax, faceted by scenario
  p_scatter <- ggplot(val_out$daily,
                      aes(x = Delta_obs, y = Delta_sim, colour = scenario)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
    facet_wrap(~scenario, ncol = 4) +
    scale_colour_manual(values = pal_sc) +
    labs(title    = "Chapter 3 — Obs vs sim ΔTmax (all scenarios)",
         subtitle = "Summer 2021, Blois HOBO sensors",
         x = "Observed ΔTmax (°C)", y = "Simulated ΔTmax (°C)", colour = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "none")

  # Seasonal bias boxplot
  df_seasonal <- val_out$daily %>%
    mutate(
      month   = as.integer(format(date, "%m")),
      season  = dplyr::case_when(
        month %in% 3:5  ~ "Spring (MAM)",
        month %in% 6:8  ~ "Summer (JJA)",
        month %in% 9:11 ~ "Autumn (SON)",
        TRUE            ~ "Other"
      ),
      bias = Delta_sim - Delta_obs,
      mode = dplyr::if_else(grepl("^DYN", scenario), "Dynamic", "Static")
    )

  p_seasonal <- ggplot(df_seasonal, aes(x = season, y = bias, fill = scenario)) +
    geom_boxplot(alpha = 0.6, outlier.size = 0.6,
                 position = position_dodge(0.85), width = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~mode, ncol = 2) +
    scale_fill_manual(values = pal_sc) +
    labs(title = "Chapter 3 — Seasonal bias by scenario",
         x = NULL, y = "Bias sim − obs ΔTmax (°C)", fill = NULL) +
    theme_bw(base_size = 10) + theme(legend.position = "right")

  # Head-to-head: best static vs best dynamic (by RMSE)
  best_static  <- val_out$metrics %>%
    filter(!grepl("^DYN", scenario)) %>% slice_min(rmse, n = 1) %>% pull(scenario)
  best_dynamic <- val_out$metrics %>%
    filter(grepl("^DYN", scenario))  %>% slice_min(rmse, n = 1) %>% pull(scenario)

  p_h2h <- ggplot(
    val_out$daily %>% filter(scenario %in% c(best_static, best_dynamic)),
    aes(x = Delta_obs, y = Delta_sim, colour = scenario)
  ) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    scale_colour_manual(
      values = c(pal_sc[best_static], pal_sc[best_dynamic]),
      labels = c(sprintf("Best static: %s", best_static),
                 sprintf("Best dynamic: %s", best_dynamic))
    ) +
    labs(title    = "Chapter 3 — Best static vs best dynamic scenario",
         subtitle = sprintf("Static: %s | Dynamic: %s", best_static, best_dynamic),
         x = "Observed ΔTmax (°C)", y = "Simulated ΔTmax (°C)", colour = NULL) +
    theme_bw(base_size = 11)

  # ---- 7a-extra. Absolute Tmax scatter + micro vs macro -----------------------

  p_abs_scatter <- plot_tmax_absolute_scatter(val_out$daily, ncol = 5)
  p_micro_macro <- plot_tmax_micro_vs_macro(val_out$daily,   ncol = 5)

  if (FLAGS_C3$SAVE_FIGURES) {
    save_plot(p_rmse,        file.path(figs_dir, "c3_rmse_all.png"),
              width = 11, height = 7)
    save_plot(p_scatter,     file.path(figs_dir, "c3_scatter_matrix.png"),
              width = 18, height = 14, dpi = 300)
    save_plot(p_seasonal,    file.path(figs_dir, "c3_seasonal_bias.png"),
              width = 14, height = 6, dpi = 300)
    save_plot(p_h2h,         file.path(figs_dir, "c3_best_h2h.png"),
              width = 8, height = 6)
    save_plot(p_abs_scatter, file.path(figs_dir, "c3_tmax_absolute_scatter.png"),
              width = 18, height = 12, dpi = 200)
    save_plot(p_micro_macro, file.path(figs_dir, "c3_tmax_micro_vs_macro.png"),
              width = 18, height = 12, dpi = 200)
    cat("  Validation figures saved to", figs_dir, "\n")
  }

  # ---- 7b. Heatwave temporal shift (Bonus 2 — Figure F6) ---------------------

  if (FLAGS_C3$RUN_TIMING_SHIFT) {
    cat("\n== Heatwave temporal shift analysis ==\n")
    shift_res <- compute_tmax_shift_scenarios(
      scenario_names    = names(scenarios_c3),
      scenario_nc_root  = file.path(CFG_C3$out_dir, "nc"),
      csv_file          = CFG_C3$hobo_temp_csv,
      forcing_file      = CFG_C3$forcing_file,
      date_seq          = CFG_C3$date_seq,
      ids_to_remove     = CFG_C3$ids_to_remove,
      heatwave_threshold = 30
    )

    if (!is.null(shift_res)) {
      df_shift_summary <- summarise_tmax_shift(shift_res)
      cat("\n  Timing shift summary (hours vs ERA5):\n")
      print(df_shift_summary, n = Inf)
      write.csv(df_shift_summary,
                file.path(CFG_C3$tables_dir, "c3_tmax_shift_summary.csv"),
                row.names = FALSE)

      p_shift <- plot_tmax_shift(shift_res, pal_sc = pal_sc, heatwave_threshold = 30)
      if (FLAGS_C3$SAVE_FIGURES)
        save_plot(p_shift, file.path(figs_dir, "c3_tmax_shift_heatwave.png"),
                  width = 8, height = 5)
      cat(sprintf("  F6 saved: c3_tmax_shift_heatwave.png\n"))
    }
  }

  # ---- 8. Spatial mapping: domain-wide ΔTmax ----------------------------------

  if (FLAGS_C3$RUN_SPATIAL_MAP) {
    cat("== Spatial mapping — domain-wide ΔTmax ==\n")

    df_macro <- extract_macro_daily(CFG_C3$forcing_file, CFG_C3$date_seq)
    hw_days  <- pick_heatwave_days(df_macro, n = 3L)
    cat(sprintf("  Heatwave reference days: %s\n", paste(hw_days, collapse = ", ")))

    dom_rast <- load_domain_rasters(CFG_C3$in_dir)

    # Identify best scenario and its LAI column
    best_sc_name <- val_out$metrics %>% slice_min(rmse, n = 1) %>% pull(scenario)
    lai_col_best <- dplyr::case_when(
      grepl("ALS_DOPT",   best_sc_name) ~ "LAI_ALS_DOPT",
      grepl("ALS",        best_sc_name) ~ "LAI_ALS",
      grepl("S2_DOPT",    best_sc_name) ~ "LAI_S2_DOPT",
      grepl("S2_RESCALED",best_sc_name) ~ "LAI_S2_RESCALED",
      grepl("RF",         best_sc_name) ~ "LAI_RF",
      TRUE                              ~ "LAI_S2_ATBD"
    )

    df_for_emul <- val_out$daily %>%
      filter(scenario == best_sc_name) %>%
      left_join(df_hobo %>% dplyr::select(x, y, all_of(lai_col_best), Hmax, fCover),
                by = c("x", "y")) %>%
      rename(LAI = all_of(lai_col_best))

    cat(sprintf("  Fitting GAM emulator on %s (%d plot-days)...\n",
                best_sc_name, nrow(df_for_emul)))
    gam_emul <- fit_microclimate_emulator(df_for_emul, lai_col = "LAI", k = 5)
    cat(sprintf("  Emulator R² = %.3f (training)\n", summary(gam_emul)$r.sq))

    r_lai_ref  <- load_static_lai(
      file.path(CFG_C3$in_dir, "s2lai_summer_atbd_res_10_m.tif"), "LAI_S2_ATBD"
    )
    lai_path_key <- switch(lai_col_best,
      LAI_ALS       = "lai_als",
      LAI_ALS_DOPT  = "lai_als_dopt",
      LAI_S2_DOPT   = "lai_s2_dopt",
      "lai_s2_atbd"
    )
    r_lai_best <- load_static_lai(make_lai_paths(CFG_C3)[[lai_path_key]], lai_col_best)

    for (k in seq_along(hw_days)) {
      tmax_k   <- df_macro$Tmax_macro[df_macro$date == hw_days[k]]
      r_impact <- map_correction_impact(
        gam_emul, r_lai_ref, r_lai_best,
        dom_rast$hmax, dom_rast$fcover,
        tmax_day   = tmax_k,
        agg_factor = CFG_C3$agg_factor
      )
      r_pred_best <- predict_dtmax_domain(
        gam_emul, r_lai_best,
        dom_rast$hmax, dom_rast$fcover,
        tmax_day   = tmax_k,
        agg_factor = CFG_C3$agg_factor
      )

      p_map_impact <- plot_correction_impact(
        r_impact,
        title      = sprintf("ΔTmax impact: %s vs S2_ATBD (best = %s)",
                             lai_col_best, best_sc_name),
        date_label = sprintf("%s — ERA5 Tmax = %.1f°C", hw_days[k], tmax_k)
      )
      p_map_pred <- ggplot(
        as.data.frame(r_pred_best, xy = TRUE) %>% setNames(c("x", "y", "dtmax")),
        aes(x = x, y = y, fill = dtmax)
      ) +
        geom_raster() +
        scale_fill_viridis_c(option = "B", na.value = "grey90",
                              name = "ΔTmax\n(°C)") +
        coord_equal() +
        labs(title    = sprintf("Predicted ΔTmax — %s", best_sc_name),
             subtitle = sprintf("%s — ERA5 Tmax = %.1f°C", hw_days[k], tmax_k),
             x = "Easting (m)", y = "Northing (m)") +
        theme_bw(base_size = 10) +
        theme(axis.text = element_text(size = 7))

      p_maps <- p_map_pred + p_map_impact +
        patchwork::plot_annotation(
          title    = "Chapter 3 — Domain-wide microclimatic mapping",
          subtitle = sprintf("GAM emulator trained on HOBO validation | %d plots",
                             length(unique(df_for_emul$id_plot)))
        )

      if (FLAGS_C3$SAVE_FIGURES) {
        save_plot(p_maps,
                  file.path(figs_dir, sprintf("c3_map_%s.png", format(hw_days[k]))),
                  width = 14, height = 7, dpi = 300)
        terra::writeRaster(
          r_pred_best,
          file.path(figs_dir, sprintf("c3_dtmax_pred_%s_%s.tif",
                                       best_sc_name, hw_days[k])),
          overwrite = TRUE
        )
        terra::writeRaster(
          r_impact,
          file.path(figs_dir, sprintf("c3_correction_impact_%s.tif", hw_days[k])),
          overwrite = TRUE
        )
      }
    }
    cat(sprintf("  Spatial maps written (%d heatwave days)\n", length(hw_days)))
  }
}

cat("\n== Chapter 3 Step 2 complete ==\n")
cat("  Outputs in:", CFG_C3$out_dir, "\n")
