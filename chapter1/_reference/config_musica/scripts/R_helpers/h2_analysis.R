# ==============================================================================
# Chapter 1 — H2 analysis: uniform vs real LAD
# ==============================================================================

#' Pivot H2 results to wide format and compute per-plot, per-day difference.
#'
#' @param df_h2 Stacked H2 scenario results (from extract_all_scenarios).
#' @return Wide tibble with columns: x, y, date, Real, Uniform, diff.
summarise_h2 <- function(df_h2) {
  df_h2 %>% pivot_wider(id_cols = c(x, y, date), names_from = scenario, values_from = Delta_Tmax) %>%
    rename(Real = matches("real_LAD"), Uniform = matches("uniform_LAD")) %>% mutate(diff = Real - Uniform)
}

#' Density plot of ΔTmax distributions for real vs uniform LAD.
#'
#' @param df_h2 Stacked H2 scenario results.
#' @return ggplot object.
plot_h2_distributions <- function(df_h2) {
  ggplot(df_h2, aes(x = Delta_Tmax, fill = scenario, colour = scenario)) + geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values  = c("#31688e", "#d8576b")) + scale_colour_manual(values = c("#31688e", "#d8576b")) +
    labs(title = "H2 — ΔTmax distributions: real vs uniform LAD", x = expression(Delta * T[max] ~ (degree*C)), y = "Density")
}

#' Paired scatter plot of ΔTmax for two scenarios (hexbin + stats annotation).
#'
#' @param df_h2_wide  Wide H2 tibble (x, y, date, + two scenario columns).
#' @param col_x       Name of the x-axis scenario column (default "Real").
#' @param col_y       Name of the y-axis scenario column (default "Uniform").
#' @param title       Plot title.
#' @param subtitle    Plot subtitle.
#' @return ggplot object.
plot_h2_paired <- function(df_h2_wide, col_x = "Real", col_y = "Uniform",
                            title = "H2 — Per-plot, per-day ΔTmax: real vs uniform LAD",
                            subtitle = "Paired comparison: impact of vertical profile SHAPE on cooling") {
  x_vals <- df_h2_wide[[col_x]]; y_vals <- df_h2_wide[[col_y]]
  r2_val      <- cor(x_vals, y_vals, use = "complete.obs")^2
  rmse_val    <- sqrt(mean((y_vals - x_vals)^2, na.rm = TRUE))
  bias_val    <- mean(y_vals - x_vals, na.rm = TRUE)
  stats_label <- sprintf("R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", r2_val, rmse_val, bias_val)
  axis_lims   <- range(c(x_vals, y_vals), na.rm = TRUE)

  df_plot <- df_h2_wide; df_plot$.x <- x_vals; df_plot$.y <- y_vals
  ggplot(df_plot, aes(x = .x, y = .y)) +
    geom_hex(bins = 100) +
    scale_fill_gradient(low = "grey80", high = "midnightblue", trans = "log10", name = "Count\n(log)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", colour = "#d73027", se = FALSE, linewidth = 1.2) +
    coord_fixed(xlim = axis_lims, ylim = axis_lims) +
    annotate("text", x = -Inf, y = Inf, label = stats_label,
             hjust = -0.1, vjust = 1.2, fontface = "bold", size = 5) +
    labs(title = title, subtitle = subtitle,
         x = bquote(.(col_x) ~ Delta * T[max] ~ (degree*C)),
         y = bquote(.(col_y) ~ Delta * T[max] ~ (degree*C))) +
    theme_bw(base_size = 14) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
}

#' Split H2 differences by macroclimate regime (heatwave vs normal).
#'
#' @param df_h2_wide Wide H2 tibble with diff column.
#' @param df_macro   Daily macroclimate.
#' @param threshold  Tmax_macro threshold for heatwave (default 30°C).
#' @return Invisible ggplot (also prints summary table).
analyse_h2_distribution <- function(df_h2_wide, df_macro, threshold = 30) {
  df <- df_h2_wide %>% left_join(df_macro, by = "date") %>% mutate(period = ifelse(Tmax_macro >= threshold, "heatwave", "normal"))
  summ <- df %>% group_by(period) %>% summarise(n = n(), mean_diff = mean(diff, na.rm = TRUE), sd_diff = sd(diff, na.rm = TRUE), .groups = "drop")
  print(summ)
  p <- ggplot(df, aes(x = diff, fill = period)) + geom_histogram(bins = 80, alpha = 0.6, position = "identity") + geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c(normal = "#31688e", heatwave = "#d8576b")) + labs(title = "H2 — distribution of (Real - Uniform) ΔTmax", subtitle = "Split by macroclimate regime", x = "ΔTmax difference (°C)", y = "Count")
  print(p)
  invisible(p)
}

#' Seasonal stability of the H2 LAD shape effect (diff by month and regime).
#'
#' Answers: is the ±1°C real-vs-uniform effect stable across spring/summer/canicule,
#' or is it driven by a single extreme period?
#'
#' @param df_h2_wide Wide H2 tibble with diff column (Real - Uniform).
#' @param df_macro   Daily macroclimate dataframe (date, Tmax_macro).
#' @param threshold  Heatwave Tmax_macro threshold in °C (default 30).
#' @param out_path   PNG output path.
#' @return Invisible patchwork ggplot.
plot_h2_seasonal <- function(df_h2_wide, df_macro, threshold = 30,
                              out_path = "outputs/figures/annex/A17_h2_seasonal.png") {
  df <- df_h2_wide %>%
    inner_join(df_macro %>% dplyr::select(date, Tmax_macro), by = "date") %>%
    mutate(
      month  = factor(format(date, "%b"), levels = c("Jun", "Jul", "Aug", "Sep")),
      regime = factor(ifelse(Tmax_macro >= threshold, "Heatwave", "Normal"),
                      levels = c("Normal", "Heatwave"))
    )

  # Panel A: boxplot by month with mean annotations
  month_sum <- df %>%
    group_by(month) %>%
    summarise(
      mean_diff = mean(diff, na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    )

  pA <- ggplot(df, aes(x = month, y = diff, fill = month)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_text(
      data = month_sum,
      aes(x = month, y = 1.6,
          label = sprintf("mean=%.2f°C\nn=%d", mean_diff, n)),
      vjust = 0, size = 3.2, fontface = "bold", inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_fill_viridis_d(option = "mako", end = 0.85) +
    labs(
      title    = "A. LAD effect by month (outliers beyond ±2°C hidden)",
      x        = NULL,
      y        = "Diff Real − Uniform ΔTmax (°C)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  # Panel B: density by heatwave regime
  regime_sum <- df %>%
    group_by(regime) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE), n = dplyr::n(), .groups = "drop")

  pB <- ggplot(df, aes(x = diff, fill = regime, colour = regime)) +
    geom_density(alpha = 0.3, linewidth = 0.9) +
    geom_vline(data = regime_sum,
               aes(xintercept = mean_diff, colour = regime),
               linetype = "dashed", linewidth = 0.9) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.5) +
    geom_label(
      data = regime_sum %>%
        mutate(y_pos = ifelse(regime == "Heatwave", 3.5, 5.5)),
      aes(x = mean_diff, y = y_pos,
          label = sprintf("%s\n%.2f°C  n=%d", regime, mean_diff, n),
          colour = regime),
      fill = "white", label.size = 0.3, size = 3.2, fontface = "bold",
      inherit.aes = FALSE
    ) +
    coord_cartesian(xlim = c(-2, 2)) +
    scale_fill_manual(values   = c(Normal = "#31688e", Heatwave = "#d8576b")) +
    scale_colour_manual(values = c(Normal = "#31688e", Heatwave = "#d8576b")) +
    labs(
      title  = "B. Heatwave vs Normal distribution (±2°C shown)",
      x      = "Diff Real − Uniform ΔTmax (°C)",
      y      = "Density",
      fill   = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  subtitle_str <- if (all(c("Heatwave", "Normal") %in% regime_sum$regime)) {
    sprintf("Mean effect: Heatwave = %.2f°C | Normal = %.2f°C",
            regime_sum$mean_diff[regime_sum$regime == "Heatwave"],
            regime_sum$mean_diff[regime_sum$regime == "Normal"])
  } else {
    paste(sprintf("%s = %.2f°C", regime_sum$regime, regime_sum$mean_diff), collapse = " | ")
  }

  fig <- pA + pB +
    patchwork::plot_annotation(
      title    = "H2 — Seasonal stability of the LAD profile effect (Real − Uniform)",
      subtitle = subtitle_str,
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(colour = "grey30"))
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(fig, out_path, width = 13, height = 6)
  cat(sprintf("  [h2_seasonal] PNG : %s\n", out_path))
  invisible(fig)
}

#' Scatter H2 Real vs Uniform LAD, faceted by structural cluster.
#'
#' Shows that the H2 null result (negligible mean difference) holds within each
#' structural archetype. Rules out compensating effects between clusters masking
#' a real H2 signal.
#'
#' @param df_h2_wide Wide H2 tibble (x, y, date, Real, Uniform, diff).
#' @param df_sample  cLHS sample with Cluster column (x, y, Cluster).
#' @param out_path   PNG output path.
#' @return Invisible ggplot.
plot_h2_cluster_conditional <- function(df_h2_wide, df_sample,
                                         out_path = "outputs/figures/annex/A20_h2_cluster_conditional.png") {
  df <- df_h2_wide %>%
    inner_join(
      df_sample %>% dplyr::select(x, y, Cluster) %>% mutate(Cluster = as.factor(Cluster)),
      by = c("x", "y")
    )

  if (nrow(df) == 0) {
    warning("[plot_h2_cluster_conditional] join produced 0 rows — coordinate mismatch")
    return(invisible(NULL))
  }

  lims <- range(c(df$Real, df$Uniform), na.rm = TRUE)

  cl_stats <- df %>%
    group_by(Cluster) %>%
    summarise(
      n    = n(),
      r2   = cor(Real, Uniform, use = "complete.obs")^2,
      rmse = sqrt(mean((Uniform - Real)^2, na.rm = TRUE)),
      bias = mean(Uniform - Real, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("R²=%.2f\nRMSE=%.2f°C\nBias=%+.2f°C\nn=%d", r2, rmse, bias, n),
      x_pos = lims[1],
      y_pos = lims[2]
    )

  p <- ggplot(df, aes(x = Real, y = Uniform)) +
    geom_hex(bins = 60) +
    scale_fill_gradient(low = "grey85", high = "midnightblue",
                        trans = "log10", name = "Count\n(log)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.7) +
    geom_text(data = cl_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = -0.05, vjust = 1.1, fontface = "bold", size = 3.2,
              inherit.aes = FALSE) +
    coord_fixed(xlim = lims, ylim = lims) +
    facet_wrap(~ Cluster, labeller = label_both) +
    labs(
      title    = "H2 — Real vs Uniform LAD conditionnel par cluster structural",
      subtitle = "Le resultat nul H2 est verifie au sein de chaque archettype — pas d'effet compensatoire inter-clusters",
      x        = expression("Real LAD" ~ Delta * T[max] ~ (degree*C)),
      y        = expression("Uniform LAD" ~ Delta * T[max] ~ (degree*C))
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold"),
      legend.position = "right"
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 13, height = 10)
  cat(sprintf("  [h2_cluster_conditional] PNG : %s\n", out_path))
  invisible(p)
}

#' Compute paired validation metrics between HOBO observations and simulations.
#'
#' @param hobo_res List with element $daily (id_plot, date, scenario, Delta_obs, Delta_sim).
#' @return Tibble of per-scenario metrics: n, r2, rmse, mae, bias.
validate_hobo_paired <- function(hobo_res) {
  df <- hobo_res$daily
  key_counts <- df %>% count(id_plot, date, name = "n_sc")
  common_keys <- key_counts %>% filter(n_sc == length(unique(df$scenario)))
  df %>% inner_join(common_keys %>% select(id_plot, date), by = c("id_plot", "date")) %>% group_by(scenario) %>%
    summarise(n = n(), r2 = cor(Delta_obs, Delta_sim, use = "complete.obs")^2, rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)),
              mae  = mean(abs(Delta_obs - Delta_sim), na.rm = TRUE), bias = mean(Delta_sim - Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(rmse)
}

#' Analyse la structure de la diff (Real - Uniform) par archétype et métriques structurelles.
#'
#' Produit 4 livrables :
#'   1. Boxplot horizontal diff par archétype (PNG)
#'   2. Table heatwave × archétype (CSV + console)
#'   3. Scatter corrélation diff vs métriques structurelles (PNG + CSV)
#'   4. Histogramme amplitude temporelle de l'effet LAD par plot (PNG)
#'
#' @param df_h2_w      Dataframe wide H2 (x, y, date, Real, Uniform, diff).
#' @param df_sample    Dataframe cLHS (x, y, Cluster, LAI, Hmax, fCover, ...).
#' @param df_macro     Dataframe macroclimat journalier (date, Tmax_macro, ...).
#' @param h_median_vec Vecteur numérique H_median par plot (même ordre que
#'                     df_sample) ; NULL si non disponible.
#' @param fpc_scores   Matrice nx3 des scores FPC1/2/3 ; NULL si non disponible.
#' @param out_dir      Répertoire de sortie pour PNG et CSV.
#' @return Invisible list $per_plot, $cor_tbl, $tbl_cluster.
analyse_h2_structure <- function(df_h2_w, df_sample, df_macro,
                                  h_median_vec = NULL, fpc_scores = NULL,
                                  out_dir = "outputs/audit") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df_meta <- df_sample %>%
    dplyr::select(x, y, Cluster, LAI, Hmax, fCover) %>%
    mutate(Cluster = as.factor(Cluster))

  df_per_plot <- df_h2_w %>%
    group_by(x, y) %>%
    summarise(
      diff_mean = mean(diff, na.rm = TRUE),
      diff_q05  = quantile(diff, 0.05, na.rm = TRUE),
      diff_q95  = quantile(diff, 0.95, na.rm = TRUE),
      n_dates   = n(),
      .groups   = "drop"
    ) %>%
    inner_join(df_meta, by = c("x", "y")) %>%
    mutate(amplitude = diff_q95 - diff_q05)

  # Livrable 1 — Boxplot diff par cluster
  cl_stats <- df_per_plot %>%
    group_by(Cluster) %>%
    summarise(n = n(), median_diff = median(diff_mean, na.rm = TRUE), .groups = "drop")

  subtitle_1 <- paste(
    sprintf("Cl.%s : n=%d, méd.=%+.3f°C",
            cl_stats$Cluster, cl_stats$n, cl_stats$median_diff),
    collapse = " | "
  )

  p1 <- ggplot(df_per_plot,
               aes(x = diff_mean, y = Cluster, fill = Cluster, colour = Cluster)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(height = 0.2, alpha = 0.3, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    scale_fill_viridis_d(option = "turbo") +
    scale_colour_viridis_d(option = "turbo") +
    labs(
      title    = "H2 — Structure de l'effet LAD par cluster (per-plot mean)",
      subtitle = subtitle_1,
      x        = "diff Real - Uniform ΔTmax (°C)",
      y        = "Cluster"
    ) +
    theme(legend.position = "none")

  print(p1)
  ggsave(file.path(out_dir, "h2_boxplot_cluster.png"), p1,
         width = 10, height = 6, dpi = 150)

  # Livrable 2 — Table heatwave × cluster
  df_dates <- df_h2_w %>%
    left_join(df_macro %>% dplyr::select(date, Tmax_macro), by = "date") %>%
    mutate(period = ifelse(Tmax_macro >= 30, "heatwave", "normal")) %>%
    inner_join(df_meta %>% dplyr::select(x, y, Cluster), by = c("x", "y"))

  tbl2 <- df_dates %>%
    group_by(Cluster, period) %>%
    summarise(
      n           = n(),
      median_diff = median(diff, na.rm = TRUE),
      mean_diff   = mean(diff, na.rm = TRUE),
      sd_diff     = sd(diff, na.rm = TRUE),
      q05         = quantile(diff, 0.05, na.rm = TRUE),
      q95         = quantile(diff, 0.95, na.rm = TRUE),
      .groups     = "drop"
    ) %>%
    arrange(Cluster, period)

  write.csv(tbl2, file.path(out_dir, "h2_structure_by_cluster.csv"), row.names = FALSE)
  cat("\n── H2 Structure : diff par cluster × période ──\n")
  for (cl in levels(tbl2$Cluster)) {
    sub <- tbl2 %>% filter(Cluster == cl)
    cat(sprintf("  Cluster %s:\n", cl))
    for (i in seq_len(nrow(sub))) {
      cat(sprintf(
        "    %-10s  n=%5d  méd=%+.3f°C  moy=%+.3f°C  sd=%.3f°C  [q05=%+.3f ; q95=%+.3f]\n",
        sub$period[i], sub$n[i], sub$median_diff[i],
        sub$mean_diff[i], sub$sd_diff[i], sub$q05[i], sub$q95[i]
      ))
    }
  }
  cat(sprintf("  → CSV : %s/h2_structure_by_cluster.csv\n", out_dir))

  # Livrable 3 — Corrélation diff ↔ métriques structurelles
  df_corr <- df_per_plot %>% dplyr::select(x, y, diff_mean, LAI, Hmax, fCover)

  if (!is.null(h_median_vec) && length(h_median_vec) == nrow(df_sample)) {
    df_hmed <- df_sample %>%
      dplyr::select(x, y) %>%
      mutate(H_median = h_median_vec)
    df_corr <- df_corr %>% left_join(df_hmed, by = c("x", "y"))
  }
  if (!is.null(fpc_scores) && nrow(fpc_scores) == nrow(df_sample)) {
    df_fpc <- df_sample %>%
      dplyr::select(x, y) %>%
      mutate(FPC1 = fpc_scores[, 1],
             FPC2 = fpc_scores[, 2],
             FPC3 = fpc_scores[, 3])
    df_corr <- df_corr %>% left_join(df_fpc, by = c("x", "y"))
  }

  metric_cols <- setdiff(names(df_corr), c("x", "y", "diff_mean"))

  cor_rows <- lapply(metric_cols, function(m) {
    x_v <- df_corr[[m]]
    y_v <- df_corr$diff_mean
    ok  <- complete.cases(x_v, y_v)
    if (sum(ok) < 5) return(NULL)
    data.frame(
      metric   = m,
      pearson  = cor(x_v[ok], y_v[ok], method = "pearson"),
      spearman = cor(x_v[ok], y_v[ok], method = "spearman"),
      n        = sum(ok)
    )
  })
  df_cor_tbl <- bind_rows(cor_rows)
  write.csv(df_cor_tbl, file.path(out_dir, "h2_diff_correlations.csv"), row.names = FALSE)

  annot_df <- df_cor_tbl %>% mutate(label = sprintf("r = %.2f", pearson))

  df_long <- df_corr %>%
    pivot_longer(cols = all_of(metric_cols), names_to = "metric", values_to = "metric_val") %>%
    left_join(annot_df %>% dplyr::select(metric, label), by = "metric")

  p3 <- ggplot(df_long, aes(x = metric_val, y = diff_mean)) +
    geom_point(alpha = 0.4, size = 1.5, colour = "#31688e") +
    geom_smooth(method = "lm", se = TRUE, colour = "#d8576b",
                fill = "#d8576b", alpha = 0.15) +
    geom_text(
      data = df_long %>% group_by(metric, label) %>% slice(1),
      aes(x = -Inf, y = Inf, label = label),
      hjust = -0.1, vjust = 1.3, fontface = "bold", size = 3.5
    ) +
    facet_wrap(~ metric, scales = "free_x") +
    labs(
      title = "H2 — LAD diff (real-uniform) vs structural metrics per plot",
      x     = "Metric value",
      y     = "Mean diff ΔTmax (°C)"
    )

  print(p3)
  ggsave(file.path(out_dir, "h2_scatter_metrics.png"), p3,
         width = 12, height = 8, dpi = 150)
  cat(sprintf("  → CSV corrélations : %s/h2_diff_correlations.csv\n", out_dir))

  # Livrable 4 — Distribution amplitude temporelle
  pct_gt1 <- 100 * mean(df_per_plot$amplitude > 1, na.rm = TRUE)
  med_amp  <- median(df_per_plot$amplitude, na.rm = TRUE)

  p4 <- ggplot(df_per_plot, aes(x = amplitude, fill = Cluster)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
    geom_vline(xintercept = med_amp, linetype = "dashed",
               colour = "black", linewidth = 1) +
    scale_fill_viridis_d(option = "turbo") +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("%.1f%%\nplots amplitude > 1°C", pct_gt1),
             hjust = 1.1, vjust = 1.3, fontface = "bold", size = 4) +
    labs(
      title    = "H2 — Temporal amplitude of the LAD effect per plot (q95 - q05)",
      subtitle = sprintf("Global median = %.2f°C | %.1f%% plots with amplitude > 1°C",
                         med_amp, pct_gt1),
      x        = "LAD effect amplitude (q95 - q05) (°C)",
      y        = "Number of plots"
    )

  print(p4)
  ggsave(file.path(out_dir, "h2_amplitude_histogram.png"), p4,
         width = 10, height = 6, dpi = 150)

  cat(sprintf("\n[H2 structure] 5 livrables produits dans %s/\n", out_dir))
  invisible(list(per_plot = df_per_plot, cor_tbl = df_cor_tbl,
                 tbl_cluster = tbl2))
}
