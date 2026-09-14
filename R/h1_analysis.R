# ==============================================================================
# Chapter 1 — H1 analysis: forward inclusion, LOO, and Shapley attribution
# ==============================================================================

#' Score scenarios by their divergence from a reference (full-real) scenario.
#'
#' @param df_all            Stacked scenario results (x, y, date, Delta_Tmax, scenario).
#' @param ref_scenario_name Name of the reference scenario.
#' @return Tibble sorted by RMSE with columns: scenario, n, mean_diff, mae, rmse.
score_scenarios_vs_reference <- function(df_all, ref_scenario_name) {
  df_ref <- df_all %>% filter(scenario == ref_scenario_name) %>% dplyr::select(x, y, date, Delta_ref = Delta_Tmax)
  df_all %>% inner_join(df_ref, by = c("x", "y", "date")) %>% mutate(diff = Delta_Tmax - Delta_ref) %>% group_by(scenario) %>%
    summarise(n = n(), mean_diff = mean(diff, na.rm = TRUE), mae = mean(abs(diff), na.rm = TRUE),
              rmse = sqrt(mean(diff^2, na.rm = TRUE)),
              r2   = cor(Delta_Tmax, Delta_ref, use = "complete.obs")^2,
              .groups = "drop") %>% arrange(rmse)
}

#' Bar chart of scenario RMSE vs reference, coloured by mean bias.
#'
#' @param df_scores         Tibble from score_scenarios_vs_reference.
#' @param ref_scenario_name Name of the reference scenario to exclude from the plot.
#' @return ggplot object.
plot_scenario_hierarchy <- function(df_scores, ref_scenario_name) {
  # Position-aware decoder : convertit le code H1F_<L>_<H>_<F>_<D> en label
  # explicite "Fact L:X H:X F:X D:X" sans ambiguïté.
  # Convention (cf. scenarios.R) :
  #   LAI    pos1 : m=mean,  r=real
  #   Hmax   pos2 : m=mean,  r=real
  #   fCover pos3 : a=absent (=1, canopée fermée),  r=real
  #   LAD    pos4 : m=mean-shape,  r=real,  u=uniform
  readable <- function(scenarios) {
    decode_factorial <- function(s) {
      parts <- strsplit(sub("^H1F_", "", s), "_", fixed = TRUE)[[1]]
      if (length(parts) != 4) return(s)
      lookup <- c(m = "μ", r = "r", u = "u", a = "∅")
      keys   <- c("L", "H", "F", "D")
      sprintf("Fact %s",
              paste(sprintf("%s:%s", keys, lookup[parts]), collapse = " "))
    }
    decode_loo <- function(s) {
      var <- sub("^H1l_drop([A-Za-z]+)_.*", "\\1", s)
      val <- sub("^H1l_drop[A-Za-z]+_(.*)$", "\\1", s)
      sprintf("LOO −%s%s", var,
              if (val == "mean") "" else sprintf(" (%s)", val))
    }
    decode_one <- function(s) {
      if (grepl("^H1f_0_Null_baseline$",    s)) return("Null baseline")
      if (grepl("^H1f_1_LAI_only$",         s)) return("Fwd1 LAI")
      if (grepl("^H1f_2_LAI_Hmax$",         s)) return("Fwd2 +Hmax")
      if (grepl("^H1f_3_LAI_Hmax_fCover$",  s)) return("Fwd3 +fCover")
      if (grepl("^H1f_4_Full_real$",        s)) return("Fwd4 +LAD")
      if (grepl("^H1F_",     s)) return(decode_factorial(s))
      if (grepl("^H1l_drop", s)) return(decode_loo(s))
      gsub("_", " ", s)
    }
    vapply(scenarios, decode_one, character(1))
  }

  df_plot <- df_scores %>%
    filter(scenario != ref_scenario_name) %>%
    mutate(label = readable(scenario))
  ggplot(df_plot, aes(x = reorder(label, rmse), y = rmse, fill = mean_diff)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "#31688e", mid = "white", high = "#d8576b",
                         midpoint = 0, name = "Mean bias\n(°C)") +
    labs(
      title   = "H1 — divergence from full-real reference",
      caption = paste0(
        "Fact = factorial 2^4 | Fwd = forward inclusion | LOO = leave-one-out\n",
        "L = LAI, H = Hmax, F = fCover, D = LAD profile\n",
        "μ = mean (baseline), r = real (LiDAR), u = uniform LAD, ∅ = absent (fCover = 1)"
      ),
      x = NULL, y = "RMSE vs reference (°C)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y  = element_text(size = 8, family = "mono"),
      plot.caption = element_text(size = 7.5, colour = "grey30", hjust = 0,
                                  lineheight = 1.3)
    )
}

#' Dual-panel scenario hierarchy : RMSE (left) vs R² (right).
#'
#' Même ordre Y des deux côtés (tri par RMSE) → permet de lire ligne par ligne
#' la cohérence entre les deux métriques. Si une ligne a RMSE bas mais R² bas
#' aussi, le modèle est précis en magnitude mais mauvais en pattern (et vice
#' versa). Lecture : un scénario peut être "bon RMSE / mauvais R²" si la sim
#' est proche en moyenne mais ne discrimine pas les plots entre eux.
#'
#' @param df_scores         Tibble from score_scenarios_vs_reference (avec colonnes
#'   rmse, r2, mean_diff).
#' @param ref_scenario_name Nom du scénario de référence à exclure.
#' @param out_path          Chemin PNG de sortie.
#' @return Invisible patchwork ggplot.
plot_scenario_hierarchy_dual <- function(df_scores, ref_scenario_name,
                                          out_path = "outputs/figures/annex/A05b_h1_hierarchy_dual.png") {
  # Reuse le decoder readable() local (copie courte)
  readable <- function(scenarios) {
    decode_factorial <- function(s) {
      parts <- strsplit(sub("^H1F_", "", s), "_", fixed = TRUE)[[1]]
      if (length(parts) != 4) return(s)
      lookup <- c(m = "μ", r = "r", u = "u", a = "∅")
      keys   <- c("L", "H", "F", "D")
      sprintf("Fact %s", paste(sprintf("%s:%s", keys, lookup[parts]), collapse = " "))
    }
    decode_loo <- function(s) {
      var <- sub("^H1l_drop([A-Za-z]+)_.*", "\\1", s)
      val <- sub("^H1l_drop[A-Za-z]+_(.*)$", "\\1", s)
      sprintf("LOO −%s%s", var, if (val == "mean") "" else sprintf(" (%s)", val))
    }
    decode_one <- function(s) {
      if (grepl("^H1f_0_Null_baseline$",   s)) return("Null baseline")
      if (grepl("^H1f_1_LAI_only$",        s)) return("Fwd1 LAI")
      if (grepl("^H1f_2_LAI_Hmax$",        s)) return("Fwd2 +Hmax")
      if (grepl("^H1f_3_LAI_Hmax_fCover$", s)) return("Fwd3 +fCover")
      if (grepl("^H1f_4_Full_real$",       s)) return("Fwd4 +LAD")
      if (grepl("^H1F_",     s)) return(decode_factorial(s))
      if (grepl("^H1l_drop", s)) return(decode_loo(s))
      gsub("_", " ", s)
    }
    vapply(scenarios, decode_one, character(1))
  }

  if (!"r2" %in% names(df_scores))
    stop("plot_scenario_hierarchy_dual: df_scores must contain an 'r2' column")

  df <- df_scores %>%
    dplyr::filter(scenario != ref_scenario_name) %>%
    dplyr::mutate(label = readable(scenario))

  # Ordre Y commun : tri par RMSE croissant → meilleur RMSE en bas après coord_flip
  order_labels <- df %>% dplyr::arrange(rmse) %>% dplyr::pull(label)
  df$label_f <- factor(df$label, levels = order_labels)

  bias_range <- max(abs(df$mean_diff), na.rm = TRUE)
  fill_scale <- ggplot2::scale_fill_gradient2(
    low = "#31688e", mid = "white", high = "#d8576b",
    midpoint = 0, limits = c(-bias_range, bias_range),
    name = "Mean bias\n(°C)"
  )

  base_theme <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8, family = "mono"),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )

  p_rmse <- ggplot2::ggplot(df, ggplot2::aes(x = label_f, y = rmse, fill = mean_diff)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    fill_scale +
    ggplot2::labs(title = "RMSE (lower = better)", x = NULL,
                  y = "RMSE vs reference (°C)") +
    base_theme

  p_r2 <- ggplot2::ggplot(df, ggplot2::aes(x = label_f, y = r2, fill = mean_diff)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    fill_scale +
    ggplot2::scale_y_continuous(limits = c(0, max(1, max(df$r2, na.rm = TRUE) * 1.05))) +
    ggplot2::labs(title = "R² (higher = better)", x = NULL,
                  y = "R² vs reference") +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())   # libre vu l'ordre commun

  fig <- (p_rmse | p_r2) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title    = "H1 — scenario hierarchy : RMSE vs R² (same Y-order, sorted by RMSE)",
      caption  = paste0(
        "Fact = factorial 2^4 | Fwd = forward inclusion | LOO = leave-one-out\n",
        "L = LAI, H = Hmax, F = fCover, D = LAD profile\n",
        "μ = mean (baseline), r = real (LiDAR), u = uniform LAD, ∅ = absent (fCover = 1)"
      ),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", size = 13),
        plot.caption  = ggplot2::element_text(size = 7.5, colour = "grey30",
                                              hjust = 0, lineheight = 1.3)
      )
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(fig, out_path, width = 14, height = 8)
  cat(sprintf("  [scenario_hierarchy_dual] PNG : %s\n", out_path))
  invisible(fig)
}

#' Line plot of RMSE along the forward inclusion sequence.
#'
#' @param df_scores    Tibble from score_scenarios_vs_reference.
#' @param forward_order Ordered character vector of scenario names.
#' @return ggplot object.
plot_forward_curve <- function(df_scores, forward_order) {
  short_labels <- c(
    "H1f_0_Null_baseline"    = "Null",
    "H1f_1_LAI_only"         = "+LAI",
    "H1f_2_LAI_Hmax"         = "+Hmax",
    "H1f_3_LAI_Hmax_fCover"  = "+fCover",
    "H1f_4_Full_real"        = "+LAD profile"
  )
  x_labels <- ifelse(forward_order %in% names(short_labels),
                     short_labels[forward_order], forward_order)

  df_plot <- df_scores %>%
    filter(scenario %in% forward_order) %>%
    mutate(scenario = factor(scenario, levels = forward_order),
           step     = as.integer(scenario))
  ggplot(df_plot, aes(x = step, y = rmse)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3, colour = "#31688e") +
    geom_text(aes(label = sprintf("%.2f°C", rmse)),
              vjust = -0.9, size = 3.5, colour = "grey30") +
    scale_x_continuous(breaks = seq_along(forward_order), labels = x_labels) +
    labs(title    = "H1 — Progressive inclusion of structural descriptors",
         subtitle = "Reference = Full-real scenario (all real LiDAR descriptors)",
         x = NULL, y = "RMSE vs full-real reference (°C)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = margin(t = 5, r = 10, b = 10, l = 5, unit = "pt"))
}

#' Scatter plots ΔTmax scénario vs référence — un panneau par étape H1 forward.
#'
#' Demandé lors de la réunion du 17/04 : produire un scatter par paramètre
#' (pas seulement pour le profil complet). Chaque panneau montre la relation
#' entre le ΔTmax du scénario et celui de la référence full-real, avec R²,
#' RMSE et biais annotés.
#'
#' @param df_all           Dataframe empilé (x, y, date, Delta_Tmax, scenario).
#' @param ref_scenario_name Nom du scénario de référence (ex. "REF_all_real").
#' @param forward_names    Vecteur ordonné des noms de scénarios H1 forward à
#'   inclure (sans ref ni Null_baseline). Si NULL, tous sauf la référence.
#' @param out_path         Chemin PNG de sortie.
#' @return Invisible ggplot.
plot_h1_scenario_scatters <- function(df_all, ref_scenario_name,
                                       forward_names = NULL,
                                       out_path = "outputs/figures/04_h1_scatters_vs_ref.png") {
  label_map <- c(
    H1f_1_LAI_only        = "+ LAI",
    H1f_2_LAI_Hmax        = "+ Hmax",
    H1f_3_LAI_Hmax_fCover = "+ fCover",
    H1f_4_Full_real       = "+ Profil LAD"
  )

  df_ref <- df_all %>%
    filter(scenario == ref_scenario_name) %>%
    dplyr::select(x, y, date, Delta_ref = Delta_Tmax)

  sc_names <- if (!is.null(forward_names)) forward_names else
    setdiff(unique(df_all$scenario), c(ref_scenario_name, "H1f_0_Null_baseline"))

  df_joined <- df_all %>%
    filter(scenario %in% sc_names) %>%
    inner_join(df_ref, by = c("x", "y", "date")) %>%
    mutate(
      scenario_lab = dplyr::recode(scenario, !!!label_map, .default = scenario),
      scenario_lab = factor(scenario_lab, levels = label_map[label_map %in% unique(scenario_lab)])
    )

  lims <- range(c(df_joined$Delta_Tmax, df_joined$Delta_ref), na.rm = TRUE)
  x_ann <- lims[1] + 0.05 * diff(lims)
  y_ann <- lims[2] - 0.05 * diff(lims)

  df_stats <- df_joined %>%
    group_by(scenario_lab) %>%
    summarise(
      r2   = cor(Delta_Tmax, Delta_ref, use = "complete.obs")^2,
      rmse = sqrt(mean((Delta_Tmax - Delta_ref)^2, na.rm = TRUE)),
      bias = mean(Delta_Tmax - Delta_ref, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("R²=%.2f\nRMSE=%.2f°C\nBias=%+.2f°C", r2, rmse, bias),
      x_pos = x_ann, y_pos = y_ann
    )

  p <- ggplot(df_joined, aes(x = Delta_ref, y = Delta_Tmax)) +
    geom_hex(bins = 60) +
    scale_fill_gradient(low = "grey85", high = "midnightblue",
                        trans = "log10", name = "Count\n(log)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = FALSE, colour = "#d73027", linewidth = 1.1) +
    geom_text(
      data = df_stats,
      aes(x = x_pos, y = y_pos, label = label),
      hjust = 0, vjust = 1, fontface = "bold", size = 3.2,
      inherit.aes = FALSE
    ) +
    coord_fixed(xlim = lims, ylim = lims) +
    facet_wrap(~ scenario_lab, nrow = 1) +
    labs(
      title    = "H1 — ΔTmax scenario vs reference (full-real), by parameter",
      subtitle = "Progressive addition: LAI → Hmax → fCover → complete LAD profile",
      x        = expression("Reference (full-real)" ~ Delta * T[max] ~ (degree*C)),
      y        = expression("Scenario" ~ Delta * T[max] ~ (degree*C))
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position  = "right",
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold"))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 16, height = 6)
  cat(sprintf("  [h1_scatters] PNG : %s\n", out_path))
  invisible(p)
}

#' Comparaison Forward vs LOO pour valider la robustesse du ranking des variables.
#'
#' Si l'ordre des barres Forward et LOO est cohérent (LAI > Hmax > fCover > profil
#' LAD ou inverse), le ranking n'est pas un artefact de l'ordre d'insertion.
#' Les deux scénarios LOO pour le profil LAD (uniform vs mean-shape) illustrent
#' la sensibilité au choix de la baseline.
#'
#' @param df_scores     Dataframe issu de score_scenarios_vs_reference().
#' @param forward_order Vecteur ordonné des noms de scénarios H1 forward.
#' @param out_path      Chemin PNG de sortie.
#' @return Invisible ggplot.
compare_h1_rankings <- function(df_scores,
                                 forward_order,
                                 out_path = "outputs/figures/02_h1_ranking.png") {

  fwd_labels <- c(
    H1f_1_LAI_only        = "LAI",
    H1f_2_LAI_Hmax        = "Hmax",
    H1f_3_LAI_Hmax_fCover = "fCover",
    H1f_4_Full_real       = "Profil LAD"
  )
  loo_labels <- c(
    H1l_dropLAI_mean      = "LAI",
    H1l_dropHmax_mean     = "Hmax",
    H1l_dropfCover_mean   = "fCover",
    H1l_dropLAD_uniform   = "Profil LAD",
    H1l_dropLAD_meanShape = "Profil LAD"
  )
  loo_approach <- c(
    H1l_dropLAI_mean      = "LOO",
    H1l_dropHmax_mean     = "LOO",
    H1l_dropfCover_mean   = "LOO",
    H1l_dropLAD_uniform   = "LOO (base. uniforme)",
    H1l_dropLAD_meanShape = "LOO (base. moy-forme)"
  )
  var_order <- c("LAI", "Hmax", "fCover", "Profil LAD")
  metrics   <- c("rmse", "mae", "mean_diff")

  fwd_scs <- forward_order[forward_order %in% df_scores$scenario]
  df_fwd_scores <- df_scores %>%
    filter(scenario %in% fwd_scs) %>%
    mutate(step = match(scenario, forward_order)) %>%
    arrange(step)

  df_fwd_long <- purrr::map_dfr(metrics, function(m) {
    vals  <- setNames(df_fwd_scores[[m]], df_fwd_scores$scenario)
    steps <- fwd_scs
    purrr::map_dfr(seq_along(steps)[-1], function(i) {
      tibble(variable = fwd_labels[steps[i]], approach = "Forward",
             metric = m, delta = vals[steps[i - 1]] - vals[steps[i]])
    })
  }) %>% filter(!is.na(variable))

  loo_scs_present <- intersect(names(loo_labels), df_scores$scenario)
  df_loo_long     <- tibble()
  if (length(loo_scs_present) > 0) {
    full_row <- df_scores %>% filter(scenario == "H1f_4_Full_real")
    df_loo_long <- purrr::map_dfr(metrics, function(m) {
      full_val <- if (nrow(full_row) > 0) full_row[[m]] else 0
      df_scores %>%
        filter(scenario %in% loo_scs_present) %>%
        mutate(variable = loo_labels[scenario],
               approach = loo_approach[scenario],
               metric   = m,
               delta    = !!sym(m) - full_val) %>%
        dplyr::select(variable, approach, metric, delta)
    })
  }

  approach_levels <- c("Forward", "LOO", "LOO (base. uniforme)", "LOO (base. moy-forme)")
  approach_colors <- c(
    "Forward"              = "#31688e",
    "LOO"                  = "#d8576b",
    "LOO (base. uniforme)" = "#f0a500",
    "LOO (base. moy-forme)"= "#7b2d8b"
  )
  metric_labels <- c(
    rmse      = "ΔRMSE (°C)",
    mae       = "ΔMAE (°C)",
    mean_diff = "ΔBias (°C, signé)"
  )

  df_plot <- bind_rows(df_fwd_long, df_loo_long) %>%
    mutate(
      variable   = factor(variable, levels = var_order),
      approach   = factor(approach, levels = approach_levels),
      metric_lab = factor(metric_labels[metric], levels = metric_labels)
    )

  p <- ggplot(df_plot, aes(x = variable, y = delta, fill = approach)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_manual(values = approach_colors, name = NULL, drop = TRUE) +
    facet_wrap(~ metric_lab, scales = "free_y", ncol = 1) +
    labs(
      title    = "H1 — Ranking robustness: Forward vs LOO",
      subtitle = "Forward/LOO agreement = order-independent ranking | 2 LOO bars for LAD profile = sensitivity to baseline",
      x        = NULL,
      y        = "Error reduction (°C, positive = improvement)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 20, hjust = 1),
      strip.text      = element_text(face = "bold"),
      plot.title      = element_text(face = "bold"),
      legend.position = "top"
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 9, height = 12)
  cat(sprintf("  [h1_rankings] PNG : %s\n", out_path))
  invisible(p)
}

#' Scatter 2x2 : ΔTmax journalier LiDAR-réel (y) vs scénario avec une seule
#' variable dégradée (x). Chaque panneau = une variable remplacée par sa valeur
#' absente (LAI moyen, Hmax moyen, fCover=1, LAD uniforme).
#'
#' @param df_all  Dataframe empilé (x, y, date, Delta_Tmax, scenario).
#' @param ref_scenario_name Scénario de référence "full real".
#' @param out_path Chemin PNG de sortie.
#' @return Invisible ggplot.
plot_loo_scatters_2x2 <- function(df_all,
                                   ref_scenario_name = "H1f_4_Full_real",
                                   out_path = "outputs/figures/03_h1_loo_scatters_2x2.png") {
  loo_map <- c(
    H1l_dropLAI_mean    = "LAI → moyen",
    H1l_dropHmax_mean   = "Hmax → moyen",
    H1l_dropfCover_mean = "fCover → 1",
    H1l_dropLAD_uniform = "LAD → uniforme"
  )
  panel_order <- unname(loo_map)

  df_ref <- df_all %>%
    filter(scenario == ref_scenario_name) %>%
    dplyr::select(x, y, date, Delta_real = Delta_Tmax)

  loo_scs <- intersect(names(loo_map), unique(df_all$scenario))
  if (length(loo_scs) == 0) {
    warning("[plot_loo_scatters_2x2] aucun scénario LOO trouvé dans df_all")
    return(invisible(NULL))
  }

  df_joined <- df_all %>%
    filter(scenario %in% loo_scs) %>%
    inner_join(df_ref, by = c("x", "y", "date")) %>%
    mutate(panel = factor(loo_map[scenario], levels = panel_order))

  lims <- range(c(df_joined$Delta_Tmax, df_joined$Delta_real), na.rm = TRUE)
  x_ann <- lims[1] + 0.05 * diff(lims)
  y_ann <- lims[2] - 0.05 * diff(lims)

  df_stats <- df_joined %>%
    group_by(panel) %>%
    summarise(
      r2   = cor(Delta_Tmax, Delta_real, use = "complete.obs")^2,
      rmse = sqrt(mean((Delta_Tmax - Delta_real)^2, na.rm = TRUE)),
      bias = mean(Delta_Tmax - Delta_real, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("R²=%.2f\nRMSE=%.2f°C\nBias=%+.2f°C",
                      r2, rmse, bias),
      x_pos = x_ann, y_pos = y_ann
    )

  p <- ggplot(df_joined, aes(x = Delta_Tmax, y = Delta_real)) +
    geom_hex(bins = 50) +
    scale_fill_gradient(low = "grey88", high = "midnightblue",
                        trans = "log10", name = "N (log)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.7) +
    geom_smooth(method = "lm", se = FALSE, colour = "#d73027", linewidth = 1) +
    geom_text(data = df_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = 0, vjust = 1, size = 3, fontface = "bold",
              inherit.aes = FALSE) +
    coord_fixed(xlim = lims, ylim = lims) +
    facet_wrap(~ panel, nrow = 2, ncol = 2) +
    labs(
      title    = "Impact of a structural variable — daily ΔTmax",
      subtitle = paste0("x = scenario with degraded variable | y = real LiDAR (Full-real)\n",
                        "Cross-read Fig. 4: LAI-only (Fwd, R²≈0.70) > all-but-LAI (LOO, R²≈0.57)",
                        " — LAI spatial variability dominates alone"),
      caption  = paste0("Negative bias on fCover (LOO): replacing real fCover with mean fCover",
                        " artificially closes the canopy on the most open plots",
                        " — MuSICA then simulates more cooling, producing a systematic negative bias."),
      x        = expression(Delta * T[max] ~ "degraded scenario (°C)"),
      y        = expression(Delta * T[max] ~ "real LiDAR (°C)")
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"),
          plot.caption    = element_text(size = 8, colour = "grey35", hjust = 0),
          legend.position = "right")

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 10, height = 10)
  cat(sprintf("  [loo_scatters_2x2] PNG : %s\n", out_path))
  invisible(p)
}

# ---- Shapley exact attribution (2^4 lattice) ---------------------------------

#' Extract metric values for the 16 coalitions from df_scores.
#'
#' The coalition→scenario mapping is generated dynamically when
#' `factorial_scenarios` is provided (preferred — avoids silent breakage if
#' scenario naming conventions change). Falls back to the hardcoded map when
#' NULL, but validates all 16 keys against df_scores$scenario.
#'
#' Metric handling :
#'   - `"rmse"` / `"mae"` : loss-type (lower = better). Retournés tels quels.
#'   - `"r2"` : score-type (higher = better). Renvoyés en forme "loss" (1 − R²)
#'              pour que la machinerie Shapley (v(S) = v_empty − v(S))
#'              s'applique sans changement. v(N) sera alors R²_full − R²_null.
#'
#' @param df_scores           Tibble from score_scenarios_vs_reference.
#' @param factorial_scenarios Optional named list from scenarios_h1_factorial().
#' @param metric              "rmse" (défaut), "mae" ou "r2".
#' @return Named numeric vector (16 valeurs) avec attributs :
#'   `metric` (nom de la métrique) et `raw_values` (valeurs originales avant
#'   conversion loss pour le cas R²).
extract_shapley_coalitions <- function(df_scores, factorial_scenarios = NULL,
                                        metric = c("rmse", "mae", "r2")) {
  metric <- match.arg(metric)

  if (!is.null(factorial_scenarios)) {
    # Build map from the scenario objects: each has $grid_row with
    # columns LAI/Hmax/fCover/LAD encoding real vs degraded.
    # setNames(values, names): values = scenario names, names = coalition keys.
    coalition_to_scenario <- setNames(
      vapply(factorial_scenarios, function(sc) sc$name, character(1)),
      vapply(factorial_scenarios, function(sc) {
        g <- sc$grid_row
        paste0(as.integer(g$LAI    == "real"),
               as.integer(g$Hmax   == "real"),
               as.integer(g$fCover == "real"),
               as.integer(g$LAD    == "real"))
      }, character(1))
    )
  } else {
    coalition_to_scenario <- c(
      "0000" = "H1F_m_m_a_u", "1000" = "H1F_r_m_a_u",
      "0100" = "H1F_m_r_a_u", "0010" = "H1F_m_m_r_u",
      "0001" = "H1F_m_m_a_r", "1100" = "H1F_r_r_a_u",
      "1010" = "H1F_r_m_r_u", "1001" = "H1F_r_m_a_r",
      "0110" = "H1F_m_r_r_u", "0101" = "H1F_m_r_a_r",
      "0011" = "H1F_m_m_r_r", "1110" = "H1F_r_r_r_u",
      "1101" = "H1F_r_r_a_r", "1011" = "H1F_r_m_r_r",
      "0111" = "H1F_m_r_r_r", "1111" = "H1F_r_r_r_r"
    )
    # Validate all 16 scenario names exist in df_scores
    missing_sc <- setdiff(coalition_to_scenario, df_scores$scenario)
    if (length(missing_sc) > 0)
      stop(sprintf("Treillis Shapley — scénarios absents de df_scores : %s",
                   paste(missing_sc, collapse = ", ")))
    unexpected_sc <- setdiff(df_scores$scenario[grepl("^H1F_", df_scores$scenario)],
                             coalition_to_scenario)
    if (length(unexpected_sc) > 0)
      warning(sprintf("Scénarios H1F_* dans df_scores non inclus dans la carte Shapley : %s",
                      paste(unexpected_sc, collapse = ", ")))
  }

  if (!metric %in% names(df_scores))
    stop(sprintf("Colonne '%s' absente de df_scores (cols dispo : %s)",
                 metric, paste(names(df_scores), collapse = ", ")))

  raw_values <- sapply(coalition_to_scenario, function(sc) {
    val <- df_scores[[metric]][df_scores$scenario == sc]
    if (length(val) == 0) {
      warning(sprintf("Scenario manquant pour treillis Shapley : %s", sc))
      return(NA_real_)
    }
    val[1]
  })
  names(raw_values) <- names(coalition_to_scenario)
  if (any(is.na(raw_values)))
    stop("Treillis Shapley incomplet : certains scenarios factorial manquent.")

  # Convertir R² (score) en forme loss pour reutiliser compute_shapley_exact
  loss_values <- if (metric == "r2") 1 - raw_values else raw_values
  names(loss_values) <- names(raw_values)
  attr(loss_values, "metric")     <- metric
  attr(loss_values, "raw_values") <- raw_values
  loss_values
}

#' Compute exact Shapley values from the 2^4 factorial loss lattice.
#'
#' Generic en metric : accepte n'importe quelle vecteur "loss" (lower=better)
#' via la convention v(S) = v_empty − loss(S). Pour R², l'appelant
#' (extract_shapley_coalitions) fournit déjà 1−R² en entrée.
#'
#' @param loss_values Named numeric vector (16 coalitions) from
#'   extract_shapley_coalitions. Peut porter un attribut `metric` ("rmse"|"r2"|"mae").
#' @return List : $shapley, $rmse_values (legacy alias), $v_empty,
#'   $total_effect, $total_sum, $check_diff, $metric, $units.
compute_shapley_exact <- function(loss_values) {
  metric <- attr(loss_values, "metric")
  if (is.null(metric)) metric <- "rmse"
  units  <- switch(metric, "rmse" = "°C", "mae" = "°C", "r2" = "", "")

  v_empty <- as.numeric(loss_values["0000"])
  v <- function(b) v_empty - as.numeric(loss_values[b])

  make_bin <- function(bits_on, n = 4) {
    b <- rep("0", n); b[bits_on] <- "1"; paste(b, collapse = "")
  }
  variables <- c("LAI", "Hmax", "fCover", "LAD")

  shapley_vals <- sapply(seq_along(variables), function(var_bit) {
    others <- setdiff(seq_len(4), var_bit)
    total  <- 0
    for (s in 0:3) {
      subsets <- if (s == 0) list(integer(0)) else combn(others, s, simplify = FALSE)
      w <- factorial(s) * factorial(4 - s - 1) / factorial(4)
      for (S in subsets)
        total <- total + w * (v(make_bin(c(S, var_bit))) - v(make_bin(S)))
    }
    total
  })
  names(shapley_vals) <- variables

  total_shapley <- sum(shapley_vals)
  total_effect  <- v("1111")

  unit_suffix <- if (nzchar(units)) sprintf(" %s", units) else ""
  cat(sprintf("── Valeurs Shapley exactes (metric=%s) ──\n", metric))
  for (i in seq_along(variables))
    cat(sprintf("  phi(%-6s) = %+.4f%s  (%+.1f%% de v(N))\n",
                variables[i], shapley_vals[i], unit_suffix,
                100 * shapley_vals[i] / total_effect))
  cat(sprintf("\n  Somme phi = %+.4f  |  v(N) = %+.4f  |  ecart = %.2e\n\n",
              total_shapley, total_effect, abs(total_shapley - total_effect)))

  list(shapley = shapley_vals, rmse_values = loss_values,   # legacy name
       loss_values = loss_values,
       v_empty = v_empty, total_effect = total_effect,
       total_sum = total_shapley,
       check_diff = abs(total_shapley - total_effect),
       metric = metric, units = units)
}

#' Pearson correlation matrix of structural LiDAR metrics on the cLHS sample.
#'
#' Documents the input correlation structure that motivates Shapley-Owen over
#' classical Sobol indices (which assume independent inputs).
#'
#' @param df_sample cLHS sample dataframe with structural columns.
#' @param vars      Column names to include (intersected with names(df_sample)).
#' @param out_path  PNG output path.
#' @return Invisible list: $cor_matrix, $n, $plot.
plot_structural_correlations <- function(df_sample,
                                          vars = c("LAI", "Hmax", "fCover", "FPC1"),
                                          out_path = "outputs/figures/annex/A22_structural_correlations.png") {
  vars_present <- intersect(vars, names(df_sample))
  if (length(vars_present) < 2) {
    warning("[plot_structural_correlations] fewer than 2 vars found in df_sample")
    return(invisible(NULL))
  }

  df_mat <- df_sample %>% dplyr::select(all_of(vars_present)) %>% drop_na()
  n_obs  <- nrow(df_mat)
  cor_mat <- cor(df_mat, method = "pearson")

  cat(sprintf("\n── Pearson correlation matrix (cLHS sample, n=%d) ──\n", n_obs))
  print(round(cor_mat, 3))

  cor_long <- as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("var1") %>%
    tidyr::pivot_longer(-var1, names_to = "var2", values_to = "r") %>%
    mutate(
      var1  = factor(var1, levels = rev(vars_present)),
      var2  = factor(var2, levels = vars_present),
      label = sprintf("%.2f", r)
    )

  p <- ggplot(cor_long, aes(x = var2, y = var1, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label, colour = abs(r) > 0.55),
              size = 4.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#d73027",
      midpoint = 0, limits = c(-1, 1), name = "Pearson r"
    ) +
    scale_colour_manual(values = c("FALSE" = "grey20", "TRUE" = "white"),
                        guide  = "none") +
    scale_x_discrete(position = "top") +
    coord_fixed() +
    labs(
      title    = "Structural correlation matrix — cLHS sample",
      subtitle = sprintf(
        paste0("n = %d plots | Pearson r | High correlations (e.g. LAI-fCover)",
               " justify Shapley-Owen over classical Sobol (which assumes independence)"),
        n_obs
      ),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(colour = "grey30", size = 9),
      axis.text       = element_text(face = "bold", size = 11),
      panel.grid      = element_blank(),
      legend.position = "right"
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 7, height = 6)
  cat(sprintf("  [structural_corr] PNG : %s\n", out_path))
  invisible(list(cor_matrix = cor_mat, n = n_obs, plot = p))
}

# ---- ANOVA bootstrap (stratified resampling of plots) ------------------------

#' Bootstrap 95% CI on SS_princ and SS_tot from the 2^4 ANOVA factorial design.
#'
#' Resamples plots (with replacement, stratified by cluster) and recomputes
#' the 16-coalition RMSE table for each iteration. Each iteration refits the
#' saturated 2^4 lm and extracts SS_princ / SS_tot per variable, plus the
#' scalar sigma_princ (sum of main-effect SS fractions).
#'
#' Exact when all plots have equal date counts (verified inside; stops otherwise).
#' RMSE_boot(S) = sqrt( mean_i_boot( MSE_i(S) ) )
#'
#' @param df_all        Stacked scenario results (x, y, date, Delta_Tmax, scenario).
#' @param coalition_map Named character vector: names = coalition keys "0000"…"1111",
#'   values = corresponding scenario names in df_all.
#' @param ref_name      Reference scenario name (coalition "1111", i.e. full-real).
#' @param cluster_vec   Named integer/character vector: names = "x y" plot keys,
#'   values = cluster labels. Build as setNames(df_sample$Cluster, paste(x, y)).
#'   NULL falls back to single-cluster (unstratified) bootstrap.
#' @param n_boot        Number of bootstrap iterations (default 50).
#' @return Long tibble: boot_id, variable, SS_princ, SS_tot, sigma_princ, int_pct.
bootstrap_anova_factorial <- function(df_all, coalition_map, ref_name,
                                       cluster_vec = NULL, n_boot = 50) {
  all_scs <- unique(as.character(coalition_map))

  df_ref <- df_all %>%
    filter(scenario == ref_name) %>%
    dplyr::select(x, y, date, Delta_ref = Delta_Tmax)

  df_fac <- df_all %>%
    filter(scenario %in% all_scs) %>%
    inner_join(df_ref, by = c("x", "y", "date")) %>%
    mutate(sq_err = (Delta_Tmax - Delta_ref)^2)

  # Balance check — formula exact only when all plots share the same date count
  n_dates_check <- df_fac %>%
    filter(scenario == all_scs[1]) %>%
    count(x, y, name = "n_dates")
  if (dplyr::n_distinct(n_dates_check$n_dates) > 1)
    stop(sprintf(
      "[ANOVA bootstrap] Unbalanced plots (date counts %d–%d); RMSE formula not exact.",
      min(n_dates_check$n_dates), max(n_dates_check$n_dates)
    ))

  # Pre-aggregate to per-plot MSE per coalition (avoids re-expanding in each iteration)
  df_mse <- df_fac %>%
    mutate(plot_key = paste(x, y)) %>%
    group_by(plot_key, scenario) %>%
    summarise(mse = mean(sq_err, na.rm = TRUE), .groups = "drop")

  # Attach cluster labels and verify full match
  plot_keys_all <- unique(df_mse$plot_key)
  if (!is.null(cluster_vec)) {
    n_matched <- sum(plot_keys_all %in% names(cluster_vec))
    if (n_matched == 0)
      stop(paste0(
        "[ANOVA bootstrap] cluster_vec keys don't match any df_mse plot keys.\n",
        "  Expected: setNames(df_sample$Cluster, paste(round(x), round(y)))"
      ))
    if (n_matched < length(plot_keys_all)) {
      unmatched <- plot_keys_all[!plot_keys_all %in% names(cluster_vec)]
      stop(sprintf(
        "[ANOVA bootstrap] Partial cluster mismatch: %d/%d plots unmatched.\n  First 5: %s",
        length(unmatched), length(plot_keys_all),
        paste(head(unmatched, 5), collapse = ", ")
      ))
    }
    df_plots <- tibble(
      plot_key = plot_keys_all,
      cluster  = cluster_vec[plot_keys_all]
    )
    cat(sprintf("  [ANOVA bootstrap] %d plots matched to cluster_vec (%d clusters)\n",
                nrow(df_plots), dplyr::n_distinct(df_plots$cluster)))
  } else {
    df_plots <- tibble(plot_key = plot_keys_all, cluster = "all")
  }

  set.seed(43)  # use 43 (not 42) to avoid collision with bootstrap_shapley's seed
  var_map   <- c("LAI" = "L", "Hmax" = "H", "fCover" = "F", "LAD" = "D")
  boot_list <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    # Stratified with-replacement bootstrap: prop=1 with replace=TRUE is
    # equivalent to slice_sample(n=n()) but works across all dplyr versions.
    sampled_keys <- df_plots %>%
      group_by(cluster) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      pull(plot_key)

    # left_join preserves duplicate keys → correctly reflects bootstrap resampling
    df_boot_mse <- dplyr::left_join(
      tibble(plot_key = sampled_keys), df_mse, by = "plot_key"
    )

    rmse_tbl <- df_boot_mse %>%
      group_by(scenario) %>%
      summarise(rmse_b = sqrt(mean(mse, na.rm = TRUE)), .groups = "drop")

    rmse_map   <- setNames(rmse_tbl$rmse_b, rmse_tbl$scenario)
    rmse_b_vec <- vapply(coalition_map, function(sc) as.numeric(rmse_map[sc]), numeric(1))
    names(rmse_b_vec) <- names(coalition_map)

    if (any(is.na(rmse_b_vec))) {
      warning(sprintf("[ANOVA bootstrap] boot %d: %d NA coalitions — skipping",
                      b, sum(is.na(rmse_b_vec))))
      next
    }

    v_empty_b <- as.numeric(rmse_b_vec["0000"])
    df_b <- tibble(
      key = names(rmse_b_vec),
      v_S = v_empty_b - as.numeric(rmse_b_vec)
    ) %>%
      mutate(
        L = 2L * as.integer(substr(key, 1, 1)) - 1L,
        H = 2L * as.integer(substr(key, 2, 2)) - 1L,
        F = 2L * as.integer(substr(key, 3, 3)) - 1L,
        D = 2L * as.integer(substr(key, 4, 4)) - 1L
      )

    mod_b   <- lm(v_S ~ L * H * F * D, data = df_b)
    aov_b   <- anova(mod_b)
    ss_b    <- setNames(aov_b[["Sum Sq"]], rownames(aov_b))
    ss_b    <- ss_b[!grepl("Residuals", names(ss_b))]
    ss_tot_b <- sum(ss_b)
    term_nb  <- names(ss_b)

    df_princ_b <- purrr::imap_dfr(var_map, function(code, nm) {
      tibble(variable = nm, SS_princ = ss_b[[code]] / ss_tot_b)
    })
    df_tot_b <- purrr::imap_dfr(var_map, function(code, nm) {
      inv <- term_nb[sapply(strsplit(term_nb, ":", fixed = TRUE),
                            function(pts) code %in% pts)]
      tibble(variable = nm, SS_tot = sum(ss_b[inv]) / ss_tot_b)
    })

    sigma_b <- sum(df_princ_b$SS_princ)
    int_b   <- 100 * sum(ss_b[term_nb[stringr::str_count(term_nb, ":") >= 1]]) / ss_tot_b

    boot_list[[b]] <- tibble(
      boot_id     = b,
      variable    = df_princ_b$variable,
      SS_princ    = df_princ_b$SS_princ,
      SS_tot      = df_tot_b$SS_tot,
      sigma_princ = sigma_b,
      int_pct     = int_b
    )

    if (b %% 10 == 0 || b == n_boot)
      cat(sprintf("  [ANOVA bootstrap] boot %d/%d done\n", b, n_boot))
  }
  dplyr::bind_rows(boot_list)
}

#' ANOVA factorial decomposition from the 2^4 design (cross-check for Shapley).
#'
#' Decomposes the variance of v(S) = RMSE_null - RMSE(S) across the 16 coalitions
#' using a saturated 2^4 ANOVA with effect coding. Yields SS_principal / SS_total
#' (analogous to Sobol S1) and SS_total_i / SS_total (analogous to Sobol ST) for
#' each variable.
#'
#' NOTE: this is NOT Sobol stricto sensu. Sobol indices require the continuous
#' empirical distribution of each input (Saltelli sampling). This computation
#' binarises each input (real vs baseline) and therefore estimates factorial
#' effects, not conditional variance over the full input space. The correct label
#' is "ANOVA factorial decomposition (Sobol-like, 2-level design)".
#'
#' Despite the limitation, convergence of rankings between Shapley and the ANOVA
#' decomposition is a valid robustness argument (two independent decompositions
#' of the same design variance agree on the hierarchy).
#'
#' @param rmse_values Named numeric vector from extract_shapley_coalitions
#'   (names = "0000"..."1111", values = RMSE vs full-real reference).
#' @param shap_res    Output of compute_shapley_exact ($shapley, $total_effect).
#'   NULL to omit the Shapley comparison bars.
#' @param out_path    PNG output path.
#' @return Invisible list: $df_princ, $df_tot, $aov_table, $plot.
compute_anova_factorial <- function(rmse_values, shap_res = NULL,
                                     rho_lai_fcover = NULL,
                                     out_path = "outputs/figures/annex/A21_anova_vs_shapley.png",
                                     df_all = NULL,
                                     coalition_map = NULL,
                                     cluster_vec = NULL,
                                     n_boot = 50) {
  v_empty <- as.numeric(rmse_values["0000"])

  # Convert RMSE to improvement over null (higher = better, consistent with v(S))
  df <- tibble(
    key = names(rmse_values),
    v_S = v_empty - as.numeric(rmse_values)
  ) %>%
    mutate(
      LAI    = as.integer(substr(key, 1, 1)),
      Hmax   = as.integer(substr(key, 2, 2)),
      fCover = as.integer(substr(key, 3, 3)),
      LAD    = as.integer(substr(key, 4, 4)),
      L = 2L * LAI    - 1L,   # effect coding: baseline -> -1, real -> +1
      H = 2L * Hmax   - 1L,
      F = 2L * fCover - 1L,
      D = 2L * LAD    - 1L
    )

  # Saturated 2^4 ANOVA (16 obs, 15 effects + intercept -> 0 residual df, exact fit)
  # No p-values possible (df_residual = 0) — decomposition is purely algebraic.
  mod     <- lm(v_S ~ L * H * F * D, data = df)
  aov_tbl <- anova(mod)
  terms   <- rownames(aov_tbl)
  ss_vec  <- setNames(aov_tbl[["Sum Sq"]], terms)
  ss_vec  <- ss_vec[!grepl("Residuals", names(ss_vec))]
  ss_tot  <- sum(ss_vec)

  var_map    <- c("LAI" = "L", "Hmax" = "H", "fCover" = "F", "LAD" = "D")
  term_names <- names(ss_vec)

  # SS_principal: main effect only (analogous to Sobol S1)
  df_princ <- purrr::imap_dfr(var_map, function(code, nm) {
    tibble(variable = nm, SS_princ = ss_vec[[code]] / ss_tot)
  })

  # SS_total_i: all terms involving variable i (analogous to Sobol ST)
  # Use strsplit to compare exact term components — avoids substring false matches
  # (e.g. grepl("L", "L:H") would also match if any future code contains "L" as substring).
  df_tot <- purrr::imap_dfr(var_map, function(code, nm) {
    involves <- term_names[
      sapply(strsplit(term_names, ":", fixed = TRUE),
             function(parts) code %in% parts)
    ]
    tibble(variable = nm, SS_tot = sum(ss_vec[involves]) / ss_tot)
  })

  df_int <- tibble(term = term_names, SS = ss_vec) %>%
    mutate(
      n_vars = stringr::str_count(term, ":") + 1L,
      type   = dplyr::case_when(
        n_vars == 1 ~ "Main", n_vars == 2 ~ "2-way",
        n_vars == 3 ~ "3-way",   TRUE        ~ "4-way"
      ),
      pct = round(100 * SS / ss_tot, 2)
    )

  sigma_princ <- round(sum(df_princ$SS_princ), 3)
  # "Main" is the correct level — df_int$type is "Main"/"2-way"/"3-way"/"4-way";
  # filtering on "Principal" (old wrong value) always returns TRUE -> int_pct=100%.
  int_pct     <- round(100 * sum(df_int$SS[df_int$type != "Main"]) / ss_tot, 1)

  cat(sprintf("  [ANOVA sanity] sigma_princ = %.3f | int_pct = %.1f%% | sum = %.1f%%\n",
              sigma_princ, int_pct, sigma_princ * 100 + int_pct))

  # ---- Optional bootstrap CI on sigma_princ ------------------------------------
  df_boot_anova <- NULL
  ci_sigma_lo   <- NA_real_
  ci_sigma_hi   <- NA_real_

  if (!is.null(df_all) && !is.null(coalition_map)) {
    cat(sprintf("  [ANOVA bootstrap] Running %d iterations (stratified by cluster)...\n", n_boot))
    ref_name_b <- unname(coalition_map["1111"])
    df_boot_anova <- bootstrap_anova_factorial(
      df_all        = df_all,
      coalition_map = coalition_map,
      ref_name      = ref_name_b,
      cluster_vec   = cluster_vec,
      n_boot        = n_boot
    )
    sigma_ci <- df_boot_anova %>%
      dplyr::distinct(boot_id, sigma_princ) %>%
      dplyr::summarise(
        ci_lo = quantile(sigma_princ, 0.025, na.rm = TRUE),
        ci_hi = quantile(sigma_princ, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
    ci_sigma_lo <- sigma_ci$ci_lo
    ci_sigma_hi <- sigma_ci$ci_hi
    cat(sprintf("  [ANOVA bootstrap] sigma_princ = %.3f [%.3f, %.3f] (CI width = %.3f)\n",
                sigma_princ, ci_sigma_lo, ci_sigma_hi, ci_sigma_hi - ci_sigma_lo))
    dir.create("outputs/audit", recursive = TRUE, showWarnings = FALSE)
    write.csv(df_boot_anova, "outputs/audit/anova_bootstrap.csv", row.names = FALSE)
    cat("  [ANOVA bootstrap] CSV: outputs/audit/anova_bootstrap.csv\n")
  }

  rho_str <- if (!is.null(rho_lai_fcover))
    sprintf("rho(LAI,fCover)=%.2f — Shapley robust to correlations", rho_lai_fcover)
  else
    "Shapley robust to correlations (provide rho_lai_fcover to quantify)"

  sigma_ci_str <- if (!is.na(ci_sigma_lo))
    sprintf("Sigma-SS_princ = %.3f [%.3f, %.3f] 95%%CI (n_boot=%d) | interactions = %.1f%%",
            sigma_princ, ci_sigma_lo, ci_sigma_hi, n_boot, int_pct)
  else
    sprintf("Sigma-SS_princ = %.3f (interactions = %.1f%% of SS_total)", sigma_princ, int_pct)

  subtitle_str <- sprintf(
    "Same 16 coalitions (2^4 factorial) | Exact ANOVA by main effects + interactions\n%s | %s",
    sigma_ci_str, rho_str
  )

  cat("\n── ANOVA factorial decomposition (2^4, exact) ──\n")
  cat("  NOTE: 2 niveaux / facteur -> fraction de SS, pas indice Sobol strict\n")
  cat(sprintf("  Total SS = %.6f\n\n", ss_tot))
  for (i in seq_len(nrow(df_int)))
    cat(sprintf("  %-20s [%s]  SS=%.6f  (%5.2f%%)\n",
                df_int$term[i], df_int$type[i], df_int$SS[i], df_int$pct[i]))
  cat(sprintf("\n  SS_princ: %s  Sigma=%.3f\n",
              paste(sprintf("%s=%.3f", df_princ$variable, df_princ$SS_princ),
                    collapse=" | "),
              sigma_princ))
  cat(sprintf("  SS_tot_i: %s\n\n",
              paste(sprintf("%s=%.3f", df_tot$variable, df_tot$SS_tot),
                    collapse=" | ")))

  # ---- Comparison plot -------------------------------------------------------
  df_plot <- df_princ %>%
    rename(value = SS_princ) %>%
    mutate(index = "SS main / SS_total") %>%
    bind_rows(
      df_tot %>% rename(value = SS_tot) %>%
        mutate(index = "SS total (incl. interactions) / SS_total")
    )

  if (!is.null(shap_res)) {
    phi    <- shap_res$shapley
    v_N    <- shap_res$total_effect
    df_phi <- tibble(
      variable = names(phi),
      value    = as.numeric(phi) / v_N,   # phi_i / v(N), Σ = 1 by efficiency axiom
      index    = "Shapley phi_i / v(N)"
    ) %>% filter(variable %in% names(var_map))
    df_plot <- bind_rows(df_phi, df_plot)
  }

  var_order <- c("LAI", "fCover", "Hmax", "LAD")
  index_levels <- c("Shapley phi_i / v(N)",
                    "SS main / SS_total",
                    "SS total (incl. interactions) / SS_total")
  df_plot <- df_plot %>%
    mutate(
      variable = factor(variable, levels = var_order),
      index    = factor(index, levels = intersect(index_levels, unique(index)))
    )

  # Attach bootstrap CI — NA for Shapley rows (kept so position_dodge is uniform)
  df_plot$ci_lo <- NA_real_
  df_plot$ci_hi <- NA_real_

  if (!is.null(df_boot_anova)) {
    ci_var <- df_boot_anova %>%
      dplyr::group_by(variable) %>%
      dplyr::summarise(
        ci_lo_princ = quantile(SS_princ, 0.025, na.rm = TRUE),
        ci_hi_princ = quantile(SS_princ, 0.975, na.rm = TRUE),
        ci_lo_tot   = quantile(SS_tot,   0.025, na.rm = TRUE),
        ci_hi_tot   = quantile(SS_tot,   0.975, na.rm = TRUE),
        .groups     = "drop"
      )
    for (i in seq_len(nrow(df_plot))) {
      vr  <- as.character(df_plot$variable[i])
      idx <- as.character(df_plot$index[i])
      ci_row <- ci_var[ci_var$variable == vr, ]
      if (nrow(ci_row) == 0) next
      if (idx == "SS main / SS_total") {
        df_plot$ci_lo[i] <- ci_row$ci_lo_princ
        df_plot$ci_hi[i] <- ci_row$ci_hi_princ
      } else if (idx == "SS total (incl. interactions) / SS_total") {
        df_plot$ci_lo[i] <- ci_row$ci_lo_tot
        df_plot$ci_hi[i] <- ci_row$ci_hi_tot
      }
    }
  }

  p <- ggplot(df_plot, aes(x = variable, y = value, fill = index)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             colour = "grey25", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = 0.75),
                  width = 0.25, linewidth = 0.7, na.rm = TRUE) +
    scale_fill_manual(
      values = c(
        "Shapley phi_i / v(N)"                       = "#31688e",
        "SS main / SS_total"                         = "#35b779",
        "SS total (incl. interactions) / SS_total"   = "#fca50a"
      ),
      name = NULL
    ) +
    labs(
      title    = "H1 — Triangulation: Shapley vs 2^4 factorial ANOVA",
      subtitle = subtitle_str,
      caption  = paste0(
        "2^4 factorial ANOVA (2 levels/factor: real vs baseline) — this is NOT Sobol stricto sensu\n",
        "(Sobol requires sampling the continuous empirical distribution; here inputs are binarised).\n",
        "Convergence of Shapley / ANOVA rankings provides method-independent robustness evidence."
      ),
      x = NULL,
      y = "Variance fraction (SS_i / SS_total or phi_i / v(N))"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(colour = "grey30", size = 9),
      plot.caption    = element_text(size = 7.5, colour = "grey35", hjust = 0,
                                     lineheight = 1.3),
      legend.position = "bottom"
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 12, height = 7)
  cat(sprintf("  [anova_factorial] PNG : %s\n", out_path))

  write.csv(df_int, sub("\\.png$", "_anova_table.csv", out_path), row.names = FALSE)

  invisible(list(df_princ = df_princ, df_tot = df_tot, aov_table = df_int,
                 plot = p, df_boot = df_boot_anova))
}

#' Three-panel figure: Shapley bars, Forward/LOO/Shapley comparison, coalition lattice.
#'
#' @param shap_res  Result of compute_shapley_exact.
#' @param df_scores Tibble from score_scenarios_vs_reference.
#' @param out_path  PNG output path.
#' @return Invisible patchwork ggplot.
plot_shapley_attribution <- function(shap_res, df_scores,
                                     out_path = "outputs/figures/02_shapley_attribution.png",
                                     panels = "ABC") {
  variables <- c("LAI", "Hmax", "fCover", "LAD")
  var_order <- c("LAI", "fCover", "Hmax", "LAD")
  rv        <- shap_res$rmse_values
  metric    <- if (is.null(shap_res$metric)) "rmse" else shap_res$metric
  units     <- if (is.null(shap_res$units))  "°C"  else shap_res$units
  metric_lbl <- toupper(metric)
  unit_paren <- if (nzchar(units)) sprintf(" (%s)", units) else ""

  # as.numeric() strips inner names to avoid "LAI.0000" concatenation bug
  forward_gain <- c(
    LAI    = as.numeric(rv["0000"] - rv["1000"]),
    Hmax   = as.numeric(rv["1000"] - rv["1100"]),
    fCover = as.numeric(rv["1100"] - rv["1110"]),
    LAD    = as.numeric(rv["1110"] - rv["1111"])
  )
  loo_loss <- c(
    LAI    = as.numeric(rv["0111"] - rv["1111"]),
    Hmax   = as.numeric(rv["1011"] - rv["1111"]),
    fCover = as.numeric(rv["1101"] - rv["1111"]),
    LAD    = as.numeric(rv["1110"] - rv["1111"])
  )

  # Panel A: Shapley bars
  df_A <- data.frame(
    variable = factor(variables, levels = var_order),
    phi      = as.numeric(shap_res$shapley),
    label    = sprintf("%+.3f\n(%+.0f%%)",
                       shap_res$shapley,
                       100 * shap_res$shapley / shap_res$total_effect)
  )
  pal_vars <- c(LAI = "#440154", fCover = "#31688e", Hmax = "#35b779", LAD = "#fde725")
  pA <- ggplot(df_A, aes(x = variable, y = phi, fill = variable)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_text(aes(label = label, vjust = ifelse(phi >= 0, -0.3, 1.3)),
              size = 4, fontface = "bold") +
    scale_fill_manual(values = pal_vars) +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.20))) +
    labs(title    = sprintf("A. Shapley attribution of canopy structural drivers — metric: %s",
                              metric_lbl),
         subtitle = sprintf("v(N) = %.3f%s | Σφ = %.4f",
                            shap_res$total_effect, unit_paren, shap_res$total_sum),
         x = NULL, y = sprintf("phi%s", unit_paren)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", size = 12))

  # Panel B: Forward / LOO / Shapley comparison
  df_B <- data.frame(
    variable = rep(variables, 3),
    method   = rep(c("Forward", "LOO", "Shapley"), each = 4),
    value    = c(as.numeric(forward_gain[variables]),
                 as.numeric(loo_loss[variables]),
                 as.numeric(shap_res$shapley))
  ) %>% mutate(variable = factor(variable, levels = var_order),
               method   = factor(method, levels = c("Forward", "LOO", "Shapley")))

  pB <- ggplot(df_B, aes(x = variable, y = value, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.7,
             colour = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = c(Forward = "#31688e", LOO = "#d8576b", Shapley = "#5ec962"),
                      name = NULL) +
    labs(title    = "B. Method comparison : Forward / LOO / Shapley",
         subtitle = "Shapley averages over all permutations — resolves Forward/LOO divergence",
         x = NULL, y = sprintf("Attribution%s", unit_paren)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          axis.text.x = element_text(face = "bold", size = 12))

  # Panel C: coalition lattice
  df_C <- data.frame(
    coalition = names(rv), rmse = as.numeric(rv),
    v_S = shap_res$v_empty - as.numeric(rv)
  ) %>% mutate(
    size  = factor(sapply(coalition, function(s) sum(as.integer(strsplit(s,"")[[1]]))),
                   levels = 0:4, labels = paste0("|S|=", 0:4)),
    label = sapply(coalition, function(s) {
      bits <- as.integer(strsplit(s,"")[[1]])
      lbl  <- paste(c("L","H","F","D")[bits == 1], collapse="")
      if (lbl == "") "∅" else lbl
    })
  )
  pC <- ggplot(df_C, aes(x = size, y = v_S)) +
    geom_jitter(aes(colour = size), width = 0.15, size = 3.5, alpha = 0.85) +
    geom_text(aes(label = label), size = 2.8, vjust = -1.1,
              colour = "grey30", fontface = "bold") +
    geom_hline(yintercept = shap_res$total_effect, linetype = "dashed",
               colour = "red", linewidth = 0.6) +
    annotate("text", x = 0.6, y = shap_res$total_effect,
             label = sprintf("v(N)=%.3f%s", shap_res$total_effect, unit_paren),
             hjust = 0, vjust = -0.5, colour = "red", fontface = "bold", size = 3.5) +
    scale_colour_viridis_d(option = "plasma", end = 0.85) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(title    = "C. The 2^4 = 16 coalitions of the Shapley lattice",
         subtitle = sprintf("v(S) = %s(∅) − %s(S) | L=LAI, H=Hmax, F=fCover, D=LAD",
                             metric_lbl, metric_lbl),
         x = "Coalition size |S|", y = sprintf("v(S)%s", unit_paren)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))

  ann <- patchwork::plot_annotation(
    title    = sprintf("Shapley attribution of canopy structural drivers on ΔTmax (metric: %s)",
                        metric_lbl),
    subtitle = "Blois oak forest | MuSICA simulations | Summer 2021-2022",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(colour = "grey30"))
  )
  if (panels == "A") {
    fig <- pA + ann
    fig_w <- 7; fig_h <- 6
  } else if (panels == "AB") {
    fig <- (pA | pB) + ann
    fig_w <- 14; fig_h <- 6
  } else if (panels == "C") {
    fig <- pC + ann
    fig_w <- 8; fig_h <- 6
  } else {
    fig <- (pA + pB) / pC + ann +
      patchwork::plot_layout(heights = c(1.2, 1))
    fig_w <- 14; fig_h <- 11
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(fig, out_path, width = fig_w, height = fig_h)

  csv_path <- sub("\\.png$", ".csv", out_path)
  write.csv(data.frame(
    variable     = variables,
    shapley      = as.numeric(shap_res$shapley),
    shapley_pct  = 100 * shap_res$shapley / shap_res$total_effect,
    forward_gain = as.numeric(forward_gain[variables]),
    loo_loss     = as.numeric(loo_loss[variables])
  ), csv_path, row.names = FALSE)
  cat(sprintf("  [Shapley] PNG : %s\n  [Shapley] CSV : %s\n", out_path, csv_path))
  invisible(fig)
}

# ---- Shapley a 3 joueurs : {LAI, fCover, Structure3D=(Hmax,LAD)} -------------
# Motivation : dans .lad_rescale(), le profil LAD est co-defini avec Hmax (sa
# forme normalisee z/Hmax_real est etiree au Hmax_target). Les traiter comme
# 2 joueurs independants dans le 2^4 cree des coalitions structurellement
# incoherentes (e.g. LAD reel de Hmax=39m force sur Hmax=21m). Cette analyse
# 2^3 = 8 coalitions respecte le couplage : Hmax et LAD ne varient qu'ensemble.

#' Extract 8 coalitions for the 3-player Shapley {LAI, fCover, Structure3D}.
#'
#' Reutilise 8 des 16 simulations factorielles existantes — pas de nouveau run.
#' Structure3D prend 2 valeurs :
#'   - baseline = (Hmax=mean, LAD=uniform)
#'   - real     = (Hmax=real, LAD=real)
#' Les coalitions mixtes (Hmax real + LAD uniform, ou Hmax mean + LAD real),
#' qui correspondent aux cas structurellement incoherents, sont ECARTEES.
#'
#' @param df_scores           Tibble from score_scenarios_vs_reference.
#' @param factorial_scenarios Liste de scenarios factoriels (pour mapping dynamique).
#' @param metric              "rmse" (defaut), "mae" ou "r2".
#' @return Named numeric vector de 8 valeurs (loss-form), avec attributs.
extract_shapley_coalitions_3p <- function(df_scores, factorial_scenarios = NULL,
                                           metric = c("rmse", "mae", "r2")) {
  metric <- match.arg(metric)

  # Mapping 3p coalition (LAI, fCover, Structure3D) -> cle 4p (LAI, Hmax, fCover, LAD)
  coal3p_to_4p <- c(
    "000" = "0000",  # all baseline
    "100" = "1000",  # +LAI
    "010" = "0010",  # +fCover
    "001" = "0101",  # +Structure3D = +Hmax & +LAD ensemble
    "110" = "1010",  # +LAI +fCover
    "101" = "1101",  # +LAI +Structure3D
    "011" = "0111",  # +fCover +Structure3D
    "111" = "1111"   # full real
  )

  if (!is.null(factorial_scenarios)) {
    coal4p_to_scenario <- setNames(
      vapply(factorial_scenarios, function(sc) sc$name, character(1)),
      vapply(factorial_scenarios, function(sc) {
        g <- sc$grid_row
        paste0(as.integer(g$LAI    == "real"),
               as.integer(g$Hmax   == "real"),
               as.integer(g$fCover == "real"),
               as.integer(g$LAD    == "real"))
      }, character(1))
    )
    coal3p_to_scenario <- coal4p_to_scenario[unname(coal3p_to_4p)]
    names(coal3p_to_scenario) <- names(coal3p_to_4p)
  } else {
    fallback <- c("0000" = "H1F_m_m_a_u", "1000" = "H1F_r_m_a_u",
                   "0010" = "H1F_m_m_r_u", "0101" = "H1F_m_r_a_r",
                   "1010" = "H1F_r_m_r_u", "1101" = "H1F_r_r_a_r",
                   "0111" = "H1F_m_r_r_r", "1111" = "H1F_r_r_r_r")
    coal3p_to_scenario <- fallback[unname(coal3p_to_4p)]
    names(coal3p_to_scenario) <- names(coal3p_to_4p)
  }

  if (!metric %in% names(df_scores))
    stop(sprintf("Colonne '%s' absente de df_scores", metric))

  raw_values <- sapply(coal3p_to_scenario, function(sc) {
    val <- df_scores[[metric]][df_scores$scenario == sc]
    if (length(val) == 0) {
      warning(sprintf("3p Shapley: scenario manquant : %s", sc))
      return(NA_real_)
    }
    val[1]
  })
  names(raw_values) <- names(coal3p_to_scenario)
  if (any(is.na(raw_values)))
    stop("3p Shapley: certaines coalitions manquent.")

  loss_values <- if (metric == "r2") 1 - raw_values else raw_values
  names(loss_values) <- names(raw_values)
  attr(loss_values, "metric")     <- metric
  attr(loss_values, "raw_values") <- raw_values
  attr(loss_values, "variables")  <- c("LAI", "fCover", "Structure3D")
  loss_values
}

#' Compute exact Shapley values for the 3-player coalition {LAI, fCover, Structure3D}.
#'
#' @param loss_values Output of extract_shapley_coalitions_3p().
#' @return List : $shapley, $loss_values, $v_empty, $total_effect, $total_sum,
#'   $check_diff, $metric, $units, $variables.
compute_shapley_exact_3p <- function(loss_values) {
  variables <- attr(loss_values, "variables")
  if (is.null(variables)) variables <- c("LAI", "fCover", "Structure3D")
  n <- length(variables)
  metric <- attr(loss_values, "metric"); if (is.null(metric)) metric <- "rmse"
  units  <- switch(metric, "rmse" = "°C", "mae" = "°C", "r2" = "", "")

  v_empty <- as.numeric(loss_values["000"])
  v <- function(b) v_empty - as.numeric(loss_values[b])
  make_bin <- function(bits_on, n = 3) {
    b <- rep("0", n); b[bits_on] <- "1"; paste(b, collapse = "")
  }

  shapley_vals <- sapply(seq_along(variables), function(var_bit) {
    others <- setdiff(seq_len(n), var_bit)
    total  <- 0
    for (s in 0:(n - 1)) {
      subsets <- if (s == 0) list(integer(0)) else combn(others, s, simplify = FALSE)
      w <- factorial(s) * factorial(n - s - 1) / factorial(n)
      for (S in subsets)
        total <- total + w * (v(make_bin(c(S, var_bit), n)) - v(make_bin(S, n)))
    }
    total
  })
  names(shapley_vals) <- variables

  total_shapley <- sum(shapley_vals)
  total_effect  <- v(paste(rep("1", n), collapse = ""))
  unit_suffix   <- if (nzchar(units)) sprintf(" %s", units) else ""

  cat(sprintf("── Valeurs Shapley 3-joueurs exactes (metric=%s) ──\n", metric))
  for (i in seq_along(variables))
    cat(sprintf("  phi(%-12s) = %+.4f%s  (%+.1f%% de v(N))\n",
                variables[i], shapley_vals[i], unit_suffix,
                100 * shapley_vals[i] / total_effect))
  cat(sprintf("\n  Somme phi = %+.4f  |  v(N) = %+.4f  |  ecart = %.2e\n\n",
              total_shapley, total_effect, abs(total_shapley - total_effect)))

  list(shapley = shapley_vals, loss_values = loss_values, rmse_values = loss_values,
       v_empty = v_empty, total_effect = total_effect,
       total_sum = total_shapley,
       check_diff = abs(total_shapley - total_effect),
       metric = metric, units = units, variables = variables)
}

#' Bar plot of 3-player Shapley attribution.
#'
#' @param shap_res Result of compute_shapley_exact_3p().
#' @param out_path PNG output path.
#' @return Invisible ggplot.
plot_shapley_3p <- function(shap_res,
                             out_path = "outputs/figures/02_shapley_3players.png") {
  variables  <- shap_res$variables
  metric     <- shap_res$metric
  units      <- shap_res$units
  metric_lbl <- toupper(metric)
  unit_paren <- if (nzchar(units)) sprintf(" (%s)", units) else ""
  unit_phi   <- if (nzchar(units)) sprintf(" %s", units) else ""

  df <- data.frame(
    variable = factor(variables, levels = variables),
    phi      = as.numeric(shap_res$shapley),
    label    = sprintf("%+.3f%s\n(%+.0f%%)",
                       as.numeric(shap_res$shapley), unit_phi,
                       100 * as.numeric(shap_res$shapley) / shap_res$total_effect)
  )
  pal <- c("LAI" = "#440154", "fCover" = "#31688e", "Structure3D" = "#fde725")

  p <- ggplot(df, aes(x = variable, y = phi, fill = variable)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    geom_text(aes(label = label, vjust = ifelse(phi >= 0, -0.3, 1.3)),
              size = 4.5, fontface = "bold") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.20))) +
    labs(
      title    = sprintf("Shapley attribution 3-joueurs — metric: %s", metric_lbl),
      subtitle = sprintf(
        "Structure3D = (Hmax, LAD) traites en bloc inseparable | v(N) = %.3f%s | Sigma phi = %.4f",
        shap_res$total_effect, unit_paren, shap_res$total_sum),
      caption = paste0(
        "8 coalitions (2^3) au lieu de 16 (2^4) : respecte le couplage Hmax-LAD impose par .lad_rescale().\n",
        "Variante methodologique du Shapley 4-joueurs principal pour controler l'artefact d'orthogonalite forcee."
      ),
      x = NULL, y = sprintf("phi%s", unit_paren)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold"),
          plot.caption    = element_text(size = 8.5, colour = "grey30", hjust = 0,
                                          lineheight = 1.3),
          axis.text.x     = element_text(face = "bold", size = 12))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 8, height = 6)
  cat(sprintf("  [Shapley 3p] PNG : %s\n", out_path))
  invisible(p)
}

# ---- Shapley bootstrap (stability over plot sub-samples) ---------------------

#' Bootstrap Shapley values over stratified sub-samples of plots.
#'
#' Reuses existing extracted scenario results — no MuSICA re-runs.
#' Each iteration draws prop=0.5 of plots *within each Cluster* (stratified,
#' without replacement), preserving the cLHS balance of the original design.
#'
#' @param df_all_scenarios  Stacked scenario results (must include all 16 H1F_*).
#' @param df_sample         cLHS sample dataframe with a Cluster column.
#' @param ref_scenario_name Reference scenario name.
#' @param n_boot            Number of bootstrap iterations (default 50).
#' @param factorial_scenarios Optional — passed to extract_shapley_coalitions for
#'   dynamic coalition mapping (see extract_shapley_coalitions()).
#' @return Long tibble with columns: boot_id, variable, phi, phi_pct.
bootstrap_shapley <- function(df_all_scenarios, df_sample,
                               ref_scenario_name,
                               n_boot              = 50,
                               factorial_scenarios = NULL,
                               metric              = c("rmse", "mae", "r2")) {
  metric <- match.arg(metric)
  set.seed(42)
  if (!"row_id" %in% names(df_sample))
    df_sample$row_id <- seq_len(nrow(df_sample))

  variables <- c("LAI", "Hmax", "fCover", "LAD")
  boot_list <- vector("list", n_boot)

  for (b in seq_len(n_boot)) {
    # Stratified: draw half of plots per cluster (preserves cLHS structure)
    idx <- df_sample %>%
      group_by(Cluster) %>%
      slice_sample(prop = 0.5) %>%
      pull(row_id)

    sub_xy <- df_sample[idx, c("x", "y")]

    df_sub <- df_all_scenarios %>%
      inner_join(sub_xy, by = c("x", "y"))

    shap_b <- tryCatch({
      df_scores_b  <- score_scenarios_vs_reference(df_sub, ref_scenario_name)
      coalitions_b <- extract_shapley_coalitions(df_scores_b, factorial_scenarios,
                                                   metric = metric)
      capture.output(res_b <- compute_shapley_exact(coalitions_b))
      tibble(
        boot_id  = b,
        variable = variables,
        phi      = as.numeric(res_b$shapley),
        phi_pct  = 100 * res_b$shapley / res_b$total_effect,
        metric   = metric
      )
    }, error = function(e) {
      warning(sprintf("[Shapley bootstrap] boot %d failed: %s", b, e$message))
      NULL
    })

    boot_list[[b]] <- shap_b
    if (b %% 10 == 0 || b == n_boot)
      cat(sprintf("  [Shapley bootstrap %s] boot %d/%d done\n", metric, b, n_boot))
  }
  bind_rows(boot_list)
}

#' Bar + jitter plot of bootstrap Shapley values with mean ± sd.
#'
#' @param df_boot    Result of bootstrap_shapley.
#' @param shap_exact Exact Shapley result (compute_shapley_exact) — overlaid as diamonds.
#' @param out_path   PNG output path.
#' @return Invisible ggplot.
plot_shapley_bootstrap <- function(df_boot, shap_exact = NULL,
                                    out_path = "outputs/figures/annex/A18_shapley_bootstrap.png") {
  variables <- c("LAI", "Hmax", "fCover", "LAD")
  pal_vars  <- c(LAI = "#440154", fCover = "#31688e", Hmax = "#35b779", LAD = "#fde725")

  # Detection metric : depuis df_boot$metric ou depuis shap_exact, defaut rmse
  metric <- if ("metric" %in% names(df_boot) && length(unique(df_boot$metric)) == 1)
              unique(df_boot$metric)
            else if (!is.null(shap_exact) && !is.null(shap_exact$metric))
              shap_exact$metric
            else "rmse"
  units      <- switch(metric, "rmse" = "°C", "mae" = "°C", "r2" = "", "")
  unit_paren <- if (nzchar(units)) sprintf(" (%s)", units) else ""
  metric_lbl <- toupper(metric)

  df_sum <- df_boot %>%
    group_by(variable) %>%
    summarise(
      mean_phi = mean(phi, na.rm = TRUE),
      sd_phi   = sd(phi,  na.rm = TRUE),
      ci_lo    = quantile(phi, 0.025, na.rm = TRUE),
      ci_hi    = quantile(phi, 0.975, na.rm = TRUE),
      n_boot   = dplyr::n(),
      .groups  = "drop"
    ) %>%
    mutate(variable = factor(variable, levels = variables))

  df_pts <- df_boot %>% mutate(variable = factor(variable, levels = variables))

  p <- ggplot(df_sum, aes(x = variable, y = mean_phi, fill = variable)) +
    geom_col(width = 0.6, colour = "grey20", linewidth = 0.3, alpha = 0.85) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.8) +
    geom_jitter(data = df_pts, aes(x = variable, y = phi),
                width = 0.12, size = 2.8, alpha = 0.75, colour = "grey25",
                inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = pal_vars) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.20)))

  if (!is.null(shap_exact)) {
    df_exact <- tibble(
      variable  = factor(names(shap_exact$shapley), levels = variables),
      phi_exact = as.numeric(shap_exact$shapley)
    )
    p <- p + geom_point(data = df_exact,
                        aes(x = variable, y = phi_exact),
                        shape = 18, size = 6, colour = "#d73027",
                        inherit.aes = FALSE)
  }

  lad_ci <- df_sum %>% filter(variable == "LAD")
  lad_ci_str <- if (nrow(lad_ci) == 1)
    sprintf("95%% CI LAD = [%+.3f ; %+.3f]%s",
            lad_ci$ci_lo, lad_ci$ci_hi, unit_paren)
  else ""

  n_b <- dplyr::n_distinct(df_boot$boot_id)
  p <- p +
    labs(
      title    = sprintf("Shapley bootstrap (metric: %s) — stability over %d sub-samples",
                          metric_lbl, n_b),
      subtitle = paste0("Bars = mean | Line = 95% empirical CI | points = draws | ◆ = exact Shapley\n",
                        lad_ci_str),
      x        = NULL,
      y        = sprintf("phi%s", unit_paren)
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title      = element_text(face = "bold"),
          axis.text.x     = element_text(face = "bold", size = 12))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 8, height = 6)
  cat(sprintf("  [Shapley bootstrap] PNG : %s\n", out_path))
  invisible(p)
}

#' Wilcoxon signed-rank test on the marginal gain of the last Forward step.
#'
#' Tests whether adding the full LAD profile (H1f_4) gives significantly smaller
#' absolute errors than the fCover-only scenario (H1f_3), and reports Cohen's d.
#'
#' @param df_all          Stacked scenario results (x, y, date, Delta_Tmax, scenario).
#' @param ref_name        Reference scenario name.
#' @param sc_reduced      Name of the reduced scenario (default H1f_3_LAI_Hmax_fCover).
#' @param sc_full         Name of the full scenario (default H1f_4_Full_real).
#' @return Named list: p_value, median_gain_degC, cohen_d, n_pairs.
test_forward_last_step <- function(df_all, ref_name,
                                    sc_reduced = "H1f_3_LAI_Hmax_fCover",
                                    sc_full    = "H1f_4_Full_real") {
  df_ref  <- df_all %>% filter(scenario == ref_name) %>%
    dplyr::select(x, y, date, Delta_ref = Delta_Tmax)

  get_err <- function(sc) {
    df_all %>% filter(scenario == sc) %>%
      inner_join(df_ref, by = c("x", "y", "date")) %>%
      mutate(abs_err = abs(Delta_Tmax - Delta_ref)) %>%
      pull(abs_err)
  }
  e3 <- get_err(sc_reduced)
  e4 <- get_err(sc_full)

  n_pairs <- min(length(e3), length(e4))
  e3 <- e3[seq_len(n_pairs)]; e4 <- e4[seq_len(n_pairs)]

  wt       <- wilcox.test(e3, e4, paired = TRUE, alternative = "greater")
  med_gain <- median(e3 - e4, na.rm = TRUE)
  pooled_sd <- sd(c(e3, e4), na.rm = TRUE)
  cohen_d  <- (mean(e3, na.rm = TRUE) - mean(e4, na.rm = TRUE)) / pooled_sd

  cat(sprintf(
    "\n[Wilcoxon H1f_3 vs H1f_4]\n  n=%d paires | W=%g | p=%.3e\n  Gain médian = %.4f°C | Cohen d = %.4f\n",
    n_pairs, wt$statistic, wt$p.value, med_gain, cohen_d
  ))
  invisible(list(p_value = wt$p.value, median_gain_degC = med_gain,
                 cohen_d = cohen_d, n_pairs = n_pairs))
}

#' Horizontal RMSE bar chart for the parsimony ladder (figure A24).
#'
#' Displays 5 key coalitions from Null baseline to Full LiDAR as a ladder,
#' labelled with RMSE and R². Requires score_scenarios_vs_reference() output
#' with an r2 column.
#'
#' @param df_scores  Tibble from score_scenarios_vs_reference (must contain r2).
#' @return ggplot object, or NULL if any required scenario is absent.
plot_parsimony_ladder <- function(df_scores) {
  LADDER <- data.frame(
    scenario = c("H1f_0_Null_baseline", "H1f_1_LAI_only",
                 "H1F_r_m_r_u",         "H1f_3_LAI_Hmax_fCover",
                 "H1f_4_Full_real"),
    label    = c("Null baseline",        "LAI only (S2+GEDI)",
                 "LAI + fCover (S2+GEDI)", "LAI + Hmax + fCover",
                 "Full LiDAR"),
    stringsAsFactors = FALSE
  )

  present <- LADDER$scenario %in% df_scores$scenario
  if (!all(present)) {
    missing_sc <- LADDER$scenario[!present]
    message(sprintf("[parsimony_ladder] Missing: %s — partial ladder only",
                    paste(missing_sc, collapse = ", ")))
    LADDER <- LADDER[present, ]
  }
  if (nrow(LADDER) < 2) return(NULL)

  df_plot <- df_scores %>%
    dplyr::filter(scenario %in% LADDER$scenario) %>%
    dplyr::left_join(LADDER, by = "scenario") %>%
    dplyr::mutate(label = factor(label, levels = LADDER$label))

  pal <- c(
    "Null baseline"           = "#aaaaaa",
    "LAI only (S2+GEDI)"      = "#fdae61",
    "LAI + fCover (S2+GEDI)"  = "#f46d43",
    "LAI + Hmax + fCover"     = "#d73027",
    "Full LiDAR"              = "#4575b4"
  )

  x_max <- max(df_plot$rmse, na.rm = TRUE) * 1.45
  has_r2 <- "r2" %in% names(df_plot) && !all(is.na(df_plot$r2))
  df_plot$bar_lbl <- if (has_r2)
    sprintf("RMSE = %.3f°C | R² = %.2f", df_plot$rmse, df_plot$r2)
  else
    sprintf("RMSE = %.3f°C", df_plot$rmse)

  ggplot2::ggplot(df_plot, ggplot2::aes(x = rmse, y = label, fill = label)) +
    ggplot2::geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = bar_lbl),
                       hjust = -0.06, size = 3.6) +
    ggplot2::scale_fill_manual(values = pal, na.value = "grey70") +
    ggplot2::scale_x_continuous(
      limits = c(0, x_max),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title    = "Parsimony ladder — RMSE vs full-LiDAR reference",
      subtitle = "LAI + fCover (Sentinel-2 + GEDI) captures >86% of LiDAR signal",
      x        = "RMSE vs REF_all_real (°C)",
      y        = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.major.y = ggplot2::element_blank())
}
