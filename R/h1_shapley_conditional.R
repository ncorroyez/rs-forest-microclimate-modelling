# ==============================================================================
# Conditional Shapley analysis for LAD stratified by Hmax level.
#
# Justification : dans le 2^4 factoriel, LAD et Hmax sont co-definis
# (LAD = fonction sur [0, Hmax]). Les coalitions "LAD=real, Hmax=mu" sont
# structurellement incoherentes — .lad_rescale() comprime le profil reel sur
# un domaine artificiel. Cette analyse decompose phi_LAD selon que Hmax est
# coherent (real) ou incoherent (mean) dans la coalition.
#
# phi_LAD_global = mean(phi_LAD_given_Hmax_real, phi_LAD_given_Hmax_baseline)
# ==============================================================================

# Mapping table des 8 paires (S sans LAD, S avec LAD)
# Code bit order : LAI | Hmax | fCover | LAD
# Group A = Hmax=real in S (coherent domain for LAD)
# Group B = Hmax=mean in S (incoherent domain — .lad_rescale compresses)
.LAD_PAIRS <- list(
  list(S = "0000", S_LAD = "0001", w = 1/4,  group = "B"),
  list(S = "1000", S_LAD = "1001", w = 1/12, group = "B"),
  list(S = "0100", S_LAD = "0101", w = 1/12, group = "A"),
  list(S = "0010", S_LAD = "0011", w = 1/12, group = "B"),
  list(S = "1100", S_LAD = "1101", w = 1/12, group = "A"),
  list(S = "1010", S_LAD = "1011", w = 1/12, group = "B"),
  list(S = "0110", S_LAD = "0111", w = 1/12, group = "A"),
  list(S = "1110", S_LAD = "1111", w = 1/4,  group = "A")
)

#' Compute conditional Shapley decomposition of phi_LAD by Hmax level.
#'
#' @param loss_values Named numeric vector of 16 loss values from
#'   extract_shapley_coalitions() (names = "0000"..."1111"). For R² metric,
#'   pass (1 - r2) values as usual.
#' @return List : phi_global, phi_A, phi_B, phi_A_raw, phi_B_raw,
#'   deltas (data.frame), metric, check_diff.
compute_shapley_conditional_LAD <- function(loss_values) {
  metric <- attr(loss_values, "metric")
  if (is.null(metric)) metric <- "rmse"

  # delta_S = loss(S) - loss(S ∪ {LAD})  (positive = LAD reduces loss)
  deltas <- vapply(.LAD_PAIRS, function(pp) {
    as.numeric(loss_values[pp$S]) - as.numeric(loss_values[pp$S_LAD])
  }, numeric(1))

  groups  <- vapply(.LAD_PAIRS, function(pp) pp$group, character(1))
  weights <- vapply(.LAD_PAIRS, function(pp) pp$w, numeric(1))

  phi_A_raw  <- sum(weights[groups == "A"] * deltas[groups == "A"])
  phi_B_raw  <- sum(weights[groups == "B"] * deltas[groups == "B"])
  phi_global <- phi_A_raw + phi_B_raw

  # Conditional Shapleys : normalize by group weight (1/2 each)
  phi_A <- 2 * phi_A_raw
  phi_B <- 2 * phi_B_raw

  # Sanity check : (phi_A + phi_B) / 2 == phi_global
  check <- abs((phi_A + phi_B) / 2 - phi_global)
  if (check > 1e-9)
    stop(sprintf("[Conditional Shapley] sanity check failed: |delta| = %.2e", check))

  list(
    phi_global = phi_global,
    phi_A      = phi_A,
    phi_B      = phi_B,
    phi_A_raw  = phi_A_raw,
    phi_B_raw  = phi_B_raw,
    deltas     = data.frame(
      pair  = seq_along(.LAD_PAIRS),
      S     = vapply(.LAD_PAIRS, `[[`, "", "S"),
      S_LAD = vapply(.LAD_PAIRS, `[[`, "", "S_LAD"),
      group = groups,
      w     = weights,
      delta = deltas,
      stringsAsFactors = FALSE
    ),
    metric     = metric,
    check_diff = check
  )
}

#' Bootstrap conditional Shapley for LAD (stratified by cluster).
#'
#' Réplique la logique de bootstrap_shapley() en stockant phi_A, phi_B et
#' phi_global par itération.
bootstrap_shapley_conditional_LAD <- function(df_all_scenarios, df_sample,
                                               ref_scenario_name,
                                               n_boot              = 50,
                                               factorial_scenarios = NULL,
                                               metric              = "rmse") {
  set.seed(44)   # different from bootstrap_shapley (42) and ANOVA bootstrap (43)
  if (!"row_id" %in% names(df_sample))
    df_sample$row_id <- seq_len(nrow(df_sample))

  out <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    idx <- df_sample %>%
      group_by(Cluster) %>%
      slice_sample(prop = 0.5) %>%
      pull(row_id)
    sub_xy <- df_sample[idx, c("x", "y")]
    df_sub <- df_all_scenarios %>% inner_join(sub_xy, by = c("x", "y"))

    res <- tryCatch({
      df_scores_b  <- score_scenarios_vs_reference(df_sub, ref_scenario_name)
      coalitions_b <- extract_shapley_coalitions(df_scores_b, factorial_scenarios,
                                                   metric = metric)
      compute_shapley_conditional_LAD(coalitions_b)
    }, error = function(e) {
      warning(sprintf("[Conditional Shapley boot %d] %s", b, e$message))
      NULL
    })

    if (!is.null(res))
      out[[b]] <- tibble(
        boot_id    = b,
        phi_global = res$phi_global,
        phi_A      = res$phi_A,
        phi_B      = res$phi_B,
        metric     = metric
      )

    if (b %% 10 == 0 || b == n_boot)
      cat(sprintf("  [Conditional Shapley %s] boot %d/%d done\n", metric, b, n_boot))
  }
  dplyr::bind_rows(out)
}

#' Plot conditional Shapley LAD with bootstrap CI (3 horizontal bars).
plot_shapley_conditional_LAD <- function(shap_cond, df_boot,
                                          out_path,
                                          metric = "rmse") {
  units      <- switch(metric, "rmse" = "°C", "mae" = "°C", "r2" = "", "")
  unit_paren <- if (nzchar(units)) sprintf(" (%s)", units) else ""
  unit_phi   <- if (nzchar(units)) sprintf(" %s", units) else ""

  ci <- dplyr::summarise(df_boot,
    phi_global_lo = quantile(phi_global, 0.025, na.rm = TRUE),
    phi_global_hi = quantile(phi_global, 0.975, na.rm = TRUE),
    phi_A_lo      = quantile(phi_A,      0.025, na.rm = TRUE),
    phi_A_hi      = quantile(phi_A,      0.975, na.rm = TRUE),
    phi_B_lo      = quantile(phi_B,      0.025, na.rm = TRUE),
    phi_B_hi      = quantile(phi_B,      0.975, na.rm = TRUE)
  )

  df_plot <- data.frame(
    label = c("phi_LAD  (global)",
               "phi_LAD  |  Hmax = real (coherent)",
               "phi_LAD  |  Hmax = baseline (mu)"),
    phi   = c(shap_cond$phi_global, shap_cond$phi_A, shap_cond$phi_B),
    ci_lo = c(ci$phi_global_lo, ci$phi_A_lo, ci$phi_B_lo),
    ci_hi = c(ci$phi_global_hi, ci$phi_A_hi, ci$phi_B_hi)
  )
  df_plot$label <- factor(df_plot$label, levels = rev(df_plot$label))
  df_plot$sig   <- ifelse((df_plot$ci_lo > 0) | (df_plot$ci_hi < 0), "*", "")

  p <- ggplot(df_plot, ggplot2::aes(x = phi, y = label, fill = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                         colour = "grey50", linewidth = 0.6) +
    ggplot2::geom_col(width = 0.55, colour = "grey20", linewidth = 0.4) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                             height = 0.18, linewidth = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%+.3f%s%s\n[%+.3f, %+.3f]",
                                    phi, unit_phi, sig, ci_lo, ci_hi),
                    x = ifelse(phi >= 0, ci_hi, ci_lo)),
      hjust = ifelse(df_plot$phi >= 0, -0.1, 1.1),
      size = 3.7, fontface = "bold", lineheight = 1.0
    ) +
    ggplot2::scale_fill_manual(values = c(
      "phi_LAD  (global)"                    = "#888888",
      "phi_LAD  |  Hmax = real (coherent)"   = "#440154",
      "phi_LAD  |  Hmax = baseline (mu)"     = "#fde725"
    ), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.45, 0.45))) +
    ggplot2::labs(
      title    = "Conditional Shapley analysis for LAD stratified by Hmax",
      subtitle = paste0(
        "phi_LAD_global = average of two conditionals.  ",
        "Group A : LAD acts on its real Hmax domain.  ",
        "Group B : LAD compressed onto baseline Hmax_mean by .lad_rescale()."
      ),
      caption  = sprintf("Metric : %s | %d bootstrap iterations stratified by cluster | 95%% empirical CI | * = CI excludes 0",
                         toupper(metric), nrow(df_boot)),
      x = sprintf("phi_LAD%s", unit_paren),
      y = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey30",
                                             lineheight = 1.2),
      plot.caption  = ggplot2::element_text(size = 8, colour = "grey40",
                                             hjust = 0),
      axis.text.y   = ggplot2::element_text(size = 11, face = "bold")
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 12, height = 6)
  cat(sprintf("  [Shapley conditional LAD] PNG : %s\n", out_path))
  invisible(p)
}
