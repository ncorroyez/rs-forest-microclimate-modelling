# ==============================================================================
# Chapter 1 — GAMM emulator (structure + macroclimate)
# ==============================================================================

#' Prepare scaled inputs for the GAMM, integrating macroclimate and shape drivers.
#'
#' @param df_daily      Daily scenario results for the reference scenario.
#' @param df_sample     cLHS sample dataframe.
#' @param fpc_scores    Optional FPC score matrix (n_plots × n_fpc).
#' @param h_median_vec  Optional H_median vector (length = nrow(df_sample)).
#' @return Dataframe ready for fit_reference_gamm, with scaled predictors.
prepare_gamm_data <- function(df_daily, df_sample, fpc_scores = NULL, h_median_vec = NULL) {
  df <- df_daily %>%
    mutate(plot_id = as.factor(paste0("X", x, "_Y", y)),
           date_factor = as.factor(date))

  if (!is.null(h_median_vec)) {
    df_sample$H_median <- h_median_vec
    df <- df %>% left_join(
      df_sample %>% mutate(x = round(x), y = round(y)) %>% dplyr::select(x, y, H_median),
      by = c("x", "y"))
  }

  if ("Cluster" %in% names(df)) df$Cluster <- as.factor(df$Cluster)

  scale_cols <- c("LAI", "Hmax", "fCover", "H_median", "Tmax_macro",
                  "Wind_mean", "Rad_mean", "VPD_mean")
  cols_to_scale <- intersect(scale_cols, names(df))
  for (nm in cols_to_scale) df[[paste0(nm, "_sc")]] <- as.numeric(scale(df[[nm]]))

  if (!is.null(fpc_scores)) {
    fpc_df <- as.data.frame(fpc_scores)
    names(fpc_df) <- paste0("FPC", seq_len(ncol(fpc_df)), "_sc")
    fpc_df[] <- lapply(fpc_df, function(v) as.numeric(scale(v)))
    fpc_df$x <- round(df_sample$x); fpc_df$y <- round(df_sample$y)
    df <- df %>% left_join(fpc_df, by = c("x", "y"))
    n_na_fpc <- sum(is.na(df$FPC1_sc))
    if (n_na_fpc > 0)
      warning(sprintf(
        "[prepare_gamm_data] FPC join produced %d NA rows (%.0f%% of %d) — coordinate mismatch between df_daily and fpc_scores. bam() will silently drop these rows.",
        n_na_fpc, 100 * n_na_fpc / nrow(df), nrow(df)
      ))
  }
  cat(sprintf("  [GAMM] df_gamm: %d rows, %d unique plots, %d unique dates\n",
              nrow(df), n_distinct(df$plot_id), n_distinct(df$date_factor)))
  return(df)
}

#' Build a safe s() term: cap k at n_unique-1 to prevent mgcv construction errors.
#'
#' mgcv raises "fewer unique covariate combinations than specified maximum degrees
#' of freedom" when k > n_unique for any smooth term. This is common for
#' near-constant predictors in summer (e.g. precipitation).
#'
#' @param var  Column name in data.
#' @param k    Desired k (will be silently capped).
#' @param data Dataframe containing var.
#' @return Character string like "s(var, k=N)".
safe_smooth_term <- function(var, k, data) {
  n_u    <- length(unique(data[[var]][!is.na(data[[var]])]))
  k_safe <- max(3L, min(as.integer(k), n_u - 1L))
  sprintf("s(%s, k=%d)", var, k_safe)
}

#' Fit the reference GAMM with toggles for macroclimate and shape predictor.
#'
#' Uses bam() with scat() family and fREML for efficiency on large datasets.
#'
#' @param df_gamm    Dataframe from prepare_gamm_data.
#' @param shape_type Shape predictor: "CLUSTER", "FPC", "H_MEDIAN", or "NONE".
#' @param clim_vars  Character vector of ERA5 base names to include as climate
#'   smooths (e.g. c("Tmax_macro","Wind_mean","Rad_mean","VPD_mean")).
#'   The function appends "_sc" automatically. Pass character(0) to exclude all.
#' @return Fitted bam model object.
fit_reference_gamm <- function(df_gamm, shape_type = "NONE",
                                clim_vars = c("Tmax_macro", "Wind_mean",
                                              "Rad_mean", "VPD_mean")) {
  base_terms <- c(
    safe_smooth_term("LAI_sc",    7L, df_gamm),
    safe_smooth_term("Hmax_sc",   7L, df_gamm),
    safe_smooth_term("fCover_sc", 7L, df_gamm)
  )

  clim_terms <- if (length(clim_vars) > 0) {
    cat(sprintf("  [GAMM] Climate drivers: %s\n", paste(clim_vars, collapse = ", ")))
    vapply(paste0(clim_vars, "_sc"), function(v) safe_smooth_term(v, 7L, df_gamm),
           character(1))
  } else {
    cat("  [GAMM] Structural drivers only (no climate covariates)...\n")
    character(0)
  }

  shape_terms <- switch(shape_type,
                        "CLUSTER" = {
                          cat("  [GAMM] Including profile Cluster type as categorical predictor...\n")
                          c("Cluster")
                        },
                        "FPC" = {
                          cat("  [GAMM] Including FPC scores as shape predictors...\n")
                          c(safe_smooth_term("FPC1_sc", 7L, df_gamm),
                            safe_smooth_term("FPC2_sc", 7L, df_gamm),
                            safe_smooth_term("FPC3_sc", 7L, df_gamm))
                        },
                        "H_MEDIAN" = {
                          cat("  [GAMM] Including structural shape (H_median) as predictor...\n")
                          c(safe_smooth_term("H_median_sc", 7L, df_gamm))
                        },
                        {
                          cat("  [GAMM] Structural macro-metrics only (no specific shape variable)...\n")
                          character(0)
                        }
  )

  # date_factor as RE absorbs unmodeled daily climate variance not fully captured
  # by ERA5 smooths (local cloud cover, phenological state, measurement noise).
  # RE keeps df minimal and prevents per-day residuals from biasing fixed structural
  # effects. plot_id as RE controls spatial pseudo-replication (~120 obs per plot).
  re_terms <- c("s(date_factor, bs='re')", "s(plot_id, bs='re')")
  rhs  <- paste(c(base_terms, clim_terms, shape_terms, re_terms), collapse = " + ")
  form <- as.formula(paste("Delta_Tmax ~", rhs))

  # discrete=FALSE: avoids mgcv pre-binning that fails when df_gamm contains
  # many extra numeric columns (FPC*_sc, H_median_sc) not in the formula.
  # ~36k rows — fREML without discrete is fast enough.
  #
  # scat() (scaled-t) is used over Gaussian because DeltaTmax >= 0 has a heavy
  # right tail (large cooling events) and variance that increases with the fitted
  # value (heteroscedasticity visible in scale-location plot A11). scat() handles
  # both without a log transform that would change scale interpretation.
  # AIC reduction vs Gaussian family is ~50% (see A11 diagnostics).
  bam(form, data = df_gamm, family = scat(), method = "fREML", discrete = FALSE)
}

#' Plot marginal effects of the GAMM predictors.
#'
#' @param gam_model  Fitted bam/gam model.
#' @param shape_type Shape predictor type used (see fit_reference_gamm).
#' @param clim_vars  ERA5 base names included in the model (appends "_sc").
#' @return patchwork ggplot.
plot_gamm_marginal_effects <- function(gam_model, shape_type = "NONE",
                                        clim_vars = c("Tmax_macro", "Wind_mean",
                                                      "Rad_mean", "VPD_mean")) {
  terms_to_plot <- c("LAI_sc", "Hmax_sc")

  if (shape_type == "CLUSTER")  terms_to_plot <- c(terms_to_plot, "Cluster")
  if (shape_type == "FPC")      terms_to_plot <- c(terms_to_plot, "FPC1_sc")
  if (shape_type == "H_MEDIAN") terms_to_plot <- c(terms_to_plot, "H_median_sc")
  if (length(clim_vars) > 0)
    terms_to_plot <- c(terms_to_plot, intersect(paste0(clim_vars, "_sc"),
                                                names(model.frame(gam_model))))

  plots <- lapply(terms_to_plot, function(term) {
    pred <- ggpredict(gam_model, terms = term)
    plot(pred) +
      coord_cartesian(ylim = c(-3, 3)) +
      labs(title = paste("Effect of", gsub("_sc", "", term)),
           x = if (term == "Cluster") "Profile cluster type" else paste(gsub("_sc", "", term), "(Scaled)"),
           y = "Predicted ΔTmax (°C)") +
      theme_bw() + theme(plot.title = element_text(size = 10))
  })

  wrap_plots(plots, ncol = 2) + plot_annotation(title = "GAMM Marginal Effects", subtitle = "Holding all other variables at their mean")
}

#' Residual diagnostics for the reference GAMM (QQ, histogram, fitted vs residuals, AIC, ACF).
#'
#' @param gam_scat Fitted scat() bam model.
#' @param df_gamm  Data used to fit the model.
#' @param plot     Whether to produce base-graphics plots (default TRUE).
#' @return Dataframe with diagnostic statistics.
diagnose_residuals <- function(gam_scat, df_gamm, plot = TRUE) {
  gam_gauss <- bam(formula(gam_scat), data = df_gamm, family = gaussian(), method = "fREML", discrete = TRUE)
  resid_g <- residuals(gam_gauss, type = "deviance")
  resid_r <- residuals(gam_scat,  type = "response")
  idx_sw  <- sample(length(resid_g), min(5000, length(resid_g))); sw_test <- shapiro.test(resid_g[idx_sw])
  excess_kurt <- mean((resid_g - mean(resid_g))^4) / sd(resid_g)^4 - 3
  aic_g <- AIC(gam_gauss); aic_s <- AIC(gam_scat)

  if (plot) {
    op <- par(mfrow = c(2, 3)); on.exit(par(op))
    qqnorm(resid_g, main = "QQ-plot (Gaussian residuals)", pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3)); qqline(resid_g, col = "red", lwd = 2)
    hist(resid_g, breaks = 80, freq = FALSE, main = "Residual distribution", xlab = "Residual", col = "lightblue", border = "white"); curve(dnorm(x, mean = mean(resid_g), sd = sd(resid_g)), add = TRUE, col = "red", lwd = 2)
    plot(fitted(gam_gauss), resid_g, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3), xlab = "Fitted", ylab = "Residual", main = "Residuals vs fitted"); abline(h = 0, col = "red", lwd = 2); lines(lowess(fitted(gam_gauss), resid_g), col = "blue", lwd = 2)
    barplot(c(Gaussian = aic_g, `scat()` = aic_s), col = c("lightblue", "salmon"), main = "AIC comparison", ylab = "AIC")
    # ACF on scat() response residuals sorted by (plot_id, date) within each plot.
    # Sorting within-plot ensures lags are temporal (not cross-plot mixing).
    # Residual order follows df_gamm row order; sort by plot×date to make lags
    # meaningful within-plot, then compute mean per-lag across plots.
    if ("plot_id" %in% names(df_gamm) && "date_factor" %in% names(df_gamm)) {
      df_acf <- data.frame(plot_id = df_gamm$plot_id,
                           date    = df_gamm$date_factor,
                           r       = resid_r) %>%
        dplyr::arrange(plot_id, date)
      acf(df_acf$r, lag.max = 14,
          main = "ACF scat() residuals (sorted by plot × date)",
          col  = "#31688e", lwd = 2)
    } else {
      acf(resid_r, lag.max = 14, main = "ACF scat() response residuals",
          col = "#31688e", lwd = 2)
    }
    # Variance vs fitted (heteroscedasticity)
    plot(fitted(gam_scat), abs(resid_r), pch = 16, cex = 0.3,
         col = rgb(0, 0, 0.5, 0.3), xlab = "Fitted (scat)", ylab = "|Residual|",
         main = "Scale-Location (scat)")
    lines(lowess(fitted(gam_scat), abs(resid_r)), col = "red", lwd = 2)
  }

  data.frame(n_obs = length(resid_g), shapiro_W = unname(sw_test$statistic), shapiro_p = sw_test$p.value, excess_kurt = excess_kurt, aic_gauss = aic_g, aic_scat = aic_s, delta_aic = aic_g - aic_s, scat_preferred = (aic_g - aic_s) > 10)
}

#' Fit and compare 3 GAMM variants for the shape-predictor choice.
#'
#' Variants compared:
#'   1. NONE       — structural macro-metrics only (LAI, Hmax, fCover + clim + RE)
#'   2. H_MEDIAN   — adds s(H_median_sc)
#'   3. FPC        — adds s(FPC1_sc) + s(FPC2_sc) + s(FPC3_sc)
#'
#' The df_gamm dataframe must already contain H_median_sc and FPC*_sc columns
#' (prepare_gamm_data must be called with both fpc_scores and h_median_vec).
#'
#' @param df_gamm   Dataframe from prepare_gamm_data.
#' @param clim_vars ERA5 base names to include (see fit_reference_gamm). Default: full set.
#' @param out_path  Optional PNG path for the comparison bar chart.
#' @return Invisible list: $models (named list of 3 bam objects), $comparison (tibble).
compare_gamm_variants <- function(df_gamm,
                                   clim_vars = c("Tmax_macro", "Wind_mean",
                                                 "Rad_mean", "VPD_mean"),
                                   out_path = NULL) {
  cat("  [GAMM compare] Fitting 3 shape-predictor variants...\n")

  variants <- c("NONE", "H_MEDIAN", "FPC")
  models <- setNames(
    lapply(variants, function(vt) {
      cat(sprintf("    → shape_type = %s\n", vt))
      tryCatch(
        fit_reference_gamm(df_gamm, shape_type = vt, clim_vars = clim_vars),
        error = function(e) {
          warning(sprintf("[compare_gamm_variants] %s failed: %s", vt, e$message))
          NULL
        }
      )
    }),
    variants
  )

  # Keep only variants that converged
  models <- Filter(Negate(is.null), models)
  if (length(models) == 0) stop("[compare_gamm_variants] all 3 variants failed")

  df_cmp <- purrr::imap_dfr(models, function(m, nm) {
    s <- summary(m)
    tibble(
      variant  = nm,
      dev_expl = round(s$dev.expl * 100, 2),
      aic      = round(AIC(m), 1),
      n_obs    = nrow(model.frame(m))
    )
  }) %>%
    mutate(delta_aic = round(aic - min(aic), 1))

  cat("\n── GAMM variant comparison ──\n")
  print(as.data.frame(df_cmp))
  cat("\n")

  if (!is.null(out_path)) {
    df_long <- df_cmp %>%
      pivot_longer(c(dev_expl, delta_aic), names_to = "metric", values_to = "value") %>%
      mutate(
        metric  = dplyr::recode(metric,
                                dev_expl  = "Deviance explained (%)",
                                delta_aic = "ΔAIC (lower = better)"),
        variant = factor(variant, levels = c("NONE", "H_MEDIAN", "FPC"))
      )

    p <- ggplot(df_long, aes(x = variant, y = value, fill = variant)) +
      geom_col(width = 0.6, colour = "grey30", linewidth = 0.3) +
      geom_text(aes(label = round(value, 1)),
                vjust = -0.4, size = 4, fontface = "bold") +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
      facet_wrap(~ metric, scales = "free_y") +
      scale_fill_manual(values = c(NONE = "#aaaaaa", H_MEDIAN = "#2a9d8f", FPC = "#e76f51")) +
      labs(
        title    = "GAMM — comparison of 3 shape-descriptor variants",
        subtitle = "NONE: integrated metrics | H_MEDIAN: LAD centre of gravity | FPC: 3 first harmonics",
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", strip.text = element_text(face = "bold"))

    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    save_plot(p, out_path, width = 10, height = 5)
    cat(sprintf("  [GAMM compare] PNG: %s\n", out_path))
  }

  invisible(list(models = models, comparison = df_cmp))
}

#' Fit the FPC GAMM without s(Hmax_sc) — robustness test for FPC1/Hmax collinearity.
#'
#' Identical to fit_reference_gamm(shape_type="FPC") except s(Hmax_sc) is dropped
#' from the base structural terms. All other smooths are preserved:
#' s(LAI_sc) + s(fCover_sc) + s(FPC1_sc) + s(FPC2_sc) + s(FPC3_sc)
#' + climate smooths + s(date_factor, bs="re") + s(plot_id, bs="re").
#' Same scat() family and fREML method.
#'
#' @param df_gamm   Dataframe from prepare_gamm_data (must contain FPC*_sc columns).
#' @param clim_vars ERA5 base names to include as climate smooths.
#' @return Fitted bam model object.
fit_gamm_without_hmax <- function(df_gamm,
                                   clim_vars = c("Tmax_macro", "Wind_mean", "Rad_mean")) {
  base_terms <- c(
    safe_smooth_term("LAI_sc",    7L, df_gamm),
    safe_smooth_term("fCover_sc", 7L, df_gamm)
    # s(Hmax_sc) deliberately omitted — sensitivity test
  )

  clim_terms <- if (length(clim_vars) > 0) {
    cat(sprintf("  [GAMM no-Hmax] Climate drivers: %s\n", paste(clim_vars, collapse = ", ")))
    vapply(paste0(clim_vars, "_sc"), function(v) safe_smooth_term(v, 7L, df_gamm),
           character(1))
  } else {
    character(0)
  }

  shape_terms <- c(
    safe_smooth_term("FPC1_sc", 7L, df_gamm),
    safe_smooth_term("FPC2_sc", 7L, df_gamm),
    safe_smooth_term("FPC3_sc", 7L, df_gamm)
  )

  re_terms <- c("s(date_factor, bs='re')", "s(plot_id, bs='re')")
  rhs  <- paste(c(base_terms, clim_terms, shape_terms, re_terms), collapse = " + ")
  form <- as.formula(paste("Delta_Tmax ~", rhs))
  cat("  [GAMM no-Hmax] Fitting:", deparse(form), "\n")
  bam(form, data = df_gamm, family = scat(), method = "fREML", discrete = FALSE)
}

#' Two-panel comparison of FPC1 marginal effect with vs without Hmax in the model.
#'
#' Extracts the s(FPC1_sc) smooth from both models, plots them side by side on a
#' shared y-axis, and returns the smooth data and amplitudes for export.
#'
#' @param gam_with_hmax    Fitted bam from fit_reference_gamm(shape_type="FPC").
#' @param gam_without_hmax Fitted bam from fit_gamm_without_hmax().
#' @return Invisible list: $plot (patchwork ggplot), $df_smooth (data.frame),
#'   $amp_with (°C), $amp_without (°C).
plot_fpc1_robustness_to_hmax <- function(gam_with_hmax, gam_without_hmax) {
  get_fpc1_stats <- function(m) {
    st  <- summary(m)$s.table
    idx <- grep("FPC1_sc", rownames(st))
    list(edf  = round(st[idx, "edf"], 2),
         pval = st[idx, "p-value"])
  }

  stats_w  <- get_fpc1_stats(gam_with_hmax)
  stats_wo <- get_fpc1_stats(gam_without_hmax)

  pred_w  <- as.data.frame(ggpredict(gam_with_hmax,    terms = "FPC1_sc"))
  pred_wo <- as.data.frame(ggpredict(gam_without_hmax, terms = "FPC1_sc"))

  amp_w  <- max(pred_w$predicted)  - min(pred_w$predicted)
  amp_wo <- max(pred_wo$predicted) - min(pred_wo$predicted)

  y_lim <- range(c(pred_w$conf.low,  pred_w$conf.high,
                   pred_wo$conf.low, pred_wo$conf.high))

  make_panel <- function(pred, stats, amp, panel_title) {
    pval_str <- if (stats$pval < 0.001) "p<0.001" else sprintf("p=%.3f", stats$pval)
    ggplot2::ggplot(pred, ggplot2::aes(x = x, y = predicted)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                           alpha = 0.2, fill = "#e76f51") +
      ggplot2::geom_line(colour = "#e76f51", linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::coord_cartesian(ylim = y_lim) +
      ggplot2::labs(
        title    = panel_title,
        subtitle = sprintf("edf = %.2f | %s | amplitude = %.3f°C",
                           stats$edf, pval_str, amp),
        x = "FPC1 (scaled)",
        y = "Marginal effect on ΔTmax (°C)"
      ) +
      ggplot2::theme_bw(base_size = 11)
  }

  pA <- make_panel(pred_w,  stats_w,  amp_w,  "With Hmax")
  pB <- make_panel(pred_wo, stats_wo, amp_wo, "Without Hmax")

  fig <- (pA | pB) +
    patchwork::plot_annotation(
      title   = "FPC1 marginal effect: with vs without Hmax in the model",
      caption = "ρ(FPC1, Hmax) = 0.65 | GAMM concurvity s(FPC1_sc)–s(Hmax_sc) = 0.56"
    )

  pred_w$model  <- "with_Hmax"
  pred_wo$model <- "without_Hmax"
  df_smooth <- rbind(pred_w, pred_wo)

  invisible(list(plot = fig, df_smooth = df_smooth, amp_with = amp_w, amp_without = amp_wo))
}

#' Fit the NONE GAMM without s(fCover_sc) — robustness test for LAI/fCover collinearity.
#'
#' Identical to fit_reference_gamm(shape_type="NONE") except s(fCover_sc) is dropped.
#' Structural terms kept: s(LAI_sc) + s(Hmax_sc) + climate smooths + RE.
#' Same scat() family and fREML method.
#'
#' @param df_gamm   Dataframe from prepare_gamm_data.
#' @param clim_vars ERA5 base names to include as climate smooths.
#' @return Fitted bam model object.
fit_gamm_without_fcover <- function(df_gamm,
                                     clim_vars = c("Tmax_macro", "Wind_mean", "Rad_mean")) {
  base_terms <- c(
    safe_smooth_term("LAI_sc",  7L, df_gamm),
    safe_smooth_term("Hmax_sc", 7L, df_gamm)
    # s(fCover_sc) deliberately omitted — sensitivity test
  )

  clim_terms <- if (length(clim_vars) > 0) {
    cat(sprintf("  [GAMM no-fCover] Climate drivers: %s\n", paste(clim_vars, collapse = ", ")))
    vapply(paste0(clim_vars, "_sc"), function(v) safe_smooth_term(v, 7L, df_gamm),
           character(1))
  } else {
    character(0)
  }

  re_terms <- c("s(date_factor, bs='re')", "s(plot_id, bs='re')")
  rhs  <- paste(c(base_terms, clim_terms, re_terms), collapse = " + ")
  form <- as.formula(paste("Delta_Tmax ~", rhs))
  cat("  [GAMM no-fCover] Fitting:", deparse(form), "\n")
  bam(form, data = df_gamm, family = scat(), method = "fREML", discrete = FALSE)
}

#' Two-panel comparison of LAI marginal effect with vs without fCover in the model.
#'
#' Extracts the s(LAI_sc) smooth from both models, plots side by side on a shared
#' y-axis, and returns smooth data and amplitudes for export.
#'
#' @param gam_with_fcover    Fitted bam from fit_reference_gamm(shape_type="NONE").
#' @param gam_without_fcover Fitted bam from fit_gamm_without_fcover().
#' @return Invisible list: $plot, $df_smooth, $amp_with (°C), $amp_without (°C).
plot_lai_robustness_to_fcover <- function(gam_with_fcover, gam_without_fcover) {
  get_lai_stats <- function(m) {
    st  <- summary(m)$s.table
    idx <- grep("LAI_sc", rownames(st))
    list(edf  = round(st[idx, "edf"], 2),
         pval = st[idx, "p-value"])
  }

  stats_w  <- get_lai_stats(gam_with_fcover)
  stats_wo <- get_lai_stats(gam_without_fcover)

  pred_w  <- as.data.frame(ggpredict(gam_with_fcover,    terms = "LAI_sc"))
  pred_wo <- as.data.frame(ggpredict(gam_without_fcover, terms = "LAI_sc"))

  amp_w  <- max(pred_w$predicted)  - min(pred_w$predicted)
  amp_wo <- max(pred_wo$predicted) - min(pred_wo$predicted)

  y_lim <- range(c(pred_w$conf.low,  pred_w$conf.high,
                   pred_wo$conf.low, pred_wo$conf.high))

  make_panel <- function(pred, stats, amp, panel_title) {
    pval_str <- if (stats$pval < 0.001) "p<0.001" else sprintf("p=%.3f", stats$pval)
    ggplot2::ggplot(pred, ggplot2::aes(x = x, y = predicted)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                           alpha = 0.2, fill = "#31688e") +
      ggplot2::geom_line(colour = "#31688e", linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      ggplot2::coord_cartesian(ylim = y_lim) +
      ggplot2::labs(
        title    = panel_title,
        subtitle = sprintf("edf = %.2f | %s | amplitude = %.3f°C",
                           stats$edf, pval_str, amp),
        x = "LAI (scaled)",
        y = "Marginal effect on ΔTmax (°C)"
      ) +
      ggplot2::theme_bw(base_size = 11)
  }

  pA <- make_panel(pred_w,  stats_w,  amp_w,  "With fCover")
  pB <- make_panel(pred_wo, stats_wo, amp_wo, "Without fCover")

  fig <- (pA | pB) +
    patchwork::plot_annotation(
      title   = "LAI marginal effect: with vs without fCover in the model",
      caption = "ρ(LAI, fCover) = 0.842 | GAMM concurvity s(LAI_sc)–s(fCover_sc) = 0.842"
    )

  pred_w$model  <- "with_fCover"
  pred_wo$model <- "without_fCover"
  df_smooth <- rbind(pred_w, pred_wo)

  invisible(list(plot = fig, df_smooth = df_smooth, amp_with = amp_w, amp_without = amp_wo))
}

#' Diagnostic de concurvité pour le GAMM de référence.
#'
#' Teste la concurvité globale (full=TRUE) et paire-à-paire (full=FALSE) pour
#' les termes de forçage macroclimatique (Tmax, Rad, VPD, Wind) physiquement
#' corrélés.  Lève une alerte console si le worst pairwise > 0.8 sur ce groupe.
#' Sauvegarde deux CSV dans out_dir.
#'
#' @param gam_model  Modèle bam/gam ajusté (sorti de fit_reference_gamm).
#' @param out_dir    Répertoire de sortie pour les CSV.
#' @return Invisible list $full et $pairwise.
run_concurvity_audit <- function(gam_model, out_dir = "outputs/audit") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  conc_full <- mgcv::concurvity(gam_model, full = TRUE)
  conc_pair <- mgcv::concurvity(gam_model, full = FALSE)

  cat("\n── Concurvity audit (full = TRUE) — worst-case par terme smooth ──\n")
  print(round(conc_full, 3))

  cat("\n── Concurvity audit (full = FALSE) — matrice pairwise worst ──\n")
  worst_mat <- conc_pair[["worst"]]
  print(round(worst_mat, 3))

  clim_terms  <- c("s(Tmax_macro_sc)", "s(Rad_mean_sc)", "s(VPD_mean_sc)", "s(Wind_mean_sc)")
  terms_found <- intersect(clim_terms, rownames(worst_mat))

  if (length(terms_found) >= 2) {
    sub_mat      <- worst_mat[terms_found, terms_found, drop = FALSE]
    diag(sub_mat) <- NA
    max_conc     <- max(sub_mat, na.rm = TRUE)

    if (max_conc > 0.8) {
      cat(sprintf(
        "\n[!] CONCURVITE ELEVEE : worst pairwise macroclimat = %.3f > 0.8\n",
        max_conc
      ))
      cat("    Termes concernés :", paste(terms_found, collapse = ", "), "\n")
      cat("    → Envisager de retirer un terme redondant ou d'utiliser",
          "un tenseur produit ti() pour contrôler la collinéarité.\n\n")
    } else {
      cat(sprintf(
        "\n[ok] Concurvité acceptable : worst pairwise macroclimat = %.3f ≤ 0.8\n\n",
        max_conc
      ))
    }
  }

  write.csv(as.data.frame(conc_full),
            file.path(out_dir, "concurvity_full.csv"), row.names = TRUE)
  write.csv(as.data.frame(worst_mat),
            file.path(out_dir, "concurvity_pairwise_worst.csv"), row.names = TRUE)
  cat(sprintf("  [concurvity] CSV sauvegardés → %s/\n", out_dir))

  invisible(list(full = conc_full, pairwise = conc_pair))
}
