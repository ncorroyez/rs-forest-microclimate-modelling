# ==============================================================================
# Chapter 1 — Unit tests for pipeline invariants
# ==============================================================================

#' Test LAI conservation across all make_lad_* constructors.
#'
#' @param df_sample cLHS sample dataframe.
#' @param eps       Absolute tolerance on sum(density) vs LAI (default 1e-6).
#' @return Invisible logical: TRUE if all checks pass.
test_lad_lai_conservation <- function(df_sample, eps = 1e-6) {
  cat("── [TEST] LAI conservation in make_lad_* ──\n")
  n_test <- min(20, nrow(df_sample))
  idx    <- round(seq(1, nrow(df_sample), length.out = n_test))

  lad_mean_fn <- make_lad_mean_factory(df_sample)
  lad_ct_fn   <- make_lad_cluster_type_factory(df_sample)

  fns <- list(
    real         = make_lad_real,
    uniform      = make_lad_uniform,
    mean_shape   = lad_mean_fn,
    cluster_type = lad_ct_fn
  )

  ok_all <- TRUE
  for (fn_name in names(fns)) {
    fn <- fns[[fn_name]]
    for (i in idx) {
      row <- df_sample[i, ]
      lai <- as.numeric(row$LAI)
      # Pass lai explicitly so lai_tgt is well-defined — tests conservation of
      # whatever value the pipeline will actually pass (plot_row$LAI via lai_fn).
      prof <- tryCatch(fn(row, lai = lai), error = function(e) {
        cat(sprintf("  [FAIL] %s plot %d error: %s\n", fn_name, i, e$message))
        NULL
      })
      if (is.null(prof)) { ok_all <- FALSE; next }
      s <- sum(prof$density, na.rm = TRUE)
      if (abs(s - lai) > eps) {
        cat(sprintf("  [FAIL] %s plot %d: sum=%.8f  LAI=%.8f  diff=%.2e\n",
                    fn_name, i, s, lai, abs(s - lai)))
        ok_all <- FALSE
      }
    }
  }

  if (ok_all) cat(sprintf("  [PASS] All %d LAD constructors × %d plots conserve LAI (tol=%.0e)\n",
                           length(fns), n_test, eps))
  invisible(ok_all)
}

#' Test Shapley efficiency axiom: Σφ_i = v(N) - v(∅) = total_effect.
#'
#' @param shap_res Result of compute_shapley_exact().
#' @param eps      Absolute tolerance in °C (default 1e-6).
#' @return Invisible logical.
test_shapley_efficiency <- function(shap_res, eps = 1e-6) {
  cat("── [TEST] Shapley efficiency axiom (Σφ = v(N)) ──\n")
  diff <- abs(shap_res$check_diff)
  if (is.na(diff)) {
    cat("  [FAIL] check_diff is NA — Shapley computation produced NA values (check coalition mapping)\n")
    return(invisible(FALSE))
  }
  if (diff < eps) {
    cat(sprintf("  [PASS] |Σφ - v(N)| = %.2e°C  (tol=%.0e)\n", diff, eps))
    invisible(TRUE)
  } else {
    cat(sprintf("  [FAIL] |Σφ - v(N)| = %.2e°C  ≥ tol %.0e\n", diff, eps))
    invisible(FALSE)
  }
}

#' Test REF_all_real ↔ H2_real_LAD identity.
#'
#' Both scenarios use identical full-LiDAR inputs; their ΔTmax must match
#' to floating-point precision. Stops the pipeline (via stop()) if violated.
#'
#' @param df_all_scenarios Stacked scenario results.
#' @param df_h2            H2 scenario results (from extract_all_scenarios).
#' @param ref_name         Reference scenario name (default "REF_all_real").
#' @param eps              RMSE threshold in °C (default 1e-6).
#' @return Invisible logical; calls stop() on failure.
test_ref_h2_reconciliation <- function(df_all_scenarios, df_h2,
                                        ref_name = "REF_all_real",
                                        eps      = 1e-6) {
  cat("── [TEST] REF_all_real ↔ H2_real_LAD reconciliation ──\n")

  df_ref  <- df_all_scenarios %>%
    filter(scenario == ref_name) %>%
    dplyr::select(x, y, date, Delta_REF = Delta_Tmax)
  df_real <- df_h2 %>%
    filter(scenario == "H2_real_LAD") %>%
    dplyr::select(x, y, date, Delta_H2 = Delta_Tmax)

  df_chk <- inner_join(df_ref, df_real, by = c("x", "y", "date"))

  if (nrow(df_chk) == 0) {
    warning("[TEST] REF/H2_real reconciliation: no overlapping rows — check pipeline order")
    return(invisible(FALSE))
  }

  rmse_chk <- sqrt(mean((df_chk$Delta_REF - df_chk$Delta_H2)^2, na.rm = TRUE))

  if (rmse_chk < eps) {
    cat(sprintf("  [PASS] RMSE(REF, H2_real) = %.2e°C  (n=%d, tol=%.0e)\n",
                rmse_chk, nrow(df_chk), eps))
    invisible(TRUE)
  } else {
    stop(sprintf(
      "[TEST FAIL] REF_all_real ↔ H2_real_LAD RMSE = %.6f°C (> %.0e). Pipeline broken — scenario definitions diverged.",
      rmse_chk, eps
    ))
  }
}

#' Run all pipeline unit tests.
#'
#' Calls available tests depending on which objects are provided.
#'
#' @param df_sample  cLHS sample dataframe (required for LAD tests).
#' @param shap_res   compute_shapley_exact() result, or NULL.
#' @param df_all     df_all_scenarios, or NULL.
#' @param df_h2      H2 scenario results, or NULL.
#' @param ref_name   Reference scenario name.
#' @return Invisible logical: TRUE if all run tests pass.
run_unit_tests <- function(df_sample,
                            shap_res = NULL,
                            df_all   = NULL,
                            df_h2    = NULL,
                            ref_name = "REF_all_real") {
  cat("\n══════════════════════════════════════════════\n")
  cat("  UNIT TESTS — pipeline invariants\n")
  cat("══════════════════════════════════════════════\n")

  results <- list()
  results[["lad_lai"]] <- test_lad_lai_conservation(df_sample)

  if (!is.null(shap_res))
    results[["shapley"]] <- test_shapley_efficiency(shap_res)

  if (!is.null(df_all) && !is.null(df_h2))
    results[["ref_h2"]] <- tryCatch(
      test_ref_h2_reconciliation(df_all, df_h2, ref_name),
      error = function(e) { cat(sprintf("  %s\n", e$message)); FALSE }
    )

  results_vec <- unlist(results)
  cat("──────────────────────────────────────────────\n")
  cat(sprintf("  %d / %d tests passed\n", sum(results_vec), length(results_vec)))
  cat("══════════════════════════════════════════════\n\n")
  invisible(all(results_vec))
}
