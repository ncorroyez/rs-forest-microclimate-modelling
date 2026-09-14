# ==============================================================================
# LOVB analysis — Module 5 : tables (CSV + markdown to console)
#
# Five tables :
#   Tab 1 — Global cLHS metrics : 4 vars x 2 metrics (RMSE, Bias, pct>0.5C, rank)
#   Tab 2 — Per-cluster cLHS metrics : 4 clusters x 4 vars x 2 metrics
#   Tab 3 — Ranking convergence (Shapley vs LOVB vs LVA on cLHS), Kendall tau
#   Tab 4 — Raw archetype values : 4 arch x 4 vars x 2 metrics (REF, LOVB, Delta)
#   Tab 5 — Per-archetype ranking convergence + Kendall tau
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
})

.VAR_ORDER <- c("LAI", "fCover", "Hmax", "LAD")

# ---- helper : ranking by absolute value (or RMSE), 1 = largest -----------------
.lovb_rank <- function(x) rank(-abs(x), ties.method = "average")

# Pretty markdown printer for data.table
.lovb_print_md <- function(DT, caption = NULL) {
  if (!is.null(caption)) cat(sprintf("\n### %s\n\n", caption))
  if (nrow(DT) == 0) { cat("(empty)\n"); return(invisible(NULL)) }
  # Header
  cat("|", paste(names(DT), collapse = " | "), "|\n")
  cat("|", paste(rep("---", ncol(DT)), collapse = " | "), "|\n")
  for (i in seq_len(nrow(DT))) {
    vals <- vapply(DT[i], function(v) {
      if (is.numeric(v)) sprintf("%.3f", v) else as.character(v)
    }, character(1))
    cat("|", paste(vals, collapse = " | "), "|\n")
  }
  invisible(NULL)
}

# ==============================================================================
# Tab 1 — Global cLHS metrics
# ==============================================================================
lovb_tab1_global <- function(DT_contrib_clhs,
                              out_csv = here("outputs/lovb/tables/tab_LOVB_global.csv")) {
  stopifnot("x" %in% names(DT_contrib_clhs))
  rows <- list()
  for (m in c("mean", "P90")) {
    ref_col <- paste0("Tmax_", m, "_REF")
    for (v in .VAR_ORDER) {
      lovb_col <- paste0("Tmax_", m, "_LOVB_", v)
      d <- DT_contrib_clhs[[ref_col]] - DT_contrib_clhs[[lovb_col]]
      rows[[length(rows) + 1L]] <- data.table(
        Variable     = v,
        Metric       = m,
        RMSE_global  = sqrt(mean(d^2, na.rm = TRUE)),
        Bias_global  = mean(d, na.rm = TRUE),
        pct_above_0.5C = 100 * mean(abs(d) > 0.5, na.rm = TRUE),
        N            = sum(!is.na(d))
      )
    }
  }
  DT <- rbindlist(rows)
  DT[, Ranking := .lovb_rank(RMSE_global), by = Metric]
  setorder(DT, Metric, Ranking)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT, out_csv)
  cli_alert_success("Saved {.path {out_csv}}")
  .lovb_print_md(DT, "Tab 1 — Global cLHS LOVB metrics (sorted by Metric then RMSE rank)")
  DT
}

# ==============================================================================
# Tab 2 — Per-cluster cLHS metrics (32 rows)
# ==============================================================================
lovb_tab2_per_archetype <- function(DT_contrib_clhs,
                                     out_csv = here("outputs/lovb/tables/tab_LOVB_per_archetype.csv")) {
  stopifnot("Cluster" %in% names(DT_contrib_clhs))
  cluster_labs <- c("1" = "C1 Open", "2" = "C2 Dense",
                     "3" = "C3 Bas couvert", "4" = "C4 Inter")
  rows <- list()
  for (cl in c("1", "2", "3", "4")) {
    DT_cl <- DT_contrib_clhs[as.character(Cluster) == cl]
    if (nrow(DT_cl) == 0) next
    for (m in c("mean", "P90")) {
      ref_col <- paste0("Tmax_", m, "_REF")
      for (v in .VAR_ORDER) {
        lovb_col <- paste0("Tmax_", m, "_LOVB_", v)
        d <- DT_cl[[ref_col]] - DT_cl[[lovb_col]]
        rows[[length(rows) + 1L]] <- data.table(
          Archetype       = cluster_labs[cl],
          Variable        = v,
          Metric          = m,
          RMSE            = sqrt(mean(d^2, na.rm = TRUE)),
          Bias            = mean(d, na.rm = TRUE),
          pct_above_0.5C  = 100 * mean(abs(d) > 0.5, na.rm = TRUE),
          N               = sum(!is.na(d))
        )
      }
    }
  }
  DT <- rbindlist(rows)
  setorder(DT, Archetype, Metric, -RMSE)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT, out_csv)
  cli_alert_success("Saved {.path {out_csv}}")
  .lovb_print_md(DT, "Tab 2 — Per-cluster cLHS LOVB metrics (32 rows)")
  DT
}

# ==============================================================================
# Tab 3 — Ranking convergence cLHS (Shapley vs LOVB vs LVA)
# ==============================================================================
lovb_tab3_convergence <- function(DT_contrib_clhs,
                                   shapley_csv = here("outputs/figures/02_shapley_attribution.csv"),
                                   hobo_gam_csv = here("outputs/lovb/tables/tab_HOBO_convergence.csv"),
                                   out_csv     = here("outputs/lovb/tables/tab_convergence_rankings.csv")) {
  # Shapley rank (cLHS principal)
  df_sh <- fread(shapley_csv)
  shap_phi <- setNames(df_sh$shapley, df_sh$variable)
  rank_shap <- .lovb_rank(shap_phi[.VAR_ORDER])
  names(rank_shap) <- .VAR_ORDER

  # LOVB ranks (mean + P90) — by global RMSE
  ranks_lovb <- list()
  for (m in c("mean", "P90")) {
    ref_col <- paste0("Tmax_", m, "_REF")
    rmses <- vapply(.VAR_ORDER, function(v) {
      d <- DT_contrib_clhs[[ref_col]] - DT_contrib_clhs[[paste0("Tmax_", m, "_LOVB_", v)]]
      sqrt(mean(d^2, na.rm = TRUE))
    }, numeric(1))
    ranks_lovb[[m]] <- rank(-rmses, ties.method = "average")
  }

  # LVA ranks (mean + P90) — by global RMSE between LVA_v and NULL
  ranks_lva <- list()
  for (m in c("mean", "P90")) {
    null_col <- paste0("Tmax_", m, "_NULL")
    rmses <- vapply(.VAR_ORDER, function(v) {
      d <- DT_contrib_clhs[[paste0("Tmax_", m, "_LVA_", v)]] - DT_contrib_clhs[[null_col]]
      sqrt(mean(d^2, na.rm = TRUE))
    }, numeric(1))
    ranks_lva[[m]] <- rank(-rmses, ties.method = "average")
  }

  # GAM observed rank (multivariate scat on 53 HOBO sites) — partial_R2
  rank_gam_obs <- setNames(rep(NA_real_, length(.VAR_ORDER)), .VAR_ORDER)
  if (file.exists(hobo_gam_csv)) {
    df_gam <- fread(hobo_gam_csv)
    if ("Rank_obs" %in% names(df_gam) && all(.VAR_ORDER %in% df_gam$Variable)) {
      rank_gam_obs <- setNames(df_gam$Rank_obs[match(.VAR_ORDER, df_gam$Variable)],
                                .VAR_ORDER)
      cli_alert("HOBO GAM ranks loaded from {.path {hobo_gam_csv}}")
    }
  } else {
    cli_alert_warning("No HOBO GAM convergence file — Rank_GAM_obs will be NA")
  }

  DT <- data.table(
    Variable          = .VAR_ORDER,
    Rank_Shapley      = rank_shap,
    Rank_LOVB_mean    = ranks_lovb$mean,
    Rank_LOVB_P90     = ranks_lovb$P90,
    Rank_LVA_mean     = ranks_lva$mean,
    Rank_LVA_P90      = ranks_lva$P90,
    Rank_GAM_obs_HOBO = rank_gam_obs
  )

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT, out_csv)
  cli_alert_success("Saved {.path {out_csv}}")
  .lovb_print_md(DT, "Tab 3 — Convergence rankings (cLHS sim + HOBO obs GAM)")

  # Kendall tau matrix (NA if GAM_obs missing)
  rank_cols <- c("Rank_Shapley", "Rank_LOVB_mean", "Rank_LOVB_P90",
                  "Rank_LVA_mean",  "Rank_LVA_P90",
                  "Rank_GAM_obs_HOBO")
  K <- matrix(NA_real_, length(rank_cols), length(rank_cols),
              dimnames = list(rank_cols, rank_cols))
  for (i in seq_along(rank_cols)) for (j in seq_along(rank_cols)) {
    a <- DT[[rank_cols[i]]]; b <- DT[[rank_cols[j]]]
    K[i, j] <- if (all(is.na(a)) || all(is.na(b))) NA_real_
                else cor(a, b, method = "kendall", use = "complete.obs")
  }
  cat("\n### Kendall tau — Shapley / LOVB / LVA / GAM obs HOBO (all rankings)\n\n")
  print(round(K, 3))

  invisible(list(table = DT, kendall = K))
}

# ==============================================================================
# Tab 4 — Raw archetype values (32 rows : 4 arch x 4 var x 2 metrics)
# ==============================================================================
lovb_tab4_archetypes <- function(DT_contrib_archetypes,
                                  out_csv = here("outputs/lovb/tables/tab_LOVB_archetypes.csv")) {
  stopifnot("archetype" %in% names(DT_contrib_archetypes))
  rows <- list()
  for (i in seq_len(nrow(DT_contrib_archetypes))) {
    arch <- DT_contrib_archetypes$archetype[i]
    for (m in c("mean", "P90")) {
      ref <- DT_contrib_archetypes[[paste0("Tmax_", m, "_REF")]][i]
      for (v in .VAR_ORDER) {
        lovb_val <- DT_contrib_archetypes[[paste0("Tmax_", m, "_LOVB_", v)]][i]
        rows[[length(rows) + 1L]] <- data.table(
          Archetype       = arch,
          Variable        = v,
          Metric          = m,
          DeltaTmax_REF   = ref,
          DeltaTmax_LOVB  = lovb_val,
          Delta_v         = ref - lovb_val
        )
      }
    }
  }
  DT <- rbindlist(rows)
  setorder(DT, Archetype, Metric, -Delta_v)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT, out_csv)
  cli_alert_success("Saved {.path {out_csv}}")
  .lovb_print_md(DT, "Tab 4 — Raw archetype LOVB values (32 rows)")
  DT
}

# ==============================================================================
# Tab 5 — Per-archetype convergence Shapley vs LOVB (mean, P90)
# ==============================================================================
lovb_tab5_convergence_archetypes <- function(DT_contrib_archetypes,
                                              shapley_arch_csv = here("outputs/h1/shapley_archetypes_values.csv"),
                                              out_csv = here("outputs/lovb/tables/tab_convergence_archetypes.csv")) {
  if (!file.exists(shapley_arch_csv))
    stop("Missing Shapley archetypes CSV : ", shapley_arch_csv)
  df_sh <- fread(shapley_arch_csv)   # archetype, variable, phi, pct, ...

  archs <- unique(DT_contrib_archetypes$archetype)
  rows <- list()
  for (arch in archs) {
    # Shapley ranks per archetype
    sh_a <- df_sh[archetype == arch & variable %in% .VAR_ORDER]
    phi  <- setNames(sh_a$phi, sh_a$variable)
    if (length(phi) < 4L) next
    rank_shap <- .lovb_rank(phi[.VAR_ORDER])

    # LOVB ranks per archetype
    DT_a <- DT_contrib_archetypes[archetype == arch]
    rank_lovb <- list()
    for (m in c("mean", "P90")) {
      ref <- DT_a[[paste0("Tmax_", m, "_REF")]]
      deltas <- vapply(.VAR_ORDER, function(v) {
        ref - DT_a[[paste0("Tmax_", m, "_LOVB_", v)]]
      }, numeric(1))
      rank_lovb[[m]] <- .lovb_rank(deltas)
    }

    tau_mean <- cor(rank_shap, rank_lovb$mean, method = "kendall")
    tau_P90  <- cor(rank_shap, rank_lovb$P90,  method = "kendall")

    for (v in .VAR_ORDER) {
      rows[[length(rows) + 1L]] <- data.table(
        Archetype       = arch,
        Variable        = v,
        Rank_Shapley    = rank_shap[v],
        Rank_LOVB_mean  = rank_lovb$mean[v],
        Rank_LOVB_P90   = rank_lovb$P90[v],
        Kendall_tau_mean = tau_mean,
        Kendall_tau_P90  = tau_P90
      )
    }
  }
  DT <- rbindlist(rows)

  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  fwrite(DT, out_csv)
  cli_alert_success("Saved {.path {out_csv}}")
  .lovb_print_md(DT, "Tab 5 — Per-archetype ranking convergence")

  cat("\n### Kendall tau summary per archetype\n\n")
  tau_summary <- unique(DT[, .(Archetype,
                                Kendall_tau_Shapley_vs_LOVBmean = Kendall_tau_mean,
                                Kendall_tau_Shapley_vs_LOVBP90  = Kendall_tau_P90)])
  .lovb_print_md(tau_summary)

  invisible(list(table = DT, tau_summary = tau_summary))
}
