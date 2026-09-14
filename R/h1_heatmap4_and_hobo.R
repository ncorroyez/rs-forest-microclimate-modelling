# ==============================================================================
# H1 MEB deck — simplified 4-column heatmap + non-parametric HOBO validation
#
# Replaces the previous 6-column heatmap (Shapley | LOVB mean | LOVB P90 |
# LVA mean | LVA P90 | GAM HOBO) by a 4-column triangulation
# (Shapley | LOVB mean | LVA mean | Consensus) and validates the
# observational ranking via :
#   (A) Spearman correlations + bootstrap CI95 between observed slope and
#       simulated LOVB delta per HOBO sensor
#   (B) Jonckheere-Terpstra monotonic trend test on slope ~ binned trait
#
# No GAM. No dplyr / purrr / tidyr.
# Outputs under outputs/lovb/ (per user choice).
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
})

HOBO_CACHE_PATH <- here::here("outputs/lovb/data/DT_HOBO_scalars.rds")

# ---- Helper : compute slope_obs + LAD_L2 + traits at 53 HOBO sensors ---------
# Calls lovb_observed_buffering() to get slope + traits + FPC1, then computes
# LAD_L2 = ||lad_profile - uniform_profile||_2 per sensor from raw LAD layers.
load_hobo_scalars <- function(force = FALSE) {
  if (!force && file.exists(HOBO_CACHE_PATH)) {
    cli_alert("Reusing cached HOBO scalars : {.path {HOBO_CACHE_PATH}}")
    return(readRDS(HOBO_CACHE_PATH))
  }
  cli_h2("Compute slope_obs + LAD_L2 at 53 HOBO sensors")

  # Source heavy helpers only here (not at file top) to keep step 1 cheap
  suppressMessages({
    library(ncdf4); library(tidyverse); library(sf); library(terra)
    library(musica.tools); library(rmusica)
  })
  source(here::here("R/config.R"))
  source(here::here("R/io.R"))
  source(here::here("R/lad.R"))
  source(here::here("R/forest.R"))
  source(here::here("R/scenarios.R"))
  source(here::here("R/musica.R"))
  source(here::here("R/validation.R"))
  source(here::here("R/fpca.R"))
  source(here::here("R/lovb_07_hobo_analysis.R"))

  DT_obs <- as.data.table(lovb_observed_buffering())   # id_plot, slope, intercept, LAI, Hmax, fCover, FPC1
  setnames(DT_obs, "slope", "slope_obs")

  # LAD_L2 : need raw LAD profile  ->  reload inputs
  rasters        <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
  df_hobo_inputs <- as.data.table(build_hobo_inputs(CFG$hobo_geojson, rasters$stack,
                                                      CFG$ids_to_remove))
  lad_cols <- grep("^LAD_Layer_", names(df_hobo_inputs), value = TRUE)
  z_breaks <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  dz <- median(diff(sort(unique(z_breaks))))   # uniform layer thickness

  l2_uniform <- function(lad_vec, hmax, z, dz) {
    lad_vec[is.na(lad_vec)] <- 0
    in_canopy <- z <= hmax
    n_in <- sum(in_canopy)
    if (n_in == 0L) return(NA_real_)
    lai_tot <- sum(lad_vec) * dz
    unif <- rep(0, length(z))
    unif[in_canopy] <- lai_tot / (n_in * dz)
    sqrt(sum((lad_vec - unif)^2) * dz)
  }

  lad_mat <- as.matrix(df_hobo_inputs[, ..lad_cols])
  hmax_v  <- as.numeric(df_hobo_inputs$Hmax)
  lad_l2  <- vapply(seq_len(nrow(lad_mat)),
                     function(i) l2_uniform(lad_mat[i, ], hmax_v[i], z_breaks, dz),
                     numeric(1L))
  DT_l2 <- data.table(id_plot = df_hobo_inputs$id_plot, LAD_L2 = lad_l2)

  DT_out <- merge(DT_obs, DT_l2, by = "id_plot")
  saveRDS(DT_out, HOBO_CACHE_PATH)
  cli_alert_success("Saved {.path {HOBO_CACHE_PATH}}")
  DT_out[]
}

# Step 1 : heatmap 4 cols + consensus ------------------------------------------
build_heatmap_4col <- function() {
  cli_h1("Step 1 — heatmap 4 cols + consensus")
  DT <- fread(here::here("outputs/lovb/tables/tab_convergence_rankings.csv"))

  # Keep only Shapley | LOVB_mean | LVA_mean columns ; compute consensus
  methods <- c("Rank_Shapley", "Rank_LOVB_mean", "Rank_LVA_mean")
  labels  <- c("Shapley",      "LOVB (mean)",    "LVA (mean)")
  DT[, mean_rank := rowMeans(.SD), .SDcols = methods]
  DT[, consensus_rank := frank(mean_rank, ties.method = "average")]
  setorder(DT, consensus_rank)
  cli_alert("Variable order by consensus :")
  print(DT[, .(Variable, Rank_Shapley, Rank_LOVB_mean, Rank_LVA_mean,
                mean_rank = round(mean_rank, 2), consensus_rank)])

  # Long format
  long <- melt(DT, id.vars = "Variable", measure.vars = methods,
                variable.name = "method_raw", value.name = "Rank")
  long[, method := factor(labels[match(method_raw, methods)],
                           levels = c(labels, "Consensus\n(mean rank)"))]
  long[, Rank_int := Rank]
  cons_long <- DT[, .(Variable    = Variable,
                       method      = factor("Consensus\n(mean rank)",
                                            levels = c(labels, "Consensus\n(mean rank)")),
                       Rank        = round(mean_rank, 1),
                       Rank_int    = consensus_rank)]
  all_long <- rbind(long[, .(Variable, method, Rank, Rank_int)], cons_long)
  all_long[, Variable := factor(Variable, levels = DT$Variable)]
  all_long[, fill_rank := as.character(round(Rank_int))]

  pal <- c("1" = "#440154", "2" = "#3B528B",
            "3" = "#5DC863", "4" = "#FDE725")

  p <- ggplot(all_long, aes(x = method, y = Variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = ifelse(grepl("Consensus", method),
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(all_long$fill_rank) <= 2, "white", "grey15"),
               fontface = "bold", size = 6) +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_y_discrete(limits = rev) +
    labs(
      title    = "Attribution ranking (4-method triangulation)",
      subtitle = "Dark = rank 1 (strongest attributed effect)   |   light = rank 4 (weakest)",
      x = NULL, y = NULL,
      caption = paste0(
        "Shapley / LOVB / LVA : MuSICA simulations (400 cLHS plots x 16 coalitions). ",
        "Mean over 122 summer days.\n",
        "P90 rankings invariant to substitution (see supplementary tab_convergence_rankings)."
      )
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 15),
           plot.subtitle = element_text(colour = "grey25", size = 11),
           plot.caption  = element_text(size = 8, colour = "grey35", hjust = 0,
                                         lineheight = 1.3),
           axis.text.x   = element_text(face = "bold", size = 11),
           axis.text.y   = element_text(face = "bold", size = 13),
           panel.grid    = element_blank())

  fig_path <- here::here("outputs/lovb/figures/fig_heatmap_4col_simulation.png")
  ggsave(fig_path, p, width = 9.5, height = 4.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")

  # Consensus table CSV
  tab_out <- DT[, .(Variable,
                     Shapley_rank   = Rank_Shapley,
                     LOVB_rank      = Rank_LOVB_mean,
                     LVA_rank       = Rank_LVA_mean,
                     Consensus_mean = round(mean_rank, 2),
                     Consensus_rank = consensus_rank)]
  tab_path <- here::here("outputs/lovb/tables/tab_consensus_4methods.csv")
  fwrite(tab_out, tab_path)
  cli_alert_success("Saved {.path {tab_path}}")
  invisible(list(fig = fig_path, tab = tab_path, data = tab_out))
}

# ---- Step 2A : Spearman + bootstrap CI95 --------------------------------------
# x = slope_obs (53 HOBO sensors)
# y = Delta_v_LOVB_mean simulated at HOBO (from DT_contrib_HOBO)
spearman_boot_ci <- function(x, y, B = 1000L, seed = 42L) {
  ok <- !is.na(x) & !is.na(y)
  x  <- x[ok]; y <- y[ok]
  n  <- length(x)
  rho_hat <- suppressWarnings(cor(x, y, method = "spearman"))
  p_val   <- suppressWarnings(cor.test(x, y, method = "spearman",
                                          exact = FALSE)$p.value)
  set.seed(seed)
  rho_boot <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    rho_boot[b] <- suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  }
  ci <- quantile(rho_boot, c(0.025, 0.975), na.rm = TRUE, type = 7)
  list(rho = rho_hat, p = p_val, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
}

run_spearman_validation <- function() {
  cli_h1("Step 2A — Spearman correlations slope_obs vs simulated LOVB delta")
  DT_h    <- load_hobo_scalars()
  DT_sim  <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  DT      <- merge(DT_h[, .(id_plot, slope_obs, LAI, Hmax, fCover, FPC1, LAD_L2)],
                    DT_sim[, .(id_plot,
                                Delta_LAI_mean, Delta_Hmax_mean,
                                Delta_fCover_mean, Delta_LAD_mean)],
                    by = "id_plot")
  cli_alert("n = {nrow(DT)} HOBO sensors with slope_obs + simulated Delta_v")

  vars   <- c("LAI", "Hmax", "fCover", "LAD")
  ycols  <- c("Delta_LAI_mean", "Delta_Hmax_mean",
                "Delta_fCover_mean", "Delta_LAD_mean")
  rows   <- vector("list", length(vars))
  for (i in seq_along(vars)) {
    r <- spearman_boot_ci(DT$slope_obs, DT[[ycols[i]]])
    rows[[i]] <- data.table(Variable = vars[i],
                              rho_S    = r$rho,
                              p_value  = r$p,
                              ci_lo    = r$ci_lo,
                              ci_hi    = r$ci_hi)
  }
  tab <- rbindlist(rows)
  cli_alert("Spearman + bootstrap CI95 (1000 reps) :")
  print(tab)

  # Save table
  tab_path <- here::here("outputs/lovb/tables/tab_HOBO_spearman.csv")
  fwrite(tab, tab_path)
  cli_alert_success("Saved {.path {tab_path}}")

  # Long data for plotting
  long <- rbindlist(lapply(seq_along(vars), function(i)
    data.table(Variable = vars[i],
                slope_obs = DT$slope_obs,
                Delta_v   = DT[[ycols[i]]])))
  long[, Variable := factor(Variable, levels = vars)]

  ann <- copy(tab)
  ann[, label := sprintf("rho = %+.2f  p = %s\n[%+.2f ; %+.2f]",
                           rho_S,
                           ifelse(p_value < 1e-3, "<0.001",
                                   sprintf("%.3f", p_value)),
                           ci_lo, ci_hi)]
  ann[, Variable := factor(Variable, levels = vars)]

  p <- ggplot(long, aes(x = slope_obs, y = Delta_v)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
    geom_point(colour = "#2C5F2D", alpha = 0.65, size = 1.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                  colour = "#FFB400", linewidth = 0.8, alpha = 0.7) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.2, size = 3.1, fontface = "bold",
               inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = "HOBO validation : observed slope vs simulated LOVB delta",
      subtitle = "n = 53 sensors  |  Spearman rho + bootstrap CI95 (1000 reps)",
      x = "Observed buffering slope (Tmax_obs ~ Tmax_macro)",
      y = "Simulated LOVB delta (K, mean over 122 d)",
      caption = "Non-parametric rank correlation. LM line shown for visual reference only."
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 11),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
           panel.grid.minor = element_blank())
  fig_path <- here::here("outputs/lovb/figures/fig_HOBO_validation_spearman.png")
  ggsave(fig_path, p, width = 11.4, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(tab)
}

# ---- Step 2B : Jonckheere-Terpstra trend test --------------------------------
# 3-bin (low/mid/high) by trait, JT test for monotonic trend on slope_obs.
# Manual implementation (DescTools unavailable). Asymptotic normal p-value.
jt_test <- function(x, group, alternative = c("two.sided", "increasing", "decreasing")) {
  alternative <- match.arg(alternative)
  group <- as.factor(group)
  g_levels <- levels(group)
  k <- length(g_levels)
  J <- 0
  for (i in seq_len(k - 1L)) {
    for (j in seq.int(i + 1L, k)) {
      x_i <- x[group == g_levels[i]]
      x_j <- x[group == g_levels[j]]
      for (a in x_i) {
        J <- J + sum(x_j > a) + 0.5 * sum(x_j == a)
      }
    }
  }
  n   <- length(x)
  n_g <- tabulate(group, nbins = k)
  mu  <- (n^2 - sum(n_g^2)) / 4
  sigma2 <- (n^2 * (2 * n + 3) - sum(n_g^2 * (2 * n_g + 3))) / 72
  z <- (J - mu) / sqrt(sigma2)
  p <- switch(alternative,
                "two.sided"  = 2 * (1 - pnorm(abs(z))),
                "increasing" = 1 - pnorm(z),
                "decreasing" = pnorm(z))
  list(statistic = J, z = z, p.value = p, mu = mu, var = sigma2,
        alternative = alternative)
}

run_jonckheere <- function() {
  cli_h1("Step 2B — Jonckheere-Terpstra monotonic trend test")
  DT <- load_hobo_scalars()
  vars <- c("LAI", "Hmax", "fCover", "LAD_FPC1", "LAD_L2")
  trait_cols <- c(LAI = "LAI", Hmax = "Hmax", fCover = "fCover",
                   LAD_FPC1 = "FPC1", LAD_L2 = "LAD_L2")
  alt_lookup <- c(LAI = "decreasing", Hmax = "decreasing",
                   fCover = "decreasing",
                   LAD_FPC1 = "two.sided", LAD_L2 = "two.sided")

  rows <- vector("list", length(vars)); bins_long <- list()
  for (v in vars) {
    tc <- trait_cols[v]
    trait <- DT[[tc]]
    slope <- DT$slope_obs
    keep  <- !is.na(trait) & !is.na(slope)
    trait <- trait[keep]; slope <- slope[keep]

    # 3 equal-N bins via rank (robust to ties, e.g. fCover saturated at 1)
    n   <- length(trait)
    rnk <- rank(trait, ties.method = "first")
    bin <- cut(rnk, breaks = c(0, n/3, 2*n/3, n),
                include.lowest = TRUE, labels = c("low", "mid", "high"))
    bin <- factor(bin, levels = c("low", "mid", "high"))

    res <- jt_test(slope, bin, alternative = alt_lookup[v])

    med  <- tapply(slope, bin, median)
    iqr  <- tapply(slope, bin, IQR)
    nbin <- tapply(slope, bin, length)

    trend <- if (res$p.value < 0.05) {
      if (res$z > 0) "Monotonic increase" else "Monotonic decrease"
    } else "No monotonic trend"

    rows[[length(rows) + 1L]] <- data.table(
      Variable      = v,
      Alternative   = alt_lookup[v],
      n_low         = nbin["low"],  n_mid = nbin["mid"], n_high = nbin["high"],
      med_low       = round(med["low"],  3),
      med_mid       = round(med["mid"],  3),
      med_high      = round(med["high"], 3),
      iqr_low       = round(iqr["low"],  3),
      iqr_mid       = round(iqr["mid"],  3),
      iqr_high      = round(iqr["high"], 3),
      JT_statistic  = res$statistic,
      JT_z          = round(res$z, 3),
      JT_pvalue     = signif(res$p.value, 4),
      Interpretation = trend
    )
    bins_long[[v]] <- data.table(Variable = v, bin = bin, slope = slope,
                                    p_value = res$p.value, alt = alt_lookup[v])
  }
  tab <- rbindlist(rows)
  cli_alert("JT results :")
  print(tab[, .(Variable, Alternative, JT_z, JT_pvalue, Interpretation)])
  tab_path <- here::here("outputs/lovb/tables/tab_HOBO_jonckheere.csv")
  fwrite(tab, tab_path)
  cli_alert_success("Saved {.path {tab_path}}")

  # Companion figure : boxplots by bin
  long <- rbindlist(bins_long)
  long[, Variable := factor(Variable, levels = vars)]
  ann <- tab[, .(Variable, label = sprintf("JT p = %s",
                                              ifelse(JT_pvalue < 1e-3, "<0.001",
                                                      sprintf("%.3f", JT_pvalue))))]
  ann[, Variable := factor(Variable, levels = vars)]

  p <- ggplot(long, aes(x = bin, y = slope)) +
    geom_boxplot(fill = "#2C5F2D", alpha = 0.55, colour = "grey20",
                  outlier.size = 1.3) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1, colour = "grey30") +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.4, size = 3.2, fontface = "bold",
               inherit.aes = FALSE) +
    facet_wrap(~ Variable, nrow = 1L) +
    labs(
      title    = "Observed buffering slope by trait bin (low / mid / high)",
      subtitle = "Jonckheere-Terpstra non-parametric monotonic trend test  |  n = 53 HOBO",
      x = NULL, y = "Observed slope (Tmax_obs ~ Tmax_macro)",
      caption = paste0("Alternative = decreasing for LAI/Hmax/fCover ",
                        "(physical expectation : trait ↑  ->  buffering ↑  ->  slope ↓).\n",
                        "Alternative = two.sided for LAD_FPC1 and LAD_L2 (no a priori direction).")
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold"),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0,
                                            lineheight = 1.3),
           panel.grid.minor = element_blank())
  fig_path <- here::here("outputs/lovb/figures/fig_HOBO_jonckheere_bins.png")
  ggsave(fig_path, p, width = 13, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(tab)
}

# ---- Step 3 : console synthesis ----------------------------------------------
print_synthesis <- function() {
  cli_h1("Step 3 — Synthesis")

  sim <- fread(here::here("outputs/lovb/tables/tab_consensus_4methods.csv"))
  sp  <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman.csv"))
  jt  <- fread(here::here("outputs/lovb/tables/tab_HOBO_jonckheere.csv"))

  cat("\n=====================================================\n")
  cat(" SIMULATION RANKING (4-method triangulation, n=400 cLHS)\n")
  cat("=====================================================\n")
  setorder(sim, Consensus_rank)
  cat(sprintf("%-8s  %-7s  %-5s  %-5s  %-9s\n",
                "Variable", "Shapley", "LOVB", "LVA", "Consensus"))
  for (i in seq_len(nrow(sim))) {
    r <- sim[i]
    cat(sprintf("%-8s  %-7d  %-5d  %-5d  %.2f  (#%d)\n",
                  r$Variable, r$Shapley_rank, r$LOVB_rank, r$LVA_rank,
                  r$Consensus_mean, r$Consensus_rank))
  }

  cat("\n=====================================================\n")
  cat(" OBSERVATIONAL VALIDATION (n=53 HOBO, non-parametric)\n")
  cat("=====================================================\n")
  cat(sprintf("%-10s  %-13s  %-12s  %-13s  %s\n",
                "Variable", "Spearman_rho", "p_value", "Jonckheere_p", "Bulk_effect"))

  # Map JT rows : LAI/Hmax/fCover from JT, LAD aggregated from JT LAD_FPC1 + LAD_L2
  vars_disp <- c("LAI", "Hmax", "fCover", "LAD")
  for (v in vars_disp) {
    sp_row <- sp[Variable == v]
    if (v == "LAD") {
      # report both JT scalars
      jt_fpc1 <- jt[Variable == "LAD_FPC1"]
      jt_l2   <- jt[Variable == "LAD_L2"]
      cat(sprintf("%-10s  %+10.2f     %-12s  FPC1: %-7s  %s\n",
                    v, sp_row$rho_S,
                    ifelse(sp_row$p_value < 1e-3, "<0.001",
                            sprintf("%.3f", sp_row$p_value)),
                    sprintf("%.3f", jt_fpc1$JT_pvalue),
                    "NULL (H2 OK)"))
      cat(sprintf("%-10s  %-13s  %-12s  L2  : %-7s  %s\n",
                    "  (alt scalar)", "", "",
                    sprintf("%.3f", jt_l2$JT_pvalue),
                    "NULL concordant"))
    } else {
      jt_row <- jt[Variable == v]
      cat(sprintf("%-10s  %+10.2f     %-12s  %-13s  %s\n",
                    v, sp_row$rho_S,
                    ifelse(sp_row$p_value < 1e-3, "<0.001",
                            sprintf("%.3f", sp_row$p_value)),
                    ifelse(jt_row$JT_pvalue < 1e-3, "<0.001",
                            sprintf("%.3f", jt_row$JT_pvalue)),
                    "CONFIRMED"))
    }
  }
  cat("\nFiles : outputs/lovb/figures/fig_heatmap_4col_simulation.png\n")
  cat("        outputs/lovb/figures/fig_HOBO_validation_spearman.png\n")
  cat("        outputs/lovb/figures/fig_HOBO_jonckheere_bins.png\n")
  cat("        outputs/lovb/tables/tab_{consensus_4methods,HOBO_spearman,HOBO_jonckheere}.csv\n")
}

# Step 1 entry -----------------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  step <- if (length(args) >= 1L) args[[1L]] else "1"

  if (step == "1") {
    build_heatmap_4col()
  } else if (step == "2A") {
    run_spearman_validation()
  } else if (step == "2B") {
    run_jonckheere()
  } else if (step == "3") {
    print_synthesis()
  } else {
    stop("Unknown step : ", step)
  }
}
