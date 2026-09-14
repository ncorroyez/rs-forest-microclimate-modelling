# ==============================================================================
# H1 — 2-level attribution heatmap + non-parametric HOBO validation
#
# Level 1 (archetype) : exact attribution at 4 K-means centroid plots
#                        Shapley_arch | LOVB_arch | LVA_arch
# Level 2 (cLHS)      : mean attribution over 400 cLHS plots
#                        Shapley_cLHS | LOVB_cLHS | LVA_cLHS
# Plus Consensus column (mean rank over the 6 methods).
#
# Non-parametric observational validation at 53 HOBO sensors :
#   Spearman + bootstrap CI95  (refresh from cache, LOESS smoother)
#   Jonckheere-Terpstra + bootstrap CI95 on z-statistic  (NEW)
#
# All cluster sizes equal in cLHS sample (n_C{1..4} = 100), so pooled
# weighted = pooled unweighted -> use simple mean.
#
# No mgcv. No dplyr / purrr / tidyr.
# Outputs under outputs/lovb/.
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
})

# ---- Step 1 : archetype-level attribution ------------------------------------
build_attribution_archetypes <- function() {
  cli_h1("Step 1 — archetype-level attribution (pooled across 4 K-means centroids)")
  vars <- c("LAI", "Hmax", "fCover", "LAD")

  # 1.1 Shapley : use pre-computed pooled_4archetypes row in CSV
  sh <- fread(here::here("outputs/h1/shapley_archetypes_values.csv"))
  sh_pool <- sh[archetype == "pooled_4archetypes",
                 .(variable, phi = phi, phi_abs = abs(phi),
                    ci95_lo, ci95_hi)]
  setnames(sh_pool, "phi", "Shapley_arch_value")
  sh_pool[, Shapley_arch_rank := frank(-phi_abs, ties.method = "min")]

  # 1.2 LOVB archetype : mean of |Delta_v| over 4 archetypes, Metric == mean
  lovb_a <- fread(here::here("outputs/lovb/tables/tab_LOVB_archetypes.csv"))
  lovb_pool <- lovb_a[Metric == "mean",
                       .(LOVB_arch_value = mean(abs(Delta_v))),
                       by = .(variable = Variable)]
  lovb_pool[, LOVB_arch_rank := frank(-LOVB_arch_value, ties.method = "min")]

  # 1.3 LVA archetype : same logic from DT_contrib_archetypes LVA_v_mean
  DT_a <- readRDS(here::here("outputs/lovb/data/DT_contrib_archetypes.rds"))
  lva_rows <- data.table(
    variable        = vars,
    LVA_arch_value  = c(mean(abs(DT_a$LVA_LAI_mean)),
                         mean(abs(DT_a$LVA_Hmax_mean)),
                         mean(abs(DT_a$LVA_fCover_mean)),
                         mean(abs(DT_a$LVA_LAD_mean)))
  )
  lva_rows[, LVA_arch_rank := frank(-LVA_arch_value, ties.method = "min")]

  # Merge
  tab <- merge(merge(sh_pool[, .(variable, Shapley_arch_value, Shapley_arch_rank)],
                       lovb_pool, by = "variable"),
                lva_rows, by = "variable")
  tab[, variable := factor(variable, levels = vars)]
  setorder(tab, variable)

  cli_alert("Archetype-pooled attribution :")
  print(tab)
  out_path <- here::here("outputs/lovb/tables/tab_attribution_archetypes_pooled.csv")
  fwrite(tab, out_path)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(tab)
}

# ---- Step 2 : cLHS-level attribution (cache reuse) ---------------------------
build_attribution_cLHS <- function() {
  cli_h1("Step 2 — cLHS-level attribution (mean over 400 plots)")
  vars <- c("LAI", "Hmax", "fCover", "LAD")

  # 2.1 Shapley cLHS : 02_shapley_attribution.csv
  sh <- fread(here::here("outputs/figures/02_shapley_attribution.csv"))
  sh <- sh[, .(variable, Shapley_cLHS_value = shapley,
                phi_abs = abs(shapley))]
  sh[, Shapley_cLHS_rank := frank(-phi_abs, ties.method = "min")]

  # 2.2 LOVB cLHS : tab_LOVB_global.csv, Metric == mean (RMSE_global)
  lovb_g <- fread(here::here("outputs/lovb/tables/tab_LOVB_global.csv"))
  lovb_pool <- lovb_g[Metric == "mean",
                       .(variable = Variable,
                          LOVB_cLHS_value = RMSE_global)]
  lovb_pool[, LOVB_cLHS_rank := frank(-LOVB_cLHS_value, ties.method = "min")]

  # 2.3 LVA cLHS : DT_contrib_cLHS$LVA_v_mean, mean of |LVA| over 400 plots
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  lva_rows <- data.table(
    variable = vars,
    LVA_cLHS_value = c(mean(abs(DT_c$LVA_LAI_mean),    na.rm = TRUE),
                        mean(abs(DT_c$LVA_Hmax_mean),   na.rm = TRUE),
                        mean(abs(DT_c$LVA_fCover_mean), na.rm = TRUE),
                        mean(abs(DT_c$LVA_LAD_mean),    na.rm = TRUE))
  )
  lva_rows[, LVA_cLHS_rank := frank(-LVA_cLHS_value, ties.method = "min")]

  tab <- merge(merge(sh[, .(variable, Shapley_cLHS_value, Shapley_cLHS_rank)],
                       lovb_pool, by = "variable"),
                lva_rows, by = "variable")
  tab[, variable := factor(variable, levels = vars)]
  setorder(tab, variable)

  cli_alert("cLHS-level attribution :")
  print(tab)
  out_path <- here::here("outputs/lovb/tables/tab_attribution_cLHS_mean.csv")
  fwrite(tab, out_path)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(tab)
}

# ---- Step 3 : 7-column 2-level heatmap ---------------------------------------
build_heatmap_2levels <- function() {
  cli_h1("Step 3 — 7-column 2-level heatmap")

  arch <- fread(here::here("outputs/lovb/tables/tab_attribution_archetypes_pooled.csv"))
  clhs <- fread(here::here("outputs/lovb/tables/tab_attribution_cLHS_mean.csv"))
  vars <- c("LAI", "fCover", "Hmax", "LAD")
  arch[, variable := factor(variable, levels = vars)]
  clhs[, variable := factor(variable, levels = vars)]

  # Build long table with 6 methods + consensus
  methods <- c("Shapley_arch", "LOVB_arch", "LVA_arch",
                "Shapley_cLHS", "LOVB_cLHS", "LVA_cLHS")
  labels  <- c("Shapley", "LOVB (mean)", "LVA (mean)",
                "Shapley", "LOVB (mean)", "LVA (mean)")
  groups  <- c("Level 1 — 4 K-means archetypes", "Level 1 — 4 K-means archetypes",
                "Level 1 — 4 K-means archetypes",
                "Level 2 — 400 cLHS plots", "Level 2 — 400 cLHS plots",
                "Level 2 — 400 cLHS plots")

  long_list <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    m <- methods[i]
    src <- if (grepl("_arch$", m)) arch else clhs
    rank_col <- paste0(m, "_rank")
    long_list[[i]] <- src[, .(variable, method = labels[i],
                                group = groups[i], col_idx = i,
                                Rank = get(rank_col))]
  }
  long <- rbindlist(long_list)

  # Consensus = mean rank over the 6 methods
  cons <- long[, .(Rank = mean(Rank), col_idx = 7L), by = variable]
  cons[, method := "Consensus"]
  cons[, group  := "Consensus"]
  long_full <- rbind(long[, .(variable, method, group, col_idx, Rank)],
                      cons[, .(variable, method, group, col_idx, Rank)])

  # X positions (numeric for geom_tile, then label via scale_x_continuous)
  long_full[, fill_rank := as.character(round(Rank))]

  # Vertical separator positions : between col 3 / 4 and 6 / 7
  sep_xs <- c(3.5, 6.5)

  pal <- c("1" = "#440154", "2" = "#3B528B",
            "3" = "#5DC863", "4" = "#FDE725")

  x_breaks <- 1:7
  x_labels <- c("Shapley", "LOVB", "LVA",
                 "Shapley", "LOVB", "LVA", "Consensus")
  top_labels_df <- data.table(
    x = c(2, 5, 7),
    label = c("Level 1 — archetype basis (4 K-means centroids)",
                "Level 2 — cLHS population (400 plots)",
                "Consensus")
  )

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(method == "Consensus",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2, "white", "grey15"),
               fontface = "bold", size = 5.5) +
    geom_vline(xintercept = sep_xs, colour = "white", linewidth = 3) +
    geom_text(data = top_labels_df,
               aes(x = x, y = 4.85, label = label),
               inherit.aes = FALSE, fontface = "bold", size = 3.3,
               colour = "grey20") +
    scale_fill_manual(values = pal, name = "Rank") +
    scale_x_continuous(breaks = x_breaks, labels = x_labels,
                        expand = c(0, 0)) +
    scale_y_discrete(limits = rev, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.5, 5.1), clip = "off") +
    labs(
      title    = "Attribution ranking at 2 aggregation levels (7-column triangulation)",
      subtitle = "Dark = rank 1 (strongest attributed effect)   |   light = rank 4 (weakest)",
      x = NULL, y = NULL,
      caption = paste0(
        "Level 1 — exact attribution at 4 K-means archetype centroids ",
        "(pooled by equal cluster size, n=100 each).\n",
        "Level 2 — mean attribution over 400 cLHS plots.\n",
        "LAI consistently rank 1-2 ; LAD consistently rank 3-4 across all 6 methods and both aggregation scales."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 14),
           plot.subtitle = element_text(colour = "grey25", size = 11),
           plot.caption  = element_text(size = 8, colour = "grey35", hjust = 0,
                                         lineheight = 1.3),
           axis.text.x   = element_text(face = "bold", size = 10),
           axis.text.y   = element_text(face = "bold", size = 13),
           plot.margin   = margin(15, 10, 10, 10),
           panel.grid    = element_blank())

  fig_path <- here::here("outputs/lovb/figures/fig_heatmap_2levels.png")
  ggsave(fig_path, p, width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Helper : load all HOBO scalars (slope_obs + dTmax_obs + traits) ---------
load_hobo_scalars_full <- function() {
  DT_h <- readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds"))   # slope_obs, LAI, Hmax, fCover, FPC1, LAD_L2
  cg   <- readRDS(here::here("outputs/lovb/data/DT_cross_gam_targets.rds"))
  hobo_obs <- as.data.table(cg$HOBO)[, .(id_plot = id, dTmax_obs = dTmax_mean)]
  merge(DT_h, hobo_obs, by = "id_plot")
}

# ---- Spearman + bootstrap helper (generic) -----------------------------------
spearman_table <- function(x_obs, ymat, var_names, B = 1000L, seed = 42L) {
  ok_x <- !is.na(x_obs)
  rows <- vector("list", length(var_names))
  for (i in seq_along(var_names)) {
    y  <- ymat[[i]]
    ok <- ok_x & !is.na(y)
    rho_hat <- suppressWarnings(cor(x_obs[ok], y[ok], method = "spearman"))
    p_val   <- suppressWarnings(cor.test(x_obs[ok], y[ok], method = "spearman",
                                            exact = FALSE)$p.value)
    set.seed(seed + i)
    n_ok <- sum(ok); xs <- x_obs[ok]; ys <- y[ok]
    rb <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n_ok, n_ok, replace = TRUE)
      rb[b] <- suppressWarnings(cor(xs[idx], ys[idx], method = "spearman"))
    }
    ci <- quantile(rb, c(0.025, 0.975), na.rm = TRUE, type = 7)
    rows[[i]] <- data.table(Variable = var_names[i],
                              rho_S = rho_hat, p_value = p_val,
                              ci_lo = unname(ci[1]), ci_hi = unname(ci[2]))
  }
  rbindlist(rows)
}

# ---- Step 4A : Spearman with LOESS, both slope_obs and dTmax_obs -------------
run_spearman_dual <- function() {
  cli_h1("Step 4A — Spearman dual (slope_obs + dTmax_obs) with LOESS")

  DT_h <- load_hobo_scalars_full()
  DT_sim <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  DT <- merge(DT_h[, .(id_plot, slope_obs, dTmax_obs)],
                DT_sim[, .(id_plot,
                            Delta_LAI_mean, Delta_Hmax_mean,
                            Delta_fCover_mean, Delta_LAD_mean)],
                by = "id_plot")
  cli_alert("n = {nrow(DT)} HOBO sensors  |  slope_obs and dTmax_obs joined")
  vars  <- c("LAI", "Hmax", "fCover", "LAD")
  ycols <- c("Delta_LAI_mean", "Delta_Hmax_mean",
              "Delta_fCover_mean", "Delta_LAD_mean")

  ymat <- DT[, ..ycols]
  tab_slope <- spearman_table(DT$slope_obs,  ymat, vars)
  tab_dTmax <- spearman_table(DT$dTmax_obs,  ymat, vars)
  tab_slope[, x_metric := "slope_obs"]
  tab_dTmax[, x_metric := "dTmax_obs"]

  fwrite(tab_slope, here::here("outputs/lovb/tables/tab_HOBO_spearman.csv"))
  fwrite(tab_dTmax, here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))
  cli_alert("Spearman (slope_obs) :"); print(tab_slope)
  cli_alert("Spearman (dTmax_obs) :"); print(tab_dTmax)

  # 2 separate figures (independent x axes per metric, cleaner)
  make_fig <- function(x_vec, tab_ann, x_label, x_axis_title, suffix) {
    long <- rbindlist(lapply(seq_along(vars), function(i)
      data.table(Variable = vars[i],
                  x_obs   = x_vec,
                  Delta_v = ymat[[i]])))
    long[, Variable := factor(Variable, levels = vars)]
    ann <- copy(tab_ann)
    ann[, label := sprintf("rho = %+.2f  p = %s\n[%+.2f ; %+.2f]",
                             rho_S,
                             ifelse(p_value < 1e-3, "<0.001",
                                     sprintf("%.3f", p_value)),
                             ci_lo, ci_hi)]
    ann[, Variable := factor(Variable, levels = vars)]

    p <- ggplot(long, aes(x = x_obs, y = Delta_v)) +
      geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
      geom_point(colour = "#2C5F2D", alpha = 0.65, size = 1.8) +
      geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                    colour = "#FFB400", fill = "#FFB400",
                    linewidth = 0.8, alpha = 0.18, span = 0.9) +
      geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
                 hjust = -0.05, vjust = 1.2, size = 3.1, fontface = "bold",
                 inherit.aes = FALSE, colour = "grey15") +
      facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
      labs(
        title    = paste0("HOBO validation (", x_label,
                            ") : observed vs simulated LOVB delta"),
        subtitle = "n = 53 sensors  |  Spearman rho + bootstrap CI95 (1000 reps)  |  LOESS span = 0.9",
        x = x_axis_title,
        y = "Simulated LOVB delta (K, mean over 122 d)",
        caption = "Non-parametric rank correlation. LOESS shown for visual reference only."
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey92", colour = NA),
             strip.text       = element_text(face = "bold", size = 11),
             plot.title       = element_text(face = "bold"),
             plot.subtitle    = element_text(colour = "grey25"),
             plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
             panel.grid.minor = element_blank())
    fig_path <- here::here(paste0("outputs/lovb/figures/hobo/fig_HOBO_LOVB_spearman_",
                                     suffix, ".png"))
    ggsave(fig_path, p, width = 11.4, height = 4.3, dpi = 300)
    cli_alert_success("Saved {.path {fig_path}}")
  }
  make_fig(DT$slope_obs, tab_slope, "slope_obs",
            "Observed buffering slope (Tmax_obs ~ Tmax_macro)", "slope")
  make_fig(DT$dTmax_obs, tab_dTmax, "dTmax_obs",
            "Observed buffering mean  Tmax_macro − Tmax_obs  (K)", "dTmax")
  invisible(list(slope = tab_slope, dTmax = tab_dTmax))
}

# ---- (legacy) Step 4A : Spearman figure (LOESS smoother) ---------------------
refresh_spearman_figure <- function() {
  cli_h1("Step 4A (legacy) — Spearman figure refresh (LOESS smoother)")

  DT_h <- readRDS(here::here("outputs/lovb/data/DT_HOBO_scalars.rds"))
  DT_sim <- readRDS(here::here("outputs/lovb/data/DT_contrib_HOBO.rds"))
  DT <- merge(DT_h[, .(id_plot, slope_obs)],
                DT_sim[, .(id_plot,
                            Delta_LAI_mean, Delta_Hmax_mean,
                            Delta_fCover_mean, Delta_LAD_mean)],
                by = "id_plot")
  vars  <- c("LAI", "Hmax", "fCover", "LAD")
  ycols <- c("Delta_LAI_mean", "Delta_Hmax_mean",
              "Delta_fCover_mean", "Delta_LAD_mean")

  long <- rbindlist(lapply(seq_along(vars), function(i)
    data.table(Variable = vars[i],
                slope_obs = DT$slope_obs,
                Delta_v   = DT[[ycols[i]]])))
  long[, Variable := factor(Variable, levels = vars)]

  tab <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman.csv"))
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
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                  colour = "#FFB400", fill = "#FFB400",
                  linewidth = 0.8, alpha = 0.18, span = 0.9) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.2, size = 3.1, fontface = "bold",
               inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
    labs(
      title    = "HOBO validation : observed slope vs simulated LOVB delta",
      subtitle = "n = 53 sensors  |  Spearman rho + bootstrap CI95 (1000 reps)  |  LOESS smoother",
      x = "Observed buffering slope (Tmax_obs ~ Tmax_macro)",
      y = "Simulated LOVB delta (K, mean over 122 d)",
      caption = "Non-parametric rank correlation. LOESS shown for visual reference only (span = 0.9)."
    ) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey92", colour = NA),
           strip.text       = element_text(face = "bold", size = 11),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
           panel.grid.minor = element_blank())
  fig_path <- here::here("outputs/lovb/figures/hobo/fig_HOBO_LOVB_spearman_legacy.png")
  ggsave(fig_path, p, width = 11.4, height = 4.3, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Jonckheere-Terpstra with bootstrap CI95 ---------------------------------
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
      for (a in x_i) J <- J + sum(x_j > a) + 0.5 * sum(x_j == a)
    }
  }
  n   <- length(x); n_g <- tabulate(group, nbins = k)
  mu  <- (n^2 - sum(n_g^2)) / 4
  sigma2 <- (n^2 * (2*n + 3) - sum(n_g^2 * (2*n_g + 3))) / 72
  z <- (J - mu) / sqrt(sigma2)
  p <- switch(alternative,
                "two.sided"  = 2 * (1 - pnorm(abs(z))),
                "increasing" = 1 - pnorm(z),
                "decreasing" = pnorm(z))
  list(statistic = J, z = z, p.value = p, alternative = alternative)
}

jt_boot_ci <- function(x, group, alternative, B = 1000L, seed = 42L) {
  set.seed(seed)
  n <- length(x)
  zb <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    zb[b] <- tryCatch(jt_test(x[idx], group[idx], alternative)$z,
                       error = function(e) NA_real_)
  }
  quantile(zb, c(0.025, 0.975), na.rm = TRUE, type = 7)
}

run_jonckheere_dual <- function(B = 1000L) {
  cli_h1("Step 4B — JT dual (slope_obs + dTmax_obs) with bootstrap CI95")
  DT <- load_hobo_scalars_full()
  vars <- c("LAI", "Hmax", "fCover", "LAD_FPC1", "LAD_L2")
  trait_cols <- c(LAI = "LAI", Hmax = "Hmax", fCover = "fCover",
                   LAD_FPC1 = "FPC1", LAD_L2 = "LAD_L2")
  # Physical direction depends on metric :
  #   slope_obs : trait ↑  ->  slope ↓     (decreasing)
  #   dTmax_obs : trait ↑  ->  buffering ↑ (increasing)
  alt_slope <- c(LAI = "decreasing", Hmax = "decreasing",
                  fCover = "decreasing",
                  LAD_FPC1 = "two.sided", LAD_L2 = "two.sided")
  alt_dTmax <- c(LAI = "increasing", Hmax = "increasing",
                  fCover = "increasing",
                  LAD_FPC1 = "two.sided", LAD_L2 = "two.sided")

  run_one_metric <- function(y_col, suffix, y_label, alt_lookup) {
    rows <- list(); bins_long <- list()
    for (v in vars) {
      tc <- trait_cols[v]
      trait <- DT[[tc]]; y <- DT[[y_col]]
      keep <- !is.na(trait) & !is.na(y)
      trait <- trait[keep]; y <- y[keep]
      n <- length(trait)
      rnk <- rank(trait, ties.method = "first")
      bin <- cut(rnk, breaks = c(0, n/3, 2*n/3, n),
                  include.lowest = TRUE, labels = c("low", "mid", "high"))
      bin <- factor(bin, levels = c("low", "mid", "high"))

      res <- jt_test(y, bin, alternative = alt_lookup[v])
      ci  <- jt_boot_ci(y, bin, alt_lookup[v], B = B,
                          seed = 42L + which(vars == v))

      med <- tapply(y, bin, median); iqr <- tapply(y, bin, IQR)
      trend <- if (res$p.value < 0.05) {
        if (res$z > 0) "Monotonic increase" else "Monotonic decrease"
      } else "No monotonic trend"

      rows[[length(rows) + 1L]] <- data.table(
        Variable       = v,
        Metric         = suffix,
        Alternative    = alt_lookup[v],
        med_low = round(med["low"], 3), med_mid = round(med["mid"], 3),
        med_high = round(med["high"], 3),
        iqr_low = round(iqr["low"], 3), iqr_mid = round(iqr["mid"], 3),
        iqr_high = round(iqr["high"], 3),
        JT_statistic = res$statistic,
        JT_z         = round(res$z, 3),
        JT_pvalue    = signif(res$p.value, 4),
        JT_z_ci_lo   = round(ci[1], 3),
        JT_z_ci_hi   = round(ci[2], 3),
        Interpretation = trend
      )
      bins_long[[v]] <- data.table(Variable = v, bin = bin, y = y,
                                     p_value = res$p.value, z = res$z,
                                     ci_lo = ci[1], ci_hi = ci[2])
    }
    tab <- rbindlist(rows)
    tab_path <- here::here(paste0("outputs/lovb/tables/tab_HOBO_jonckheere_",
                                     suffix, ".csv"))
    fwrite(tab, tab_path)
    cli_alert("JT ({suffix}) :")
    print(tab[, .(Variable, JT_z, JT_z_ci_lo, JT_z_ci_hi,
                    JT_pvalue, Interpretation)])
    cli_alert_success("Saved {.path {tab_path}}")

    # Figure
    long <- rbindlist(bins_long)
    long[, Variable := factor(Variable, levels = vars)]
    ann <- tab[, .(Variable, label = sprintf("JT p = %s\nz = %+.2f [%+.2f ; %+.2f]",
                                                ifelse(JT_pvalue < 1e-3, "<0.001",
                                                        sprintf("%.3f", JT_pvalue)),
                                                JT_z, JT_z_ci_lo, JT_z_ci_hi))]
    ann[, Variable := factor(Variable, levels = vars)]
    p <- ggplot(long, aes(x = bin, y = y)) +
      geom_boxplot(fill = "#2C5F2D", alpha = 0.55, colour = "grey20",
                    outlier.size = 1.3) +
      geom_jitter(width = 0.12, alpha = 0.55, size = 1, colour = "grey30") +
      geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
                 hjust = -0.05, vjust = 1.3, size = 2.9, fontface = "bold",
                 inherit.aes = FALSE) +
      facet_wrap(~ Variable, nrow = 1L, scales = "free_y") +
      labs(
        title    = paste0("Jonckheere-Terpstra trend test (", suffix, ")"),
        subtitle = paste0("Tertile binning by trait  |  n = 53 HOBO  |  ",
                            "JT z bootstrap CI95 (", B, " reps)"),
        x = NULL, y = y_label,
        caption = paste0("Alternative direction = ",
                          if (suffix == "slope") "decreasing" else "increasing",
                          " for LAI/Hmax/fCover  (physical : trait ↑  ->  buffering ↑),  ",
                          "two.sided for LAD (FPC1 and L2).")
      ) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey92", colour = NA),
             strip.text       = element_text(face = "bold"),
             plot.title       = element_text(face = "bold"),
             plot.subtitle    = element_text(colour = "grey25"),
             plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0),
             panel.grid.minor = element_blank())
    fig_path <- here::here(paste0("outputs/lovb/figures/hobo/fig_HOBO_LOVB_jonckheere_",
                                     suffix, ".png"))
    ggsave(fig_path, p, width = 13, height = 4.3, dpi = 300)
    cli_alert_success("Saved {.path {fig_path}}")
    tab
  }

  tab_slope <- run_one_metric("slope_obs", "slope",
                                "Observed slope (Tmax_obs ~ Tmax_macro)",
                                alt_slope)
  tab_dTmax <- run_one_metric("dTmax_obs", "dTmax",
                                "Observed buffering mean Tmax_macro − Tmax_obs (K)",
                                alt_dTmax)
  invisible(list(slope = tab_slope, dTmax = tab_dTmax))
}

# ---- Step 5 : synthesis console + interpretation table ----------------------
classify_interpretation <- function(rho, jt_p, jt_z) {
  bulk_dir <- (jt_z < 0 && rho > 0) ||   # slope metric : rho>0 + JT decreasing
                (jt_z > 0 && rho < 0)       # dTmax metric : rho<0 + JT increasing
  rho_abs <- abs(rho)
  if (rho_abs > 0.3 && jt_p < 0.05 && bulk_dir) "BULK CONFIRMED"
  else if (rho_abs > 0.3 && jt_p < 0.05 && !bulk_dir) "INCONSISTENT"
  else if ((rho_abs >= 0.2 && rho_abs <= 0.3) ||
            (jt_p >= 0.05 && jt_p <= 0.10))  "BULK MARGINAL"
  else if (rho_abs < 0.2 && jt_p > 0.10)     "NULL (H2 OK)"
  else                                         "AMBIGUOUS"
}

build_validation_table <- function() {
  sp_s <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman.csv"))
  sp_d <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))
  jt_s <- fread(here::here("outputs/lovb/tables/tab_HOBO_jonckheere_slope.csv"))
  jt_d <- fread(here::here("outputs/lovb/tables/tab_HOBO_jonckheere_dTmax.csv"))

  # JT for LAD : prefer FPC1 (consistent with rest of pipeline)
  jt_s_main <- jt_s[Variable %in% c("LAI", "Hmax", "fCover", "LAD_FPC1")]
  jt_d_main <- jt_d[Variable %in% c("LAI", "Hmax", "fCover", "LAD_FPC1")]
  jt_s_main[, Variable := gsub("_FPC1$", "", Variable)]
  jt_d_main[, Variable := gsub("_FPC1$", "", Variable)]

  build_one <- function(sp, jt, metric) {
    DT <- merge(sp[, .(Variable, rho_S, p_value, ci_lo, ci_hi)],
                  jt[, .(Variable, JT_z, JT_pvalue,
                          JT_z_ci_lo, JT_z_ci_hi)],
                  by = "Variable")
    DT[, x_metric := metric]
    DT[, Interpretation := mapply(classify_interpretation,
                                    rho_S, JT_pvalue, JT_z)]
    setcolorder(DT, c("Variable", "x_metric"))
    DT
  }
  tab <- rbind(build_one(sp_s, jt_s_main, "slope_obs"),
                 build_one(sp_d, jt_d_main, "dTmax_obs"))
  tab[, Variable := factor(Variable, levels = c("LAI", "Hmax", "fCover", "LAD"))]
  setorder(tab, x_metric, Variable)

  out_path <- here::here("outputs/lovb/tables/tab_HOBO_validation_nonparam.csv")
  fwrite(tab, out_path)
  cli_alert_success("Saved {.path {out_path}}")
  tab
}

print_synthesis <- function() {
  cli_h1("Step 5 — Synthesis")

  # 2-level attribution table
  arch <- fread(here::here("outputs/lovb/tables/tab_attribution_archetypes_pooled.csv"))
  clhs <- fread(here::here("outputs/lovb/tables/tab_attribution_cLHS_mean.csv"))
  vars <- c("LAI", "fCover", "Hmax", "LAD")
  arch[, variable := factor(variable, levels = vars)]
  clhs[, variable := factor(variable, levels = vars)]
  setorder(arch, variable); setorder(clhs, variable)

  cat("\n=====================================================\n")
  cat(" ATTRIBUTION RANKING — 2 LEVELS, 6 METHODS + CONSENSUS\n")
  cat("=====================================================\n")
  cat(sprintf("%-8s  | %-19s | %-19s | %s\n",
                "Variable",
                "  ARCHETYPES (n=4) ", "    cLHS (n=400)   ",
                "Consensus"))
  cat(sprintf("%-8s  | %-5s %-5s %-5s | %-5s %-5s %-5s |\n",
                "        ",
                "Shap","LOVB","LVA", "Shap","LOVB","LVA"))
  cons_vals <- numeric(length(vars))
  for (k in seq_along(vars)) {
    v <- vars[k]
    ra <- arch[variable == v]; rc <- clhs[variable == v]
    cons_vals[k] <- mean(c(ra$Shapley_arch_rank, ra$LOVB_arch_rank, ra$LVA_arch_rank,
                            rc$Shapley_cLHS_rank, rc$LOVB_cLHS_rank, rc$LVA_cLHS_rank))
  }
  cons_ranks <- frank(cons_vals, ties.method = "min")
  for (k in seq_along(vars)) {
    v <- vars[k]
    ra <- arch[variable == v]; rc <- clhs[variable == v]
    cat(sprintf("%-8s  | %-5d %-5d %-5d | %-5d %-5d %-5d | %.1f  (#%d)\n",
                  v,
                  ra$Shapley_arch_rank, ra$LOVB_arch_rank, ra$LVA_arch_rank,
                  rc$Shapley_cLHS_rank, rc$LOVB_cLHS_rank, rc$LVA_cLHS_rank,
                  cons_vals[k], cons_ranks[k]))
  }

  # Observational validation
  tab <- build_validation_table()

  cat("\n=====================================================\n")
  cat(" OBSERVATIONAL VALIDATION (n=53 HOBO, non-parametric, dual y-metric)\n")
  cat("=====================================================\n")
  cat(sprintf("%-8s  | %-9s | %-20s | %-9s %-7s | %s\n",
                "Variable", "y-metric", "Spearman rho [CI95]",
                "JT p   ", "JT z  ", "Interpretation"))
  for (m in c("slope_obs", "dTmax_obs")) {
    for (v in vars) {
      r <- tab[Variable == v & x_metric == m]
      cat(sprintf("%-8s  | %-9s | %+5.2f [%+5.2f ; %+5.2f] | %-9s %+5.2f | %s\n",
                    v, m, r$rho_S, r$ci_lo, r$ci_hi,
                    ifelse(r$JT_pvalue < 1e-3, "<0.001",
                            sprintf("%.3f", r$JT_pvalue)),
                    r$JT_z, r$Interpretation))
    }
    cat("\n")
  }
}

# ---- Step 4C : daily-level Spearman with block bootstrap by sensor -----------
# Builds the daily observed dTmax panel (sensor x day) and the daily simulated
# Delta_v_LOVB panel from DT_daily_HOBO. Pooled Spearman across all daily pairs,
# block bootstrap resampling SENSORS (not days) to handle autocorrelation.

DAILY_OBS_CACHE <- here::here("outputs/lovb/data/DT_HOBO_daily_obs.rds")

build_daily_obs_panel <- function(force = FALSE) {
  if (!force && file.exists(DAILY_OBS_CACHE)) {
    cli_alert("Reusing cached daily observed panel : {.path {DAILY_OBS_CACHE}}")
    return(readRDS(DAILY_OBS_CACHE))
  }
  cli_h2("Build daily observed dTmax panel (sensor x day)")
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

  df_macro       <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
  df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq,
                                     df_macro, CFG$ids_to_remove)
  m <- as.data.table(df_hobo_daily)
  if (!("Tmax_macro" %in% names(m)))
    m <- merge(m, as.data.table(df_macro)[, .(date, Tmax_macro)], by = "date")
  m <- m[, .(id_plot, date, Tmax_obs, Tmax_macro,
               dTmax_obs = Tmax_macro - Tmax_obs)]
  m <- m[!is.na(dTmax_obs)]
  saveRDS(m, DAILY_OBS_CACHE)
  cli_alert_success("Saved {.path {DAILY_OBS_CACHE}}  ({nrow(m)} sensor-day rows)")
  m
}

# Pooled daily Spearman with block bootstrap by sensor
run_daily_spearman_block <- function(B = 1000L, seed = 42L) {
  cli_h1("Step 4C — daily Spearman + block bootstrap by sensor")
  obs <- build_daily_obs_panel()
  sim <- as.data.table(readRDS(here::here("outputs/lovb/data/DT_daily_HOBO.rds")))

  # Build daily Delta_v_LOVB sim per variable :
  # Delta_v_LOVB(s,d) = DeltaTmax_sim(s,d | 1111) - DeltaTmax_sim(s,d | LOVB_v)
  bit_map <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110",
                REF = "1111")
  vars <- c("LAI", "Hmax", "fCover", "LAD")

  sim_ref <- sim[bit_code == bit_map[["REF"]],
                   .(id_plot, date, dTmax_REF = Delta_Tmax)]
  sim_w   <- merge(obs[, .(id_plot, date, dTmax_obs)],
                    sim_ref, by = c("id_plot", "date"))
  for (v in vars) {
    sim_v <- sim[bit_code == bit_map[[v]],
                  .(id_plot, date, dTmax_LOVB = Delta_Tmax)]
    setnames(sim_v, "dTmax_LOVB", paste0("dTmax_LOVB_", v))
    sim_w <- merge(sim_w, sim_v, by = c("id_plot", "date"))
  }
  for (v in vars) {
    col_lovb <- paste0("dTmax_LOVB_", v)
    sim_w[, paste0("Delta_", v, "_LOVB") := dTmax_REF - get(col_lovb)]
  }
  cli_alert("Daily panel : {nrow(sim_w)} pairs across {length(unique(sim_w$id_plot))} sensors")

  # Pooled Spearman estimate (all daily pairs)
  spear_pooled <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

  sensors <- unique(sim_w$id_plot)
  set.seed(seed)

  results <- vector("list", length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    delta_col <- paste0("Delta_", v, "_LOVB")
    rho_hat <- spear_pooled(sim_w$dTmax_obs, sim_w[[delta_col]])

    # Block bootstrap : resample sensors with replacement
    rho_b <- numeric(B)
    for (b in seq_len(B)) {
      idx_sens <- sample(sensors, length(sensors), replace = TRUE)
      sim_b <- sim_w[id_plot %in% idx_sens]
      rho_b[b] <- spear_pooled(sim_b$dTmax_obs, sim_b[[delta_col]])
    }
    ci <- quantile(rho_b, c(0.025, 0.975), na.rm = TRUE, type = 7)
    results[[i]] <- data.table(
      Variable        = v,
      rho_daily       = rho_hat,
      ci_lo_block     = unname(ci[1]),
      ci_hi_block     = unname(ci[2]),
      n_pairs         = nrow(sim_w),
      n_sensors       = length(sensors),
      B_boot          = B
    )
  }
  tab <- rbindlist(results)
  cli_alert("Daily Spearman with block bootstrap CI95 :")
  print(tab)

  # Compare with sensor-level results
  sp_sensor <- fread(here::here("outputs/lovb/tables/tab_HOBO_spearman_dTmax.csv"))
  compare <- merge(tab,
                    sp_sensor[, .(Variable, rho_sensor = rho_S,
                                    sensor_ci_lo = ci_lo, sensor_ci_hi = ci_hi)],
                    by = "Variable")
  compare[, Variable := factor(Variable, levels = vars)]
  setorder(compare, Variable)
  cli_alert("Sensor-level (n=53) vs daily-level (block-boot) comparison :")
  print(compare[, .(Variable, rho_sensor, sensor_ci_lo, sensor_ci_hi,
                     rho_daily, ci_lo_block, ci_hi_block)])

  out_path <- here::here("outputs/lovb/tables/tab_HOBO_spearman_daily_block.csv")
  fwrite(compare, out_path)
  cli_alert_success("Saved {.path {out_path}}")
  invisible(compare)
}

# ---- script entry ------------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  step <- if (length(args) >= 1L) args[[1L]] else "1"

  if (step == "1") {
    build_attribution_archetypes()
  } else if (step == "2") {
    build_attribution_cLHS()
  } else if (step == "3") {
    build_heatmap_2levels()
  } else if (step == "4A") {
    run_spearman_dual()
  } else if (step == "4B") {
    run_jonckheere_dual()
  } else if (step == "4C") {
    run_daily_spearman_block()
  } else if (step == "5") {
    print_synthesis()
  } else {
    stop("Unknown step : ", step)
  }
}
