# ==============================================================================
# H1 — Extra Type B figures (LVA + Shapley) and reordered heatmap.
#
# Reuses sign conventions from R/lovb_03_contrib.R as validated :
#   LOVB Delta_v = Tmax_REF - Tmax_LOVB_v    (positive = v contributes)
#   LVA  Delta_v = Tmax_LVA_v - Tmax_NULL    (positive = v contributes)
#   Shapley phi_v = exact attribution on 2^4 = 16 coalitions
#
# Figure 1 (LOVB) already exists - DO NOT regenerate.
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
})

# Same constants as R/lovb_04_figures.R for visual consistency
.VAR_ORDER    <- c("LAI", "fCover", "Hmax", "LAD")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")

# ---- Figure 2 : LVA Type B ---------------------------------------------------
build_typeB_LVA <- function(metric = "mean") {
  cli_h1("Figure 2 — Type B LVA (mean over 122 d)")
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))

  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")
  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lva_col   <- paste0("LVA_", v, "_", metric)
    trait_col <- trait_map[v]
    data.table(
      Cluster   = DT_c$Cluster,
      variable  = v,
      trait     = DT_c[[trait_col]],
      LVA       = DT_c[[lva_col]]
    )
  }))
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  cli_alert("LVA Type B : {nrow(long) / length(.VAR_ORDER)} plots per panel")

  p <- ggplot(long, aes(x = trait, y = LVA)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                linewidth = 0.6) +
    geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.6) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE,
                  colour = "grey20", fill = "grey70", alpha = 0.3) +
    facet_wrap(~ variable, nrow = 1, scales = "free_x",
                labeller = labeller(variable = function(x)
                                       paste("Contribution of", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    labs(
      title    = sprintf("Contribution Delta_v vs trait (%s) — LVA",
                          toupper(metric)),
      subtitle = "x-axis = real value of v in the plot (FPC1 for LAD shape) | LOESS span=0.75 with 95%CI",
      x = "Trait value (real per plot)",
      y = sprintf("LVA_v   Tmax (%s)  [degree C]", metric)
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
           plot.title      = element_text(face = "bold"),
           legend.position = "bottom")

  fig_path <- here::here(sprintf("outputs/lovb/figures/fig_typeB_LVA_%s.png",
                                    metric))
  # Same dimensions as fig_LOVB_typeB_mean.png : 14 x 5.5 in @ 300 dpi
  ggsave(fig_path, p, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Figure 4 : reordered heatmap v3 -----------------------------------------
build_heatmap_v3 <- function() {
  cli_h1("Figure 4 — reordered heatmap (LOVB -> LVA -> Shapley)")

  arch <- fread(here::here("outputs/lovb/tables/tab_attribution_archetypes_pooled.csv"))
  clhs <- fread(here::here("outputs/lovb/tables/tab_attribution_cLHS_mean.csv"))
  vars <- c("LAI", "fCover", "Hmax", "LAD")
  arch[, variable := factor(variable, levels = vars)]
  clhs[, variable := factor(variable, levels = vars)]

  # New column order : LOVB | LVA | Shapley  (instead of Shapley | LOVB | LVA)
  methods <- c("LOVB_arch",   "LVA_arch",   "Shapley_arch",
                "LOVB_cLHS",   "LVA_cLHS",   "Shapley_cLHS")
  labels  <- c("LOVB (mean)", "LVA (mean)", "Shapley",
                "LOVB (mean)", "LVA (mean)", "Shapley")
  groups  <- c(rep("Level 1 — 4 K-means archetypes", 3L),
                rep("Level 2 — 400 cLHS plots", 3L))

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

  cons <- long[, .(Rank = mean(Rank), col_idx = 7L), by = variable]
  cons[, method := "Consensus"]; cons[, group := "Consensus"]
  long_full <- rbind(long[, .(variable, method, group, col_idx, Rank)],
                      cons[, .(variable, method, group, col_idx, Rank)])
  long_full[, fill_rank := as.character(round(Rank))]

  pal <- c("1" = "#440154", "2" = "#3B528B",
            "3" = "#5DC863", "4" = "#FDE725")

  x_breaks <- 1:7
  x_labels <- c("LOVB", "LVA", "Shapley",
                 "LOVB", "LVA", "Shapley", "Consensus")
  top_labels_df <- data.table(
    x     = c(2, 5, 7),
    label = c("Level 1 — archetype basis (4 K-means centroids)",
                "Level 2 — cLHS population (400 plots)",
                "Consensus")
  )
  sep_xs <- c(3.5, 6.5)

  p <- ggplot(long_full, aes(x = col_idx, y = variable, fill = fill_rank)) +
    geom_tile(colour = "white", linewidth = 1.0) +
    geom_text(aes(label = ifelse(method == "Consensus",
                                   sprintf("%.1f", Rank),
                                   sprintf("%g", Rank))),
               colour = ifelse(as.numeric(long_full$fill_rank) <= 2,
                                 "white", "grey15"),
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
      title    = "Attribution ranking at 2 aggregation levels (column order : LOVB -> LVA -> Shapley)",
      subtitle = "Dark = rank 1 (strongest attributed effect)  |  light = rank 4 (weakest)",
      x = NULL, y = NULL,
      caption = paste0(
        "Column order reflects narrative flow : ",
        "LOVB (remove one trait from full canopy) -> LVA (add one trait alone vs null) -> Shapley (exact average over 16 coalitions).\n",
        "Level 1 = exact attribution at 4 K-means archetypes (pooled by equal cluster size, n=100 each).  ",
        "Level 2 = mean attribution over 400 cLHS plots."
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

  fig_path <- here::here("outputs/lovb/figures/fig_heatmap_2levels_v3.png")
  ggsave(fig_path, p, width = 12, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Figure 5 : Type A LVA strip-violin --------------------------------------
# Tmax_NULL is constant across the 400 plots (same baseline canopy) so a normal
# REF-vs-LOVB scatter degenerates. We show distribution of Tmax_LVA_v as
# violin+jitter per variable, with horizontal reference line at Tmax_NULL.
build_typeA_LVA <- function(metric = "mean") {
  cli_h1("Figure 5 — Type A LVA (strip-violin, Tmax_NULL constant)")
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  null_val <- unique(DT_c[[paste0("Tmax_", metric, "_NULL")]])
  cli_alert("Tmax_NULL ({metric}) = {round(null_val, 4)}  (constant on {nrow(DT_c)} plots)")

  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lva_col <- paste0("Tmax_", metric, "_LVA_", v)
    data.table(Cluster = DT_c$Cluster, variable = v,
                Tmax_LVA = DT_c[[lva_col]])
  }))
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  p <- ggplot(long, aes(x = Cluster_lab, y = Tmax_LVA)) +
    geom_hline(yintercept = null_val, linetype = "dashed", colour = "red",
                linewidth = 0.6) +
    geom_violin(aes(fill = Cluster_lab), alpha = 0.35,
                  colour = "grey30", linewidth = 0.4) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.55,
                  size = 1.1) +
    facet_wrap(~ variable, nrow = 1, scales = "free_y",
                labeller = labeller(variable = function(x)
                                       paste0("Tmax_LVA_", x, " (v alone real)"))) +
    scale_fill_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                        name = NULL) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    labs(
      title    = sprintf("Type A LVA — distribution of Tmax_LVA_v per cluster (%s)",
                          toupper(metric)),
      subtitle = sprintf("Red dashed line = Tmax_NULL = %.3f K (constant across all 400 plots, identical baseline canopy)",
                          null_val),
      x = NULL,
      y = sprintf("Tmax LVA_v  (%s)  [degree C]", metric),
      caption = "Strip-violin (not scatter) : Tmax_NULL is identical on every plot because the baseline canopy (LAI=mu, Hmax=mu, fCover=1, LAD=uniform) is plot-invariant. Vertical spread = effect of activating only v on the plot."
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text       = element_text(face = "bold"),
           plot.title       = element_text(face = "bold"),
           plot.subtitle    = element_text(colour = "grey25"),
           plot.caption     = element_text(size = 8, colour = "grey35", hjust = 0,
                                            lineheight = 1.3),
           legend.position  = "bottom",
           axis.text.x      = element_text(angle = 25, hjust = 1))

  fig_path <- here::here(sprintf("outputs/lovb/figures/fig_typeA_LVA_%s.png",
                                    metric))
  ggsave(fig_path, p, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

# ---- Figure 6 : Type A Shapley consistency vs LOVB ---------------------------
# Scatter phi_v (per-plot Shapley, 16 coalitions) vs Delta_v_LOVB (per-plot,
# max-order approximation). Annotated with Spearman rho per panel.
build_typeA_Shapley_consistency <- function() {
  cli_h1("Figure 6 — Type A Shapley consistency (phi_v vs Delta_v_LOVB)")
  DT_phi <- readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
  DT_c   <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))

  # Reshape DT_c long for Delta_v_LOVB per (plot, variable)
  long_lovb <- rbindlist(lapply(.VAR_ORDER, function(v) {
    data.table(x = DT_c$x, y = DT_c$y,
                variable = v,
                Delta_LOVB = DT_c[[paste0("Delta_", v, "_mean")]])
  }))
  # Wide-style merge with phi
  long <- merge(DT_phi[, .(x, y, Cluster, variable, phi)],
                  long_lovb, by = c("x", "y", "variable"))
  long[, variable    := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  ann <- long[, .(rho   = suppressWarnings(cor(phi, Delta_LOVB, method = "spearman")),
                    p_val = suppressWarnings(cor.test(phi, Delta_LOVB,
                                                        method = "spearman",
                                                        exact = FALSE)$p.value)),
                by = variable]
  ann[, label := sprintf("rho_S = %+.2f\np = %s",
                           rho,
                           ifelse(p_val < 1e-3, "<0.001",
                                   sprintf("%.3f", p_val)))]
  cli_alert("Spearman per variable :")
  print(ann)

  p <- ggplot(long, aes(x = Delta_LOVB, y = phi, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_text(data = ann, aes(x = -Inf, y = Inf, label = label),
               hjust = -0.05, vjust = 1.3, size = 3.1, fontface = "bold",
               inherit.aes = FALSE, colour = "grey15") +
    facet_wrap(~ variable, nrow = 1, scales = "free",
                labeller = labeller(variable = function(x)
                                       paste("Contribution of", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    labs(
      title    = "Type A Shapley — consistency : phi_v (exact, 16 coalitions) vs Delta_v_LOVB (max-order)",
      subtitle = "1 point per plot (n=400). Dashed line = y=x identity. Spearman rho = rank consistency.",
      x = "Delta_v_LOVB  (Tmax_REF - Tmax_LOVB_v)  [K]",
      y = "phi_v  Shapley  (exact, mean over 122 d)  [K]",
      caption = "Consistency check inter-methodes. Strong rho = Shapley and LOVB attribute similar plot-level effects ; low rho = methods diverge on individual plots."
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
           plot.title      = element_text(face = "bold"),
           plot.subtitle   = element_text(colour = "grey25"),
           plot.caption    = element_text(size = 8, colour = "grey35", hjust = 0),
           legend.position = "bottom")

  fig_path <- here::here("outputs/lovb/figures/fig_typeA_Shapley_consistency.png")
  ggsave(fig_path, p, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(list(fig = fig_path, ann = ann))
}

# ---- script entry ------------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  step <- if (length(args) >= 1L) args[[1L]] else "2"
  if (step == "2") {
    build_typeB_LVA("mean")
  } else if (step == "4") {
    build_heatmap_v3()
  } else if (step == "5") {
    build_typeA_LVA("mean")
  } else if (step == "6") {
    build_typeA_Shapley_consistency()
  } else {
    stop("Unknown step : ", step)
  }
}
