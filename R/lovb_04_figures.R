# ==============================================================================
# LOVB analysis — Module 4 : figures
#
# Six functions (one per figure type, each handles both metrics mean + P90) :
#   lovb_fig1_typeA            : 4 panels REF vs LOVB scatter, color=Cluster (cLHS)
#   lovb_fig2_typeB            : 4 panels Delta_v vs trait, color=Cluster (cLHS)
#   lovb_fig3_lectureD         : 4x4 grid (Cluster x Variable) REF vs LOVB (cLHS)
#   lovb_fig4_archetypes       : barplot Delta_v per archetype x variable
#   lovb_fig5_archetypes_scat  : 4 panels REF vs LOVB scatter, 4 archetypes (markers)
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
  library(ggplot2)
})

# Color palettes (cohérent avec figures Shapley existantes du projet)
.PAL_VAR <- c(LAI = "#440154", fCover = "#31688e",
              Hmax = "#35b779", LAD = "#fde725")

# Clusters using viridis 'turbo' family — distinct des couleurs variables
.PAL_CLUSTER <- c("1" = "#440154", "2" = "#31688e", "3" = "#35b779", "4" = "#fde725")
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Bas couvert", "4" = "C4 Inter")

.VAR_ORDER <- c("LAI", "fCover", "Hmax", "LAD")

# ---- helper : metrics aggregation per variable (RMSE/bias/pct>0.5C) ----------
.lovb_metrics_one <- function(DT_contrib, metric = c("mean", "P90"), var) {
  metric  <- match.arg(metric)
  ref_col <- paste0("Tmax_", metric, "_REF")
  lovb_col <- paste0("Tmax_", metric, "_LOVB_", var)
  diff <- DT_contrib[[ref_col]] - DT_contrib[[lovb_col]]
  list(
    rmse    = sqrt(mean(diff^2, na.rm = TRUE)),
    bias    = mean(diff, na.rm = TRUE),
    pct_05  = 100 * mean(abs(diff) > 0.5, na.rm = TRUE),
    n       = sum(!is.na(diff))
  )
}

# ==============================================================================
# Figure 4 — barplot Delta_v par archetype x variable
# ==============================================================================
lovb_fig4_archetypes <- function(DT_contrib, metric = c("mean", "P90"),
                                   out_path = NULL) {
  metric <- match.arg(metric)
  stopifnot("archetype" %in% names(DT_contrib))

  # Long format pour ggplot
  cols <- paste0("Delta_", .VAR_ORDER, "_", metric)
  long <- melt(DT_contrib[, c("archetype", "Cluster", cols), with = FALSE],
                id.vars = c("archetype", "Cluster"),
                measure.vars = cols,
                variable.name = "var_metric", value.name = "Delta")
  long[, variable := sub("^Delta_(.*)_(mean|P90)$", "\\1", var_metric)]
  long[, variable := factor(variable, levels = .VAR_ORDER)]
  long[, archetype_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_LABS)]

  ylab_metric <- if (metric == "mean") "summer mean" else "summer P90"

  p <- ggplot(long, aes(x = variable, y = Delta, fill = variable)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.5) +
    geom_col(width = 0.7, colour = "grey25", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%+.2f", Delta),
                   vjust = ifelse(Delta >= 0, -0.5, 1.3)),
              size = 3.2, fontface = "bold") +
    facet_wrap(~ archetype_lab, nrow = 1) +
    scale_fill_manual(values = .PAL_VAR, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.20, 0.20))) +
    labs(
      title = sprintf("LOVB contribution by archetype  (%s)", toupper(metric)),
      subtitle = sprintf(
        "Delta_v  =  Tmax_REF_%s  -  Tmax_LOVB_v_%s   |   positive = removing v reduces sub-canopy warming",
        metric, metric),
      caption = "Bit code REF=1111  |  LOVB_LAI=0111  LOVB_Hmax=1011  LOVB_fCover=1101  LOVB_LAD=1110",
      x = NULL,
      y = sprintf("Delta_v   Tmax (%s)  [degree C]", ylab_metric)
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text   = element_text(face = "bold"),
          plot.title   = element_text(face = "bold"),
          plot.caption = element_text(size = 8, colour = "grey35"))

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p, width = 12, height = 5.5, dpi = 300)
    cli_alert_success("Saved {.path {out_path}}")
  }
  invisible(p)
}

# ==============================================================================
# Figure 5 — REF vs LOVB scatter pour archetypes (4 panels par variable)
# ==============================================================================
lovb_fig5_archetypes_scat <- function(DT_contrib, metric = c("mean", "P90"),
                                        out_path = NULL) {
  metric <- match.arg(metric)
  stopifnot("archetype" %in% names(DT_contrib))

  ref_col <- paste0("Tmax_", metric, "_REF")
  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lovb_col <- paste0("Tmax_", metric, "_LOVB_", v)
    data.table(
      archetype   = DT_contrib$archetype,
      Cluster     = DT_contrib$Cluster,
      variable    = v,
      Tmax_REF    = DT_contrib[[ref_col]],
      Tmax_LOVB   = DT_contrib[[lovb_col]]
    )
  }))
  long[, variable := factor(variable, levels = .VAR_ORDER)]
  long[, archetype_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, archetype_lab := factor(archetype_lab, levels = .CLUSTER_LABS)]

  # Stats par variable (4 points each)
  stats <- long[, .(rmse = sqrt(mean((Tmax_REF - Tmax_LOVB)^2, na.rm = TRUE)),
                     bias = mean(Tmax_REF - Tmax_LOVB,         na.rm = TRUE)),
                 by = variable]

  lims <- range(c(long$Tmax_REF, long$Tmax_LOVB), na.rm = TRUE)
  lims <- c(min(lims) - 0.2, max(lims) + 0.2)

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB)) +
    geom_abline(slope = 1, intercept = 0,
                 linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_point(aes(colour = archetype_lab, shape = archetype_lab),
                size = 4.2, stroke = 1.1) +
    geom_text(data = stats,
               aes(x = lims[1] + 0.05 * diff(lims),
                    y = lims[2] - 0.05 * diff(lims),
                    label = sprintf("RMSE=%.2f\nbias=%+.2f", rmse, bias)),
               hjust = 0, vjust = 1, fontface = "bold", size = 3.2,
               inherit.aes = FALSE) +
    facet_wrap(~ variable, nrow = 1,
                labeller = labeller(variable = function(x) paste("Removing", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    scale_shape_manual(values = c(15, 16, 17, 18), name = NULL) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(
      title    = sprintf("REF vs LOVB at each archetype  (%s)", toupper(metric)),
      subtitle = "Points off the y=x line indicate a structural contribution to sub-canopy warming",
      x = sprintf("Tmax (REF, %s)  [degree C]", metric),
      y = sprintf("Tmax (LOVB, %s)  [degree C]", metric)
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"),
          legend.position = "bottom")

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p, width = 13, height = 5.5, dpi = 300)
    cli_alert_success("Saved {.path {out_path}}")
  }
  invisible(p)
}

# ==============================================================================
# Figure 1 — REF vs LOVB scatter (4 panels par variable) sur 400 plots cLHS
# ==============================================================================
lovb_fig1_typeA <- function(DT_contrib, metric = c("mean", "P90"),
                              out_path = NULL) {
  metric <- match.arg(metric)
  stopifnot("x" %in% names(DT_contrib))

  ref_col <- paste0("Tmax_", metric, "_REF")
  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lovb_col <- paste0("Tmax_", metric, "_LOVB_", v)
    data.table(x = DT_contrib$x, y = DT_contrib$y,
               Cluster = DT_contrib$Cluster, variable = v,
               Tmax_REF = DT_contrib[[ref_col]],
               Tmax_LOVB = DT_contrib[[lovb_col]])
  }))
  long[, variable := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  stats <- long[, .(rmse   = sqrt(mean((Tmax_REF - Tmax_LOVB)^2, na.rm = TRUE)),
                     bias   = mean(Tmax_REF - Tmax_LOVB, na.rm = TRUE),
                     pct_05 = 100 * mean(abs(Tmax_REF - Tmax_LOVB) > 0.5, na.rm = TRUE)),
                 by = variable]

  lims <- range(c(long$Tmax_REF, long$Tmax_LOVB), na.rm = TRUE)
  lims <- c(floor(lims[1]), ceiling(lims[2]))

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                 linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_point(alpha = 0.55, size = 1.8) +
    geom_text(data = stats,
               aes(x = lims[2] - 0.05 * diff(lims),
                    y = lims[1] + 0.05 * diff(lims),
                    label = sprintf("RMSE=%.2f\nbias=%+.2f\nN>0.5=%.0f%%",
                                     rmse, bias, pct_05)),
               hjust = 1, vjust = 0, fontface = "bold", size = 3.0,
               inherit.aes = FALSE) +
    facet_wrap(~ variable, nrow = 1,
                labeller = labeller(variable = function(x) paste("Removing", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(
      title    = sprintf("REF vs LOVB across 400 cLHS plots  (%s)", toupper(metric)),
      subtitle = "Each point = one plot. Deviation from y=x quantifies the contribution of removing variable v.",
      x = sprintf("Tmax (REF, %s)  [degree C]", metric),
      y = sprintf("Tmax (LOVB, %s)  [degree C]", metric)
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"),
          legend.position = "bottom")

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p, width = 14, height = 5.5, dpi = 300)
    cli_alert_success("Saved {.path {out_path}}")
  }
  invisible(p)
}

# ==============================================================================
# Figure 2 — Delta_v vs trait, LOESS, 4 panels (cLHS)
# ==============================================================================
lovb_fig2_typeB <- function(DT_contrib, metric = c("mean", "P90"),
                              out_path = NULL) {
  metric <- match.arg(metric)
  stopifnot(all(c("x", "LAI", "Hmax", "fCover", "FPC1") %in% names(DT_contrib)))

  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")
  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    delta_col <- paste0("Delta_", v, "_", metric)
    trait_col <- trait_map[v]
    data.table(
      Cluster   = DT_contrib$Cluster,
      variable  = v,
      trait     = DT_contrib[[trait_col]],
      Delta     = DT_contrib[[delta_col]]
    )
  }))
  long[, variable := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  p <- ggplot(long, aes(x = trait, y = Delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.6) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE,
                 colour = "grey20", fill = "grey70", alpha = 0.3) +
    facet_wrap(~ variable, nrow = 1, scales = "free_x",
                labeller = labeller(variable = function(x) paste("Contribution of", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         name = NULL) +
    labs(
      title    = sprintf("Contribution Delta_v vs trait  (%s)", toupper(metric)),
      subtitle = "x-axis = real value of v in the plot (FPC1 for LAD shape) | LOESS span=0.75 with 95%CI",
      x = "Trait value (real per plot)",
      y = sprintf("Delta_v   Tmax (%s)  [degree C]", metric)
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"),
          legend.position = "bottom")

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p, width = 14, height = 5.5, dpi = 300)
    cli_alert_success("Saved {.path {out_path}}")
  }
  invisible(p)
}

# ==============================================================================
# Figure 3 — Lecture D : 4x4 grid (Cluster x Variable) REF vs LOVB (cLHS)
# ==============================================================================
lovb_fig3_lectureD <- function(DT_contrib, metric = c("mean", "P90"),
                                 out_path = NULL) {
  metric <- match.arg(metric)
  stopifnot(all(c("x", "Cluster") %in% names(DT_contrib)))

  ref_col <- paste0("Tmax_", metric, "_REF")
  long <- rbindlist(lapply(.VAR_ORDER, function(v) {
    lovb_col <- paste0("Tmax_", metric, "_LOVB_", v)
    data.table(x = DT_contrib$x, y = DT_contrib$y,
               Cluster = DT_contrib$Cluster, variable = v,
               Tmax_REF = DT_contrib[[ref_col]],
               Tmax_LOVB = DT_contrib[[lovb_col]])
  }))
  long[, variable := factor(variable, levels = .VAR_ORDER)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]

  stats <- long[, .(rmse = sqrt(mean((Tmax_REF - Tmax_LOVB)^2, na.rm = TRUE))),
                 by = .(Cluster_lab, variable)]

  lims <- range(c(long$Tmax_REF, long$Tmax_LOVB), na.rm = TRUE)
  lims <- c(floor(lims[1]), ceiling(lims[2]))

  p <- ggplot(long, aes(x = Tmax_REF, y = Tmax_LOVB, colour = Cluster_lab)) +
    geom_abline(slope = 1, intercept = 0,
                 linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_text(data = stats,
               aes(x = lims[2] - 0.05 * diff(lims),
                    y = lims[1] + 0.05 * diff(lims),
                    label = sprintf("RMSE=%.2f", rmse)),
               hjust = 1, vjust = 0, fontface = "bold", size = 3.0,
               inherit.aes = FALSE) +
    facet_grid(Cluster_lab ~ variable, switch = "y",
                labeller = labeller(variable = function(x) paste("Removing", x))) +
    scale_colour_manual(values = setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS),
                         guide = "none") +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(
      title    = sprintf("Per-cluster LOVB scatter  (%s)", toupper(metric)),
      subtitle = "Rows = cluster (Open / Dense / Bas / Inter) | Cols = removed variable",
      x = sprintf("Tmax (REF, %s)  [degree C]", metric),
      y = sprintf("Tmax (LOVB, %s)  [degree C]", metric)
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text       = element_text(face = "bold"),
          strip.placement  = "outside",
          plot.title       = element_text(face = "bold"))

  if (!is.null(out_path)) {
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(out_path, p, width = 12, height = 12, dpi = 300)
    cli_alert_success("Saved {.path {out_path}}")
  }
  invisible(p)
}
