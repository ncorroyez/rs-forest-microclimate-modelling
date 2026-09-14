# ==============================================================================
# Silvilaser v3ter — Type B figures (LOVB, LVA, Shapley) with LAD as violin
# instead of FPC1 scatter. Same delta values, only LAD panel visualisation changes.
#
# Layout : 3 scatter panels (LAI, fCover, Hmax) + 1 violin panel (LAD) side-by-side.
# Y-axis shared across all 4 panels for direct visual comparison.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ggplot2)
  library(patchwork)
})

OUT_DIR <- here::here("outputs/figs_silvilaser_v3ter")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Color palette matching existing pipeline
.CLUSTER_LABS <- c("1" = "C1 Open", "2" = "C2 Dense",
                    "3" = "C3 Low cover", "4" = "C4 Inter")
.PAL_CLUSTER  <- c("C1 Open" = "#E69F00", "C2 Dense" = "#0072B2",
                    "C3 Low cover" = "#009E73", "C4 Inter" = "#CC79A7")

# ---- Build long format data for a given method --------------------------------
# method : "LOVB", "LVA", or "Shapley"
# Returns data.table with cols : Cluster, Cluster_lab, variable, trait, delta
build_long_data <- function(method) {
  DT_c <- readRDS(here::here("outputs/lovb/data/DT_contrib_cLHS.rds"))
  vars <- c("LAI", "fCover", "Hmax", "LAD")
  trait_map <- c(LAI = "LAI", fCover = "fCover", Hmax = "Hmax", LAD = "FPC1")

  if (method == "LOVB") {
    rows <- lapply(vars, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("Delta_", v, "_mean")]]
    ))
  } else if (method == "LVA") {
    rows <- lapply(vars, function(v) data.table(
      Cluster = DT_c$Cluster, variable = v,
      trait   = DT_c[[trait_map[v]]],
      delta   = DT_c[[paste0("LVA_", v, "_mean")]]
    ))
  } else if (method == "Shapley") {
    DT_phi <- readRDS(here::here("outputs/lovb/data/DT_shapley_per_plot.rds"))
    # phi already in long format : x, y, Cluster, LAI, Hmax, fCover, FPC1, variable, phi
    DT_phi[, trait := NA_real_]
    for (v in vars) {
      DT_phi[variable == v, trait := get(trait_map[v])]
    }
    rows <- list(DT_phi[, .(Cluster, variable, trait, delta = phi)])
  } else stop("Unknown method : ", method)

  long <- rbindlist(rows)
  long[, variable    := factor(variable, levels = vars)]
  long[, Cluster_lab := .CLUSTER_LABS[as.character(Cluster)]]
  long[, Cluster_lab := factor(Cluster_lab, levels = .CLUSTER_LABS)]
  long[]
}

# ---- Build the 4-panel composite figure ---------------------------------------
build_typeB_v3ter <- function(method, y_label, title, fig_name) {
  cli_h1(sprintf("Type B v3ter — %s", method))
  long <- build_long_data(method)

  # Shared Y range across all 4 panels
  y_range <- range(long$delta, na.rm = TRUE)
  y_lim <- c(min(y_range[1], 0) * 1.05, max(y_range[2], 0) * 1.05)
  cli_alert("Shared Y range : [{round(y_lim[1],2)}, {round(y_lim[2],2)}]")

  # LAD distribution summary for caption
  lad_vals <- long[variable == "LAD", delta]
  cli_alert("LAD distribution :  range [{round(min(lad_vals,na.rm=TRUE),3)}, {round(max(lad_vals,na.rm=TRUE),3)}]  median {round(median(lad_vals,na.rm=TRUE),3)}")

  pal <- setNames(unname(.PAL_CLUSTER), .CLUSTER_LABS)

  # Scatter panels (LAI, fCover, Hmax)
  make_scatter <- function(v) {
    sub <- long[variable == v]
    ggplot(sub, aes(x = trait, y = delta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                  linewidth = 0.6) +
      geom_point(aes(colour = Cluster_lab), alpha = 0.55, size = 1.6) +
      geom_smooth(method = "loess", span = 0.75, se = TRUE,
                    colour = "grey20", fill = "grey70", alpha = 0.3) +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
      scale_y_continuous(limits = y_lim) +
      labs(x = "Trait value", y = y_label,
            title = paste("Contribution of", v)) +
      theme_bw(base_size = 11) +
      theme(plot.title  = element_text(face = "bold", hjust = 0.5, size = 11),
             legend.position = "none",
             panel.grid.minor = element_blank())
  }

  # Violin panel (LAD)
  sub_lad <- long[variable == "LAD"]
  p_lad <- ggplot(sub_lad, aes(x = Cluster_lab, y = delta, fill = Cluster_lab)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40",
                linewidth = 0.6) +
    geom_violin(alpha = 0.45, colour = "grey30", linewidth = 0.4,
                  scale = "width", trim = FALSE) +
    geom_jitter(aes(colour = Cluster_lab), width = 0.18, alpha = 0.6,
                  size = 0.9) +
    geom_boxplot(width = 0.12, alpha = 0.7, colour = "grey20",
                  fill = "white", outlier.shape = NA, linewidth = 0.35) +
    scale_fill_manual(values = pal, name = NULL, drop = FALSE) +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = y_lim) +
    labs(x = NULL, y = NULL, title = "Contribution of LAD (distribution)") +
    theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
           axis.text.x      = element_text(angle = 25, hjust = 1),
           axis.text.y      = element_blank(),
           axis.ticks.y     = element_blank(),
           legend.position  = "none",
           panel.grid.minor = element_blank(),
           panel.grid.major.x = element_blank())

  p1 <- make_scatter("LAI")
  p2 <- make_scatter("fCover") + theme(axis.title.y = element_blank(),
                                          axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank())
  p3 <- make_scatter("Hmax")   + theme(axis.title.y = element_blank(),
                                          axis.text.y = element_blank(),
                                          axis.ticks.y = element_blank())

  # Combined with patchwork ; common legend at the bottom
  combined <- (p1 + p2 + p3 + p_lad) +
    plot_layout(nrow = 1, widths = c(1, 1, 1, 1)) +
    plot_annotation(
      title    = title,
      subtitle = "x-axis = real value of v in the plot. LAD shown as distribution (functional object, no scalar projection). | LOESS span = 0.75 with 95%CI",
      caption  = sprintf("4 panels share Y axis [%.2f ; %.2f]. Violin: outer = density, dots = 400 plots colored by archetype, white box = IQR + median.",
                          y_lim[1], y_lim[2]),
      theme = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(colour = "grey25", size = 10),
                     plot.caption  = element_text(size = 8, colour = "grey35",
                                                    hjust = 0))
    )

  # Add a shared bottom legend (from one of the scatter panels)
  legend_plot <- make_scatter("LAI") +
    theme(legend.position = "bottom", legend.text = element_text(size = 10))
  legend <- cowplot::get_legend(legend_plot)
  final <- cowplot::plot_grid(combined, legend, ncol = 1,
                                rel_heights = c(1, 0.07))

  fig_path <- file.path(OUT_DIR, fig_name)
  # 4200 x 1650 px @ 300 dpi = 14 x 5.5 in
  ggsave(fig_path, final, width = 14, height = 5.5, dpi = 300)
  cli_alert_success("Saved {.path {fig_path}}")
  invisible(fig_path)
}

if (sys.nframe() == 0) {
  suppressMessages(library(cowplot))
  build_typeB_v3ter("LOVB",
                      "Delta_v Tmax (mean) [degree C]",
                      "Contribution Delta_v vs trait (MEAN) — LOVB",
                      "fig_typeB_LOVB_v3ter.png")
  build_typeB_v3ter("LVA",
                      "Delta_v Tmax (mean) [degree C]",
                      "Contribution Delta_v vs trait (MEAN) — LVA",
                      "fig_typeB_LVA_v3ter.png")
  build_typeB_v3ter("Shapley",
                      "phi_v Tmax (mean) [degree C]",
                      "Contribution phi_v vs trait (MEAN) — Shapley",
                      "fig_typeB_Shapley_v3ter.png")
  cli_alert_success("All 3 figures saved in {.path {OUT_DIR}}")
}
