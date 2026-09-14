# ==============================================================================
# Forward-validation curves: r (zoomed to 0.5-1) and RMSE (sim vs obs across each
# cluster's HOBO loggers) as LiDAR variables are added cumulatively in each
# cluster's own importance order. One line per cluster (+ All); each point
# labelled with the variable just added; the LAI step is a diamond.
# Reads the table written by c1_forward_percluster.R (no recompute).
#   Rscript c1_forward_rcurve.R
# Out: out_files/Chapter1/figures/FigSh_forward_rcurve.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork); source("R/cluster_relabel.R")
})
d <- fread("outputs/figures_pipeline_z05/tables/tab_forward_percluster.csv")
d[, trait := sub(" \\(REF\\)", "", sub("^\\+", "", add))]
d[add == "Baseline", trait := ""]
d[, is_lai := trait == "LAI"]
d[, rowlab := factor(rowlab, levels = c("P1","P2","P3","P4","All"))]
pal <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850", All="grey30")

mk <- function(metric, ylab, ylim = NULL, drop_baseline = FALSE) {
  dd <- copy(d); dd[, val := get(metric)]
  if (drop_baseline) dd <- dd[is.finite(val)]
  set.seed(1)
  g <- ggplot(dd, aes(step, val, colour = rowlab, group = rowlab)) +
    geom_line(linewidth = 0.9) +
    geom_point(data = dd[!(is_lai)], size = 2) +
    geom_point(data = dd[(is_lai)], shape = 18, size = 4.2) +
    geom_text_repel(aes(label = trait), size = 2.8, fontface = "bold", show.legend = FALSE,
                    max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2, box.padding = 0.25) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = 0:4, labels = c("Baseline","step 1","step 2","step 3","step 4"),
                       expand = expansion(mult = c(0.04, 0.04))) +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  g + coord_cartesian(xlim = c(-0.15, 4.15), ylim = ylim)   # keep Baseline tick even if r empty there
}

pr <- mk("r", "r (sim vs obs)", ylim = c(0.5, 1.0), drop_baseline = TRUE) +
  scale_y_continuous(breaks = seq(0.5, 1.0, 0.1)) +
  annotate("text", x = 1.0, y = 0.52, hjust = 0, size = 2.7, colour = "#1A9850",
           label = "P4: +LAD alone r = -0.07 (off scale); LAI rescues it")
prm <- mk("RMSE", "RMSE (°C)")
p <- (pr / prm) + plot_layout(guides = "collect") +
  plot_annotation(caption = "Variables added cumulatively in each cluster's own importance order. Diamond = LAI step.") &
  theme(legend.position = "right", plot.caption = element_text(size = 9, colour = "grey40"))

ggsave("out_files/Chapter1/figures/FigSh_forward_rcurve.png", p, width = 9.5, height = 7.2, dpi = 200, bg = "white")
cat("DONE -> out_files/Chapter1/figures/FigSh_forward_rcurve.png\n")
