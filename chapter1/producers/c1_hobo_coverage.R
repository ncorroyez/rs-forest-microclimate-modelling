# ==============================================================================
# Coverage figure: where do the 53 HOBO loggers fall in the cLHS trait space?
# Shows the validation loggers span the full structural gradient INCLUDING the
# dense end (answers the Reviewer-2 question raised when the logger→archetype
# units bug was found). cLHS = 400 design plots (background, by archetype);
# HOBO = 53 loggers (canonical clusters.rds assignment, foreground).
# LAI put on a common ONE-SIDED scale (cLHS column is two-sided -> /2).
#   Rscript c1_hobo_coverage.R
# Out: FigSh_hobo_coverage.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); source("R/cluster_relabel.R")
})
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, P := relabel_cluster(Cluster)]
samp[, LAI1 := LAI / 2]                                   # two-sided -> one-sided
H <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")  # canonical P, LAI one-sided
setnames(H, "LAI", "LAI1")

PAL <- c(P1 = "#D9A441", P2 = "#7FBC41", P3 = "#3690C0", P4 = "#1B7837")
base <- function(xvar, xlab, yvar, ylab) {
  ggplot() +
    geom_point(data = samp, aes(.data[[xvar]], .data[[yvar]], colour = P), alpha = .18, size = .9) +
    geom_point(data = H, aes(.data[[xvar]], .data[[yvar]], fill = P), shape = 21, colour = "black",
               size = 2.4, stroke = .4) +
    scale_colour_manual(values = PAL, name = "cLHS archetype") +
    scale_fill_manual(values = PAL, name = "HOBO (53)") +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
}
pA <- base("LAI1", "LAI (one-sided)", "Hmax", expression(H[max]~(m))) + ggtitle("LAI × height")
pB <- base("LAI1", "LAI (one-sided)", "fCover", "fCover")             + ggtitle("LAI × cover")
pC <- base("fCover", "fCover", "Hmax", expression(H[max]~(m)))        + ggtitle("cover × height")

# annotate dense coverage
nd <- nrow(H[LAI1 >= 4]); cLHSp4_lai <- round(mean(samp[P=="P4", LAI1]), 1)
sub <- sprintf("The 53 loggers span the gradient incl. the dense end: %d have one-sided LAI ≥ 4 (cLHS-P4 mean ≈ %.1f). Validation is not confined to open stands.",
               nd, cLHSp4_lai)
fig <- (pA + pB + pC) + plot_layout(guides = "collect") +
  plot_annotation(title = "Validation loggers cover the full structural gradient",
                  subtitle = sub,
                  theme = theme(plot.subtitle = element_text(size = 8.5, colour = "grey35"),
                                plot.title = element_text(face = "bold"))) &
  theme(legend.position = "bottom")
ggsave("out_files/Chapter1/figures/FigSh_hobo_coverage.png", fig,
       width = 15, height = 4.8, dpi = 200, bg = "white")
cat(sprintf("DONE -> FigSh_hobo_coverage.png | HOBO with one-sided LAI>=4: %d/53 ; cLHS-P4 mean LAI %.1f\n",
            nd, cLHSp4_lai))
