# ==============================================================================
# c1_radius_clhs_figure.R
#
# Appendix figure: how far the clipping radius moves the simulated buffering of
# a plot, read archetype by archetype on the 400 cLHS design plots.
#
# Each radius is compared to the 20 m reference plot by plot, so the paired
# difference removes the between-plot variance and leaves the radius effect
# alone. The grey band gives the magnitude of the leaf-quantity steps reported
# in the chapter (0.15 to 0.51 degC), which is the scale the radius effect has
# to be read against.
#
# Input : out_files/Chapter1/tables/radius_clhs_by_archetype.csv
# Output: out_files/Chapter1/figures/FigAnnex_radius_clhs.png
#   Rscript c1_radius_clhs_figure.R
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork) })
source("R/cluster_relabel.R")

tdir <- "out_files/Chapter1/tables"
fdir <- "out_files/Chapter1/figures"
dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(tdir, "radius_clhs_by_archetype.csv"))
d <- d[group != "ALL"]
d[, group := factor(group, levels = paste0("P", 1:4))]

REF      <- 20
LQ_LO    <- 0.15   # smallest leaf-quantity step reported in the chapter
LQ_HI    <- 0.51   # largest
# Sur un axe log, position_dodge() ne peut pas ecarter des intervalles inegaux :
# on decale x d'un facteur multiplicatif par archetype.
d[, xr := radius * 1.035^(as.integer(group) - 2.5)]

panel <- function(y, lo, hi, ylab, band = FALSE) {
  g <- ggplot(d, aes(xr, get(y), colour = group, group = group))
  if (band) {
    g <- g +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -LQ_HI, ymax = -LQ_LO,
               fill = "grey85", alpha = 0.5) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = LQ_LO, ymax = LQ_HI,
               fill = "grey85", alpha = 0.5) +
      annotate("text", x = Inf, y = (LQ_LO + LQ_HI) / 2, hjust = 1.05, vjust = 0.5,
               label = "leaf-quantity steps", size = 3.1, colour = "grey30")
  }
  g +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.4) +
    geom_vline(xintercept = REF, linetype = "dashed", colour = "grey55") +
    geom_linerange(aes(ymin = get(lo), ymax = get(hi)),
                   linewidth = 0.6, alpha = 0.85) +
    geom_line(linewidth = 0.6, alpha = 0.85) +
    geom_point(size = 2.1) +
    scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
    # 5 to 50 m on a linear axis crowds every label below 20 m and isolates 50 m;
    # a log axis spaces the sweep as it was designed.
    scale_x_log10(breaks = sort(unique(c(d$radius, REF))),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    labs(x = "Clipping radius (m)", y = ylab) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

g <- panel("d_dtmax", "dt_lo", "dt_hi",
           expression(paste("Paired change in ", Delta, italic(T)[max], " (°C)")),
           band = TRUE) +
     panel("d_slope", "sl_lo", "sl_hi", "Paired change in coupling slope") +
     plot_layout(ncol = 2, guides = "collect") &
     theme(legend.position = "top")

out <- file.path(fdir, "FigAnnex_radius_clhs.png")
ggsave(out, g, width = 10.5, height = 4.4, dpi = 300, bg = "white")
cat("DONE ->", out, "\n")
