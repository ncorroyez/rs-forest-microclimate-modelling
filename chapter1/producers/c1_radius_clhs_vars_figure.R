# ==============================================================================
# c1_radius_clhs_vars_figure.R
#
# Companion to the radius sweep: how the extracted canopy variables move with the
# clipping radius on the 400 cLHS plots, read archetype by archetype. Shows that
# the structural inputs (LAI, Hmax, fCover) are near-invariant to the footprint,
# which is why the simulated buffering is too (FigAnnex_radius_clhs).
#
# Input : out_files/Chapter1/tables/radius_clhs_vars.csv (c1_radius_clhs_vars.R)
# Output: out_files/Chapter1/figures/FigAnnex_radius_clhs_vars.png
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("R/cluster_relabel.R")

V <- fread("out_files/Chapter1/tables/radius_clhs_vars.csv")
V[, P := factor(P, levels = paste0("P", 1:4))]
REF <- 20

summ <- function(var) V[, .(med = median(get(var), na.rm = TRUE),
                            lo  = quantile(get(var), .1, na.rm = TRUE),
                            hi  = quantile(get(var), .9, na.rm = TRUE)),
                        by = .(radius, P)][order(P, radius)]

panel <- function(var, ylab) {
  d <- summ(var)
  ggplot(d, aes(radius, med, colour = P, fill = P, group = P)) +
    geom_vline(xintercept = REF, linetype = "dashed", colour = "grey55") +
    geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.12) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) +
    scale_colour_manual(values = PAL_CLUSTER, name = NULL) +
    scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
    scale_x_log10(breaks = sort(unique(V$radius)),
                  labels = function(x) format(x, drop0trailing = TRUE)) +
    labs(x = "Clipping radius (m)", y = ylab) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top")
}

g <- panel("LAI", "One-sided LAI") +
     panel("Hmax", expression(italic(H)[max]~"(m)")) +
     panel("fCover", "Fractional cover") +
     plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "top")

out <- "out_files/Chapter1/figures/FigAnnex_radius_clhs_vars.png"
ggsave(out, g, width = 12, height = 4.2, dpi = 300, bg = "white")
cat("DONE ->", out, "\n")

# median shift vs 20 m reference, for the caption
R0 <- V[radius == REF, .(id_plot, LAI0 = LAI, H0 = Hmax, f0 = fCover)]
D  <- merge(V[radius != REF], R0, by = "id_plot")
cat("\nmedian |shift| vs 20 m, pooled:\n")
cat(sprintf("  LAI  %.3f   Hmax %.2f m   fCover %.3f\n",
            median(abs(D$LAI - D$LAI0), na.rm = TRUE),
            median(abs(D$Hmax - D$H0), na.rm = TRUE),
            median(abs(D$fCover - D$f0), na.rm = TRUE)))
