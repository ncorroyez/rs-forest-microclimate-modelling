# ==============================================================================
# Appendix D figure: within-archetype dispersion of the three scalar traits and
# the height-density DECOUPLING at the open end. Companion to Table D1.
#
# Shows (top) the realized per-archetype distribution of one-sided LAI, fCover
# and Hmax across the 100 cLHS plots of each archetype, and (bottom) Hmax vs
# fCover by archetype: in the open archetype P1 height spans the whole forest
# range (4-40 m) independently of cover, whereas in the dense archetypes height
# and cover are pinned together. The cLHS draw imposes no external bounds -- the
# bounds are the empirical within-cluster marginal ranges, faithfully reproduced.
# Uses the floored sample (fCover<0.5->0.5) actually passed to MuSICA, matching
# Table D1 and Appendix E.
#   Rscript c1_clhs_trait_dispersion.R
# Out: out_files/Chapter1/figures/FigAnnex_clhs_trait_dispersion.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); source("R/cluster_relabel.R")
})
FIG <- "out_files/Chapter1/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
PAL <- c(P1 = "#D9A441", P2 = "#7FBC41", P3 = "#3690C0", P4 = "#1B7837")

s <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
s[, P := relabel_cluster(Cluster)]
s[, LAI1 := LAI / 2]                                   # two-sided -> one-sided

# per-archetype within-cluster Hmax-fCover correlation (decoupling readout)
cors <- s[, .(r = cor(Hmax, fCover), n = .N, hsd = sd(Hmax)), by = P][order(P)]

viol <- function(yvar, ylab) {
  ggplot(s, aes(P, .data[[yvar]], fill = P)) +
    geom_violin(colour = NA, alpha = .55, scale = "width", width = .85) +
    geom_boxplot(width = .16, outlier.size = .5, fill = "white", alpha = .9) +
    scale_fill_manual(values = PAL, guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
}
pL <- viol("LAI1", "LAI (one-sided)")    + ggtitle("Leaf quantity")
pF <- viol("fCover", "fCover")           + ggtitle("Cover (floored at 0.5)")
pH <- viol("Hmax", expression(H[max]~(m))) + ggtitle("Canopy height")

# bottom: Hmax vs fCover -> decoupling at the open end
lab <- cors[, sprintf("%s: r = %+.2f (SD %.1f m)", P, r, hsd)]
pS <- ggplot(s, aes(fCover, Hmax, colour = P)) +
  geom_point(alpha = .5, size = 1.1) +
  scale_colour_manual(values = PAL, name = "Archetype") +
  labs(x = "fCover", y = expression(H[max]~(m)),
       title = "Height is decoupled from cover at the open end (P1), pinned in the dense end (P4)") +
  annotate("text", x = 0.50, y = c(39, 36.5, 34, 31.5), hjust = 0, vjust = 1,
           label = lab, size = 3, colour = PAL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

sub <- "cLHS imposes no external bounds: per-trait ranges are the empirical within-cluster marginals, reproduced by the design. LAI dispersion is near-constant across the gradient; fCover saturates; Hmax dispersion is far largest in the open archetype P1, where it spans the whole forest height range (4-40 m) at near-zero correlation with cover."

fig <- (pL | pF | pH) / pS +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(
    title = "Within-archetype trait dispersion and height-density decoupling",
    subtitle = sub,
    theme = theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(size = 8.3, colour = "grey35")))

ggsave(file.path(FIG, "FigAnnex_clhs_trait_dispersion.png"), fig,
       width = 11, height = 8.2, dpi = 200, bg = "white")
cat("DONE -> FigAnnex_clhs_trait_dispersion.png\n")
print(cors)
