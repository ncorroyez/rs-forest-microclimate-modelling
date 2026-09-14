# ==============================================================================
# Compact 2-panel figure for the methods-justification tests (SupMat).
#  (a) Test 1: baseline ΔTmax shift when using MEDIAN instead of MEAN central
#      values (per cluster). Near zero except open P1 (fCover 0.5 floor).
#  (b) Test 2: Jensen gap f(mean archetype) - mean_i f(real_i) per cluster.
#      Near zero except open P1 -> justifies cLHS-distribution as MAIN.
#   Rscript c1_methodo_figure.R
# Out: out_files/Chapter1/figures/FigSh_methodo_tests.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); source("R/cluster_relabel.R")
})
t1 <- fread("out_files/Chapter1/tables/tab_methodo_test1_median_vs_mean.csv")[, .(Cluster, val = dTmax)]
t2 <- fread("out_files/Chapter1/tables/tab_methodo_test2_clhs_vs_archetype.csv")[, .(Cluster, val = gap_Tmax)]
t1[, P := relabel_cluster(Cluster)]; t2[, P := relabel_cluster(Cluster)]

bar <- function(d, ttl, ylab) ggplot(d, aes(P, val, fill = P)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_col(width = 0.66, colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%+.2f", val), vjust = ifelse(val >= 0, -0.4, 1.3)), size = 3.4) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  labs(title = ttl, x = NULL, y = ylab) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11, face = "bold"))

p <- bar(t1, "(a) Baseline: median - mean central values", expression(Delta*T[max]~shift~("°C"))) +
     bar(t2, "(b) Jensen: f(mean archetype) - mean f(plots)", expression(Delta*T[max]~gap~("°C")))
p <- p + plot_annotation(
  caption = "Immaterial in saturated dense clusters; only open, steep-response P1 is sensitive. Keep mean baseline (flag P1); cLHS distribution as MAIN.") &
  theme(plot.caption = element_text(size = 9, colour = "grey40", hjust = 0))
ggsave("out_files/Chapter1/figures/FigSh_methodo_tests.png", p, width = 10, height = 4.6, dpi = 200, bg = "white")
cat("DONE -> FigSh_methodo_tests.png\n")
