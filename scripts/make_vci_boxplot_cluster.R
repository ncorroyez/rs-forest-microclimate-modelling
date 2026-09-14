# ==============================================================================
# Boxplot of VCI (vertical complexity index) per archetype cluster P1->P4, on the
# article-canonical cLHS sample (400 plots, 100/cluster). VCI measures the
# vertical inhomogeneity of the LAD profile (NOT pure shape): the supervisors
# want to see whether the 4 LAD archetypes differ in their degree of vertical
# structuring. One box per cluster, jittered points, PAL_CLUSTER.
#   Rscript scripts/make_vci_boxplot_cluster.R
# Out: out_files/Chapter1/figures/FigSh_vci_boxplot_cluster.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })

d <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
d <- d[, .(VCI, Cluster)]
d[, P := relabel_cluster(Cluster)]

# per-cluster summary (median, IQR) for the caption / table
summ <- d[, .(n = .N, median = median(VCI, na.rm = TRUE),
              q25 = quantile(VCI, 0.25, na.rm = TRUE),
              q75 = quantile(VCI, 0.75, na.rm = TRUE),
              mean = mean(VCI, na.rm = TRUE)), by = P][order(P)]
cat("=== VCI per cluster ===\n"); print(summ)
fwrite(summ, "out_files/Chapter1/tables/tab_vci_by_cluster.csv")

set.seed(1)
p <- ggplot(d, aes(P, VCI, fill = P, colour = P)) +
  geom_boxplot(width = 0.55, alpha = 0.35, outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.12, size = 1.1, alpha = 0.45) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  scale_colour_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = NULL, y = "VCI (vertical complexity index)") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("out_files/Chapter1/figures/FigSh_vci_boxplot_cluster.png",
       p, width = 6.2, height = 5.0, dpi = 200, bg = "white")
cat("DONE -> out_files/Chapter1/figures/FigSh_vci_boxplot_cluster.png\n")
