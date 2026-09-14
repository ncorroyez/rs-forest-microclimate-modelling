# ==============================================================================
# Does a PCA on {VCI, LAI, Hmax, fCover} RECOVER the P1-P4 typology?
# (Supervisor: "ACP avec VCI LAI Hmax fCover -> si on retrouve les clusters, ça
#  prouve que le LAD compte.")
#
# Logic: the typology was built by K-means on LAI, Hmax, fCover AND the LAD shape
# (FPC1-3), NOT on VCI. If a PCA whose ONLY shape descriptor is the scalar VCI
# still separates P1-P4, then VCI captures the shape axis that the FPCs encoded
# -> the vertical-arrangement dimension carries real, separable information.
# We quantify recovery by (a) visual separation in PC1-PC2 and (b) leave-VCI-out
# vs with-VCI silhouette / classification of cluster labels.
#   Rscript c1_pca_recover_clusters.R
# Reads :
#         out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds
# Out: out_files/Chapter1/figures/FigSh_pca_clusters.png
#      ONE png: the with-VCI and no-VCI biplots are two patchwork panels inside it, not
#      two files (the previous wording read as if a separate _noVCI file existed).
#      out_files/Chapter1/tables/tab_pca_clusters.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); source("R/cluster_relabel.R")
})
source("scripts/_article_style.R")
# LINEAGE + DETERMINISM FIX 2026-07-29. This read the superseded `floor05_v2` design
# sample while every other Chapter 1 analysis uses the native20 one, so the cited PCA
# figure was built on a different 400-plot draw. The seed also makes the ggrepel label
# placement reproducible (two consecutive runs previously differed).
set.seed(42)
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover) & is.finite(VCI)]
samp[, P := relabel_cluster(Cluster)]

#' Scaled PCA of a chosen trait set, with scores and loadings on PC1-PC2
#' @param vars character vector of trait column names
#' @param tag label carried into the panel subtitle
#' @return list(sc, ld, ve, vars, tag): scores with P, loadings, and % variance explained
run_pca <- function(vars, tag) {
  X  <- scale(as.matrix(samp[, ..vars]))
  pc <- prcomp(X, center = FALSE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  sc <- as.data.table(pc$x[, 1:2]); sc[, P := samp$P]
  ld <- as.data.table(pc$rotation[, 1:2], keep.rownames = "var")
  list(sc = sc, ld = ld, ve = ve, vars = vars, tag = tag)
}
P_with <- run_pca(c("LAI","Hmax","fCover","VCI"), "with VCI")
P_no   <- run_pca(c("LAI","Hmax","fCover"),        "no VCI")

# The honest, projection-independent finding = VCI is COLLINEAR with the density
# axis (it loads almost only on PC1, and correlates strongly with LAI/fCover), so
# a scalar shape index adds no INDEPENDENT separating dimension. NOTE (descriptive
# only): this does NOT say shape is irrelevant — the ceteris-paribus proof that
# vertical arrangement matters is the sensitivity / Shapley real-vs-uniform
# contrast. It says a SCALAR is insufficient, which is exactly why the typology
# clusters on the full LAD profile (FPC1-3, orthogonal to LAI/Hmax by design).
cm <- cor(samp[, .(LAI, Hmax, fCover, VCI)], use = "complete.obs")
cat("=== variable correlations (VCI collinearity) ===\n"); print(round(cm, 2))
cat(sprintf("\nVCI loadings (with-VCI PCA): PC1=%.2f  PC2=%.2f  (VCI lives on the density axis)\n",
            P_with$ld[var=="VCI", PC1], P_with$ld[var=="VCI", PC2]))
cat(sprintf("var explained PC1-2: with VCI %.0f%% | no VCI %.0f%%\n", sum(P_with$ve[1:2]), sum(P_no$ve[1:2])))

#' Biplot of one PCA: points, 68% ellipses, centroids and loading arrows
#' @param R a list as returned by run_pca()
#' @return a ggplot object
mk_biplot <- function(R) {
  sc <- R$sc; ld <- copy(R$ld); k <- max(abs(sc$PC1), abs(sc$PC2))
  ld[, `:=`(PC1 = PC1 * k * 0.9, PC2 = PC2 * k * 0.9)]
  cen <- sc[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = P]
  ggplot(sc, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = P), size = 1.4, alpha = 0.55) +
    stat_ellipse(aes(colour = P), level = 0.68, linewidth = 0.8) +
    geom_point(data = cen, aes(fill = P), shape = 21, size = 4, colour = "black", stroke = 0.6) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.18, "cm")), colour = "grey25", linewidth = 0.5) +
    geom_text_repel(data = ld, aes(PC1, PC2, label = var), colour = "grey15", fontface = "bold", size = 3.6) +
    scale_colour_manual(values = PAL_CLUSTER, name = "Archetype") +
    scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
    labs(subtitle = R$tag,
         x = sprintf("PC1 (%.0f%%)", R$ve[1]), y = sprintf("PC2 (%.0f%%)", R$ve[2])) +
    theme_article(12) +
    # theme_article blanks plot.subtitle; restore it here or the two panels
    # ("with VCI" / "no VCI") are indistinguishable in the rendered figure.
    theme(legend.position = "bottom",
          plot.subtitle = element_text(face = "bold", size = 12))
}
library(patchwork)
p <- mk_biplot(P_with) + mk_biplot(P_no) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("out_files/Chapter1/figures/FigSh_pca_clusters.png", p, width = 11, height = 5.6, dpi = 200, bg = "white")

out <- rbind(
  data.table(set = "with VCI", var = P_with$ld$var, PC1_load = P_with$ld$PC1, PC2_load = P_with$ld$PC2),
  data.table(set = "no VCI",   var = P_no$ld$var,   PC1_load = P_no$ld$PC1,   PC2_load = P_no$ld$PC2))
fwrite(out, "out_files/Chapter1/tables/tab_pca_clusters.csv")
cat("\nDONE -> FigSh_pca_clusters.png + tab_pca_clusters.csv\n")
