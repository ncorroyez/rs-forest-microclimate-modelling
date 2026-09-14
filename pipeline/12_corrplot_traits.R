# ==============================================================================
# ANNEX — pairwise correlation of the canopy traits, per archetype (P1..P4) and
# over ALL 400 cLHS points. Shows that the trait collinearity (notably LAI-fCover)
# is real and varies with canopy density. Variables: LAI, fCover, Hmax (scalars)
# and VCI (vertical complexity index = LAD-shape scalar).
#   Rscript pipeline/12_corrplot_traits.R
# Out: out_files/Chapter1/figures/FigAnnex_corrplot_traits.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(corrplot); source("R/cluster_relabel.R") })
FIG <- "out_files/Chapter1/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
VARS <- c("LAI","fCover","Hmax","VCI")

s <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
s[, P := relabel_cluster(Cluster)]
groups <- c(paste0("P",1:4), "All")
sub_of <- function(g) if (g == "All") s else s[P == g]

png(file.path(FIG, "FigAnnex_corrplot_traits.png"), width = 1500, height = 1050, res = 150)
op <- par(mfrow = c(2, 3), mar = c(1, 1, 2.5, 1), oma = c(0, 0, 0, 0))
for (g in groups) {
  d <- sub_of(g)[, ..VARS]
  M <- cor(d, use = "pairwise.complete.obs")
  corrplot(M, method = "color", type = "upper", diag = FALSE,
           addCoef.col = "black", number.cex = 0.95, tl.col = "black", tl.cex = 1.0,
           col = colorRampPalette(c("#2166AC","white","#B2182B"))(200), cl.pos = "n",
           mar = c(0, 0, 1.5, 0), title = sprintf("%s  (n = %d)", g, nrow(d)))
}
plot.new()  # 6th cell: shared note
text(0.5, 0.6, "Pearson correlations\namong canopy traits", cex = 1.1, font = 2)
text(0.5, 0.32, "VCI = vertical complexity\nindex (LAD-shape scalar)", cex = 0.95, col = "grey30")
par(op); dev.off()
cat("DONE -> out_files/Chapter1/figures/FigAnnex_corrplot_traits.png\n")
# also print the matrices for the record
for (g in groups) { cat(sprintf("\n=== %s (n=%d) ===\n", g, nrow(sub_of(g))));
  print(round(cor(sub_of(g)[, ..VARS], use="pairwise.complete.obs"), 2)) }
