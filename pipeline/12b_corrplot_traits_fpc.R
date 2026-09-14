# ==============================================================================
# ANNEX (extended) — trait correlations PLUS the 3 LAD-shape FPCs, per archetype
# (P1..P4) and over ALL 400 cLHS points. FPCs are computed once on the 400 cLHS
# LAD profiles with the same double-normalised FPCA as the typology (R/fpca.R::
# compute_fpca), then correlated against the scalar traits within each group.
#   Rscript pipeline/12b_corrplot_traits_fpc.R
# Out: out_files/Chapter1/figures/FigAnnex_corrplot_traits_fpc.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(corrplot); library(fda); library(ggplot2)
  source("R/fpca.R"); source("R/cluster_relabel.R")
})
FIG <- "out_files/Chapter1/figures"; dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

s <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
s[, P := relabel_cluster(Cluster)]

# ---- FPCA on the 400 cLHS LAD profiles (double-normalised, like the typology) --
ladcols  <- grep("^LAD_Layer_", names(s), value = TRUE)
z_breaks <- as.numeric(sub("LAD_Layer_", "", ladcols))
ord      <- order(z_breaks); ladcols <- ladcols[ord]; z_breaks <- z_breaks[ord]
mat_lad  <- as.matrix(s[, ..ladcols])
fp <- compute_fpca(mat_lad, hmax_vec = s$Hmax, z_breaks = z_breaks, n_harm = 3)
sc <- fp$fpca$scores[, 1:3]; colnames(sc) <- c("FPC1","FPC2","FPC3")
s[, c("FPC1","FPC2","FPC3") := as.data.table(sc)]
varprop <- round(100 * fp$varprop[1:3], 1)
cat(sprintf("FPCA variance explained: FPC1=%.1f%% FPC2=%.1f%% FPC3=%.1f%%\n", varprop[1], varprop[2], varprop[3]))

VARS <- c("LAI","fCover","Hmax","VCI","FPC1","FPC2","FPC3")
groups <- c(paste0("P",1:4), "All")
sub_of <- function(g) if (g == "All") s else s[P == g]

png(file.path(FIG, "FigAnnex_corrplot_traits_fpc.png"), width = 1650, height = 1150, res = 150)
op <- par(mfrow = c(2, 3), mar = c(1, 1, 2.5, 1))
for (g in groups) {
  d <- sub_of(g)[, ..VARS]
  M <- cor(d, use = "pairwise.complete.obs")
  corrplot(M, method = "color", type = "upper", diag = FALSE,
           addCoef.col = "black", number.cex = 0.7, tl.col = "black", tl.cex = 0.95,
           col = colorRampPalette(c("#2166AC","white","#B2182B"))(200), cl.pos = "n",
           mar = c(0,0,1.5,0), title = sprintf("%s  (n = %d)", g, nrow(d)))
}
plot.new()
text(0.5, 0.72, "Traits + 3 LAD-shape FPCs", cex = 1.1, font = 2)
text(0.5, 0.46, sprintf("FPC1 %.0f%% · FPC2 %.0f%% · FPC3 %.0f%%\nof LAD-profile shape variance",
                        varprop[1], varprop[2], varprop[3]), cex = 0.9, col = "grey30")
text(0.5, 0.18, "VCI = vertical complexity index", cex = 0.85, col = "grey45")
par(op); dev.off()
cat("DONE -> out_files/Chapter1/figures/FigAnnex_corrplot_traits_fpc.png\n")
for (g in groups) { cat(sprintf("\n=== %s (n=%d) ===\n", g, nrow(sub_of(g))));
  print(round(cor(sub_of(g)[, ..VARS], use="pairwise.complete.obs"), 2)) }
