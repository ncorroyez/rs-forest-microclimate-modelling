# ==============================================================================
# Apply the corrected (normalized-height) VCI raster to the cLHS sample + figures.
# Re-extracts VCI from in_files/vci_res_10_m_norm.tif at the 400 cLHS pixels,
# reports the per-cluster shift (old buggy vs new), saves a corrected sample, and
# regenerates the #14 VCI-per-cluster deliverable + the VCI~LAI collinearity check.
# VCI is a DESCRIPTOR (not a MuSICA input) → no sims change.
#   Rscript c1_vci_apply_recompute.R
# Out: clhs_sample_floor05_v2_vcinorm.rds, tab_vci_by_cluster_norm.csv,
#      FigStation_trait_variation_cluster.png (corrected), printed correlation check
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(terra); library(ggplot2); source("R/cluster_relabel.R")
})
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")
stopifnot(file.exists("in_files/vci_res_10_m_norm.tif"))
r <- rast("in_files/vci_res_10_m_norm.tif")

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, VCI_old := VCI]
samp[, VCI := terra::extract(r, cbind(x, y))[, 1]]          # corrected VCI (height-normalized)
samp[, P := relabel_cluster(Cluster)]

cat("=== VCI per cluster: OLD (buggy, abs elevation) vs NEW (normalized) ===\n")
cmp <- samp[is.finite(VCI) & is.finite(VCI_old),
            .(n=.N, VCI_old=round(median(VCI_old),3), VCI_new=round(median(VCI),3),
              ratio=round(median(VCI_old)/median(VCI),3)), by=P][order(P)]
print(cmp)
cat(sprintf("ALL: old %.3f -> new %.3f (ratio %.3f); n missing new = %d\n",
            median(samp$VCI_old,na.rm=TRUE), median(samp$VCI,na.rm=TRUE),
            median(samp$VCI_old,na.rm=TRUE)/median(samp$VCI,na.rm=TRUE), sum(!is.finite(samp$VCI))))

saveRDS(samp, "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2_vcinorm.rds")
fwrite(samp[is.finite(VCI), .(VCI_mean=round(mean(VCI),3), VCI_sd=round(sd(VCI),3),
            VCI_median=round(median(VCI),3), n=.N), by=P][order(P)],
       "out_files/Chapter1/tables/tab_vci_by_cluster_norm.csv")

# ---- collinearity / conclusion check with corrected VCI ----------------------
sp <- function(a,b) cor(a,b,method="spearman",use="complete.obs")
S <- samp[is.finite(VCI)&is.finite(LAI)]
cat(sprintf("\nVCI~LAI Spearman: OLD %+.2f | NEW %+.2f  (collinearity conclusion should hold)\n",
            sp(samp$VCI_old, samp$LAI), sp(S$VCI, S$LAI)))

# ---- regenerate #14 figure with corrected VCI (same style as c1_trait_variation_cluster) ----
S <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
ladcols <- grep("^LAD_Layer_", names(S), value=TRUE); z <- as.numeric(gsub("LAD_Layer_","",ladcols))
mat <- as.matrix(S[, ..ladcols]); mat[is.na(mat)] <- 0
S[, LAD_com := (mat %*% z)[,1] / pmax(rowSums(mat),1e-9) / Hmax]
TRAITS <- c(LAI="LAI", Hmax="Hmax", fCover="fCover", VCI="VCI (struct. vert.)", LAD_com="LAD top-heaviness")
L <- melt(S[, c("P", names(TRAITS)), with=FALSE], id.vars="P", variable.name="trait", value.name="val")
L <- L[is.finite(val)]; L[, zsc := (val-mean(val))/sd(val), by=trait]
L[, trait := factor(TRAITS[as.character(trait)], levels=unname(TRAITS))]
p <- ggplot(L, aes(P, zsc, fill=P, colour=P)) +
  geom_hline(yintercept=0, linetype=2, colour="grey45") +
  geom_violin(alpha=0.18, colour=NA, scale="width", width=0.85) +
  geom_boxplot(width=0.18, alpha=0.35, outlier.shape=NA) +
  facet_wrap(~trait, nrow=1) +
  scale_fill_manual(values=PAL, guide="none") + scale_colour_manual(values=PAL, guide="none") +
  labs(x=NULL, y="Variation normalisée (z-score, 0 = moyenne globale)",
       title="Variation des traits par cluster (VCI corrigé, hauteurs normalisées)",
       subtitle="400 pixels LiDAR cLHS. VCI recalculé sur hauteurs au-dessus du sol (van Ewijk 2011) ; médiane ~0,86, comparable à Eva Gril") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=7.8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_trait_variation_cluster.png", p, width=12, height=4.6, dpi=200, bg="white")
cat("DONE -> corrected sample + tab_vci_by_cluster_norm.csv + FigStation_trait_variation_cluster.png\n")
