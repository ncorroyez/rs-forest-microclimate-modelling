# ==============================================================================
# Finalize the corrected (normalized-height) VCI: merge the per-pixel results from
# c1_vci_points_tilewise.R (vci_points_partial.csv) into the cLHS sample, write a
# corrected sample, the per-cluster VCI table, and regenerate the #14 figure.
# Runs AFTER c1_vci_points_tilewise.R has completed the 400 pixels.
#   Rscript c1_vci_finalize.R
# Out: clhs_sample_floor05_v2_vcinorm.rds, tab_vci_by_cluster_norm.csv,
#      FigStation_trait_variation_cluster.png (corrected)
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
PAL  <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")
PART <- "out_files/Chapter1/tables/vci_points_partial.csv"
stopifnot(file.exists(PART))
vci  <- unique(fread(PART), by="pid")                       # pid, vci_norm, tile

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, pid := .I][, VCI_old := VCI]
samp <- merge(samp, vci[, .(pid, vci_norm)], by="pid", all.x=TRUE)
nmiss <- sum(!is.finite(samp$vci_norm))
if (nmiss > 0) cat(sprintf("WARNING: %d/%d pixels sans VCI recalculé — recompute incomplet (relancer c1_vci_points_tilewise.R)\n",
                           nmiss, nrow(samp)))
samp[is.finite(vci_norm), VCI := vci_norm][, P := relabel_cluster(Cluster)]

cat("=== VCI par cluster : OLD (buggy) -> NEW (normalisé, r=25, tous retours) ===\n")
print(samp[is.finite(VCI)&is.finite(VCI_old),
           .(n=.N, VCI_old=round(median(VCI_old),3), VCI_new=round(median(VCI),3),
             ratio=round(median(VCI_old)/median(VCI),3)), by=P][order(P)])

saveRDS(samp, "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2_vcinorm.rds")
fwrite(samp[is.finite(VCI), .(VCI_mean=round(mean(VCI),3), VCI_sd=round(sd(VCI),3),
            VCI_median=round(median(VCI),3), n=.N), by=P][order(P)],
       "out_files/Chapter1/tables/tab_vci_by_cluster_norm.csv")

# ---- #14 figure with corrected VCI ------------------------------------------
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
       title="Variation des traits par cluster (VCI corrigé : hauteurs normalisées, van Ewijk 2011)",
       subtitle="400 pixels LiDAR cLHS. VCI recalculé sur hauteurs au-dessus du sol ; médiane ~0,86 (niveau Eva Gril, même site)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=7.6,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_trait_variation_cluster.png", p, width=12, height=4.6, dpi=200, bg="white")
cat("DONE -> sample corrigé + tab_vci_by_cluster_norm.csv + FigStation_trait_variation_cluster.png\n")
