# ==============================================================================
# FAST corrected VCI at the 400 cLHS pixels (bypasses the slow wall-to-wall raster).
# Clip r=25 m at each cLHS pixel (SAME radius as the LAD/LAI extraction → consistent),
# normalize_height(tin()), VCI on ALL returns above ground (van Ewijk 2011).
# Replaces the buggy absolute-elevation VCI; VCI is a DESCRIPTOR, not a MuSICA input.
#   Rscript c1_vci_points_recompute.R
# Out: clhs_sample_floor05_v2_vcinorm.rds, tab_vci_by_cluster_norm.csv,
#      FigStation_trait_variation_cluster.png (corrected), per-cluster old->new printout
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(data.table); library(parallel); library(ggplot2); source("R/cluster_relabel.R")
})
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")
CTG <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
ctg <- readLAScatalog(CTG, select = "xyzc"); opt_progress(ctg) <- FALSE   # ALL returns (no return filter)

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, pid := .I]
cat(sprintf("clipping r=25 m at %d cLHS pixels ...\n", nrow(samp)))
clips <- clip_circle(ctg, samp$x, samp$y, 25)
if (inherits(clips, "LAS")) clips <- list(clips)

vci_one <- function(pt) {
  if (is.null(pt) || npoints(pt) < 200) return(NA_real_)
  ptn <- tryCatch(normalize_height(pt, tin()), error = function(e) NULL)
  if (is.null(ptn)) return(NA_real_)
  z <- ptn$Z; z <- z[is.finite(z) & z >= 0]
  if (length(z) < 200) return(NA_real_)
  zmax <- max(z); if (zmax < 2) return(NA_real_)
  tryCatch(VCI(z, zmax = zmax), error = function(e) NA_real_)
}
vci_norm <- unlist(mclapply(clips, vci_one, mc.cores = 4))
samp[, VCI_old := VCI][, VCI := vci_norm[pid]][, P := relabel_cluster(Cluster)]

cat("\n=== VCI par cluster : OLD (buggy abs elevation) -> NEW (normalisé, r=25, tous retours) ===\n")
cmp <- samp[is.finite(VCI) & is.finite(VCI_old),
            .(n=.N, VCI_old=round(median(VCI_old),3), VCI_new=round(median(VCI),3),
              ratio=round(median(VCI_old)/median(VCI),3)), by=P][order(P)]
print(cmp)
cat(sprintf("ALL: old %.3f -> new %.3f (ratio %.3f) ; n manquants new = %d / %d\n",
            median(samp$VCI_old,na.rm=TRUE), median(samp$VCI,na.rm=TRUE),
            median(samp$VCI_old,na.rm=TRUE)/median(samp$VCI,na.rm=TRUE),
            sum(!is.finite(samp$VCI)), nrow(samp)))

saveRDS(samp, "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2_vcinorm.rds")
fwrite(samp[is.finite(VCI), .(VCI_mean=round(mean(VCI),3), VCI_sd=round(sd(VCI),3),
            VCI_median=round(median(VCI),3), n=.N), by=P][order(P)],
       "out_files/Chapter1/tables/tab_vci_by_cluster_norm.csv")

sp <- function(a,b) cor(a,b,method="spearman",use="complete.obs")
S0 <- samp[is.finite(VCI)&is.finite(LAI)]
cat(sprintf("\nVCI~LAI Spearman: OLD %+.2f | NEW %+.2f (colinéarité ~conservée)\n",
            sp(samp$VCI_old, samp$LAI), sp(S0$VCI, S0$LAI)))

# ---- regenerate #14 figure with corrected VCI -------------------------------
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
       title="Variation des traits par cluster (VCI corrigé : hauteurs normalisées, r=25 m, tous retours)",
       subtitle="400 pixels LiDAR cLHS. VCI van Ewijk 2011 sur hauteurs au-dessus du sol ; médiane ~0,86, comparable à Eva Gril (même site)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=7.6,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_trait_variation_cluster.png", p, width=12, height=4.6, dpi=200, bg="white")
cat("DONE -> corrected sample + tab_vci_by_cluster_norm.csv + FigStation_trait_variation_cluster.png\n")
