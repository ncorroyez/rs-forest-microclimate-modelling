# ==============================================================================
# Réu 2026-06-26 #14 — Variation des traits par cluster : normalisée (z-score) +
# en différence (écart à la moyenne globale), signe positif/négatif selon cluster.
# Inclut le VCI (proxy structure verticale, même s'il n'entre pas dans MuSICA) et
# un scalaire LAD (hauteur du centre de masse relatif = top-heaviness).
#   Rscript c1_trait_variation_cluster.R
# Out: tab_trait_variation_cluster.csv + FigStation_trait_variation_cluster.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  source("R/cluster_relabel.R")
})
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
samp[, P := relabel_cluster(Cluster)]

# scalaire LAD : hauteur du centre de masse / Hmax (top-heaviness, 0=bas, 1=haut)
ladcols <- grep("^LAD_Layer_", names(samp), value=TRUE)
z <- as.numeric(gsub("LAD_Layer_", "", ladcols))
mat <- as.matrix(samp[, ..ladcols]); mat[is.na(mat)] <- 0
samp[, LAD_com := (mat %*% z)[,1] / pmax(rowSums(mat), 1e-9) / Hmax]

TRAITS <- c(LAI="LAI", Hmax="Hmax", fCover="fCover", VCI="VCI (struct. vert.)", LAD_com="LAD top-heaviness")
L <- melt(samp[, c("P", names(TRAITS)), with=FALSE], id.vars="P", variable.name="trait", value.name="val")
L <- L[is.finite(val)]
L[, z := (val - mean(val)) / sd(val), by=trait]                       # normalisé (z-score global)
L[, trait := factor(TRAITS[as.character(trait)], levels=unname(TRAITS))]

# table : moyenne native, SD, écart à la globale (différence), z-moyen par cluster
tab <- L[, .(mean=mean(val), sd=sd(val), z_mean=mean(z)), by=.(trait, P)]
glob <- L[, .(g=mean(val)), by=trait]
tab <- merge(tab, glob, by="trait"); tab[, diff_glob := mean - g]
fwrite(tab[order(trait,P)], "out_files/Chapter1/tables/tab_trait_variation_cluster.csv")
cat("=== écart à la moyenne globale (différence native) et z par cluster ===\n")
print(tab[order(trait,P), .(trait,P, mean=round(mean,3), sd=round(sd,3),
                            diff_glob=round(diff_glob,3), z_mean=round(z_mean,2))])

p <- ggplot(L, aes(P, z, fill=P, colour=P)) +
  geom_hline(yintercept=0, linetype=2, colour="grey45") +
  geom_violin(alpha=0.18, colour=NA, scale="width", width=0.85) +
  geom_boxplot(width=0.18, alpha=0.35, outlier.shape=NA) +
  facet_wrap(~trait, nrow=1) +
  scale_fill_manual(values=PAL, guide="none") + scale_colour_manual(values=PAL, guide="none") +
  labs(x=NULL, y="Variation normalisée (z-score, 0 = moyenne globale)",
       title="Variation des traits par cluster : normalisée + signe (écart à la globale)",
       subtitle="400 pixels LiDAR cLHS. La variation diffère par cluster (positive ou négative) ; VCI inclus comme proxy de structure verticale (hors-modèle)") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=7.8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_trait_variation_cluster.png", p, width=12, height=4.6, dpi=200, bg="white")
cat("DONE -> tab_trait_variation_cluster.csv + FigStation_trait_variation_cluster.png\n")
