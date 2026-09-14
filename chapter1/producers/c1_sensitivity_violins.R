# ==============================================================================
# Réu 2026-06-26 #8 — Remplacer les barplots d'importance par des VIOLONS de la
# distribution de sensibilité par-plot (400 pixels cLHS). Pour chaque trait, la
# distribution de la réponse ΔTmax par unité native, par cluster (médiane = ce
# que le barplot montrait ; le violon montre en plus la dispersion + le signe).
#   Rscript c1_sensitivity_violins.R
# Out: FigStation_sensitivity_violins.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)

S <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot","part_.*csv$",full.names=TRUE), fread), fill=TRUE)
S[, P := relabel_cluster(Cluster)]
# per-unit ΔTmax response to a +increment (native units, MÊME convention que #11/#12) :
# LAI per +1 stored-LAI (LAI_add = pas +2 → /2), Hmax per +1 m, fCover per +0.1.
# Vue par-plot, +incrément unilatéral, LAD réel — complémentaire de la pente ±SD
# centroïde (#11/#12) ; l'écart en P1 reflète la non-linéarité + LAD réel vs moyen.
S[, `:=`(LAI = LAI_add/2, Hmax = Hmax_add/5, fCover = fCov_add)]
L <- melt(S[, .(P, LAI, Hmax, fCover)], id.vars="P", variable.name="trait", value.name="s")
L <- L[is.finite(s)]
UNI <- c(LAI="ΔTmax / +1 LAI", Hmax="ΔTmax / +1 m Hmax", fCover="ΔTmax / +0,1 fCover")
L[, trait := factor(UNI[as.character(trait)], levels=unname(UNI))]

med <- L[, .(m=median(s)), by=.(trait,P)]
cat("=== médiane de sensibilité par-plot (°C/unité) — ce que résumait le barplot ===\n")
print(dcast(med, P~trait, value.var="m"))

p <- ggplot(L, aes(P, s, fill=P, colour=P)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_violin(alpha=0.18, colour=NA, scale="width", width=0.9) +
  geom_boxplot(width=0.16, alpha=0.4, outlier.shape=NA) +
  facet_wrap(~trait, scales="free_y", nrow=1) +
  scale_fill_manual(values=PAL, guide="none") + scale_colour_manual(values=PAL, guide="none") +
  labs(x=NULL, y="Sensibilité ΔTmax par-plot (°C / unité native)",
       title="Distribution de sensibilité ΔTmax par trait et cluster (violons, pas barplots)",
       subtitle="400 pixels cLHS, +incrément par-plot (LAD réel) ; même unité que #11/#12 (pente ±SD centroïde). Médiane = ancien barplot. Négatif = refroidit") +
  theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold",size=9),
       plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_sensitivity_violins.png", p, width=12, height=4.6, dpi=200, bg="white")
cat("DONE -> FigStation_sensitivity_violins.png\n")
