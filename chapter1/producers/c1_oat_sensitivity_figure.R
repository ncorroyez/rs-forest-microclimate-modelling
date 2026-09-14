# ==============================================================================
# Réu 2026-06-26 #11/#13 — Sensibilité de ΔTmax à ±1 SD intra-cluster par trait.
# Depuis tab_oat_sd_cluster.csv : pour chaque archétype, ΔTmax du centroïde
# (point) et son excursion quand UN trait varie de ± la SD propre au cluster
# (barre d'erreur). Remplace le "plot moyen" par une plage de sensibilité.
#   (a) excursion ΔTmax par trait × cluster ; (b) span (=2·SD) = importance locale
#   Rscript c1_oat_sensitivity_figure.R
# Out: tab_oat_sensitivity.csv + FigStation_oat_sensitivity.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
PALc <- c(P1="#D7191C", P2="#FDAE61", P3="#74C476", P4="#1A9850")  # thermal PAL_CLUSTER (article std)
PALt <- c(LAI="#1B7837", Hmax="#BDBDBD", fCover="#7FBC41")
TR <- c("LAI","Hmax","fCover")

OFFSET <- 0.214   # ERA5 − CHS41 (cf. c1_bias_origin_chs41) : exprime ΔTmax vs station 1,5 m
oat <- fread("out_files/Chapter1/tables/tab_oat_sd_cluster.csv")
oat[, dTmax := dTmax + OFFSET]   # bases vs station (spans invariants à l'offset)
sd  <- fread("out_files/Chapter1/tables/tab_cluster_meansd.csv")
base <- oat[tag=="base", .(P, base=dTmax)]

S <- rbindlist(lapply(TR, function(t) {
  pl <- oat[tag==paste0(t,"_p"), .(P, plus=dTmax)]; mi <- oat[tag==paste0(t,"_m"), .(P, minus=dTmax)]
  m <- Reduce(function(a,b) merge(a,b,by="P"), list(base, pl, mi))
  sdv <- sd[, .(P, sdv=get(paste0(t,"_sd")))]
  m <- merge(m, sdv, by="P")
  m[, `:=`(trait=t, lo=pmin(plus,minus), hi=pmax(plus,minus),
           span=abs(plus-minus), per_unit=(plus-minus)/(2*sdv))]   # ΔTmax / unité native (#12)
  m[]
}))
S[, P := factor(P, levels=names(PALc))][, trait := factor(trait, levels=TR)]
fwrite(S[order(P,trait), .(P,trait,base,minus,plus,sd=sdv,lo,hi,span,per_unit)],
       "out_files/Chapter1/tables/tab_oat_sensitivity.csv")
cat("=== sensibilité ΔTmax à ±1 SD intra-cluster ===\n")
print(S[order(P,trait), .(P,trait, base=round(base,2), span=round(span,3), per_unit=round(per_unit,3))])

# (a) excursion ΔTmax par trait, facet cluster
bl <- base[, P:=factor(P,levels=names(PALc))]
pA <- ggplot(S, aes(trait, colour=trait)) +
  geom_hline(data=bl, aes(yintercept=base), linetype=3, colour="grey55") +
  geom_linerange(aes(ymin=lo, ymax=hi), linewidth=1.4, alpha=0.9) +
  geom_point(aes(y=base), shape=21, fill="white", size=2.4, stroke=1) +
  facet_wrap(~P, nrow=1) +
  scale_colour_manual(values=PALt, guide="none") +
  labs(x=NULL, y="ΔTmax (°C) — point = centroïde, barre = ±1 SD du trait",
       title="(a) Sensibilité de ΔTmax à ±1 SD intra-cluster, par trait et archétype",
       subtitle="ΔTmax vs station CHS 41 (1,5 m). Pointillé = centroïde (P1 amplifie → P4 tamponne). LAI = principal levier, s'érode vers le dense") +
  theme_bw(base_size=11) + theme(strip.text=element_text(face="bold"),
       plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))

# (b) span = amplitude de réponse à ±1 SD (importance locale, sans barplot)
pB <- ggplot(S, aes(P, span, colour=trait, group=trait)) +
  geom_line(linewidth=0.9, alpha=0.8) + geom_point(size=2.6) +
  scale_colour_manual(values=PALt, name=NULL) +
  labs(x=NULL, y="Excursion ΔTmax sur ±1 SD (°C)",
       title="(b) Amplitude de réponse par trait le long du gradient",
       subtitle="LAI domine et sature ; Hmax quasi-inerte ; fCover surtout en P1 ouvert") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))

ggsave("out_files/Chapter1/figures/FigStation_oat_sensitivity.png",
       pA + pB + plot_layout(widths=c(1.5,1)), width=13, height=5, dpi=200, bg="white")
cat("DONE -> tab_oat_sensitivity.csv + FigStation_oat_sensitivity.png\n")
