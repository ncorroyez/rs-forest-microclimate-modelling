# ==============================================================================
# Réu 2026-06-26 #14 — Delta/delta NORMALISÉ par SD intra-groupe, à côté du
# delta/delta par unité native. Deux lectures complémentaires :
#   (a) ΔTmax / unité NATIVE  (pente physique ; PAS comparable entre traits)
#   (b) ΔTmax / ±1 SD intra-groupe  (signé, sans dimension d'échelle de trait →
#       DIRECTEMENT comparable LAI vs Hmax vs fCover, ET entre clusters ;
#       intègre la variation réellement disponible par cluster — le point #14).
# Dérivé des tables existantes (pas de MuSICA) :
#   - clusters : (plus−minus)/2 depuis tab_oat_sensitivity.csv
#   - global   : per_unit(All) × SD_global depuis tab_deltadelta.csv + cLHS
#   Rscript c1_deltadelta_normalized.R
# Out: tab_deltadelta_normalized.csv + FigStation_deltadelta_normalized.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
PALt <- c(LAI="#1B7837", Hmax="#9E9E9E", fCover="#7FBC41"); TR <- c("LAI","Hmax","fCover")

oat <- fread("out_files/Chapter1/tables/tab_oat_sensitivity.csv")          # P,trait,minus,plus,sd,per_unit
dd  <- fread("out_files/Chapter1/tables/tab_deltadelta.csv")               # has P=="All" per_unit
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
sd_g <- c(LAI=sd(samp$LAI), Hmax=sd(samp$Hmax), fCover=sd(samp$fCover))    # pooled global SD (matches deltadelta)

# (b) per ±1 SD, SIGNED. clusters: (plus−minus)/2 ; global: per_unit_All × SD_global
cl <- oat[, .(P=as.character(P), trait, per_unit, sd_used=sd, dT_per_SD=(plus-minus)/2)]
gl <- dd[P=="All", .(P="All", trait=as.character(trait), per_unit,
                     sd_used=sd_g[as.character(trait)], dT_per_SD=per_unit*sd_g[as.character(trait)])]
N  <- rbind(gl, cl)
N[, P := factor(P, levels=c("All","P1","P2","P3","P4"))][, trait := factor(trait, levels=TR)]
fwrite(N[order(trait,P), .(P,trait,per_unit,sd_used,dT_per_SD)],
       "out_files/Chapter1/tables/tab_deltadelta_normalized.csv")
cat("=== (b) ΔTmax par ±1 SD intra-groupe (signé, °C) ===\n")
print(dcast(N, P~trait, value.var="dT_per_SD")[order(P)])

# ---- (a) per native unit (rappel du delta/delta) -----------------------------
A <- copy(N); A[, pu := ifelse(trait=="fCover", per_unit*0.1, per_unit)]
UNI <- c(LAI="par +1 LAI", Hmax="par +1 m", fCover="par +0,1")
pA <- ggplot(A, aes(P, pu, colour=trait)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_segment(aes(xend=P, yend=0), linewidth=0.9) +
  geom_point(aes(shape=P=="All"), size=3.2) +
  scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18), guide="none") +
  facet_wrap(~trait, scales="free_y", nrow=1,
     labeller=labeller(trait=function(x) paste0(x,"\n(", UNI[x], ")"))) +
  scale_colour_manual(values=PALt, guide="none") +
  labs(x=NULL, y="ΔTmax / Δtrait (°C, unité native)",
       title="(a) Par unité native — pente physique (NON comparable entre traits)",
       subtitle="fCover paraît énorme car +1 fCover couvre tout [0,1] ; sert à comparer un MÊME trait entre clusters") +
  theme_bw(base_size=11) + theme(strip.text=element_text(face="bold",size=8.5),
       plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=7.8,colour="grey35"))

# ---- (b) per ±1 SD, shared scale → comparable across traits ------------------
pB <- ggplot(N, aes(P, dT_per_SD, colour=trait, group=trait)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_line(linewidth=0.7, alpha=0.6) +
  geom_point(aes(shape=P=="All"), size=3.2) +
  scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18), guide="none") +
  scale_colour_manual(values=PALt, name=NULL) +
  labs(x=NULL, y="ΔTmax par +1 SD intra-groupe (°C, signé)",
       title="(b) Normalisé par SD intra-groupe — COMPARABLE entre traits et clusters",
       subtitle="LAI domine et s'érode P1→P4 ; fCover 2e, surtout à l'ouvert ; Hmax ≈ 0 et change de signe. Diamant = global") +
  theme_bw(base_size=11) + theme(legend.position="bottom",
       plot.title=element_text(size=10,face="bold"), plot.subtitle=element_text(size=7.8,colour="grey35"))

ggsave("out_files/Chapter1/figures/FigStation_deltadelta_normalized.png",
       pA / pB, width=10, height=8, dpi=200, bg="white")
cat("DONE -> tab_deltadelta_normalized.csv + FigStation_deltadelta_normalized.png\n")
