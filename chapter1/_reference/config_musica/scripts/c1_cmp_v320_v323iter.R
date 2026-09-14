# ==============================================================================
# Comparison figures: legacy v3.2.0 (no yoyo, validated) vs v3.2.3 + ABL_flag='iter'
# (yoyo, real MERRA-2 PBLH). Two figures:
#   (1) HOBO validation: simulated vs observed ΔTmax per plot, per version.
#   (2) Per-unit trait-importance per archetype, per version.
# Inputs: tab_hobo_iter_validation.csv, tab_version_iter_compare.csv (already built).
#   Rscript c1_cmp_v320_v323iter.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
FIG <- "out_files/Chapter3/figures"; VERS <- c("v3.2.0", "v3.2.3 iter")

# ---- Figure 1: HOBO validation (sim vs obs) ----------------------------------
H <- fread("out_files/Chapter3/tables/tab_hobo_iter_validation.csv")
L <- melt(H, id.vars=c("id_plot","obs"), measure.vars=VERS, variable.name="version", value.name="sim")
L <- L[is.finite(sim) & is.finite(obs)]
lab <- L[, .(r=cor(sim,obs), bias=mean(sim-obs), rmse=sqrt(mean((sim-obs)^2)), n=.N), by=version]
lab[, txt := sprintf("r = %.2f\nbias = %+.2f °C\nRMSE = %.2f\nn = %d", r, bias, rmse, n)]
rng <- range(c(L$sim, L$obs), na.rm=TRUE)
p1 <- ggplot(L, aes(obs, sim)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey55") +
  geom_smooth(method="lm", se=FALSE, colour="#D7791B", linewidth=0.8) +
  geom_point(alpha=0.7, size=1.6, colour="#2C3E50") +
  geom_text(data=lab, aes(x=rng[1], y=rng[2], label=txt), hjust=0, vjust=1, size=3, colour="grey20") +
  facet_wrap(~version, nrow=1) + coord_equal(xlim=rng, ylim=rng) +
  labs(x="Observed ΔT[max] (HOBO, °C)", y="Simulated ΔT[max] (1 m, °C)",
       title="HOBO validation: legacy v3.2.0 (no yoyo) vs v3.2.3 yoyo/iter (real BLH)",
       subtitle="53 loggers, REF (real z05 canopy). The validated v3.2.0 keeps the between-plot ranking (r=0.93); the yoyo binary does not recover it.") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11,face="bold"), plot.subtitle=element_text(size=8.3,colour="grey35"))
ggsave(file.path(FIG,"FigCmp_hobo_validation_v320_v323iter.png"), p1, width=9, height=5, dpi=200, bg="white")

# ---- Figure 2: per-unit trait importance ------------------------------------
I <- fread("out_files/Chapter3/tables/tab_version_iter_compare.csv")
I[version == "v3.2.3 iter (realBLH)", version := "v3.2.3 iter"]   # match the comparison label
I <- I[version %in% VERS & P != "All"]
M <- melt(I, id.vars=c("P","version"), variable.name="trait", value.name="imp")
M[, trait := factor(trait, levels=c("LAI","fCover","LAD","Hmax"))]
M[, version := factor(version, levels=VERS)]
PAL <- c(LAI="#1B7837", fCover="#7FBC41", LAD="#762A83", Hmax="#BDBDBD")
p2 <- ggplot(M, aes(P, imp, fill=trait)) +
  geom_col(position=position_dodge(0.8), width=0.74) +
  geom_text(aes(label=sprintf("%.2f",imp)), position=position_dodge(0.8), vjust=-0.3, size=2.3, colour="grey30") +
  facet_wrap(~version, nrow=1) + scale_fill_manual(values=PAL, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.10))) +
  labs(x=NULL, y="ΔT[max] importance (°C)",
       title="Per-unit trait importance per archetype: v3.2.0 vs v3.2.3 yoyo/iter",
       subtitle="Importance = |median per-unit sensitivity| × within-archetype dispersion (LAD = |real-vs-uniform|). 400 cLHS plots.") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), legend.position="bottom",
        strip.text=element_text(face="bold"), plot.title=element_text(size=11,face="bold"),
        plot.subtitle=element_text(size=8.3,colour="grey35"))
ggsave(file.path(FIG,"FigCmp_importance_v320_v323iter.png"), p2, width=9.5, height=5, dpi=200, bg="white")

cat("=== validation metrics ===\n"); print(lab[,.(version,r=round(r,3),bias=round(bias,2),rmse=round(rmse,2),n)])
cat("\n=== importance (per archetype) ===\n"); print(I)
cat("DONE -> FigCmp_hobo_validation_v320_v323iter.png + FigCmp_importance_v320_v323iter.png\n")
