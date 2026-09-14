# ==============================================================================
# Illustrate the v3.2.0-vs-v3.2.3 trade-off on two axes, for the case for v3.2.3:
#   AMPLITUDE  : v3.2.0 recovers more of the observed spread; v3.2.3 compresses
#                (distributions of the metric across plots, SD-recovery ratio).
#   RANKING    : v3.2.3 orders plots as well or better (rank-rank + Spearman),
#                so the Pearson "collapse" is amplitude/outlier-driven, not rank.
# Reads recap_JJAS_plotdata.csv (no recompute). Metrics: Gril slope, dTmax summer.
#   Rscript fig_v323_amplitude_vs_ranking.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
OUT_FIG <- "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/figures"
D <- fread("/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/comparaison_versions/tables/recap_JJAS_plotdata.csv")

PAL <- c("Observed"="grey30","MuSICA v3.2.0"="#0072B2","MuSICA v3.2.3"="#D55E00")
METS <- list(
  list(key="slope", obs="slope_obs", v0="slope_v320", v3="slope_v323",
       lab="Gril slope  (dT_micro/dT_macro)", ref=1, reflab="no buffering"),
  list(key="dTmaxS", obs="dTmaxS_obs", v0="dTmaxS_v320", v3="dTmaxS_v323",
       lab="dTmax summer (degC)", ref=0, reflab="= free air"))

panels <- list()
for (M in METS){
  o<-D[[M$obs]]; v0<-D[[M$v0]]; v3<-D[[M$v3]]
  # ---- amplitude: distributions + SD-recovery ratio --------------------------
  L <- rbind(data.table(src="Observed",v=o), data.table(src="MuSICA v3.2.0",v=v0), data.table(src="MuSICA v3.2.3",v=v3))
  L[, src:=factor(src, levels=names(PAL))]
  sdo<-sd(o); rec<-data.table(src=factor(c("MuSICA v3.2.0","MuSICA v3.2.3"),levels=names(PAL)),
                              r=c(sd(v0)/sdo, sd(v3)/sdo))
  pa <- ggplot(L, aes(src, v, colour=src, fill=src)) +
    { if(is.finite(M$ref)) geom_hline(yintercept=M$ref, linetype="dotted", colour="grey55") } +
    geom_violin(alpha=.12, colour=NA) +
    geom_boxplot(width=.22, alpha=0, outlier.shape=NA, linewidth=.5) +
    geom_jitter(width=.08, size=.9, alpha=.5) +
    geom_text(data=rec, aes(y=Inf, label=sprintf("SD recov.\n%.0f%%", 100*r)), vjust=1.2, size=2.8, show.legend=FALSE) +
    scale_colour_manual(values=PAL, guide="none") + scale_fill_manual(values=PAL, guide="none") +
    labs(x=NULL, y=M$lab, subtitle=sprintf("Amplitude — %s", M$key)) +
    theme_article(10) + theme(axis.text.x=element_text(angle=12, hjust=1))
  # ---- ranking: rank-rank + Spearman -----------------------------------------
  R <- rbind(data.table(model="MuSICA v3.2.0", ro=rank(o), rp=rank(v0)),
             data.table(model="MuSICA v3.2.3", ro=rank(o), rp=rank(v3)))
  R[, model:=factor(model, levels=c("MuSICA v3.2.0","MuSICA v3.2.3"))]
  sp <- data.table(model=factor(c("MuSICA v3.2.0","MuSICA v3.2.3"),levels=levels(R$model)),
                   rho=c(cor(o,v0,method="spearman"), cor(o,v3,method="spearman")),
                   r  =c(cor(o,v0), cor(o,v3)))
  pr <- ggplot(R, aes(ro, rp, colour=model)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60") +
    geom_point(size=1.5, alpha=.8) +
    geom_text(data=sp, aes(x=-Inf, y=Inf, label=sprintf("%s: rho=%.2f  (r=%.2f)", model, rho, r), colour=model),
              hjust=-0.05, vjust=c(1.4,3.0), size=2.8, show.legend=FALSE) +
    scale_colour_manual(values=PAL[2:3], name=NULL) +
    labs(x="Rank (observed)", y="Rank (predicted)", subtitle=sprintf("Ranking — %s", M$key)) +
    coord_fixed() + theme_article(10) + theme(legend.position="none")
  panels[[M$key]] <- list(a=pa, r=pr)
}

g <- (panels$slope$a | panels$slope$r) / (panels$dTmaxS$a | panels$dTmaxS$r) +
  plot_annotation(title="The v3.2.3 trade-off: loses amplitude, keeps (or improves) the ranking",
                  subtitle="Left = spread of the metric across 53 plots (v3.2.0 recovers more of the observed SD; v3.2.3 compresses). Right = rank-rank vs observed (v3.2.3 Spearman >= v3.2.0; the Pearson r drop is amplitude/outlier-driven).",
                  theme=theme(plot.subtitle=element_text(size=9, colour="grey35")))
ggsave_article(file.path(OUT_FIG,"Fig_v323_amplitude_vs_ranking"), g, 11, 8.4)
cat("DONE -> Fig_v323_amplitude_vs_ranking\n")
# quick numbers to console
for (M in METS){ o<-D[[M$obs]];v0<-D[[M$v0]];v3<-D[[M$v3]]
  cat(sprintf("[%s] SD obs=%.3f | v3.2.0=%.3f (%.0f%%) | v3.2.3=%.3f (%.0f%%) | Spearman v0=%.2f v3=%.2f | Pearson v0=%.2f v3=%.2f\n",
    M$key, sd(o), sd(v0),100*sd(v0)/sd(o), sd(v3),100*sd(v3)/sd(o),
    cor(o,v0,method="spearman"),cor(o,v3,method="spearman"), cor(o,v0),cor(o,v3))) }
