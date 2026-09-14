# ==============================================================================
# Chapter 3 — Evolution of the buffering / amplifying plot composition ACROSS
# SCENARIOS and months. For each scenario (+ observed HOBO), counts how many of the
# 53 plots buffer (slope-and-equilibrium slope < 1) vs amplify (slope >= 1) per month
# (Apr–Nov). Reads the per-plot slopes saved by boxplots_monthly.R (Table 4p).
# Two panels: (a) monthly composition faceted by scenario; (b) pooled Apr–Nov bar.
# Run from z_Example root:  Rscript fig_ampbuff_scenarios.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")

D <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4p_monthly_boxplot_data.csv")
D <- D[metric=="slope" & is.finite(value)]
# Matched comparison: per month, keep only the plots that have an OBSERVED slope that
# month. April observations cover only 30 of 53 loggers (progressive deployment / gaps),
# so the simulations are restricted to the same 30 plots for a like-for-like composition.
obs_keys <- D[source=="Observed (HOBO)", paste(id_plot, mo)]
D <- D[paste(id_plot, mo) %in% obs_keys]
src_lvl <- c("Observed (HOBO)","LiDAR full","LiDAR d_opt","S2t·ALS/ATBD","S2t·ALS/opt","S2 ATBD","S2 opt",
             "S2 ATBD ×ratio→full","S2 opt ×ratio→d_opt")
short <- c("Observed (HOBO)"="Obs","LiDAR full"="LiDAR full","LiDAR d_opt"="LiDAR d_opt",
           "S2t·ALS/ATBD"="S2t·ALS/ATBD","S2t·ALS/opt"="S2t·ALS/opt","S2 ATBD"="ATBD","S2 opt"="opt",
           "S2 ATBD ×ratio→full"="ATBD→full","S2 opt ×ratio→d_opt"="opt→d_opt")
D[, source := factor(source, levels=src_lvl)]
D[, sx := factor(short[as.character(source)], levels=unname(short))]
D[, mo := factor(mo, levels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))]
D[, regime := fifelse(value < 1, "buffering (slope < 1)", "amplifying (slope ≥ 1)")]
D[, regime := factor(regime, levels=c("amplifying (slope ≥ 1)","buffering (slope < 1)"))]

pal <- c("buffering (slope < 1)"=unname(PAL_GRP["buffering"]),
         "amplifying (slope ≥ 1)"=unname(PAL_GRP["amplifying"]))

# (a) monthly composition, faceted by scenario
ca <- D[, .N, by=.(sx, mo, regime)]
ga <- ggplot(ca, aes(mo, N, fill=regime)) +
  geom_col(width=0.8) + facet_wrap(~sx, ncol=4) +
  scale_fill_manual(values=pal, name=NULL) +
  labs(x=NULL, y="number of plots (of 53)", tag="(a)",
       subtitle="Monthly buffering / amplifying composition per scenario (slope-and-equilibrium slope < 1 = buffering)") +
  theme_article(9) + theme(legend.position="top", axis.text.x=element_text(angle=45,hjust=1,size=7))

# (b)/(c) pooled % buffering plot-months per scenario, two windows
mk_pool <- function(sub, tag, sub_lab){
  cb <- sub[, .(pct_buff = 100*mean(value < 1), n = .N), by=sx]
  ggplot(cb, aes(sx, pct_buff)) +
    geom_col(width=0.7, fill=unname(PAL_GRP["buffering"])) +
    geom_hline(yintercept=cb[sx=="Obs"]$pct_buff, linetype=2, colour="grey30") +
    geom_text(aes(label=sprintf("%.0f%%", pct_buff)), vjust=-0.3, size=2.7) +
    scale_y_continuous(limits=c(0,100)) +
    labs(x=NULL, y="% buffering (plot-months)", tag=tag, subtitle=sub_lab) +
    theme_article(9) + theme(axis.text.x=element_text(angle=45,hjust=1))
}
gb <- mk_pool(D,                                       "(b)", "Pooled Apr–Nov; dashed = observed")
gc <- mk_pool(D[mo %in% c("Jun","Jul","Aug","Sep")],   "(c)", "Pooled Jun–Sep (summer plateau)")

g <- ga / (gb | gc) + plot_layout(heights=c(2.0,1.1))
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_ampbuff_scenarios"), g, 9.0, 8.0)
ggsave_article("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/Fig_ampbuff_scenarios", g, 9.0, 8.0)
cat("DONE -> Fig_ampbuff_scenarios\n")
print(D[, .(AprNov=100*mean(value<1)), by=sx])
print(D[mo %in% c("Jun","Jul","Aug","Sep"), .(JunSep=100*mean(value<1)), by=sx])
