# ==============================================================================
# Chapter 3 — month-by-month boxplots, FACETED BY MONTH (April–November), with the
# six scenarios plus observed HOBO on the x-axis. Per-plot ΔTmax and slope (Tair at
# 1 m, no shift). Reads the per-plot data saved by boxplots_monthly.R (Table 4p).
# Run from z_Example root:  Rscript boxplots_monthly_byMonth.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
D <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4p_monthly_boxplot_data.csv")

src_lvl <- c("Observed (HOBO)","LiDAR full","LiDAR d_opt","S2t·ALS/ATBD","S2t·ALS/opt","S2 ATBD","S2 opt",
             "S2 ATBD ×ratio→full","S2 opt ×ratio→d_opt")
short  <- c("Observed (HOBO)"="Obs","LiDAR full"="LiDAR full","LiDAR d_opt"="LiDAR d_opt",
            "S2t·ALS/ATBD"="S2t·ALS/ATBD","S2t·ALS/opt"="S2t·ALS/opt","S2 ATBD"="ATBD","S2 opt"="opt",
            "S2 ATBD ×ratio→full"="ATBD→full","S2 opt ×ratio→d_opt"="opt→d_opt")
D[, source := factor(source, levels=src_lvl)]
D[, sx := factor(short[as.character(source)], levels=unname(short))]
D[, mo := factor(mo, levels=c("Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))]
pal <- c(Observed="grey40", LiDAR=unname(PAL_SENSOR["LiDAR"]), "Sentinel-2"=unname(PAL_SENSOR["Sentinel-2"]),
         Fusion=unname(PAL_SENSOR["Fusion"]))

mk <- function(met, ylab, href){
  ggplot(D[metric==met], aes(sx, value, fill=typ)) +
    { if(!is.na(href)) geom_hline(yintercept=href, linetype=3, colour="grey50") } +
    geom_boxplot(outlier.size=0.4, linewidth=0.3) +
    facet_wrap(~mo, ncol=4) +
    scale_fill_manual(values=pal, name=NULL) +
    labs(x=NULL, y=ylab) + theme_article(9) +
    theme(legend.position="top", axis.text.x=element_text(angle=45,hjust=1,size=7))
}
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_boxplot_byMonth_dTmax"), mk("dTmax","ΔTmax (sub-canopy − macro, °C)", 0), 10, 5.5)
ggsave_article(file.path(outdir,"Fig_boxplot_byMonth_slope"), mk("slope","slope-and-equilibrium slope", 1), 10, 5.5)
cat("DONE -> 2 by-month boxplot figures\n")
