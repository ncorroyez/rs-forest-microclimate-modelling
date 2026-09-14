# ==============================================================================
# Chapter 3 — Figure 6: two-component diagnostic of the Sentinel-2 deficit for
# forcing the sub-canopy microclimate. Plots split by LiDAR LAI (open = low half,
# dense = high half). Top row: ΔTmax RMSE (magnitude accuracy); bottom row:
# ΔTmax R² (between-plot ranking). Reads Table4k_stratified.csv.
# Journal specs (target RSE): double-column figure, width ≈ 6.9 in; vector PDF;
# fonts legible at print size; colour-blind-safe sensor palette (PAL_SENSOR).
# LAI shown one-sided (LiDAR reference LAI). Run from z_Example root.
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
R <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4k_stratified.csv")
R <- R[var=="LAI_ALS"]
ord <- c("LiDAR full","LiDAR d_opt","S2 ATBD (dyn)","S2 opt (dyn)","S2 ATBD ×ratio→full","S2 opt ×ratio→d_opt")
R[, scenario:=factor(scenario, levels=rev(ord))]
R[, typ:=ifelse(grepl("LiDAR",scenario),"LiDAR","Sentinel-2")]

ml <- melt(R, id.vars=c("scenario","grp","typ"),
           measure.vars=c("dT_rmse","dT_r2"), variable.name="metric", value.name="val")
ml[, panel := fcase(
  metric=="dT_rmse" & grp=="low",  "(a) open canopy (n=26)  —  ΔTmax RMSE (°C)",
  metric=="dT_rmse" & grp=="high", "(b) dense canopy (n=27)  —  ΔTmax RMSE (°C)",
  metric=="dT_r2"  & grp=="low",   "(c) open canopy (n=26)  —  ΔTmax R²",
  metric=="dT_r2"  & grp=="high",  "(d) dense canopy (n=27)  —  ΔTmax R²")]
ml[, panel := factor(panel, levels=c(
  "(a) open canopy (n=26)  —  ΔTmax RMSE (°C)",
  "(b) dense canopy (n=27)  —  ΔTmax RMSE (°C)",
  "(c) open canopy (n=26)  —  ΔTmax R²",
  "(d) dense canopy (n=27)  —  ΔTmax R²"))]

g <- ggplot(ml, aes(scenario, val, fill=typ)) +
  geom_col(width=0.7) +
  geom_text(aes(label=sprintf("%.2f",val)), hjust=-0.15, size=2.5) +
  facet_wrap(~panel, scales="free_x", ncol=2) +
  coord_flip() +
  scale_fill_manual(values=PAL_SENSOR, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.20))) +
  labs(x=NULL, y=NULL) +
  theme_article(10) +
  theme(legend.position="top", panel.spacing=unit(0.6,"lines"),
        axis.text.y=element_text(size=8))
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_stratified_diagnostic"), g, 6.9, 5.6)
cat("DONE -> Fig_stratified_diagnostic\n")
