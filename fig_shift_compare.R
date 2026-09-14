# ==============================================================================
# Chapter 3 — figure: effect of the -2h sim time shift on the two validation
# metrics (ΔTmax & slope-and-equilibrium slope), across the 6 scenarios and the
# period set. Shows ΔTmax is shift-invariant (daily max) while the slope is
# shift-sensitive. Reads Table4j_full_shift_periods.csv.
# Run from z_Example root:  Rscript fig_shift_compare.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
R <- fread("/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4j_full_shift_periods.csv")
R[, period:=factor(period, levels=c("spring","summer","autumn","whole","summer_hot","whole_hot"))]

ml <- melt(R, id.vars=c("scenario","period","shift"),
           measure.vars=c("dT_r2","sl_r2"), variable.name="metric", value.name="R2")
ml[, metric:=ifelse(metric=="dT_r2","ΔT[max]~R^2","slope~R^2")]
ml[, shift:=factor(shift, levels=c("no shift","-2h shift"))]
ml[, typ:=ifelse(grepl("LiDAR",scenario),"LiDAR","S2")]

# paired dumbbell: no-shift vs -2h per scenario/period/metric
w <- dcast(ml, scenario+period+metric+typ ~ shift, value.var="R2")
setnames(w, c("no shift","-2h shift"), c("noshift","shift2"))

g <- ggplot(w, aes(y=scenario)) +
  geom_segment(aes(x=noshift, xend=shift2, yend=scenario, colour=typ), linewidth=0.5, alpha=0.6) +
  geom_point(aes(x=noshift, shape="no shift", colour=typ), size=1.6) +
  geom_point(aes(x=shift2,  shape="-2h shift", colour=typ), size=1.6) +
  facet_grid(period ~ metric, labeller=label_parsed) +
  scale_colour_manual(values=c(LiDAR="#1A9850", S2="#D7191C"), name=NULL) +
  scale_shape_manual(values=c("no shift"=1, "-2h shift"=16), name=NULL) +
  labs(x=expression(R^2~"(sim vs obs, 53 HOBO)"), y=NULL) +
  coord_cartesian(xlim=c(0,1)) +
  theme_article(9) +
  theme(legend.position="bottom", axis.text.y=element_text(size=7),
        panel.spacing=unit(0.4,"lines"))
outdir <- "/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_shift_compare"), g, 9, 11)
cat("DONE -> Fig_shift_compare\n")
