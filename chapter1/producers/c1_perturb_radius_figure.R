# ==============================================================================
# c1_perturb_radius_figure.R — footprint robustness of the Ch1 attribution.
# The per-unit levers (leaf quantity, cover, height, real-vs-uniform profile) from
# the 8-sim perturbation, read across the 7 clipping radii, on the plots present at
# ALL radii (constant sample). Shows the attribution ordering is invariant to the
# footprint. Input: sensitivity_radius_clhs/r*.csv. Output: FigAnnex_perturb_radius.png
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
source("R/cluster_relabel.R")
d <- "out_files/Chapter1/tables/sensitivity_radius_clhs"
RAD <- c(5,10,12.5,15,20,25,50); REF <- 20
fs <- setNames(file.path(d, sprintf("r%s.csv", gsub("\\.","p",RAD))), RAD)
L <- rbindlist(lapply(names(fs), function(r){ x<-fread(fs[[r]]); x[,radius:=as.numeric(r)]; x}), fill=TRUE)
# constant sample: plots present at every radius
common <- L[, .N, by=id_plot][N==length(RAD), id_plot]
L <- L[id_plot %in% common]; L[, P := factor(P, levels=paste0("P",1:4))]
cat("common plots across all", length(RAD), "radii:", length(common), "\n")

lev <- c(LAI_up="Leaf-quantity lever (°C per 0.5 LAI)",
         fCov_up="Cover lever (°C per 10 pts)",
         Hmax_up="Height lever (°C per m)",
         dT_LAD="Profile: real − uniform (°C)")
panel <- function(col, ylab){
  s <- L[, .(med=median(get(col),na.rm=TRUE),
             lo=quantile(get(col),.25,na.rm=TRUE), hi=quantile(get(col),.75,na.rm=TRUE)),
         by=.(radius,P)][order(P,radius)]
  ggplot(s, aes(radius, med, colour=P, fill=P, group=P)) +
    geom_hline(yintercept=0, colour="grey55", linewidth=.3) +
    geom_vline(xintercept=REF, linetype="dashed", colour="grey55") +
    geom_ribbon(aes(ymin=lo,ymax=hi), colour=NA, alpha=.12) +
    geom_line(linewidth=.6) + geom_point(size=1.7) +
    scale_colour_manual(values=PAL_CLUSTER,name=NULL) + scale_fill_manual(values=PAL_CLUSTER,guide="none") +
    scale_x_log10(breaks=RAD, labels=function(x)format(x,drop0trailing=TRUE)) +
    labs(x="Clipping radius (m)", y=ylab) +
    theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="top")
}
g <- panel("LAI_up",lev["LAI_up"]) + panel("fCov_up",lev["fCov_up"]) +
     panel("Hmax_up",lev["Hmax_up"]) + panel("dT_LAD",lev["dT_LAD"]) +
     plot_layout(ncol=2, guides="collect") & theme(legend.position="top")
out <- "out_files/Chapter1/figures/FigAnnex_perturb_radius.png"
ggsave(out, g, width=10, height=7, dpi=300, bg="white")
cat("DONE ->", out, "\n")

# paired median |shift| vs 20 m, pooled, per lever
R0 <- L[radius==REF, .(id_plot, l0=NA_real_)]
for(col in names(lev)){
  m <- merge(L[radius!=REF,.(id_plot,radius,v=get(col))], L[radius==REF,.(id_plot,v0=get(col))], by="id_plot")
  cat(sprintf("  %-8s median |lever(r) − lever(20m)| = %.3f\n", col, median(abs(m$v-m$v0),na.rm=TRUE)))
}
