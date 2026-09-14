# ==============================================================================
# Does each lever's strength track the plot's own operating point, or only its archetype?
# Per-plot sensitivity against the plot's baseline trait, with the fCover floor filtered out.
#
# Reads : out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv   (stage A3)
# Writes: outputs/figures_chap1/fig_operating_point_main.{png,pdf}
#   Rscript scripts/c1_operating_point_main.R
# ==============================================================================
# Main-text operating-point figure (tight 2x2). Committee rebuttal: trait importance
# is a smooth function of the OPERATING POINT, not the clustering. Per trait: the
# per-1-SD ΔTmax lever (scalars) or the real-vs-uniform swap (profile) vs its OWN
# baseline; points = 400 plots coloured by cluster; robust-loess smooth; each panel
# annotated with variance partition R² (operating point vs cluster label).
# fCover floor (0.5) filtered from the fCover panel. No sims. Out: fig_operating_point_main.{png,pdf}
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(mgcv)
  library(patchwork);source("R/cluster_relabel.R");source("scripts/_article_style.R")})
# CANONICAL SOURCE 2026-07-29: read the v2 table (2 metrics x 2 periods) that Figure 4
# also uses, not the v1 parts. The DeltaTmax columns are identical to 5e-15, but the v1
# parts would silently diverge the moment v2 is regenerated. Column names are mapped to
# the v1 names the rest of this script expects.
M <- fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
setnames(M, c("LAI_up_dt_all","LAI_dn_dt_all","Hmax_up_dt_all","Hmax_dn_dt_all",
              "fCov_up_dt_all","fCov_dn_dt_all","LAD_dt_all","base_dt_all"),
            c("LAI_up","LAI_dn","Hmax_up","Hmax_dn","fCov_up","fCov_dn","dT_LAD","base"))
M[, `:=`(step_LAI=0.5, step_Hmax=1.0, step_fCov=0.10)]
# native-unit design: plot the change actually simulated for the + step of each trait
M[, `:=`(stepLp=step_LAI, stepHp=step_Hmax, stepFp=pmin(fCover+step_fCov,1)-fCover)]
M[, `:=`(sLAI=LAI_up*stepLp, sHm=Hmax_up*stepHp, sfC=fCov_up*stepFp)]
M[, metric:="Tmax_all"]
M[,P:=factor(relabel_cluster(Cluster),levels=c("P1","P2","P3","P4"))]; M[,LAI1:=LAI]   # sample LAI is ALREADY one-sided (matches Table 1); no /2
#' Variance partition of one lever: smooth operating point against cluster label
#' @param d data.table of the 400 design plots
#' @param yc name of the lever column
#' @param xc name of the operating-point column
#' @return named numeric c(op, cl): GAM R-squared on x, and linear R-squared on the archetype
r2 <- function(d,yc,xc){ y<-d[[yc]]; x<-d[[xc]]
  c(op=summary(gam(y~s(x,k=5)))$r.sq, cl=summary(lm(y~d$P))$r.squared) }
#' One operating-point panel: lever against its own baseline, coloured by archetype
#' @param dd data.table of plots to draw
#' @param yc lever column name
#' @param xc operating-point column name
#' @param xlab x axis label
#' @param ylab y axis label
#' @param ttl panel subtitle
#' @return a ggplot object
pan <- function(dd,yc,xc,xlab,ylab,ttl){ d<-dd[is.finite(get(yc))&is.finite(get(xc))]
  rr<-r2(d,yc,xc); lab<-sprintf("R²: operating point %.2f | cluster %.2f", rr["op"], rr["cl"])
  ggplot(d,aes(get(xc),get(yc)))+geom_hline(yintercept=0,colour="grey80",linewidth=0.3)+
    geom_point(aes(colour=P),size=0.8,alpha=0.5)+
    geom_smooth(method="loess",span=0.8,method.args=list(family="symmetric"),colour="grey12",fill="grey65",alpha=0.22,linewidth=0.95)+
    annotate("text",x=Inf,y=Inf,label=lab,hjust=1.03,vjust=1.4,size=2.9,colour="grey25")+
    scale_colour_manual(values=PAL_CLUSTER,name="Archetype")+labs(x=xlab,y=ylab,subtitle=ttl)+theme_article(11) }
yl<-expression(Delta*T[max]~"sensitivity for the step (°C)")
pa<-pan(M,"sLAI","LAI1","baseline one-sided LAI (open → dense)",yl,"Leaf quantity")
pb<-pan(M,"sHm","Hmax",expression("baseline"~italic(H)[max]~"(m)"),yl,"Height")
pc<-pan(M[fCover>0.5001],"sfC","fCover","baseline fCover (open → closed)",yl,"Cover")
pd<-pan(M,"dT_LAD","LAI1","baseline one-sided LAI (open → dense)",expression(Delta*T[max]~"real − uniform profile contrast (°C)"),"Vertical profile, real vs uniform")
fig<-(pa|pb)/(pc|pd)+plot_layout(guides="collect")+plot_annotation(tag_levels="a",
  caption="Each lever vs its own operating point; points = 400 plots by archetype; robust-loess smooth. fCover panel excludes the 0.5 floor.")&
  theme(legend.position="bottom",plot.tag=element_text(face="bold"))
ggsave_article("outputs/figures_chap1/fig_operating_point_main",fig,9.5,8)
cat("DONE\n")
