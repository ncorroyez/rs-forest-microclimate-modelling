# ==============================================================================
# Figures 3/4 (and their hot-day companion S2): per-archetype attribution in native units.
#
# ONE script, TWO stages: PERIOD=all METRIC=dt is stage B2, PERIOD=hot METRIC=dt is stage B18.
#
# Reads : out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv   (stage A3)
# Writes: out_files/Chapter1/figures/Fig3_attribution_units{,_hot,_slope,_slope_hot}.{png,pdf}
#   Rscript scripts/c1_fig3_units.R                        # stage B2
#   PERIOD=hot METRIC=dt Rscript scripts/c1_fig3_units.R   # stage B18
# ==============================================================================
# Fig 3, NATIVE-UNIT design (committee 2026-07-28): the change in ΔTmax produced by
# a fixed, interpretable step of each trait — LAI ±0.5 (one-sided), Hmax ±1 m,
# fCover ±10 points — with the + and the − side shown SEPARATELY rather than averaged,
# so asymmetry (curvature) is visible. The vertical profile keeps its own footing:
# the full real-versus-uniform swap, which has no per-unit equivalent.
#
# PERIOD/METRIC selector (added 2026-07-29). The 3200 design NetCDFs were re-read by
# c1_sensitivity_units_extract2.R, which recovered the buffering slope and the hot-day
# subset the first extraction discarded. Same design, same steps, same convention:
#   PERIOD=all|hot   METRIC=dt|sl
# PERIOD=hot METRIC=dt produces the supplementary hot-day companion (Fig. S2).
suppressPackageStartupMessages({library(data.table);library(ggplot2);source("scripts/_article_style.R");source("R/cluster_relabel.R")})
LV<-c("P1","P2","P3","P4")
PERIOD<-Sys.getenv("PERIOD","all"); METRIC<-Sys.getenv("METRIC","dt")
stopifnot(PERIOD%in%c("all","hot"), METRIC%in%c("dt","sl"))
MT<-paste0(METRIC,"_",PERIOD)
U<-fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
U[,P:=factor(relabel_cluster(Cluster),levels=LV)]
# map the selected metric onto the column names the rest of the script expects
setnames(U,c(sprintf("LAI_up_%s",MT),sprintf("LAI_dn_%s",MT),sprintf("Hmax_up_%s",MT),
             sprintf("Hmax_dn_%s",MT),sprintf("fCov_up_%s",MT),sprintf("fCov_dn_%s",MT),
             sprintf("LAD_%s",MT)),
         c("LAI_up","LAI_dn","Hmax_up","Hmax_dn","fCov_up","fCov_dn","dT_LAD"))
U[,`:=`(step_LAI=0.5, step_Hmax=1.0, step_fCov=0.10)]
# per-unit slopes -> degrees for the ACTUAL step taken
# ACTUAL simulated change, not extrapolated to the nominal step: near the cover
# bounds (0.5, 1) and the LAI/Hmax floors the realised step is smaller, and
# rescaling it back to the nominal step would inflate the lever (in P4 the cover
# lever would read -0.33 instead of the -0.27 actually simulated). We therefore
# plot what was run; where a bound bites, the lever is if anything understated.
U[,`:=`(stepLp=pmin(LAI+step_LAI,Inf)-LAI,               stepLm=LAI-pmax(LAI-step_LAI,0.1),
        stepHp=step_Hmax,                                 stepHm=Hmax-pmax(Hmax-step_Hmax,3),
        stepFp=pmin(fCover+step_fCov,1)-fCover,           stepFm=fCover-pmax(fCover-step_fCov,0.5))]
U[,`:=`(LAI_plus =LAI_up *stepLp,  LAI_minus =LAI_dn *(-stepLm),
        Hmax_plus=Hmax_up*stepHp,  Hmax_minus=Hmax_dn*(-stepHm),
        fCov_plus=fCov_up*stepFp,  fCov_minus=fCov_dn*(-stepFm))]
L<-melt(U,id.vars=c("pid","P"),
        measure.vars=c("LAI_plus","LAI_minus","fCov_plus","fCov_minus","Hmax_plus","Hmax_minus","dT_LAD"),
        variable.name="lev",value.name="d")
L[,dir:=fifelse(grepl("_plus$",lev),"increase",fifelse(grepl("_minus$",lev),"decrease","swap"))]
L[,trait:=sub("_(plus|minus)$","",lev)][trait=="dT_LAD",trait:="profile"]
L[,trait:=factor(trait,levels=c("LAI","fCov","Hmax","profile"),
   labels=c("LAI  (±0.5)","fCover  (±10 pts)","Hmax  (±1 m)","profile  (real vs uniform)"))]
L[,dir:=factor(dir,levels=c("increase","decrease","swap"),
   labels=c("increase","decrease","real vs uniform"))]
p<-ggplot(L[is.finite(d)],aes(P,d,fill=P,alpha=dir))+
  geom_hline(yintercept=0,colour="grey65",linewidth=0.4)+
  geom_violin(aes(group=interaction(P,dir)),scale="width",width=0.8,colour=NA,
              position=position_dodge(width=0.85))+
  geom_boxplot(aes(group=interaction(P,dir)),width=0.13,fill="white",alpha=0.9,
               outlier.size=0.2,linewidth=0.25,position=position_dodge(width=0.85),show.legend=FALSE)+
  facet_wrap(~trait,nrow=1,scales="free_y")+
  scale_fill_manual(values=PAL_CLUSTER,guide="none")+
  scale_alpha_manual(values=c(increase=1,decrease=0.45,`real vs uniform`=0.8),name=NULL,
                     guide=guide_legend(override.aes=list(fill="grey25",colour=NA)))+
  labs(x=NULL,y=if(METRIC=="dt") expression("simulated change in "*Delta*T[max]~"(°C)") else "simulated change in buffering slope")+
  theme_article(11)+theme(legend.position="bottom")
OUT<-sprintf("out_files/Chapter1/figures/Fig3_attribution_units%s%s",
             if(METRIC=="sl") "_slope" else "", if(PERIOD=="hot") "_hot" else "")
ggsave_article(OUT,p,10,4.8)
cat(sprintf("metric=%s period=%s -> %s\n",METRIC,PERIOD,basename(OUT)))
S<-L[is.finite(d),.(med=round(median(d),3)),by=.(trait,dir,P)]
print(dcast(S,trait+dir~P,value.var="med"))
cat("DONE\n")
