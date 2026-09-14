# Fig S1 (hot-day attribution, native-step ΔTmax per archetype), CHS41-Rmerge/no-wind.
# RE-AUTHORED 2026-09-14: the shipped FigS1_attribution_hot_chs41.png was plotted interactively
# on 2026-09-08 and its script was lost (see PRODUCER_FIGURE_MAP_2026-09-14.md). This mirrors
# c1_fig4_fig5_chs41.R (the full-summer Fig 4) on the hot-day table c1_hotdays_chs41.R writes.
# Data are faithful; validate the aesthetic against the shipped png (md5 will not match).
# Reads : out_files/Chapter1/tables/hotdays_chs41.csv   (P, LAI_up, LAI_dn, fCov_up, Hmax_up, dT_LAD)
# Writes: out_files/Chapter1/figures/FigS1_attribution_hot_chs41.png  (hand-copy to chapter1/manuscript/figures/)
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("R/cluster_relabel.R")
FIG<-Sys.getenv("CH1_CHS41_FIGDIR","out_files/Chapter1/figures"); dir.create(FIG,showWarnings=FALSE,recursive=TRUE)
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="top")
T<-fread("out_files/Chapter1/tables/hotdays_chs41.csv"); T[,P:=factor(P,levels=paste0("P",1:4))]

# native-step ΔTmax by archetype, per variable. Hot-day table carries both directions for
# leaf area only; cover, height and profile are the single native "up" step (as in Fig 4).
L<-rbind(
  T[,.(P,var="Leaf area (±0.5)",dir="up",val=LAI_up*0.5)], T[,.(P,var="Leaf area (±0.5)",dir="down",val=LAI_dn*0.5)],
  T[,.(P,var="Cover (+10 pts)",dir="up",val=fCov_up*0.10)],
  T[,.(P,var="Height (+1 m)",dir="up",val=Hmax_up)],
  T[,.(P,var="Profile (real−uniform)",dir="up",val=dT_LAD)])
L[,var:=factor(var,levels=c("Leaf area (±0.5)","Cover (+10 pts)","Height (+1 m)","Profile (real−uniform)"))]
L<-L[is.finite(val)]
g<-ggplot(L,aes(P,val,fill=P,alpha=dir))+
  geom_hline(yintercept=0,colour="grey60",linewidth=.3)+
  geom_violin(aes(group=interaction(P,dir)),position=position_dodge(.7),width=.9,linewidth=.2,colour="grey40",draw_quantiles=c(.25,.5,.75))+
  scale_fill_manual(values=PAL_CLUSTER,guide="none")+scale_alpha_manual(values=c(up=.9,down=.4),name=NULL)+
  facet_wrap(~var,scales="free_y",nrow=1)+labs(x=NULL,y=expression("hot-day "*Delta*italic(T)[max]*" (°C)"))+th
ggsave(file.path(FIG,"FigS1_attribution_hot_chs41.png"),g,width=13,height=3.6,dpi=300,bg="white")

cat("=== FigS1 hot-day native-step ΔTmax: median by archetype ===\n")
for(p in c("P1","P2","P3","P4")){d<-T[P==p]
 cat(sprintf("%s: LAI+0.5 %+.3f | fCov+10 %+.3f | Hmax+1 %+.4f | profile %+.3f\n",
   p,median(d$LAI_up,na.rm=TRUE)*0.5,median(d$fCov_up,na.rm=TRUE)*0.10,
   median(d$Hmax_up,na.rm=TRUE),median(d$dT_LAD,na.rm=TRUE)))}
cat("DONE -> FigS1_attribution_hot_chs41.png\n")
