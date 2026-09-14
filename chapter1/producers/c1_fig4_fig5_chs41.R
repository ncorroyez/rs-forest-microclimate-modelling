# Fig 4 (attribution, native-step ΔTmax per archetype) + Fig 5 (each lever vs its
# baseline operating point) on CHS41-Rmerge/no-wind. Also prints Table G1 (slope).
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
source("R/cluster_relabel.R")
FIG<-"out_files/Chapter1/figures"; dir.create(FIG,showWarnings=FALSE,recursive=TRUE)
th<-theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position="top")
T<-fread("out_files/Chapter1/tables/perturb_chs41_nowind.csv"); T[,P:=factor(P,levels=paste0("P",1:4))]

# ---- Fig 4: native-step ΔTmax by archetype, per variable (up=dark, down=pale) ----
L<-rbind(
  T[,.(P,var="Leaf area (±0.5)",dir="up",val=LAI_up*0.5)], T[,.(P,var="Leaf area (±0.5)",dir="down",val=LAI_dn*0.5)],
  T[,.(P,var="Cover (±10 pts)",dir="up",val=fCov_up*0.10)], T[,.(P,var="Cover (±10 pts)",dir="down",val=fCov_dn*0.10)],
  T[,.(P,var="Height (±1 m)",dir="up",val=Hmax_up)], T[,.(P,var="Height (±1 m)",dir="down",val=Hmax_dn)],
  T[,.(P,var="Profile (real−uniform)",dir="up",val=dT_LAD)])
L[,var:=factor(var,levels=c("Leaf area (±0.5)","Cover (±10 pts)","Height (±1 m)","Profile (real−uniform)"))]
g4<-ggplot(L,aes(P,val,fill=P,alpha=dir))+
  geom_hline(yintercept=0,colour="grey60",linewidth=.3)+
  geom_violin(aes(group=interaction(P,dir)),position=position_dodge(.7),width=.9,linewidth=.2,colour="grey40",draw_quantiles=c(.25,.5,.75))+
  scale_fill_manual(values=PAL_CLUSTER,guide="none")+scale_alpha_manual(values=c(up=.9,down=.4),name=NULL)+
  facet_wrap(~var,scales="free_y",nrow=1)+labs(x=NULL,y=expression(Delta*italic(T)[max]*" (°C)"))+th
ggsave(file.path(FIG,"Fig4_attribution_chs41.png"),g4,width=13,height=4,dpi=300,bg="white")

# ---- Fig 5: each lever vs baseline LAI operating point, by archetype ----
op<-rbind(
  T[,.(P,LAI,lever="Leaf area",y=LAI_up)],T[,.(P,LAI,lever="Cover",y=fCov_up)],
  T[,.(P,LAI,lever="Height",y=Hmax_up)],T[,.(P,LAI,lever="Profile",y=dT_LAD)])
op[,lever:=factor(lever,levels=c("Leaf area","Cover","Height","Profile"))]
r2<-op[,.(r2=round(summary(lm(y~LAI))$r.squared,2),r2lab=round(summary(lm(y~P))$r.squared,2)),by=lever]
g5<-ggplot(op,aes(LAI,y,colour=P))+geom_hline(yintercept=0,colour="grey70",linewidth=.3)+
  geom_point(size=.9,alpha=.5)+geom_smooth(aes(group=1),method="loess",se=FALSE,colour="grey25",linewidth=.6)+
  scale_colour_manual(values=PAL_CLUSTER,name=NULL)+
  geom_text(data=r2,aes(x=Inf,y=Inf,label=paste0("R²(op)=",sprintf('%.2f',r2),"  R²(lab)=",sprintf('%.2f',r2lab))),inherit.aes=FALSE,hjust=1.05,vjust=1.4,size=2.8)+
  facet_wrap(~lever,scales="free_y",nrow=1)+labs(x="Baseline one-sided LAI",y="lever (°C per unit)")+th
ggsave(file.path(FIG,"Fig5_operating_point_chs41.png"),g5,width=13,height=4,dpi=300,bg="white")

# ---- Table G1 (slope metric) per-step effects + CI ----
S<-fread("out_files/Chapter1/tables/perturb_chs41_nowind_SLOPE.csv")
set.seed(1);B<-2000
ci<-function(x){x<-x[is.finite(x)];m<-median(x);b<-replicate(B,median(sample(x,length(x),TRUE)));sprintf("%+.4f (%+.4f, %+.4f)",m,quantile(b,.025),quantile(b,.975))}
cat("=== Table E1 (SLOPE, CHS41): per-step MEDIAN slope effect by archetype ===\n")
for(p in c("P1","P2","P3","P4")){d<-S[P==p]
 cat(sprintf("%s: LAI+0.5 %s | fCov+10 %s | Hmax+1 %s | profile %s\n",p,ci(d$LAI_up*0.5),ci(d$fCov_up*0.10),ci(d$Hmax_up),ci(d$dT_LAD)))}
cat("DONE -> Fig4, Fig5 (chs41)\n")
