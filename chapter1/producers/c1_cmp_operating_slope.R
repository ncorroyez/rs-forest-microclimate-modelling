# ==============================================================================
# Operating point on the BUFFERING SLOPE (REF-canopy micro/macro slope) —
# ERA5 vs FR-Blo. Companion to cmp03/cmp04 (which are on base ΔTmax).
# Design-independent (REF sim only, no ±SD perturbation).
#   cmp03b_operating_point_slope       (median by cluster)
#   cmp04b_operating_point_slope_perplot (per-plot scatter)
# → out_files/Chapter1/figures/cmp_era5_frblo/
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2)
  source("R/cluster_relabel.R"); source("scripts/_article_style.R") })
OUT<-"out_files/Chapter1/figures/cmp_era5_frblo"; palF<-c("ERA5 (original)"="#7570B3","FR-Blo (new)"="#E66101")
ERA5_VER<-Sys.getenv("ERA5_VER","iter")
bslope<-function(ver,lab){m<-rbindlist(lapply(list.files(sprintf("out_files/Chapter1/tables/metrics6_%s",ver),"part.*csv",full.names=T),fread),fill=T)
  m<-m[metric=="slope_all"]; m[,P:=as.character(relabel_cluster(Cluster))]; m[,.(pid,P,slope=base)][,src:=lab][]}
E<-bslope(ERA5_VER,"ERA5 (original)"); F<-bslope("frblo","FR-Blo (new)")

# cmp03b: median base slope by cluster
B<-rbind(E[,.(slope=median(slope)),by=P][,src:="ERA5 (original)"],F[,.(slope=median(slope)),by=P][,src:="FR-Blo (new)"]); B[,P:=factor(P,levels=c("P1","P2","P3","P4"))]
p3<-ggplot(B,aes(P,slope,fill=src))+geom_hline(yintercept=1,colour="grey55",linetype=2)+geom_col(position=position_dodge(0.7),width=0.62)+
  scale_fill_manual(values=palF,name=NULL)+coord_cartesian(ylim=c(min(B$slope)*0.98,1.02))+
  labs(x=NULL,y="REF buffering slope (micro/macro)",subtitle="Operating point — REF-canopy buffering slope by archetype (<1 = buffering). Design-independent.")+
  theme_article(12)+theme(legend.position="bottom")
ggsave_article(file.path(OUT,"cmp03b_operating_point_slope"),p3,6.5,5)

# cmp04b: per-plot scatter
W<-merge(E[,.(pid,P,se=slope)],F[,.(pid,sf=slope)],by="pid"); r<-cor(W$se,W$sf,use="complete.obs"); lim<-range(c(W$se,W$sf,1),na.rm=T)
p4<-ggplot(W,aes(se,sf,colour=P))+geom_abline(slope=1,intercept=0,linetype=2,colour="grey55")+
  geom_hline(yintercept=1,linetype=3,colour="grey80")+geom_vline(xintercept=1,linetype=3,colour="grey80")+
  geom_point(size=1.6,alpha=0.75)+scale_colour_manual(values=PAL_CLUSTER,name=NULL)+coord_equal(xlim=lim,ylim=lim)+
  labs(x="REF slope ERA5",y="REF slope FR-Blo",subtitle=sprintf("REF buffering slope per plot: ERA5 vs FR-Blo (r=%.3f). Design-independent.",r))+
  theme_article(12)+theme(legend.position="bottom")
ggsave_article(file.path(OUT,"cmp04b_operating_point_slope_perplot"),p4,6,6.3)
cat(sprintf("median REF slope by cluster:\n"));print(dcast(B,P~src,value.var="slope"))
cat(sprintf("per-plot REF slope agreement r=%.3f\n",r))
cat("DONE -> cmp03b + cmp04b\n")
