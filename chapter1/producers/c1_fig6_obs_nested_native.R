# ==============================================================================
# Model-free nested R2 on the 53 loggers: what does each LiDAR variable add to OBSERVED
# buffering, with no simulation involved? Profile shape enters last, and adds nothing.
#
# Reads : in_files/data_Blois_utm31n.geojson / in_files_native20/*.tif  (logger traits)
#         in_files/lad_z05/Blois_lad_z05_r25.csv                        (profile shape)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv        (stage A2: dTmax_obs)
# Writes: out_files/Chapter1/figures/Fig6_obs_nested_native.{png,pdf}
#   Rscript scripts/c1_fig6_obs_nested_native.R
# ==============================================================================
# Fig 6 (model-free nested regression) — NATIVE logger traits.
# Committee 2026-07-28: show a TRUE NESTED PROGRESSION rather than four unrelated
# arms: LAI -> LAI+fCover -> LAI+fCover+Hmax -> full (+VCI), so each bar adds one
# trait to the previous model and the increment is readable directly.
# The structure-only arm is kept as a reference point (grey) because the text
# quotes it, but it is set apart from the nested sequence.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggtext);library(terra);library(sf);source("scripts/_article_style.R")})
g<-st_read("in_files/data_Blois_utm31n.geojson",quiet=TRUE);r0<-rast("in_files_native20/lai_z1_res_10_m.tif");g<-st_transform(g,crs(r0));v<-vect(g)
D<-data.table(id_plot=g$id_plot,
  LAI=terra::extract(rast("in_files_native20/lai_z1_res_10_m.tif"),v)[,2],
  fCover=terra::extract(rast("in_files_native20/fCover_res_10_m.tif"),v)[,2],
  Hmax=terra::extract(rast("in_files_native20/max_res_10_m.tif"),v)[,2],
  VCI=terra::extract(rast("in_files_native20/vci_res_10_m.tif"),v)[,2])
# vertical-profile SHAPE at each logger: relative height of the LAD centroid
# (same variable as Fig. F3). VCI is deliberately NOT used here: it is a density
# descriptor, not a measure of profile shape, so a VCI-based arm would not test
# the vertical dimension at all (committee 2026-07-28).
DZ<-0.5
lad<-fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
lyr<-grep("^LAD_Layer_",names(lad),value=TRUE); hh<-as.numeric(sub("LAD_Layer_","",lyr))
Lm<-as.matrix(lad[,..lyr]); Lm[!is.finite(Lm)]<-0
lad[,cen:=(Lm%*%hh)[,1]/pmax(rowSums(Lm),1e-9)]
lad[,profile:=cen/pmax(Hmax,1e-9)]                       # 0-1, higher = top-heavy
oc<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot,dTmax_obs)]
oc<-merge(oc,lad[,.(id_plot,profile)],by="id_plot")
D<-merge(D,oc,by="id_plot");D<-D[is.finite(dTmax_obs)&is.finite(LAI)&is.finite(fCover)&is.finite(Hmax)&is.finite(profile)];D[,y:=-dTmax_obs]
#' R-squared of a linear model on the merged 53-logger table D
#' @param f model formula
#' @return numeric R-squared
r2<-function(f)summary(lm(f,D))$r.squared
#' Adjusted R-squared of a linear model on D
#' @param f model formula
#' @return numeric adjusted R-squared
adj<-function(f)summary(lm(f,D))$adj.r.squared
# ---- nested sequence ---------------------------------------------------------
r_l   <-r2(y~LAI)
r_lf  <-r2(y~LAI+fCover)
r_lfh <-r2(y~LAI+fCover+Hmax)
r_full<-r2(y~LAI+fCover+Hmax+profile)
r_str <-r2(y~Hmax+profile)                      # structure alone, reference only
cat(sprintf("nested: LAI %.3f | +fCover %.3f | +Hmax %.3f | +VCI(full) %.3f | structure-only %.3f | adjR2 full %.3f\n",
            r_l,r_lf,r_lfh,r_full,r_str,adj(y~LAI+fCover+Hmax+profile)))
#' Format the R-squared increment between two nested models
#' @param a R-squared of the smaller model
#' @param b R-squared of the larger model
#' @return character, signed to two decimals
inc<-function(a,b) sprintf("%+.2f",b-a)
M<-data.table(
  model=c("LAI","LAI + fCover","LAI + fCover<br>+ *H*<sub>max</sub>","all four<br>(+ profile)","structure only<br>(*H*<sub>max</sub> + profile)"),
  kind =c("nested","nested","nested","nested","reference"),
  R2   =c(r_l,r_lf,r_lfh,r_full,r_str))
M[,model:=factor(model,levels=model)]
PAL<-c(nested="#1A9850",reference="#999999")
# increment labels between consecutive nested bars
seg<-data.table(x=1:3+0.5, lab=c(inc(r_l,r_lf),inc(r_lf,r_lfh),inc(r_lfh,r_full)))
p<-ggplot(M,aes(model,R2,fill=kind))+geom_col(width=0.68,alpha=0.9)+
  geom_text(aes(label=sprintf("%.2f",R2)),vjust=-0.6,size=4,fontface=2)+
  geom_text(data=seg,aes(x=x,y=max(M$R2)+0.10,label=paste0("Δ",lab)),inherit.aes=FALSE,
            size=3.4,colour="grey25")+
  scale_fill_manual(values=PAL,guide="none")+coord_cartesian(ylim=c(0,1.05))+
  labs(x=NULL,y=expression("observed-buffering "*italic(R)^2*" (n = 53)"))+
  theme_article(12)+theme(axis.text.x=ggtext::element_markdown())
ggsave_article("out_files/Chapter1/figures/Fig6_obs_nested_native",p,7.4,4.6)
cat("DONE\n")
