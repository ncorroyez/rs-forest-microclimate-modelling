# ==============================================================================
# Supplementary hot-day companions S3-S5, from the hot-day extraction of stage A4:
#   S3  glass-ceiling dumbbell, observed vs full-canopy simulation
#   S4  forward inclusion, cumulative LiDAR variables
#   S5  model-free nested R2 on hot-day observed buffering
# Pure re-plotting; no MuSICA run.
#
# Reads : out_files/Chapter1/hot_extract.rds                      (stage A4)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv  (stage A2: P labels)
#         in_files/data_Blois_utm31n.geojson                      (logger footprints)
#         in_files_native20/{lai_z1,fCover,max}_res_10_m.tif      (traits)
#         in_files/lad_z05/Blois_lad_z05_r25.csv                  (LAD -> profile shape)
# Writes: outputs/figures_chap1/FigS_obs_vs_sim_hot.{png,pdf}
#         outputs/figures_chap1/FigS_forward_hot.{png,pdf}
#         outputs/figures_chap1/FigS_obs_nested_hot.{png,pdf}
#   Rscript c1_hot_figures.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggtext);source("scripts/_article_style.R")})
H<-readRDS("out_files/Chapter1/hot_extract.rds"); OBS<-H$obs; CO<-H$coal
clu<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
clab<-clu[,.(id_plot,P=factor(P,levels=c("P1","P2","P3","P4")))]
steps<-c("+LAI"="1000","+H<sub>max</sub>"="1100","+fCover"="1110","+LAD"="1111")
# ---- Fig 4 hot: glass ceiling dumbbell (full canopy 1111) ----
V<-merge(OBS,CO[bit=="1111",.(id_plot=id,sim_dt=dt_hot,sim_sl=sl_hot)],by="id_plot")
V<-merge(V,clab,by="id_plot")[is.finite(obs_dt_hot)&is.finite(sim_dt)&is.finite(obs_sl_hot)&is.finite(sim_sl)]
setorder(V,P,obs_dt_hot); V[,rank:=.I]
#' Hot-day dumbbell panel pairing an observed and a simulated column
#' @param oc name of the observed column in V
#' @param sc name of the simulated column in V
#' @param ref x value of the dashed no-effect reference line
#' @param xl x axis label
#' @return a ggplot object
dum<-function(oc,sc,ref,xl){L<-melt(V,id.vars=c("rank","P"),measure.vars=c(oc,sc),variable.name="s",value.name="v");L[,s:=factor(s,levels=c(oc,sc),labels=c("observed","simulated"))]
  ggplot()+geom_vline(xintercept=ref,linetype="dashed",colour="grey55")+geom_segment(data=V,aes(x=get(oc),xend=get(sc),y=rank,yend=rank,colour=P),linewidth=0.5,alpha=0.55)+
  geom_point(data=L,aes(v,rank,colour=P,shape=s),size=2)+scale_colour_manual(values=PAL_CLUSTER,name=NULL)+scale_shape_manual(values=c(observed=16,simulated=17),name=NULL)+
  guides(colour=guide_legend(order=1,override.aes=list(shape=15,size=3.5)),shape=guide_legend(order=2,override.aes=list(colour="grey30")))+labs(x=xl,y=NULL)+theme_article(11)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())}
f4<-(dum("obs_dt_hot","sim_dt",0,expression(Delta*T[max]~"(°C)"))|dum("obs_sl_hot","sim_sl",1,"micro-macro buffering slope"))+plot_layout(guides="collect")+plot_annotation(tag_levels="a")&theme(legend.position="bottom",plot.tag=element_text(face="bold",size=14))
ggsave_article("outputs/figures_chap1/FigS_obs_vs_sim_hot",f4,9,6.5)
# ---- Fig 5 hot: forward ----
#' Cross-plot r at each forward step, per archetype and pooled
#' @param valcol name of the simulated hot-day column in the coalition table
#' @param obscol name of the matching observed hot-day column
#' @return data.table(rowlab, step, r) with rowlab in P1..P4 and "All"
fwd<-function(valcol,obscol){D<-merge(CO[,.(id_plot=id,bit,val=get(valcol))],OBS[,.(id_plot,obs=get(obscol))],by="id_plot");D<-merge(D,clab,by="id_plot")
  rc<-function(i2)sapply(steps,function(b){m<-D[bit==b&id_plot%in%i2];if(sd(m$val,na.rm=T)>0)cor(m$val,m$obs,use="complete.obs") else NA})
  rbindlist(c(lapply(c("P1","P2","P3","P4"),function(cl)data.table(rowlab=cl,step=names(steps),r=rc(clab[P==cl,id_plot]))),list(data.table(rowlab="All",step=names(steps),r=rc(clab$id_plot)))))}
#' One forward-inclusion panel from a table produced by fwd()
#' @param R data.table(rowlab, step, r)
#' @param yl y axis label
#' @return a ggplot object
mkf<-function(R,yl){R[,step:=factor(step,levels=names(steps))];R[,rowlab:=factor(rowlab,levels=c("All","P1","P2","P3","P4"))]
  ggplot(R,aes(step,r,colour=rowlab,group=rowlab))+geom_line(aes(linewidth=rowlab=="All"))+geom_point(aes(size=rowlab=="All"))+scale_colour_manual(values=c(All="grey15",PAL_CLUSTER),name=NULL)+scale_linewidth_manual(values=c(0.7,1.6),guide="none")+scale_size_manual(values=c(2,3.4),guide="none")+labs(x="LiDAR variable added (cumulative)",y=yl)+theme_article(11)+theme(axis.text.x=ggtext::element_markdown())}
f5<-(mkf(fwd("dt_hot","obs_dt_hot"),expression("cross-plot "*italic(r)*" ("*Delta*T[max]*")"))|mkf(fwd("sl_hot","obs_sl_hot"),expression("cross-plot "*italic(r)*" (slope)")))+plot_layout(guides="collect")+plot_annotation(tag_levels="a")&theme(legend.position="bottom",plot.tag=element_text(face="bold",size=14))
ggsave_article("outputs/figures_chap1/FigS_forward_hot",f5,9.5,4.8)
# ---- Fig 8 hot: obs-only NESTED R2 --------------------------------------------
# Mirrors scripts/c1_fig6_obs_nested_native.R exactly (native20 traits, nested
# progression, vertical-profile shape as the fourth term) but on hot-day
# observed buffering. VCI is deliberately NOT used: it is a density descriptor,
# not a measure of profile shape, so a VCI arm would not test the vertical
# dimension at all (committee 2026-07-28). The previous version of this block
# still used VCI and the superseded four-arm layout, so it did not match the
# Figure 8 it claims to reproduce.
suppressPackageStartupMessages({library(terra);library(sf)})
g<-st_read("in_files/data_Blois_utm31n.geojson",quiet=TRUE)
r0<-rast("in_files_native20/lai_z1_res_10_m.tif");g<-st_transform(g,crs(r0));v<-vect(g)
D6<-data.table(id_plot=g$id_plot,
  LAI   =terra::extract(rast("in_files_native20/lai_z1_res_10_m.tif"),v)[,2],
  fCover=terra::extract(rast("in_files_native20/fCover_res_10_m.tif"),v)[,2],
  Hmax  =terra::extract(rast("in_files_native20/max_res_10_m.tif"),v)[,2])
lad<-fread("in_files/lad_z05/Blois_lad_z05_r25.csv")
lyr<-grep("^LAD_Layer_",names(lad),value=TRUE); hh<-as.numeric(sub("LAD_Layer_","",lyr))
Lm<-as.matrix(lad[,..lyr]); Lm[!is.finite(Lm)]<-0
lad[,cen:=(Lm%*%hh)[,1]/pmax(rowSums(Lm),1e-9)]
lad[,profile:=cen/pmax(Hmax,1e-9)]                    # 0-1, higher = top-heavy
D6<-merge(merge(D6,lad[,.(id_plot,profile)],by="id_plot"),OBS[,.(id_plot,obs_dt_hot)],by="id_plot")
D6<-D6[is.finite(obs_dt_hot)&is.finite(LAI)&is.finite(fCover)&is.finite(Hmax)&is.finite(profile)]
D6[,y:=-obs_dt_hot]
#' R-squared of a linear model on the hot-day observation table D6
#' @param f model formula
#' @return numeric R-squared
r2<-function(f)summary(lm(f,D6))$r.squared
r_l<-r2(y~LAI); r_lf<-r2(y~LAI+fCover); r_lfh<-r2(y~LAI+fCover+Hmax)
r_full<-r2(y~LAI+fCover+Hmax+profile); r_str<-r2(y~Hmax+profile)
M<-data.table(
  model=c("LAI","LAI + fCover","LAI + fCover<br>+ *H*<sub>max</sub>","all four<br>(+ profile)","structure only<br>(*H*<sub>max</sub> + profile)"),
  kind =c("nested","nested","nested","nested","reference"),
  R2   =c(r_l,r_lf,r_lfh,r_full,r_str))
M[,model:=factor(model,levels=model)]
#' Format the R-squared increment between two nested models
#' @param a R-squared of the smaller model
#' @param b R-squared of the larger model
#' @return character, signed to two decimals
inc<-function(a,b) sprintf("%+.2f",b-a)
seg<-data.table(x=1:3+0.5, lab=c(inc(r_l,r_lf),inc(r_lf,r_lfh),inc(r_lfh,r_full)))
f6<-ggplot(M,aes(model,R2,fill=kind))+geom_col(width=0.68,alpha=0.9)+
  geom_text(aes(label=sprintf("%.2f",R2)),vjust=-0.6,size=4,fontface=2)+
  geom_text(data=seg,aes(x=x,y=max(M$R2)+0.10,label=paste0("Δ",lab)),inherit.aes=FALSE,size=3.4,colour="grey25")+
  scale_fill_manual(values=c(nested="#1A9850",reference="#999999"),guide="none")+
  coord_cartesian(ylim=c(0,1.05))+
  labs(x=NULL,y=bquote("observed hot-day buffering "*italic(R)^2*" (n = "*.(nrow(D6))*")"))+
  theme_article(12)+theme(axis.text.x=ggtext::element_markdown())
ggsave_article("outputs/figures_chap1/FigS_obs_nested_hot",f6,7.4,4.6)
cat(sprintf("S3hot: dt r=%.2f sl r=%.2f | S5hot nested: LAI %.2f | +fCover %.2f | +Hmax %.2f | +profile %.2f | structure-only %.2f\nDONE\n",
  cor(V$sim_dt,V$obs_dt_hot),cor(V$sim_sl,V$obs_sl_hot),r_l,r_lf,r_lfh,r_full,r_str))
