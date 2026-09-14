# ==============================================================================
# Fig 5 (glass-ceiling dumbbell) + Fig J1 (forward inclusion), regenerated under
# the canonical time-matched ΔTmax convention (R/dtmax_convention.R) and on the
# NATIVE20 lineage throughout. Fig 5 previously read the z05/station variant and
# Fig J1 read cached convention-A metrics; both are re-derived here from
# out_files/musica_native20_forward/ so the whole chapter shares one lineage.
# No new MuSICA runs.
# Reads :
#         out_files/musica_native20_forward/{0000,1000,1100,1110,1111}/*.nc
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv  (stage A1)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv      (stage A2)
# Writes: out_files/Chapter1/figures/Fig5_obs_vs_sim_dumbbell.{png,pdf}
#         out_files/Chapter1/figures/FigJ1_forward_inclusion.{png,pdf}
#         out_files/Chapter1/tables/tab_forward_inclusion_convB.csv
#   Rscript scripts/c1_fig5_j1_convB.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)
  library(ggplot2);library(patchwork);library(ggtext)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("scripts/_article_style.R"); source("R/cluster_relabel.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1
FORC<-"in_files/FR-Blo_2021_v2.nc"; MREF<-macro_ref(FORC,ds); MOBS<-macro_ref_obs(FORC,ds)
macH<-{nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time")
  t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
  d<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
  d[as.Date(time)%in%ds][,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]}
LV<-c("P1","P2","P3","P4")
clu<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot,P=factor(P,levels=LV))]
V<-fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")   # deja convention B

# ---------------- Fig 5 : dumbbell obs vs sim, par archetype ------------------
D5<-merge(V,clu,by="id_plot")[is.finite(obs_dt)&is.finite(sim_dt)]
#' Dumbbell panel pairing an observed and a simulated column, one row per logger
#' @param D per-logger table carrying both columns and P
#' @param oc name of the observed column
#' @param sc name of the simulated column
#' @param xl x axis label
#' @param ref x value of the dashed no-effect reference line
#' @param ordby column fixing the vertical order. Both panels are ordered by the
#'   observed DeltaTmax (committee 2026-08-09) so that a logger sits on the same
#'   row in (a) and (b), and the two orderings can be read against each other.
#' @return a ggplot object
mk_dumb<-function(D,oc,sc,xl,ref,ordby=oc){
  d<-copy(D); setorderv(d,ordby); d[,rank:=.I]
  L<-melt(d,id.vars=c("id_plot","rank","P"),measure.vars=c(oc,sc),variable.name="src",value.name="v")
  L[,src:=factor(src,levels=c(oc,sc),labels=c("observed","MuSICA (simulated)"))]
  ggplot()+geom_vline(xintercept=ref,linetype="dashed",colour="grey55")+
    geom_segment(data=d,aes(x=get(oc),xend=get(sc),y=rank,yend=rank,colour=P),linewidth=0.5,alpha=0.6)+
    geom_point(data=L,aes(v,rank,colour=P,shape=src),size=2.2)+
    scale_colour_manual(values=PAL_CLUSTER,name=NULL)+
    scale_shape_manual(values=c("observed"=16,"MuSICA (simulated)"=17),name=NULL)+
    labs(x=xl,y=expression("sub-canopy plots, ordered by observed "*Delta*T[max]))+theme_article(11)+
    theme(legend.position="bottom",axis.text.y=element_blank(),axis.ticks.y=element_blank())}
pa<-mk_dumb(D5,"obs_dt","sim_dt",expression(Delta*T[max]~"(°C)"),0)
pb<-mk_dumb(D5[is.finite(obs_sl)&is.finite(sim_sl)],"obs_sl","sim_sl","micro-macro buffering slope",1,
            ordby="obs_dt")
f5<-(pa|pb)+plot_layout(guides="collect")+plot_annotation(tag_levels="a")&
     theme(legend.position="bottom",plot.tag=element_text(face="bold",size=13))
ggsave_article("out_files/Chapter1/figures/Fig5_obs_vs_sim_dumbbell",f5,9.5,6.4)
cat(sprintf("Fig5: n=%d | dTmax r=%.3f biais=%+.2f | slope r=%.3f\n",nrow(D5),
  cor(D5$obs_dt,D5$sim_dt),mean(D5$sim_dt-D5$obs_dt),
  cor(D5$obs_sl,D5$sim_sl,use="complete.obs")))

# ---------------- Fig J1 : forward inclusion, coalitions native20 -------------
FWD<-"out_files/musica_native20_forward"
steps<-c("uniform"="0000","+LAI"="1000","+Hmax"="1100","+fCover"="1110","+LAD"="1111")
ids<-V$id_plot
#' DeltaTmax and buffering slope for one coalition at one logger
#' @param bit 4-bit coalition code, for instance "1110"
#' @param id logger id
#' @return named numeric c(dt, sl); both NA when the run is unreadable
ext<-function(bit,id){ m<-micro_hourly_at(file.path(FWD,bit,sprintf("musica_out_HOBO_%s.nc",id)),Z)
  if(is.null(m)) return(c(dt=NA_real_,sl=NA_real_))
  mm<-merge(m,macH,by="time")
  c(dt=delta_tmax_mean(m,MREF,ds),
    sl=if(nrow(mm)>50) as.numeric(coef(lm(Tmic~Tmac,mm))[2]) else NA_real_)}
COAL<-rbindlist(lapply(names(steps),function(nm) rbindlist(lapply(ids,function(id){
  v<-ext(steps[[nm]],id); data.table(step=nm,id_plot=id,sim_dt=v["dt"],sim_sl=v["sl"])}))))
COAL<-merge(COAL,V[,.(id_plot,obs_dt,obs_sl)],by="id_plot")
COAL<-merge(COAL,clu,by="id_plot")
#' Cross-plot agreement of a simulated column against its observed counterpart
#' @param D table of plots for one forward step and group
#' @param sv name of the simulated column
#' @param ov name of the observed column
#' @return named numeric c(r, rmse); NA when fewer than 4 pairs or zero simulated variance
sc<-function(D,sv,ov) if(sum(is.finite(D[[sv]])&is.finite(D[[ov]]))>3 && sd(D[[sv]],na.rm=TRUE)>0)
  c(r=cor(D[[sv]],D[[ov]],use="complete.obs"),
    rmse=sqrt(mean((D[[sv]]-D[[ov]])^2,na.rm=TRUE))) else c(r=NA,rmse=NA)
R<-rbindlist(lapply(names(steps),function(nm){
  d<-COAL[step==nm]
  rbind(data.table(step=nm,grp="All",metric="dTmax",t(sc(d,"sim_dt","obs_dt"))),
        data.table(step=nm,grp="All",metric="slope",t(sc(d,"sim_sl","obs_sl"))),
        rbindlist(lapply(LV,function(p){dd<-d[P==p]
          rbind(data.table(step=nm,grp=p,metric="dTmax",t(sc(dd,"sim_dt","obs_dt"))),
                data.table(step=nm,grp=p,metric="slope",t(sc(dd,"sim_sl","obs_sl"))))})))}))
R[,step:=factor(step,levels=names(steps))][,grp:=factor(grp,levels=c("All",LV))]
fwrite(R,"out_files/Chapter1/tables/tab_forward_inclusion_convB.csv")
#' One forward-inclusion panel: a metric across the cumulative variable additions
#' @param met metric to draw, "dTmax" or "slope"
#' @param yv column to plot on y, "r" or "rmse"
#' @param ylab y axis label
#' @return a ggplot object
mkp<-function(met,yv,ylab){ ggplot(R[metric==met&is.finite(get(yv))],aes(step,get(yv),colour=grp,group=grp))+
  geom_line(aes(linewidth=grp=="All"))+geom_point(aes(size=grp=="All"))+
  scale_colour_manual(values=c(All="grey15",PAL_CLUSTER),name=NULL)+
  scale_linewidth_manual(values=c(0.7,1.6),guide="none")+scale_size_manual(values=c(1.9,3.2),guide="none")+
  labs(x="LiDAR variable added (cumulative)",y=ylab)+theme_article(11)+
  theme(axis.text.x=element_markdown(angle=20,hjust=1))}
j1<-(mkp("dTmax","r",expression("cross-plot "*italic(r)*" ("*Delta*T[max]*")")) |
     mkp("slope","r",expression("cross-plot "*italic(r)*" (slope)")))/
    (mkp("dTmax","rmse",expression("RMSE ("*Delta*T[max]*", °C)")) |
     mkp("slope","rmse","RMSE (slope)"))+
    plot_layout(guides="collect")+plot_annotation(tag_levels="a")&
    theme(legend.position="bottom",plot.tag=element_text(face="bold",size=13))
ggsave_article("out_files/Chapter1/figures/FigJ1_forward_inclusion",j1,9.5,8.4)
cat("\n=== forward inclusion (All) ===\n"); print(R[grp=="All"][,.(step,metric,r=round(r,3),rmse=round(rmse,2))])
cat("DONE\n")
