# ==============================================================================
# Extra comparison axes between MuSICA v3.2.0 and v3.2.3-iter (STATIC_ALS), to
# round out the advantages/disadvantages ledger. JJAS 2021, 53 plots, vs HOBO.
#   1. DAY vs NIGHT buffering slope   (ABL coupling acts on daytime convection)
#   2. Diurnal temperature range DTR = Tmax - Tmin   (direct coupling signature)
#   3. dTmin on cold nights (Tmin - Tmin_macro)      (the untested tail)
#   4. within-plot temporal tracking r(sim,obs) hourly   (operational fidelity)
#   5. mean bias (sim - obs)                          (calibration)
# Reads the existing nc; writes a tidy ledger + a 4-panel figure.
#   Rscript compare_binaries_extra_axes.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(ggplot2)
  library(dplyr); library(patchwork); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
YEAR<-2021L; MM<-6:9; DAY<-10:15; NIGHT<-c(22,23,0,1,2,3)
OUT_TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
OUT_FIG<-"/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
NCDIRS<-c("v3.2.0"=file.path(CFG_C3$out_dir,"nc","STATIC_ALS"),
          "v3.2.3 iter"=file.path(CFG_C3$out_dir,"nc_v323iter","STATIC_ALS"))
slp<-function(y,x){k<-is.finite(x)&is.finite(y); if(sum(k)<24)return(NA_real_); coef(lm(y[k]~x[k]))[2]}

# macro hourly
ncf<-nc_open(CFG_C3$forcing_file)
mh<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),T_macro=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
mh<-mh[year(time)==YEAR & month(time)%in%MM][,`:=`(hr=floor_date(time,"hour"),H=hour(time),d=as.Date(time))]
mday<-mh[,.(Tmac_max=max(T_macro),Tmac_min=min(T_macro)),by=d]

prep<-load_lai_prep(CFG_C3); df<-as.data.table(prep$df_plots)

# observed hourly
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&year(time)==YEAR&month(time)%in%MM,
       .(id_plot=as.character(id_plot),T=t_hobo,hr=floor_date(time,"hour"))]
hb<-merge(hb,mh[,.(hr,T_macro,H,d)],by="hr")

metrics_from<-function(r){ # r: data.table(T, T_macro, H, d)
  rd<-r[,.(Tmax=max(T),Tmin=min(T)),by=d]; rd<-merge(rd,mday,by="d")
  list(slope_day=slp(r[H%in%DAY]$T, r[H%in%DAY]$T_macro),
       slope_night=slp(r[H%in%NIGHT]$T, r[H%in%NIGHT]$T_macro),
       DTR=mean(rd$Tmax-rd$Tmin,na.rm=TRUE),
       dTmin=mean(rd$Tmin-rd$Tmac_min,na.rm=TRUE)) }

P<-rbindlist(lapply(df$id_plot,function(id){
  ro<-hb[id_plot==id,.(T,T_macro,H,d)]; mo<-metrics_from(ro)
  out<-data.table(id_plot=id, src="Observed", slope_day=mo$slope_day, slope_night=mo$slope_night,
                  DTR=mo$DTR, dTmin=mo$dTmin, temp_r=NA_real_, bias=NA_real_)
  for(v in names(NCDIRS)){ f<-file.path(NCDIRS[v],sprintf("musica_out_HOBO_%s.nc",id)); if(!file.exists(f))next
    s<-as.data.table(get_tair_at_z(nc_open(f),1.0));s<-s[year(time)==YEAR&month(time)%in%MM][,hr:=floor_date(time,"hour")]
    s<-merge(s,mh[,.(hr,T_macro,H,d)],by="hr"); rr<-s[,.(T=Tair_sim,T_macro,H,d)]; m<-metrics_from(rr)
    so<-merge(s[,.(hr,Tsim=Tair_sim)], hb[id_plot==id,.(hr,Tobs=T)], by="hr")
    out<-rbind(out, data.table(id_plot=id, src=v, slope_day=m$slope_day, slope_night=m$slope_night,
               DTR=m$DTR, dTmin=m$dTmin, temp_r=cor(so$Tsim,so$Tobs,use="complete.obs"),
               bias=mean(so$Tsim-so$Tobs,na.rm=TRUE))) }
  out }))
P[,src:=factor(src,levels=c("Observed","v3.2.0","v3.2.3 iter"))]
fwrite(P, file.path(OUT_TAB,"compare_binaries_extra_axes.csv"))

# ---- ledger: per axis, sim-vs-obs ranking(cor)/bias/amplitude ----------------
O<-dcast(P, id_plot~src, value.var=c("slope_day","slope_night","DTR","dTmin"))
cat("=== LEDGER (JJAS, n=53) — sim vs observed ===\n")
ax<-function(name,oc,c0,c3,unit=""){ o<-O[[oc]]
  cat(sprintf("%-14s | v3.2.0: r=%.2f rho=%.2f SDrec=%3.0f%% bias=%+.2f | v3.2.3: r=%.2f rho=%.2f SDrec=%3.0f%% bias=%+.2f %s\n",
    name, cor(o,O[[c0]],use="c"),cor(o,O[[c0]],method="spearman",use="c"),100*sd(O[[c0]],na.rm=T)/sd(o,na.rm=T),mean(O[[c0]]-o,na.rm=T),
          cor(o,O[[c3]],use="c"),cor(o,O[[c3]],method="spearman",use="c"),100*sd(O[[c3]],na.rm=T)/sd(o,na.rm=T),mean(O[[c3]]-o,na.rm=T),unit)) }
ax("slope DAY",   "slope_day_Observed","slope_day_v3.2.0","slope_day_v3.2.3 iter")
ax("slope NIGHT", "slope_night_Observed","slope_night_v3.2.0","slope_night_v3.2.3 iter")
ax("DTR (degC)",  "DTR_Observed","DTR_v3.2.0","DTR_v3.2.3 iter","degC")
ax("dTmin (degC)","dTmin_Observed","dTmin_v3.2.0","dTmin_v3.2.3 iter","degC")
cat("\n--- within-plot temporal tracking & bias (mean over 53 plots) ---\n")
print(P[src!="Observed",.(temporal_r=round(mean(temp_r),3), temporal_r_sd=round(sd(temp_r),3),
                          bias_degC=round(mean(bias),3)), by=src])

# ---- figure: 4 panels --------------------------------------------------------
PAL<-c("Observed"="grey30","v3.2.0"="#0072B2","v3.2.3 iter"="#D55E00")
# (a) day vs night slope
SL<-melt(P,id.vars=c("id_plot","src"),measure.vars=c("slope_day","slope_night"),variable.name="period",value.name="slope")
SL[,period:=factor(period,labels=c("Day 10-15h","Night 22-03h"))]
pa<-ggplot(SL,aes(src,slope,colour=src,fill=src))+geom_hline(yintercept=1,linetype="dotted",colour="grey55")+
  geom_violin(alpha=.12,colour=NA)+geom_boxplot(width=.25,alpha=0,outlier.shape=NA,linewidth=.5)+geom_jitter(width=.08,size=.7,alpha=.45)+
  facet_wrap(~period)+scale_colour_manual(values=PAL,guide="none")+scale_fill_manual(values=PAL,guide="none")+
  labs(x=NULL,y="Buffering slope",subtitle="(a) Day vs night: where ABL coupling acts")+
  theme_article(9)+theme(axis.text.x=element_text(angle=12,hjust=1))
# (b) DTR 1:1
DT<-dcast(P,id_plot~src,value.var="DTR"); limd<-range(c(DT$Observed,DT$`v3.2.0`,DT$`v3.2.3 iter`),na.rm=T)
DTl<-rbind(data.table(model="v3.2.0",obs=DT$Observed,pred=DT$`v3.2.0`),data.table(model="v3.2.3 iter",obs=DT$Observed,pred=DT$`v3.2.3 iter`))
pb<-ggplot(DTl,aes(pred,obs,colour=model))+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+
  geom_point(size=1.3,alpha=.8)+scale_colour_manual(values=PAL[2:3],name=NULL)+coord_fixed(xlim=limd,ylim=limd)+
  labs(x="Predicted DTR (degC)",y="Observed DTR",subtitle="(b) Diurnal range Tmax-Tmin")+theme_article(9)+theme(legend.position="top")
# (c) dTmin ranking (rank-rank)
DM<-dcast(P,id_plot~src,value.var="dTmin")
DMl<-rbind(data.table(model="v3.2.0",ro=rank(DM$Observed),rp=rank(DM$`v3.2.0`),rho=cor(DM$Observed,DM$`v3.2.0`,method="spearman")),
           data.table(model="v3.2.3 iter",ro=rank(DM$Observed),rp=rank(DM$`v3.2.3 iter`),rho=cor(DM$Observed,DM$`v3.2.3 iter`,method="spearman")))
pc<-ggplot(DMl,aes(ro,rp,colour=model))+geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey60")+geom_point(size=1.3,alpha=.8)+
  scale_colour_manual(values=PAL[2:3],name=NULL)+coord_fixed()+labs(x="Rank obs",y="Rank pred",subtitle="(c) Cold-night dTmin ranking")+
  theme_article(9)+theme(legend.position="top")
# (d) temporal r + bias
pd<-ggplot(P[src!="Observed"],aes(src,temp_r,colour=src))+geom_boxplot(width=.4,outlier.shape=NA)+geom_jitter(width=.1,size=.8,alpha=.5)+
  scale_colour_manual(values=PAL[2:3],guide="none")+labs(x=NULL,y="Within-plot hourly r(sim,obs)",subtitle="(d) Temporal tracking fidelity")+
  theme_article(9)
g<-(pa|pb)/(pc|pd)+plot_annotation(title="Extra comparison axes: MuSICA v3.2.0 vs v3.2.3-iter (JJAS, n=53)",
  theme=theme(plot.title=element_text(face="bold")))
ggsave_article(file.path(OUT_FIG,"Fig_compare_binaries_extra_axes"),g,11,8.4)
cat("\nDONE -> Fig_compare_binaries_extra_axes + compare_binaries_extra_axes.csv\n")
