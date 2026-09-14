# ==============================================================================
# Comparaison STATISTIQUE détaillée v3.2.0 vs v3.2.3 iter (vs HOBO) sur les
# températures et métriques associées : T horaire 1 m, ΔTmax & slope (tous / 10% chauds).
# Par plot (53). Sorties : tableau (r/R²/RMSE/MAE/biais + IC bootstrap),
# tests appariés (Wilcoxon sur |erreur| par plot ; bootstrap sur Δr), 2 figures.
#   Rscript c1_cmp_stats.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R","\\.R$",full.names=TRUE); src<-src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source)); source("pipeline/00_config.R")
set.seed(1); NB <- 3000; Z<-1; SH<-2L; HOTQ<-0.90; ds<-CFG$date_seq
dm <- as.data.table(extract_macro_daily(CFG$forcing_file, ds)); setnames(dm,"Tmax_macro","Tmac")
hot <- dm[Tmac>=quantile(Tmac,HOTQ,na.rm=TRUE),date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))

# --- HOBO horaire par plot ---
hb <- as.data.table(read.csv(CFG$hobo_temp_csv))[position_sensor=="a" & !(id_plot %in% CFG$ids_to_remove)]
hb[, datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]; hb<-hb[as.Date(datetime)%in%ds]
hb[, time:=floor_date(datetime,"hour")]; hb<-hb[,.(t_obs=mean(t_hobo,na.rm=TRUE)),by=.(id_plot,time)]; hb[,date:=as.Date(time)]
ids <- sort(unique(hb$id_plot))
# obs métriques par plot
obs <- hb[, { d<-.SD[,.(Tx=max(t_obs)),by=date]; d<-merge(d,dm,by="date")
  .(dTmax_all=mean(d$Tx-d$Tmac,na.rm=TRUE), dTmax_hot=mean(d[date%in%hot,Tx-Tmac],na.rm=TRUE),
    slope_all=coef(lm(t_obs~Tair_era5, merge(.SD,era5,by="time")))[2],
    slope_hot={m<-merge(.SD[date%in%hot],era5,by="time"); if(nrow(m)>10) coef(lm(t_obs~Tair_era5,m))[2] else NA_real_}) }, by=id_plot]

# --- sim par plot/version : Tc horaire 1 m, métriques + erreur horaire ---
extract <- function(nc_file) {
  if(!file.exists(nc_file)||file.size(nc_file)<1e5) return(NULL)
  nc<-tryCatch(nc_open(nc_file),error=function(e)NULL); if(is.null(nc))return(NULL); on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var)))return(NULL)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
  vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh
  if(Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15; tt<-floor_date(t0+dhours(th),"hour")
  data.table(time=tt, Tc=Tc)
}
VERS <- c("v3.2.0"="out_files/musica_hobo_z05/1111","v3.2.3 iter"="out_files/musica_hobo_z05_iter/1111")
perplot <- rbindlist(lapply(names(VERS),function(v) rbindlist(lapply(ids,function(id){
  s<-extract(file.path(VERS[[v]],sprintf("musica_out_HOBO_%s.nc",id))); if(is.null(s))return(NULL)
  # ΔTmax (-2h, daily max)
  sd<-copy(s); sd[,date:=as.Date(time-lubridate::hours(SH))]; dd<-sd[date%in%ds,.(Tx=max(Tc)),by=date]; dd<-merge(dd,dm,by="date")
  # slope (no shift) + horaire vs obs
  m<-merge(s, era5, by="time"); mo<-merge(s, hb[id_plot==id,.(time,t_obs)], by="time")
  data.table(id_plot=id, version=v,
    dTmax_all=mean(dd$Tx-dd$Tmac,na.rm=TRUE), dTmax_hot=mean(dd[date%in%hot,Tx-Tmac],na.rm=TRUE),
    slope_all=coef(lm(Tc~Tair_era5,m))[2], slope_hot={mh<-m[as.Date(time)%in%hot]; if(nrow(mh)>10)coef(lm(Tc~Tair_era5,mh))[2] else NA_real_},
    Th_rmse=sqrt(mean((mo$Tc-mo$t_obs)^2,na.rm=TRUE)), Th_r=cor(mo$Tc,mo$t_obs,use="complete.obs"), Th_n=nrow(mo))
}))))

# --- assemblage sim+obs par métrique (between-plot) ---
MM <- c("dTmax_all","dTmax_hot","slope_all","slope_hot")
D <- melt(perplot, id.vars=c("id_plot","version"), measure.vars=MM, variable.name="metric", value.name="sim")
O <- melt(obs, id.vars="id_plot", measure.vars=MM, variable.name="metric", value.name="obs")
D <- merge(D,O,by=c("id_plot","metric")); D<-D[is.finite(sim)&is.finite(obs)]; D[,err:=sim-obs]

# --- stats + IC bootstrap + tests ---
bootci <- function(x,y,f,nb=NB){ n<-length(x); v<-replicate(nb,{i<-sample(n,n,TRUE);f(x[i],y[i])}); quantile(v,c(.025,.975),na.rm=TRUE) }
ST <- D[, { ci<-bootci(obs,sim,function(a,b)cor(a,b))
  .(r=cor(obs,sim), R2=cor(obs,sim)^2, RMSE=sqrt(mean(err^2)), MAE=mean(abs(err)), bias=mean(err),
    r_lo=ci[1], r_hi=ci[2], n=.N) }, by=.(metric,version)]
# tests appariés v3.2.0 vs v3.2.3 (mêmes plots)
tests <- rbindlist(lapply(MM, function(mt){
  w<-dcast(D[metric==mt],id_plot~version,value.var="err"); w<-w[complete.cases(w)]
  a<-abs(w[["v3.2.0"]]); b<-abs(w[["v3.2.3 iter"]])
  wp<-wilcox.test(a,b,paired=TRUE)$p.value
  # bootstrap Δr (v3.2.0 − v3.2.3) sur plots appariés
  s<-dcast(D[metric==mt],id_plot~version,value.var="sim"); o<-D[metric==mt & version=="v3.2.0",.(id_plot,obs)]
  m<-merge(s,o,by="id_plot"); m<-m[complete.cases(m)]; n<-nrow(m)
  dr<-replicate(NB,{i<-sample(n,n,TRUE); cor(m$obs[i],m[["v3.2.0"]][i])-cor(m$obs[i],m[["v3.2.3 iter"]][i])})
  data.table(metric=mt, wilcox_absErr_p=wp, dr=cor(m$obs,m[["v3.2.0"]])-cor(m$obs,m[["v3.2.3 iter"]]),
             dr_lo=quantile(dr,.025), dr_hi=quantile(dr,.975), dr_p=mean(dr<=0)) }))
# horaire : test apparié sur RMSE par plot
wh <- dcast(perplot, id_plot~version, value.var="Th_rmse"); wh<-wh[complete.cases(wh)]
hourly_wilcox_p <- wilcox.test(wh[["v3.2.0"]], wh[["v3.2.3 iter"]], paired=TRUE)$p.value
hourly_tab <- perplot[, .(Th_RMSE_mean=mean(Th_rmse,na.rm=TRUE), Th_r_mean=mean(Th_r,na.rm=TRUE)), by=version]

fwrite(ST, "out_files/Chapter3/tables/tab_stats_v320_v323.csv")
fwrite(tests, "out_files/Chapter3/tables/tab_stats_tests_v320_v323.csv")
cat("=== skill par métrique×version (r, R², RMSE, MAE, biais, IC r) ===\n"); print(ST[order(metric,version)])
cat("\n=== T horaire (1 m) ===\n"); print(hourly_tab); cat(sprintf("Wilcoxon apparié RMSE horaire v3.2.0<v3.2.3 : p=%.2g\n", hourly_wilcox_p))
cat("\n=== tests appariés (v3.2.0 vs v3.2.3) ===\n"); print(tests)

# ---- Figure 1 : forest r ± IC bootstrap 95% ----
MPLAB <- c(dTmax_all="ΔTmax · tous", dTmax_hot="ΔTmax · chauds", slope_all="slope · tous", slope_hot="slope · chauds")
ST[, metric:=factor(metric, levels=names(MPLAB), labels=MPLAB)]
f1 <- ggplot(ST, aes(r, metric, colour=version)) +
  geom_errorbarh(aes(xmin=r_lo,xmax=r_hi), height=.25, position=position_dodge(.5)) +
  geom_point(size=2.6, position=position_dodge(.5)) +
  scale_colour_manual(values=c("v3.2.0"="#1B7837","v3.2.3 iter"="#9970AB"), name=NULL) +
  scale_x_continuous(limits=c(0,1)) +
  labs(x="r (corrélation entre-plots sim vs obs) ± IC bootstrap 95%", y=NULL,
       title="Skill v3.2.0 vs v3.2.3 iter — r ± IC 95% (bootstrap, 53 plots)",
       subtitle="IC disjoints = différence significative. Δr (v3.2.0−v3.2.3) > 0 sur les 4 métriques.") +
  theme_bw(base_size=12)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.title=element_text(size=11,face="bold"),plot.subtitle=element_text(size=8.5,colour="grey35"))
ggsave("out_files/Chapter3/figures/FigCmp_stats_skill_v320_v323iter.png", f1, width=8.5,height=4.5,dpi=200,bg="white")

# ---- Figure 2 : |erreur| par plot appariée + Wilcoxon p ----
D2 <- copy(D); D2[, metric:=factor(metric, levels=names(MPLAB), labels=MPLAB)]
medd <- dcast(D[, .(m=median(abs(err))), by=.(metric,version)], metric~version, value.var="m")
medd[, lower := ifelse(`v3.2.0` < `v3.2.3 iter`, "v3.2.0", "v3.2.3")]
pv <- merge(tests[,.(metric, wilcox_absErr_p)], medd[,.(metric=as.character(metric),lower)], by="metric")
pv[, metric := factor(MPLAB[metric], levels=MPLAB)]
pv[, lab := sprintf("Wilcoxon apparié p=%.2g\n(plus bas : %s)", wilcox_absErr_p, lower)]
f2 <- ggplot(D2, aes(version, abs(err), fill=version)) +
  geom_boxplot(width=.55, outlier.size=.6, alpha=.85) +
  geom_text(data=pv, aes(x=1.5, y=Inf, label=lab), vjust=1.1, size=2.5, inherit.aes=FALSE, colour="grey20") +
  facet_wrap(~metric, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("v3.2.0"="#1B7837","v3.2.3 iter"="#9970AB"), guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0.02,0.28))) +
  labs(x=NULL, y="|erreur| par plot = |sim − obs|",
       title="Erreur absolue (NIVEAU) par plot, appariée v3.2.0 vs v3.2.3 iter (53 plots)",
       subtitle="L'erreur absolue favorise plutôt v3.2.3 (biais chaud moindre) — significatif pour slope & ΔTmax-chauds. MAIS c'est le NIVEAU :\nle RANKING entre-plots (r, fig. skill) favorise nettement v3.2.0. Les deux figures = structure (v3.2.0) vs niveau (v3.2.3).") +
  theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face="bold",size=9),
    axis.text.x=element_text(angle=20,hjust=1),plot.title=element_text(size=11,face="bold"),
    plot.subtitle=element_text(size=8.5,colour="grey35"))
ggsave("out_files/Chapter3/figures/FigCmp_stats_error_v320_v323iter.png", f2, width=10,height=4.8,dpi=200,bg="white")
for (f in c("FigCmp_stats_skill_v320_v323iter.png","FigCmp_stats_error_v320_v323iter.png"))
  file.copy(file.path("out_files/Chapter3/figures",f), "Chapitre1/comparaison_versions/figures/", overwrite=TRUE)
cat("DONE -> tab_stats* + FigCmp_stats_{skill,error}_v320_v323iter.png\n")
