# ==============================================================================
# Chapter 3 — "beyond summer": validate ΔTmax & slope (1 m, -2h, all temps) on the
# 53 HOBO over THREE windows — spring leaf-out shoulder, summer, autumn senescence
# shoulder. Tests whether the DYNAMIC S2 forcing catches up to the STATIC LiDAR at
# the shoulders, where the canopy leaf area changes and a single summer value is wrong.
# Run from z_Example root:  Rscript extend_seasonal.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(data.table); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
source("/home/corroyez/Documents/NC_Full/scripts/_article_style.R")
ncp <- file.path(CFG_C3$out_dir, "nc")
scen <- c(DYN_S2_ATBD_NM="S2 ATBD (dyn)", DYN_S2_OPT_NM="S2 opt (dyn)",
          DYN_ATBD_rfull_NM="S2 ATBD ×ratio→full (dyn)", DYN_OPT_rdopt_NM="S2 opt ×ratio→d_opt (dyn)",
          STATIC_ALS="LiDAR full (static)", STATIC_ALS_DOPT="LiDAR d_opt (static)")
windows <- list(
  spring = c("2021-04-15","2021-05-31"),
  summer = c("2021-06-01","2021-09-30"),
  autumn = c("2021-10-01","2021-11-15"))

# hourly macro (whole year), shifted/used as in the summer extractor
ncf <- nc_open(CFG_C3$forcing_file)
macro_h <- data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"), Tmacro=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
macro_h[, date:=as.Date(time)]
hb0 <- as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb0[, time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb0 <- hb0[position_sensor=="a" & !(id_plot %in% CFG_C3$ids_to_remove), .(id_plot=as.character(id_plot), time, Tobs=t_hobo)]
hb0[, hr:=floor_date(time,"hour")]
slp <- function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
read_sim_1m <- function(f){ nc<-tryCatch(nc_open(f),error=function(e)NULL); if(is.null(nc))return(NULL)
  res<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL); nc_close(nc); if(is.null(res)||!nrow(res))return(NULL)
  r<-as.data.table(res); r[,time:=time-hours(2)]; r[,.(time,Tsim=Tair_sim)] }

one_window <- function(wn, rng){
  ds <- seq(as.Date(rng[1]), as.Date(rng[2]), by="day")
  dmd <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds))
  hd <- as.data.table(read_hobo_daily(CFG_C3$hobo_temp_csv, ds, dmd, CFG_C3$ids_to_remove)); hd[,id_plot:=as.character(id_plot)]
  mh <- macro_h[date %in% ds]; mhr <- mh[,.(hr=floor_date(time,"hour"),Tmacro)]
  hb <- merge(hb0[as.Date(time)%in%ds], mhr, by="hr")
  OBSdT <- hd[, .(do=mean(Delta_obs,na.rm=TRUE)), by=id_plot]
  OBSsl <- hb[, .(so=slp(Tobs,Tmacro)), by=id_plot]
  OBS <- merge(OBSdT,OBSsl,by="id_plot")
  rbindlist(lapply(names(scen), function(scn){
    dir<-file.path(ncp,scn); fs<-list.files(dir,pattern="\\.nc$",full.names=TRUE); if(!length(fs))return(NULL)
    dT <- as.data.table(extract_deltatmax_hobo_scenario(dir,dmd,ds))[, .(ds=mean(Delta_Tmax,na.rm=TRUE)), by=.(id_plot=as.character(id_plot))]
    sl <- rbindlist(lapply(fs,function(f){ id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)"); s<-read_sim_1m(f); if(is.null(s))return(NULL)
      s<-merge(s, mh[,.(time,Tmacro)], by="time"); data.table(id_plot=id, ss=slp(s$Tsim,s$Tmacro)) }),fill=TRUE)
    d <- merge(merge(dT,sl,by="id_plot"), OBS, by="id_plot")
    dd<-d[is.finite(ds)&is.finite(do)]; sd<-d[is.finite(ss)&is.finite(so)]
    data.table(window=wn, scenario=scen[[scn]], n=nrow(dd),
      dT_r2=cor(dd$ds,dd$do)^2, sl_r2=cor(sd$ss,sd$so)^2)
  }),fill=TRUE)
}
res <- rbindlist(lapply(names(windows), function(w) one_window(w, windows[[w]])), fill=TRUE)
res[, window:=factor(window, levels=c("spring","summer","autumn"),
       labels=c("spring (leaf-out)","summer","autumn (senescence)"))]
fwrite(res, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4i_seasonal_windows.csv")
cat("\n=== ΔTmax & slope R² by window ===\n"); print(res[order(window,-dT_r2)])
ml <- melt(res, id.vars=c("window","scenario"), measure.vars=c("dT_r2","sl_r2"), variable.name="metric", value.name="R2")
ml[, metric:=ifelse(metric=="dT_r2","ΔT_max","slope")]
ml[, typ:=fifelse(grepl("LiDAR",scenario),"LiDAR","S2")]
g <- ggplot(ml, aes(window, R2, group=scenario, colour=typ)) +
  geom_line(aes(linetype=scenario), linewidth=0.8) + geom_point(size=1.8) +
  facet_wrap(~metric) + scale_colour_manual(values=c(LiDAR="#1A9850", S2="#D7191C"), name=NULL) +
  scale_linetype_manual(values=c(1,2,3,4,1,2), guide=guide_legend(ncol=2)) +
  labs(x=NULL, y=expression(R^2~"(sim vs obs, 53 HOBO)"), linetype=NULL) +
  coord_cartesian(ylim=c(0,1)) + theme_article(11) + theme(legend.position="bottom", axis.text.x=element_text(angle=15,hjust=1))
outdir<-"/home/corroyez/Documents/NC_Full/outputs/figures_article_ch3"
ggsave_article(file.path(outdir,"Fig_seasonal_windows"), g, 9, 5)
cat("\nDONE seasonal windows\n")
