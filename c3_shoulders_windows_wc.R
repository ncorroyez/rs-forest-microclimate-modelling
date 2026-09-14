# ==============================================================================
# Chapter 3 §3.5 — does variable (S2-driven) LAI beat a fixed / parametric LAI at
# the phenological SHOULDERS? Per-plot ΔTmax vs HOBO over leaf-out / summer /
# autumn, genuine-53 v3.2.3. Contrast: CONST (LAI flat all year) vs STATIC
# (parametric seasonal curve) vs DYN (S2 phenology) vs FUSION_H (layered).
# Metrics per window: R², SD-recovery, bias, RMSE, n (HOBO deployment varies).
#   Rscript c3_shoulders_windows.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
d3<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); d1<-"out_files/Chapter3/nc_genuine53_windcorr"
# scenario -> (dir,label,type)
SC<-list(
  CONST_ALS       =list(d=d3,lab="CONST (LAI fixe an)",    typ="fixe"),
  STATIC_ALS      =list(d=d3,lab="STATIC (pheno param.)", typ="param"),
  DYN_ALS         =list(d=d3,lab="DYN LiDAR-mag x S2-timing",typ="dyn"),
  DYN_S2_ATBD     =list(d=d3,lab="DYN S2 (mag+timing)",   typ="dyn"),
  DYN_S2_ANNUAL   =list(d=d3,lab="DYN S2 annual",         typ="dyn"),
  FUSION_H        =list(d=d1,lab="FUSION par couche",     typ="fusion"))
wins<-list(leafout=c("2021-04-01","2021-05-31"), summer=c("2021-06-01","2021-09-30"), autumn=c("2021-10-01","2021-11-30"))

ncf<-nc_open(CFG_C3$forcing_file); MA<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA[,date:=as.Date(time)]; mday_all<-MA[,.(Tmax_macro=max(Tm,na.rm=TRUE)),by=date]
hb0<-as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb0[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb0<-hb0[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove),.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time))]
read_days<-function(dir,scn){ fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
  rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
    s<-as.data.table(r)[,.(date=as.Date(time),Tsim=Tair_sim)]; s[,id_plot:=id]; s[,.(Tmax_sim=max(Tsim,na.rm=TRUE)),by=.(id_plot,date)]}),fill=TRUE) }
SIMS<-lapply(names(SC),function(s) read_days(SC[[s]]$d,s)); names(SIMS)<-names(SC)

res<-list()
for(wn in names(wins)){ ds<-seq(as.Date(wins[[wn]][1]),as.Date(wins[[wn]][2]),by="day"); md<-mday_all[date%in%ds]
  ob<-merge(hb0[date%in%ds,.(Tmax_obs=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],md,by="date")[,.(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)),by=id_plot]
  ob<-ob[is.finite(do)]; sdo<-sd(ob$do)
  for(s in names(SC)){ sday<-SIMS[[s]][date%in%ds]
    pp<-merge(sday[md,on="date",nomatch=0L][,.(dsim=mean(Tmax_sim-Tmax_macro,na.rm=TRUE)),by=id_plot], ob, by="id_plot")
    pp<-pp[is.finite(dsim)&is.finite(do)]
    res[[length(res)+1]]<-data.table(window=wn, scenario=SC[[s]]$lab, typ=SC[[s]]$typ, n=nrow(pp),
      R2=round(cor(pp$dsim,pp$do)^2,3), SDrec=round(sd(pp$dsim)/sdo,2), bias=round(mean(pp$dsim-pp$do),2), RMSE=round(sqrt(mean((pp$dsim-pp$do)^2)),2)) } }
R<-rbindlist(res); R[,window:=factor(window,levels=names(wins))]
fwrite(R,"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table5_shoulders_windows_genuine.csv")
for(wn in names(wins)){ cat(sprintf("\n===== %s =====\n",toupper(wn))); print(R[window==wn,.(scenario,typ,n,R2,SDrec,bias,RMSE)]) }
cat("\nDONE -> Table5_shoulders_windows_genuine.csv\n")
