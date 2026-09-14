# ==============================================================================
# Chapter 3 — dig, SLOPE version. Same as c3_dig_ladopt.R but the thermal-coupling
# slope (per-plot coef lm(Tsim~Tmacro), hourly, summer) as well as ΔTmax, each
# with R² AND SD-recovery (SD_sim/SD_obs), bias, RMSE. Tests whether the slope
# R² of DYN_ALS_LADOPT is genuine or another collapsed-amplitude artifact.
#   Rscript c3_dig_ladopt_slope.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
d3 <- file.path(CFG_C3$out_dir, "nc_genuine53")
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
ncf<-nc_open(CFG_C3$forcing_file); macro<-data.table(time=force_utc_nc(CFG_C3$forcing_file,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
macro<-macro[as.Date(time)%in%ds]; macro[,hr:=floor_date(time,"hour")]; mday<-macro[,.(Tmax_macro=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA_real_
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv)); hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,macro[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmax_obs=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmax_obs-Tmax_macro,na.rm=TRUE)),by=id_plot],
           hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn,dir=d3){ fs<-list.files(file.path(dir,scn),pattern="\\.nc$",full.names=TRUE)
  sim<-rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
    r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
    s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))]
    sm<-merge(s,macro[,.(hr,Tm)],by="hr"); sd<-merge(s[,.(Tmax_sim=max(Tsim,na.rm=TRUE)),by=date],mday,by="date")
    data.table(id_plot=id, ds=mean(sd$Tmax_sim-sd$Tmax_macro,na.rm=TRUE), ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)
  merge(sim,OBS,by="id_plot") }
R2<-function(x,y) cor(x,y,use="complete.obs")^2; sdo_d<-sd(OBS$do); sdo_s<-sd(OBS$so,na.rm=TRUE)
row<-function(s,d){ data.table(scenario=s,n=nrow(d),
  dT_R2=round(R2(d$ds,d$do),3), dT_SDrec=round(sd(d$ds)/sdo_d,2), dT_bias=round(mean(d$ds-d$do),2), dT_RMSE=round(sqrt(mean((d$ds-d$do)^2)),2),
  sl_R2=round(R2(d$ss,d$so),3), sl_SDrec=round(sd(d$ss,na.rm=TRUE)/sdo_s,2), sl_bias=round(mean(d$ss-d$so,na.rm=TRUE),3), sl_RMSE=round(sqrt(mean((d$ss-d$so)^2,na.rm=TRUE)),3)) }
cat(sprintf("obs SD: ΔTmax=%.2f °C | slope=%.3f\n\n", sdo_d, sdo_s))
cat("=== LADOPT effect on the SLOPE (full vs top-d_opt LAD) ===\n")
for(p in list(c("DYN_ALS","DYN_ALS_LADOPT"),c("STATIC_S2_ATBD","STATIC_S2_ATBD_LADOPT"),c("DYN_S2_DOPT","DYN_S2_DOPT_LADOPT"))){
  a<-row(p[1],getpp(p[1])); b<-row(p[2],getpp(p[2]))
  cat(sprintf("  %-22s sl_R2=%.3f sl_SDrec=%.2f  ->  %-24s sl_R2=%.3f sl_SDrec=%.2f\n",p[1],a$sl_R2,a$sl_SDrec,p[2],b$sl_R2,b$sl_SDrec)) }
cat("\n=== full table: DYN_ALS_LADOPT vs DYN_S2 family (both metrics) ===\n")
fam<-c("DYN_ALS_LADOPT","DYN_ALS","STATIC_ALS","DYN_S2_ATBD","STATIC_S2_ATBD","DYN_S2_DOPT","DYN_S2_DOPT_LADOPT","DYN_S2_RESCALED")
T<-rbindlist(lapply(fam,function(s) row(s,getpp(s))))
print(T[order(-sl_R2)])
fwrite(T,"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table4_dig_ladopt_slope.csv")
cat("\nDONE -> Table4_dig_ladopt_slope.csv\n")
