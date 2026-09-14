# ==============================================================================
# Footprint-radius sweep: which LiDAR extraction radius best validates MuSICA
# v3.2.3 iter + FR-Blo against the 53 HOBO loggers. For each radius, extract the
# logger canopy traits (LAI/Hmax/fCover/LAD) within that circular buffer, run the
# REF (real-canopy) sim, and score ΔTmax + slope vs observed (station 1.5 m ref).
#   Radii: 5, 10, 12.5, 15, 20, 25, 50 m. 53 × 7 = 371 sims (REF only).
#   Rscript c1_radius_sweep_validation.R [radius]   (one radius, or all)
# Out: tab_radius_sweep_validation.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(sf); library(terra)
  library(parallel); library(rmusica); library(musica.tools)
  src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("pipeline/00_config.R") })
FORC<-"in_files/FR-Blo_2021_v2.nc"; MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
RADII<-c(5,10,12.5,15,20,25,50); arg<-commandArgs(trailingOnly=TRUE); if(length(arg)>0) RADII<-as.numeric(arg)
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; SH<-2L
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=F];d$Julian_day<-366;rbind(ph,d)}

# raw 10 m rasters (buffer25 source)
r_lai<-rast("in_files/lai_z1_res_10_m.tif"); r_hmax<-rast("in_files/max_res_10_m.tif")
r_fcov<-rast("in_files/fCover_res_10_m.tif"); r_lad<-rast("in_files/lad_profiles_z1_res_10_m.tif")
hpts<-st_read(CFG$hobo_geojson,quiet=TRUE)%>%filter(!id_plot%in%CFG$ids_to_remove)

extract_radius<-function(rad){
  bv<-vect(st_buffer(hpts,rad))
  LAI<-terra::extract(r_lai,bv,fun="mean",na.rm=T)[,2]; Hmax<-terra::extract(r_hmax,bv,fun="max",na.rm=T)[,2]
  fC<-terra::extract(r_fcov,bv,fun="mean",na.rm=T)[,2]; lad<-terra::extract(r_lad,bv,fun="mean",na.rm=T)
  d<-data.frame(id_plot=hpts$id_plot,LAI=LAI,Hmax=Hmax,fCover=pmax(fC,0.5))
  ladm<-lad[,-1,drop=F]; zc<-as.numeric(gsub("[^0-9.]","",names(ladm)))
  for(k in seq_along(zc)) d[[sprintf("LAD_Layer_%s",zc[k])]]<-ladm[,k]
  d[is.finite(d$LAI)&is.finite(d$Hmax)&is.finite(d$fCover),] }

# macro (station) + obs
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc);macH<-macH[as.Date(time)%in%ds]
macD<-macH[,.(Tmax_mac=max(Tmac)),by=.(date=as.Date(time))]
hobo<-as.data.table(read.csv(CFG$hobo_temp_csv));hobo[,datetime:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hobo<-hobo[position_sensor=="a"&as.Date(datetime)%in%ds&!id_plot%in%CFG$ids_to_remove];hobo[,time:=floor_date(datetime,"hour")]
oh<-hobo[,.(Tmic=mean(t_hobo,na.rm=T)),by=.(id_plot,time)]; od<-hobo[,.(Tmax_mic=max(t_hobo,na.rm=T)),by=.(id_plot,date=as.Date(datetime))]
obs_dt<-merge(od,macD,by="date")[,.(obs_dt=mean(Tmax_mic-Tmax_mac,na.rm=T)),by=id_plot]
obs_sl<-merge(oh,macH,by="time")[,.(obs_sl=coef(lm(Tmic~Tmac))[2]),by=id_plot]
ext<-function(p){if(!file.exists(p)||file.size(p)<1e5)return(list(dt=NA,sl=NA));nc<-try(nc_open(p),silent=T);if(inherits(nc,"try-error"))return(list(dt=NA,sl=NA));on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var)))return(list(dt=NA,sl=NA))
  tu<-ncatt_get(nc,"time","units")$value;t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC");th<-ncvar_get(nc,"time");Tk<-ncvar_get(nc,"Tair_z");rh<-ncvar_get(nc,"relative_height");vh<-median(ncvar_get(nc,"veget_height_top"),na.rm=T);zl<-rh*vh
  if(Z<=zl[1]){il<-1;ih<-1;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1;w<-(Z-zl[il])/(zl[ih]-zl[il])};Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
  tv<-t0+dhours(th)-hours(SH);dd<-data.table(date=as.Date(floor_date(tv,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date];dt<-merge(dd,macD,by="date")[,mean(Tmax-Tmax_mac,na.rm=T)]
  mm<-merge(data.table(time=floor_date(t0+dhours(th),"hour"),Tmic=Tc),macH,by="time");sl<-if(nrow(mm)>50)coef(lm(Tmic~Tmac,mm))[2] else NA;list(dt=dt,sl=sl)}

RES<-list()
for(rad in RADII){
  df<-extract_radius(rad); OUT<-sprintf("out_files/musica_hobo_radius/r%s",gsub("\\.","p",rad)); dir.create(OUT,recursive=T,showWarnings=F)
  invisible(mclapply(seq_len(nrow(df)),function(i){pr<-df[i,,drop=F];nc<-file.path(OUT,sprintf("musica_out_HOBO_%s.nc",pr$id_plot));if(file.exists(nc)&&file.size(nc)>1e6)return(NULL)
    sc<-list(lai_fn=function(p)pr$LAI,hmax_fn=function(p)pr$Hmax,fcover_fn=function(p)pr$fCover,lad_fn=make_lad_real,phenology_fn=mk_phen(pr$LAI))
    tryCatch(run_musica_one(pr,sc,nc,FORC,MB,extra_setup=ABL),error=function(e)cat(sprintf("ERR r%s %s\n",rad,pr$id_plot)))},mc.cores=4,mc.preschedule=F))
  E<-rbindlist(lapply(df$id_plot,function(id){e<-ext(file.path(OUT,sprintf("musica_out_HOBO_%s.nc",id)));data.table(id_plot=id,sim_dt=e$dt,sim_sl=e$sl)}))
  V<-merge(merge(E,obs_dt,by="id_plot"),obs_sl,by="id_plot"); ok<-is.finite(V$sim_dt)&is.finite(V$obs_dt);oks<-is.finite(V$sim_sl)&is.finite(V$obs_sl)
  RES[[as.character(rad)]]<-data.table(radius=rad, n=sum(ok), LAI_med=median(df$LAI), LAI_max=max(df$LAI),
    r=cor(V$sim_dt[ok],V$obs_dt[ok]), R2=cor(V$sim_dt[ok],V$obs_dt[ok])^2,
    RMSE=sqrt(mean((V$sim_dt[ok]-V$obs_dt[ok])^2)), MAE=mean(abs(V$sim_dt[ok]-V$obs_dt[ok])),
    bias=mean(V$sim_dt[ok]-V$obs_dt[ok]), slope_r=cor(V$sim_sl[oks],V$obs_sl[oks]))
  cat(sprintf("r=%4sm : LAI med %.1f max %.1f | ΔTmax r=%.3f R2=%.2f RMSE=%.2f MAE=%.2f bias=%+.2f | slope_r=%.3f\n",
    rad,median(df$LAI),max(df$LAI),RES[[as.character(rad)]]$r,RES[[as.character(rad)]]$R2,RES[[as.character(rad)]]$RMSE,RES[[as.character(rad)]]$MAE,RES[[as.character(rad)]]$bias,RES[[as.character(rad)]]$slope_r))
}
tab<-rbindlist(RES); fwrite(tab,"out_files/Chapter1/tables/tab_radius_sweep_validation.csv")
cat("\nDONE -> tab_radius_sweep_validation.csv\n")
