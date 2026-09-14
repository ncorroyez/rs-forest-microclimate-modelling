# ==============================================================================
# HOBO validation on the NATIVE-20 m raster, single CONTAINING pixel per logger
# (fixed grid, off-center ~7.5 m) — isolates the CENTERING effect vs the centered
# 20 m square clip (c1_hobo_native20_centered.R), same 20 m / 400 m2 footprint.
# Reads in_files_native20 rasters directly (NOT build_hobo_inputs, which reads in_files).
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(data.table)
  library(rmusica); library(musica.tools); library(parallel) })
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R"); source("R/wind_correction.R")
IN<-"in_files_native20"; FORC<-"in_files/FR-Blo_2021_v2.nc"; MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
lai<-rast(file.path(IN,"lai_z1_res_10_m.tif")); hmax<-rast(file.path(IN,"max_res_10_m.tif"))
fcov<-rast(file.path(IN,"fCover_res_10_m.tif")); lad<-rast(file.path(IN,"lad_profiles_z1_res_10_m.tif")); sec<-rast(file.path(IN,"sec_theta_res_10_m.tif"))
names(lad)<-sprintf("LAD_Layer_%.1f",seq(1.5,39.5,by=1))
g<-st_read(CFG$hobo_geojson,quiet=TRUE); g<-st_transform(g,crs(lai)); v<-vect(g)
E<-cbind(terra::extract(lai,v)[,2,drop=FALSE], Hmax=terra::extract(hmax,v)[,2], fCover=terra::extract(fcov,v)[,2],
         sec=terra::extract(sec,v)[,2], terra::extract(lad,v)[,-1], crds(v))
names(E)[1]<-"LAI"; E$id_plot<-g$id_plot; E<-as.data.table(E)
E<-E[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
gsec<-mean(E$sec,na.rm=TRUE); E[,LAI:=LAI/gsec]; E[fCover<0.5,fCover:=0.5]        # global scan-angle (as native cLHS)
cat(sprintf("extracted %d loggers ; global sec=%.3f ; LAI mean=%.2f\n",nrow(E),gsec,mean(E$LAI)))
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
OUT<-"out_files/hobo_native20_pixel"; dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
invisible(mclapply(seq_len(nrow(E)),function(i){pr<-as.data.frame(E[i]);id<-pr$id_plot
  nc<-file.path(OUT,sprintf("musica_out_HOBO_%s.nc",id)); if(file.exists(nc)&&file.size(nc)>1e6)return(NULL)
  ff<-windcorr_forcing(pr$Hmax)
  sc<-list(lai_fn=function(p)pr$LAI,hmax_fn=function(p)pr$Hmax,fcover_fn=function(p)pr$fCover,lad_fn=make_lad_real,phenology_fn=mk_phen(pr$LAI))
  tryCatch(run_musica_one(pr,sc,nc,ff,MB,extra_setup=ABL),error=function(e)message("ERR ",id,": ",e$message)); NULL
},mc.cores=6,mc.preschedule=FALSE))
cat(sprintf("PIXEL nc: %d/%d\n",length(list.files(OUT,"\\.nc$")),nrow(E)))
