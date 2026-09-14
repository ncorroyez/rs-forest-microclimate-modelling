# Native-20m correction chain: uncorrected (raw station wind, no scan-angle) and
# wind-only (windcorr, no scan-angle), for the §4.2 correction-chain numbers.
# Baseline native (windcorr + scan-angle) = c1_hobo_native20_pixel.R (r=0.42).
#   Rscript c1_hobo_native20_chain.R <uncorr|windonly>
suppressPackageStartupMessages({ library(terra); library(sf); library(ncdf4); library(lubridate); library(data.table)
  library(rmusica); library(musica.tools); library(parallel) })
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R"); source("R/wind_correction.R")
MODE<-commandArgs(trailingOnly=TRUE)[1]; IN<-"in_files_native20"
FORC<-"in_files/FR-Blo_2021_v2.nc"; MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
lai<-rast(file.path(IN,"lai_z1_res_10_m.tif"));hmax<-rast(file.path(IN,"max_res_10_m.tif"));fcov<-rast(file.path(IN,"fCover_res_10_m.tif"))
lad<-rast(file.path(IN,"lad_profiles_z1_res_10_m.tif"));sec<-rast(file.path(IN,"sec_theta_res_10_m.tif"));names(lad)<-sprintf("LAD_Layer_%.1f",seq(1.5,39.5,by=1))
g<-st_read(CFG$hobo_geojson,quiet=TRUE);g<-st_transform(g,crs(lai));v<-vect(g)
E<-cbind(LAI=terra::extract(lai,v)[,2],Hmax=terra::extract(hmax,v)[,2],fCover=terra::extract(fcov,v)[,2],sec=terra::extract(sec,v)[,2],terra::extract(lad,v)[,-1],crds(v))
E<-as.data.table(E);E$id_plot<-g$id_plot;E<-E[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
E[fCover<0.5,fCover:=0.5]                                       # NO scan-angle correction for the chain (uncorr/windonly)
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
OUT<-sprintf("out_files/hobo_native20_%s",MODE);dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
invisible(mclapply(seq_len(nrow(E)),function(i){pr<-as.data.frame(E[i]);id<-pr$id_plot
  nc<-file.path(OUT,sprintf("musica_out_HOBO_%s.nc",id));if(file.exists(nc)&&file.size(nc)>1e6)return(NULL)
  ff<-if(MODE=="windonly") windcorr_forcing(pr$Hmax) else FORC   # uncorr = raw station wind
  sc<-list(lai_fn=function(p)pr$LAI,hmax_fn=function(p)pr$Hmax,fcover_fn=function(p)pr$fCover,lad_fn=make_lad_real,phenology_fn=mk_phen(pr$LAI))
  tryCatch(run_musica_one(pr,sc,nc,ff,MB,extra_setup=ABL),error=function(e)message("ERR ",id));NULL
},mc.cores=8,mc.preschedule=FALSE))
cat(sprintf("%s nc: %d/%d\n",MODE,length(list.files(OUT,"\\.nc$")),nrow(E)))
