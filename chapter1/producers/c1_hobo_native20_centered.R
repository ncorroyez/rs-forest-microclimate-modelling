# ==============================================================================
# HOBO validation with a 20 m SQUARE clipped EXACTLY on each logger (centered),
# traits computed natively from the LAS with the native-raster recipe
# (MacArthur-Horn dz=1, k=0.5, z0=1) + scan-angle + wind correction + iter.
# Tests whether the fixed-grid off-centering (loggers ~7.5 m off pixel center)
# and the 20 vs 25 m footprint matter vs the raster-buffer25 native validation.
# Out: out_files/hobo_native20_centered/musica_out_HOBO_<id>.nc
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(sf); library(data.table); library(lubridate)
  library(rmusica); library(musica.tools); library(parallel) })
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R"); options(lidR.progress=FALSE,lidR.verbose=FALSE)
CTG<-"/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
DZ<-1; Z0<-1; KEXT<-0.5; HALF<-10                         # 20 m square = +/-10 m; native recipe
ZB<-seq(1.5,39.5,by=1)
ctg<-readLAScatalog(CTG); opt_select(ctg)<-"xyzca"; opt_progress(ctg)<-FALSE
g<-st_read("in_files/data_Blois_utm31n.geojson",quiet=TRUE); g<-st_transform(g,32631)
xy<-st_coordinates(g); ids<-g$id_plot
OUT<-"out_files/hobo_native20_centered"; dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
res<-mclapply(seq_len(nrow(xy)),function(i){id<-ids[i]
  nc<-file.path(OUT,sprintf("musica_out_HOBO_%s.nc",id)); if(file.exists(nc)&&file.size(nc)>1e6)return(NULL)
  las<-tryCatch(suppressMessages(clip_rectangle(ctg,xy[i,1]-HALF,xy[i,2]-HALF,xy[i,1]+HALF,xy[i,2]+HALF)),error=function(e)NULL)
  if(is.null(las)||is.empty(las)||nrow(las@data)<50)return(NULL)
  las<-tryCatch(filter_duplicates(las),error=function(e)las)
  las<-tryCatch(suppressMessages(normalize_height(las,tin())),error=function(e)NULL); if(is.null(las))return(NULL)
  fin<-is.finite(las@data$Z)&las@data$Z>=0; Z<-las@data$Z[fin]; if(length(Z)<50)return(NULL)
  sa<-if("ScanAngle"%in%names(las@data)) las@data$ScanAngle[fin] else NULL
  sec<-if(!is.null(sa)&&any(is.finite(sa))) mean(1/cos(sa[is.finite(sa)]*pi/180)) else 1
  d<-tryCatch(lidR::LAD(Z,dz=DZ,k=KEXT,z0=Z0),error=function(e)NULL); if(is.null(d)||nrow(d)==0)return(NULL)
  lai1<-sum(d$lad,na.rm=TRUE)*DZ/sec                        # scan-angle corrected LAI
  hmax<-min(max(Z,na.rm=TRUE),40)
  fcov<-max(mean(Z>2,na.rm=TRUE),0.5)
  ladf<-data.frame(height=d$z,density=d$lad/sec)
  ff<-windcorr_forcing(hmax)
  sc<-list(lai_fn=function(p)lai1,hmax_fn=function(p)hmax,fcover_fn=function(p)fcov,
           lad_fn=function(plot_row,hmax=NULL,lai=NULL)ladf,phenology_fn=mk_phen(lai1))
  tryCatch(run_musica_one(data.frame(Cluster=1,x=xy[i,1],y=xy[i,2],id_plot=id),sc,nc,ff,MB,extra_setup=ABL),
           error=function(e)message("ERR ",id,": ",e$message))
  if(file.exists(nc)) message(sprintf("done %s (LAI=%.2f Hmax=%.1f sec=%.3f)",id,lai1,hmax,sec)); NULL
},mc.cores=6,mc.preschedule=FALSE)
cat(sprintf("CENTERED SQUARE nc: %d/%d\n",length(list.files(OUT,"\\.nc$")),nrow(xy)))
