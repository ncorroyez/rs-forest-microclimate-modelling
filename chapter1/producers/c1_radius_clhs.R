# ==============================================================================
# Radius sensitivity on the cLHS design (station forcing, NO wind correction,
# + scan-angle-corrected LAI + v3.2.3 iter), so it matches the reported model.
# Per cLHS plot: clip_circle(radius) -> LAD/LAI/Hmax/fCover + <sec theta> -> MuSICA.
#   Rscript c1_radius_corr.R <radius_m> [nmax]
# Out: out_files/radius_test_corr/<r>m/musica_out_HOBO_<id>.nc
# ==============================================================================
suppressPackageStartupMessages({
  library(lidR); library(terra); library(sf); library(data.table); library(lubridate)
  library(rmusica); library(musica.tools); library(parallel) })
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("pipeline/00_config.R")
options(lidR.progress=FALSE)
args<-commandArgs(trailingOnly=TRUE); radius<-as.numeric(args[1]); nmax<-if(length(args)>=2) as.integer(args[2]) else Inf
CTG<-"/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
DZ<-0.5; Z0<-0.5; KEXT<-0.5
ctg<-readLAScatalog(CTG); opt_select(ctg)<-"*"; opt_progress(ctg)<-FALSE
# cLHS design instead of the instrumented plots: Chapter 1 rests on the 400 cLHS
# samples only. Plot ids are the row index, prefixed so the file names stay unique.
CL <- readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds")
if (nzchar(Sys.getenv("NPER"))) {                 # optional: N per archetype
  n_per <- as.integer(Sys.getenv("NPER")); set.seed(42)
  CL <- do.call(rbind, lapply(split(CL, CL$Cluster),
                              function(d) d[sample(seq_len(nrow(d)), min(n_per, nrow(d))), ]))
}
xy  <- as.matrix(CL[, c("x", "y")])
ids <- sprintf("clhs_%03d", seq_len(nrow(CL)))
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
OUT<-file.path("out_files/radius_test_clhs",paste0(radius,"m")); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
np<-min(nrow(xy),nmax); cat(sprintf("RADIUS %s m (corrected) — %d plots\n",radius,np))
invisible(mclapply(seq_len(np), function(i){id<-ids[i]
  nc<-file.path(OUT,sprintf("musica_out_%s.nc",id)); if(file.exists(nc)&&file.size(nc)>1e6) return(NULL)
  las<-tryCatch(suppressMessages(clip_circle(ctg,xy[i,1],xy[i,2],radius)),error=function(e)NULL)
  if(is.null(las)||is.empty(las)||nrow(las@data)<50){message(" skip ",id);return(NULL)}
  las<-tryCatch(filter_duplicates(las),error=function(e)las)
  las<-tryCatch(suppressMessages(normalize_height(las,knnidw())),error=function(e)NULL); if(is.null(las))return(NULL)
  fin<-is.finite(las@data$Z)&las@data$Z>=0; Z<-las@data$Z[fin]; if(length(Z)<50)return(NULL)
  sa<-if("ScanAngle"%in%names(las@data)) las@data$ScanAngle[fin] else NULL
  sec<-if(!is.null(sa)&&any(is.finite(sa))) mean(1/cos(sa[is.finite(sa)]*pi/180)) else 1
  d<-tryCatch(lidR::LAD(Z,dz=DZ,k=KEXT,z0=Z0),error=function(e)NULL); if(is.null(d)||nrow(d)==0)return(NULL)
  lai1<-sum(d$lad,na.rm=TRUE)*DZ/sec                    # one-sided LAI, scan-angle corrected
  hmax<-min(max(Z,na.rm=TRUE),40)
  fcov<-max(mean(Z>2,na.rm=TRUE),0.5)                   # canopy-return fraction (floored 0.5)
  ladf<-data.frame(height=d$z, density=d$lad/sec)       # corrected LAD profile (renormalised in R/lad.R to lai1)
  ff<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"   # no wind correction (Ch1/Ch3 aligned)
  sc<-list(lai_fn=function(p)lai1,hmax_fn=function(p)hmax,fcover_fn=function(p)fcov,
           lad_fn=function(plot_row,hmax=NULL,lai=NULL)ladf,phenology_fn=mk_phen(lai1))
  tryCatch(run_musica_one(data.frame(Cluster=1,x=xy[i,1],y=xy[i,2],id_plot=id),sc,nc,ff,MB,extra_setup=ABL),
           error=function(e)message(" ERR ",id,": ",e$message))
  if(file.exists(nc)) message(sprintf(" done %s r%s (LAI=%.2f Hmax=%.1f sec=%.3f)",id,radius,lai1,hmax,sec))
  NULL
}, mc.cores=10, mc.preschedule=FALSE))
cat(sprintf("RADIUS %s DONE\n",radius))
