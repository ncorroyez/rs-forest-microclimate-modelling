# ==============================================================================
# c1_h2_gaussian_sweep_chs41.R — PROTOTYPE: does vertical DIFFUSENESS (sigma) move
# ΔTmax, independent of peak height (mu)? Gaussian LAD parametrisation after
# Fujiwara et al. (2026, IJRS): peak position mu and width sigma (fractions of Hmax),
# systematically swept, renormalised to preserve LAI. CHS41-Rmerge, NO wind.
# Config as the H2 controlled test: Hmax=25 m, cover=0.87. ΔTmax time-matched.
#   Rscript c1_h2_gaussian_sweep_chs41.R "5"        # dry-run: LAI 5 only
#   Rscript c1_h2_gaussian_sweep_chs41.R "3,5"      # full prototype
# Caveat: MuSICA re-bins the profile onto ~10 veg layers, so sigma below ~0.10*Hmax
# (~2.5 m here) is not resolvable; grid stays at/above that.
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
args<-commandArgs(trailingOnly=TRUE)
LAI_LEVELS<-if(length(args)) as.numeric(strsplit(args[1],",")[[1]]) else c(3,5)
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; MREF<-macro_ref(FORC,ds)
HMAX<-25; FCOVER<-0.87
MU_GRID<-c(0.35,0.50,0.65)        # peak height as fraction of Hmax (bottom / mid / top)
SIGMA_GRID<-c(0.12,0.22,0.35)     # vertical diffuseness as fraction of Hmax (narrow / med / broad)
make_lad_gaussian<-function(mu,sigma,canopy_base=1){force(mu);force(sigma);force(canopy_base)
  function(plot_row,hmax=NULL,lai=NULL){
    if(is.null(hmax))hmax<-as.numeric(plot_row$Hmax);if(is.null(lai))lai<-as.numeric(plot_row$LAI)
    ml<-ceiling(hmax);h<-seq_len(ml);z<-(h-0.5)/hmax
    w<-ifelse(h>=canopy_base, exp(-((z-mu)^2)/(2*sigma^2)), 0);w[!is.finite(w)]<-0
    dens<-if(sum(w)>0) lai*w/sum(w) else rep(0,ml);data.frame(height=h,density=dens)}}
com_frac<-function(lad_fn,hmax=HMAX){p<-lad_fn(NULL,hmax=hmax,lai=1);sum(p$density*p$height)/(sum(p$density)*hmax)}
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
NCDIR<-"out_files/H2_gaussian_chs41"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
metric_of<-function(nc){m<-micro_hourly_at(nc,Z);if(is.null(m))return(NA_real_);delta_tmax_mean(m,MREF,ds)}
grid<-CJ(LAI=LAI_LEVELS,mu=MU_GRID,sigma=SIGMA_GRID)
R<-rbindlist(lapply(seq_len(nrow(grid)),function(i){g<-grid[i]
  ladf<-make_lad_gaussian(g$mu,g$sigma); prow<-data.frame(id_plot=sprintf("mu%g_s%g",g$mu,g$sigma),Hmax=HMAX,LAI=g$LAI,fCover=FCOVER)
  nc<-file.path(NCDIR,sprintf("LAI%g_mu%03d_s%03d.nc",g$LAI,round(g$mu*100),round(g$sigma*100)))
  sc<-list(lai_fn=function(p)g$LAI,hmax_fn=function(p)HMAX,fcover_fn=function(p)FCOVER,lad_fn=ladf,phenology_fn=mk_phen(g$LAI))
  if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(prow,sc,nc,FORC,MB,extra_setup=ABL),error=function(e)cat("ERR",g$LAI,g$mu,g$sigma,conditionMessage(e),"\n"))
  data.table(LAI=g$LAI,mu=g$mu,sigma=g$sigma,com=round(com_frac(ladf),3),dtmax=metric_of(nc))}),fill=TRUE)
fwrite(R,"out_files/Chapter1/tables/h2_gaussian_chs41.csv")
cat("nc NA:",sum(is.na(R$dtmax)),"/",nrow(R),"\n\n")
for(l in LAI_LEVELS){cat(sprintf("=== LAI %g ===\n",l))
  cat("ΔTmax matrix (rows mu, cols sigma):\n")
  m<-dcast(R[LAI==l],mu~sigma,value.var="dtmax");print(m)
  se<-R[LAI==l,.(sigma_effect=round(max(dtmax)-min(dtmax),3)),by=mu]      # spread across sigma at fixed mu
  me<-R[LAI==l,.(mu_effect=round(max(dtmax)-min(dtmax),3)),by=sigma]      # spread across mu at fixed sigma
  cat(sprintf("  SIGMA effect (ΔTmax range across diffuseness at fixed mu): %.3f to %.3f °C\n",min(se$sigma_effect),max(se$sigma_effect)))
  cat(sprintf("  MU effect    (ΔTmax range across peak height at fixed sigma): %.3f to %.3f °C\n",min(me$mu_effect),max(me$mu_effect)))
  cat("  realised centroid (com) range: ",paste(range(R[LAI==l]$com),collapse=" to "),"\n\n")}
cat("DONE\n")
