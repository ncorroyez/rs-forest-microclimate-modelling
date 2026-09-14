# SUPERSEDED 2026-09-09: Fig 6 moved to the Gaussian mu/sigma grid (c1_h2_gaussian_grid_chs41.R +
# c1_fig6v2_gaussian_heatmap.R). This beta-shape centroid sweep is kept for reference only.
# ==============================================================================
# c1_h2_controlled_chs41_clean.R — Fig B1 (→ main Fig 6): controlled H2 test on the
# authoritative forcing CHS41-Rmerge, NO wind correction. Rebuilt with the working
# run harness (the legacy config-based one failed to place ./forcing.nc on CHS41).
# 7 Beta LAD shapes (bottom->top-heavy) x 6 one-sided LAI, Hmax=25 m, cover=0.87 fixed.
# ΔTmax time-matched (delta_tmax_mean), CHS41 macro. Out: h2_controlled_chs41.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; MREF<-macro_ref(FORC,ds)
HMAX<-25; FCOVER<-0.87; LAI_LEVELS<-c(1,2,3,4,5,6)
PROFILES<-list(list(lab="very_bottom",a=2,b=6),list(lab="bottom",a=2,b=4),list(lab="mid_low",a=2,b=3),
  list(lab="centered",a=2,b=2),list(lab="mid_high",a=3,b=2),list(lab="top",a=4,b=2),list(lab="very_top",a=6,b=2))
make_lad_beta<-function(a,b,canopy_base=1){force(a);force(b);force(canopy_base)
  function(plot_row,hmax=NULL,lai=NULL){if(is.null(hmax))hmax<-as.numeric(plot_row$Hmax);if(is.null(lai))lai<-as.numeric(plot_row$LAI)
    ml<-ceiling(hmax);h<-seq_len(ml);u<-(h-canopy_base+0.5)/(ml-canopy_base+1)
    w<-ifelse(h>=canopy_base,dbeta(pmin(pmax(u,1e-3),1-1e-3),a,b),0);w[!is.finite(w)]<-0
    dens<-if(sum(w)>0) lai*w/sum(w) else rep(0,ml);data.frame(height=h,density=dens)}}
com_frac<-function(lad_fn,hmax=HMAX){p<-lad_fn(NULL,hmax=hmax,lai=1);sum(p$density*p$height)/(sum(p$density)*hmax)}
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
NCDIR<-"out_files/H2_controlled_chs41_clean"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
metric_of<-function(nc){m<-micro_hourly_at(nc,Z);if(is.null(m))return(NA_real_);delta_tmax_mean(m,MREF,ds)}
R<-rbindlist(lapply(LAI_LEVELS,function(lai) rbindlist(lapply(PROFILES,function(pf){
  ladf<-make_lad_beta(pf$a,pf$b); prow<-data.frame(id_plot=pf$lab,Hmax=HMAX,LAI=lai,fCover=FCOVER)
  nc<-file.path(NCDIR,sprintf("LAI%g_%s.nc",lai,pf$lab))
  sc<-list(lai_fn=function(p)lai,hmax_fn=function(p)HMAX,fcover_fn=function(p)FCOVER,lad_fn=ladf,phenology_fn=mk_phen(lai))
  if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(prow,sc,nc,FORC,MB,extra_setup=ABL),error=function(e)cat("ERR",lai,pf$lab,conditionMessage(e),"\n"))
  data.table(LAI=lai,profile=pf$lab,com=round(com_frac(ladf),3),dtmax=metric_of(nc))
}))),fill=TRUE)
fwrite(R,"out_files/Chapter1/tables/h2_controlled_chs41.csv")
cat("nc NA:",sum(is.na(R$dtmax)),"/",nrow(R),"\n")
rng<-R[,.(dt_min=round(min(dtmax,na.rm=T),3),dt_max=round(max(dtmax,na.rm=T),3),
          com_effect=round(dtmax[which.max(com)]-dtmax[which.min(com)],3)),by=LAI][order(LAI)]
cat("\n=== H2 controlled (CHS41): per LAI, ΔTmax range across shapes + top-minus-bottom ===\n"); print(rng)
cat(sprintf("\nProfile-shape ΔTmax effect spans %.2f to %.2f °C across LAI\n",
    min(R[,max(dtmax)-min(dtmax),by=LAI]$V1), max(R[,max(dtmax)-min(dtmax),by=LAI]$V1)))
cat("DONE\n")
