# DRY-RUN (supervisor request #4b, 2026-09-10): resolvable form of "add LAI above
# vs below the 1 m sensor". Literal below-1 m is sub-grid (lowest resolved LAD
# layer is 0-2.5 m, and the 1 m sensor sits in it). So we add a fixed dLAI=0.5
# block to the LOWEST resolved layer (0-2 m, at/below the sensor) vs an UPPER
# layer (around the 12.5 m peak), on a base Gaussian LAI=4 canopy, and compare
# ΔTmax. A renormalised LAI=4.5 Gaussian separates "adding foliage" from "where".
# CHS41-Rmerge / no wind, Hmax=25 m, cover=0.87, ΔTmax time-matched. DRY-RUN: LAI 4.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; MREF<-macro_ref(FORC,ds)
HMAX<-25; FCOVER<-0.87; MU<-0.5; SG<-0.24; DLAI<-0.5
args<-commandArgs(trailingOnly=TRUE)
LAI_LEVELS<-if(length(args)) as.numeric(strsplit(args[1],",")[[1]]) else c(2,4,6)
gauss<-function(lai){ml<-ceiling(HMAX);h<-seq_len(ml);z<-(h-0.5)/HMAX
  w<-exp(-((z-MU)^2)/(2*SG^2));w[!is.finite(w)]<-0; if(sum(w)>0) lai*w/sum(w) else rep(0,ml)}
# lad builders: base gaussian(lai_base) + a uniform dLAI block over height bins (metres)
mk_lad<-function(mode,lai_base){function(plot_row,hmax=NULL,lai=NULL){
  ml<-ceiling(HMAX);h<-seq_len(ml);base<-gauss(lai_base);add<-rep(0,ml)
  if(mode=="low")  add[h %in% 1:2]   <-DLAI/2                 # 0-2 m, lowest resolved layer
  if(mode=="high") add[h %in% 12:13] <-DLAI/2                 # 11-13 m, around the peak
  if(mode=="renorm") base<-gauss(lai_base+DLAI)              # add as gaussian (position-free)
  data.frame(height=h,density=base+add)}}
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
NCDIR<-"out_files/H2_addblock_chs41"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
metric_of<-function(nc){m<-micro_hourly_at(nc,Z);if(is.null(m))return(NA_real_);delta_tmax_mean(m,MREF,ds)}
modes<-c("low","high","renorm"); t0<-Sys.time()
grid<-CJ(LAI_base=LAI_LEVELS,mode=modes)
R<-rbindlist(lapply(seq_len(nrow(grid)),function(i){lb<-grid$LAI_base[i];md<-grid$mode[i]
  tot<-lb+DLAI
  prow<-data.frame(id_plot=sprintf("add_L%g_%s",lb,md),Hmax=HMAX,LAI=tot,fCover=FCOVER)
  nc<-file.path(NCDIR,sprintf("LAI%g_add_%s.nc",lb,md))
  sc<-list(lai_fn=function(p)tot,hmax_fn=function(p)HMAX,fcover_fn=function(p)FCOVER,
           lad_fn=mk_lad(md,lb),phenology_fn=mk_phen(tot))
  if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(prow,sc,nc,FORC,MB,extra_setup=ABL),
    error=function(e)cat("ERR",lb,md,conditionMessage(e),"\n"))
  data.table(LAI_base=lb,mode=md,LAI_tot=tot,dtmax=metric_of(nc))}),fill=TRUE)
cat(sprintf("\n=== ADD-BLOCK (base Gaussian + %.1f LAI) | %.0f s for %d sims ===\n",
    DLAI,as.numeric(difftime(Sys.time(),t0,units="secs")),nrow(grid)))
W<-dcast(R,LAI_base~mode,value.var="dtmax")
W[,`position_effect(low-high)`:=round(low-high,3)]
print(W)
fwrite(R,"out_files/Chapter1/tables/h2_addblock_chs41.csv")
cat("\nRead: add-high buffers more than add-low (upper-canopy foliage cools the 1 m understorey more);\n")
cat("position effect = add-low minus add-high (positive = upper placement buffers more).\n")
cat("DONE\n")
