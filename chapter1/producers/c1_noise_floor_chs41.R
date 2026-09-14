# ==============================================================================
# Numerical noise floor (±0.01 LAI control) on the AUTHORITATIVE forcing:
# CHS41-Rmerge, NO wind correction (baseline config of every Chapter 1 result).
# Perturb LAI by a physically negligible ±0.01 and measure the ΔTmax change; it
# sets the threshold below which a lever cannot be called real.
#   Rscript c1_noise_floor_chs41.R [n_plots]
# Out: out_files/Chapter1/tables/tab_noise_floor_chs41.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table)
  library(parallel);library(rmusica);library(musica.tools)})
args<-commandArgs(trailingOnly=TRUE); NP<-if(length(args)) as.integer(args[1]) else 40
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE)
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"; ABL<-list("abl_flag"='"iter"')
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); MREF<-macro_ref(FORC,ds)
samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub<-samp[,.SD[sample(.N,min(.N,100))],by=Cluster]; sub[,pid:=sprintf("S%04d",.I)]
set.seed(11); pick<-sub[,.SD[sample(.N,ceiling(NP/4))],by=Cluster]
NCDIR<-"out_files/Chapter1/nc_noise_floor_chs41"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,
  leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l))
  d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
one<-function(i){
  p<-as.data.frame(pick[i])                                   # NO wind correction: drive with FORC
  g<-function(tag,lai){nc<-file.path(NCDIR,sprintf("%s_%s.nc",p$pid,tag))
    sc<-list(lai_fn=function(x)lai,hmax_fn=function(x)p$Hmax,fcover_fn=function(x)p$fCover,
             lad_fn=make_lad_real,phenology_fn=mk_phen(lai))
    if(!file.exists(nc)||file.size(nc)<1000) run_musica_one(p,sc,nc,FORC,MB,extra_setup=ABL)
    m<-micro_hourly_at(nc,1.0); if(is.null(m)) NA_real_ else delta_tmax_mean(m,MREF,ds)}
  data.table(pid=p$pid,Cluster=p$Cluster,LAI=p$LAI,P=as.character(relabel_cluster(p$Cluster)),
             base=g("base",p$LAI), up=g("eps_p",p$LAI+0.01), dn=g("eps_m",p$LAI-0.01))}
NCORES<-if(length(args)>=2) as.integer(args[2]) else 4
R<-rbindlist(mclapply(seq_len(nrow(pick)),function(i) tryCatch(one(i),error=function(e){cat("ERR",i,conditionMessage(e),"\n");NULL}),mc.cores=NCORES),fill=TRUE)
R<-R[is.finite(up)&is.finite(dn)][,noise:=up-dn]
fwrite(R,"out_files/Chapter1/tables/tab_noise_floor_chs41.csv")
cat(sprintf("\n=== NOISE FLOOR CHS41/no-wind (LAI +-0.01, n=%d plots) ===\n",nrow(R)))
cat(sprintf("  |dTmax(+0.01)-dTmax(-0.01)| : median %.5f | q95 %.5f | max %.5f degC\n",
  median(abs(R$noise)),quantile(abs(R$noise),0.95),max(abs(R$noise))))
cat(sprintf("  median |response| per 0.02 step: %.5f (linear extrap of +-0.5 step ~ 0.02/0.5 * P4-lever)\n",median(abs(R$noise))))
for(pp in c("P1","P2","P3","P4")){d<-R[P==pp]; if(nrow(d)) cat(sprintf("  %s: median noise %.5f (n=%d)\n",pp,median(d$noise),nrow(d)))}
cat(sprintf("  P4 leaf-area lever (0.5 step, median) = 0.205 degC ; ratio lever/noise ~ %.0fx\n",0.205/median(abs(R$noise))))
