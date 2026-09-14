# ==============================================================================
# c1_perturb_chs41_nowind.R
#
# Fig 4 / Table G attribution, rerun on the authoritative forcing: CHS41-Rmerge
# (h_sbl median 550 m, MERRA-2), NO wind correction (native20 config; wind becomes
# an appendix sensitivity). The shipped harness (c1_sensitivity_perplot_units.R)
# used FR-Blo_2021_v2.nc (h_sbl 48.8 m) + wind correction; verified that the low
# h_sbl flattens the mid-density LAI lever (P2 −0.105 -> −0.146; dense P3/P4
# unchanged). Native-unit 8-sim design, identical steps/bounds to the shipped one.
#
# Native cLHS profiles (from the rds, not re-clipped), keyed clhs_%03d, full 400.
# Shell-parallel by chunk (no mclapply/fork). Resumable (skip existing nc).
#   Rscript c1_perturb_chs41_nowind.R <chunk> <K>
# Out: out_files/Chapter1/tables/perturb_chs41_nowind_part<chunk>of<K>.csv
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
source("Chapter3_config.R"); CFG_C3$musica_cmd<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE)
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"; ABL<-list("abl_flag"='"iter"')
.forc<-FORC                                     # NO wind correction (native20; wind = appendix)
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z_FIX<-1.0
MREF<-macro_ref(FORC,ds)

args<-commandArgs(trailingOnly=TRUE)
CHUNK<-if(length(args)>=1) as.integer(args[1]) else 1L
K    <-if(length(args)>=2) as.integer(args[2]) else 1L
NMAX <-{e<-Sys.getenv("DRY_NMAX"); if(nzchar(e)) as.integer(e) else Inf}

D_LAI<-0.5; D_HMAX<-1.0; D_FCOV<-0.10
LAI_FLOOR<-0.1; HMAX_FLOOR<-2; HMAX_CEIL<-40; FCOV_MIN<-0.5; FCOV_MAX<-1
mk_phen<-function(lai1)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,
  leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=lai1))
  d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}

samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[,id_plot:=sprintf("clhs_%03d",seq_len(.N))]; samp[,P:=as.character(relabel_cluster(Cluster))]
samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
np<-min(nrow(samp),NMAX)
NCDIR<-"out_files/Chapter1/nc_perturb_chs41_nowind"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
metric_of<-function(nc){m<-micro_hourly_at(nc,Z_FIX);if(is.null(m))return(NA_real_);delta_tmax_mean(m,MREF,ds)}
per_unit<-function(sim,base,step) if(!is.finite(sim)||!is.finite(base)||abs(step)<1e-9) NA_real_ else (sim-base)/step
run_get<-function(prow,tag,lai,hmax,fcov,ladf=make_lad_real){
  nc<-file.path(NCDIR,sprintf("%s_%s.nc",prow$id_plot,tag))
  sc<-list(lai_fn=function(p)lai,hmax_fn=function(p)hmax,fcover_fn=function(p)fcov,lad_fn=ladf,phenology_fn=mk_phen(lai))
  if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(prow,sc,nc,.forc,CFG_C3$musica_cmd,extra_setup=ABL),error=function(e)NULL)
  metric_of(nc)}

my<-which(((seq_len(np)-1L)%%K)==(CHUNK-1L))
cat(sprintf("Fig4 rerun CHS41/no-wind | chunk %d/%d | %d plots | 8 sims/plot\n",CHUNK,K,length(my)))
R<-rbindlist(lapply(my,function(i){prow<-as.data.frame(samp[i])
  r<-tryCatch({
    base<-run_get(prow,"base",prow$LAI,prow$Hmax,prow$fCover)
    unif<-run_get(prow,"unifLAD",prow$LAI,prow$Hmax,prow$fCover,ladf=make_lad_uniform)
    laip_v<-prow$LAI+D_LAI; laim_v<-max(prow$LAI-D_LAI,LAI_FLOOR)
    hmp_v<-min(prow$Hmax+D_HMAX,HMAX_CEIL); hmm_v<-max(prow$Hmax-D_HMAX,HMAX_FLOOR)
    fcp_v<-min(prow$fCover+D_FCOV,FCOV_MAX); fcm_v<-max(prow$fCover-D_FCOV,FCOV_MIN)
    laip<-run_get(prow,"LAIp",laip_v,prow$Hmax,prow$fCover); laim<-run_get(prow,"LAIm",laim_v,prow$Hmax,prow$fCover)
    hmp<-run_get(prow,"Hmaxp",prow$LAI,hmp_v,prow$fCover);   hmm<-run_get(prow,"Hmaxm",prow$LAI,hmm_v,prow$fCover)
    fcp<-run_get(prow,"fCovp",prow$LAI,prow$Hmax,fcp_v);     fcm<-run_get(prow,"fCovm",prow$LAI,prow$Hmax,fcm_v)
    data.table(id_plot=prow$id_plot,P=prow$P,LAI=prow$LAI,Hmax=prow$Hmax,fCover=prow$fCover,base=base,
      LAI_up=per_unit(laip,base,laip_v-prow$LAI), LAI_dn=per_unit(laim,base,laim_v-prow$LAI),
      Hmax_up=per_unit(hmp,base,hmp_v-prow$Hmax), Hmax_dn=per_unit(hmm,base,hmm_v-prow$Hmax),
      fCov_up=per_unit(fcp,base,fcp_v-prow$fCover), fCov_dn=per_unit(fcm,base,fcm_v-prow$fCover),
      dT_LAD=if(is.finite(base)&&is.finite(unif)) base-unif else NA_real_)
  },error=function(e){cat(sprintf(" ERR %s: %s\n",prow$id_plot,conditionMessage(e)));NULL})
  if(!is.null(r)) cat(sprintf(" %s done\n",prow$id_plot)); r}),fill=TRUE)
odir<-"out_files/Chapter1/tables"; dir.create(odir,recursive=TRUE,showWarnings=FALSE)
of<-if(K==1L) file.path(odir,"perturb_chs41_nowind.csv") else file.path(odir,sprintf("perturb_chs41_nowind_part%02dof%02d.csv",CHUNK,K))
fwrite(R,of); cat("wrote",of,"\n")
