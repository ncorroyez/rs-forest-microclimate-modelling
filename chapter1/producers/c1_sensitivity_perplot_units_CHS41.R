# ==============================================================================
# VARIANTE CHS41 de c1_sensitivity_perplot_units.R (original intact).
# Deux changements et deux seulement : forcage = hybride CHS41 (h_sbl MERRA-2,
# mediane 550 m, contre 49 m dans FR-Blo_2021_v2.nc) et correction de vent
# desactivee, pour aligner l attribution sur la convention retenue pour les runs
# CHS41 : validation aux 53 capteurs (chapitre 3) et balayage de rayon (chapitre 1).
# Sorties dans des repertoires distincts ; rien n est ecrase.
#
# NATIVE-UNIT sensitivity design (committee 2026-07-28): replaces the ±1 SD steps
# by fixed, interpretable steps so effects read as degrees per unit:
#     LAI  ± 0.5 (one-sided)      Hmax ± 1 m       fCover ± 0.10 (10 points)
# The + and the − side are stored separately (never averaged), and the profile is
# still the real-vs-uniform swap. ΔTmax uses the canonical time-matched convention
# (R/dtmax_convention.R): sub-canopy read at the hour of the macro daily maximum,
# NO clock shift.
#   Rscript c1_sensitivity_perplot_units.R <chunk> <K>
# Out: out_files/Chapter1/tables/sensitivity_perplot_units_CHS41/part_<chunk>.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools) })
args<-commandArgs(trailingOnly=TRUE); chunk<-as.integer(args[1]); K<-as.integer(args[2])
src<-list.files("R",pattern="\\.R$",full.names=TRUE); src<-src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source)); source("Chapter3_config.R")
CFG_C3$musica_cmd<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE)
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"; ABL<-list("abl_flag"='"iter"')
WC<-FALSE; source("R/wind_correction.R"); .forc<-FORC   # sans correction de vent, convention CHS41
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z_FIX<-1.0
MREF<-macro_ref(FORC,ds)                       # canonical macro: daily max + its hour

samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub<-samp[,.SD[sample(.N,min(.N,100))],by=Cluster]; sub[,pid:=sprintf("S%04d",.I)]
my<-which(((seq_len(nrow(sub))-1)%%K)==(chunk-1))
cat(sprintf("units chunk %d/%d : %d plots\n",chunk,K,length(my)))

# --- fixed native steps -------------------------------------------------------
D_LAI<-0.5; D_HMAX<-1.0; D_FCOV<-0.10
# Physical bounds on the perturbed traits. A step that would cross one is truncated,
# and the effect reported is the change ACTUALLY simulated, never rescaled to the
# nominal step. HMAX_FLOOR corrected 3 -> 2 m and HMAX_CEIL added (2026-07-30, author):
# the shipped 3200 runs used a 3 m floor and no ceiling, which over-truncated 2 plots
# and let 2 others reach 41 m. With the height lever at or below 0.004 degC per metre
# that is immaterial to every reported median, so those runs were NOT redone.
LAI_FLOOR<-0.1; HMAX_FLOOR<-2; HMAX_CEIL<-40; FCOV_MIN<-0.5; FCOV_MAX<-1
mk_phen<-function(lai1)function(p){
  ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,
      leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=lai1))
  d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE]; d$Julian_day<-366; rbind(ph,d)}
NCDIR<-"out_files/Chapter1/nc_sensitivity_perplot_units_CHS41"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
metric_of<-function(nc){ m<-micro_hourly_at(nc,Z_FIX); if(is.null(m)) return(NA_real_)
  delta_tmax_mean(m,MREF,ds) }
run_get<-function(prow,tag,lai,hmax,fcov,ladf=make_lad_real){
  nc<-file.path(NCDIR,sprintf("%s_%s.nc",prow$pid,tag))
  sc<-list(lai_fn=function(p)lai,hmax_fn=function(p)hmax,fcover_fn=function(p)fcov,
           lad_fn=ladf,phenology_fn=mk_phen(lai))
  if(!file.exists(nc)||file.size(nc)<1000) run_musica_one(prow,sc,nc,.forc,CFG_C3$musica_cmd,extra_setup=ABL)
  metric_of(nc)}
per_unit<-function(sim,base,step){ if(!is.finite(sim)||!is.finite(base)||!is.finite(step)||abs(step)<1e-9) return(NA_real_); (sim-base)/step }

pdir<-"out_files/Chapter1/tables/sensitivity_perplot_units_CHS41"; dir.create(pdir,recursive=TRUE,showWarnings=FALSE)
pf<-file.path(pdir,sprintf("part_%02d.csv",chunk)); if(file.exists(pf)) file.remove(pf)
nok<-0L
for(j in seq_along(my)){
  prow<-as.data.frame(sub[my[j]])
  .forc<<-if(WC) windcorr_forcing(prow$Hmax) else FORC
  r<-tryCatch({
    base<-run_get(prow,"base",prow$LAI,prow$Hmax,prow$fCover)
    unif<-run_get(prow,"unifLAD",prow$LAI,prow$Hmax,prow$fCover,ladf=make_lad_uniform)
    laip_v<-prow$LAI+D_LAI;               laim_v<-max(prow$LAI-D_LAI,LAI_FLOOR)
    hmp_v <-min(prow$Hmax+D_HMAX,HMAX_CEIL); hmm_v <-max(prow$Hmax-D_HMAX,HMAX_FLOOR)
    fcp_v <-min(prow$fCover+D_FCOV,FCOV_MAX); fcm_v<-max(prow$fCover-D_FCOV,FCOV_MIN)
    laip<-run_get(prow,"LAIp",laip_v,prow$Hmax,prow$fCover)
    laim<-run_get(prow,"LAIm",laim_v,prow$Hmax,prow$fCover)
    hmp <-run_get(prow,"Hmaxp",prow$LAI,hmp_v,prow$fCover)
    hmm <-run_get(prow,"Hmaxm",prow$LAI,hmm_v,prow$fCover)
    fcp <-run_get(prow,"fCovp",prow$LAI,prow$Hmax,fcp_v)
    fcm <-run_get(prow,"fCovm",prow$LAI,prow$Hmax,fcm_v)
    data.table(pid=prow$pid,Cluster=prow$Cluster,LAI=prow$LAI,Hmax=prow$Hmax,fCover=prow$fCover,
      base=base,
      # degrees per unit, + side and - side kept SEPARATE
      LAI_up  =per_unit(laip,base, laip_v-prow$LAI),  LAI_dn  =per_unit(laim,base, laim_v-prow$LAI),
      Hmax_up =per_unit(hmp ,base, hmp_v -prow$Hmax), Hmax_dn =per_unit(hmm ,base, hmm_v -prow$Hmax),
      fCov_up =per_unit(fcp ,base, fcp_v -prow$fCover),fCov_dn=per_unit(fcm ,base, fcm_v -prow$fCover),
      step_LAI=D_LAI, step_Hmax=D_HMAX, step_fCov=D_FCOV,
      dT_LAD  =base-unif)
  },error=function(e) NULL)
  if(!is.null(r)){ fwrite(r,pf,append=file.exists(pf)); nok<-nok+1L }
  if(j%%5==0) cat(sprintf("  units chunk %d: %d/%d\n",chunk,j,length(my)))
}
cat(sprintf("units chunk %d DONE (%d/%d)\n",chunk,nok,length(my)))
