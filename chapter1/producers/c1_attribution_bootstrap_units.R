# ==============================================================================
# Bootstrap medians and 95% intervals of the native-unit sensitivities (Table G2).
#
# Reads : out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv   (stage A3)
# Writes: out_files/Chapter1/tables/tab_attribution_units_bootstrap.csv
#   Rscript c1_attribution_bootstrap_units.R
# ==============================================================================
# Bootstrap CI on the per-archetype median trait effect, NATIVE-UNIT design.
# Effects are the change actually simulated for a fixed step (LAI ±0.5, cover ±10
# points, height ±1 m); the + and the − side are kept separate. Replaces the former
# per-SD bootstrap, which is no longer used anywhere.
suppressPackageStartupMessages({library(data.table);source("R/cluster_relabel.R")})
LV<-c("P1","P2","P3","P4"); set.seed(1)
# CANONICAL SOURCE 2026-07-30: read the v2 table (the one Table G2, Fig. 4 and Fig. 5
# use) instead of the v1 parts. They agree today, but a pipeline re-run of v2 would have
# left Table G1's intervals matching nothing. Column names are mapped to the v1 names.
U<-fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
setnames(U,c("LAI_up_dt_all","LAI_dn_dt_all","Hmax_up_dt_all","Hmax_dn_dt_all",
             "fCov_up_dt_all","fCov_dn_dt_all","LAD_dt_all","base_dt_all"),
           c("LAI_up","LAI_dn","Hmax_up","Hmax_dn","fCov_up","fCov_dn","dT_LAD","base"))
U[,`:=`(step_LAI=0.5, step_Hmax=1.0, step_fCov=0.10)]
U[,P:=factor(relabel_cluster(Cluster),levels=LV)]
U[,`:=`(stepLp=step_LAI, stepLm=LAI-pmax(LAI-step_LAI,0.1),
        stepHp=step_Hmax, stepHm=Hmax-pmax(Hmax-step_Hmax,3),
        stepFp=pmin(fCover+step_fCov,1)-fCover, stepFm=fCover-pmax(fCover-step_fCov,0.5))]
U[,`:=`(LAI_plus=LAI_up*stepLp, LAI_minus=LAI_dn*(-stepLm),
        fCov_plus=fCov_up*stepFp, fCov_minus=fCov_dn*(-stepFm),
        Hmax_plus=Hmax_up*stepHp, Hmax_minus=Hmax_dn*(-stepHm))]
#' Median and 95% bootstrap interval of a per-plot effect vector
#' @param x numeric vector of effects, non-finite values dropped
#' @param B number of bootstrap replicates (default 4000, the count Table G1's caption states)
#' @return numeric c(median, lo, hi); all NA when fewer than 5 finite values
bmed<-function(x,B=4000){x<-x[is.finite(x)];if(length(x)<5)return(c(NA,NA,NA))
  b<-replicate(B,median(sample(x,length(x),TRUE)));c(median(x),quantile(b,.025),quantile(b,.975))}
VARS<-c("LAI_plus","LAI_minus","fCov_plus","fCov_minus","Hmax_plus","Hmax_minus","dT_LAD")
LAB<-c("LAI +0.5","LAI −0.5","fCover +10 pts","fCover −10 pts","Hmax +1 m","Hmax −1 m","profile (full swap)")
R<-rbindlist(lapply(LV,function(p) rbindlist(lapply(seq_along(VARS),function(i){
  v<-bmed(U[P==p][[VARS[i]]]); data.table(P=p,effect=LAB[i],median=v[1],lo=v[2],hi=v[3])}))))
fwrite(R,"out_files/Chapter1/tables/tab_attribution_units_bootstrap.csv")
cat("=== effet median (degC) et IC 95% bootstrap, par archetype ===\n")
for(p in LV){cat(sprintf("\n%s:\n",p))
  for(i in seq_along(VARS)){r<-R[P==p&effect==LAB[i]]
    cat(sprintf("  %-20s %+.3f  [%+.3f, %+.3f]\n",LAB[i],r$median,r$lo,r$hi))}}
# separation LAI vs profil (cote +)
cat("\n=== LAI(+) vs profil : IC disjoints ? ===\n")
for(p in LV){a<-R[P==p&effect=="LAI +0.5"];b<-R[P==p&effect=="profile (full swap)"]
  sep<-(a$hi<b$lo)||(b$hi<a$lo); cat(sprintf("  %s : %s\n",p,ifelse(sep,"SEPARES","chevauchement")))}
