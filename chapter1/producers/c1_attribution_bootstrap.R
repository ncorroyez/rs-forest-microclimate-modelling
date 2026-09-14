# Audit 2.1: bootstrap 95% CI on the per-archetype median trait sensitivity (ΔTmax),
# to test whether LAI's lever is significantly separated from the vertical profile.
# Replicates the per-plot sensitivity of c1_deltadelta_perplot_frblo.R (NORM=sd, ΔTmax).
suppressPackageStartupMessages({library(data.table);source("R/cluster_relabel.R")})
LV<-c("P1","P2","P3","P4")
VER<-Sys.getenv("VER","frblo")
M<-rbindlist(lapply(list.files(sprintf("out_files/Chapter1/tables/metrics6_%s",VER),"part_.*csv$",full.names=TRUE),fread),fill=TRUE)
M<-M[metric=="Tmax_all"]; M[,P:=factor(relabel_cluster(Cluster),levels=LV)]
.clhs<-if(VER=="native20")"out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds" else "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"
vci<-as.data.table(readRDS(.clhs))[,.(pid=sprintf("S%04d",.I),VCI)]
M<-merge(M,vci,by="pid",all.x=TRUE)
M[,`:=`(LAI_s=(LAI_add-LAI_rem)/2, fCover_s=(fCov_add-fCov_rem)/2, Hmax_s=(Hmax_add-Hmax_rem)/2, LAD_s=dT_LAD)]
# LAD_s stays = dT_LAD (raw real-vs-uniform contrast, full swap / upper bound); NO (1-VCI) rescaling
# (committee 2026-07: a VCI maps to infinitely many profiles; the swap is what we actually compute).
vars<-c(LAI="LAI_s",fCover="fCover_s",Hmax="Hmax_s",LAD="LAD_s")
set.seed(1); B<-4000
boot_med<-function(x){x<-x[is.finite(x)];if(length(x)<3)return(c(NA,NA,NA));qs<-quantile(replicate(B,median(sample(x,replace=TRUE))),c(.025,.5,.975));c(qs)}
cat("=== Per-archetype median ΔTmax sensitivity (°C per +1 SD) with 95% bootstrap CI ===\n")
res<-rbindlist(lapply(LV,function(p)rbindlist(lapply(names(vars),function(v){
  ci<-boot_med(M[P==p][[vars[v]]]); data.table(P=p,trait=v,lo=ci[1],med=ci[2],hi=ci[3])}))))
res[,`:=`(lo=round(lo,3),med=round(med,3),hi=round(hi,3))]
for(p in LV){cat(sprintf("\n%s:\n",p));d<-res[P==p];for(i in 1:nrow(d))cat(sprintf("  %-7s median %+.3f  [%+.3f, %+.3f]\n",d$trait[i],d$med[i],d$lo[i],d$hi[i]))
  lai<-d[trait=="LAI"];lad<-d[trait=="LAD"]
  sep<-if(lai$hi<lad$lo|lad$hi<lai$lo)"SEPARATED (non-overlapping CIs)" else "overlap"
  cat(sprintf("  -> LAI vs LAD 95%% CIs: %s\n",sep))}
# per-plot paired dominance: fraction of plots where |LAI| > |LAD|
cat("\n=== per-plot paired: fraction with |LAI sensitivity| > |profile| ===\n")
for(p in LV){m<-M[P==p];fr<-mean(abs(m$LAI_s)>abs(m$LAD_s),na.rm=TRUE);cat(sprintf("  %s: %.0f%% (n=%d)\n",p,100*fr,nrow(m)))}
# RAW dT_LAD (real-minus-uniform ΔTmax, VCI-free) — protects the "small lever" claim from any VCI-normalisation artifact
cat("\n=== RAW dT_LAD per archetype (°C, real-vs-uniform, NO VCI normalisation) ===\n")
for(p in LV){x<-M[P==p]$dT_LAD;x<-x[is.finite(x)];cat(sprintf("  %s: median %+.3f  mean %+.3f  [%+.3f, %+.3f]  n=%d\n",p,median(x),mean(x),quantile(x,.05),quantile(x,.95),length(x)))}
fwrite(res, if(VER=="native20") "out_files/Chapter1/tables/tab_attribution_bootstrap_ci_native20.csv" else "out_files/Chapter1/tables/tab_attribution_bootstrap_ci.csv")
cat("DONE\n")
