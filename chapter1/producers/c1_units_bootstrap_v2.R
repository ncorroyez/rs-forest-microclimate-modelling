# ==============================================================================
# Bootstrap of the native-unit sensitivities (Table G1)
# Bootstrap confidence intervals on the per-archetype median trait effects of the
# native-unit design, for BOTH metrics and BOTH periods (sensitivity_perplot_units_v2).
# Resamples plots within each archetype; effects are the change ACTUALLY simulated
# (truncated steps kept as run, never rescaled to the nominal step).
# Reads :
#         out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv  (stage A3)
# Out: out_files/Chapter1/tables/tab_attribution_units_bootstrap_v2.csv
# ==============================================================================
suppressPackageStartupMessages({library(data.table);source("R/cluster_relabel.R")})
X<-fread("out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv")
X[,P:=factor(relabel_cluster(Cluster),levels=c("P1","P2","P3","P4"))]
X[,`:=`(sLp=0.5, sLm=LAI-pmax(LAI-0.5,0.1),
        sHp=1.0, sHm=Hmax-pmax(Hmax-1,3),
        sFp=pmin(fCover+0.10,1)-fCover, sFm=fCover-pmax(fCover-0.10,0.5))]
set.seed(1); NB<-4000   # matches Table G1 (c1_attribution_bootstrap_units.R), whose caption states 4000
#' Bootstrap median and 95% percentile interval of one effect vector
#' @param v numeric vector of per-plot effects, non-finite values dropped
#' @return list(median, lo, hi, n); all NA with n = 0 when fewer than 5 finite values
bs<-function(v){v<-v[is.finite(v)]; if(length(v)<5) return(list(NA_real_,NA_real_,NA_real_,0L))
  m<-replicate(NB,median(sample(v,length(v),TRUE)))
  list(median(v),unname(quantile(m,0.025)),unname(quantile(m,0.975)),length(v))}
EFF<-list(
  "LAI +0.5"      =function(D,mt) D[[sprintf("LAI_up_%s",mt)]] * D$sLp,
  "LAI -0.5"      =function(D,mt) D[[sprintf("LAI_dn_%s",mt)]] * (-D$sLm),
  "fCover +10 pts"=function(D,mt) D[[sprintf("fCov_up_%s",mt)]]* D$sFp,
  "fCover -10 pts"=function(D,mt) D[[sprintf("fCov_dn_%s",mt)]]* (-D$sFm),
  "Hmax +1 m"     =function(D,mt) D[[sprintf("Hmax_up_%s",mt)]]* D$sHp,
  "Hmax -1 m"     =function(D,mt) D[[sprintf("Hmax_dn_%s",mt)]]* (-D$sHm),
  "profile (full swap)"=function(D,mt) D[[sprintf("LAD_%s",mt)]])
out<-rbindlist(lapply(c("dt_all","dt_hot","sl_all","sl_hot"),function(mt)
  rbindlist(lapply(levels(X$P),function(p){ D<-X[P==p]
    rbindlist(lapply(names(EFF),function(e){ r<-bs(EFF[[e]](D,mt))
      data.table(metric=sub("_.*","",mt), period=sub(".*_","",mt), P=p, effect=e,
                 median=r[[1]], lo=r[[2]], hi=r[[3]], n=r[[4]]) }))}))))
fwrite(out,"out_files/Chapter1/tables/tab_attribution_units_bootstrap_v2.csv")
cat("=== pente, ete complet : les trois leviers sont-ils separables ? ===\n")
print(out[metric=="sl"&period=="all"&effect%in%c("LAI +0.5","fCover +10 pts","profile (full swap)"),
          .(P,effect,median=round(median,4),lo=round(lo,4),hi=round(hi,4))])
cat(sprintf("\n%d lignes ecrites (4 metriques/periodes x 4 archetypes x %d effets)\nDONE\n",nrow(out),length(EFF)))
