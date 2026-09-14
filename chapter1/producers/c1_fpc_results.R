# ==============================================================================
# FPC-typology RESULTS (run once the 3200 sims finish): per-archetype per-SD
# sensitivity (ΔTmax) + bootstrap co-lead + figure, on the REAL balanced FPC
# cLHS sims (sensitivity_perplot_fpc). Mirrors the article deltadelta/bootstrap.
#   Rscript c1_fpc_results.R
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("scripts/_article_style.R") })
set.seed(1); NB <- 3000
M <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot_fpc","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
stopifnot(nrow(M) > 0); cat(sprintf("FPC plots: %d (par archétype: %s)\n", nrow(M), paste(table(M$P),collapse="/")))
M[, fCov_add:=ifelse(fCover>=0.95,0,fCov_add)][, fCov_rem:=ifelse(fCover<=0.5001,NA_real_,fCov_rem)]

# per +1 within-archetype SD (signed), incl LAD = real-vs-uniform per SD of (1-VCI)
sens <- M[, {
  sLAI<-sd(LAI/2,na.rm=T); sFC<-sd(fCover,na.rm=T)/0.1; sHM<-sd(Hmax,na.rm=T)/5; sVCI<-sd(VCI,na.rm=T)
  .(LAI=round(median((LAI_add-LAI_rem)/2*sLAI,na.rm=T),3),
    fCover=round(median((fCov_add-fCov_rem)/2*sFC,na.rm=T),3),
    Hmax=round(median((Hmax_add-Hmax_rem)/2*sHM,na.rm=T),3),
    LAD=round(median(dT_LAD/pmax(1-VCI,0.05)*sVCI,na.rm=T),3), n=.N)
}, by=P][order(P)]
cat("\n=== FPC typology — ΔTmax sensitivity per +1 SD (signed) ===\n"); print(sens)
fwrite(sens, "out_files/Chapter1/tables/tab_fpc_results_sensitivity.csv")

# bootstrap LAI - LAD importance difference per archetype
imp_of <- function(d){ c(LAI=median(abs((d$LAI_add-d$LAI_rem)/2),na.rm=T)*sd(d$LAI,na.rm=T)/2, LAD=abs(median(d$dT_LAD,na.rm=T))) }
ci <- function(v) c(med=median(v), lo=quantile(v,.025,names=F), hi=quantile(v,.975,names=F))
res <- rbindlist(lapply(sort(unique(M$P)), function(p){
  d<-M[P==p]; n<-nrow(d)
  B<-replicate(NB,{ i<-sample(n,n,TRUE); v<-imp_of(d[i]); unname(v["LAI"]-v["LAD"]) })
  c0<-ci(B); pt<-imp_of(d)
  data.table(P=p,n=n,impLAI=round(pt["LAI"],3),impLAD=round(pt["LAD"],3),
             LAI_minus_LAD=round(c0[1],3),lo=round(c0[2],3),hi=round(c0[3],3),straddle0=c0[2]<=0&c0[3]>=0)
}))
cat("\n=== bootstrap co-lead (LAI-LAD, 95% CI) ===\n"); print(res)
fwrite(res, "out_files/Chapter1/tables/tab_fpc_results_bootstrap.csv")
cat("\nFrozen ref: P1 LAI≈0.57→P4≈0.12, co-lead (straddle 0) ONLY in P4. FPC verdict above.\n")

# figure: per-SD sensitivity by archetype (Fig2-equivalent)
L <- melt(sens[, .(P,LAI,fCover,Hmax,LAD)], id.vars="P", variable.name="trait", value.name="dd")
p <- ggplot(L, aes(P, dd, colour=trait, group=trait)) +
  geom_hline(yintercept=0, colour="grey55", linewidth=0.4) +
  geom_line(linewidth=0.6) + geom_point(size=2.8) +
  scale_colour_manual(values=c(LAI="#1B9E77",fCover="#D95F02",Hmax="#7570B3",LAD="#E7298A"), name=NULL) +
  labs(x=NULL, y=expression(Delta*T[max]~"sensitivity per +1 SD ("*degree*"C)"),
       subtitle="FPC typology {LAI,Hmax,fCover,FPC1-3} — real balanced cLHS (v3.2.0/ERA5)") +
  theme_article(12) + theme(legend.position="bottom")
ggsave_article("out_files/Chapter1/figures/Fig_fpc_sensitivity_per_archetype", p, 7, 5)
cat("DONE -> tab_fpc_results_* + Fig_fpc_sensitivity_per_archetype.{png,pdf}\n")
