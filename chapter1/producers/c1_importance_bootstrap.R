# ==============================================================================
# Test du co-lead (durcissement de la Fig 2) : bootstrap des DIFFÉRENCES d'importance
# par archétype (rééch. des plots intra-cluster, 3000×, binaire legacy v3.2.0).
#  - imp_LAI − imp_LAD  (P4 : straddle 0 → co-lead/égalité réelle)
#  - imp_LAD − imp_fCover, imp_LAD − imp_Hmax  (P3/P4 : IC>0 → le LAD passe #2)
# Même formule d'importance que c1_importance_sensitivity.R. NE modifie PAS la Fig 2.
#   Rscript c1_importance_bootstrap.R
# Out: tab_importance_bootstrap.csv + FigSh_importance_bootstrap_colead.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R"); source("scripts/_article_style.R") })
set.seed(1); NB <- 3000
S <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
LAD <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")[,.(pid,dT_LAD)]
S <- merge(S, LAD, by="pid", all.x=TRUE); S[, P := relabel_cluster(Cluster)]

imp_of <- function(d) {  # importances d'un sous-échantillon (mêmes formules que Fig 2)
  s_LAI <- median(abs((d$LAI_add - d$LAI_rem)/2), na.rm=TRUE)
  s_fC  <- median(abs((d$fCov_add - d$fCov_rem)/2), na.rm=TRUE)/0.1
  s_Hm  <- median(abs((d$Hmax_add - d$Hmax_rem)/2), na.rm=TRUE)/5
  c(LAI = s_LAI*sd(d$LAI,na.rm=TRUE)/2, fCover = s_fC*sd(d$fCover,na.rm=TRUE),
    Hmax = s_Hm*sd(d$Hmax,na.rm=TRUE), LAD = abs(median(d$dT_LAD,na.rm=TRUE)))
}
ci <- function(v) c(est=median(v), lo=quantile(v,.025,names=FALSE), hi=quantile(v,.975,names=FALSE))

res <- rbindlist(lapply(c("P1","P2","P3","P4"), function(p){
  d <- S[P==p]; n <- nrow(d); pt <- imp_of(d)
  B <- t(replicate(NB, { i<-sample(n,n,TRUE); imp_of(d[i]) }))
  dLL <- B[,"LAI"]-B[,"LAD"]; dLf <- B[,"LAD"]-B[,"fCover"]; dLh <- B[,"LAD"]-B[,"Hmax"]
  data.table(P=p, n=n,
    impLAI=round(pt["LAI"],3), impLAD=round(pt["LAD"],3), impfC=round(pt["fCover"],3), impHm=round(pt["Hmax"],3),
    LAI_minus_LAD=round(ci(dLL)[1],3), LL_lo=round(ci(dLL)[2],3), LL_hi=round(ci(dLL)[3],3),
    LAD_minus_fCover=round(ci(dLf)[1],3), Lf_lo=round(ci(dLf)[2],3), Lf_hi=round(ci(dLf)[3],3), Lf_p=round(mean(dLf<=0),3),
    LAD_minus_Hmax=round(ci(dLh)[1],3), Lh_lo=round(ci(dLh)[2],3), Lh_hi=round(ci(dLh)[3],3), Lh_p=round(mean(dLh<=0),3))
}))
fwrite(res, "out_files/Chapter1/tables/tab_importance_bootstrap.csv")
cat("=== bootstrap des différences d'importance (IC 95%) ===\n"); print(res)
cat("\nLecture : LAI−LAD straddle 0 en P4 ? (co-lead) ; LAD−fCover et LAD−Hmax >0 en P3/P4 ? (LAD #2)\n")

# figure : différences clés avec IC (P3, P4)
G <- rbind(
  res[, .(P, comp="LAI − LAD", est=LAI_minus_LAD, lo=LL_lo, hi=LL_hi)],
  res[, .(P, comp="LAD − fCover", est=LAD_minus_fCover, lo=Lf_lo, hi=Lf_hi)],
  res[, .(P, comp="LAD − Hmax", est=LAD_minus_Hmax, lo=Lh_lo, hi=Lh_hi)])
G[, comp:=factor(comp, levels=c("LAI − LAD","LAD − fCover","LAD − Hmax"))]
G[, P := factor(P, levels=c("P4","P3","P2","P1"))]   # P1 top -> P4 bottom
p <- ggplot(G, aes(est, P, colour=comp)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.22, position=position_dodge(.55), linewidth=0.5) +
  geom_point(size=2.4, position=position_dodge(.55)) +
  scale_colour_manual(values=c("LAI − LAD"="#1B9E77","LAD − fCover"="#D95F02","LAD − Hmax"="#7570B3"), name=NULL) +
  labs(x="Importance difference (°C) ± bootstrap 95% CI", y=NULL,
       subtitle="LAI − LAD straddles 0 in dense P4 = real co-lead; LAD − fCover / LAD − Hmax > 0 in P3/P4 = profile overtakes cover & height.") +
  theme_article(11) + theme(legend.position="bottom")
ggsave_article("out_files/Chapter1/figures/FigSh_importance_bootstrap_colead", p, 8.5, 4.5)
cat("DONE -> tab_importance_bootstrap.csv + FigSh_importance_bootstrap_colead (png+pdf)\n")
