# ==============================================================================
# #3 — LAI vs fCover separability. They are collinear (r=0.76), so the FIELD
# cannot tell them apart (#2: LMG fCover 0.36 ≈ LAI 0.28, RF ties them). The MODEL
# can: the on-manifold sensitivity perturbs each ORTHOGONALLY (LAI±1 at fixed
# fCover, fCover±0.1 at fixed LAI). This script bootstraps the model-side
# importance DIFFERENCE LAI−fCover (and LAI−Hmax) per archetype to test whether,
# once separation is enforced, the LAI lead over cover is statistically resolved.
# Same importance formula as Fig 2 / c1_importance_bootstrap.R.
#   Rscript c1_lai_fcover_separability.R
# Out: tab_lai_fcover_separability.csv + FigSh_lai_fcover_separability.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
set.seed(1); NB <- 3000
S <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
LAD <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")[,.(pid,dT_LAD)]
S <- merge(S, LAD, by="pid", all.x=TRUE); S[, P := relabel_cluster(Cluster)]

imp_of <- function(d) {  # identical to Fig 2 importances
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
  dLf <- B[,"LAI"]-B[,"fCover"]; dLh <- B[,"LAI"]-B[,"Hmax"]
  data.table(P=p, n=n,
    impLAI=round(pt["LAI"],3), impfCover=round(pt["fCover"],3), impHmax=round(pt["Hmax"],3),
    LAI_minus_fCover=round(ci(dLf)[1],3), Lf_lo=round(ci(dLf)[2],3), Lf_hi=round(ci(dLf)[3],3),
    Lf_p_le0=round(mean(dLf<=0),3),
    LAI_minus_Hmax=round(ci(dLh)[1],3), Lh_lo=round(ci(dLh)[2],3), Lh_hi=round(ci(dLh)[3],3),
    Lh_p_le0=round(mean(dLh<=0),3))
}))
fwrite(res, "out_files/Chapter1/tables/tab_lai_fcover_separability.csv")
cat("\n=== #3 LAI vs fCover separability (MODEL, orthogonal perturbation, bootstrap) ===\n")
print(res)
cat("\nLAI−fCover > 0 with CI excluding 0 => model resolves LAI ahead of cover in that archetype.\n")
cat("p_le0 = bootstrap share of replicates where LAI <= fCover (one-sided).\n\n")
for (p in c("P1","P2","P3","P4")) {
  r <- res[P==p]
  verdict <- if (r$Lf_lo > 0) "LAI > fCover (resolved)" else if (r$Lf_hi < 0) "fCover > LAI (resolved)" else "LAI ≈ fCover (not separable even in model)"
  cat(sprintf("  %s: LAI−fCover = %+.3f [%.3f, %.3f], p(LAI<=fC)=%.3f -> %s\n",
              p, r$LAI_minus_fCover, r$Lf_lo, r$Lf_hi, r$Lf_p_le0, verdict))
}

# figure
G <- rbind(
  res[, .(P, comp="LAI − fCover", est=LAI_minus_fCover, lo=Lf_lo, hi=Lf_hi)],
  res[, .(P, comp="LAI − Hmax",   est=LAI_minus_Hmax,   lo=Lh_lo, hi=Lh_hi)])
G[, comp:=factor(comp, levels=c("LAI − fCover","LAI − Hmax"))]
p <- ggplot(G, aes(est, P, colour=comp)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=.22, position=position_dodge(.55)) +
  geom_point(size=2.2, position=position_dodge(.55)) +
  scale_colour_manual(values=c("LAI − fCover"="#1B7837","LAI − Hmax"="#BDBDBD"), name=NULL) +
  labs(x="Model importance difference (°C) ± bootstrap 95% CI", y=NULL,
       title="Model's ceteris-paribus attribution favours LAI over cover/height in every archetype",
       subtitle="Orthogonal on-manifold perturbation (robust to normalization). NB: observed association favours fCover (#2) — the traits co-vary (r=0.76); which facet of quantity 'leads' is method-dependent.") +
  theme_bw(base_size=12)+theme(legend.position="bottom",panel.grid.minor=element_blank(),
    plot.title=element_text(size=11,face="bold"),plot.subtitle=element_text(size=8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigSh_lai_fcover_separability.png", p, width=8.5,height=4.2,dpi=200,bg="white")
cat("\nDONE -> tab_lai_fcover_separability.csv + FigSh_lai_fcover_separability.png\n")
