# ==============================================================================
# TEST: make_lad_real vs make_lad_uniform — visualise the exact contrast the sim
# sees for dT_LAD = base − unifLAD. One representative plot per archetype (the
# plot closest to the cluster median dT_LAD). Overlays density(height) real vs
# uniform, same total LAI. Diagnoses whether "profile anti-buffering" is a
# top-heavy-vs-uniform shape effect or something else.
#   → out_files/Chapter1/figures/Fig_test_lad_real_vs_uniform.{png,pdf}
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2)
  source("R/cluster_relabel.R"); source("scripts/_article_style.R"); source("R/lad.R") })
s <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
s[, pid := sprintf("S%04d", .I)]; s[, P := as.character(relabel_cluster(Cluster))]
m <- rbindlist(lapply(list.files("out_files/Chapter1/tables/metrics6_frblo","part.*csv",full.names=T),fread),fill=T)
m <- m[metric=="Tmax_all",.(pid,dT_LAD)]
s <- merge(s, m, by="pid")
# representative plot per cluster = closest to cluster-median dT_LAD
reps <- s[, .SD[which.min(abs(dT_LAD - median(dT_LAD)))], by=P]

prof <- rbindlist(lapply(seq_len(nrow(reps)), function(i){
  pr <- as.data.frame(reps[i])
  lr <- make_lad_real(pr);  lr$type <- "real (LiDAR)"
  lu <- make_lad_uniform(pr); lu$type <- "uniform"
  d  <- rbind(lr, lu); d$P <- reps$P[i]
  d$lab <- sprintf("%s  (LAI=%.1f, Hmax=%.0fm, dT_LAD=%+.2f)", reps$P[i], reps$LAI[i]/2, reps$Hmax[i], reps$dT_LAD[i])
  as.data.table(d) }))
prof[, P := factor(P, levels=c("P1","P2","P3","P4"))]
setorder(prof, P); prof[, lab := factor(lab, levels=unique(lab))]

p <- ggplot(prof, aes(density, height, colour=type)) +
  geom_hline(yintercept=1, linetype=3, colour="grey70") +               # understory sensor 1 m
  geom_path(linewidth=0.9) +
  facet_wrap(~lab, nrow=1, scales="free") +
  scale_colour_manual(values=c("real (LiDAR)"="#1A9850","uniform"="#7570B3"), name=NULL) +
  labs(x=expression("LAD (m"^2*" m"^{-3}*")"), y="Height (m)",
       subtitle="make_lad_real vs make_lad_uniform (same total LAI). Dotted = 1 m understory. Positive dT_LAD = real warms vs uniform.") +
  theme_article(11) + theme(legend.position="bottom")
ggsave_article("out_files/Chapter1/figures/Fig_test_lad_real_vs_uniform", p, 11, 4.8)

# quantitative: leaf fraction in the LOW canopy (<1/3 Hmax) real vs uniform
cat("=== leaf-area fraction below 1/3 Hmax (low canopy) — more low = more understory shading ===\n")
low<-prof[, {h3<-max(height)/3; .(frac_low=sum(density[height<=h3])/sum(density))}, by=.(P,type)]
print(dcast(low,P~type,value.var="frac_low"))
cat("DONE -> Fig_test_lad_real_vs_uniform\n")
