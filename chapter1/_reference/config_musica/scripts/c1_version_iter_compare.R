# ==============================================================================
# 3-way per-unit sensitivity importance comparison: legacy v3.2.0, v3.2.3 (none),
# v3.2.3 (iter, real MERRA-2 PBLH). Same importance formula as Fig 2
# (c1_importance_sensitivity.R) for all three. Answers: does the yoyo (iter) change
# the trait ranking / density-dependent pattern vs none?
#   Rscript c1_version_iter_compare.R
# Out: tab_version_iter_compare.csv + FigSh_version_iter_compare.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })

imp_from <- function(M) {
  M[, P := relabel_cluster(Cluster)]
  blk <- function(Ss) {
    a <- Ss[, .(s_LAI=median(abs((LAI_add-LAI_rem)/2),na.rm=TRUE),
                s_fC =median(abs((fCov_add-fCov_rem)/2),na.rm=TRUE)/0.1,
                s_Hm =median(abs((Hmax_add-Hmax_rem)/2),na.rm=TRUE)/5,
                sd_LAI=sd(LAI,na.rm=TRUE)/2, sd_fC=sd(fCover,na.rm=TRUE), sd_Hm=sd(Hmax,na.rm=TRUE),
                LADc=abs(median(dT_LAD,na.rm=TRUE)))]
    data.table(LAI=a$s_LAI*a$sd_LAI, fCover=a$s_fC*a$sd_fC, Hmax=a$s_Hm*a$sd_Hm, LAD=a$LADc)
  }
  R <- rbind(M[, blk(.SD), by=P], cbind(P="All", blk(M)))
  R[, c("LAI","fCover","Hmax","LAD") := lapply(.SD, round, 3), .SDcols=c("LAI","fCover","Hmax","LAD")]
  R[order(P)]
}
load_parts <- function(dir) rbindlist(lapply(list.files(dir,"part_.*csv$",full.names=TRUE),fread),fill=TRUE)

# v3.2.0 (legacy): per-unit parts + dT_LAD from separate table
L <- load_parts("out_files/Chapter3/tables/sensitivity_perplot")
L <- merge(L, fread("out_files/Chapter3/tables/tab_sensitivity_perplot_lad.csv")[,.(pid,dT_LAD)], by="pid", all.x=TRUE)
imp <- rbind(
  imp_from(L)[, version:="v3.2.0"],
  imp_from(load_parts("out_files/Chapter3/tables/sensitivity_perplot_v323"))[, version:="v3.2.3 none"],
  imp_from(load_parts("out_files/Chapter3/tables/sensitivity_perplot_iter"))[, version:="v3.2.3 iter (realBLH)"]
)
fwrite(imp, "out_files/Chapter3/tables/tab_version_iter_compare.csv")
cat("=== per-unit importance (°C) by version ===\n"); print(imp)
cat("\n=== rank order per cluster ===\n")
for (cl in c("P1","P2","P3","P4","All")) for (v in unique(imp$version)) {
  o <- names(sort(unlist(imp[P==cl & version==v, .(LAI,fCover,Hmax,LAD)]), decreasing=TRUE))
  cat(sprintf("  %-3s %-22s : %s\n", cl, v, paste(o,collapse=">")))
}
# does iter differ from none? (per-cluster max abs importance diff)
n <- imp[version=="v3.2.3 none"]; it <- imp[version=="v3.2.3 iter (realBLH)"]
m <- merge(n, it, by="P", suffixes=c(".none",".iter"))
m[, maxdiff := pmax(abs(LAI.none-LAI.iter),abs(fCover.none-fCover.iter),abs(Hmax.none-Hmax.iter),abs(LAD.none-LAD.iter))]
cat("\n=== max |importance(iter) - importance(none)| per cluster (°C) ===\n"); print(m[,.(P,maxdiff=round(maxdiff,3))])

L2 <- melt(imp[P!="All"], id.vars=c("P","version"), variable.name="trait", value.name="imp")
L2[, trait := factor(trait, levels=c("LAI","fCover","LAD","Hmax"))]
L2[, version := factor(version, levels=c("v3.2.0","v3.2.3 none","v3.2.3 iter (realBLH)"))]
PAL <- c(LAI="#1B7837", fCover="#7FBC41", LAD="#762A83", Hmax="#BDBDBD")
p <- ggplot(L2, aes(P, imp, fill=trait)) +
  geom_col(position=position_dodge(0.8), width=0.74) + facet_wrap(~version, nrow=1) +
  scale_fill_manual(values=PAL, name=NULL) + scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL, y=expression(Delta*T[max]~importance~(degree*C)),
       title="Per-unit sensitivity importance: legacy v3.2.0 vs v3.2.3 (none) vs v3.2.3 iter/yoyo (real MERRA-2 PBLH)",
       subtitle="Does turning on the ABL yoyo (iter) with real boundary-layer height change which traits matter? Compare middle vs right panel.") +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), legend.position="bottom",
        strip.text=element_text(face="bold"), plot.title=element_text(size=10,face="bold"),
        plot.subtitle=element_text(size=8,colour="grey35"))
ggsave("out_files/Chapter3/figures/FigSh_version_iter_compare.png", p, width=12, height=5, dpi=200, bg="white")
cat("\nDONE -> tab_version_iter_compare.csv + FigSh_version_iter_compare.png\n")
