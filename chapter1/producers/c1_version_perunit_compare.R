# ==============================================================================
# Compare the ARTICLE per-unit sensitivity importance (ΔLAI/ΔfCover/ΔHmax + LAD
# real-vs-uniform) between legacy v3.2.0 and v3.2.3, on the SAME 400 cLHS plots.
# Importance computed IDENTICALLY to c1_importance_sensitivity.R for both versions:
#   scalar = |median per-unit slope| × within-cluster SD (one-sided LAI)
#   LAD    = |median real-vs-uniform contrast|
# Legacy reads sensitivity_perplot/part_*.csv + tab_sensitivity_perplot_lad.csv;
# v3.2.3 reads sensitivity_perplot_v323/part_*.csv (has dT_LAD column inline).
#   Rscript c1_version_perunit_compare.R
# Out: tab_version_perunit_compare.csv + FigSh_version_perunit_compare.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })

imp_from <- function(M) {   # M has per-unit cols AND a dT_LAD column
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

# legacy: per-unit parts + dT_LAD from the separate LAD table
Sleg <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
Llad <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")[,.(pid,dT_LAD)]
Mleg <- merge(Sleg, Llad, by="pid", all.x=TRUE)
imp_leg <- imp_from(Mleg)[, version:="v3.2.0"]
# v323: dT_LAD already inline in the parts
Sv3 <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot_v323","part_.*csv$",full.names=TRUE),fread),fill=TRUE)
imp_v3 <- imp_from(Sv3)[, version:="v3.2.3"]

OUT <- rbind(imp_leg, imp_v3)
fwrite(OUT, "out_files/Chapter1/tables/tab_version_perunit_compare.csv")
cat("=== per-unit sensitivity importance, legacy vs v3.2.3 (°C) ===\n"); print(OUT)

# ranking order per cluster per version
cat("\n=== rank order per cluster ===\n")
for (cl in c("P1","P2","P3","P4","All")) {
  o_l <- names(sort(unlist(imp_leg[P==cl, .(LAI,fCover,Hmax,LAD)]), decreasing=TRUE))
  o_v <- names(sort(unlist(imp_v3 [P==cl, .(LAI,fCover,Hmax,LAD)]), decreasing=TRUE))
  cat(sprintf("  %-3s  v3.2.0: %-22s | v3.2.3: %-22s | %s\n", cl,
              paste(o_l,collapse=">"), paste(o_v,collapse=">"),
              if(identical(o_l,o_v)) "MATCH" else "DIFFER"))
}

# figure: per-cluster grouped bars, facet by version
L <- melt(OUT[P!="All"], id.vars=c("P","version"), variable.name="trait", value.name="imp")
L[, trait := factor(trait, levels=c("LAI","fCover","LAD","Hmax"))]
PAL <- c(LAI="#1B7837", fCover="#7FBC41", LAD="#762A83", Hmax="#BDBDBD")
p <- ggplot(L, aes(P, imp, fill=trait)) +
  geom_col(position=position_dodge(0.8), width=0.74) +
  facet_wrap(~version, nrow=1) +
  scale_fill_manual(values=PAL, name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x=NULL, y=expression(Delta*T[max]~importance~(degree*C)),
       title="Per-unit sensitivity: the density-dependent SHIFT is robust to version, its MAGNITUDE is not (400 cLHS)",
       subtitle=paste0("Pooled trait order identical (LAI>fCover>LAD>Hmax) and in both, LAI's role erodes while LAD's rises toward the dense end. ",
                       "But v3.2.3\namplifies the vertical effect — LAD LEADS LAI in P3-P4 (and the open P1 collapses) — whereas the validated v3.2.0 gives the ",
                       "conservative co-lead.\nv3.2.3 validates worse (between-plot r 0.93 -> 0.63), so the field does not support its stronger vertical role.")) +
  theme_bw(base_size=12) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position="bottom", strip.text=element_text(face="bold"),
        plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigSh_version_perunit_compare.png", p, width=10, height=5, dpi=200, bg="white")
cat("\nDONE -> tab_version_perunit_compare.csv + FigSh_version_perunit_compare.png\n")
