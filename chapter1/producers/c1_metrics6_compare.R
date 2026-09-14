# ==============================================================================
# Compare v3.2.0 vs v3.2.3-iter ACROSS the 6 metric×period combos (ΔTmax, micro/
# macro slope, ΔVPDmax × all-days / hot-10%). Per-unit importance (same formula as
# Fig 2) per metric × period × cluster × version → heatmap of the leading trait,
# annotated with LAD's relative share. Inputs: metrics6_{v320,iter}/part_*.csv.
#   Rscript c1_metrics6_compare.R
# Out: tab_metrics6_compare.csv + FigCmp_metrics6_v320_v323iter.png
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })

imp_block <- function(M) {   # M: rows for one (ver,metric); has P + per-unit + LAI/Hmax/fCover + dT_LAD
  blk <- function(S) {
    a <- S[, .(s_LAI=median(abs((LAI_add-LAI_rem)/2),na.rm=TRUE),
               s_fC =median(abs((fCov_add-fCov_rem)/2),na.rm=TRUE)/0.1,
               s_Hm =median(abs((Hmax_add-Hmax_rem)/2),na.rm=TRUE)/5,
               sd_LAI=sd(LAI,na.rm=TRUE)/2, sd_fC=sd(fCover,na.rm=TRUE), sd_Hm=sd(Hmax,na.rm=TRUE),
               LADc=abs(median(dT_LAD,na.rm=TRUE)))]
    data.table(LAI=a$s_LAI*a$sd_LAI, fCover=a$s_fC*a$sd_fC, Hmax=a$s_Hm*a$sd_Hm, LAD=a$LADc)
  }
  M[, blk(.SD), by=P]
}
load_ver <- function(ver) rbindlist(lapply(list.files(sprintf("out_files/Chapter1/tables/metrics6_%s",ver),
                                                       "part_.*csv$", full.names=TRUE), fread), fill=TRUE)
D <- rbind(load_ver("v320"), load_ver("iter"))
D[, P := relabel_cluster(Cluster)]
VLAB <- c(v320="v3.2.0", iter="v3.2.3 iter")

# importance per (ver, metric, cluster)
imp <- D[, imp_block(.SD), by=.(ver, metric)]
imp[, tot := LAI+fCover+Hmax+LAD]
imp[, `:=`(lead = c("LAI","fCover","Hmax","LAD")[max.col(.SD)], LAD_share = LAD/tot),
    .SDcols=c("LAI","fCover","Hmax","LAD")]
fwrite(imp, "out_files/Chapter1/tables/tab_metrics6_compare.csv")
cat("=== leading trait & LAD share per metric×period×cluster×version ===\n")
print(imp[, .(ver, metric, P, lead, LAD_share=round(LAD_share,2),
              LAI=round(LAI,3), LAD=round(LAD,3))][order(ver, metric, P)])

# ---- heatmap: leading trait (fill) + LAD share (label), rows=metric×period ----
MLAB <- c(Tmax_all="ΔTmax · all", Tmax_hot="ΔTmax · hot10%",
          slope_all="slope · all", slope_hot="slope · hot10%",
          VPD_all="ΔVPDmax · all", VPD_hot="ΔVPDmax · hot10%")
imp[, metric := factor(metric, levels=rev(names(MLAB)), labels=rev(MLAB))]
imp[, verlab := factor(VLAB[ver], levels=VLAB)]
imp[, lead := factor(lead, levels=c("LAI","fCover","LAD","Hmax"))]
PAL <- c(LAI="#1B7837", fCover="#7FBC41", LAD="#762A83", Hmax="#BDBDBD")
p <- ggplot(imp, aes(P, metric, fill=lead)) +
  geom_tile(colour="white", linewidth=1.2) +
  geom_text(aes(label=sprintf("%s\nLAD %.0f%%", lead, 100*LAD_share)), size=2.6,
            colour=ifelse(imp$lead=="LAD","white","grey15")) +
  facet_wrap(~verlab, nrow=1) +
  scale_fill_manual(values=PAL, name="trait dominant", drop=FALSE) +
  labs(x=NULL, y=NULL,
       title="Trait dominant par métrique × période × archétype : v3.2.0 vs v3.2.3 iter (yoyo)",
       subtitle="Importance per-unit (|sensibilité médiane| × dispersion intra-cluster ; LAD = |réel-vs-uniforme|). Couleur = trait dominant ; étiquette = part du LAD.") +
  theme_bw(base_size=12) +
  theme(panel.grid=element_blank(), strip.text=element_text(face="bold"),
        plot.title=element_text(size=11,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"),
        legend.position="bottom")
ggsave("out_files/Chapter1/figures/FigCmp_metrics6_v320_v323iter.png", p, width=10, height=6, dpi=200, bg="white")
cat("DONE -> tab_metrics6_compare.csv + FigCmp_metrics6_v320_v323iter.png\n")
