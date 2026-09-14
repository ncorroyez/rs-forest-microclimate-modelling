# ==============================================================================
# Additional LAD vertical-distribution descriptors (supervisor: "Lmom, Hmed etc,
# variabilité autour de la moyenne"). Treats each plot's LAD(z) profile as a leaf-
# area distribution over RELATIVE height (z/Hmax, pure shape) and computes:
#   Hmed_rel  : relative height holding 50% of leaf area (median leaf height)
#   Hmean_rel : leaf-area centroid (Σ z_rel·LAD / ΣLAD)
#   Lskew (t3): L-skewness of the leaf-height distribution (asymmetry; +=top-heavy)
#   Lcv       : L-CV = L-scale / L-mean (dispersion of leaf height)
# Per cluster: boxplots + mean table. Complements VCI by separating WHERE leaf
# area sits (Hmed) and HOW asymmetric/dispersed it is, and shows within-cluster
# variability. No MuSICA — reads the static cLHS sample.
#   Rscript c1_lmom_hmed.R
# Out: out_files/Chapter1/figures/FigSh_lmom_hmed.png + tab_lmom_hmed.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(lmom); source("R/cluster_relabel.R")
})
df <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
df <- df[is.finite(Hmax) & Hmax > 0]
lad_cols <- grep("^LAD_Layer_", names(df), value = TRUE)
z_real   <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

descr_one <- function(d, hmax) {
  d[is.na(d)] <- 0; if (sum(d) <= 0) return(c(Hmed_rel=NA, Hmean_rel=NA, Lskew=NA, Lcv=NA))
  zr <- z_real / hmax; keep <- zr <= 1.05 & d > 0; zr <- zr[keep]; w <- d[keep]
  if (length(zr) < 3) return(c(Hmed_rel=NA, Hmean_rel=NA, Lskew=NA, Lcv=NA))
  cum <- cumsum(w) / sum(w)
  Hmed <- approx(cum, zr, 0.5, rule = 2)$y
  Hmean <- sum(zr * w) / sum(w)
  samp <- rep(zr, times = pmax(1L, round(w / max(w) * 200)))   # weighted resample for L-moments
  lm <- tryCatch(samlmu(samp, nmom = 3), error = function(e) c(NA, NA, NA))
  Lcv <- if (is.finite(lm[1]) && lm[1] != 0) lm[2] / lm[1] else NA
  c(Hmed_rel = Hmed, Hmean_rel = Hmean, Lskew = unname(lm[3]), Lcv = Lcv)
}
Mat <- as.matrix(df[, ..lad_cols])
R <- as.data.table(t(vapply(seq_len(nrow(df)),
       function(i) descr_one(Mat[i, ], df$Hmax[i]),
       FUN.VALUE = c(Hmed_rel=0, Hmean_rel=0, Lskew=0, Lcv=0))))
R[, Cluster := df$Cluster]; R[, P := relabel_cluster(Cluster)]

summ <- R[, .(Hmed_rel = round(mean(Hmed_rel, na.rm=TRUE),3), Hmean_rel = round(mean(Hmean_rel, na.rm=TRUE),3),
              Lskew = round(mean(Lskew, na.rm=TRUE),3), Lcv = round(mean(Lcv, na.rm=TRUE),3)), by = P][order(P)]
cat("=== mean LAD distribution descriptors per cluster ===\n"); print(summ)
fwrite(summ, "out_files/Chapter1/tables/tab_lmom_hmed.csv")

LAB <- c(Hmed_rel="Median~leaf~height~(z/H[max])", Hmean_rel="Leaf-area~centroid~(z/H[max])",
         Lskew="L-skewness~(asymmetry)", Lcv="L-CV~(dispersion)")
ml <- melt(R, id.vars = "P", measure.vars = names(LAB), variable.name = "metric", value.name = "val")
ml[, metric := factor(metric, levels = names(LAB), labels = LAB)]
ml <- ml[is.finite(val)]
set.seed(1)
p <- ggplot(ml, aes(P, val, fill = P, colour = P)) +
  geom_boxplot(width = 0.6, alpha = 0.35, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.12, size = 0.7, alpha = 0.35) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1, labeller = label_parsed) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  scale_colour_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"))
ggsave("out_files/Chapter1/figures/FigSh_lmom_hmed.png", p, width = 12, height = 4.2, dpi = 200, bg = "white")
cat("DONE -> FigSh_lmom_hmed.png + tab_lmom_hmed.csv\n")
