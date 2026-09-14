# ==============================================================================
# Per-cluster sensitivity figures from tab_sensitivity_percluster.csv.
#
# HEADLINE (local) : metric response to a LOCAL change of each variable at the
#   cluster's own operating point (within-cluster p10-p90 finite difference),
#   normalized to the variable's full-sample range. Shows SATURATION: ΔTmax
#   sensitivity to LAI collapses from open (P1) to dense (P4) canopies. This is
#   the physical precondition for LAD to LEAD in the per-cluster Shapley.
# COMPLEMENT (global) : intrinsic potency = response across each variable's full
#   range (others at cluster mean). Cluster-invariant: LAI is the most potent
#   lever everywhere; only its within-cluster room-to-move is local — which is
#   why per-cluster Shapley reorders without LAI losing its physical power.
#   Rscript c1_sensitivity_figure.R
# Out: out_files/Chapter1/figures/FigSh_sensitivity_local.png
#      out_files/Chapter1/figures/FigSh_sensitivity_potency.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R")
})
s <- fread("out_files/Chapter1/tables/tab_sensitivity_percluster.csv")
s[, P := relabel_cluster(Cluster)]
s[, variable := factor(variable, levels = c("LAI","fCover","Hmax","LAD"))]

LAB <- c(Tmax_all = "Delta*T[max]~('°C')",
         VPD_all  = "VPD[max]~('kPa')",
         slope_all = "Micro/macro~slope")

mk_fig <- function(scope_sel, ylab, out, w = 11) {
  d <- melt(s[scope == scope_sel], id.vars = c("P","variable"),
            measure.vars = names(LAB), variable.name = "metric", value.name = "val")
  d[, metric := factor(metric, levels = names(LAB), labels = LAB)]
  p <- ggplot(d, aes(variable, val, fill = P)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_col(position = position_dodge(width = 0.8), width = 0.72,
             colour = "grey25", linewidth = 0.2) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1, labeller = label_parsed) +
    scale_fill_manual(values = PAL_CLUSTER, name = "Cluster") +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(out, p, width = w, height = 5.2, dpi = 200, bg = "white")
  cat("DONE ->", out, "\n")
}

mk_fig("local",
       "Local sensitivity (metric change per full range, at operating point)",
       "out_files/Chapter1/figures/FigSh_sensitivity_local.png")
mk_fig("global",
       "Intrinsic potency (metric change across full range)",
       "out_files/Chapter1/figures/FigSh_sensitivity_potency.png")
