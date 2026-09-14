# ==============================================================================
# Figure 7: does the ALS footprint radius change the validation? Pooled and per archetype,
# on the aligned clock. The per-group n in the legend (53/8/12/13/20) is the frozen split.
#
# Reads : out_files/Chapter1/tables/tab_radius_sweep_aligned.csv   (stage A5)
# Writes: out_files/Chapter1/figures/Fig7_footprint_radius.{png,pdf}
#   Rscript scripts/c1_fig7_footprint_radius.R
# ==============================================================================
# Figure 7 — footprint-radius sensitivity of the between-plot DeltaTmax fit to the
# 53 loggers, pooled and within each archetype, on the corrected baseline.
#
# WRITTEN 2026-07-29. Figure 7 previously had NO generating script in the repo: it
# was produced ad hoc and its underlying table (tab_radius_by_cluster_convB.csv)
# predated the clock alignment, so its values ran about 0.03 high. This script
# rebuilds it from tab_radius_sweep_aligned.csv (c1_radius_rescore.R), which scores
# the `radius_test_corr` lineage (canopy wind corrected, u = 1.233, matching the
# reported baseline) on the aligned clock.
suppressPackageStartupMessages({library(data.table);library(ggplot2);source("scripts/_article_style.R")})
R <- fread("out_files/Chapter1/tables/tab_radius_sweep_aligned.csv")
L <- melt(R, id.vars = "radius", measure.vars = c("r","P1","P2","P3","P4"),
          variable.name = "grp", value.name = "corr")
L[, grp := factor(ifelse(grp == "r", "All", as.character(grp)),
                  levels = c("All","P1","P2","P3","P4"))]
NLAB <- c(All = 53, P1 = 8, P2 = 12, P3 = 13, P4 = 20)
L[, lab := paste0(grp, " (n = ", NLAB[as.character(grp)], ")")]
L[, lab := factor(lab, levels = unique(lab[order(grp)]))]
PAL <- c(setNames("grey25", levels(L$lab)[1]),
         setNames(unname(PAL_CLUSTER[c("P1","P2","P3","P4")]), levels(L$lab)[2:5]))
p <- ggplot(L, aes(radius, corr, colour = lab, group = lab)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45") +
  geom_line(aes(linewidth = grp == "All")) +
  geom_point(aes(size = grp == "All")) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_linewidth_manual(values = c(0.7, 1.6), guide = "none") +
  scale_size_manual(values = c(2, 3.2), guide = "none") +
  scale_x_continuous(breaks = R$radius) +
  labs(x = "footprint radius (m)",
       y = expression("between-plot "*italic(r)*" against the sub-canopy measurements ("*Delta*T[max]*")")) +
  theme_article(12)
ggsave_article("out_files/Chapter1/figures/Fig7_footprint_radius", p, 8, 5)
print(R[, .(radius, pooled = round(r,3), P1 = round(P1,2), P2 = round(P2,2),
            P3 = round(P3,2), P4 = round(P4,2))])
cat("DONE\n")
