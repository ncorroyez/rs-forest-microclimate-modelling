# ==============================================================================
# Fig 5 (compact) — forward-inclusion as a CUMULATIVE-r CURVE (replaces the
# scatter grid FigSh_forward_percluster). One line per archetype: validation r
# (vs the 53 HOBO) as LiDAR traits are added in that archetype's own importance
# order. Reads the existing forward table (NO recompute / NO sims).
#   Rscript c1_forward_curve.R
# Out: out_files/Chapter1/figures/FigSh_forward_curve.{png,pdf}
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("scripts/_article_style.R") })
a <- fread("outputs/figures_pipeline_z05/tables/tab_forward_percluster.csv")
a <- a[step >= 1 & is.finite(r)]                                   # drop baseline (no r)
a[, trait := sub("^\\+", "", sub(" \\(REF\\)", "", add))]          # trait added at this step
a[, rowlab := factor(rowlab, levels = c("P1","P2","P3","P4","All"))]
PAL <- c(PAL_CLUSTER, All = "grey35")

p <- ggplot(a, aes(step, r, colour = rowlab, group = rowlab)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  geom_text(aes(label = trait), size = 2.5, vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_x_continuous(breaks = 1:4, labels = c("+1st","+2nd","+3rd","+4th"),
                     expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(x = "Traits added (each archetype's own importance order)",
       y = expression("Validation "*italic(r)*" vs 53 loggers"),
       subtitle = "Leaf quantity (LAI) added first lifts the fit most everywhere; the 2nd trait differs by density — fractional cover in the open (P1, P2), the vertical profile (LAD) in the dense archetypes (P3, P4).") +
  theme_article(11) + theme(legend.position = "bottom",
       plot.subtitle = element_text(size = 8, colour = "grey35"))
ggsave_article("out_files/Chapter1/figures/FigSh_forward_curve", p, 8.5, 4.8)
cat("DONE -> FigSh_forward_curve.{png,pdf}\n")
