# ==============================================================================
# log(slope) vs structure on the doubled-grid lineage, read near the ground.
# Companion to scripts/c1_airres1_lineage.R; see that header for the two-step
# reading (grid effect first, sampling height second) and the caveats.
# Three series: the loggers, the chapter's baseline, and the 20-layer grid read
# at its lowest node (0.02 to 0.26 m).
# Reads : out_files/Chapter1/tables/tab_airres1_compare.csv
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv
# Writes: out_files/Chapter1/figures/Fig_logslope_airres1.{png,pdf}
#   Rscript scripts/c1_logslope_airres1.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(grid)})
source("scripts/_article_style.R")

D <- fread("out_files/Chapter1/tables/tab_airres1_compare.csv")
D <- merge(D, fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, LAI, VCI)], by = "id_plot")
D[, `:=`(obs = log(obs_sl), base = log(g10_sl1m), low = log(g20_sllow))]

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c(obs  = "Field sensors (HOBO, 1 m)",
         base = "MuSICA, 10 layers, read at 1 m",
         low  = "MuSICA, 20 layers, read at 0.02 to 0.26 m")
L <- rbindlist(lapply(names(MET), function(m) rbindlist(lapply(names(SRC), function(k)
  data.table(metric = MET[[m]], x = D[[m]], y = D[[k]], src = SRC[[k]])))))
L[, `:=`(metric = factor(metric, levels = unname(MET)), src = factor(src, levels = unname(SRC)))]

PAL <- setNames(c("#128d84", "#E69F00", "#5E3C99"), SRC)
p <- ggplot(L, aes(x, y, colour = src)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.9) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +
  scale_colour_manual(values = PAL, name = NULL) +
  facet_wrap(~ metric, nrow = 1, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = "log(slope)") + theme_article(12) +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.background = element_blank(),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(4, "pt")) +
  guides(colour = guide_legend(nrow = 1))
ggsave_article("out_files/Chapter1/figures/Fig_logslope_airres1", p, 10.5, 4.2)
cat(sprintf("log(slope) span: obs %.3f | baseline %.3f | 20 layers low %.3f\n",
            diff(range(D$obs)), diff(range(D$base)), diff(range(D$low))))
cat(sprintf("plots with slope > 1: obs %d | baseline %d | 20 layers low %d\n",
            sum(D$obs_sl > 1), sum(D$g10_sl1m > 1), sum(D$g20_sllow > 1)))
cat("DONE\n")
