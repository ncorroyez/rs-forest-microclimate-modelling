# ==============================================================================
# The log(slope) figure read at the model's LOWEST air node instead of 1 m.
# Companion to scripts/c1_sampling_height_test.R; see that script's header for
# why a fixed 0.1 m is not available (no air node below 0.0285 x Hmax, and
# micro_hourly_at clamps rather than extrapolates) and for the Hmax confound.
# Reads : out_files/Chapter1/tables/tab_sampling_height_test.csv
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv
# Writes: out_files/Chapter1/figures/Fig_logslope_sampling_height.{png,pdf}
#   Rscript scripts/c1_logslope_sampling_height.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(grid)})
source("scripts/_article_style.R")

D <- fread("out_files/Chapter1/tables/tab_sampling_height_test.csv")
C <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, LAI, VCI)]
D <- merge(D, C, by = "id_plot")
D[, `:=`(obs = log(obs_sl), s1 = log(sl_1m), s0 = log(sl_low))]

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c(obs = "Field sensors (HOBO, 1 m)",
         s1  = "MuSICA read at 1 m",
         s0  = "MuSICA read at its lowest node (0.09 to 1.08 m)")
L <- rbindlist(lapply(names(MET), function(m) rbindlist(lapply(names(SRC), function(k)
  data.table(metric = MET[[m]], x = D[[m]], y = D[[k]], src = SRC[[k]])))))
L[, `:=`(metric = factor(metric, levels = unname(MET)), src = factor(src, levels = unname(SRC)))]

PAL <- setNames(c("#128d84", "#E69F00", "#8B2500"), SRC)
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
ggsave_article("out_files/Chapter1/figures/Fig_logslope_sampling_height", p, 10.5, 4.2)
cat("DONE\n")
