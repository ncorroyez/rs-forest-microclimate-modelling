# ==============================================================================
# The figure JO asked for on 2026-08-21: log(slope) against Hmax, VCI and PAI,
# with one model version at 1 m WITH the ABL iteration and one at 1 m WITHOUT,
# everything else identical. Loggers kept as the reference cloud.
# Reads : out_files/Chapter1/tables/tab_noiter_compare.csv (scripts/c1_noiter_lineage.R)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv
# Writes: out_files/Chapter1/figures/Fig_logslope_iter_noiter.{png,pdf}
#   Rscript scripts/c1_logslope_iter_noiter.R
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(grid)})
source("scripts/_article_style.R")

D <- fread("out_files/Chapter1/tables/tab_noiter_compare.csv")
D <- merge(D, fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, LAI, VCI, Hmax)],
           by = "id_plot", suffixes = c("", ".c"))
D[, `:=`(obs = log(obs_sl), it = log(it_sl1m), no = log(no_sl1m))]

MET <- c(Hmax = "Maximum height (m)", VCI = "Vertical complexity index",
         LAI = "Leaf area index (one-sided)")
SRC <- c(obs = "Field sensors (HOBO, 1 m)",
         it  = "MuSICA at 1 m, with ABL iteration",
         no  = "MuSICA at 1 m, without ABL iteration")
L <- rbindlist(lapply(names(MET), function(m) rbindlist(lapply(names(SRC), function(k)
  data.table(metric = MET[[m]], x = D[[m]], y = D[[k]], src = SRC[[k]])))))
L[, `:=`(metric = factor(metric, levels = unname(MET)), src = factor(src, levels = unname(SRC)))]

PAL <- setNames(c("#128d84", "#E69F00", "#0072B2"), SRC)
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
ggsave_article("out_files/Chapter1/figures/Fig_logslope_iter_noiter", p, 10.5, 4.2)
cat("DONE\n")
