# ==============================================================================
# PIPELINE STAGE 06b — Metric (trait) RANKING deduced from the Shapley boxplots,
# for the HOBO validation set. Rank within each (metric × period × cluster) by
# |median φ| (1 = largest contribution). Robustness from the bootstrap CI (05b).
# Outputs: tables/tab_shapley_ranking.csv + fig_shapley_ranking.{png,pdf}
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
cli_h1("STAGE 06b — trait ranking (from Shapley boxplots)")

ci <- fread(file.path(PIPE$OUT_TAB, "tab_shapley_ci.csv"))   # from 05b
ci[, abs_med := abs(median)]
ci[, rank := frank(-abs_med, ties.method = "min"), by = .(metric, period, Cluster)]
setorder(ci, metric, period, Cluster, rank)
fwrite(ci[, .(metric, period, Cluster, rank, trait, median, ci_lo, ci_hi,
              robust = excludes_zero, direction)],
       file.path(PIPE$OUT_TAB, "tab_shapley_ranking.csv"))

cli_h2("Ranking — All (pooled HOBO set)")
print(ci[Cluster=="All", .(metric, period, rank, trait, median,
        CI=sprintf("[%.3f,%.3f]",ci_lo,ci_hi), direction)][order(metric,period,rank)])

# ---- article ranking figure: ordered Cleveland dot plot + bootstrap CI ------
# Headline = pooled HOBO set (All). Traits ordered by |median φ| within each
# metric×period panel (manual reorder_within; tidytext not available).
mlabs <- c(Tmax="Delta*T[max]~(degree*C)", VPD="Delta*VPD[max]~(kPa)")
tl <- c(LAI="LAI", Hmax="italic(H)[max]", fCover="fCover", LAD="LAD")
pl_lab <- c("All period"="All summer", "10% hottest"="10% hottest days")

pan_lvls <- c(paste0(mlabs["Tmax"],"*' — ", pl_lab["All period"], "'"),
              paste0(mlabs["Tmax"],"*' — ", pl_lab["10% hottest"], "'"),
              paste0(mlabs["VPD"], "*' — ", pl_lab["All period"], "'"),
              paste0(mlabs["VPD"], "*' — ", pl_lab["10% hottest"], "'"))
d <- ci[Cluster == "All"]
d[, panel := factor(paste0(mlabs[metric], "*' — ", pl_lab[period], "'"), levels = pan_lvls)]
d[, key := paste(trait, metric, period, sep = "___")]
d <- d[order(panel, abs_med)]                              # ascending → largest on top
d[, key := factor(key, levels = key)]
d[, dir := factor(direction, levels = c("buffers","amplifies","ns"))]
dir_col <- c(buffers="#2166AC", amplifies="#B2182B", ns="grey60")

p <- ggplot(d, aes(median, key, colour = dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22, linewidth = 0.8) +
  geom_point(size = 3.4) +
  geom_text(aes(label = paste0("#", rank)), colour = "grey30", size = 3.6,
            hjust = -0.0, vjust = -1.0) +
  facet_wrap(~ panel, scales = "free", ncol = 2, labeller = label_parsed) +
  scale_y_discrete(labels = function(k) parse(text = tl[sub("___.*", "", k)])) +
  scale_colour_manual(values = dir_col, name = NULL,
                      labels = c(buffers="buffers (φ<0)", amplifies="amplifies (φ>0)", ns="CI spans 0")) +
  labs(x = expression("Shapley "*varphi[v]*"  (median ± bootstrap 95% CI)"), y = NULL,
       title = "Trait ranking — pooled HOBO set",
       caption = "Ordered by |median φ| within each panel. #1 = largest contribution. Per-cluster ranking: see fig_shapley_by_cluster.") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 14), panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.caption = element_text(size = 10.5, colour = "grey40"))
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_ranking.png"), p, width = 11, height = 7.5, dpi = 300, bg = "white")
ggsave(file.path(PIPE$OUT_FIG,"fig_shapley_ranking.pdf"), p, width = 11, height = 7.5, device = cairo_pdf)
cli_alert_success("Saved fig_shapley_ranking (Cleveland dot + CI) + tab_shapley_ranking.csv")
