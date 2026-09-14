# ==============================================================================
# Figure F2: what each correction level does to the validation, per logger. The scan-angle
# and the k = 0.65 corrections are DIFFERENT corrections; merging them is how this got
# confused originally, so they stay separate columns throughout.
#
# Reads : out_files/Chapter1/tables/tab_correction_levels{,_perplot}.csv   (stage A6)
# Writes: out_files/Chapter1/figures/F2_lai_correction_compare.{png,pdf}
#   Rscript scripts/fig_lai_correction_compare.R
# ==============================================================================
# Fig. F2 — effect of the leaf-area-retrieval corrections on the microclimate
# model, scored against the 53 HOBO loggers.
#
# REBUILT 2026-07-29. The previous version was wrong twice over: it read the
# `frblo` lineage (not native20) and it HARD-CODED the middle level at
# (r = 0.47, bias = 0.72, amp = 11%), which are the k = 0.65 endpoint's own
# values, so panel (b) drew a three-point line in which two points were the same
# run. Both panels now come from tab_correction_levels*.csv (native20, aligned
# clock, time-matched convention), produced by c1_correction_levels_rescore.R.
#
# The panel shows ONE axis, leaf-area retrieval, all three levels sharing the
# canopy wind correction:
#   wind only (k = 0.5)  ->  + scan angle (k = 0.5, the reported baseline)
#                        ->  + optimized k = 0.65 (on top of the scan angle).
# The wind axis (uncorrected wind) is deliberately NOT placed on this line; it
# is a different correction and merging the two is how the original got confused.
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("scripts/_article_style.R")
TAB <- "out_files/Chapter1/tables"
R  <- fread(file.path(TAB, "tab_correction_levels.csv"))
PP <- fread(file.path(TAB, "tab_correction_levels_perplot.csv"))
LV <- c("wind only\n(k = 0.5)", "+ scan angle\n(k = 0.5)", "+ optimized k\n(0.65)")
key <- c(windonly = LV[1], full = LV[2], k065 = LV[3])

# --- panel a: sim vs obs, reported baseline vs optimized k -------------------
sc <- PP[variant %in% c("full", "k065") & is.finite(obs_dt) & is.finite(sim_dt)]
sc[, set := ifelse(variant == "full", "baseline (k = 0.5 + scan angle)",
                                      "optimized k = 0.65")]
lim <- range(c(sc$obs_dt, sc$sim_dt), na.rm = TRUE)
pa <- ggplot(sc, aes(obs_dt, sim_dt, colour = set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = c("baseline (k = 0.5 + scan angle)" = "grey30",
                                 "optimized k = 0.65" = "#D7191C"), name = NULL) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = expression("Observed  "*Delta*T[max]*"  (°C)"),
       y = expression("Simulated  "*Delta*T[max]*"  (°C)"),
       subtitle = "The optimized k collapses the simulated range toward zero") +
  theme_article(12) + theme(legend.position = "bottom")

# --- panel b: metrics across the three leaf-area correction levels -----------
M <- R[variant %in% names(key)][, .(level = factor(key[variant], levels = LV),
                                    r = r, bias = bias, amp = amp)]
stopifnot(nrow(M) == 3L, uniqueN(round(M$r, 3)) == 3L)   # guard: no duplicated level
Ml <- melt(M, id.vars = "level")
Ml[, variable := factor(variable, levels = c("r", "bias", "amp"),
     labels = c("correlation r", "warm bias (°C)", "amplitude captured (%)"))]
pb <- ggplot(Ml, aes(level, value, group = 1)) +
  geom_line(colour = "#D7191C") + geom_point(colour = "#D7191C", size = 2.4) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(x = NULL, y = NULL) + theme_article(11)

p <- pa + pb + plot_layout(widths = c(1.5, 1)) + plot_annotation(tag_levels = "a")
ggsave_article("out_files/Chapter1/figures/F2_lai_correction_compare", p, 10, 6)
print(M)
cat("Saved F2_lai_correction_compare\n")
