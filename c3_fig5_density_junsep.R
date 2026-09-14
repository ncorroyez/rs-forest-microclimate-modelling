# ==============================================================================
# c3_fig5_density_junsep.R — Chapter 3, Figure 5: LiDAR leaf area against point
# density. Written 2026-09-11 because the figure in the manuscript
# (Fig5_density.png) had no producer script in either repository, so it could
# not be regenerated or checked.
#
# Source table: chapter3_S2_LAI/tables/reframe_density_pulse.csv
#   53 plots x 5 pulse-thinning fractions (100, 75, 50, 25, 10 % of native).
#   The _10m variant is the one the chapter reports: it is computed on the 10 m
#   grid the whole chapter works on. The other file is a finer-grid variant and
#   gives slightly higher rank correlations.
# Out: NC_Full/manuscripts/ch3/figures/Fig5_reframe_junsep.png
# ==============================================================================
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source("R/cluster_relabel.R")
TAB <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
FIG <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures"
th  <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())

D <- fread(file.path(TAB, "reframe_density_pulse_10m.csv"))
D[, P := factor(P, levels = paste0("P", 1:4))]
D[, pct := 100 * frac]

# (b) between-plot rank correlation against the full-density leaf area
full <- D[frac == 1, .(id_plot, LAI_full = LAI)]
R <- merge(D, full, by = "id_plot")[, .(rho = cor(LAI, LAI_full, method = "spearman"),
                                        dmed = median(abs(LAI - LAI_full))), by = pct][order(pct)]
print(R)
stopifnot(nrow(R) == 5, all(R$rho >= 0.99))          # the claim the figure makes
stopifnot(abs(R[pct == 10, dmed] - 0.071) < 0.002)   # the value the chapter reports

pa <- ggplot(D, aes(pct, LAI, colour = P, fill = P)) +
  stat_summary(fun = median, geom = "line", linewidth = .7) +
  stat_summary(fun = median, geom = "point", size = 2) +
  stat_summary(fun.min = function(x) quantile(x, .25),
               fun.max = function(x) quantile(x, .75),
               geom = "ribbon", alpha = .18, colour = NA) +
  scale_colour_manual(values = PAL_CLUSTER, name = "Archetype") +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  labs(x = "Point density retained (% of native)",
       y = expression("LiDAR LAI (m"^2*" m"^-2*")"), tag = "(a)") + th

pb <- ggplot(R, aes(pct, rho)) +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey55") +
  geom_line(linewidth = .7) + geom_point(size = 2) +
  labs(x = "Point density retained (% of native)",
       y = "Spearman ρ (vs 100 %)", tag = "(b)") + th

ggsave(file.path(FIG, "Fig5_reframe_junsep.png"), pa + pb, width = 9, height = 4.1,
       dpi = 300, bg = "white")
cat("DONE: Fig5_reframe_junsep.png\n")
