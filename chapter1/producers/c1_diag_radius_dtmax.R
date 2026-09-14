# Diagnostic (supervisor request #5/#6, 2026-09-10): footprint-radius sensitivity
# of ΔTmax on the 400 cLHS plots, as the paired difference against the 20 m
# reference (NOT a correlation), per archetype. Pure re-plot of an existing
# table produced under CHS41-Rmerge / no wind correction.
# In : out_files/Chapter1/tables/radius_clhs_by_archetype.csv
# Out: out_files/Chapter1/figures/diag_2026-09-10/diag_radius_dtmax_clhs.png
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
R <- fread("out_files/Chapter1/tables/radius_clhs_by_archetype.csv")
R[, group := factor(group, levels = c("P1","P2","P3","P4","ALL"))]
PAL <- c(P1="#D7191C", P2="#FDAE61", P3="#A6D96A", P4="#1A9850", ALL="grey30")
# 20 m is the reference (d_dtmax = 0 by construction); add it for the eye
ref <- R[, .(radius = 20, d_dtmax = 0, dt_lo = 0, dt_hi = 0), by = group]
Rp <- rbind(R[, .(group, radius, d_dtmax, dt_lo, dt_hi)], ref)[order(group, radius)]
g <- ggplot(Rp, aes(radius, d_dtmax, colour = group, fill = group)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = "grey60") +
  geom_vline(xintercept = 20, linetype = "dotted", colour = "grey60") +
  geom_ribbon(aes(ymin = dt_lo, ymax = dt_hi), alpha = .12, colour = NA) +
  geom_line(linewidth = .6) + geom_point(size = 1.6) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_fill_manual(values = PAL, name = NULL) +
  facet_wrap(~group, nrow = 1) +
  labs(x = "footprint radius (m)",
       y = expression(Delta*italic(T)[max]~"vs 20 m reference (°C)")) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
        strip.text = element_text(face = "bold"))
ggsave("out_files/Chapter1/figures/diag_2026-09-10/diag_radius_dtmax_clhs.png",
       g, width = 11, height = 3, dpi = 300, bg = "white")
# console summary: max |deviation| from the 20 m reference per archetype
s <- R[, .(max_abs_dev = max(abs(d_dtmax)),
           at_radius = radius[which.max(abs(d_dtmax))],
           dev_5m = d_dtmax[radius == 5][1]), by = group]
print(s)
cat("\nn plots per radius (ALL):\n"); print(R[group=="ALL", .(radius, n)])
