# Diagnostic (supervisor request #2, 2026-09-10; corrected 2026-09-11): normalise
# each native-step effect by the plot's own baseline ΔTmax magnitude (`base`), so
# the relative weight of each lever is comparable across archetypes.
# IMPORTANT sign structure: base ΔTmax is POSITIVE (amplifying) in the open/inter-
# mediate archetypes (P1 +0.185, P2 +0.194, P3 +0.144 median) and NEGATIVE
# (buffering) in the dense P4 (-0.295). So the normalised fraction changes meaning
# across the amplify->buffer transition; read within an archetype, not across it.
# Only the 29 plots within |base| <= 0.05 °C of the sign change are dropped
# (0 in P1, 2 in P2, 14 in P3, 13 in P4); P1 is kept (its base is well away from 0).
# In : out_files/Chapter1/tables/perturb_chs41_nowind.csv
# Out: out_files/Chapter1/figures/diag_2026-09-10/diag_levers_normalised_by_buffering.png
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
source("R/cluster_relabel.R")
T <- fread("out_files/Chapter1/tables/perturb_chs41_nowind.csv")
T[, P := factor(P, levels = paste0("P",1:4))]
# native-step effects (per-unit columns x native step), and the profile contrast
E <- T[, .(id_plot, P, base,
           `Leaf area` = LAI_up*0.5, Cover = fCov_up*0.10,
           Height = Hmax_up, Profile = dT_LAD)]
Lm <- melt(E, id.vars = c("id_plot","P","base"), variable.name = "lever", value.name = "eff")
Lm[, frac := eff / abs(base)]                 # effect as a fraction of the plot's |ΔTmax|
Lm[, base_ok := abs(base) > 0.05]             # drop only plots straddling the amplify/buffer sign change
cat("=== median native-step effect as a fraction of the plot's own |ΔTmax| ===\n")
print(Lm[base_ok == TRUE, .(median_frac = round(median(frac, na.rm = TRUE), 3),
                            n = .N), by = .(P, lever)][order(lever, P)])
cat("\nplots dropped (|base| <= 0.05 °C, near the sign change), per archetype:\n")
print(Lm[lever=="Leaf area" & base_ok==FALSE, .N, by = P])
D <- Lm[base_ok == TRUE]
D[, lever := factor(lever, levels = c("Leaf area","Cover","Height","Profile"))]
g <- ggplot(D, aes(P, frac, fill = P)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = .3) +
  geom_violin(width = .9, linewidth = .2, colour = "grey40", draw_quantiles = c(.25,.5,.75)) +
  scale_fill_manual(values = PAL_CLUSTER, guide = "none") +
  facet_wrap(~lever, nrow = 1) +
  labs(x = NULL, y = "native-step effect / plot's own |ΔT_max|",
       subtitle = "Each lever / the plot's baseline ΔTmax. Base amplifies (+) in P1-P3, buffers (-) in P4; read within an archetype.") +
  coord_cartesian(ylim = c(-0.9, 0.9)) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                   strip.text = element_text(face = "bold"))
ggsave("out_files/Chapter1/figures/diag_2026-09-10/diag_levers_normalised_by_buffering.png",
       g, width = 11, height = 3.4, dpi = 300, bg = "white")
cat("DONE -> diag_levers_normalised_by_buffering.png\n")
