# ==============================================================================
# Sensitivity-based trait IMPORTANCE per archetype — the NON-Shapley replacement
# for the old per-cluster mean|phi| figure. No coalitions, no global-mean baseline,
# no off-manifold canopies: every perturbation is a small one-at-a-time move at
# each plot's OWN real operating point (on the data manifold), the readout asked
# for by the supervisors ("sensibilite au sein de chaque cluster").
#
# Importance of a trait within an archetype = how much ΔTmax that trait moves over
# its naturally realized variation in that archetype:
#   - continuous traits (LAI, fCover, Hmax): |median local per-unit slope| ×
#       within-archetype dispersion of the trait (SD, native units; LAI one-sided)
#   - vertical profile (LAD), categorical real-vs-uniform: |median(real − uniform)|
#       contrast at each plot's own real LAI/Hmax/fCover
# These are two genuinely distinct quantities (a slope×spread vs a shape contrast),
# reported side by side exactly as the analysis frames them.
#
# Reads the per-plot sensitivity (c1_sensitivity_perplot_chunk.R) and the per-plot
# LAD contrast (c1_sensitivity_perplot_lad.R) — NO sims.
#   Rscript c1_importance_sensitivity.R
# Out: tab_importance_sensitivity.csv + FigSh_importance_sensitivity.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R")
})

S <- rbindlist(lapply(list.files("out_files/Chapter1/tables/sensitivity_perplot",
                                 "part_.*csv$", full.names = TRUE), fread), fill = TRUE)
stopifnot(nrow(S) > 0)
S[, P := relabel_cluster(Cluster)]
LAD <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")
LAD[, P := relabel_cluster(Cluster)]

# per-plot local per-unit effects are already ΔTmax per native unit:
#   LAI_*  = per 1 ONE-SIDED LAI unit   | the cLHS LAI column is two-sided -> SD/2
#   fCov_* = per 0.1 fCover             | Hmax_* = per 5 m
# compute importance grouped (per archetype) AND pooled ("All") in one pass
imp_block <- function(Ssub, LADsub) {
  a <- Ssub[, .(
    s_LAI = median(abs((LAI_add  - LAI_rem )/2), na.rm = TRUE),       # |per 1 one-sided LAI unit|
    s_fC  = median(abs((fCov_add - fCov_rem)/2), na.rm = TRUE) / 0.1, # per fCover unit
    s_Hm  = median(abs((Hmax_add - Hmax_rem)/2), na.rm = TRUE) / 5,   # per m
    sd_LAI = sd(LAI, na.rm = TRUE) / 2,                               # one-sided dispersion
    sd_fC  = sd(fCover, na.rm = TRUE),
    sd_Hm  = sd(Hmax, na.rm = TRUE))]
  data.table(LAI = a$s_LAI * a$sd_LAI, fCover = a$s_fC * a$sd_fC,
             Hmax = a$s_Hm * a$sd_Hm, LAD = abs(median(LADsub$dT_LAD, na.rm = TRUE)))
}
IMP <- rbind(
  S[, imp_block(.SD, LAD[P == .BY$P]), by = P],
  cbind(P = "All", imp_block(S, LAD))
)
IMP[, c("LAI","fCover","Hmax","LAD") := lapply(.SD, round, 3), .SDcols = c("LAI","fCover","Hmax","LAD")]
IMP <- IMP[order(P)]
fwrite(IMP, "out_files/Chapter1/tables/tab_importance_sensitivity.csv")
cat("=== sensitivity-based importance (°C of ΔTmax generated within archetype) ===\n")
print(IMP)
cat("\nP4 co-lead check: LAI", IMP[P=="P4", LAI], "vs LAD", IMP[P=="P4", LAD], "\n")

# ---- bootstrap 95% CIs per trait per archetype (point estimates above unchanged) ----
# resample plots within each archetype, recompute the SAME importance formulas; the
# overlap of the LAI and LAD intervals in P4 is the visual proof of the near-tie co-lead.
set.seed(1); NB <- 2000
Sb <- merge(S, LAD[, .(pid, dT_LAD)], by = "pid", all.x = TRUE)
imp_boot <- function(d) {
  s_LAI <- median(abs((d$LAI_add - d$LAI_rem)/2), na.rm = TRUE)
  s_fC  <- median(abs((d$fCov_add - d$fCov_rem)/2), na.rm = TRUE) / 0.1
  s_Hm  <- median(abs((d$Hmax_add - d$Hmax_rem)/2), na.rm = TRUE) / 5
  c(LAI = s_LAI*sd(d$LAI,na.rm=TRUE)/2, fCover = s_fC*sd(d$fCover,na.rm=TRUE),
    Hmax = s_Hm*sd(d$Hmax,na.rm=TRUE), LAD = abs(median(d$dT_LAD,na.rm=TRUE)))
}
CI <- rbindlist(lapply(c("P1","P2","P3","P4"), function(p){
  d <- Sb[P == p]; n <- nrow(d)
  B <- t(replicate(NB, imp_boot(d[sample(n, n, TRUE)])))
  data.table(P = p, trait = factor(c("LAI","fCover","Hmax","LAD")),
             lo = apply(B, 2, quantile, .025, names = FALSE),
             hi = apply(B, 2, quantile, .975, names = FALSE))
}))

# ---- figure: grouped bars, archetype × trait (P1–P4 only; "All" is for ordering) ----
L <- melt(IMP[P != "All"], id.vars = "P", variable.name = "trait", value.name = "imp")
L <- merge(L, CI, by = c("P", "trait"))
TLAB <- c(LAI = "LAI (leaf quantity)", fCover = "fCover", Hmax = "H[max]",
          LAD = "LAD (vertical profile)")
L[, trait := factor(trait, levels = names(TLAB))]
PAL <- c(LAI = "#1B7837", fCover = "#7FBC41", Hmax = "#BDBDBD", LAD = "#762A83")

p <- ggplot(L, aes(P, imp, fill = trait)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8),
                width = 0.2, colour = "grey30", linewidth = 0.35) +
  geom_text(aes(y = hi, label = sprintf("%.2f", imp)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2.3, colour = "grey25") +
  scale_fill_manual(values = PAL, labels = TLAB[names(PAL)], name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = expression(Delta*T[max]~importance~(degree*C)),
       title = "Sensitivity-based trait importance per archetype (non-Shapley)",
       subtitle = paste0("Leaf quantity dominates throughout; its lead erodes as the canopy saturates ",
                         "(P1→P4) until, in dense P4, the\nvertical profile becomes comparable to (co-leads with) leaf quantity ",
                         "(overlapping 95% CIs = a near-tie),\nhaving already overtaken cover and height. Error bars: 2000-resample bootstrap within archetype.")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"))
ggsave("out_files/Chapter1/figures/FigSh_importance_sensitivity.png", p,
       width = 8.5, height = 5.2, dpi = 200, bg = "white")
cat("DONE -> FigSh_importance_sensitivity.png + tab_importance_sensitivity.csv\n")
