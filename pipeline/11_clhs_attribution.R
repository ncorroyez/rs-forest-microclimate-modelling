# ==============================================================================
# STAGE 11 — cLHS Shapley ATTRIBUTION figures (article Fig 2 & Fig 3).
# Canonical, unambiguous, LOCAL. Reads the per-point Shapley parts computed on the
# ARTICLE-canonical sample clhs_sample_floor05_v2.rds (NOT floor05):
#   per-cluster baseline -> out_files/Chapter1/tables/shapley_parts/        (Fig 2)
#   global  baseline     -> out_files/Chapter1/tables/shapley_parts_global/ (Fig 3)
# (Compute side = c3_shapley_chunk.R / c1_shapley_global_chunk.R over the cached
#  6400 nc in out_files/Chapter1/nc_shapley2x[_global]; no MuSICA re-run here.)
# Writes the EXACT names the manuscript cites, to the LOCAL z_Example figures dir.
#   Rscript pipeline/11_clhs_attribution.R
# Out: out_files/Chapter1/figures/FigSh_shapley_clhs_cluster.png  (Fig 2)
#      out_files/Chapter1/figures/FigSh_shapley_clhs_global.png   (Fig 3)
#      + tables tab_shapley_clhs_{cluster_importance,global_signed}.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); source("R/cluster_relabel.R") })
FIG <- "out_files/Chapter1/figures"; TAB <- "out_files/Chapter1/tables"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
Fv  <- c("LAI","fCover","LAD","Hmax")          # display order: quantity, cover, profile, height
read_parts <- function(sub) {
  f <- list.files(file.path(TAB, sub), pattern = "part_.*\\.csv$", full.names = TRUE)
  stopifnot(length(f) > 0)
  S <- rbindlist(lapply(f, fread), fill = TRUE)
  S[is.finite(dTmax_full) & is.finite(dTmax_base)]
}

# ---- Fig 2 : per-cluster baseline, within-archetype importance = mean |phi| ----
SC <- read_parts("shapley_parts")
SC[, P := factor(relabel_cluster(Cluster), levels = c("P1","P2","P3","P4"))]
impC <- melt(SC, id.vars = c("pid","P"), measure.vars = Fv, variable.name = "trait", value.name = "phi")[
  , .(imp = mean(abs(phi), na.rm = TRUE)), by = .(P, trait)]
impC[, trait := factor(trait, levels = Fv)]
fwrite(dcast(impC, P ~ trait, value.var = "imp"), file.path(TAB, "tab_shapley_clhs_cluster_importance.csv"))
cat("=== within-archetype importance mean|phi| (ΔTmax, °C) ===\n"); print(dcast(impC, P ~ trait, value.var = "imp"))

g2 <- ggplot(impC, aes(trait, imp, fill = trait)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", imp)), vjust = -0.3, size = 3.0) +
  facet_wrap(~ P, nrow = 1) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression("mean |"*phi*"|  (°C on "*Delta*T[max]*")")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 0))
ggsave(file.path(FIG, "FigSh_shapley_clhs_cluster.png"), g2, width = 9.2, height = 3.4, dpi = 200, bg = "white")
cat("DONE -> FigSh_shapley_clhs_cluster.png\n")

# ---- Fig 3 : common (global) baseline, landscape importance = mean |phi| -------
# (Importance, like Fig 2. The SIGNED phi straddles zero here because a single
#  global baseline imposes off-manifold canopies on atypical stands — see
#  manuscript caveat — so the directional readout is not robust; the importance
#  ranking is. We report both: bars = mean|phi|, table also carries signed stats.)
SG <- read_parts("shapley_parts_global")
mlG <- melt(SG, id.vars = "pid", measure.vars = Fv, variable.name = "trait", value.name = "phi")
mlG[, trait := factor(trait, levels = Fv)]
sumG <- mlG[, .(imp = mean(abs(phi), na.rm = TRUE), median_signed = median(phi, na.rm = TRUE),
                mean_signed = mean(phi, na.rm = TRUE), frac_neg = mean(phi < 0, na.rm = TRUE)), by = trait]
fwrite(sumG, file.path(TAB, "tab_shapley_clhs_global.csv"))
cat("\n=== global-baseline landscape importance + signed stats ===\n"); print(sumG)

g3 <- ggplot(sumG, aes(trait, imp, fill = trait)) +
  geom_col(width = 0.66, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", imp)), vjust = -0.3, size = 3.2) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = expression("mean |"*phi*"|  (°C on "*Delta*T[max]*")")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(FIG, "FigSh_shapley_clhs_global.png"), g3, width = 6.6, height = 3.6, dpi = 200, bg = "white")
cat("DONE -> FigSh_shapley_clhs_global.png\n")
