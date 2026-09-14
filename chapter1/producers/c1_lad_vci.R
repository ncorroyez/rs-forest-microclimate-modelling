# ==============================================================================
# Scalar interpretation of the vertical-profile (LAD) effect via the VCI.
# The real-vs-uniform contrast dT_LAD IS conceptually "VCI_real -> VCI=1 (uniform,
# maximal vertical evenness)". Here we show the contrast scales with the profile's
# departure from uniform, (1 - VCI), WITHIN each archetype (i.e. density-controlled),
# which the bulk centre of mass failed to capture. VCI is a DESCRIPTOR/normalizer of
# the contrast, NOT a MuSICA input (MuSICA ingests the profile, not the scalar) — and
# remains necessary-but-not-sufficient (collinear with density, one-to-many with shape).
#   Rscript c1_lad_vci.R
# Out: FigSh_lad_vci.png + tab_lad_vci.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); source("R/cluster_relabel.R")
})
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N, 100))], by = Cluster]
sub[, pid := sprintf("S%04d", .I)]; sub[, P := relabel_cluster(Cluster)]
LAD <- fread("out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")
D <- merge(sub[, .(pid, P, VCI, LAI, dep = 1 - VCI)], LAD[, .(pid, dT_LAD)], by = "pid")
D <- D[is.finite(dT_LAD) & is.finite(VCI)]

sp <- function(a, b) cor(a, b, method = "spearman", use = "complete.obs")
partial <- { r1 <- resid(lm(dT_LAD ~ LAI, D)); r2 <- resid(lm(VCI ~ LAI, D)); sp(r1, r2) }
tab <- D[, .(rho_dT_VCI = round(sp(dT_LAD, VCI), 2),
             rho_dT_dep = round(sp(dT_LAD, dep), 2),
             median_VCI = round(median(VCI), 2), n = .N), by = P][order(P)]
tab <- rbind(tab, data.table(P = "All", rho_dT_VCI = round(sp(D$dT_LAD, D$VCI), 2),
             rho_dT_dep = round(sp(D$dT_LAD, D$dep), 2),
             median_VCI = round(median(D$VCI), 2), n = nrow(D)))
fwrite(tab, "out_files/Chapter1/tables/tab_lad_vci.csv")
cat("=== dT_LAD vs VCI by archetype ===\n"); print(tab)
cat(sprintf("\nVCI~LAI collinearity rho = %+.2f ; partial rho(dT_LAD,VCI|LAI) = %+.2f\n",
            sp(D$VCI, D$LAI), partial))

PAL <- c(P1 = "#D7191C", P2 = "#FDAE61", P3 = "#A6D96A", P4 = "#1A9641")
lab <- D[, .(rho = sp(dT_LAD, dep)), by = P]
lab[, txt := sprintf("%s: rho=%+.2f", P, rho)]
p <- ggplot(D, aes(dep, dT_LAD, colour = P)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.45, size = 1.1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = PAL, name = NULL) +
  labs(x = "Departure from a uniform profile,  1 - VCI  (0 = uniform / VCI=1)",
       y = expression(dT[LAD]==Delta*T[max](real)-Delta*T[max](uniform)~~(degree*C)),
       title = "The vertical-profile effect scales with departure from uniform (VCI), within each archetype",
       subtitle = sprintf(paste0("Density-controlled: partial rho(dT_LAD, VCI | LAI) = %+.2f (more concentrated profile, lower VCI, ",
                  "larger real-vs-uniform contrast).\nThe cross-archetype growth of the contrast is the density-saturation of the ",
                  "other levers; VCI explains the within-archetype scatter (modestly), not the gradient."), partial)) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.3, colour = "grey35"))
ggsave("out_files/Chapter1/figures/FigSh_lad_vci.png", p, width = 8.5, height = 5.4, dpi = 200, bg = "white")
cat("DONE -> FigSh_lad_vci.png + tab_lad_vci.csv\n")
