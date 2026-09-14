# ==============================================================================
# Complementary analyses on the v3.2.0 vs v3.2.3 comparison (no new sims).
#   #1 version difference Δ(v323-v320) BY ARCHETYPE P1->P4 (is the yoyo effect
#      concentrated in the open canopy?)
#   #2 does each version preserve the OBSERVED P1->P4 buffering gradient?
#   #3 robustness: global ΔTmax skill WITH vs WITHOUT the 9 amplifying gaps
#   #4 Taylor diagram (correlation + normalized amplitude) per metric x version
#   #5 structural predictor of the version divergence |v323-v320|
#   Run from repo root:
#   Rscript Chapitre1/comparaison_versions/scripts/02_analyses_microclimat/cmp_v320_v323_complements.R
# Out: figures/doc_media/{image_cmp_byarchetype,image_gradient,image_taylor,image_divergence}.png
#      tables/tab_v320_v323_robustness.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(plotrix); library(ggrepel)
})
DIR <- "Chapitre1/comparaison_versions"; TBL <- file.path(DIR, "tables"); MED <- file.path(DIR, "figures/doc_media")
PAL  <- c(P1 = "#D9A441", P2 = "#7FBC41", P3 = "#3690C0", P4 = "#1B7837")
SRC  <- c(observé = "black", "v3.2.0" = "#0072B2", "v3.2.3" = "#D55E00")

cl  <- as.data.table(readRDS("outputs/figures_pipeline_z05/data/clusters.rds")); setnames(cl, "Cluster", "P")
rj  <- fread(file.path(TBL, "recap_JJAS_plotdata.csv"))
ex  <- fread(file.path(TBL, "compare_binaries_extra_axes.csv"))
ex[src == "v3.2.0", src := "v320"][src == "v3.2.3 iter", src := "v323"][src == "Observed", src := "obs"]
exw <- dcast(ex, id_plot ~ src, value.var = c("slope_night", "DTR", "dTmin"))

# unified per-plot table: obs / v320 / v323 for the 5 metrics
d <- merge(rj, exw, by = "id_plot"); d <- merge(d, cl, by = "id_plot")
MET <- data.table(
  id   = c("dTmaxS", "slope", "DTR", "dTmin", "slopeN"),
  lab  = c("ΔTmax été", "pente (tot)", "DTR", "ΔTmin", "pente (nuit)"),
  unit = c("°C", "", "°C", "°C", ""),
  obs  = c("dTmaxS_obs", "slope_obs", "DTR_obs", "dTmin_obs", "slope_night_obs"),
  v0   = c("dTmaxS_v320", "slope_v320", "DTR_v320", "dTmin_v320", "slope_night_v320"),
  v3   = c("dTmaxS_v323", "slope_v323", "DTR_v323", "dTmin_v323", "slope_night_v323"))

# ============================== #1 Δ by archetype =============================
diffL <- rbindlist(lapply(seq_len(nrow(MET)), function(i) {
  data.table(P = d$P, lab = MET$lab[i], unit = MET$unit[i],
             delta = d[[MET$v3[i]]] - d[[MET$v0[i]]])
}))
agg <- diffL[, .(m = mean(delta), se = sd(delta) / sqrt(.N)), by = .(lab, unit, P)]
agg[, facet := sprintf("%s (%s)", lab, ifelse(unit == "", "—", unit))]
f1 <- ggplot(agg, aes(P, m, fill = P)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_col(width = .7) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .25) +
  facet_wrap(~facet, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = PAL, guide = "none") +
  labs(x = NULL, y = "Δ (v3.2.3 − v3.2.0)",
       title = "1. Effet du couplage par archétype : Δ(v3.2.3 − v3.2.0), moyenne ± SE par type",
       subtitle = "Effet contrasté jour/nuit selon la densité : le refroidissement diurne (ΔTmax) est concentré dans l'ouvert (P1, ≈ −0,9 °C) ; le réchauffement nocturne (ΔTmin) et la baisse de la pente nocturne se concentrent au contraire dans le dense (P4) — signature de l'inversion nocturne captée par le couplage.") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
       plot.subtitle = element_text(size = 8.3, colour = "grey35"))
ggsave(file.path(MED, "image_cmp_byarchetype.png"), f1, width = 11, height = 3.6, dpi = 200, bg = "white")

# ============================== #2 gradient preservation ======================
grad <- rbind(
  d[, .(src = "observé", v = mean(dTmaxS_obs), se = sd(dTmaxS_obs)/sqrt(.N)), by = P],
  d[, .(src = "v3.2.0",  v = mean(dTmaxS_v320), se = sd(dTmaxS_v320)/sqrt(.N)), by = P],
  d[, .(src = "v3.2.3",  v = mean(dTmaxS_v323), se = sd(dTmaxS_v323)/sqrt(.N)), by = P])
f2 <- ggplot(grad, aes(P, v, colour = src, group = src)) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = 3) +
  geom_line(linewidth = .9) +
  geom_pointrange(aes(ymin = v - se, ymax = v + se), size = .5) +
  scale_colour_manual(values = SRC, name = NULL) +
  labs(x = "Archétype (ouvert → dense)", y = "ΔTmax été moyen (°C)",
       title = "2. Préservation du gradient observé de tampon P1 → P4",
       subtitle = "L'observation décroît de façon monotone (trouée chaude → sous-bois frais). Les deux versions retrouvent l'ordre ; v3.2.0 garde la magnitude du gradient, v3.2.3 le tasse vers 0.") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "right",
       plot.subtitle = element_text(size = 8.3, colour = "grey35"))
ggsave(file.path(MED, "image_gradient.png"), f2, width = 8.2, height = 4.6, dpi = 200, bg = "white")

# ============================== #3 robustness (no amplifiers) =================
d[, grp := ifelse(dTmaxS_obs > 0, "amplifier", "buffer")]
skill <- function(o, s) data.table(pearson = cor(o, s), spearman = cor(o, s, method = "spearman"),
                                    amp = sd(s)/sd(o))
rob <- rbindlist(lapply(c("v320", "v323"), function(v) {
  col <- paste0("dTmaxS_", v)
  rbind(
    cbind(version = v, scope = "toutes (53)",            skill(d$dTmaxS_obs, d[[col]])),
    cbind(version = v, scope = "sans trouées (44)",      d[grp == "buffer", skill(dTmaxS_obs, get(col))]))
}))
rob[, (c("pearson","spearman","amp")) := lapply(.SD, round, 2), .SDcols = c("pearson","spearman","amp")]
fwrite(rob, file.path(TBL, "tab_v320_v323_robustness.csv"))
cat("=== #3 robustness: ΔTmax skill with vs without the 9 amplifying gaps ===\n"); print(rob)

# ============================== #4 Taylor diagram =============================
png(file.path(MED, "image_taylor.png"), width = 1500, height = 1500, res = 200)
par(mar = c(4, 4, 4, 2))
pchs <- c(15, 16, 17, 18, 8); names(pchs) <- MET$lab
first <- TRUE
for (i in seq_len(nrow(MET))) {
  o <- d[[MET$obs[i]]]
  for (v in c("v0", "v3")) {
    s   <- d[[MET[[v]][i]]]
    col <- if (v == "v0") SRC[["v3.2.0"]] else SRC[["v3.2.3"]]
    taylor.diagram(o, s, add = !first, normalize = TRUE, col = col, pch = pchs[i],
                   pcex = 1.6, main = "4. Diagramme de Taylor (normalisé) — v3.2.0 vs v3.2.3",
                   sd.arcs = TRUE, ref.sd = TRUE)
    first <- FALSE
  }
}
legend("topright", legend = MET$lab, pch = pchs, pt.cex = 1.3, bty = "n", title = "métrique", cex = .9)
legend("right", legend = c("v3.2.0", "v3.2.3"), col = c(SRC[["v3.2.0"]], SRC[["v3.2.3"]]),
       pch = 19, bty = "n", title = "version", cex = .9)
dev.off()
cat(sprintf("DONE -> image_taylor.png\n"))

# ============================== #5 divergence predictor =======================
d[, absdiff := abs(dTmaxS_v323 - dTmaxS_v320)]
PRED <- c("LAI_ALS", "fCover", "Hmax", "VCI")
z <- copy(d); z[, (PRED) := lapply(.SD, scale), .SDcols = PRED]
lm_std <- lm(reformulate(PRED, "absdiff"), data = z)
co <- summary(lm_std)$coefficients
beta <- data.table(trait = rownames(co), b = co[, 1], p = co[, 4])[trait != "(Intercept)"]
beta[, trait := factor(trait, levels = trait[order(abs(b))])]
cat("\n=== #5 standardized predictors of |Δ dTmaxS| (R2 =", round(summary(lm_std)$r.squared, 2), ") ===\n"); print(beta)
top <- as.character(beta[which.max(abs(b)), trait])

pa <- ggplot(beta, aes(b, trait, fill = b > 0)) +
  geom_col(width = .6) + geom_vline(xintercept = 0, colour = "grey50") +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
  labs(x = "coefficient standardisé", y = NULL,
       title = sprintf("Prédicteurs de la divergence (R² = %.2f)", summary(lm_std)$r.squared)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
pb <- ggplot(d, aes(.data[[top]], absdiff)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = .7) +
  geom_point(aes(fill = P), shape = 21, colour = "black", size = 2.3, stroke = .3) +
  scale_fill_manual(values = PAL, name = "Archétype") +
  labs(x = top, y = "|ΔTmax v3.2.3 − v3.2.0| (°C)", title = sprintf("Divergence vs %s", top)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
f5 <- (pa | pb) +
  plot_annotation(title = "5. Où les deux binaires divergent : prédicteurs structuraux de |v3.2.3 − v3.2.0|",
                  subtitle = "La divergence se concentre dans les peuplements ouverts / peu couvrants — la signature des trouées où le couplage agit.",
                  theme = theme(plot.title = element_text(face = "bold"),
                                plot.subtitle = element_text(size = 8.3, colour = "grey35")))
ggsave(file.path(MED, "image_divergence.png"), f5, width = 11, height = 4.4, dpi = 200, bg = "white")
cat("\nALL DONE — 4 figures + 1 table written.\n")
