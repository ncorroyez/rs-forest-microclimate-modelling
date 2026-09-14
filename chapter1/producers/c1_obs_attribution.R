# ==============================================================================
# #2 — Observational attribution: does LAI dominate the BUFFERING in the DATA,
# independently of MuSICA? Regress OBSERVED ΔTmax (53 HOBO) on the 4 LiDAR traits.
# This is the model-free triangulation the committee asked for (no black box).
#
# Method (collinearity-aware): (a) LMG relative importance (Lindeman-Merenda-Gold,
# averages R² over all predictor orderings — the standard with correlated traits);
# (b) random-forest permutation importance as a non-linear cross-check;
# (c) simple + partial Pearson correlations; (d) open (P1+P2) vs closed (P3) split.
#
# HARD CAVEATS (printed): n=53; severe collinearity (all traits r=0.75–0.86, VCI
# 0.82–0.86 with everything); and NO P4 plots in the field (P1=7,P2=15,P3=31,P4=0)
# => the model's dense-P4 LAD co-lead is NOT field-testable here. Triangulation,
# not a clean attribution.
#   Rscript c1_obs_attribution.R
# Out: tab_obs_attribution.csv + FigSh_obs_attribution.png
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(relaimpo); library(randomForest)
})
set.seed(1)
D <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")
TR <- c("LAI", "fCover", "Hmax", "VCI")
# sign convention: dTmax_obs < 0 = buffered; analyze buffering magnitude = -dTmax_obs
D[, buffering := -dTmax_obs]

# ---- (a) LMG relative importance (share of model R²), full sample --------------
lmg_of <- function(d) {
  if (nrow(d) < length(TR) + 3 || any(sapply(d[, ..TR], function(x) sd(x) == 0))) return(NULL)
  f <- as.formula(paste("buffering ~", paste(TR, collapse = " + ")))
  m <- lm(f, data = d)
  ri <- calc.relimp(m, type = "lmg", rela = TRUE)   # rela=TRUE -> shares sum to 1
  c(R2 = summary(m)$r.squared, ri@lmg[TR])
}
full <- lmg_of(D)

# bootstrap CI of LMG shares (full sample)
B <- 2000
bs <- t(replicate(B, { i <- sample(nrow(D), nrow(D), TRUE); v <- lmg_of(D[i]); if (is.null(v)) rep(NA, 5) else v }))
colnames(bs) <- c("R2", TR)
ci <- function(col) quantile(bs[, col], c(.025, .975), na.rm = TRUE)

# ---- (b) RF permutation importance (non-linear cross-check) --------------------
rf <- randomForest(x = as.data.frame(D[, ..TR]), y = D$buffering,
                   importance = TRUE, ntree = 2000)
rfimp <- importance(rf, type = 1)[, 1]                # %IncMSE
rfimp <- rfimp / sum(pmax(rfimp, 0))                  # normalize to shares

# ---- (c) simple & partial Pearson correlations vs buffering --------------------
simple <- sapply(TR, function(v) cor(D$buffering, D[[v]]))
partial <- sapply(TR, function(v) {           # partial r controlling the other 3
  others <- setdiff(TR, v)
  ry <- resid(lm(reformulate(others, "buffering"), data = D))
  rx <- resid(lm(reformulate(others, v),          data = D))
  cor(ry, rx)
})

# NB: per-archetype / gradient LMG is DELIBERATELY NOT reported. At ~14 loggers
# per stratum with predictors collinear at r=0.82–0.86, a 4-predictor LMG split is
# degenerate (returns are unstable, e.g. exact 0.00 = variance fully absorbed by a
# collinear partner). The defensible observational statement is full-sample only.
# The field DOES span dense stands (one-sided LAI up to ~6.4, matching cLHS-P4 ~5.4),
# but trait collinearity prevents isolating the profile's independent role at ANY
# density — so the model's dense co-lead is neither confirmed nor refuted here.

# ---- assemble & report ---------------------------------------------------------
res <- data.table(
  trait        = TR,
  simple_r     = round(simple, 3),
  partial_r    = round(partial, 3),
  LMG_share    = round(full[TR], 3),
  LMG_lo       = round(sapply(TR, function(t) ci(t)[1]), 3),
  LMG_hi       = round(sapply(TR, function(t) ci(t)[2]), 3),
  RF_share     = round(rfimp[TR], 3))
setorder(res, -LMG_share)
fwrite(res, "out_files/Chapter1/tables/tab_obs_attribution.csv")

cat("\n=== #2 OBSERVATIONAL ATTRIBUTION (53 HOBO, model-free) ===\n")
cat(sprintf("Full-sample lm R² = %.2f (n=%d). LMG shares sum to 1.\n", full["R2"], nrow(D)))
print(res)
cat(sprintf("\nLeader (LMG): %s  (%.0f%% of explained variance); quantity (LAI+fCover) = %.0f%%\n",
            res$trait[1], 100 * res$LMG_share[1], 100 * sum(res[trait %in% c("LAI","fCover"), LMG_share])))
cat("\n*** CAVEATS ***\n")
cat(" - n=53; severe collinearity (traits r=0.75–0.86; VCI 0.82–0.86 with all) -> partials unstable,\n")
cat("   the profile's INDEPENDENT role cannot be isolated observationally at any density.\n")
cat(" - field DOES span dense stands (one-sided LAI up to ~6.4 ~ cLHS-P4 ~5.4); per-archetype LMG is\n")
cat("   NOT reported (degenerate at ~14/stratum). Triangulation of quantity-dominance, not attribution.\n")

# ---- figure --------------------------------------------------------------------
res[, trait := factor(trait, levels = rev(trait))]
G <- melt(res[, .(trait, LMG = LMG_share, RF = RF_share)], id.vars = "trait")
pA <- ggplot(res, aes(LMG_share, trait)) +
  geom_col(fill = "#3690C0", width = .6) +
  geom_errorbarh(aes(xmin = LMG_lo, xmax = LMG_hi), height = .2, colour = "grey25") +
  labs(x = "LMG share of explained variance (±95% boot)", y = NULL,
       title = sprintf("Observed buffering ~ LiDAR traits (53 HOBO, R²=%.2f)", full["R2"]),
       subtitle = "Model-free check: LMG relative importance. Quantity (LAI+fCover) dominates. CAVEAT n=53, collinearity 0.75–0.86 limits trait separability.") +
  theme_bw(base_size = 11) + theme(plot.subtitle = element_text(size = 8, colour = "grey35"))
pB <- ggplot(G, aes(value, trait, fill = variable)) +
  geom_col(position = position_dodge(.65), width = .6) +
  scale_fill_manual(values = c(LMG = "#3690C0", RF = "#A6611A"), name = NULL) +
  labs(x = "Importance share", y = NULL, title = "LMG vs random-forest (cross-check)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("out_files/Chapter1/figures/FigSh_obs_attribution.png",
       pA, width = 7.5, height = 3.6, dpi = 200, bg = "white")
ggsave("out_files/Chapter1/figures/FigSh_obs_attribution_crosscheck.png",
       pB, width = 6, height = 3.8, dpi = 200, bg = "white")
cat("DONE -> tab_obs_attribution.csv + FigSh_obs_attribution.png (+_crosscheck)\n")
