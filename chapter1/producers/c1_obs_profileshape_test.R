# ==============================================================================
# Does OBSERVED buffering respond to profile top-heaviness once leaf area and cover are
# controlled? Partial correlation and added R2, overall and within the dense P4.
#
# Reads : in_files/lad_z05/Blois_lad_z05_r25.csv                   (top-heaviness)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv   (stage A2)
# Writes: outputs/figures_pipeline/annex/fig_obs_profileshape.{png,pdf}
#         out_files/Chapter1/tables/tab_obs_profileshape.csv
#   Rscript scripts/c1_obs_profileshape_test.R
# ==============================================================================
# Model-free test that the VERTICAL PROFILE SHAPE adds no buffering information beyond
# leaf quantity, using an INTERPRETABLE shape metric (not VCI, not FPC scores):
#   top-heaviness = relative height of the LAD centroid = (sum z*LAD / sum LAD) / Hmax.
# Higher = foliage concentrated high (the H2 "top-heavy" case).
# Test: does top-heaviness explain observed buffering (ΔTmax) AFTER controlling for LAI
# and fractional cover? Partial correlation + added-R², overall and within the dense P4.
suppressMessages({ library(data.table); library(ggplot2) })
# partial cor(y, x | Z) = cor of residuals of y~Z and x~Z (no ppcor dependency)
#' Partial correlation of y and x given Z, as the correlation of their residuals.
#' Kept manual so the script carries no ppcor dependency.
#' @param y response vector
#' @param x predictor of interest
#' @param Z data.frame of variables to control for
#' @return named numeric c(est, p)
pcor_manual <- function(y, x, Z){ ry <- residuals(lm(y ~ ., data=cbind(y=y, Z))); rx <- residuals(lm(x ~ ., data=cbind(x=x, Z))); ct <- cor.test(ry, rx); c(est=unname(ct$estimate), p=ct$p.value) }
source("scripts/_article_style.R")
DZ <- 0.5
lad <- fread("in_files/lad_z05/Blois_lad_z05_r25.csv")           # lidR profiles (shape unaffected by scan-angle scalar)
lyr <- grep("^LAD_Layer_", names(lad), value = TRUE); h <- as.numeric(sub("LAD_Layer_", "", lyr))
Lmat <- as.matrix(lad[, ..lyr]); Lmat[!is.finite(Lmat)] <- 0
lad[, LAIp   := rowSums(Lmat) * DZ]
lad[, cen    := (Lmat %*% h)[,1] / pmax(rowSums(Lmat), 1e-9)]     # LAD centroid height (m)
lad[, topheavy := cen / pmax(Hmax, 1e-9)]                          # relative centroid height (0-1)
# fraction of leaf area above mid-canopy (per plot)
lad[, frac_high := vapply(seq_len(.N), function(i){ hm <- Hmax[i]; w <- h >= 0.5*hm; s <- sum(Lmat[i,]); if(s<=0) NA_real_ else sum(Lmat[i, w])/s }, numeric(1))]

oc <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")  # dTmax_obs, LAI, fCover, cluster
D <- merge(lad[, .(id_plot, topheavy, frac_high, LAIp, Hmax)], oc, by = "id_plot")
D <- D[is.finite(dTmax_obs) & is.finite(topheavy) & is.finite(LAI) & is.finite(fCover)]
cat(sprintf("n = %d plots ; top-heaviness range [%.2f, %.2f]\n", nrow(D), min(D$topheavy), max(D$topheavy)))

#' Print the added R-squared, partial correlation and coefficient for one subset
#' @param dt subset of the 53-logger table
#' @param tag label printed in the heading
#' @return invisible NULL; called for its printed output
report <- function(dt, tag){
  if (nrow(dt) < 6) { cat(sprintf("\n[%s] n=%d too small\n", tag, nrow(dt))); return(invisible()) }
  m0 <- lm(dTmax_obs ~ LAI + fCover, dt); m1 <- lm(dTmax_obs ~ LAI + fCover + topheavy, dt)
  r0 <- summary(m0)$r.squared; r1 <- summary(m1)$r.squared
  pc <- tryCatch(pcor_manual(dt$dTmax_obs, dt$topheavy, dt[, .(LAI, fCover)]), error=function(e) c(est=NA,p=NA))
  cat(sprintf("\n=== %s (n=%d) ===\n", tag, nrow(dt)))
  cat(sprintf("  simple cor(buffering ΔTmax, top-heaviness) = %.3f (p=%.3f)\n",
              cor(dt$dTmax_obs, dt$topheavy), cor.test(dt$dTmax_obs, dt$topheavy)$p.value))
  cat(sprintf("  R²(ΔTmax~LAI+fCover) = %.3f ; +top-heaviness = %.3f ; ΔR² = %.3f\n", r0, r1, r1-r0))
  cat(sprintf("  PARTIAL cor(ΔTmax, top-heaviness | LAI, fCover) = %.3f (p=%.3f)\n", pc["est"], pc["p"]))
  cat(sprintf("  top-heaviness coef in full model: %+.3f (p=%.3f)\n",
              coef(summary(m1))["topheavy","Estimate"], coef(summary(m1))["topheavy","Pr(>|t|)"]))
}
report(D, "ALL plots")
# dense subset (top-third LAI) where leaf quantity varies least = the sharpest shape test
report(D[LAI >= quantile(LAI, 2/3)], "dense third (high LAI)")

# --- added-variable figure: buffering residual vs top-heaviness residual (after LAI+fCover) ---
D[, res_buf := residuals(lm(dTmax_obs ~ LAI + fCover, D))]
D[, res_top := residuals(lm(topheavy ~ LAI + fCover, D))]
pcr <- cor(D$res_buf, D$res_top)
p <- ggplot(D, aes(res_top, res_buf)) +
  geom_hline(yintercept=0, colour="grey85") + geom_vline(xintercept=0, colour="grey85") +
  geom_smooth(method="lm", se=TRUE, colour="grey30", fill="grey85", linewidth=0.7) +
  geom_point(size=2.4, colour="#1A9850") +
  labs(x="Top-heaviness residual  (relative LAD-centroid height | LAI, fCover)",
       y="Buffering residual  (observed ΔTmax | LAI, fCover)",
       subtitle=sprintf("partial r = %.2f (n = %d) after controlling for LAI and fCover", pcr, nrow(D))) +
  theme_article(12)
OUT <- "outputs/figures_pipeline/annex"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(OUT,"fig_obs_profileshape.png"), p, width=7.5, height=5.5, dpi=300, bg="white")
ggsave(file.path(OUT,"fig_obs_profileshape.pdf"), p, width=7.5, height=5.5, device=cairo_pdf)
fwrite(D[, .(id_plot, dTmax_obs, LAI, fCover, topheavy, frac_high)], "out_files/Chapter1/tables/tab_obs_profileshape.csv")
cat("\nSaved fig_obs_profileshape + table\n")
