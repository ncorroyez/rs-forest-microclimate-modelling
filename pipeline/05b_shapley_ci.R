# ==============================================================================
# PIPELINE STAGE 05b — Bootstrap robustness of the Shapley attribution (TABLE).
# 95% percentile CI of the MEDIAN φ per (metric × period × Cluster × trait),
# resampling PLOTS with replacement (B=2000, seeded). Frames the attribution as
# "robust to plot resampling" (CI excludes 0), NOT a population significance test
# — φ are deterministic model outputs over non-independent, stratified plots, and
# the 4 traits are additivity-constrained (Σφ = REF − Base). No figure changes.
# Output: tables/tab_shapley_ci.csv
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
cli_h1("STAGE 05b — Shapley bootstrap CI (table)")

phi <- as.data.table(readRDS(file.path(PIPE$OUT_DATA, "shapley_perplot.rds")))
phi[, Cluster := as.character(Cluster)]
phiC <- rbind(phi, copy(phi)[, Cluster := "All"])          # add pooled "All"
phiC[, Cluster := factor(Cluster, levels = c(paste0("P",1:4), "All"))]

B <- 2000L; SEED <- 42L
boot_median_ci <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 2L) return(list(median = if (n) median(x) else NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, n = n))
  set.seed(SEED)
  meds <- vapply(seq_len(B), function(b) median(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  ci <- quantile(meds, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  list(median = median(x), ci_lo = ci[1], ci_hi = ci[2], n = n)
}

ci <- phiC[, { r <- boot_median_ci(phi)
               .(n = r$n, median = r$median, ci_lo = r$ci_lo, ci_hi = r$ci_hi) },
           by = .(metric, period, Cluster, trait)]
ci[, excludes_zero := is.finite(ci_lo) & (ci_lo > 0 | ci_hi < 0)]
ci[, direction := fifelse(!excludes_zero, "ns",
                   fifelse(median < 0, "buffers", "amplifies"))]
setorder(ci, metric, period, Cluster, median)

# round for the table
num <- c("median","ci_lo","ci_hi")
ci_out <- copy(ci)[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(ci_out, file.path(PIPE$OUT_TAB, "tab_shapley_ci.csv"))

cli_alert_success("Saved tab_shapley_ci.csv ({nrow(ci)} cells, B={B})")
cli_h2("ΔTmax — robust (CI excludes 0) attributions per cluster")
print(ci_out[metric=="Tmax" & excludes_zero==TRUE,
             .(period, Cluster, trait, median, CI = sprintf("[%.3f, %.3f]", ci_lo, ci_hi), direction)])
cli_alert("Cells with CI excluding 0: {sum(ci$excludes_zero)}/{nrow(ci)} (Tmax {sum(ci$excludes_zero & ci$metric=='Tmax')}, VPD {sum(ci$excludes_zero & ci$metric=='VPD')})")
