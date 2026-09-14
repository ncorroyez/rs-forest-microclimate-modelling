# ==============================================================================
# c1_windcorr_diagnostic.R
#
# Why does removing the per-plot wind correction improve the open-canopy ranking
# while leaving the dense archetypes untouched?
#
# The correction scales the forcing wind by wind_factor(Hmax), a neutral
# log-profile ratio. That factor is nearly constant across tall canopies and
# spreads widely across short ones, so it injects between-plot variance almost
# only where the canopy is low. This script tests whether the plots whose
# simulated dTmax moves most between the two lineages are the plots whose wind
# factor departs most from the sample median.
#
# A positive association supports reading the correction as propagated structural
# uncertainty rather than as a physical signal MuSICA fails to reproduce.
#
# Inputs : /tmp/val_wind.csv, /tmp/val_nowind.csv  (canonical extractor, both lineages)
#          out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv  (Hmax, archetype)
# Output : out_files/Chapter1/tables/tab_windcorr_diagnostic.csv
#   Rscript c1_windcorr_diagnostic.R
# ==============================================================================

suppressPackageStartupMessages(library(data.table))
source("R/wind_correction.R")

W  <- fread("/tmp/val_wind.csv")[,   .(id_plot, sim_wind = sim_dt, obs_dt)]
NW <- fread("/tmp/val_nowind.csv")[, .(id_plot, sim_nowind = sim_dt)]
CL <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,
            .(id_plot = as.character(id_plot), Hmax, P)]

D <- merge(merge(W, NW, by = "id_plot"), CL, by = "id_plot")
stopifnot(nrow(D) == 53)

D[, f := wind_factor(pmax(Hmax, 2))]
D[, f_dev := abs(f - median(f))]              # departure from the sample median
D[, d_sim := sim_nowind - sim_wind]           # how far the simulated value moved
D[, err_wind   := abs(sim_wind   - obs_dt)]   # per-plot absolute error, each lineage
D[, err_nowind := abs(sim_nowind - obs_dt)]
D[, d_err := err_wind - err_nowind]           # > 0 means removing the correction helped

cat(sprintf("\nwind factor: %.2f to %.2f (median %.2f)\n", min(D$f), max(D$f), median(D$f)))

sp <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  sprintf("rho = %+.2f, p = %.4f", ct$estimate, ct$p.value)
}

cat("\n== Does the simulated value move most where the factor is most extreme? ==\n")
cat("  |d_sim| vs |f - median(f)|, all 53 : ", sp(D$f_dev, abs(D$d_sim)), "\n")
cat("\n== Does removing the correction help most where the factor is most extreme? ==\n")
cat("  d_err vs |f - median(f)|, all 53   : ", sp(D$f_dev, D$d_err), "\n")

cat("\n== Per archetype ==\n")
for (p in sort(unique(D$P))) {
  d <- D[P == p]
  cat(sprintf("  %s (n=%2d) f %.2f-%.2f | mean |d_sim| = %.3f degC | mean d_err = %+.3f degC\n",
              p, nrow(d), min(d$f), max(d$f), mean(abs(d$d_sim)), mean(d$d_err)))
}

# Sparse and dense as the validation reports them (median Hmax split is not used;
# the archetypes are the chapter's own grouping).
D[, stratum := ifelse(P == "P4", "dense-P4", "other")]
cat("\n== Movement by canopy height ==\n")
cat("  correlation Hmax vs |d_sim| : ", sp(D$Hmax, abs(D$d_sim)), "\n")

fwrite(D, "out_files/Chapter1/tables/tab_windcorr_diagnostic.csv")
cat("\nwrote out_files/Chapter1/tables/tab_windcorr_diagnostic.csv\n")
