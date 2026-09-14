# ==============================================================================
# Canopy-ratio control (Section 4.3): does the canopy ratio of Starck et al. add
# anything to observed buffering beyond leaf quantity?
#
# WHY THIS EXISTS. The manuscript claimed "no raw association with buffering
# (Spearman rho = 0.04, p = 0.76)". That pair was NOT a correlation: it is the
# dR2 (0.04%) and p_CR of the VCI-baseline row of tab_canopy_ratio_added_power.csv.
# The June table's own CR row gives rho = -0.463 (p = 4.8e-4) on log(slope), and on
# the canonical ALIGNED DeltaTmax lineage the association is rho = +0.40 (p = 0.003)
# for buffering. The redundancy conclusion is unchanged and in fact strengthened:
# CR correlates, yet adds ~0.5% beyond leaf quantity and ~0% in the dense third.
#
# Reads : outputs/figs_MEB2026_final/tab_canopy_ratio_perplot.csv  (per-plot CR)
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv (aligned obs)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv     (LAI, fCover)
# Writes: out_files/Chapter1/tables/tab_canopy_ratio_check.csv
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
P <- fread("outputs/figs_MEB2026_final/tab_canopy_ratio_perplot.csv")[, .(id_plot, CR)]
V <- fread("out_files/Chapter1/tables/tab_hobo_native20_validation.csv")[, .(id_plot, obs_dt)]
C <- fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[, .(id_plot, LAI, fCover)]
D <- merge(merge(P, V, by = "id_plot"), C, by = "id_plot")
D <- D[is.finite(CR) & is.finite(obs_dt) & is.finite(LAI) & is.finite(fCover)][, y := -obs_dt]
#' R-squared of a linear model, on the merged 53-logger table by default
#' @param f model formula
#' @param d data (defaults to D, the 53-logger table)
#' @return numeric R-squared
r2 <- function(f, d = D) summary(lm(f, d))$r.squared
ct <- suppressWarnings(cor.test(D$y, D$CR, method = "spearman"))
dn <- D[LAI >= quantile(LAI, 2/3)]
out <- data.table(
  quantity = c("Spearman rho (CR vs buffering)", "p", "n",
               "dR2 beyond LAI", "dR2 beyond LAI+fCover", "dR2 beyond LAI+fCover, dense third", "n dense"),
  value = c(unname(ct$estimate), ct$p.value, nrow(D),
            r2(y ~ LAI + CR) - r2(y ~ LAI),
            r2(y ~ LAI + fCover + CR) - r2(y ~ LAI + fCover),
            r2(y ~ LAI + fCover + CR, dn) - r2(y ~ LAI + fCover, dn), nrow(dn)))
fwrite(out, "out_files/Chapter1/tables/tab_canopy_ratio_check.csv")
cat(sprintf("CR vs observed buffering: rho = %+.3f (p = %.4f, n = %d)\n", ct$estimate, ct$p.value, nrow(D)))
cat(sprintf("added R2 beyond LAI %.2f%% | beyond LAI+fCover %.2f%% | dense third %.2f%% (n = %d)\n",
  100*(r2(y~LAI+CR)-r2(y~LAI)), 100*(r2(y~LAI+fCover+CR)-r2(y~LAI+fCover)),
  100*(r2(y~LAI+fCover+CR,dn)-r2(y~LAI+fCover,dn)), nrow(dn)))
cat("DONE\n")
