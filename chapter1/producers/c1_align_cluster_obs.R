# ==============================================================================
# Align the OBSERVED columns of tab_hobo_perplot_cluster.csv onto the canonical
# clock (R/dtmax_convention.R, OBS_CLOCK_OFFSET_H = -1).
#
# WHY. That table dates from 2026-07-20, before the clock alignment. Its
# dTmax_obs is the pre-alignment observation (mean -0.730) while the headline
# validation uses the aligned one (mean -1.030); r = 0.992 but the levels differ
# by ~0.3 degC and up to 1.66 degC on one plot. The chapter was therefore running
# TWO observational conventions: the validation aligned, every model-free
# analysis (Fig. 8 nested regression, Fig. E2 top-heaviness, Appendix C residual,
# the Kruskal-Wallis by archetype) not.
#
# WHAT THIS DOES. Overwrites ONLY dTmax_obs and slope_obs, from
# tab_hobo_native20_validation.csv (obs_dt / obs_sl, produced by
# c1_hobo_native20_validation_regen.R with macro_ref_obs). Everything else is
# copied through untouched.
#
# WHY NOT REGENERATE THE WHOLE TABLE. The P column is the canonical archetype
# label, read from outputs/figures_pipeline_z05/data/clusters.rds and assigned by
# ascending LAI (R/cluster_relabel.R). Re-running the clustering risks new raw
# ids and a relabel that would cascade into every figure, the palette and every
# per-archetype number in both manuscripts. The labels stay frozen here.
# Reads :
#         out_files/Chapter1/tables/tab_hobo_native20_validation.csv (stage A1)
# Writes: out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv  (dTmax_obs/slope_obs ONLY)
#         out_files/Chapter1/tables/tab_hobo_perplot_cluster_preAlign_backup.csv (first run only)
#   Rscript c1_align_cluster_obs.R      (opt-in stage: CH1_ALIGN_OBS=TRUE)
# ==============================================================================
suppressPackageStartupMessages(library(data.table))
TAB <- "out_files/Chapter1/tables"
SRC <- file.path(TAB, "tab_hobo_perplot_cluster.csv")
BAK <- file.path(TAB, "tab_hobo_perplot_cluster_preAlign_backup.csv")

D <- fread(SRC)
if (!file.exists(BAK)) { fwrite(D, BAK); cat("backup written:", BAK, "\n") }
V <- fread(file.path(TAB, "tab_hobo_native20_validation.csv"))[, .(id_plot, obs_dt, obs_sl)]

stopifnot(nrow(D) == 53L, setequal(D$id_plot, V$id_plot), "P" %in% names(D))
P_before <- D[order(id_plot), paste(P, collapse = "|")]

M <- merge(D, V, by = "id_plot")
cat(sprintf("dTmax_obs : mean %+.3f -> %+.3f | r = %.4f | max |diff| = %.3f\n",
            mean(M$dTmax_obs), mean(M$obs_dt), cor(M$dTmax_obs, M$obs_dt), max(abs(M$obs_dt - M$dTmax_obs))))
cat(sprintf("slope_obs : mean %+.4f -> %+.4f | r = %.4f\n",
            mean(M$slope_obs, na.rm = TRUE), mean(M$obs_sl, na.rm = TRUE),
            cor(M$slope_obs, M$obs_sl, use = "complete.obs")))

D[V, on = "id_plot", `:=`(dTmax_obs = i.obs_dt, slope_obs = i.obs_sl)]
stopifnot(identical(D[order(id_plot), paste(P, collapse = "|")], P_before))   # labels frozen
fwrite(D, SRC)
cat("labels unchanged; table rewritten on the aligned clock\nDONE\n")
