# ==============================================================================
# PIPELINE — run the full spine, each stage in an ISOLATED Rscript process.
# Process isolation avoids library-masking/state accumulation across stages and
# keeps every stage independently reproducible. Stops on first failure.
#
# Usage:  Rscript pipeline/run_all.R
#         PIPE_RECLUSTER=TRUE Rscript pipeline/run_all.R   # opt-in re-clustering
#
# Optional analyses (NOT run here): pipeline/optional/{10_canopy_ratio,
# 11_vertical_profiles,12_height_robustness}.R — launch individually.
#
# ⚠️ SUPERSEDED ATTRIBUTION METHOD (2026-06-20). The Chapter 1 ARTICLE headline
# attribution is the on-manifold trait-perturbation SENSITIVITY analysis, produced
# by the c1_*.R scripts (c1_sensitivity_perplot_chunk*, c1_metrics6_perplot,
# c1_deltadelta_perplot, c1_importance_*) + the article figure set in
# chapter1/ledger/make_article_figure_set.R. The Shapley stages below (05/05b/10/10b) are
# the OLD decomposition; they are retained for provenance and still run end-to-end,
# but their outputs are NOT the article result. To reproduce the article, run the
# c1_* sensitivity chain, not these Shapley stages.
# ==============================================================================

# data/sim/validation backbone (still current) + Shapley stages (SUPERSEDED, see above)
stages <- c("01_data","02_clusters","03_musica","04_extract",
            "05_shapley","05b_shapley_ci",                 # SUPERSEDED (legacy Shapley)
            "06_boxplots","06b_ranking",
            "07_validation","07b_forward_validation","08_typology",
            "09_extra_metrics",
            "10_conditional_shapley","10b_conditional_shapley_lai_fcover")  # SUPERSEDED (legacy Shapley)

t0 <- Sys.time()
for (s in stages) {
  cat(sprintf("\n=============== RUN %s ===============\n", s))
  script <- here::here("pipeline", paste0(s, ".R"))
  code <- system2("Rscript", script, stdout = "", stderr = "")
  if (!identical(code, 0L)) stop(sprintf("Stage %s FAILED (exit %s) — pipeline halted.", s, code))
}
cat(sprintf("\nPIPELINE complete in %.1f min — figures in outputs/figures_pipeline/\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
