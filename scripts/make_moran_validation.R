# ==============================================================================
# Moran's I on the validation fields — CANONICAL native20 lineage.
# Tests whether (i) the OBSERVED per-plot buffering field is spatially
# structured, (ii) the MODELLED field reproduces that structure, and (iii) the
# validation RESIDUAL (obs - model) is spatially autocorrelated (which would
# inflate per-plot inference). Spatial weights: k-nearest-neighbour (k=8),
# row-standardized; significance by 999 Monte-Carlo permutations.
#   Rscript scripts/make_moran_validation.R
# LINEAGE FIX 2026-07-31. This read outputs/figures_pipeline_z05/data/ref_validation.rds,
# i.e. the SUPERSEDED z05 branch, so Appendix C quoted Moran statistics of a validation
# that is not the chapter's. It now reads the same native20 table the rest of Appendix C
# decomposes (tab_hobo_native20_validation.csv), which is also what make_residual_vs_topo.R
# uses, so the two Appendix C analyses finally describe one field.
# Out: out_files/Chapter1/tables/tab_moran_validation.csv
# ==============================================================================
suppressMessages({ library(here); library(data.table); library(spdep) })
OUT <- here::here("out_files/Chapter1/tables"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
KNN <- 8L; NSIM <- 999L
set.seed(42)   # moran.mc is a permutation test: without a seed its p-values are irreproducible

# ---- per-plot observed / modelled / residual + coords (native20 validation) ---
V <- fread(here::here("out_files/Chapter1/tables/tab_hobo_native20_validation.csv"))
res <- V[is.finite(sim_dt) & is.finite(obs_dt),
         .(id_plot, dTmax_sim = sim_dt, dTmax_obs = obs_dt)][
         , residual := dTmax_obs - dTmax_sim]
xy <- fread(here::here("in_files/lad_z05/Blois_lad_z05_r25.csv"))[, .(id_plot, x, y)]
D <- merge(res, xy, by = "id_plot")
stopifnot(nrow(D) >= 50)

# ---- kNN spatial weights (row-standardized) ----------------------------------
coords <- as.matrix(D[, .(x, y)])
lw <- nb2listw(knn2nb(knearneigh(coords, k = KNN)), style = "W")

mtest <- function(z) {
  mt <- moran.mc(z, lw, nsim = NSIM, zero.policy = TRUE)
  data.table(I = round(unname(mt$statistic), 3), p_value = signif(mt$p.value, 3))
}
out <- rbindlist(list(
  cbind(field = "observed ΔTmax",  mtest(D$dTmax_obs)),
  cbind(field = "modelled ΔTmax",  mtest(D$dTmax_sim)),
  cbind(field = "validation residual", mtest(D$residual))
))
out[, `:=`(n = nrow(D), knn = KNN, nsim = NSIM)]
cat(sprintf("=== Moran's I (native20, n=%d, kNN=%d, %d perms) ===\n", nrow(D), KNN, NSIM)); print(out)
fwrite(out, file.path(OUT, "tab_moran_validation.csv"))
cat("DONE ->", file.path(OUT, "tab_moran_validation.csv"), "\n")
