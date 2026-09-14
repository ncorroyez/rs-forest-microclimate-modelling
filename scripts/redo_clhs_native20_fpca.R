# ==============================================================================
# CORRECTED native-20 m cLHS: clusters + cLHS on the SIX metrics main.R intends
# (LAI, Hmax, fCover + FPC1..n of the height-normalised LAD profile), which
# redo_clhs_native20.R wrongly reduced to the 3 scalars. Writes to NEW files so
# the current 3-metric sample is untouched pending a GO/NO-GO decision.
# Output: out_files/Sensitivity_Analysis/clhs_sample_native20_fpca{,_floor05}.rds
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(cli); library(terra)
  library(sf); library(dplyr); library(tidyr); library(purrr); library(clhs)
  library(fda); library(ggplot2) })
source(here::here("R/config.R")); source(here::here("R/io.R")); source(here::here("R/forest.R"))
source(here::here("R/fpca.R")); source(here::here("R/cluster_relabel.R"))
IN <- here::here("in_files_native20"); sa <- here::here("out_files/Sensitivity_Analysis")

cli_h1("Native-20 m cLHS rebuild — 6 metrics (FPCA)")
rasters <- load_lidar_rasters(IN, agg_factor = 1)
fd <- build_forest_dataframe(rasters$stack)
df_forest <- fd$df
cli_alert("forest pixels: {nrow(df_forest)}")

# ---- FPCA on the full forest (main.R recipe) ----
fpca_res <- compute_fpca(fd$mat_lad, df_forest$Hmax, fd$z_breaks, n_harm = CFG$n_fpc_total)
n_sel    <- select_n_fpc(fpca_res$varprop, min_marginal = CFG$fpc_min_marginal)
fpc_cols <- paste0("FPC", seq_len(n_sel))
fpc_mat  <- fpca_res$fpca$scores[, seq_len(n_sel), drop = FALSE]; colnames(fpc_mat) <- fpc_cols
df_forest <- bind_cols(df_forest, as.data.frame(fpc_mat))
cluster_vars <- c("LAI", "Hmax", "fCover", fpc_cols)
cli_alert("FPC varprop: {paste(round(100*fpca_res$varprop,1),collapse=' ')} | retained: {n_sel}")

# ---- k-means on the 6 metrics, relabel by ascending cluster-mean LAI ----
set.seed(42)
clust_obj <- label_clusters(df_forest, k = 4, vars = cluster_vars)
df_cl <- as.data.table(clust_obj$df)
.ord    <- df_cl[, .(mlai = mean(LAI, na.rm = TRUE)), by = Cluster][order(mlai)]
.target <- c("3", "2", "4", "1")                          # ascending-LAI rank -> raw int that CLUSTER_RELABEL sends to P1..P4
.remap  <- setNames(.target, as.character(.ord$Cluster))
df_cl[, Cluster := factor(.remap[as.character(Cluster)], levels = c("1","2","3","4"))]
df_cl <- as.data.frame(df_cl)
cli_alert("native cluster LAI order (asc): {paste(round(.ord$mlai,2), collapse=' < ')}")
print(table(relabel_cluster(df_cl$Cluster)))

# ---- cLHS 100/cluster on the SAME 6 metrics ----
df_sample <- sample_clhs_per_cluster(df_cl, n_per_cluster = 100, vars = cluster_vars, iter = 10000)
cli_alert("cLHS: {nrow(df_sample)} plots")

# ---- per-plot scan-angle secant + global LAI scan-angle correction (as 3-metric build) ----
sec_r <- tryCatch(rast(file.path(IN, "sec_theta_res_10_m.tif")), error = function(e) NULL)
if (!is.null(sec_r)) {
  s <- terra::extract(sec_r, as.matrix(df_sample[, c("x","y")]))[,1]
  s[!is.finite(s)] <- median(s, na.rm = TRUE); df_sample$sec_theta <- s
} else df_sample$sec_theta <- 1
gsec <- mean(df_sample$sec_theta, na.rm = TRUE)
df_sample$LAI <- df_sample$LAI / gsec
cli_alert("applied global scan-angle factor {round(gsec,3)} to LAI (mean now {round(mean(df_sample$LAI),2)})")

saveRDS(df_sample, file.path(sa, "clhs_sample_native20_fpca.rds"))
df_floor <- df_sample; df_floor$fCover[df_floor$fCover < 0.5] <- 0.5
saveRDS(df_floor, file.path(sa, "clhs_sample_native20_fpca_floor05.rds"))
cli_alert_success("saved clhs_sample_native20_fpca + _floor05 rds")

# ---- compare to the 3-metric sample currently used by the article ----
old <- as.data.table(readRDS(file.path(sa, "clhs_sample_native20_floor05.rds")))
new <- as.data.table(df_floor)
cat("\nper-archetype summary (3-metric OLD | 6-metric NEW):\n")
print(old[, .(.N, LAI=round(mean(LAI),2), Hmax=round(mean(Hmax),1), fCover=round(mean(fCover),2), VCI=round(mean(VCI),2)), by=relabel_cluster(Cluster)][order(V1)])
print(new[, .(.N, LAI=round(mean(LAI),2), Hmax=round(mean(Hmax),1), fCover=round(mean(fCover),2), VCI=round(mean(VCI),2)), by=relabel_cluster(Cluster)][order(V1)])

# ---- promotion to the canonical filename -------------------------------------
# Added 2026-07-31. The chapter reads `clhs_sample_native20_floor05.rds`, and that file is
# byte-identical to the `_fpca_floor05` one written above (md5 4bfd9d8a...), i.e. the 6-metric
# FPCA sample WAS promoted over the earlier 3-metric one. No script performed the copy, so the
# canonical filename misdescribed its own recipe and the comparison printed just above reads
# "OLD vs NEW" on two identical tables. Do the promotion here, idempotently, and refuse to
# overwrite a canonical sample that differs: that would silently change every downstream run.
canon <- file.path(sa, "clhs_sample_native20_floor05.rds")
srcf  <- file.path(sa, "clhs_sample_native20_fpca_floor05.rds")
if (!file.exists(canon)) {
  file.copy(srcf, canon); cli_alert_success("promoted {basename(srcf)} -> {basename(canon)}")
} else if (identical(tools::md5sum(srcf)[[1]], tools::md5sum(canon)[[1]])) {
  cli_alert_info("canonical sample already identical to the FPCA sample, nothing to promote")
} else {
  stop("REFUSING to promote: ", basename(canon), " differs from ", basename(srcf), ".\n",
       "  Every native20 simulation on disk was run from the canonical file. Overwriting it\n",
       "  would silently invalidate them. Reconcile deliberately, then re-run this script.")
}
