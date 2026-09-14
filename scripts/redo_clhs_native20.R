# ==============================================================================
# Stage 3: rebuild the FPCA typology + k-means + cLHS design on the NATIVE 20 m
# rasters (in_files_native20/, agg_factor=1 so aggregation is a no-op). Separate
# output, nothing overwritten. Mirrors redo_clhs_v2.R exactly except the in_dir/agg.
# Output: out_files/Sensitivity_Analysis/clhs_sample_native20{,_floor05}.rds
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(cli)
  library(terra); library(sf); library(dplyr); library(tidyr); library(purrr) })
source(here::here("R/config.R")); source(here::here("R/io.R")); source(here::here("R/forest.R"))
IN <- here::here("in_files_native20"); sa <- here::here("out_files/Sensitivity_Analysis")

cli_h1("Native-20 m cLHS rebuild")
rasters <- load_lidar_rasters(IN, agg_factor = 1)                 # native 20 m, no aggregation
cli_alert("stack: {nlyr(rasters$stack)} layers @ {res(rasters$stack)[1]} m")

forest_obj <- build_forest_dataframe(rasters$stack)
df_forest  <- forest_obj$df
cli_alert("forest pixels: {nrow(df_forest)}")

set.seed(42)                                             # reproducible native k-means
clust_obj <- label_clusters(df_forest, k = 4, vars = c("LAI", "Hmax", "fCover"))
df_cl <- as.data.table(clust_obj$df)
# The k-means integers are arbitrary; CLUSTER_RELABEL (R/cluster_relabel.R) maps raw
# integers -> P by the BASELINE numbering (3->P1,2->P2,4->P3,1->P4). Remap the NATIVE
# integers by ascending cluster-mean LAI so P1=open(low LAI) .. P4=dense(high LAI),
# matching the baseline convention for an apples-to-apples per-archetype comparison.
source(here::here("R/cluster_relabel.R"))
.ord <- df_cl[, .(mlai = mean(LAI, na.rm = TRUE)), by = Cluster][order(mlai)]
.target <- c("3", "2", "4", "1")                          # ascending-LAI rank -> raw integer that CLUSTER_RELABEL sends to P1..P4
.remap <- setNames(.target, as.character(.ord$Cluster))
df_cl[, Cluster := factor(.remap[as.character(Cluster)], levels = c("1","2","3","4"))]
df_cl <- as.data.frame(df_cl)
cli_alert("native cluster LAI order (asc): {paste(round(.ord$mlai,2), collapse=' < ')}")
print(table(relabel_cluster(df_cl$Cluster)))

df_sample <- sample_clhs_per_cluster(df_cl, n_per_cluster = 100,
                                     vars = c("LAI", "Hmax", "fCover"), iter = 10000)
cli_alert("cLHS: {nrow(df_sample)} plots")

# per-plot scan-angle secant from the native sec raster (falls back to global mean)
sec_r <- tryCatch(rast(file.path(IN, "sec_theta_res_10_m.tif")), error = function(e) NULL)
if (!is.null(sec_r)) {
  s <- terra::extract(sec_r, as.matrix(df_sample[, c("x","y")]))[,1]
  s[!is.finite(s)] <- median(s, na.rm = TRUE); df_sample$sec_theta <- s
} else df_sample$sec_theta <- 1
cli_alert("sec_theta: mean={round(mean(df_sample$sec_theta),3)}")

# global scan-angle correction on LAI (matches the baseline cLHS, corrected via a
# global factor AFTER the design; a global rescale does not change the k-means partition)
gsec <- mean(df_sample$sec_theta, na.rm = TRUE)
df_sample$LAI <- df_sample$LAI / gsec
cli_alert("applied global scan-angle factor {round(gsec,3)} to LAI (mean now {round(mean(df_sample$LAI),2)})")

saveRDS(df_sample, file.path(sa, "clhs_sample_native20.rds"))
df_floor <- df_sample; df_floor$fCover[df_floor$fCover < 0.5] <- 0.5
saveRDS(df_floor, file.path(sa, "clhs_sample_native20_floor05.rds"))
cli_alert_success("saved clhs_sample_native20 + _floor05 rds")

# compare to the aggregated design
agg <- as.data.table(readRDS(file.path(sa, "clhs_sample_floor05_v2.rds")))
cat(sprintf("AGG 20m : LAI=%.2f Hmax=%.1f fCover=%.2f\n", mean(agg$LAI), mean(agg$Hmax), mean(agg$fCover)))
cat(sprintf("NAT 20m : LAI=%.2f Hmax=%.1f fCover=%.2f\n", mean(df_floor$LAI), mean(df_floor$Hmax), mean(df_floor$fCover)))
print(as.data.table(df_floor)[, .(.N, LAI=round(mean(LAI),2), Hmax=round(mean(Hmax),1), fCover=round(mean(fCover),2)), by=Cluster][order(Cluster)])
