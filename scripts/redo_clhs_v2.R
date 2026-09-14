# ==============================================================================
# Rebuild the cLHS sample with the Hmax-MAX aggregation (consistent extraction).
#
# The legacy clhs_sample.rds was built when load_lidar_rasters() used MEAN
# aggregation for ALL variables including Hmax. With the corrected
# load_lidar_rasters() that uses MAX for Hmax and MEAN for the others, the
# forest dataframe shifts in Hmax → k-means clusters re-partition → cLHS
# sampling needs to be redone.
#
# Output :
#   out_files/Sensitivity_Analysis/clhs_sample_v2.rds          (new 400 plots)
#   out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds  (same with fCover<0.5 floored)
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(terra); library(sf); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/forest.R"))

cli_h1("Rebuild cLHS sample with Hmax-MAX aggregation")

# 1. Load LiDAR with new aggregation rule
rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
cli_alert("New stack: {nlyr(rasters$stack)} layers, resolution {res(rasters$stack)[1]} m")
cli_alert("Hmax stats : mean={round(global(rasters$metrics[['Hmax']], 'mean', na.rm=TRUE)[1,1], 2)} m, "
            ~ "max={round(global(rasters$metrics[['Hmax']], 'max', na.rm=TRUE)[1,1], 2)} m")

# 2. Build forest dataframe
forest_obj <- build_forest_dataframe(rasters$stack)
df_forest <- forest_obj$df
cli_alert("Forest dataframe : {nrow(df_forest)} valid pixels")

# 3. K-means clustering (same vars and k as legacy)
clust_obj <- label_clusters(df_forest, k = 4,
                                vars = c("LAI", "Hmax", "fCover"))
df_forest_cl <- clust_obj$df
cli_alert("k-means clustering done. Cluster sizes :")
print(table(df_forest_cl$Cluster))

# 4. Cluster centroids (with new Hmax)
centroids <- df_forest_cl %>%
  group_by(Cluster) %>%
  summarise(LAI = mean(LAI, na.rm = TRUE),
             Hmax = mean(Hmax, na.rm = TRUE),
             fCover = mean(fCover, na.rm = TRUE),
             n_plots = n(), .groups = "drop")
cli_alert("New cluster centroids :")
print(centroids)

# 5. cLHS sampling (100 plots per cluster)
df_sample <- sample_clhs_per_cluster(df_forest_cl,
                                          n_per_cluster = 100,
                                          vars = c("LAI", "Hmax", "fCover"),
                                          iter = 10000)
cli_alert("cLHS sample : {nrow(df_sample)} plots (~100/cluster)")
cli_alert("Sample Hmax stats : mean={round(mean(df_sample$Hmax), 2)} m, "
            ~ "range [{round(min(df_sample$Hmax), 2)}, {round(max(df_sample$Hmax), 2)}]")

# 6. Save (don't overwrite legacy)
sa_dir <- here::here("out_files/Sensitivity_Analysis")
dir.create(sa_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(df_sample, file.path(sa_dir, "clhs_sample_v2.rds"))

# Floor05 version (cap fCover at 0.5 for plots with very low cover)
df_sample_floor <- df_sample
df_sample_floor$fCover[df_sample_floor$fCover < 0.5] <- 0.5
saveRDS(df_sample_floor, file.path(sa_dir, "clhs_sample_floor05_v2.rds"))

cli_alert_success("Saved clhs_sample_v2.rds and clhs_sample_floor05_v2.rds")

# 7. Compare with legacy quickly
legacy <- as.data.table(readRDS(file.path(sa_dir, "clhs_sample.rds")))
cli_h2("Legacy vs new sample comparison")
cat(sprintf("Legacy : Hmax mean=%.2f, LAI mean=%.2f, fCover mean=%.2f\n",
            mean(legacy$Hmax), mean(legacy$LAI), mean(legacy$fCover)))
cat(sprintf("New v2 : Hmax mean=%.2f, LAI mean=%.2f, fCover mean=%.2f\n",
            mean(df_sample$Hmax), mean(df_sample$LAI), mean(df_sample$fCover)))
