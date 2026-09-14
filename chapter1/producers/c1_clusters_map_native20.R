# ==============================================================================
# Native-20 m wall-to-wall archetype raster. Reproduces EXACTLY the partition of
# redo_clhs_native20.R (k-means on LAI/Hmax/fCover of the native forest grid,
# seed 42, relabel by ascending cluster-mean LAI) and writes it as a GeoTIFF of
# the P index (1..4). Because the recipe/seed/data are identical, the 400 cLHS
# plots inherit the same labels as clhs_sample_native20_floor05.rds -> the map is
# 100% consistent with the native archetypes. Nothing baseline is touched.
# Output: out_files/clusters/blois_clusters_P1P4_native20.tif
# ==============================================================================
suppressMessages({ library(data.table); library(here); library(terra); library(dplyr)
  library(fda); library(ggplot2) })
source(here::here("R/config.R")); source(here::here("R/io.R"))
source(here::here("R/forest.R")); source(here::here("R/fpca.R")); source(here::here("R/cluster_relabel.R"))

IN <- here::here("in_files_native20")
rasters <- load_lidar_rasters(IN, agg_factor = 1)                 # native 20 m, no aggregation
forest  <- build_forest_dataframe(rasters$stack)
df      <- as.data.table(forest$df)

# 6-metric recipe (main.R / redo_clhs_native20_fpca.R): LAI/Hmax/fCover + FPC1..n
fpca_res <- compute_fpca(forest$mat_lad, df$Hmax, forest$z_breaks, n_harm = CFG$n_fpc_total)
n_sel    <- select_n_fpc(fpca_res$varprop, min_marginal = CFG$fpc_min_marginal)
fpc_cols <- paste0("FPC", seq_len(n_sel))
fpc_mat  <- fpca_res$fpca$scores[, seq_len(n_sel), drop = FALSE]; colnames(fpc_mat) <- fpc_cols
df       <- cbind(df, as.data.table(fpc_mat))

set.seed(42)                                                      # same as redo_clhs_native20_fpca.R
clu   <- label_clusters(as.data.frame(df), k = 4, vars = c("LAI", "Hmax", "fCover", fpc_cols))
df_cl <- as.data.table(clu$df)

# identical ascending-LAI remap so the raw k-means integers map to P1..P4 as in the sample
.ord    <- df_cl[, .(mlai = mean(LAI, na.rm = TRUE)), by = Cluster][order(mlai)]
.target <- c("3", "2", "4", "1")
.remap  <- setNames(.target, as.character(.ord$Cluster))
df_cl[, Cluster := factor(.remap[as.character(Cluster)], levels = c("1","2","3","4"))]
df_cl[, P := relabel_cluster(Cluster)]
df_cl[, Pidx := as.integer(P)]
cat("native cluster LAI order (asc):", paste(round(.ord$mlai,2), collapse=" < "), "\n")
print(df_cl[!is.na(Pidx), .N, by = P][order(P)])

# write P index (1..4) to a raster on the native grid
template <- rasters$stack[[1]]
r <- terra::rast(template); terra::values(r) <- NA_integer_
cells <- terra::cellFromXY(r, as.matrix(df_cl[, .(x, y)]))
terra::values(r)[cells] <- df_cl$Pidx
names(r) <- "Profile"
dir.create(here::here("out_files/clusters"), showWarnings = FALSE, recursive = TRUE)
terra::writeRaster(r, here::here("out_files/clusters/blois_clusters_P1P4_native20.tif"),
                   overwrite = TRUE, datatype = "INT1U", NAflag = 255)
cat("saved out_files/clusters/blois_clusters_P1P4_native20.tif\n")

# ---- VALIDATION: the 400 cLHS plots must inherit the same labels as the raster --
samp <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds")))
samp[, Psamp := relabel_cluster(Cluster)]
samp[, Prast := c("P1","P2","P3","P4")[terra::extract(r, as.matrix(samp[, .(x, y)]))[,1]]]
agree <- samp[!is.na(Prast), mean(Psamp == Prast)]
cat(sprintf("VALIDATION: %.1f%% of the 400 cLHS plots match the raster label\n", 100*agree))
print(table(sample = samp$Psamp, raster = samp$Prast))
