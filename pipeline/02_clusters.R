# ==============================================================================
# PIPELINE STAGE 02 — Clusters (FROZEN) + archetypes.
# CONSTRAINT #1: the P1–P4 typology and cLHS sample are a FROZEN INPUT
# (clhs_sample_floor05_v2.rds). All validated results depend on this exact
# assignment, so this stage does NOT re-cluster by default. Set RECLUSTER=TRUE
# to re-run k-means with a pinned seed AND assert membership matches the cache.
# Output: OUT_DATA/archetypes.rds
# ==============================================================================

if (!exists("PIPE")) source(here::here("pipeline/00_config.R"))
RECLUSTER <- isTRUE(as.logical(Sys.getenv("PIPE_RECLUSTER", "FALSE")))
cli_h1("STAGE 02 — clusters (frozen) + archetypes")

stopifnot(file.exists(PIPE$CLUSTER_SAMPLE))
df_sample <- as.data.table(readRDS(PIPE$CLUSTER_SAMPLE))
pipe_assert("Cluster" %in% names(df_sample), "Frozen cLHS sample has a Cluster column")
cli_alert("Frozen sample: {nrow(df_sample)} plots ; clusters {paste(sort(unique(df_sample$Cluster)), collapse=',')}")
print(df_sample[, .(.N, LAI=round(mean(LAI),2), Hmax=round(mean(Hmax),1),
                    fCover=round(mean(fCover),2)), by=Cluster][order(Cluster)])

if (RECLUSTER) {
  cli_alert_warning("RECLUSTER=TRUE — re-running k-means (seed {PIPE$KMEANS_SEED}) and checking against cache")
  set.seed(PIPE$KMEANS_SEED)
  vars <- c("LAI","Hmax","fCover")
  sc <- scale(as.data.frame(df_sample)[, vars])
  km <- kmeans(sc, centers = length(unique(df_sample$Cluster)), nstart = 50)
  agree <- mclust::adjustedRandIndex(km$cluster, as.integer(factor(df_sample$Cluster)))
  pipe_assert(agree > 0.95, sprintf("Re-clustered membership matches cache (ARI=%.3f)", agree))
}

# archetypes (one synthetic centroid plot per cluster) — mean-shape representation
arch <- make_synthetic_archetypes(as.data.frame(df_sample))
saveRDS(arch, file.path(PIPE$OUT_DATA, "archetypes.rds"))
cli_alert("Archetypes:"); print(as.data.table(arch)[, .(Cluster, LAI=round(LAI,2), Hmax=round(Hmax,1), fCover=round(fCover,2))])
cli_alert_success("Stage 02 done (clusters frozen; archetypes cached).")
