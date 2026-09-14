# ==============================================================================
# Verify the attribution ranking is UNCHANGED after standardizing nc_shapley2x[_global]
# from the Blois binary (md5 7ec7) to the validated legacy v3.2.0 (md5 5307).
# Compares per-cluster mean|phi| (Fig 2 source = shapley_parts) and global mean|phi|
# (Fig 3 source = shapley_parts_global) between the backed-up Blois parts and the
# freshly regenerated legacy parts.
#   Rscript c1_binary_ranking_compare.R
# Out: out_files/Chapter1/tables/tab_binary_ranking_compare.csv
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); source("R/cluster_relabel.R") })
Fv <- c("LAI","Hmax","fCover","LAD")
read_parts <- function(dir) {
  f <- list.files(dir, "part_.*csv$", full.names = TRUE)
  if (!length(f)) return(NULL)
  rbindlist(lapply(f, fread), fill = TRUE)
}
meanabs <- function(S, by_cluster) {
  if (is.null(S)) return(NULL)
  if (!"Cluster" %in% names(S)) by_cluster <- FALSE
  agg <- function(d) as.list(sapply(Fv, function(v) mean(abs(d[[v]]), na.rm = TRUE)))
  if (by_cluster) S[, c(agg(.SD)), by = Cluster, .SDcols = Fv][order(Cluster)]
  else            data.table(Cluster = 0L, as.data.table(agg(S)))
}
rank_str <- function(row) paste(names(sort(unlist(row[, ..Fv]), decreasing = TRUE)), collapse = ">")

compare <- function(legacy_dir, blois_dir, by_cluster, tag) {
  L <- meanabs(read_parts(legacy_dir), by_cluster); B <- meanabs(read_parts(blois_dir), by_cluster)
  if (is.null(L) || is.null(B)) { cat(sprintf("[%s] missing parts (legacy or blois) — skip\n", tag)); return(NULL) }
  out <- rbindlist(lapply(seq_len(nrow(L)), function(i) {
    cl <- L$Cluster[i]
    data.table(set = tag, Cluster = cl,
               rank_legacy = rank_str(L[i]), rank_blois = rank_str(B[B$Cluster == cl]),
               LAI_leg = round(L[i]$LAI,3), LAI_blo = round(B[B$Cluster==cl]$LAI,3),
               fCover_leg = round(L[i]$fCover,3), fCover_blo = round(B[B$Cluster==cl]$fCover,3),
               LAD_leg = round(L[i]$LAD,3), LAD_blo = round(B[B$Cluster==cl]$LAD,3))
  }))
  out[, rank_match := rank_legacy == rank_blois]
  out
}
res <- rbind(
  compare("out_files/Chapter1/tables/shapley_parts",         "out_files/Chapter1/tables/shapley_parts_blois_bak",        TRUE,  "per-cluster"),
  compare("out_files/Chapter1/tables/shapley_parts_global",  "out_files/Chapter1/tables/shapley_parts_global_blois_bak", FALSE, "global"),
  fill = TRUE)
fwrite(res, "out_files/Chapter1/tables/tab_binary_ranking_compare.csv")
cat("\n=== ranking: legacy vs Blois (mean|phi|) ===\n"); print(res)
cat(sprintf("\nAll rankings match: %s\n", all(res$rank_match, na.rm = TRUE)))
