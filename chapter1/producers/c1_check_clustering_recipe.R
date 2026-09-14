# Which recipe produced the FROZEN typology: {LAI,Hmax,fCover} vs +FPC1-3 ?
# Recluster the forest both ways (label_clusters recipe: seed42, scale, k=4, nstart=25),
# match the 400 frozen cLHS sample by (x,y), compare to its frozen Cluster.
suppressMessages({ library(data.table); library(mclust) })
fc <- readRDS("outputs/figures_pipeline_z05/data/forest_fpca.rds")
fo <- as.data.table(fc$forest)
sco <- fc$fpca$fpca$scores                             # forest FPC scores (pca.fd, rows aligned to fo)
stopifnot(nrow(sco) == nrow(fo))
fo[, `:=`(FPC1=sco[,1], FPC2=sco[,2], FPC3=sco[,3])]

reclust <- function(vars) {
  d <- fo[complete.cases(fo[, ..vars])]
  set.seed(42); km <- kmeans(scale(as.matrix(d[, ..vars])), centers=4, iter.max=100, nstart=25)
  d[, cl := km$cluster]; d[, .(x,y,cl)]
}
r3  <- reclust(c("LAI","Hmax","fCover"))
r4  <- reclust(c("LAI","Hmax","fCover","FPC1","FPC2","FPC3"))

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))[, .(x,y,frozen=as.integer(as.character(Cluster)))]
m3 <- merge(samp, r3, by=c("x","y")); m4 <- merge(samp, r4, by=c("x","y"))
cat(sprintf("sample matched to forest: 3-var %d/400, 4-var %d/400\n", nrow(m3), nrow(m4)))

agree <- function(a, b) {                              # max agreement over label permutations + ARI
  library(combinat); best <- 0
  for (p in permn(1:4)) { mapped <- p[b]; best <- max(best, mean(mapped==a)) }
  list(pct=best, ari=adjustedRandIndex(a,b))
}
a3 <- agree(m3$frozen, m3$cl); a4 <- agree(m4$frozen, m4$cl)
cat(sprintf("\n=== agreement avec la typologie GELÉE (sur les %d placettes) ===\n", nrow(m4)))
cat(sprintf("  {LAI,Hmax,fCover}        : %.1f%% identiques  (ARI=%.3f)\n", 100*a3$pct, a3$ari))
cat(sprintf("  {LAI,Hmax,fCover,FPC1-3} : %.1f%% identiques  (ARI=%.3f)\n", 100*a4$pct, a4$ari))
cat("\n→ la recette qui reproduit le gelé = celle réellement utilisée.\n")
