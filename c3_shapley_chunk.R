# ==============================================================================
# Chapter 3 — per-point cLHS Shapley (per-cluster baseline, 2×LAI), CHUNKED for
# multi-PROCESS parallelism. NOT mclapply/fork: setup_dir() makes its tempdir with
# tempfile() in cwd, and forked workers share RNG/counter state → tempdir collisions
# would clobber each other's MuSICA runs. Separate OS processes each get independent
# tempfile state, so process-level parallelism is the safe way to use all cores.
# Baseline = per-cluster MEAN of every numeric factor (LAI, Hmax, fCover); LAD
# baseline = uniform (no vertical structure). Coalition 0000 = the cluster archetype,
# so each plot's ΔTmax is decomposed vs its own cluster archetype.
# Usage:  Rscript c3_shapley_chunk.R <chunk> <K>     # chunk in 1..K
# Each process handles rows where ((row-1) %% K == chunk-1), serially, and writes
# (incrementally, one row per plot) out_files/Chapter3/tables/shapley_parts/part_<chunk>.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(data.table)
  library(rmusica); library(musica.tools)
})
args  <- commandArgs(trailingOnly = TRUE)
chunk <- as.integer(args[1]); K <- as.integer(args[2])
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
# STANDARDIZE on the validated legacy v3.2.0 binary (md5 5307) that reproduces HOBO
# validation — NOT in_files/Blois/musica. See memory musica-binary-legacy.
CFG_C3$musica_cmd <- "/home/corroyez/Documents/musica/musica"

n_per_cluster <- 100
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))  # ARTICLE-canonical sample (LAI_b≈6.73)
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
bl <- samp[, .(LAI_b=mean(LAI,na.rm=TRUE), Hmax_b=round(mean(floor(Hmax)+1,na.rm=TRUE)),
               FC_b=mean(fCover,na.rm=TRUE)), by=Cluster]   # fCover baseline = cluster mean (was 1)
if (chunk == 1) { cat("=== per-cluster baselines (LAD baseline = uniform) ===\n"); print(bl) }
set.seed(42)
sub <- samp[, .SD[sample(.N, min(.N, n_per_cluster))], by=Cluster]
sub[, pid := sprintf("S%04d", .I)]
my_rows <- which(((seq_len(nrow(sub)) - 1) %% K) == (chunk - 1))
cat(sprintf("chunk %d/%d : %d plots\n", chunk, K, length(my_rows)))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
ncdir <- normalizePath(file.path(CFG_C3$out_dir, "nc_shapley2x"), mustWork=FALSE)
dir.create(ncdir, showWarnings=FALSE, recursive=TRUE)
forcing <- CFG_C3$forcing_file   # relative path — normalizePath breaks MuSICA's forcing read
Fv <- c("LAI","Hmax","fCover","LAD"); n <- length(Fv)
coal <- as.matrix(expand.grid(LAI=0:1, Hmax=0:1, fCover=0:1, LAD=0:1))
w <- function(s) factorial(s) * factorial(n - s - 1) / factorial(n)

metric_one <- function(prow_df, bits, b) {
  owd <- getwd(); on.exit(setwd(owd), add=TRUE)
  lai <- if (bits["LAI"])    as.numeric(prow_df$LAI)    else b$LAI_b
  hmx <- if (bits["Hmax"])   as.numeric(prow_df$Hmax)   else b$Hmax_b
  fc  <- if (bits["fCover"]) as.numeric(prow_df$fCover) else b$FC_b
  ladf<- if (bits["LAD"]) make_lad_real else make_lad_uniform
  sc <- list(name="c", lai_fn=function(p) lai, hmax_fn=function(p) hmx,
             fcover_fn=function(p) fc, lad_fn=ladf, phenology_fn=NULL)
  nc <- file.path(ncdir, sprintf("%s_%s.nc", prow_df$pid, paste(bits, collapse="")))
  run_musica_one(prow_df, sc, nc, forcing, CFG_C3$musica_cmd); setwd(owd)
  if (!file.exists(nc) || file.size(nc) < 1000) return(NA_real_)
  dT <- tryCatch(extract_deltatmax_one(nc, dm, ds), error=function(e) NULL)
  if (is.null(dT)) NA_real_ else mean(dT$Delta_Tmax, na.rm=TRUE)
}
one_plot <- function(i) {
  prow <- as.data.frame(sub[i]); b <- as.list(bl[Cluster == prow$Cluster])
  vals <- setNames(numeric(nrow(coal)), apply(coal,1,paste,collapse=""))
  for (k in seq_len(nrow(coal))) { bits <- setNames(as.integer(coal[k,]), Fv)
    vals[paste(bits,collapse="")] <- tryCatch(metric_one(prow,bits,b), error=function(e) NA_real_) }
  phi <- setNames(numeric(n), Fv)
  for (v in Fv) { others <- setdiff(Fv,v); acc <- 0
    for (m in 0:length(others)) for (S in (if(m==0) list(character(0)) else combn(others,m,simplify=FALSE))) {
      b0 <- setNames(integer(n),Fv); b0[S] <- 1L; b1 <- b0; b1[v] <- 1L
      f1 <- vals[paste(b1[Fv],collapse="")]; f0 <- vals[paste(b0[Fv],collapse="")]
      if (is.finite(f1)&&is.finite(f0)) acc <- acc + w(m)*(f1-f0) }
    phi[v] <- acc }
  data.table(pid=prow$pid, Cluster=prow$Cluster,
    dTmax_full=vals[paste(rep(1,n),collapse="")], dTmax_base=vals[paste(rep(0,n),collapse="")],
    LAI=phi["LAI"], Hmax=phi["Hmax"], fCover=phi["fCover"], LAD=phi["LAD"])
}
pdir <- file.path(CFG_C3$out_dir, "tables", "shapley_parts"); dir.create(pdir, showWarnings=FALSE, recursive=TRUE)
partfile <- file.path(pdir, sprintf("part_%02d.csv", chunk))
if (file.exists(partfile)) file.remove(partfile)            # fresh: write incrementally per plot
nok <- 0L
for (j in seq_along(my_rows)) {
  r <- tryCatch(one_plot(my_rows[j]), error=function(e) NULL)
  if (!is.null(r)) { fwrite(r, partfile, append=file.exists(partfile)); nok <- nok + 1L }
  if (j %% 5 == 0) cat(sprintf("  chunk %d: %d/%d\n", chunk, j, length(my_rows)))
}
cat(sprintf("chunk %d DONE (%d/%d rows)\n", chunk, nok, length(my_rows)))
