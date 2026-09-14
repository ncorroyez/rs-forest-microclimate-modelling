# ==============================================================================
# Chapter 3 — per-point cLHS Shapley decomposition of MuSICA buffering (ΔTmax),
# grouped by FPCA cluster, PER-CLUSTER baseline, at the corrected 2 × LAI.
# Each plot's ΔTmax is decomposed (vs its CLUSTER mean canopy) into factor
# contributions {LAI, Hmax, fCover, LAD}. Better than mean-profile archetypes:
# MuSICA is nonlinear (f(mean) biased), so per-point + aggregate is unbiased and
# yields the within-cluster distribution. 100 plots/cluster × 2^4 coalitions.
# Parallel (mclapply); setup_dir() gives each MuSICA call a unique tempfile dir.
# Run from z_Example root:  Rscript c3_shapley_clhs_2x.R
# ==============================================================================
# !!! DEPRECATED — DO NOT RUN. This script uses the WRONG cLHS sample
# (clhs_sample_floor05.rds, LAI_b~3.1) and wrote the SAME external figure filename
# as the canonical floor05_v2 run, causing the Fig 2/3 provenance ambiguity flagged
# in review/pipeline_critique.md (item 2.1). The canonical, reproducible attribution
# is now pipeline/11_clhs_attribution.R (reads floor05_v2 parts; writes local
# out_files/Chapter3/figures/). Kept only for history.
stop("DEPRECATED: use pipeline/11_clhs_attribution.R (canonical floor05_v2). See pipeline_critique.md 2.1.")

suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate)
  library(dplyr); library(tidyr); library(data.table); library(ggplot2)
  library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

n_per_cluster <- 100; MC <- 10
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
# PER-CLUSTER baselines for the numeric factors; LAD baseline = uniform.
bl <- samp[, .(LAI_b=mean(LAI,na.rm=TRUE), Hmax_b=round(mean(floor(Hmax)+1,na.rm=TRUE)),
               FC_b=1), by=Cluster]
cat("=== per-cluster baselines ===\n"); print(bl)
set.seed(42)
sub <- samp[, .SD[sample(.N, min(.N, n_per_cluster))], by=Cluster]
sub[, pid := sprintf("S%04d", .I)]
cat(sprintf("subset: %d plots (%s/cluster)\n", nrow(sub), paste(sub[,.N,by=Cluster]$N,collapse="/")))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- extract_macro_daily(CFG_C3$forcing_file, ds)
ncdir <- normalizePath(file.path(CFG_C3$out_dir, "nc_shapley2x"), mustWork=FALSE)
dir.create(ncdir, showWarnings=FALSE, recursive=TRUE)
forcing <- normalizePath(CFG_C3$forcing_file)

Fv <- c("LAI","Hmax","fCover","LAD"); n <- length(Fv)
coal <- as.matrix(expand.grid(LAI=0:1, Hmax=0:1, fCover=0:1, LAD=0:1))
w <- function(s) factorial(s) * factorial(n - s - 1) / factorial(n)

metric_one <- function(prow_df, bits, b) {
  owd <- getwd(); on.exit(setwd(owd), add=TRUE)   # callmusica setwd()s; restore cwd
  lai <- if (bits["LAI"])    as.numeric(prow_df$LAI)    else b$LAI_b
  hmx <- if (bits["Hmax"])   as.numeric(prow_df$Hmax)   else b$Hmax_b
  fc  <- if (bits["fCover"]) as.numeric(prow_df$fCover) else b$FC_b
  ladf<- if (bits["LAD"]) make_lad_real else make_lad_uniform
  sc <- list(name="c", lai_fn=function(p) lai, hmax_fn=function(p) hmx,
             fcover_fn=function(p) fc, lad_fn=ladf, phenology_fn=NULL)
  nc <- file.path(ncdir, sprintf("%s_%s.nc", prow_df$pid, paste(bits, collapse="")))
  run_musica_one(prow_df, sc, nc, forcing, CFG_C3$musica_cmd)
  setwd(owd)
  if (!file.exists(nc) || file.size(nc) < 1000) return(NA_real_)   # guard 32-byte stubs
  dT <- tryCatch(extract_deltatmax_one(nc, dm, ds), error=function(e) NULL)
  if (is.null(dT)) NA_real_ else mean(dT$Delta_Tmax, na.rm=TRUE)
}

one_plot <- function(i) {
  prow <- as.data.frame(sub[i])                       # data.frame for make_lad_real
  b <- as.list(bl[Cluster == prow$Cluster])
  vals <- setNames(numeric(nrow(coal)), apply(coal, 1, paste, collapse=""))
  for (k in seq_len(nrow(coal))) {
    bits <- setNames(as.integer(coal[k, ]), Fv)
    vals[paste(bits, collapse="")] <- tryCatch(metric_one(prow, bits, b), error=function(e) NA_real_)
  }
  phi <- setNames(numeric(n), Fv)
  for (v in Fv) {
    others <- setdiff(Fv, v); acc <- 0
    for (m in 0:length(others)) for (S in (if (m==0) list(character(0)) else combn(others, m, simplify=FALSE))) {
      b0 <- setNames(integer(n), Fv); b0[S] <- 1L; b1 <- b0; b1[v] <- 1L
      f1 <- vals[paste(b1[Fv], collapse="")]; f0 <- vals[paste(b0[Fv], collapse="")]
      if (is.finite(f1) && is.finite(f0)) acc <- acc + w(m) * (f1 - f0)
    }
    phi[v] <- acc
  }
  data.table(pid=prow$pid, Cluster=prow$Cluster,
    dTmax_full=vals[paste(rep(1,n),collapse="")], dTmax_base=vals[paste(rep(0,n),collapse="")],
    LAI=phi["LAI"], Hmax=phi["Hmax"], fCover=phi["fCover"], LAD=phi["LAD"])
}

cat(sprintf("running %d plots × 16 coalitions on %d cores...\n", nrow(sub), MC))
res <- mclapply(seq_len(nrow(sub)), function(i) tryCatch(one_plot(i), error=function(e) NULL),
                mc.cores=MC, mc.preschedule=FALSE)
S <- rbindlist(Filter(Negate(is.null), res))
fwrite(S, "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables/Table6c_shapley_clhs_2x.csv")

cat("\n=== mean Shapley (°C on ΔTmax, vs cluster baseline) per cluster ===\n")
agg <- S[, .(LAI=mean(LAI,na.rm=T), Hmax=mean(Hmax,na.rm=T), fCover=mean(fCover,na.rm=T),
             LAD=mean(LAD,na.rm=T), n=.N), by=Cluster][order(Cluster)]
print(agg, digits=3)
cat("\n=== within-cluster IMPORTANCE = mean |Shapley| per factor ===\n")
imp <- S[, .(LAI=mean(abs(LAI),na.rm=T), Hmax=mean(abs(Hmax),na.rm=T),
             fCover=mean(abs(fCover),na.rm=T), LAD=mean(abs(LAD),na.rm=T)), by=Cluster][order(Cluster)]
print(imp, digits=3)
cat("\n(additivity check: sum(Shapley) vs dTmax_full - dTmax_base)\n")
print(S[, .(sum_shap=mean(LAI+Hmax+fCover+LAD,na.rm=T), fmb=mean(dTmax_full-dTmax_base,na.rm=T)), by=Cluster], digits=3)

ml <- melt(S, id.vars=c("pid","Cluster"), measure.vars=Fv, variable.name="factor", value.name="shapley")
ml[, Cluster := paste0("Cluster ", Cluster)]
g <- ggplot(ml, aes(factor, shapley, fill=factor)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey60") +
  geom_boxplot(outlier.size=0.4, alpha=0.85) + facet_wrap(~Cluster, nrow=1) +
  scale_fill_brewer(palette="Set2", guide="none") +
  labs(title="Per-point cLHS Shapley: drivers of within-cluster ΔTmax variation (2×LAI)",
       subtitle="100 cLHS plots/cluster; PER-CLUSTER baseline (factor at cluster mean); LAD baseline = uniform. Each plot decomposed vs its cluster archetype.",
       x=NULL, y=expression(Shapley~(degree*C~"on"~Delta*T[max]))) +
  theme_minimal(base_size=11) + theme(plot.subtitle=element_text(size=8))
ggsave("/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures/FigSh_shapley_clhs_clusters.png", g, width=11, height=4, dpi=150)
cat("\nDONE\n")
