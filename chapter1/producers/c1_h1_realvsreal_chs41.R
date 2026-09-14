# ==============================================================================
# c1_h1_realvsreal.R
#
# Loggerless test of the chapter's central negative claim ("H1 not realised": the
# realistic between-plot variation in the vertical LAD profile shape barely moves
# ΔTmax). The shipped test replaced each real profile by a UNIFORM one, an extreme
# upper bound; the empirical logger-based test (ΔR²=0) was cut with the HOBO
# migration. This replaces the real profile by a REALISTIC alternative — the
# sample-mean shape rescaled to the plot's own LAI and Hmax — so the contrast is the
# effect of a plot's shape DEVIATION from the mean, not of flattening it entirely.
#
# Three arms per plot, one convention: base (real), mean (sample-mean shape),
# unif (uniform). real-vs-mean is the realistic H1 test; real-vs-unif reproduces the
# published upper bound (0.15-0.51 degC) as a cross-check.
#
# Convention identical to the shipped Fig 4 harness (c1_sensitivity_perplot_units.R):
# FR-Blo_2021_v2.nc + per-plot wind correction, v3.2.3-iter, June-September,
# ΔTmax time-matched to the macro daily max, no clock shift, T at 1 m. So the
# numbers are directly comparable to the manuscript's real-vs-uniform. If Ch1
# migrates to CHS41-Rmerge/no-wind, BOTH contrasts must be recomputed together.
#
# Keyed clhs_%03d by row index (no seeded resample), 400 cLHS plots. Shell-parallel
# by chunk (no mclapply/fork). Resumable (skip existing nc).
#   Rscript c1_h1_realvsreal.R <chunk> <K>
# Out: out_files/Chapter1/tables/h1_realvsreal_chs41_part<chunk>of<K>.csv (combine after)
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)})
src <- list.files("R","\\.R$",full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)",src)]; invisible(lapply(src,source))
source("Chapter3_config.R")
CFG_C3$musica_cmd <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC <- "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"; ABL <- list("abl_flag"='"iter"')
WC <- FALSE; source("R/wind_correction.R"); .forc <- FORC
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day"); Z_FIX <- 1.0
MREF <- macro_ref(FORC, ds)

args <- commandArgs(trailingOnly=TRUE)
CHUNK <- if (length(args)>=1) as.integer(args[1]) else 1L
K     <- if (length(args)>=2) as.integer(args[2]) else 1L

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[, id_plot := sprintf("clhs_%03d", seq_len(.N))]
samp[, P := as.character(relabel_cluster(Cluster))]
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
mean_fn <- make_lad_mean_factory(as.data.frame(samp))   # sample-mean shape, rescaled per plot
mk_phen <- function(lai1) function(p){
  ph <- as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,
        leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=lai1))
  d <- ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]; d$Julian_day <- 366; rbind(ph,d) }
NCDIR <- "out_files/Chapter1/nc_h1_realvsreal_chs41"; dir.create(NCDIR, recursive=TRUE, showWarnings=FALSE)
metric_of <- function(nc){ m <- micro_hourly_at(nc, Z_FIX); if (is.null(m)) return(NA_real_); delta_tmax_mean(m, MREF, ds) }

run_get <- function(prow, tag, ladf){
  nc <- file.path(NCDIR, sprintf("%s_%s.nc", prow$id_plot, tag))
  sc <- list(lai_fn=function(p)prow$LAI, hmax_fn=function(p)prow$Hmax, fcover_fn=function(p)prow$fCover,
             lad_fn=ladf, phenology_fn=mk_phen(prow$LAI))
  if (!file.exists(nc) || file.size(nc) < 1000)
    tryCatch(run_musica_one(prow, sc, nc, .forc, CFG_C3$musica_cmd, extra_setup=ABL), error=function(e) NULL)
  metric_of(nc) }

my <- which(((seq_len(nrow(samp)) - 1L) %% K) == (CHUNK - 1L))
cat(sprintf("H1 real-vs-real | chunk %d/%d | %d plots | FR-Blo + windcorr | 3 arms/plot\n", CHUNK, K, length(my)))
R <- rbindlist(lapply(my, function(i){
  prow <- as.data.frame(samp[i])
  .forc <<- if (WC) windcorr_forcing(prow$Hmax) else FORC
  r <- tryCatch({
    real <- run_get(prow, "real", make_lad_real)
    mean <- run_get(prow, "mean", mean_fn)
    unif <- run_get(prow, "unif", make_lad_uniform)
    data.table(id_plot=prow$id_plot, P=prow$P, LAI=prow$LAI, Hmax=prow$Hmax, fCover=prow$fCover,
               dTmax_real=real, dTmax_mean=mean, dTmax_unif=unif,
               d_real_mean=real-mean, d_real_unif=real-unif)
  }, error=function(e){ cat(sprintf(" ERR %s: %s\n", prow$id_plot, conditionMessage(e))); NULL })
  if (!is.null(r)) cat(sprintf(" %s done (real-mean %+.3f, real-unif %+.3f)\n", prow$id_plot, r$d_real_mean, r$d_real_unif))
  r
}), fill=TRUE)

odir <- "out_files/Chapter1/tables"; dir.create(odir, recursive=TRUE, showWarnings=FALSE)
of <- if (K==1L) file.path(odir,"h1_realvsreal_chs41.csv") else file.path(odir, sprintf("h1_realvsreal_chs41_part%02dof%02d.csv", CHUNK, K))
fwrite(R, of); cat("wrote", of, "\n")
