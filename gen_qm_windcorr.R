# ==============================================================================
# Chapter 3 — QUANTILE-MATCHED Sentinel-2 forcing (distribution bias correction).
# Map each plot's S2 LAI onto the LiDAR LAI at the same empirical quantile, so the
# corrected product has the LiDAR marginal DISTRIBUTION while keeping S2's between-
# plot RANK. Literature archetype F (CDF matching); grep of the 886-entry library
# and web searches found no LAI-fusion exemplar -> under-explored here.
#
# PREDICTION (stated before the run, per design): quantile mapping is a MONOTONE
# transform of LAI_S2, so it preserves S2's rank exactly (Spearman(S2,ALS)=-0.16).
# It therefore CANNOT fix the dense-canopy ranking failure; it can only correct the
# marginal distribution, i.e. warm bias and RMSE. Expected: pooled/dense ranking
# R2 ~ like raw S2 (dense ~0.02), with reduced bias. If dense R2 > ~0.05, suspect a
# pipeline bug, not a success.
#
# Same windcorr harness as gen_blend_windcorr.R (wrap_iter, per-plot wind, ABL iter).
#   DRY=1 Rscript gen_qm_windcorr.R  -> magnitude check + 3 plots -> _dry
#         Rscript gen_qm_windcorr.R  -> 53 plots -> nc_genuine53_windcorr/QUANTILE_S2
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC_BASE <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC_BASE))
ABL  <- list("abl_flag" = '"iter"'); LY <- 2020:2022; NCORES <- 10L
DRY  <- nzchar(Sys.getenv("DRY"))
NCROOT <- if(DRY) file.path(CFG_C3$out_dir,"nc_genuine53_windcorr_dry") else file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
dir.create(file.path(NCROOT,"QUANTILE_S2"), recursive=TRUE, showWarnings=FALSE)
WFDIR  <- file.path(CFG_C3$out_dir, "windcorr_forc_pblh"); dir.create(WFDIR, recursive=TRUE, showWarnings=FALSE)

windcorr_pblh <- function(hmax){
  f <- wind_factor(pmax(hmax, 2)); ff <- file.path(WFDIR, sprintf("pblh_f%.2f.nc", round(f, 2)))
  if(!file.exists(ff) || file.size(ff) < 1e5){ file.copy(FORC_BASE, ff, overwrite=TRUE)
    nc <- nc_open(ff, write=TRUE); for(v in c("Wind_E","Wind_N")) ncvar_put(nc, v, ncvar_get(nc, v) * f); nc_close(nc) }
  ff
}
inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("missing 2020/365"); d$Julian_day<-366L; rbind(ph,d) }
wrap_iter <- function(sc){ orig<-sc$phenology_fn; lai_fn<-sc$lai_fn
  sc$phenology_fn <- function(pr){ ph<-if(!is.null(orig))orig(pr) else NULL
    if(is.null(ph)) ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,
        relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)))
    inject_seed(ph) }; sc }

# ---- quantile-matched magnitude ---------------------------------------------
df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
if(any(is.na(df$LAI_S2_ATBD))) df$LAI_S2_ATBD[is.na(df$LAI_S2_ATBD)] <- median(df$LAI_S2_ATBD, na.rm=TRUE)
n <- nrow(df)
p <- rank(df$LAI_S2_ATBD, ties.method="average")/(n+1)          # S2 plotting position
df$LAI_QM <- as.numeric(quantile(df$LAI_ALS, probs=p, type=7))   # LiDAR value at that quantile

.fn <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
SC <- wrap_iter(list(name="QUANTILE_S2", lai_fn=.fn("LAI_QM"), hmax_fn=.fn("Hmax"),
                     fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=NULL))

if(DRY){ o<-order(df$LAI_S2_ATBD); idx<-unique(c(o[1],o[ceiling(n/2)],o[n])) } else idx<-seq_len(nrow(df))
cat(sprintf("Spearman rank cor(LAI_QM, LAI_S2_ATBD) = %.3f (must be 1.0: QM is monotone)\n",
            cor(df$LAI_QM, df$LAI_S2_ATBD, method="spearman")))
cat("=== quantile-matched magnitude check (selected plots) ===\n")
print(df[idx,.(id_plot, LAI_ALS=round(LAI_ALS,2), LAI_S2_ATBD=round(LAI_S2_ATBD,2),
               p=round(p[idx],2), LAI_QM=round(LAI_QM,2))])

invisible(lapply(unique(df$Hmax[idx]), windcorr_pblh))
run_one <- function(i){ prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,"QUANTILE_S2",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s",prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({ run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s", if(ok)"OK" else "FAIL", prow$id_plot) }
cat(sprintf("\nQUANTILE_S2 windcorr%s: %d plots on %d cores\n", if(DRY)" [DRY]" else "", length(idx), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(idx,run_one,mc.cores=min(NCORES,length(idx)),mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("DONE %.2f min (~%.0f min est for 53) | %s\n", dt, dt/length(idx)*53, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILS:\n");cat(fails,sep="\n");cat("\n")}
