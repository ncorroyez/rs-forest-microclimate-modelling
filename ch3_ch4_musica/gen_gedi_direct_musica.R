# ==============================================================================
# gen_gedi_direct_musica.R — Chapter 4: MuSICA scenario STATIC_GEDI_DIRECT.
# Magnitude = mean PAI of the 3 nearest power-beam GEDI footprints (GEDI_k3,
# built by c4_gedi_direct.R), used DIRECTLY (no regression). NAIVE recipe
# (FORMS-H height + uniform LAD + parametric static phenology), genuine-53
# machinery (v3.2.3-iter, per-plot wind correction, 2xLAI in run_musica_one).
# All 53 plots are simulated; scoring stratifies by distance to the footprint
# (cmp_gedi_direct.R).
#   DRY=1 Rscript gen_gedi_direct_musica.R   # 6 plots
#   Rscript gen_gedi_direct_musica.R         # full 53
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC_BASE <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC_BASE))
ABL  <- list("abl_flag" = '"iter"'); LY <- 2020:2022; NCORES <- 10L
DRY  <- nzchar(Sys.getenv("DRY"))
NCROOT <- file.path(CFG_C3$out_dir, "nc_genuine53_windcorr")
WFDIR  <- file.path(CFG_C3$out_dir, "windcorr_forc_pblh")

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

df <- readRDS(file.path(CFG_C3$out_dir, "lai_prep", "df_plots_gedi_direct.rds"))
df <- as.data.frame(df)
stopifnot(all(c("GEDI_k3", "dist_fp", "FORMS_H") %in% names(df)))
cat("GEDI_k3 fed: mean", round(mean(df$GEDI_k3), 2),
    "range", paste(round(range(df$GEDI_k3), 2), collapse="-"), "\n")

.fn <- function(col) function(pr) as.numeric(pr[[col]])
SC <- wrap_iter(list(name="STATIC_GEDI_DIRECT",
  lai_fn=.fn("GEDI_k3"), hmax_fn=.fn("FORMS_H"), fcover_fn=.fn("fCover"),
  lad_fn=make_lad_uniform, phenology_fn=NULL))
dir.create(file.path(NCROOT, "STATIC_GEDI_DIRECT"), recursive=TRUE, showWarnings=FALSE)

invisible(lapply(unique(df$Hmax), windcorr_pblh))
rows <- if (DRY) seq_len(6) else seq_len(nrow(df))
run_one <- function(i){
  prow <- df[i, ]
  nc <- file.path(NCROOT, "STATIC_GEDI_DIRECT",
                  sprintf("musica_out_HOBO_%s.nc", prow$id_plot))
  if (file.exists(nc) && file.size(nc) > 1000) return(sprintf("skip %s", prow$id_plot))
  if (file.exists(nc)) file.remove(nc)
  forc <- windcorr_pblh(prow$Hmax)
  ok <- tryCatch({ run_musica_one(prow, SC, nc, forc, BIN, extra_setup=ABL)
                   file.exists(nc) && file.size(nc) > 1000 }, error=function(e) FALSE)
  sprintf("%s %s", if (ok) "OK" else "FAIL", prow$id_plot)
}
cat(sprintf("%s: %d sims on %d cores\n", if (DRY) "DRY" else "FULL", length(rows), NCORES))
t0 <- Sys.time(); res <- unlist(mclapply(rows, run_one, mc.cores=NCORES, mc.preschedule=FALSE))
dt <- as.numeric(difftime(Sys.time(), t0, units="mins")); tab <- table(sub(" .*","",res))
cat(sprintf("DONE %.1f min | %s\n", dt, paste(names(tab), tab, sep="=", collapse=" ")))
fails <- grep("^FAIL", res, value=TRUE)
if (length(fails)) { cat("FAILS:\n"); cat(fails, sep="\n") }
