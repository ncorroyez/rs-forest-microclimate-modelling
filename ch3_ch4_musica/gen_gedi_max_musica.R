# ==============================================================================
# gen_gedi_fall_musica.R — Chapter 4: scenario DYN_ALS_GEDIMAX.
# Clone of DYN_S2_ANNUAL (LiDAR magnitude x S2 timing, real LAD, per-plot annual
# series) with ONE change: the autumn limb (DOY >= 250) is replaced by the
# GEDI-constrained leaf-fraction trajectory measured at Blois
# (c4_gedi_blois_summer.R, structure-adjusted power beams):
#   leaf fraction 1.00 @ DOY 250 (plateau holds), 0.951 @ 282 (Oct 9),
#   0.321 @ 318 (Nov 14), 0.08 @ 350; anchored on each plot's own summer
#   plateau (mean of its annual series, DOY 180-240). Leaf-out limb untouched
#   (GEDI validated the S2 timing there). Isolates the GEDI autumn correction
#   against DYN_S2_ANNUAL. genuine-53 windcorr machinery (v3.2.3-iter, 2xLAI).
#   DRY=1 Rscript gen_gedi_fall_musica.R ; Rscript gen_gedi_fall_musica.R
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

# ---- GEDI autumn leaf-fraction curve (Blois, Table14 adjusted) ---------------
f_gedi <- approxfun(c(250, 282, 318, 350), c(1.00, 0.951, 0.321, 0.08), rule = 2)

prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
ts_gedifall <- setNames(lapply(df$pid, function(p) {
  a <- tb$annual[[p]]
  if (is.null(a) || !nrow(a)) return(NULL)
  mag <- mean(a$lai[a$doy %in% 180:240], na.rm = TRUE)
  lai <- a$lai
  core <- a$doy >= 180 & a$doy < 250   # fill dips only, keep the peak
  lai[core] <- pmax(lai[core], mag)
  fall <- a$doy >= 250 & a$doy <= 350
  lai[fall] <- mag * f_gedi(a$doy[fall])
  data.frame(doy = a$doy, lai = lai)
}), df$pid)
# plots without an annual series (6 bare plots) fall back to the parametric
# phenology inside wrap_iter, as in the rest of the genuine-53 pipeline
cat("plots with GEDI-fall series:", sum(!vapply(ts_gedifall, is.null, TRUE)),
    "/", nrow(df), "(others -> parametric fallback)\n")
ts_gedifall <- ts_gedifall[!vapply(ts_gedifall, is.null, TRUE)]
ex <- ts_gedifall[["X366240_Y5269440"]]
if (!is.null(ex)) cat("example plot 41_01 lai @ doy 240/282/318:",
                      round(ex$lai[ex$doy %in% c(240, 282, 318)], 2), "\n")

.fn <- function(col) function(pr) as.numeric(pr[[col]])
mk  <- function(ts) make_phenology_fn_factory(ts, LY)
SC  <- wrap_iter(list(name="DYN_ALS_GEDIMAX",
  lai_fn=.fn("LAI_ALS"), hmax_fn=.fn("Hmax"), fcover_fn=.fn("fCover"),
  lad_fn=make_lad_real, phenology_fn=mk(ts_gedifall)))
dir.create(file.path(NCROOT, "DYN_ALS_GEDIMAX"), recursive=TRUE, showWarnings=FALSE)

invisible(lapply(unique(df$Hmax), windcorr_pblh))
rows <- if (DRY) seq_len(6) else seq_len(nrow(df))
run_one <- function(i){
  prow <- as.data.frame(df[i, ])
  nc <- file.path(NCROOT, "DYN_ALS_GEDIMAX",
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
