# ==============================================================================
# Chapter 3 — CONTINUOUS height-weighted blend fusion (soft twin of the height-
# gated switch). Per plot the leaf-area MAGNITUDE is a logistic blend of the two
# sensors on the OPERATIONAL variable FORMS-H (no per-pixel LiDAR needed for the
# weight), centred at the paper's 17.8 m gate with a FIXED scale s = 3 m:
#     w   = 1 / (1 + exp(-(FORMS_H - 17.8)/3))     (tall -> LiDAR, short -> S2)
#     LAI = w*LAI_ALS + (1-w)*LAI_S2_ATBD          (full LiDAR LAD, static pheno)
# Same windcorr recipe / harness as gen_fusionh_windcorr.R (wrap_iter, per-plot
# wind scaling, ABL iter). s is fixed, NOT swept (a swept s would be a
# hyperparameter fit on 53 plots and worthless).
#   DRY=1 Rscript gen_blend_windcorr.R  -> 3 plots (dense/open/transition), _dry
#         Rscript gen_blend_windcorr.R  -> 53 plots -> nc_genuine53_windcorr/BLEND_H
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
dir.create(file.path(NCROOT,"BLEND_H"), recursive=TRUE, showWarnings=FALSE)
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

# ---- blend magnitude (FIXED params) -----------------------------------------
H0 <- 17.8; S0 <- 3
wfun <- function(fh) 1/(1+exp(-(fh-H0)/S0))
df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
stopifnot(all(c("LAI_ALS","LAI_S2_ATBD","FORMS_H","Hmax","fCover") %in% names(df)))
if(any(is.na(df$LAI_S2_ATBD))) df$LAI_S2_ATBD[is.na(df$LAI_S2_ATBD)] <- median(df$LAI_S2_ATBD, na.rm=TRUE)
df$w <- wfun(df$FORMS_H)
df$LAI_BLEND <- df$w*df$LAI_ALS + (1-df$w)*df$LAI_S2_ATBD

.fn <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
SC <- wrap_iter(list(name="BLEND_H", lai_fn=.fn("LAI_BLEND"), hmax_fn=.fn("Hmax"),
                     fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=NULL))

if(DRY){
  idx <- unique(c(which.min(abs(df$FORMS_H-30)), which.min(abs(df$FORMS_H-8)), which.min(abs(df$FORMS_H-17.8))))
} else idx <- seq_len(nrow(df))
cat("=== blend magnitude check (selected plots) ===\n")
print(df[idx,.(id_plot, FORMS_H=round(FORMS_H,1), w=round(w,2),
               LAI_ALS=round(LAI_ALS,2), LAI_S2_ATBD=round(LAI_S2_ATBD,2), LAI_BLEND=round(LAI_BLEND,2))])

invisible(lapply(unique(df$Hmax[idx]), windcorr_pblh))
run_one <- function(i){ prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,"BLEND_H",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s",prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({ run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s", if(ok)"OK" else "FAIL", prow$id_plot) }
cat(sprintf("\nBLEND_H windcorr%s: %d plots on %d cores\n", if(DRY)" [DRY]" else "", length(idx), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(idx,run_one,mc.cores=min(NCORES,length(idx)),mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("DONE %.2f min (%.2f min/plot wall, ~%.0f min est. for 53) | %s\n",
    dt, dt/length(idx), dt/length(idx)*53, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILS:\n");cat(fails,sep="\n");cat("\n")}
cat(paste(res,collapse="\n"),"\n")
