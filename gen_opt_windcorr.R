# ==============================================================================
# Chapter 3 — CH2 ACQUIS AS A CH3 FORCING: the forest-tuned Sentinel-2 "opt"
# retrieval (atbd_optim) run through MuSICA as a static magnitude on the full
# LiDAR LAD. This Ch2 product was so far only used in LAI space (Fig A) and as an
# RF feature; never as a forcing. Per-plot magnitude = summer (DOY 152-244) mean
# of the smoothed opt series, same recipe as df$S2opt in c3_figset_data_wc.R.
#
# PREDICTION (before run): in LAI space "opt" INVERTS the between-plot ranking
# (r=+0.19 vs -0.17 for ATBD). So through MuSICA STATIC_S2_OPT should rank WORSE
# than STATIC_S2_ATBD, degrading especially the open regime where ATBD works
# (0.71). This confirms at the microclimate level the two-sided negative control
# of §4.1: a more forest-accurate optical LAI is a worse microclimate driver.
#
# Same windcorr harness as gen_blend_windcorr.R.
#   DRY=1 Rscript gen_opt_windcorr.R   -> magnitude check + 3 plots -> _dry
#         Rscript gen_opt_windcorr.R   -> 53 plots -> nc_genuine53_windcorr/STATIC_S2_OPT
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
dir.create(file.path(NCROOT,"STATIC_S2_OPT"), recursive=TRUE, showWarnings=FALSE)
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

# ---- per-plot forest-tuned "opt" magnitude (summer mean of smoothed series) ---
df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
NM <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"     # READ-ONLY
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs="EPSG:32631")
opt_files <- list.files(NM, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$", full.names=TRUE)
opt_dates <- as.Date(str_extract(basename(opt_files), "\\d{4}-\\d{2}-\\d{2}"))
vals <- as.data.frame(terra::extract(rast(opt_files), pts))[,-1,drop=FALSE]
long <- rbindlist(lapply(seq_along(opt_dates), function(i)
  data.table(plot_id=df$pid, doy=as.integer(format(opt_dates[i],"%j")), lai=pmax(vals[[i]],0))))
ts <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
df$LAI_S2_OPT <- sapply(df$pid, function(p){ a<-ts[[p]]; mean(a$lai[a$doy>=152 & a$doy<=244], na.rm=TRUE) })
if(any(is.na(df$LAI_S2_OPT))) df$LAI_S2_OPT[is.na(df$LAI_S2_OPT)] <- median(df$LAI_S2_OPT, na.rm=TRUE)

.fn <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
SC <- wrap_iter(list(name="STATIC_S2_OPT", lai_fn=.fn("LAI_S2_OPT"), hmax_fn=.fn("Hmax"),
                     fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=NULL))

if(DRY){ o<-order(df$LAI_ALS); idx<-unique(c(o[1],o[ceiling(nrow(df)/2)],o[nrow(df)])) } else idx<-seq_len(nrow(df))
cat(sprintf("LAI-space rank cor(opt, obs?) not needed here. cor(LAI_S2_OPT, LAI_S2_ATBD)=%.2f  cor(opt,LAI_ALS)=%.2f\n",
            cor(df$LAI_S2_OPT,df$LAI_S2_ATBD,use="complete.obs"), cor(df$LAI_S2_OPT,df$LAI_ALS,use="complete.obs")))
cat("=== opt magnitude check (selected plots) ===\n")
print(df[idx,.(id_plot, LAI_ALS=round(LAI_ALS,2), LAI_S2_ATBD=round(LAI_S2_ATBD,2), LAI_S2_OPT=round(LAI_S2_OPT,2))])

invisible(lapply(unique(df$Hmax[idx]), windcorr_pblh))
run_one <- function(i){ prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,"STATIC_S2_OPT",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s",prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({ run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s", if(ok)"OK" else "FAIL", prow$id_plot) }
cat(sprintf("\nSTATIC_S2_OPT windcorr%s: %d plots on %d cores\n", if(DRY)" [DRY]" else "", length(idx), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(idx,run_one,mc.cores=min(NCORES,length(idx)),mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("DONE %.2f min (~%.0f min est for 53) | %s\n", dt, dt/length(idx)*53, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILS:\n");cat(fails,sep="\n");cat("\n")}
