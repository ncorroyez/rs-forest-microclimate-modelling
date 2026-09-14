# ==============================================================================
# FUSION_GEDIMAX (layered §3.5) re-run WITH per-plot wind correction, so the fusion bar
# in figB/figI is consistent with the wind-corrected genuine-53 set. Same recipe
# as run_v323iter_ch3_extra.R (masked basis, wrap_iter), only the forcing wind is
# scaled per plot. Writes out_files/Chapter3/nc_genuine53_windcorr/FUSION_GEDIMAX/.
#   Rscript gen_fusionh_windcorr.R
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
NCROOT <- file.path(CFG_C3$out_dir, "nc_genuine53_windcorr"); dir.create(file.path(NCROOT,"FUSION_GEDIMAX"), recursive=TRUE, showWarnings=FALSE)
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

# ---- FUSION_GEDIMAX recipe (identical to run_v323iter_ch3_extra.R, masked basis) -----
prep <- load_lai_prep(CFG_C3); df <- prep$df_plots; tb <- prep$ts_by_plot; dopt <- CFG_C3$d_opt_m
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
df$LAI_ALS_DOPT[is.na(df$LAI_ALS_DOPT)] <- df$LAI_ALS[is.na(df$LAI_ALS_DOPT)]
ropt <- mean(df$LAI_S2_DOPT, na.rm=TRUE)/mean(df$LAI_S2_ATBD, na.rm=TRUE)
df$LAI_S2_DOPT[is.na(df$LAI_S2_DOPT)] <- df$LAI_S2_ATBD[is.na(df$LAI_S2_DOPT)]*ropt
.fn <- function(col) function(pr) as.numeric(pr[[col]]); lad <- make_lad_real
mk  <- function(ts) make_phenology_fn_factory(ts, LY)
get_atbd <- function(p) tb$atbd[[p]]
ropt_p   <- setNames(df$LAI_S2_DOPT/df$LAI_S2_ATBD, df$pid)
safe_ser <- function(p, scale, fb){ a<-get_atbd(p)
  if(is.null(a)||!all(c("doy","lai")%in%names(a))||nrow(a)==0) return(data.frame(doy=1:365, lai=rep(fb,365)))
  data.frame(doy=a$doy, lai=a$lai*scale) }
flat_ts  <- function(p, val){ a<-get_atbd(p); doy<-if(is.null(a)||!nrow(a))1:365 else a$doy; data.frame(doy=doy, lai=rep(val,length(doy))) }
ts_opt  <- setNames(lapply(seq_len(nrow(df)), function(i) safe_ser(df$pid[i], ropt_p[[df$pid[i]]], df$LAI_S2_DOPT[i])), df$pid)
f_gedi_fall <- approxfun(c(250,282,318,350), c(1.00,0.951,0.321,0.08), rule=2)
gedimax <- function(a){ mag<-mean(a$lai[a$doy %in% 180:240], na.rm=TRUE)
  core<-a$doy>=180 & a$doy<250; a$lai[core]<-pmax(a$lai[core], mag)
  fall<-a$doy>=250 & a$doy<=350; a$lai[fall]<-mag*f_gedi_fall(a$doy[fall]); a }
dyn_als <- function(p, val){ a<-tb$annual[[p]]; if(is.null(a)||all(is.na(a$lai))) flat_ts(p,val) else gedimax(a) }
ts_fus  <- setNames(lapply(seq_len(nrow(df)), function(i){ p<-df$pid[i]
  if(df$Hmax[i] < dopt) ts_opt[[p]] else dyn_als(p, df$LAI_ALS[i]) }), df$pid)
SC <- wrap_iter(list(name="FUSION_GEDIMAX", lai_fn=.fn("LAI_ALS"), hmax_fn=.fn("Hmax"),
                     fcover_fn=.fn("fCover"), lad_fn=lad, phenology_fn=mk(ts_fus)))

invisible(lapply(unique(df$Hmax), windcorr_pblh))
run_one <- function(i){ prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,"FUSION_GEDIMAX",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s",prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({ run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s", if(ok)"OK" else "FAIL", prow$id_plot) }
cat(sprintf("FUSION_GEDIMAX windcorr: %d plots on %d cores\n", nrow(df), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(seq_len(nrow(df)),run_one,mc.cores=NCORES,mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("DONE %.1f min | %s\n", dt, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILS:\n");cat(fails,sep="\n");cat("\n")}
