# ==============================================================================
# Chapter 3 — WIND-CORRECTION dry-run. Same genuine-53 recipe as
# gen_genuine53_v323.R, but the forcing wind is scaled per plot by the log-profile
# factor wind_factor(Hmax) (R/wind_correction.R). Base forcing is the SAME one the
# current Ch3 runs used (musica_in_Blois_pblh.nc), so ONLY the wind changes -> the
# effect is isolated. 3 reference scenarios x 53 plots -> nc_genuine53_windcorr.
#   Rscript gen_windcorr_dry.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config_CHS41.R"); source("R/wind_correction.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC_BASE <- CFG_C3$forcing_file; stopifnot(file.exists(FORC_BASE))
ABL  <- list("abl_flag"='"iter"'); LY <- 2020:2022; NCORES <- 10L
SCEN <- NULL   # set to all scenarios after all_sc is built (see below)
NCROOT <- file.path(CFG_C3$out_dir, "nc_windcorr")
WFDIR  <- file.path(CFG_C3$out_dir, "windcorr_forc_chs41")   # cache of wind-scaled pblh forcings
dir.create(WFDIR, recursive=TRUE, showWarnings=FALSE)

# per-plot wind-scaled copy of the pblh base (cached by rounded factor)
windcorr_pblh <- function(hmax){
  f  <- wind_factor(pmax(hmax, 2))
  ff <- file.path(WFDIR, sprintf("pblh_f%.2f.nc", round(f, 2)))
  if(!file.exists(ff) || file.size(ff) < 1e5){
    file.copy(FORC_BASE, ff, overwrite=TRUE)
    nc <- nc_open(ff, write=TRUE)
    for(v in c("Wind_E","Wind_N")) ncvar_put(nc, v, ncvar_get(nc, v) * f)
    nc_close(nc)
  }
  ff
}

inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("missing 2020/365"); d$Julian_day<-366L; rbind(ph,d) }

# ---- real-53 df + scenario defs (identical to genuine53 generator) ------------
tb <- load_lai_prep(CFG_C3)$ts_by_plot; dopt <- CFG_C3$d_opt_m
df <- readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")); setDT(df)
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
all_sc <- make_all_scenarios_c3(tb, LY, dopt, "full")
# Restricted to the four chapter scenarios that this builder produces. The fifth,
# STATIC_S2_OPT, has its own generator (gen_opt_windcorr_CHS41.R).
SCEN <- intersect(c("STATIC_ALS","STATIC_S2_ATBD","STATIC_S2_RESCALED","DYN_S2_ANNUAL"),
                  names(all_sc))
cat(sprintf("full run: %d scenarios: %s\n", length(SCEN), paste(SCEN, collapse=", ")))
stopifnot(all(SCEN %in% names(all_sc)))

# ---- un-masked ATBD daily series for ALL 53 (for the DYN scenario) ------------
NM <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(NM, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs=crs(rast(files[1])))
vals <- as.data.frame(terra::extract(rast(files), pts))[,-1,drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
ts_atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3); stopifnot(length(ts_atbd)==53)
atbd_summer <- sapply(df$pid, function(p){ a<-ts_atbd[[p]]; mean(a$lai[a$doy>=152 & a$doy<=244], na.rm=TRUE) })

wrap_genuine <- function(sc){
  is_dyn <- !is.null(sc$phenology_fn); lai_fn <- sc$lai_fn
  sc$phenology_fn <- function(pr){
    if(!is_dyn){
      ph <- calc_phenology(list.year=LY, nleafage=1, budburst_date=115, leaf_age_max_in=0.56,
              relative_age_firstmax=0.10, relative_age_lastmax=0.75, LAI_max_per_cohort=as.numeric(lai_fn(pr)))
    } else {
      pid <- sprintf("X%d_Y%d", round(as.numeric(pr$x)), round(as.numeric(pr$y)))
      a <- ts_atbd[[pid]]; sc_scale <- as.numeric(lai_fn(pr))/atbd_summer[[pid]]
      ph <- make_phenology_from_s2(data.frame(doy=a$doy, lai=a$lai*sc_scale), LY)
    }
    inject_seed(ph)
  }
  sc
}
SC <- lapply(all_sc[SCEN], wrap_genuine)

for(s in names(SC)) dir.create(file.path(NCROOT,s),recursive=TRUE,showWarnings=FALSE)
# pre-build the wind-factor cache serially (avoid races in mclapply)
invisible(lapply(unique(df$Hmax), windcorr_pblh))
tasks <- rbindlist(lapply(names(SC), function(s) data.table(scn=s, i=seq_len(nrow(df)))))
run_one <- function(k){ s<-tasks$scn[k]; i<-tasks$i[k]; prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,s,sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s/%s",s,prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  forc <- windcorr_pblh(prow$Hmax)
  ok<-tryCatch({ run_musica_one(prow,SC[[s]],nc,forc,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s/%s", if(ok)"OK" else "FAIL", s, prow$id_plot) }
cat(sprintf("WINDCORR full: %d runs (%d scen x 53) on %d cores\n", nrow(tasks), length(SC), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(seq_len(nrow(tasks)),run_one,mc.cores=NCORES,mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("\nDONE in %.1f min | %s\n", dt, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILURES:\n");cat(fails,sep="\n");cat("\n")}
