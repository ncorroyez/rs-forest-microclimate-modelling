# ==============================================================================
# Chapter 3 — FULLY GENUINE 53/53 generation under v3.2.3-iter. Every scenario,
# every plot, real data: real-53 magnitudes (df_plots_real53) AND real un-masked
# S2 temporal shape for the DYN phenology (53/53 from Not_Masked rasters, no
# masked-series parametric fallback). STATIC = seed-injected parametric at the
# real magnitude; DYN = un-masked ATBD shape rescaled so its summer level equals
# the scenario's real magnitude. Writes out_files/Chapter3/nc_genuine53/<sc>/.
# Modes:  DRY=1 -> 3 scenarios x the 8 formerly-failing plots -> nc_genuine53_dry
#   Rscript gen_genuine53_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC))
ABL  <- list("abl_flag"='"iter"'); LY <- 2020:2022; NCORES <- 10L
DRY  <- nzchar(Sys.getenv("DRY"))

inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("missing 2020/365"); d$Julian_day<-366L; rbind(ph,d) }

# ---- real-53 df + scenario defs ----------------------------------------------
tb <- load_lai_prep(CFG_C3)$ts_by_plot; dopt <- CFG_C3$d_opt_m
df <- readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")); setDT(df)
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
all_sc <- make_all_scenarios_c3(tb, LY, dopt, "full")

# ---- un-masked ATBD daily series for ALL 53 (run_full recipe) ----------------
NM <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(NM, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs=crs(rast(files[1])))
vals <- as.data.frame(terra::extract(rast(files), pts))[,-1,drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
ts_atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
cat(sprintf("un-masked ATBD series: %d/53 plots\n", length(ts_atbd)))
stopifnot(length(ts_atbd)==53)
atbd_summer <- sapply(df$pid, function(p){ a<-ts_atbd[[p]]; mean(a$lai[a$doy>=152 & a$doy<=244], na.rm=TRUE) })

# ---- wrap each scenario genuine: STATIC=parametric(seed), DYN=un-masked shape*mag
.mag <- function(sc, prow) as.numeric(sc$lai_fn(prow))
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
SC <- lapply(all_sc, wrap_genuine)

NCROOT <- file.path(CFG_C3$out_dir, "nc_genuine53")
prob8 <- c("41_17","41_18","41_19","41_27","41_30","41_39","41_47","41_49")
if(DRY){ NCROOT<-file.path(CFG_C3$out_dir,"nc_genuine53_dry")
  SC<-SC[intersect(c("DYN_ALS","DYN_S2_DOPT","DYN_RF"),names(SC))]; idx<-which(df$id_plot %in% prob8)
  cat(sprintf("** DRY: %d dyn scenarios x %d plots **\n",length(SC),length(idx))) } else idx<-seq_len(nrow(df))
for(s in names(SC)) dir.create(file.path(NCROOT,s),recursive=TRUE,showWarnings=FALSE)
tasks <- rbindlist(lapply(names(SC), function(s) data.table(scn=s,i=idx)))
run_one <- function(k){ s<-tasks$scn[k]; i<-tasks$i[k]; prow<-as.data.frame(df[i,])
  nc<-file.path(NCROOT,s,sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000) return(sprintf("skip %s/%s",s,prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  ok<-tryCatch({ run_musica_one(prow,SC[[s]],nc,FORC,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000 },error=function(e)FALSE)
  sprintf("%s %s/%s", if(ok)"OK" else "FAIL", s, prow$id_plot) }
cat(sprintf("launching %d runs on %d cores\n", nrow(tasks), NCORES))
t0<-Sys.time(); res<-unlist(mclapply(seq_len(nrow(tasks)),run_one,mc.cores=NCORES,mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins")); tab<-table(sub(" .*","",res))
cat(sprintf("\nDONE in %.1f min | %s\n", dt, paste(names(tab),tab,sep="=",collapse=" ")))
fails<-grep("^FAIL",res,value=TRUE); if(length(fails)){cat("FAILURES:\n");cat(fails,sep="\n");cat("\n")}
