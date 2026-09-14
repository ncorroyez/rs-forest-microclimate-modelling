# ==============================================================================
# Chapter 3 — complete the CHS41 / no-wind tree for the two dynamic Sentinel-2
# scenarios. Eight plots never produced a usable NetCDF there (six missing, two
# zero-byte): the "prob8" set of gen_genuine53_v323.R. Same recipe as that
# generator, but the CHS41-Rmerge forcing and NO wind correction, matching
# gen_opt_CHS41_nowind.R and the rest of out_files/Chapter3_CHS41/nc.
#   PLOTS=41_17 Rscript gen_fix_CHS41_dyn.R   -> one plot, timing check
#           Rscript gen_fix_CHS41_dyn.R       -> every missing/undersized plot
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config_CHS41.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork=TRUE)
FORC <- CFG_C3$forcing_file; stopifnot(file.exists(FORC))
cat("forcing:", FORC, "\n")
ABL  <- list("abl_flag"='"iter"'); LY <- 2020:2022
WANT <- c("DYN_S2_ATBD","DYN_S2_RESCALED")
MINSZ <- 1e6

inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]
  if(nrow(d)!=1L) stop("missing 2020/365"); d$Julian_day<-366L; rbind(ph,d) }

tb <- load_lai_prep(CFG_C3)$ts_by_plot; dopt <- CFG_C3$d_opt_m
df <- readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")); setDT(df)
df$pid <- sprintf("X%d_Y%d", round(df$x), round(df$y))
all_sc <- make_all_scenarios_c3(tb, LY, dopt, "full")
stopifnot(all(WANT %in% names(all_sc)))

NM <- "/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
files <- list.files(NM, pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$", full.names=TRUE)
dates <- as.Date(str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs=crs(rast(files[1])))
vals <- as.data.frame(terra::extract(rast(files), pts))[,-1,drop=FALSE]
long <- rbindlist(lapply(seq_along(dates), function(i)
  data.table(plot_id=df$pid, date=dates[i], doy=as.integer(format(dates[i],"%j")), lai=pmax(vals[[i]],0))))
ts_atbd <- smooth_s2_ts(as.data.frame(long), k=8, min_obs=3)
stopifnot(length(ts_atbd)==53)
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
SC <- lapply(all_sc[WANT], wrap_genuine)
# Read the existing tree, WRITE to a separate one: the failed runs left zero-byte
# files behind and MuSICA will not overwrite them. Nothing is deleted here.
NCROOT     <- file.path(CFG_C3$out_dir, "nc")       # what already exists
NCROOT_OUT <- file.path(CFG_C3$out_dir, "nc_fix")   # where the completions go

only <- Sys.getenv("PLOTS"); only <- if(nzchar(only)) strsplit(only,",")[[1]] else NULL
todo <- rbindlist(lapply(names(SC), function(s){
  ids <- df$id_plot
  if(!is.null(only)) ids <- ids[ids %in% only]
  keep <- vapply(ids, function(id){
    f  <- file.path(NCROOT,s,sprintf("musica_out_HOBO_%s.nc",id))
    f2 <- file.path(NCROOT_OUT,s,sprintf("musica_out_HOBO_%s.nc",id))
    ok  <- file.exists(f)  && file.size(f)  >= MINSZ
    ok2 <- file.exists(f2) && file.size(f2) >= MINSZ
    !(ok || ok2) }, logical(1))
  if(!any(keep)) return(NULL)
  data.table(scn=s, id_plot=ids[keep]) }))
if(!nrow(todo)) { cat("nothing to do\n"); quit(save="no") }
cat(sprintf("to run: %d\n", nrow(todo))); print(todo)

run_one <- function(k){
  s <- todo$scn[k]; id <- todo$id_plot[k]
  prow <- as.data.frame(df[id_plot==id,])
  nc <- file.path(NCROOT_OUT,s,sprintf("musica_out_HOBO_%s.nc",id))
  dir.create(dirname(nc), recursive=TRUE, showWarnings=FALSE)
  err <- NULL
  ok <- tryCatch({ run_musica_one(prow,SC[[s]],nc,FORC,BIN,extra_setup=ABL)
                   file.exists(nc) && file.size(nc) > MINSZ },
                 error=function(e){ err <<- conditionMessage(e); FALSE })
  sprintf("%s %s/%s (%.2f MB)%s", if(ok)"OK" else "FAIL", s, id,
          if(file.exists(nc)) file.size(nc)/1e6 else 0,
          if(is.null(err)) "" else paste0("  ERR: ", err))
}
NCORES <- min(8L, nrow(todo))
cat(sprintf("launching %d runs on %d cores\n", nrow(todo), NCORES))
t0 <- Sys.time()
res <- unlist(mclapply(seq_len(nrow(todo)), run_one, mc.cores=NCORES, mc.preschedule=FALSE))
cat(sprintf("DONE in %.1f min\n", as.numeric(difftime(Sys.time(),t0,units="mins"))))
cat(res, sep="\n"); cat("\n")
