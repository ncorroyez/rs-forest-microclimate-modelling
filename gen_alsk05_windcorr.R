# ==============================================================================
# Chapter 3 — k-robustness: re-force MuSICA with LiDAR LAI recomputed at k = 0.5
# instead of the k = 0.65 used throughout (companion/Ch2 value). LiDAR LAI from
# Beer-Lambert gap fraction is ∝ 1/k, so LAI(k=0.5) = LAI(k=0.65) × (0.65/0.5) =
# ×1.30 (a constant rescale; LAD profile shape unchanged). Tests whether the
# dense/open reversal and the switch survive the k choice Nathan flagged (Ch1
# reports k affects MuSICA absolute skill).
#   Rscript gen_alsk05_windcorr.R -> nc_genuine53_windcorr/STATIC_ALS_K05
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); FORC_BASE<-"in_files/musica_in_Blois_pblh.nc"
ABL<-list("abl_flag"='"iter"'); LY<-2020:2022; NCORES<-10L
NCROOT<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); dir.create(file.path(NCROOT,"STATIC_ALS_K05"),recursive=TRUE,showWarnings=FALSE)
WFDIR<-file.path(CFG_C3$out_dir,"windcorr_forc_pblh"); dir.create(WFDIR,recursive=TRUE,showWarnings=FALSE)
windcorr_pblh<-function(hmax){f<-wind_factor(pmax(hmax,2));ff<-file.path(WFDIR,sprintf("pblh_f%.2f.nc",round(f,2)))
  if(!file.exists(ff)||file.size(ff)<1e5){file.copy(FORC_BASE,ff,overwrite=TRUE);nc<-nc_open(ff,write=TRUE);for(v in c("Wind_E","Wind_N"))ncvar_put(nc,v,ncvar_get(nc,v)*f);nc_close(nc)};ff}
inject_seed<-function(ph){d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];if(nrow(d)!=1L)stop("seed");d$Julian_day<-366L;rbind(ph,d)}
wrap_iter<-function(sc){lai_fn<-sc$lai_fn;sc$phenology_fn<-function(pr){ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)));inject_seed(ph)};sc}
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
df$LAI_ALS_K05<-df$LAI_ALS*(0.65/0.50)   # k 0.65 -> 0.5, LAI ∝ 1/k
cat(sprintf("LAI_ALS k=0.65 median %.2f -> k=0.5 median %.2f (×1.30)\n",median(df$LAI_ALS),median(df$LAI_ALS_K05)))
.fn<-function(col)function(pr)as.numeric(pr[[col]]); lad<-make_lad_real
SC<-wrap_iter(list(name="STATIC_ALS_K05",lai_fn=.fn("LAI_ALS_K05"),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=lad,phenology_fn=NULL))
invisible(lapply(unique(df$Hmax),windcorr_pblh))
run_one<-function(i){prow<-as.data.frame(df[i,]);nc<-file.path(NCROOT,"STATIC_ALS_K05",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000)return("skip");if(file.exists(nc))file.remove(nc);forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL);file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE);if(ok)"OK" else "FAIL"}
t0<-Sys.time();res<-unlist(mclapply(seq_len(nrow(df)),run_one,mc.cores=NCORES,mc.preschedule=FALSE))
cat(sprintf("DONE %.1f min | %s\n",as.numeric(difftime(Sys.time(),t0,units="mins")),paste(names(table(res)),table(res),sep="=",collapse=" ")))
