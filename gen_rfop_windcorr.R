# ==============================================================================
# Chapter 3 — OPERATIONAL LiDAR-free empirical learner (bib archetype A, operational
# variant). Force MuSICA with the leave-one-out random-forest LAI predicted from
# Sentinel-2 + canopy-height features ONLY (no per-pixel LiDAR covariates):
#   LAIcorr_i = LOO-RF(LAI_ALS ~ S2opt + S2atbd + FORMS_H + hfrac)  [Fig G recipe]
# This is the deployable ML product (kNN/GWR/SVM/DL are interchangeable engines);
# it was so far evaluated only in LAI space (Fig G), never as a MuSICA forcing.
#
# PREDICTION: trained to recover LiDAR LAI but from operational (S2+height, LOO)
# features, it should rank like a height-driven product: usable in open canopy,
# failing in dense where S2 saturates and height alone under-ranks; i.e. behave
# between the naive S2+FORMS-H product and STATIC_RF, NOT recover the dense skill.
#
# Same windcorr harness as gen_blend_windcorr.R.
#   DRY=1 Rscript gen_rfop_windcorr.R -> check + 3 plots -> _dry
#         Rscript gen_rfop_windcorr.R -> 53 plots -> nc_genuine53_windcorr/RF_OP
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel); library(randomForest)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN  <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
FORC_BASE <- "in_files/musica_in_Blois_pblh.nc"; stopifnot(file.exists(FORC_BASE))
ABL  <- list("abl_flag" = '"iter"'); LY <- 2020:2022; NCORES <- 10L; dopt <- CFG_C3$d_opt_m
DRY  <- nzchar(Sys.getenv("DRY"))
NCROOT <- if(DRY) file.path(CFG_C3$out_dir,"nc_genuine53_windcorr_dry") else file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
dir.create(file.path(NCROOT,"RF_OP"), recursive=TRUE, showWarnings=FALSE)
WFDIR  <- file.path(CFG_C3$out_dir, "windcorr_forc_pblh"); dir.create(WFDIR, recursive=TRUE, showWarnings=FALSE)
windcorr_pblh <- function(hmax){ f<-wind_factor(pmax(hmax,2)); ff<-file.path(WFDIR,sprintf("pblh_f%.2f.nc",round(f,2)))
  if(!file.exists(ff)||file.size(ff)<1e5){file.copy(FORC_BASE,ff,overwrite=TRUE)
    nc<-nc_open(ff,write=TRUE);for(v in c("Wind_E","Wind_N"))ncvar_put(nc,v,ncvar_get(nc,v)*f);nc_close(nc)}; ff }
inject_seed<-function(ph){d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];if(nrow(d)!=1L)stop("seed");d$Julian_day<-366L;rbind(ph,d)}
wrap_iter<-function(sc){lai_fn<-sc$lai_fn;sc$phenology_fn<-function(pr){
  ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)));inject_seed(ph)};sc}

df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
# operational S2 features: summer means of opt & atbd smoothed series (Fig G recipe)
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"; pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
summ<-function(pat){f<-list.files(NM,pattern=pat,full.names=TRUE);dt<-as.integer(format(as.Date(str_extract(basename(f),"\\d{4}-\\d{2}-\\d{2}")),"%j"))
  v<-as.data.frame(terra::extract(rast(f),pts))[,-1,drop=FALSE]
  long<-rbindlist(lapply(seq_along(dt),function(i)data.table(plot_id=df$pid,doy=dt[i],lai=pmax(v[[i]],0))))
  ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3);sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
df$S2opt <-summ("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$")
df$S2atbd<-summ("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$")
df$hfrac <-pmax(df$FORMS_H,0.1)/dopt
Dl<-as.data.frame(df[,.(LAI_ALS,S2opt,S2atbd,FORMS_H,hfrac)])
set.seed(42); df$LAI_RFOP<-sapply(1:nrow(Dl),function(i)predict(randomForest(LAI_ALS~S2opt+S2atbd+FORMS_H+hfrac,Dl[-i,],ntree=400),Dl[i,]))
df$LAI_RFOP<-pmax(df$LAI_RFOP,0)
cat(sprintf("LOO operational RF: cor(LAI_RFOP,LAI_ALS)=%.2f  range[%.2f..%.2f]\n",cor(df$LAI_RFOP,df$LAI_ALS),min(df$LAI_RFOP),max(df$LAI_RFOP)))

.fn<-function(col)function(pr)as.numeric(pr[[col]]); lad<-make_lad_real
SC<-wrap_iter(list(name="RF_OP",lai_fn=.fn("LAI_RFOP"),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=lad,phenology_fn=NULL))
if(DRY){o<-order(df$LAI_ALS);idx<-unique(c(o[1],o[ceiling(nrow(df)/2)],o[nrow(df)]))} else idx<-seq_len(nrow(df))
print(df[idx,.(id_plot,LAI_ALS=round(LAI_ALS,2),FORMS_H=round(FORMS_H,1),LAI_RFOP=round(LAI_RFOP,2))])
invisible(lapply(unique(df$Hmax[idx]),windcorr_pblh))
run_one<-function(i){prow<-as.data.frame(df[i,]);nc<-file.path(NCROOT,"RF_OP",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)&&file.size(nc)>1000)return(sprintf("skip %s",prow$id_plot));if(file.exists(nc))file.remove(nc)
  forc<-windcorr_pblh(prow$Hmax)
  ok<-tryCatch({run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL);file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE)
  sprintf("%s %s",if(ok)"OK" else "FAIL",prow$id_plot)}
cat(sprintf("\nRF_OP windcorr%s: %d plots\n",if(DRY)" [DRY]" else "",length(idx)))
t0<-Sys.time();res<-unlist(mclapply(idx,run_one,mc.cores=min(NCORES,length(idx)),mc.preschedule=FALSE))
dt<-as.numeric(difftime(Sys.time(),t0,units="mins"));cat(sprintf("DONE %.2f min | %s\n",dt,paste(names(table(sub(' .*','',res))),table(sub(' .*','',res)),sep='=',collapse=' ')))
