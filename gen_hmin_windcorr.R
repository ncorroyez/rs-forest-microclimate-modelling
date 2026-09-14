# ==============================================================================
# Chapter 3 — h_min sensitivity (revision list). Recompute the LiDAR LAI and PAD
# with the vegetation floor raised to h_min = 3 m and 5 m (baseline = 2 m Bouvier),
# force MuSICA, and test whether the dense/open reversal is robust to the floor.
# LAI(h_min) = LAI_ALS × (LAD integral above h_min / total); PAD zeroed below h_min.
# NOTE: this tests the floor UPWARD only. The understory (<2 m) hypothesis for the
# open-stratum flank is NOT testable here: the derived PAD is already floored (lowest
# bin 1.5 m, ~0), so retaining sub-2 m understory would require reprocessing the raw
# point clouds. Same windcorr harness as gen_alsk05_windcorr.R.
#   DRY=1 Rscript gen_hmin_windcorr.R  -> 3 plots check
#         Rscript gen_hmin_windcorr.R  -> STATIC_ALS_HMIN3, STATIC_ALS_HMIN5
# ==============================================================================
suppressPackageStartupMessages({library(terra);library(sf);library(ncdf4);library(lubridate);library(dplyr);library(tidyr);library(stringr);library(purrr);library(data.table);library(parallel);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R");source("R/wind_correction.R")
BIN<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE);FORC_BASE<-"in_files/musica_in_Blois_pblh.nc";ABL<-list("abl_flag"='"iter"');LY<-2020:2022;NCORES<-10L;DRY<-nzchar(Sys.getenv("DRY"))
NCROOT<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr");WFDIR<-file.path(CFG_C3$out_dir,"windcorr_forc_pblh");dir.create(WFDIR,recursive=TRUE,showWarnings=FALSE)
windcorr_pblh<-function(hmax){f<-wind_factor(pmax(hmax,2));ff<-file.path(WFDIR,sprintf("pblh_f%.2f.nc",round(f,2)));if(!file.exists(ff)||file.size(ff)<1e5){file.copy(FORC_BASE,ff,overwrite=TRUE);nc<-nc_open(ff,write=TRUE);for(v in c("Wind_E","Wind_N"))ncvar_put(nc,v,ncvar_get(nc,v)*f);nc_close(nc)};ff}
inject_seed<-function(ph){d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];if(nrow(d)!=1L)stop("seed");d$Julian_day<-366L;rbind(ph,d)}
wrap_iter<-function(sc){lai_fn<-sc$lai_fn;sc$phenology_fn<-function(pr){ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)));inject_seed(ph)};sc}
df0<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")));df0$pid<-sprintf("X%d_Y%d",round(df0$x),round(df0$y))
ladc<-grep("LAD_Layer_",names(df0),value=TRUE);z<-as.numeric(gsub("LAD_Layer_","",ladc))
.fn<-function(col)function(pr)as.numeric(pr[[col]]);lad<-make_lad_real
run_hmin<-function(hm){df<-copy(df0);below<-ladc[z<hm];tot<-rowSums(as.matrix(df[,..ladc]),na.rm=TRUE)
  for(c in below) df[[c]]<-0
  frac<-rowSums(as.matrix(df[,..ladc]),na.rm=TRUE)/pmax(tot,1e-6);df$LAI_HMIN<-df0$LAI_ALS*frac
  scn<-sprintf("STATIC_ALS_HMIN%d",hm);dir.create(file.path(NCROOT,scn),recursive=TRUE,showWarnings=FALSE)
  SC<-wrap_iter(list(name=scn,lai_fn=.fn("LAI_HMIN"),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=lad,phenology_fn=NULL))
  idx<-if(DRY) c(which.max(df$LAI_ALS),which.min(abs(df$LAI_ALS-3)),which.min(df$LAI_ALS)) else seq_len(nrow(df))
  cat(sprintf("h_min=%d m: LAI median %.2f->%.2f (open plots keep %.0f%% of LAD)\n",hm,median(df0$LAI_ALS),median(df$LAI_HMIN),100*median(frac[df0$LAI_ALS<3.86])))
  invisible(lapply(unique(df$Hmax[idx]),windcorr_pblh))
  run_one<-function(i){prow<-as.data.frame(df[i,]);nc<-file.path(NCROOT,scn,sprintf("musica_out_HOBO_%s.nc",prow$id_plot));if(file.exists(nc)&&file.size(nc)>1000)return("skip");if(file.exists(nc))file.remove(nc);forc<-windcorr_pblh(prow$Hmax)
    ok<-tryCatch({run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL);file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE);if(ok)"OK" else "FAIL"}
  t0<-Sys.time();res<-unlist(mclapply(idx,run_one,mc.cores=min(NCORES,length(idx)),mc.preschedule=FALSE))
  cat(sprintf("  %s: %.1f min | %s\n",scn,as.numeric(difftime(Sys.time(),t0,units="mins")),paste(names(table(res)),table(res),sep="=",collapse=" ")))}
for(hm in c(3,5)) run_hmin(hm)
cat("DONE h_min sensitivity\n")
