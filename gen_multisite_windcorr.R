# ==============================================================================
# Chapter 3 — MULTI-SITE TRAINING of the operational LiDAR-free learner.
# Train RF (lidarlai ~ S2atbd + S2opt + FORMS-H + hfrac) on wall-to-wall pixels,
# predict LAI at the 53 Blois logger plots, force MuSICA, validate at Blois.
# Three arms, TOTAL n held fixed (6000) so only COMPOSITION differs:
#   BLOIS  : 6000 Blois pixels (plot pixels excluded)           -> RF_MS_BLOIS
#   ALL3   : 2000 each Blois+Aigoual+Mormal                     -> RF_MS_ALL3
#   LBO    : 3000 each Aigoual+Mormal (leave-Blois-out)         -> RF_MS_LBO
# Validation is Blois-only (loggers): this shows whether training composition
# improves the driver AT BLOIS; it cannot show the product works off-Blois.
# PREDICTION: S2 saturation is physical -> dense stratum stays ~0.10 in all arms;
# dense R2 > ~0.2 would signal leakage or FORMS-H carrying the signal.
#   Rscript gen_multisite_windcorr.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(sf); library(ncdf4); library(lubridate); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(data.table); library(parallel); library(randomForest)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R"); source("R/wind_correction.R")
BIN<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); FORC_BASE<-"in_files/musica_in_Blois_pblh.nc"
ABL<-list("abl_flag"='"iter"'); LY<-2020:2022; NCORES<-10L; dopt<-CFG_C3$d_opt_m
NCROOT<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr")
WFDIR<-file.path(CFG_C3$out_dir,"windcorr_forc_pblh"); dir.create(WFDIR,recursive=TRUE,showWarnings=FALSE)
RES<-"/home/corroyez/Documents/NC_Full/03_RESULTS"; DAT<-"/home/corroyez/Documents/NC_Full/01_DATA"; FR<-rast(file.path(DAT,"FORMS-H_Height_10m_cm.tif"))
windcorr_pblh<-function(hmax){f<-wind_factor(pmax(hmax,2));ff<-file.path(WFDIR,sprintf("pblh_f%.2f.nc",round(f,2)))
  if(!file.exists(ff)||file.size(ff)<1e5){file.copy(FORC_BASE,ff,overwrite=TRUE);nc<-nc_open(ff,write=TRUE);for(v in c("Wind_E","Wind_N"))ncvar_put(nc,v,ncvar_get(nc,v)*f);nc_close(nc)};ff}
inject_seed<-function(ph){d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];if(nrow(d)!=1L)stop("seed");d$Julian_day<-366L;rbind(ph,d)}
wrap_iter<-function(sc){lai_fn<-sc$lai_fn;sc$phenology_fn<-function(pr){ph<-calc_phenology(list.year=LY,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=as.numeric(lai_fn(pr)));inject_seed(ph)};sc}
smean<-function(files){doy<-as.integer(format(as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}")),"%j"));app(rast(files[doy>=152&doy<=244]),fun=function(x)mean(pmax(x,0),na.rm=TRUE))}
site_stack<-function(S){
  lidar<-clamp(rast(file.path(RES,S,"Metrics/Not_Masked/lidarlai_res_10_m.tif")),0,12,values=FALSE)  # drop only clear artefacts (>12), KEEP dense exemplars (8-12) for training
  atbd<-smean(list.files(file.path(RES,S,"Metrics/Not_Masked"),pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$",full.names=TRUE))
  opt <-smean(list.files(file.path(RES,S,"Metrics/Not_Masked"),pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$",full.names=TRUE))
  if(S=="Blois"){fh<-rast(file.path(DAT,"FORMS-H_Blois.tif"))} else {e<-project(as.polygons(ext(lidar),crs=crs(lidar)),crs(FR));fh<-crop(FR,ext(e)+200)/100}
  atbd<-project(atbd,lidar);opt<-project(opt,lidar);fh<-project(fh,lidar)  # project handles CRS+grid
  st<-c(lidar,atbd,opt,fh,fh/dopt); names(st)<-c("lidar","S2atbd","S2opt","FORMS_H","hfrac"); st}
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
cat("building site stacks...\n"); ST<-setNames(lapply(c("Blois","Aigoual","Mormal"),site_stack),c("Blois","Aigoual","Mormal"))
# sample n complete-case forest pixels per site (Blois: exclude the 53 plot cells)
samp<-function(S,n,excl_cells=NULL){st<-ST[[S]];v<-as.data.table(values(st));v[,cell:=.I];v<-v[complete.cases(v)&lidar>0.3]
  if(!is.null(excl_cells))v<-v[!cell%in%excl_cells]; set.seed(match(S,c("Blois","Aigoual","Mormal")));v[sample(.N,min(n,.N))]}
blois_cells<-unique(na.omit(terra::extract(ST[["Blois"]][[1]],terra::buffer(pts,45),cells=TRUE)$cell))  # exclude a 45 m buffer around plots (not just the cell) to kill adjacency leakage
cat(sprintf("leakage buffer: %d Blois cells excluded (45 m around 53 plots)\n",length(blois_cells)))
tr<-list(BLOIS=samp("Blois",6000,blois_cells),
         ALL3=rbind(samp("Blois",2000,blois_cells),samp("Aigoual",2000),samp("Mormal",2000)),
         LBO =rbind(samp("Aigoual",3000),samp("Mormal",3000)))
feat<-as.data.frame(terra::extract(ST[["Blois"]],pts))[,c("S2atbd","S2opt","FORMS_H","hfrac")]
for(cc in names(feat)) feat[[cc]][is.na(feat[[cc]])]<-median(feat[[cc]],na.rm=TRUE)  # impute rare cloud/edge NA at plots
fml<-lidar~S2atbd+S2opt+FORMS_H+hfrac
for(a in names(tr)){set.seed(42);rf<-randomForest(fml,data=tr[[a]],ntree=400)
  df[[paste0("LAI_",a)]]<-pmax(as.numeric(predict(rf,feat)),0)
  cat(sprintf("arm %-6s n=%d  OOB R2=%.2f  pred@plots cor(LAI_ALS)=%.2f [%.1f..%.1f]\n",a,nrow(tr[[a]]),1-rf$mse[400]/var(tr[[a]]$lidar),cor(df[[paste0("LAI_",a)]],df$LAI_ALS),min(df[[paste0("LAI_",a)]]),max(df[[paste0("LAI_",a)]])))}
# MuSICA per arm
.fn<-function(col)function(pr)as.numeric(pr[[col]]); lad<-make_lad_real
for(a in names(tr)){scn<-paste0("RF_MS2_",a); dir.create(file.path(NCROOT,scn),recursive=TRUE,showWarnings=FALSE)
  SC<-wrap_iter(list(name=scn,lai_fn=.fn(paste0("LAI_",a)),hmax_fn=.fn("Hmax"),fcover_fn=.fn("fCover"),lad_fn=lad,phenology_fn=NULL))
  invisible(lapply(unique(df$Hmax),windcorr_pblh))
  run_one<-function(i){prow<-as.data.frame(df[i,]);nc<-file.path(NCROOT,scn,sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
    if(file.exists(nc)&&file.size(nc)>1000)return("skip");if(file.exists(nc))file.remove(nc);forc<-windcorr_pblh(prow$Hmax)
    ok<-tryCatch({run_musica_one(prow,SC,nc,forc,BIN,extra_setup=ABL);file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE);if(ok)"OK" else "FAIL"}
  t0<-Sys.time();res<-unlist(mclapply(seq_len(nrow(df)),run_one,mc.cores=NCORES,mc.preschedule=FALSE))
  cat(sprintf("  %s: %.1f min | %s\n",scn,as.numeric(difftime(Sys.time(),t0,units="mins")),paste(names(table(res)),table(res),sep="=",collapse=" ")))}
cat("DONE multisite arms\n")
