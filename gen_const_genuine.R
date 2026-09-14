# ==============================================================================
# Regenerate CONST_ALS with a TRUE constant phenology (LAI_ALS flat all year =
# the naive fixed-LAI baseline), genuine-53 v3.2.3, into nc_genuine53/CONST_ALS.
# (gen_genuine53's wrap_genuine wrongly gave it the ATBD shape.)
#   Rscript gen_const_genuine.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(data.table);library(dplyr);library(stringr);library(parallel);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
BIN<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); FORC<-"in_files/musica_in_Blois_pblh.nc"; ABL<-list("abl_flag"='"iter"'); LY<-2020:2022
inject_seed<-function(ph){d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366L;rbind(ph,d)}
df<-as.data.frame(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
NCROOT<-file.path(CFG_C3$out_dir,"nc_genuine53"); dir.create(file.path(NCROOT,"CONST_ALS"),recursive=TRUE,showWarnings=FALSE)
# CONST scenario: LAI_ALS magnitude, FLAT year-round phenology, full LAD
SC<-list(name="CONST_ALS", lai_fn=function(pr)as.numeric(pr$LAI_ALS), hmax_fn=function(pr)as.numeric(pr$Hmax),
  fcover_fn=function(pr)as.numeric(pr$fCover), lad_fn=make_lad_real,
  phenology_fn=function(pr){ lai<-as.numeric(pr$LAI_ALS)
    inject_seed(make_phenology_from_s2(data.frame(doy=1:365, lai=rep(lai,365)), LY)) })
run_one<-function(i){ prow<-df[i,]; nc<-file.path(NCROOT,"CONST_ALS",sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
  if(file.exists(nc)) file.remove(nc)
  ok<-tryCatch({run_musica_one(prow,SC,nc,FORC,BIN,extra_setup=ABL); file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE)
  sprintf("%s %s", if(ok)"OK" else "FAIL", prow$id_plot) }
res<-unlist(mclapply(seq_len(nrow(df)),run_one,mc.cores=10,mc.preschedule=FALSE))
cat(sprintf("CONST_ALS regenerated: %s\n", paste(names(table(sub(" .*","",res))),table(sub(" .*","",res)),sep="=",collapse=" ")))
