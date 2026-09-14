# ==============================================================================
# Does v3.2.3 actually RANK plots via their between-plot LAD shape? Run STATIC_ALS
# with UNIFORM LAD (magnitude preserved, vertical shape removed) under both
# binaries, to compare against the existing REAL-LAD runs. Decompose:
#   - corr(real ranking, obs) vs corr(uniform ranking, obs): does LAD shape add
#     ranking skill, per binary?
#   - per-plot (real - uniform) = the LAD-shape effect on ΔTmax (c1 sensitivity);
#     is it bigger in v3.2.3? does it correlate with the obs gradient (signal) or not (noise)?
#   Rscript run_lad_uniform_test.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
prep <- load_lai_prep(CFG_C3); df <- as.data.table(prep$df_plots)
.lf <- function(col) function(pr) as.numeric(pr[[col]])
LY <- 2020:2022
inject_seed <- function(ph){ d<-ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]; d$Julian_day<-366L; rbind(ph,d) }
mk_stat <- function(col) function(pr) inject_seed(calc_phenology(list.year=LY, nleafage=1, budburst_date=115,
  leaf_age_max_in=0.56, relative_age_firstmax=0.10, relative_age_lastmax=0.75, LAI_max_per_cohort=as.numeric(pr[[col]])))

runs <- list(
  list(tag="STATIC_ALS_unifLAD", root="nc",          bin=CFG_C3$musica_cmd,
       forc=CFG_C3$forcing_file, abl=list(), pheno=NULL),                                  # v3.2.0 uniform
  list(tag="STATIC_ALS_unifLAD", root="nc_v323iter", bin=normalizePath("in_files/model-3.2.3/musica"),
       forc="in_files/musica_in_Blois_pblh.nc", abl=list("abl_flag"='"iter"'), pheno=mk_stat("LAI_ALS")))  # v3.2.3 iter uniform

for (R in runs) {
  d <- file.path(CFG_C3$out_dir, R$root, R$tag); dir.create(d, recursive=TRUE, showWarnings=FALSE)
  SC <- list(name=R$tag, lai_fn=.lf("LAI_ALS"), hmax_fn=.lf("Hmax"), fcover_fn=.lf("fCover"),
             lad_fn=make_lad_uniform, phenology_fn=R$pheno)
  run_one <- function(i){ prow<-as.data.frame(df[i,]); nc<-file.path(d,sprintf("musica_out_HOBO_%s.nc",prow$id_plot))
    if (file.exists(nc)&&file.size(nc)>1000) return("skip"); if(file.exists(nc)) file.remove(nc)
    ok<-tryCatch({run_musica_one(prow,SC,nc,R$forc,R$bin,extra_setup=R$abl); file.exists(nc)&&file.size(nc)>1000},error=function(e)FALSE)
    if(ok)"OK" else "FAIL" }
  res<-unlist(mclapply(seq_len(nrow(df)), run_one, mc.cores=10L, mc.preschedule=FALSE))
  cat(sprintf("[%s @ %s] %s\n", R$tag, R$root, paste(names(table(res)),table(res),sep="=",collapse=" ")))
}
cat("DONE uniform-LAD runs\n")
