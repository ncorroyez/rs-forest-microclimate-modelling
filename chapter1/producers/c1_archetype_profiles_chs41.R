# ==============================================================================
# c1_archetype_profiles_chs41.R — per-archetype REAL vs UNIFORM sims on CHS41-Rmerge,
# no wind, for the within-canopy vertical gradients (Fig B2 -> main Fig 7; Fig B3).
# The shipped B2/B3 read LEGACY archetype sims; this produces them on the chapter's
# authoritative forcing. Representative config per archetype = median LAI/Hmax/fCover;
# real = archetype-mean LAD shape (make_lad_cluster_type_factory), uniform = uniform LAD.
# 4 archetypes x 2 = 8 sims, full-canopy vertical output (Tair_z, wind_z, wair_z, ...).
# Out: out_files/Chapter1/nc_archetype_chs41/<P>_<real|unif>.nc
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"
MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"')
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}

samp<-as.data.frame(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp$Cluster<-samp$Cluster; samp$P<-as.character(relabel_cluster(samp$Cluster))
samp<-samp[is.finite(samp$LAI)&is.finite(samp$Hmax)&is.finite(samp$fCover),]
# archetype-mean LAD shape constructor (rescaled to each config's LAI/Hmax)
clfac<-make_lad_cluster_type_factory(samp)
NCDIR<-"out_files/Chapter1/nc_archetype_chs41"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)

for(p in paste0("P",1:4)){
  sub<-samp[samp$P==p,]
  clu<-as.numeric(names(sort(table(sub$Cluster),decreasing=TRUE))[1])   # dominant raw cluster id for this P
  prow<-data.frame(Cluster=clu, Hmax=median(sub$Hmax), LAI=median(sub$LAI), fCover=median(sub$fCover), P=p)
  cat(sprintf("%s: LAI=%.2f Hmax=%.1f fCover=%.2f (n=%d)\n",p,prow$LAI,prow$Hmax,prow$fCover,nrow(sub)))
  for(arm in c("real","unif")){
    ladf<-if(arm=="real") clfac else make_lad_uniform
    nc<-file.path(NCDIR,sprintf("%s_%s.nc",p,arm))
    sc<-list(lai_fn=function(pp)prow$LAI,hmax_fn=function(pp)prow$Hmax,fcover_fn=function(pp)prow$fCover,
             lad_fn=ladf,phenology_fn=mk_phen(prow$LAI))
    if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(prow,sc,nc,FORC,MB,extra_setup=ABL),error=function(e)cat("ERR",p,arm,conditionMessage(e),"\n"))
    cat(sprintf("   %s %s -> %s (%s)\n",p,arm,basename(nc),if(file.exists(nc)&&file.size(nc)>1000)"ok" else "FAIL"))
  }
}
cat("DONE archetype sims\n")
