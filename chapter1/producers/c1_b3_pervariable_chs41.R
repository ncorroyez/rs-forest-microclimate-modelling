# ex-B3 (Fig B1): per-variable within-canopy gradient on CHS41. Per archetype: base +
# LAI+0.5 + cover+10 + height+1 (real LAD), full vertical output. 4x4=16 sims.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc";MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE);ABL<-list("abl_flag"='"iter"')
mk_phen<-function(l)function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}
samp<-as.data.frame(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp$P<-as.character(relabel_cluster(samp$Cluster));samp<-samp[is.finite(samp$LAI)&is.finite(samp$Hmax)&is.finite(samp$fCover),]
clfac<-make_lad_cluster_type_factory(samp)
NCDIR<-"out_files/Chapter1/nc_b3_chs41";dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
for(p in paste0("P",1:4)){sub<-samp[samp$P==p,];clu<-as.numeric(names(sort(table(sub$Cluster),decreasing=TRUE))[1])
  base<-data.frame(Cluster=clu,Hmax=median(sub$Hmax),LAI=median(sub$LAI),fCover=median(sub$fCover),P=p)
  arms<-list(base=base,
    LAIp=within(base,{LAI<-LAI+0.5}),
    fCovp=within(base,{fCover<-min(fCover+0.10,1)}),
    Hmaxp=within(base,{Hmax<-min(Hmax+1,40)}))
  for(a in names(arms)){pr<-arms[[a]];nc<-file.path(NCDIR,sprintf("%s_%s.nc",p,a))
    sc<-list(lai_fn=function(x)pr$LAI,hmax_fn=function(x)pr$Hmax,fcover_fn=function(x)pr$fCover,lad_fn=clfac,phenology_fn=mk_phen(pr$LAI))
    if(!file.exists(nc)||file.size(nc)<1000) tryCatch(run_musica_one(pr,sc,nc,FORC,MB,extra_setup=ABL),error=function(e)cat("ERR",p,a,conditionMessage(e),"\n"))
    cat(sprintf("%s %s -> %s\n",p,a,if(file.exists(nc)&&file.size(nc)>1000)"ok" else "FAIL"))}}
cat("DONE b3 sims\n")
