# TEST: does a UNIFORM LAI rescale (k=0.5->0.65, x0.769) change the attribution
# ranking (LAI vs profile) or density trend (P1->P4)? BOTH metrics: dTmax + log-slope.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(rmusica);library(musica.tools)
  src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
  source("pipeline/00_config.R");source("R/cluster_relabel.R")})
FORC<-"in_files/FR-Blo_2021_v2.nc"; MB<-normalizePath("in_files/model-3.2.3/musica",mustWork=TRUE); ABL<-list("abl_flag"='"iter"'); KS<-0.5/0.65
mk_phen<-function(l){force(l);function(p){ph<-as.data.frame(calc_phenology(list.year=2020:2022,nleafage=1,budburst_date=115,leaf_age_max_in=0.56,relative_age_firstmax=0.10,relative_age_lastmax=0.75,LAI_max_per_cohort=l));d<-ph[ph$year==2020&ph$Julian_day==365,,drop=FALSE];d$Julian_day<-366;rbind(ph,d)}}
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1; SH<-2L
dm<-as.data.table(extract_macro_daily(FORC,ds))
# macro HOURLY (for slope, no shift)
fn<-nc_open(FORC);ftu<-ncatt_get(fn,"time","units")$value;ft0<-as.POSIXct(sub(".*since ","",ftu),tz="UTC");funit<-if(grepl("^seconds",ftu))"s" else "h"
fth<-ncvar_get(fn,"time");fT<-ncvar_get(fn,"Tair")-273.15;nc_close(fn)
ftime<-if(funit=="s") ft0+fth else ft0+dhours(fth)
macH<-data.table(th=floor_date(ftime,"hour"),Tmac=fT)
# extract dTmax + log-slope from one nc
mets<-function(nc){if(!file.exists(nc)||file.size(nc)<1e6)return(c(dTmax=NA,logslope=NA));n<-nc_open(nc);tu<-ncatt_get(n,"time","units")$value;t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC");th<-ncvar_get(n,"time");Tk<-ncvar_get(n,"Tair_z");rh<-ncvar_get(n,"relative_height");vh<-median(ncvar_get(n,"veget_height_top"),na.rm=T);nc_close(n)
  zl<-rh*vh;if(Z<=zl[1]){il<-1;ih<-1;w<-0}else{il<-max(which(zl<=Z));ih<-il+1;w<-(Z-zl[il])/(zl[ih]-zl[il])};Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
  tvraw<-t0+dhours(th)
  # dTmax: -2h shift, daily max
  dd<-data.table(date=as.Date(floor_date(tvraw-hours(SH),"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date];dtm<-merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=T)]
  # slope: no shift, hourly micro vs macro
  mm<-merge(data.table(th=floor_date(tvraw,"hour"),Tmic=Tc),macH,by="th");sl<-if(nrow(mm)>10)coef(lm(Tmic~Tmac,mm))[2] else NA
  c(dTmax=dtm,logslope=if(is.finite(sl)&&sl>0)log(sl) else NA)}
samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"));samp[,P:=as.character(relabel_cluster(Cluster))]
ladc<-grep("^LAD_Layer_",names(samp),value=T);hz<-as.numeric(sub("LAD_Layer_","",ladc))
LV<-c("P1","P2","P3","P4");sd_lai<-c(P1=0.70,P2=0.76,P3=0.62,P4=0.64);sd_vci<-c(P1=0.26,P2=0.06,P3=0.07,P4=0.04);vci<-c(P1=0.580,P2=0.835,P3=0.849,P4=0.864)
NCDIR<-"out_files/Chapter1/nc_ktest";dir.create(NCDIR,recursive=T,showWarnings=F)
mkl<-function(prof,lai){d<-prof;if(sum(d)>0)d<-d*lai/sum(d);data.frame(height=seq_along(d),density=d)}
run1<-function(lai,hmax,fcov,laddf,forc=FORC){sc<-list(lai_fn=function(pr)lai,hmax_fn=function(pr)hmax,fcover_fn=function(pr)max(fcov,0.5),lad_fn=function(plot_row,hmax=NULL,lai=NULL)laddf,phenology_fn=mk_phen(lai));nc<-tempfile(fileext=".nc",tmpdir=NCDIR);tryCatch(run_musica_one(data.frame(Cluster=1,x=0,y=0),sc,nc,forc,MB,extra_setup=ABL),error=function(e)NULL);m<-mets(nc);if(file.exists(nc))unlink(nc);m}
rows<-list()
for(cl in LV){cc<-samp[P==cl];L0<-mean(cc$LAI);H0<-mean(cc$Hmax);fc0<-mean(cc$fCover);ml<-colMeans(as.matrix(cc[,..ladc]),na.rm=T);ml[is.na(ml)]<-0;ml<-ml[hz<=H0]
  dep<-1-vci[cl]  # departure from uniform for the profile normalization
  fcw<-windcorr_forcing(H0)  # wind-corrected baseline forcing at the cluster base height (micro only; macro ref stays raw)
  for(ks in c(1.0,KS)){L<-L0*ks;sdL<-sd_lai[cl]*ks
    ref<-run1(L,H0,fc0,mkl(ml,L),forc=fcw);up<-run1(L+sdL,H0,fc0,mkl(ml,L+sdL),forc=fcw);dn<-run1(L-sdL,H0,fc0,mkl(ml,L-sdL),forc=fcw);unif<-run1(L,H0,fc0,data.frame(height=seq_len(ceiling(H0)),density=L/ceiling(H0)),forc=fcw)
    for(met in c("dTmax","logslope")){
      sL<-(up[met]-dn[met])/2                                  # LAI sensitivity per SD
      pr<-(ref[met]-unif[met])*sd_vci[cl]/dep                  # profile contrast per SD of (1-VCI) [FIXED]
      rows[[length(rows)+1]]<-data.table(P=cl,k=ifelse(ks==1,"k0.5","k0.65"),metric=met,LAI=round(L,2),sens_LAI=round(sL,3),prof=round(pr,3))}
    cat(sprintf("%s %s done (LAI=%.2f)\n",cl,ifelse(ks==1,"k0.5","k0.65"),L)) }}
R<-rbindlist(rows);fwrite(R,"out_files/Chapter1/tables/tab_k_rescale_test.csv")
cat("\n=== RESULT (sens_LAI vs prof, par métrique, k0.5 vs k0.65) ===\n")
for(met in c("dTmax","logslope")){cat("\n--",met,"--\n");print(dcast(R[metric==met],P~k,value.var=c("sens_LAI","prof")))}
cat("\nLAI domine le profil partout (|sensLAI|>=|prof|) ? ",all(abs(R$sens_LAI)>=abs(R$prof)-1e-6),"\n")
