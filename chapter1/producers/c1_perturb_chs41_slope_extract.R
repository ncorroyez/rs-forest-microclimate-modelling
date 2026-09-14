# Re-extract the COUPLING-SLOPE metric levers from the existing CHS41 perturbation
# NetCDFs (nc_perturb_chs41_nowind/). No MuSICA. Gives the slope-metric attribution
# (Table G1 slope) and lets us test whether the profile leads quantity on the slope.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(parallel)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day"); Z<-1
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macH<-macH[as.Date(time)%in%ds][,.(Tmac=mean(Tmac,na.rm=TRUE)),by=time]
NCDIR<-"out_files/Chapter1/nc_perturb_chs41_nowind"
slope_of<-function(tag,id){f<-file.path(NCDIR,sprintf("%s_%s.nc",id,tag))
  if(!file.exists(f))return(NA_real_); m<-tryCatch(micro_hourly_at(f,Z),error=function(e)NULL); if(is.null(m))return(NA_real_)
  mm<-merge(m,macH,by="time")[as.Date(time)%in%ds]; if(nrow(mm)<50)return(NA_real_); as.numeric(coef(lm(Tmic~Tmac,mm))[2])}
samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[,id_plot:=sprintf("clhs_%03d",seq_len(.N))]; samp[,P:=as.character(relabel_cluster(Cluster))]
samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
D_LAI<-0.5;D_HMAX<-1.0;D_FCOV<-0.10;LAI_FLOOR<-0.1;HMAX_FLOOR<-2;HMAX_CEIL<-40;FCOV_MIN<-0.5;FCOV_MAX<-1
per<-function(sim,base,step) if(!is.finite(sim)||!is.finite(base)||abs(step)<1e-9) NA_real_ else (sim-base)/step
R<-rbindlist(mclapply(seq_len(nrow(samp)),function(i){pr<-samp[i]
  b<-slope_of("base",pr$id_plot); u<-slope_of("unifLAD",pr$id_plot)
  lp<-slope_of("LAIp",pr$id_plot);lm_<-slope_of("LAIm",pr$id_plot)
  hp<-slope_of("Hmaxp",pr$id_plot);hm<-slope_of("Hmaxm",pr$id_plot)
  fp<-slope_of("fCovp",pr$id_plot);fm<-slope_of("fCovm",pr$id_plot)
  lpv<-pr$LAI+D_LAI;lmv<-max(pr$LAI-D_LAI,LAI_FLOOR);hpv<-min(pr$Hmax+D_HMAX,HMAX_CEIL);hmv<-max(pr$Hmax-D_HMAX,HMAX_FLOOR)
  fpv<-min(pr$fCover+D_FCOV,FCOV_MAX);fmv<-max(pr$fCover-D_FCOV,FCOV_MIN)
  data.table(id_plot=pr$id_plot,P=pr$P,base=b,
    LAI_up=per(lp,b,lpv-pr$LAI),LAI_dn=per(lm_,b,lmv-pr$LAI),
    Hmax_up=per(hp,b,hpv-pr$Hmax),Hmax_dn=per(hm,b,hmv-pr$Hmax),
    fCov_up=per(fp,b,fpv-pr$fCover),fCov_dn=per(fm,b,fmv-pr$fCover),
    dT_LAD=if(is.finite(b)&&is.finite(u)) b-u else NA_real_)},mc.cores=max(1,detectCores()-2)),fill=TRUE)
fwrite(R,"out_files/Chapter1/tables/perturb_chs41_nowind_SLOPE.csv")
med<-function(v) median(as.numeric(v),na.rm=TRUE)
cat("=== SLOPE metric: profile vs LAI-step vs cover-step |effect|, by archetype (CHS41) ===\n")
for(p in c("P1","P2","P3","P4")){d<-R[P==p]
 cat(sprintf("%s: pas LAI %.4f | pas couvert %.4f | profil %.4f -> profil>LAI? %s\n",
   p,abs(med(d$LAI_up))*0.5,abs(med(d$fCov_up))*0.10,abs(med(d$dT_LAD)), abs(med(d$dT_LAD))>abs(med(d$LAI_up))*0.5))}
cat("wrote perturb_chs41_nowind_SLOPE.csv\n")
