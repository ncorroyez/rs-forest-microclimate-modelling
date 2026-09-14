# Hot-days (hottest 10%) ΔTmax attribution + inline stats, CHS41, from existing perturb nc.
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(data.table);library(parallel)})
src<-list.files("R","\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source))
FORC<-"out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc";ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day");Z<-1
nc<-nc_open(FORC);tu<-ncatt_get(nc,"time","units")$value;th<-ncvar_get(nc,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
macH<-data.table(time=floor_date(t0+th*3600,"hour"),Tmac=as.numeric(ncvar_get(nc,"Tair"))-273.15);nc_close(nc)
macD<-macH[as.Date(time)%in%ds,.(Tmx=max(Tmac,na.rm=T)),by=.(date=as.Date(time))]
hot<-macD[order(-Tmx)][1:ceiling(.N*0.10)]$date
MREF<-macro_ref(FORC,hot)                                    # <-- FIX: macro ref built on hot days
NCDIR<-"out_files/Chapter1/nc_perturb_chs41_nowind"
dth<-function(tag,id){f<-file.path(NCDIR,sprintf("%s_%s.nc",id,tag));if(!file.exists(f))return(NA_real_)
  m<-tryCatch(micro_hourly_at(f,Z),error=function(e)NULL);if(is.null(m))return(NA_real_);delta_tmax_mean(m,MREF,hot,min_days=5)}
samp<-as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp[,id_plot:=sprintf("clhs_%03d",seq_len(.N))];samp[,P:=as.character(relabel_cluster(Cluster))];samp<-samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
per<-function(s,b,st) if(!is.finite(s)||!is.finite(b)||abs(st)<1e-9) NA_real_ else (s-b)/st
R<-rbindlist(mclapply(seq_len(nrow(samp)),function(i){pr<-samp[i];b<-dth("base",pr$id_plot);u<-dth("unifLAD",pr$id_plot)
  data.table(P=pr$P,LAI_up=per(dth("LAIp",pr$id_plot),b,0.5),LAI_dn=per(dth("LAIm",pr$id_plot),b,-0.5),
    fCov_up=per(dth("fCovp",pr$id_plot),b,0.10),Hmax_up=per(dth("Hmaxp",pr$id_plot),b,1),
    dT_LAD=if(is.finite(b)&&is.finite(u))b-u else NA)},mc.cores=8),fill=TRUE)
fwrite(R,"out_files/Chapter1/tables/hotdays_chs41.csv")
med<-function(v) median(as.numeric(v),na.rm=T)
cat(sprintf("hot days: %d (Tmax>=%.1f)\n",length(hot),min(macD[date%in%hot]$Tmx)))
cat("=== HOT-DAY levers per step, by archetype ===\n")
for(p in c("P1","P2","P3","P4")){d<-R[P==p]
 cat(sprintf("%s: LAI+0.5 %.3f | fCov+10 %.3f | Hmax+1 %.4f | profile %.3f\n",p,med(d$LAI_up)*0.5,med(d$fCov_up)*0.10,med(d$Hmax_up),med(d$dT_LAD)))}
# full-summer profile median for the S1 comparison (from perturb_chs41)
S<-fread("out_files/Chapter1/tables/perturb_chs41_nowind.csv")
cat("\n=== profile median: full summer vs hot days, P2/P3 ===\n")
for(p in c("P2","P3")) cat(sprintf("%s: full %.3f -> hot %.3f\n",p,med(S[P==p]$dT_LAD),med(R[P==p]$dT_LAD)))
