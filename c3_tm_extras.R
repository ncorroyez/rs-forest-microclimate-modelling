# ==============================================================================
# Chapter 3 harmonization — extras:
#  (1) figA LAI-space correlations (S2 dynamic products vs observed dTmax) old/new
#  (2) wind-correction robustness under the new convention: bias shift STATIC_ALS
#      (windcorr vs non-windcorr) + Spearman rho of scenario bias/SDrec rankings
#   Rscript c3_tm_extras.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr)
  library(terra);library(data.table);library(parallel);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
FORC<-"in_files/musica_in_Blois_pblh.nc"
PP<-fread(file.path(CFG_C3$out_dir,"tables","perplot_dtmax_conventions.csv"))
OBS<-PP[scenario=="OBS",.(id_plot,do_old=d_old,do_new=d_new)]
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))

# ---- (1) figA LAI-space -------------------------------------------------------
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
dynsum<-function(pat){files<-list.files(NM,pattern=pat,full.names=TRUE)
 dts<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
 vals<-as.data.frame(terra::extract(rast(files),pts))[,-1,drop=FALSE]
 long<-rbindlist(lapply(seq_along(dts),function(i)data.table(plot_id=df$pid,
   doy=as.integer(format(dts[i],"%j")),lai=pmax(vals[[i]],0))))
 ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3)
 sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
df$S2opt<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$")
df$S2atbd<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$")
Dl<-merge(OBS,df[,.(id_plot,LAI_ALS,S2opt,S2atbd,FORMS_H)],by="id_plot")
A<-rbindlist(lapply(c("LAI_ALS","S2opt","S2atbd","FORMS_H"),function(v)
  data.table(product=v,cor_dTmax_old=cor(Dl[[v]],Dl$do_old),cor_dTmax_new=cor(Dl[[v]],Dl$do_new))))
fwrite(A,file.path(TAB,"timematched_laispace_cor.csv"))
cat("=== figA LAI-space correlations ===\n");print(A,digits=3)

# ---- (2) wind-correction robustness -------------------------------------------
MREF<-macro_ref(FORC,ds)
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds]
mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
NWC<-file.path(CFG_C3$out_dir,"nc_genuine53")
one<-function(scn,base){fs<-list.files(file.path(base,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)")
  nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  tr<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc)
  if(is.null(tr)||!nrow(tr))return(NULL)
  m<-as.data.table(tr)[,.(time,Tmic=Tair_sim)]
  mm<-m[,.(Tmic=mean(Tmic,na.rm=TRUE)),by=.(time=floor_date(time,"hour"))][as.Date(time)%in%ds]
  dd<-merge(mm[,.(Tmx_s=max(Tmic)),by=.(date=as.Date(time))],mday,by="date")
  data.table(id_plot=id,d_old=mean(dd$Tmx_s-dd$Tmx,na.rm=TRUE),
             d_new=delta_tmax_mean(mm,MREF,ds))}),fill=TRUE)}
SCN<-list.dirs(NWC,recursive=FALSE,full.names=FALSE)
SCN<-SCN[sapply(SCN,function(s)length(list.files(file.path(NWC,s),pattern="\\.nc$"))==53)]
cat(sprintf("\nnon-windcorr scenarios: %d\n",length(SCN)))
NW<-rbindlist(mclapply(SCN,function(s){r<-one(s,NWC);r[,scenario:=s];r},mc.cores=4),fill=TRUE)
sumr<-function(X,ov){m<-merge(X,ov,by="id_plot")
  m[,.(bias_old=mean(d_old.x-d_old.y),bias_new=mean(d_new.x-d_new.y),
       SDrec_old=sd(d_old.x)/sd(d_old.y),SDrec_new=sd(d_new.x)/sd(d_new.y)),by=scenario]}
OV<-PP[scenario=="OBS",.(id_plot,d_old,d_new)]
nw<-sumr(NW,OV)
wc<-sumr(PP[scenario%in%SCN],OV)
CMP<-merge(nw,wc,by="scenario",suffixes=c("_nwc","_wc"))
fwrite(CMP,file.path(TAB,"timematched_windcorr_robustness.csv"))
cat("\n=== wind-corr robustness (scenario-level) ===\n")
cat(sprintf("STATIC_ALS bias shift (wc - nwc): old %+.3f | new %+.3f\n",
 CMP[scenario=="STATIC_ALS",bias_old_wc-bias_old_nwc],CMP[scenario=="STATIC_ALS",bias_new_wc-bias_new_nwc]))
cat(sprintf("Spearman rho bias  across %d scenarios: old %.3f | new %.3f\n",nrow(CMP),
 cor(CMP$bias_old_nwc,CMP$bias_old_wc,method="spearman"),cor(CMP$bias_new_nwc,CMP$bias_new_wc,method="spearman")))
cat(sprintf("Spearman rho SDrec across %d scenarios: old %.3f | new %.3f\n",nrow(CMP),
 cor(CMP$SDrec_old_nwc,CMP$SDrec_old_wc,method="spearman"),cor(CMP$SDrec_new_nwc,CMP$SDrec_new_wc,method="spearman")))
cat("wrote timematched_laispace_cor.csv + timematched_windcorr_robustness.csv\n")
