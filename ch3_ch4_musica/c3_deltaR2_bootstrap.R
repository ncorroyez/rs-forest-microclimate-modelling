# ==============================================================================
# Paired ΔR² bootstrap on the DENSE stratum (LAI_ALS >= 3.86, n=25), target =
# observed summer ΔTmax `do` (windcorr genuine-53). Two load-bearing paired
# differences, previously supported only by non-overlapping CIs:
#  (A) H1 dense : R²(LiDAR STATIC_ALS) − R²(S2 STATIC_S2_ATBD)   (~0.60 vs ~0.02)
#  (B) H2 nuance: R²(raw-bands RF LOO) − R²(S2 STATIC_S2_ATBD)   (~0.35 vs ~0.02)
# Protocol: `do` and dense stratum verbatim from c3_split_explore_wc.R; plot-level
# bootstrap as in c3_bootstrap_ci.R — resample the 25 dense plots with replacement,
# B=2000, percentile CI95, seed 42. Raw-bands RF preds are the LOO (53-fold)
# predictions saved by c3_rf_rawbands.R (figH2_rawbands_preds.csv); dense
# evaluation conditions on those fixed LOO predictions (no per-replicate refit).
#   Rscript c3_deltaR2_bootstrap.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"; ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"; CROSS<-3.86

## ---- target `do` (verbatim block from c3_split_explore_wc.R lines 16-24) -----
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")

## ---- simulated ΔTmax per plot (verbatim getpp from c3_split_explore_wc.R) ----
getpp<-function(scn){fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
L<-getpp("STATIC_ALS");setnames(L,c("ds","ss"),c("dL","sL")); S<-getpp("STATIC_S2_ATBD");setnames(S,c("ds","ss"),c("dS","sS"))
M<-merge(merge(merge(OBS,L,by="id_plot"),S,by="id_plot"),df[,.(id_plot,LAI_ALS)],by="id_plot")

## ---- raw-bands RF LOO predictions (from c3_rf_rawbands.R, primary variant) ---
PR<-as.data.table(fread(file.path(TAB,"figH2_rawbands_preds.csv")))[,.(id_plot=as.character(id_plot),pred)]
M<-merge(M,PR,by="id_plot")
stopifnot(nrow(M)==53)
r2<-function(a,b){ok<-is.finite(a)&is.finite(b); if(sum(ok)<4) return(NA_real_); cor(a[ok],b[ok])^2}

## ---- dense stratum + paired bootstrap ----------------------------------------
D<-M[LAI_ALS>=CROSS]; nD<-nrow(D)
cat(sprintf("dense stratum: n=%d (LAI_ALS >= %.2f)\n",nD,CROSS))
cat(sprintf("point R2 dense: LiDAR=%.3f  S2_ATBD=%.3f  rawbandsRF=%.3f\n",
            r2(D$dL,D$do),r2(D$dS,D$do),r2(D$pred,D$do)))
dA_pt<-r2(D$dL,D$do)-r2(D$dS,D$do)
dB_pt<-r2(D$pred,D$do)-r2(D$dS,D$do)
B<-2000; set.seed(42); dA<-numeric(B); dB<-numeric(B)
for(b in 1:B){ i<-sample(nD,nD,replace=TRUE)
  dA[b]<-r2(D$dL[i],D$do[i])   - r2(D$dS[i],D$do[i])
  dB[b]<-r2(D$pred[i],D$do[i]) - r2(D$dS[i],D$do[i]) }
qs<-function(v) as.numeric(quantile(v,c(.025,.975),na.rm=TRUE))
qA<-qs(dA); qB<-qs(dB)
out<-data.table(
  comparison=c("A_H1dense_LiDAR_minus_S2ATBD","B_H2_rawbandsRF_minus_S2ATBD"),
  deltaR2=c(dA_pt,dB_pt), ci_lo=c(qA[1],qB[1]), ci_hi=c(qA[2],qB[2]), n=nD)
cat("\n=== paired ΔR² bootstrap (dense, B=2000, percentile CI95, seed 42) ===\n")
print(out,digits=3)
cat(sprintf("(A) excludes zero: %s | (B) excludes zero: %s\n",
            ifelse(qA[1]>0|qA[2]<0,"YES","NO"),ifelse(qB[1]>0|qB[2]<0,"YES","NO")))
cat(sprintf("NA replicates: A=%d B=%d\n",sum(!is.finite(dA)),sum(!is.finite(dB))))
fwrite(out,file.path(TAB,"figH2_deltaR2_bootstrap.csv"))
cat(sprintf("wrote %s\n",file.path(TAB,"figH2_deltaR2_bootstrap.csv")))
