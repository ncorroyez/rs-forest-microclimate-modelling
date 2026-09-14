# ==============================================================================
# A1 — plot-level bootstrap CIs on the key summer R² claims (windcorr genuine-53).
# Resamples the 53 plots with replacement (2000 reps), recomputes each R², reports
# point estimate + percentile CI95. Addresses the "no uncertainty on R²" reviewer
# blocker. Central tests: (i) stratified reversal robust? (ii) fusion R² advantage
# over S2-alone significant or marginal? (iii) RF gain (S2 adds over height) real?
#   Rscript c3_bootstrap_ci.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(randomForest);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
WC<-file.path(CFG_C3$out_dir,"nc_genuine53_windcorr"); FORC<-"in_files/musica_in_Blois_pblh.nc"; CROSS<-3.86
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds"))); df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")
getpp<-function(scn){fs<-list.files(file.path(WC,scn),pattern="\\.nc$",full.names=TRUE)
 rbindlist(lapply(fs,function(f){id<-str_extract(basename(f),"(?<=HOBO_).*(?=\\.nc)");nc<-tryCatch(nc_open(f),error=function(e)NULL);if(is.null(nc))return(NULL)
  r<-tryCatch(get_tair_at_z(nc,1.0),error=function(e)NULL);nc_close(nc);if(is.null(r)||!nrow(r))return(NULL)
  s<-as.data.table(r)[as.Date(time)%in%ds,.(Tsim=Tair_sim,date=as.Date(time),hr=floor_date(time,"hour"))];sm<-merge(s,MA[,.(hr=time,Tm)],by="hr")
  sd_<-merge(s[,.(Tmx_s=max(Tsim,na.rm=TRUE)),by=date],mday,by="date");data.table(id_plot=id,ds=mean(sd_$Tmx_s-sd_$Tmx,na.rm=TRUE),ss=slp(sm$Tsim,sm$Tm))}),fill=TRUE)}
grab<-function(s,suf){P<-getpp(s);setnames(P,c("ds","ss"),paste0(c("d","s"),suf));P}
# dynsum for S2opt/S2atbd (for RF)
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"; pts<-terra::vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
dynsum<-function(pat){files<-list.files(NM,pattern=pat,full.names=TRUE);dts<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
 vals<-as.data.frame(terra::extract(terra::rast(files),pts))[,-1,drop=FALSE];long<-rbindlist(lapply(seq_along(dts),function(i)data.table(plot_id=df$pid,doy=as.integer(format(dts[i],"%j")),lai=pmax(vals[[i]],0))))
 ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3); sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
M<-Reduce(function(a,b)merge(a,b,by="id_plot"),list(OBS,grab("STATIC_ALS","L"),grab("STATIC_S2_ATBD","Sa")))
M<-merge(M,df[,.(id_plot,LAI_ALS,FORMS_H,pid)],by="id_plot")
M$S2opt<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$")[M$pid]
M$S2atbd<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$")[M$pid]
M$hfrac<-pmax(M$FORMS_H,0.1)/CFG_C3$d_opt_m
r2<-function(a,b){ok<-is.finite(a)&is.finite(b); if(sum(ok)<4) return(NA); cor(a[ok],b[ok])^2}
ci<-function(v) sprintf("%.3f [%.3f, %.3f]", mean(v,na.rm=T), quantile(v,.025,na.rm=T), quantile(v,.975,na.rm=T))
B<-2000; n<-nrow(M); set.seed(1)
# preallocate
st<-list(dL=numeric(B),dS=numeric(B),oL=numeric(B),oS=numeric(B),pL=numeric(B),pS=numeric(B),
         fus=numeric(B),fus_minus_S2=numeric(B),rf_h=numeric(B),rf_hs=numeric(B),rf_gain=numeric(B))
looRF<-function(D,fs){f<-as.formula(paste("LAI_ALS~",paste(fs,collapse="+")));set.seed(42)
  sapply(1:nrow(D),function(i)predict(randomForest(f,D[-i,],ntree=400),D[i,]))}
for(b in 1:B){ idx<-sample(n,n,replace=TRUE); Mb<-M[idx]; de<-Mb$LAI_ALS>=CROSS
  st$dL[b]<-r2(Mb$dL[de],Mb$do[de]); st$dS[b]<-r2(Mb$dSa[de],Mb$do[de])
  st$oL[b]<-r2(Mb$dL[!de],Mb$do[!de]); st$oS[b]<-r2(Mb$dSa[!de],Mb$do[!de])
  st$pL[b]<-r2(Mb$dL,Mb$do); st$pS[b]<-r2(Mb$dSa,Mb$do)
  fus<-ifelse(de,Mb$dL,Mb$dSa); st$fus[b]<-r2(fus,Mb$do); st$fus_minus_S2[b]<-r2(fus,Mb$do)-r2(Mb$dSa,Mb$do)
}
# RF bootstrap (cheaper: 300 reps, LOO inside)
set.seed(2); BR<-300; rfh<-numeric(BR); rfhs<-numeric(BR); rfg<-numeric(BR)
for(b in 1:BR){ idx<-sample(n,n,replace=TRUE); Db<-as.data.frame(M[idx,.(LAI_ALS,S2opt,S2atbd,FORMS_H,hfrac,do)])
  ph<-looRF(Db,c("FORMS_H","hfrac")); phs<-looRF(Db,c("S2opt","S2atbd","FORMS_H","hfrac"))
  rfh[b]<-r2(ph,Db$LAI_ALS); rfhs[b]<-r2(phs,Db$LAI_ALS); rfg[b]<-r2(phs,Db$LAI_ALS)-r2(ph,Db$LAI_ALS) }
cat("=== BOOTSTRAP CI95 (plot-level, ", B, "reps; RF ",BR," reps) ===\n",sep="")
cat("STRATIFIED ΔTmax R²:\n")
cat(sprintf("  dense  LiDAR %s | S2 %s\n", ci(st$dL), ci(st$dS)))
cat(sprintf("  open   LiDAR %s | S2 %s\n", ci(st$oL), ci(st$oS)))
cat(sprintf("  pooled LiDAR %s | S2 %s\n", ci(st$pL), ci(st$pS)))
cat("FUSION ΔTmax R²:\n")
cat(sprintf("  fusion %s\n", ci(st$fus)))
cat(sprintf("  fusion - S2alone (ΔR²) %s  [CI incl. 0 => marginal]\n", ci(st$fus_minus_S2)))
cat("RF LAI recovery:\n")
cat(sprintf("  FORMS-H only %s | S2+FORMS-H %s\n", ci(rfh), ci(rfhs)))
cat(sprintf("  RF gain (S2 adds) ΔR² %s\n", ci(rfg)))
saveRDS(st,"/tmp/claude-1001/-home-corroyez-Documents-NC-Full/d2bf1c39-bdf8-4dba-b5af-e8a63e72bc2f/scratchpad/boot.rds")
cat("\nDONE\n")
