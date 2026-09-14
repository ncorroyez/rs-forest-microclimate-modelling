# ==============================================================================
# Chapter 3 harmonization — RF-based numbers under the time-matched convention:
#  (a) raw-bands LOO RF (H2 shield): dense/pooled R2, dense CI, old vs new target
#  (b) paired ΔR2 bootstrap B (dense): rawbandsRF − S2_ATBD, old vs new
#  (c) figG ablation: LOO_R2_LAI (expected UNCHANGED) + cor_dTmax old vs new
# Reads perplot_dtmax_conventions.csv (from c3_tm_extract.R).
#   Rscript c3_tm_rf.R
# ==============================================================================
suppressPackageStartupMessages({library(terra);library(data.table);library(stringr)
  library(randomForest)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)]
invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
CROSS<-3.86; dopt<-CFG_C3$d_opt_m
PP<-fread(file.path(CFG_C3$out_dir,"tables","perplot_dtmax_conventions.csv"))
OBS<-PP[scenario=="OBS",.(id_plot,do_old=d_old,do_new=d_new)]
S2A<-PP[scenario=="STATIC_S2_ATBD",.(id_plot,dS_old=d_old,dS_new=d_new)]
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df$pid<-sprintf("X%d_Y%d",round(df$x),round(df$y))
r2<-function(a,b){ok<-is.finite(a)&is.finite(b);if(sum(ok)<4)return(NA_real_);cor(a[ok],b[ok])^2}

# ---- (a) raw bands ------------------------------------------------------------
REFL<-paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/",
  "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/",
  "L2A_T31TCN_A031222_20210614T105443_Refl")
r<-rast(REFL); pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
BND<-c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12"); names(r)<-BND
ex<-as.data.table(terra::extract(r,pts))[,-1]; df<-cbind(df,ex)
M<-merge(OBS,df,by="id_plot"); M<-merge(M,S2A,by="id_plot")
stopifnot(nrow(M)==53)
de<-M$LAI_ALS>=CROSS
loo_rf<-function(feat,y,mtry=2L,seed=42){X<-as.data.frame(M[,..feat]);set.seed(seed)
  p<-rep(NA_real_,nrow(X))
  for(i in seq_len(nrow(X))){rf<-randomForest(x=X[-i,,drop=FALSE],y=y[-i],ntree=500,mtry=mtry)
    p[i]<-as.numeric(predict(rf,X[i,,drop=FALSE]))};p}
boot_dense<-function(o,p,B=2000,seed=7){o<-o[de];p<-p[de];set.seed(seed)
  as.numeric(quantile(replicate(B,{i<-sample(length(o),replace=TRUE);r2(o[i],p[i])}),c(.025,.975),na.rm=TRUE))}
res<-list()
for(cv in c("old","new")){
  y<-if(cv=="old")M$do_old else M$do_new
  p<-loo_rf(BND,y)
  ci<-boot_dense(y,p)
  res[[cv]]<-data.table(conv=cv,pooled_r2=r2(y,p),dense_r2=r2(y[de],p[de]),
    open_r2=r2(y[!de],p[!de]),dense_ci_lo=ci[1],dense_ci_hi=ci[2])
  # ---- (b) paired ΔR2 bootstrap B: rawbandsRF − S2_ATBD on dense --------------
  s<-if(cv=="old")M$dS_old else M$dS_new
  Dp<-p[de];Do<-y[de];Ds<-s[de];nD<-sum(de)
  set.seed(42);pt<-r2(Dp,Do)-r2(Ds,Do)
  v<-replicate(2000,{i<-sample(nD,nD,replace=TRUE);r2(Dp[i],Do[i])-r2(Ds[i],Do[i])})
  res[[paste0("B_",cv)]]<-data.table(conv=cv,B_deltaR2=pt,
    B_ci_lo=quantile(v,.025,na.rm=TRUE),B_ci_hi=quantile(v,.975,na.rm=TRUE))
  cat(sprintf("[%s] rawbands: pooled=%.3f dense=%.3f [%.3f,%.3f] open=%.3f | B dR2=%.3f [%.3f,%.3f]\n",
    cv,r2(y,p),r2(y[de],p[de]),ci[1],ci[2],r2(y[!de],p[!de]),
    pt,quantile(v,.025,na.rm=TRUE),quantile(v,.975,na.rm=TRUE)))
}
RB<-rbind(res$old,res$new); BB<-rbind(res$B_old,res$B_new)
fwrite(cbind(RB,BB[,-1]),file.path(TAB,"timematched_rawbands_rf.csv"))

# ---- (c) ablation figG: LAI target (unchanged) + cor with dTmax old/new ---------
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
dynsum<-function(pat){files<-list.files(NM,pattern=pat,full.names=TRUE)
 dts<-as.Date(str_extract(basename(files),"\\d{4}-\\d{2}-\\d{2}"))
 vals<-as.data.frame(terra::extract(rast(files),pts))[,-1,drop=FALSE]
 long<-rbindlist(lapply(seq_along(dts),function(i)data.table(plot_id=df$pid,
   doy=as.integer(format(dts[i],"%j")),lai=pmax(vals[[i]],0))))
 ts<-smooth_s2_ts(as.data.frame(long),k=8,min_obs=3)
 sapply(df$pid,function(p){a<-ts[[p]];mean(a$lai[a$doy>=152&a$doy<=244],na.rm=TRUE)})}
df$S2opt<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_optim_common_res_10_m\\.tif$")
df$S2atbd<-dynsum("^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$")
df$hfrac<-pmax(df$FORMS_H,0.1)/dopt
Dl<-as.data.frame(merge(OBS,df[,.(id_plot,LAI_ALS,S2opt,S2atbd,FORMS_H,hfrac,Hmax,fCover,VCI,LCV)],by="id_plot"))
looRF<-function(fs){set.seed(42);f<-as.formula(paste("LAI_ALS~",paste(fs,collapse="+")))
  sapply(1:nrow(Dl),function(i)predict(randomForest(f,Dl[-i,],ntree=400),Dl[i,]))}
sets<-list("S2 only"=c("S2opt","S2atbd"),"FORMS-H only"=c("FORMS_H","hfrac"),
  "S2 + FORMS-H"=c("S2opt","S2atbd","FORMS_H","hfrac"),
  "S2 + LiDAR structure"=c("S2opt","S2atbd","Hmax","fCover","VCI","LCV"),
  "LiDAR structure"=c("Hmax","fCover","VCI","LCV"))
G<-rbindlist(lapply(names(sets),function(nm){p<-looRF(sets[[nm]])
  data.table(feature_set=nm,LOO_R2_LAI=cor(p,Dl$LAI_ALS)^2,
    cor_dTmax_old=cor(p,Dl$do_old),cor_dTmax_new=cor(p,Dl$do_new))}))
G<-rbind(G,data.table(feature_set="LiDAR LAI (direct)",LOO_R2_LAI=1,
  cor_dTmax_old=cor(Dl$LAI_ALS,Dl$do_old),cor_dTmax_new=cor(Dl$LAI_ALS,Dl$do_new)))
fwrite(G,file.path(TAB,"timematched_ablation.csv"))
cat("\n=== ablation (LOO_R2_LAI must equal figG_ablation.csv) ===\n");print(G,digits=3)
cat("wrote timematched_rawbands_rf.csv + timematched_ablation.csv\n")
