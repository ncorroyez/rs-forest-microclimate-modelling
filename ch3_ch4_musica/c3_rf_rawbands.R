# ==============================================================================
# H2 shield test — RAW Sentinel-2 bands vs observed summer buffering (ΔTmax).
#
# H2 claims the dense-canopy failure of S2 is a property of the S2 SIGNAL, not
# of the PROSAIL inversion. All S2 products tested so far derive from the same
# inverted LAI. Here a flexible learner (Random Forest, LOO-CV, n=53) gets the
# 10 RAW L2A bands (B02..B12, incl. red-edge/NIR/SWIR) of the SAME summer scene
# (2021-06-14 Blois) that produced the chapter's LAI_S2, and must predict the
# observed microclimate target `do` = mean(Tmax_obs - Tmax_macro) over JJAS —
# identical to Fig. 4 / figJ_stratified_cross386.csv.
#  - PRIMARY: 10 raw bands only, LOO 53 folds, RF ntree=500 mtry=2, seed 42.
#  - Robustness: + NDVI/NDRE1/NDRE2/NBR + 3x3 focal-sd texture (24 features).
#  - Metrics: pooled R2 (cor^2) + within-stratum R2 at the LAI_ALS >= 3.86
#    split (25 dense / 28 open), bootstrap CI95 on dense R2 (pair resampling).
# If dense R2 stays ~<=0.2 even for raw bands -> H2 shielded.
#   Rscript c3_rf_rawbands.R
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate); library(stringr)
  library(data.table); library(randomForest)
})
source("Chapter3_config.R")

TAB   <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
REFL  <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/",
                "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/",
                "L2A_T31TCN_A031222_20210614T105443_Refl")
FORC  <- "in_files/musica_in_Blois_pblh.nc"
CROSS <- 3.86
ds    <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")

# ---- plots + raw-band extraction (READ-ONLY raster access) -------------------
df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
r  <- rast(REFL)                                   # 10 bands B02..B12, EPSG:32631
pts <- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs="EPSG:32631")
BND <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
names(r) <- BND
ex <- as.data.table(terra::extract(r, pts))[, -1]
df <- cbind(df, ex)
cat(sprintf("extracted %d bands at %d plots; NA rows: %d\n",
            ncol(ex), nrow(df), sum(!complete.cases(ex))))

# ---- spectral indices + 3x3 focal-sd texture (robustness variant) ------------
df[, `:=`(NDVI  = (B08-B04)/(B08+B04), NDRE1 = (B8A-B05)/(B8A+B05),
          NDRE2 = (B8A-B06)/(B8A+B06), NBR   = (B08-B11)/(B08+B11))]
rt <- focal(r, w=3, fun="sd", na.rm=TRUE)          # per-band local texture
names(rt) <- paste0("sd_", BND)
ext <- as.data.table(terra::extract(rt, pts))[, -1]
df  <- cbind(df, ext)
IDX <- c("NDVI","NDRE1","NDRE2","NBR"); TEX <- names(rt)

# ---- target `do` = mean(Tmax_obs - Tmax_macro) JJAS  (verbatim block from ----
# ---- c3_split_explore_wc.R lines 17-24: identical to Fig. 4 target) ----------
ncf<-nc_open(FORC);tu<-ncatt_get(ncf,"time","units")$value;th<-ncvar_get(ncf,"time");t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr=time,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")

M <- merge(OBS, df, by="id_plot")
cat(sprintf("assembled n=%d plots (dense n=%d / open n=%d at LAI_ALS>=%.2f)\n",
            nrow(M), sum(M$LAI_ALS>=CROSS), sum(M$LAI_ALS<CROSS), CROSS))

# ---- LOO-RF (calqued on loo_cv_rf(): ntree=500, mtry=2, explicit loop) -------
r2 <- function(a,b){ ok<-is.finite(a)&is.finite(b); if(sum(ok)<4) return(NA_real_)
  cor(a[ok],b[ok])^2 }
loo_rf <- function(feat, target, mtry=2L, seed=42){
  X <- as.data.frame(M[, ..feat]); y <- M[[target]]
  set.seed(seed); p <- rep(NA_real_, nrow(X))
  for (i in seq_len(nrow(X))) {
    rf <- randomForest(x=X[-i,,drop=FALSE], y=y[-i], ntree=500, mtry=mtry)
    p[i] <- as.numeric(predict(rf, X[i,,drop=FALSE]))
  }
  p
}
boot_dense_r2 <- function(obs, pred, dense, B=2000, seed=7){
  # pair-resampling bootstrap on the dense LOO (obs,pred) pairs; conditions on
  # the fitted LOO predictions (cheap; does not refit RF per replicate)
  o <- obs[dense]; p <- pred[dense]; set.seed(seed)
  qs <- quantile(replicate(B, { i<-sample(length(o),replace=TRUE); r2(o[i],p[i]) }),
                 c(.025,.975), na.rm=TRUE)
  as.numeric(qs)
}
eval_var <- function(label, feat, target="do", mtry=2L){
  p  <- loo_rf(feat, target, mtry=mtry)
  de <- M$LAI_ALS >= CROSS
  if (label == "bands10_raw") {  # save per-plot LOO preds for downstream ΔR² bootstrap
    fwrite(data.table(id_plot=M$id_plot, pred=p, obs=M[[target]], lai_als=M$LAI_ALS),
           file.path(TAB, "figH2_rawbands_preds.csv"))
  }
  ci <- boot_dense_r2(M[[target]], p, de)
  out <- data.table(features=label, target=target, n_feat=length(feat),
    pooled_r2=r2(M[[target]],p), dense_r2=r2(M[[target]][de],p[de]),
    open_r2=r2(M[[target]][!de],p[!de]), dense_ci_lo=ci[1], dense_ci_hi=ci[2],
    n_dense=sum(de), n_open=sum(!de))
  cat(sprintf("%-28s pooled=%.3f dense=%.3f [%.3f,%.3f] open=%.3f\n",
      label, out$pooled_r2, out$dense_r2, ci[1], ci[2], out$open_r2))
  out
}

cat("\n=== LOO-RF on raw S2 signal (target: observed dTmax `do`) ===\n")
M[, do_rank := rank(do)]
res <- rbind(
  eval_var("bands10_raw",            BND),                       # PRIMARY
  eval_var("bands10_raw_mtry3",      BND, mtry=3L),               # mtry robustness
  eval_var("bands10_raw_rank",       BND, target="do_rank"),      # rank target
  eval_var("bands24_idx_texture",    c(BND,IDX,TEX))              # overfit-prone
)
fwrite(res, file.path(TAB,"figH2_rawbands_loo.csv"))

cat("\nreference figJ_stratified_cross386.csv: S2(PROSAIL) dense R2=0.017, ")
cat("LiDAR dense R2=0.595\n")
cat(sprintf("wrote %s\n", file.path(TAB,"figH2_rawbands_loo.csv")))
