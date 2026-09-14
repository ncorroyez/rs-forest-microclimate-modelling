# ==============================================================================
# Raw-bands random forest, scored on archetype P4 (replaces the retired LAI 3.86
# stratum). Same recipe as c3_rf_rawbands.R / c3_rf_bandset_test.R: leave-one-
# plot-out over the 53 plots, ntree = 500, seed 42, target = the observed summer
# daytime offset over 1 June - 30 September. References are the between-plot R2
# of that offset against each leaf-area product, computed on the same plots.
# Out: chapter3_S2_LAI/tables/TableB_rf_P4.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(terra); library(ncdf4); library(lubridate)
  library(data.table); library(randomForest)
})
source("Chapter3_config.R")
TAB  <- "/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
REFL <- paste0("/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/",
               "L2A_T31TCN_A031222_20210614T105443/Reflectance/res_10_m/",
               "L2A_T31TCN_A031222_20210614T105443_Refl")
BND  <- c("B02","B03","B04","B05","B06","B07","B08","B8A","B11","B12")
BND3 <- c("B03","B04","B08")
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")

df <- as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
r  <- rast(REFL); names(r) <- BND
pts<- vect(as.data.frame(df[,.(x,y)]), geom=c("x","y"), crs="EPSG:32631")
df <- cbind(df, as.data.table(terra::extract(r, pts))[, -1])

ncf<-nc_open("in_files/musica_in_Blois_pblh.nc")
tu<-ncatt_get(ncf,"time","units")$value; th<-ncvar_get(ncf,"time")
t0<-as.POSIXct(sub(".*since ","",tu),tz="UTC")
MA<-data.table(time=floor_date(t0+th*3600,"hour"),
               Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15); nc_close(ncf)
MA<-MA[as.Date(time)%in%ds]
mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv))
hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a" & !(id_plot%in%CFG_C3$ids_to_remove) & as.Date(time)%in%ds,
       .(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time))]
OBS<-merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[
  ,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot]
cl<-fread("out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv")[,.(id_plot=as.character(id_plot),P)]
M<-merge(merge(OBS,df,by="id_plot"),cl,by="id_plot"); stopifnot(nrow(M)==53)
p4 <- M$P=="P4"; cat("P4 n =", sum(p4), "\n")

r2<-function(a,b){ok<-is.finite(a)&is.finite(b); if(sum(ok)<4) NA_real_ else cor(a[ok],b[ok])^2}
loo_rf<-function(feat,mtry=2L,seed=42L){
  X<-as.data.frame(M[,..feat]); y<-M$do; set.seed(seed); p<-rep(NA_real_,nrow(X))
  for(i in seq_len(nrow(X))){ rf<-randomForest(x=X[-i,,drop=FALSE],y=y[-i],ntree=500,mtry=mtry)
    p[i]<-as.numeric(predict(rf,X[i,,drop=FALSE])) }; p }
boot<-function(o,p,B=2000L,seed=7L){set.seed(seed)
  as.numeric(quantile(replicate(B,{i<-sample(length(o),replace=TRUE);r2(o[i],p[i])}),
                      c(.025,.975),na.rm=TRUE))}

p10<-loo_rf(BND); p3<-loo_rf(BND3)
res<-rbindlist(list(
  data.table(predictor="Random forest, ten raw bands", type="LOO forest",
             R2=r2(M$do[p4],p10[p4]), lo=boot(M$do[p4],p10[p4])[1], hi=boot(M$do[p4],p10[p4])[2]),
  data.table(predictor="Random forest, the three inversion bands", type="LOO forest",
             R2=r2(M$do[p4],p3[p4]), lo=boot(M$do[p4],p3[p4])[1], hi=boot(M$do[p4],p3[p4])[2]),
  data.table(predictor="Retrieved Sentinel-2 LAI", type="leaf area",
             R2=r2(M$do[p4],M$LAI_S2_ATBD[p4]),
             lo=boot(M$do[p4],M$LAI_S2_ATBD[p4])[1], hi=boot(M$do[p4],M$LAI_S2_ATBD[p4])[2]),
  data.table(predictor="LiDAR LAI, full column", type="leaf area",
             R2=r2(M$do[p4],M$LAI_ALS[p4]),
             lo=boot(M$do[p4],M$LAI_ALS[p4])[1], hi=boot(M$do[p4],M$LAI_ALS[p4])[2]),
  data.table(predictor="LiDAR LAI, top d_opt layer", type="leaf area",
             R2=r2(M$do[p4],M$LAI_ALS_DOPT[p4]),
             lo=boot(M$do[p4],M$LAI_ALS_DOPT[p4])[1], hi=boot(M$do[p4],M$LAI_ALS_DOPT[p4])[2])))
res[,n:=sum(p4)]
print(res[,.(predictor,R2=round(R2,3),CI=sprintf("[%.2f, %.2f]",lo,hi),n)])
fwrite(res, file.path(TAB,"TableB_rf_P4.csv"))
cat("wrote TableB_rf_P4.csv\n")
