# ==============================================================================
# c4_gedi_direct.R — Chapter 4: DIRECT use of GEDI footprints at the HOBO plots
# ("use the LiDAR where it measures", no regression).
#
# For each of the 53 HOBO plots (Blois), take the nearest power-beam GEDI
# footprint(s) and use their PAI directly as the plot LAI. Quantify skill as a
# function of the plot-to-footprint distance, in two currencies:
#   (a) LAI space  : cor(GEDI_direct, plot LAI_ALS)
#   (b) microclimate: cor(GEDI_direct, observed summer ΔTmax and slope)
# Control that separates GEDI noise from forest spatial decorrelation: the ALS
# LAI extracted AT the footprint location (matchup LAI_lidar_cor, k-rescaled) —
# if it degrades with distance like GEDI does, the distance penalty is the
# forest, not the sensor.
# Writes NC_Full/manuscripts/ch4/tables/Table10_c4_direct_distance.csv.
#   Rscript c4_gedi_direct.R
# ==============================================================================
suppressPackageStartupMessages({library(ncdf4);library(lubridate);library(stringr);library(data.table);library(dplyr);library(rmusica);library(musica.tools)})
src<-list.files("R",pattern="\\.R$",full.names=TRUE);src<-src[!grepl("/(h1_|lovb_)",src)];invisible(lapply(src,source));source("Chapter3_config.R")
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch4/tables"
ds<-seq(as.Date("2021-06-01"),as.Date("2021-09-30"),by="day")
k_scale<-0.5/0.65

# --- observed microclimate (same prelude as c4_spacelai_gedi.R) ---------------
FMACRO<-"in_files/musica_in_Blois_pblh.nc"
ncf<-nc_open(FMACRO);MA<-data.table(time=force_utc_nc(FMACRO,"time"),Tm=as.numeric(ncvar_get(ncf,"Tair"))-273.15);nc_close(ncf)
MA<-MA[as.Date(time)%in%ds];MA[,hr:=floor_date(time,"hour")];mday<-MA[,.(Tmx=max(Tm,na.rm=TRUE)),by=.(date=as.Date(time))]
slp<-function(y,x) if(sum(is.finite(x)&is.finite(y))>=10) as.numeric(coef(lm(y~x))[2]) else NA
hb<-as.data.table(read.csv(CFG_C3$hobo_temp_csv));hb[,time:=as.POSIXct(datetime,format="%Y-%m-%d %H:%M:%S",tz="UTC")]
hb<-hb[position_sensor=="a"&!(id_plot%in%CFG_C3$ids_to_remove)&as.Date(time)%in%ds,.(id_plot=as.character(id_plot),Tobs=t_hobo,date=as.Date(time),hr=floor_date(time,"hour"))]
hb<-merge(hb,MA[,.(hr,Tm)],by="hr")
OBS<-merge(merge(hb[,.(Tmx_o=max(Tobs,na.rm=TRUE)),by=.(id_plot,date)],mday,by="date")[,.(do=mean(Tmx_o-Tmx,na.rm=TRUE)),by=id_plot],hb[,.(so=slp(Tobs,Tm)),by=id_plot],by="id_plot")

# --- plots x footprints -------------------------------------------------------
fp<-readRDS("/home/corroyez/Documents/NC_Full/output/intermediate/c4/gedi_matchup_3site.rds")
fp<-fp[site=="Blois"&power=="full"]
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
df[,id_plot:=as.character(id_plot)]
D<-sqrt(outer(df$x,fp$x,"-")^2+outer(df$y,fp$y,"-")^2)
nn<-apply(D,1,which.min)
df[,`:=`(dist_fp=D[cbind(.I,nn)],
         GEDI_direct=fp$pai_gedi[nn],
         ALS_at_fp=fp$LAI_lidar_cor[nn]*k_scale,          # distance control
         GEDI_k3=sapply(.I,function(i){o<-order(D[i,])[1:3];mean(fp$pai_gedi[o])}))]
d<-merge(df,OBS,by="id_plot")
t7<-fread(file.path(TAB,"Table7_c4_perplot_gedi_lai.csv"))[,.(id_plot=as.character(id_plot),LAI_GEDI_RF)]
d<-merge(d,t7,by="id_plot")

# --- skill vs distance (cumulative radius) ------------------------------------
prods<-c(GEDI_direct="GEDI_direct",GEDI_k3="GEDI_k3",ALS_at_fp="ALS_at_fp",
         GEDI_RF_product="LAI_GEDI_RF",S2_raw="LAI_S2_ATBD",ALS_plot="LAI_ALS")
res<-rbindlist(lapply(c(100,150,200,300,500,1300),function(r){
  s<-d[dist_fp<=r]
  rbindlist(lapply(names(prods),function(p){x<-s[[prods[[p]]]]
    data.table(radius_m=r,n=nrow(s),product=p,
      r_LAI_ALS=cor(x,s$LAI_ALS),
      r_dTmax=cor(x,s$do), r_slope=cor(x,s$so,use="complete.obs"))}))}))
num<-names(res)[sapply(res,is.numeric)];res[,(num):=lapply(.SD,round,3),.SDcols=num]
print(dcast(res,radius_m+n~product,value.var="r_LAI_ALS"),digits=3)
cat("\n--- r with observed ΔTmax ---\n")
print(dcast(res,radius_m+n~product,value.var="r_dTmax"),digits=3)
cat("\n--- r with observed slope ---\n")
print(dcast(res,radius_m+n~product,value.var="r_slope"),digits=3)
fwrite(res,file.path(TAB,"Table10_c4_direct_distance.csv"))
saveRDS(d,file.path(CFG_C3$out_dir,"lai_prep","df_plots_gedi_direct.rds"))
cat("\nwrote Table10 +",file.path(CFG_C3$out_dir,"lai_prep","df_plots_gedi_direct.rds"),"\n")
