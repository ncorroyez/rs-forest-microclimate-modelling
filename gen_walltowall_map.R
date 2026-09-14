# ==============================================================================
# Chapter 3 — WALL-TO-WALL fused LAI map over Blois (deployability demonstration).
# Applies the SAME fusion rules validated at the 53 plots, per-pixel at 10 m:
#   * hard gate : LiDAR LAI where FORMS-H >= 17.8 m, else Sentinel-2 summer LAI
#   * blend     : w*LiDAR + (1-w)*S2,  w = 1/(1+exp(-(FORMS_H-17.8)/3))  (H0,s FIXED)
# Sources (10 m, EPSG:32631):
#   LiDAR LAI = 03_RESULTS/Blois/Metrics/Not_Masked/lidarlai_res_10_m.tif (r=0.95 vs
#     the chapter's plot LAI_ALS -> a TILE-WIDE product, NOT identical to plot LAI_ALS;
#     stated as a caveat, the dense-regime magnitudes differ slightly from validation)
#   FORMS-H gate = 01_DATA/FORMS-H_Blois.tif (r=1.000 vs df$FORMS_H -> exact)
#   S2 LAI = summer (DOY 152-244) per-pixel mean of the ATBD dated rasters
# Framing: the map demonstrates DEPLOYABILITY. It does NOT extend validation; every
# skill number stays attached to the 53 loggers. No MuSICA re-run; no dTmax map.
#   Rscript gen_walltowall_map.R
# ==============================================================================
suppressPackageStartupMessages({library(terra);library(data.table);library(stringr);library(ggplot2);library(tidyterra);library(patchwork)})
setwd("/home/corroyez/Documents/z_Example_rmusica_31012025"); source("Chapter3_config.R")
OUT<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures"; dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
B<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois"; H0<-17.8; S0<-3

lidar<-rast(file.path(B,"Metrics/Not_Masked/lidarlai_res_10_m.tif")); names(lidar)<-"lidar"
lidar<-clamp(lidar,0,8,values=FALSE)   # NA implausible LAI>8 (non-forest/edge artefacts)
formsh<-rast("/home/corroyez/Documents/NC_Full/01_DATA/FORMS-H_Blois.tif"); names(formsh)<-"formsh"
# S2 summer ATBD mean
s2f<-list.files(file.path(B,"Metrics/Not_Masked"),pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$",full.names=TRUE)
s2d<-as.Date(str_extract(basename(s2f),"\\d{4}-\\d{2}-\\d{2}")); doy<-as.integer(format(s2d,"%j"))
s2summer<-s2f[doy>=152 & doy<=244]; cat(sprintf("S2 summer dates: %d\n",length(s2summer)))
s2<-app(rast(s2summer),fun=function(x)mean(pmax(x,0),na.rm=TRUE)); names(s2)<-"s2"
# align to lidar grid
formsh<-resample(formsh,lidar,method="bilinear"); s2<-resample(s2,lidar,method="bilinear")
st<-c(lidar,formsh,s2)
# fusion rules
w<-1/(1+exp(-(formsh-H0)/S0)); names(w)<-"w"
sw<-ifel(formsh>=H0,lidar,s2); names(sw)<-"switch"
bl<-w*lidar+(1-w)*s2; names(bl)<-"blend"
# mask to valid canopy (both sources present)
m<-!is.na(lidar)&!is.na(s2)&!is.na(formsh); sw<-mask(sw,m,maskvalues=FALSE); bl<-mask(bl,m,maskvalues=FALSE)
writeRaster(sw,file.path(OUT,"walltowall_LAI_switch.tif"),overwrite=TRUE)
writeRaster(bl,file.path(OUT,"walltowall_LAI_blend.tif"),overwrite=TRUE)

# ---- plot-consistency: map values at 53 plots vs the plot-scale forcing --------
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
ex<-terra::extract(c(lidar,s2,formsh,sw,bl),pts)
df$m_lidar<-ex$lidar; df$m_s2<-ex$s2; df$m_sw<-ex$switch; df$m_bl<-ex$blend
r2<-function(a,b)cor(a,b,use="complete.obs")
cat("\n=== plot consistency (map vs plot forcing) ===\n")
cat(sprintf("map LiDAR   vs df$LAI_ALS      r=%.3f\n",r2(df$m_lidar,df$LAI_ALS)))
cat(sprintf("map S2      vs df$LAI_S2_ATBD  r=%.3f\n",r2(df$m_s2,df$LAI_S2_ATBD)))
cat(sprintf("map FORMS-H vs df$FORMS_H      r=%.3f\n",r2(df$m_formsh<-terra::extract(formsh,pts)[,2],df$FORMS_H)))
# what a plot-scale ideal switch/blend would be from df values
df$p_sw<-ifelse(df$FORMS_H>=H0,df$LAI_ALS,df$LAI_S2_ATBD)
wp<-1/(1+exp(-(df$FORMS_H-H0)/S0)); df$p_bl<-wp*df$LAI_ALS+(1-wp)*df$LAI_S2_ATBD
cat(sprintf("map switch  vs plot switch     r=%.3f  RMSE=%.2f\n",r2(df$m_sw,df$p_sw),sqrt(mean((df$m_sw-df$p_sw)^2,na.rm=TRUE))))
cat(sprintf("map blend   vs plot blend      r=%.3f  RMSE=%.2f\n",r2(df$m_bl,df$p_bl),sqrt(mean((df$m_bl-df$p_bl)^2,na.rm=TRUE))))

# ---- figure: 3 panels (LiDAR, S2, blend) + weight, one shared LAI scale --------
lims<-c(0,7.5)
mk<-function(r,ttl,pal="viridis",lim=lims,nm="LAI"){ggplot()+geom_spatraster(data=r)+
  scale_fill_viridis_c(option=pal,limits=lim,na.value="transparent",name=nm)+
  labs(subtitle=ttl)+coord_sf(expand=FALSE)+theme_minimal(base_size=10)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank())}
p1<-mk(lidar,"(a) LiDAR LAI (tile-wide)"); p2<-mk(s2,"(b) Sentinel-2 summer LAI")
p3<-mk(bl,"(c) Height-blended LAI (deployable)")
p4<-mk(w,"(d) LiDAR weight w (FORMS-H gate)",pal="magma",lim=c(0,1),nm="w")
fig<-(p1|p2)/(p3|p4)+plot_annotation(title="Wall-to-wall fused leaf-area map, Blois (10 m) — deployability, validated only at 53 loggers")
ggsave(file.path(OUT,"walltowall_fused_LAI.png"),fig,width=10,height=9,dpi=200,bg="white")
ggsave(file.path(OUT,"walltowall_fused_LAI.pdf"),fig,width=10,height=9,device=cairo_pdf,bg="white")
cat(sprintf("\nwrote %s/walltowall_fused_LAI.{png,pdf} + .tif rasters\n",OUT))
