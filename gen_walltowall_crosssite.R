# ==============================================================================
# Chapter 3 — CROSS-SITE wall-to-wall fused LAI (deployability across forest types).
# Applies the SAME Blois-calibrated fusion rule (hard gate + logistic blend,
# H0=17.8 m, s=3, FIXED) per-pixel at 10 m over the three sites:
#   Blois (lowland oak), Aigoual (montane beech), Mormal (lowland oak-beech).
# Sources per site: LiDAR = Metrics/Not_Masked/lidarlai_res_10_m.tif ;
#   S2 = summer (DOY 152-244) mean of ATBD dated rasters ;
#   FORMS-H gate = FORMS-H_Blois.tif for Blois (m), else cropped from the France
#   product FORMS-H_Height_10m_cm.tif (Lambert-93, cm -> /100 = m).
# FRAMING: the 17.8 m gate is CALIBRATED ON BLOIS; its transfer to montane beech
# (Aigoual) and oak-beech (Mormal) is UNTESTED -> this is a DEPLOYABILITY demo, not
# validation (microclimate loggers exist only at Blois). Multi-site TRAINING of a
# transferable gate is the perspective this enables.
#   Rscript gen_walltowall_crosssite.R
# ==============================================================================
suppressPackageStartupMessages({library(terra);library(stringr);library(ggplot2);library(tidyterra);library(patchwork)})
RES<-"/home/corroyez/Documents/NC_Full/03_RESULTS"; DAT<-"/home/corroyez/Documents/NC_Full/01_DATA"
OUT<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/figures"; dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
H0<-17.8; S0<-3; FR<-rast(file.path(DAT,"FORMS-H_Height_10m_cm.tif"))   # France FORMS-H, cm, L93

site_layers<-function(S){
  lidar<-rast(file.path(RES,S,"Metrics/Not_Masked/lidarlai_res_10_m.tif")); names(lidar)<-"lidar"
  lidar<-clamp(lidar,0,8,values=FALSE)   # NA implausible LAI>8 (non-forest/edge artefacts: Aigoual 13.5%)
  s2f<-list.files(file.path(RES,S,"Metrics/Not_Masked"),pattern="^s2lai_\\d{4}-\\d{2}-\\d{2}_atbd_res_10_m\\.tif$",full.names=TRUE)
  doy<-as.integer(format(as.Date(str_extract(basename(s2f),"\\d{4}-\\d{2}-\\d{2}")),"%j"))
  s2<-app(rast(s2f[doy>=152&doy<=244]),fun=function(x)mean(pmax(x,0),na.rm=TRUE)); names(s2)<-"s2"
  if(S=="Blois"){ fh<-rast(file.path(DAT,"FORMS-H_Blois.tif")); fh<-resample(fh,lidar,method="bilinear") }
  else { e<-project(as.polygons(ext(lidar),crs=crs(lidar)),crs(FR))          # site extent in L93
         fc<-crop(FR,ext(e)+200)/100                                          # crop France, cm->m
         fh<-project(fc,lidar,method="bilinear") }                           # to site UTM grid
  names(fh)<-"formsh"; s2<-resample(s2,lidar,method="bilinear")
  w<-1/(1+exp(-(fh-H0)/S0)); names(w)<-"w"
  bl<-w*lidar+(1-w)*s2; names(bl)<-"blend"; sw<-ifel(fh>=H0,lidar,s2); names(sw)<-"switch"
  m<-!is.na(lidar)&!is.na(s2)&!is.na(fh); bl<-mask(bl,m,maskvalues=FALSE); sw<-mask(sw,m,maskvalues=FALSE); w<-mask(w,m,maskvalues=FALSE)
  writeRaster(bl,file.path(OUT,sprintf("walltowall_%s_blend.tif",S)),overwrite=TRUE)
  cat(sprintf("%-8s LiDAR[%.1f..%.1f] S2[%.1f..%.1f] FORMS-H[%.1f..%.1f] wmean=%.2f  frac dense(w>.5)=%.0f%%\n",
      S,minmax(lidar)[1],minmax(lidar)[2],minmax(s2)[1],minmax(s2)[2],minmax(fh)[1],minmax(fh)[2],
      global(w,"mean",na.rm=TRUE)[1,1],100*global(w>0.5,"mean",na.rm=TRUE)[1,1]))
  list(lidar=lidar,s2=s2,blend=bl,w=w,site=S)
}
mk<-function(r,ttl,lim=c(0,7.5),pal="viridis",nm="LAI"){ggplot()+geom_spatraster(data=r)+
  scale_fill_viridis_c(option=pal,limits=lim,na.value="transparent",name=nm)+labs(subtitle=ttl)+
  coord_sf(expand=FALSE)+theme_minimal(base_size=9)+
  theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),legend.key.height=unit(0.5,"cm"))}

L<-lapply(c("Blois","Aigoual","Mormal"),site_layers)
rows<-lapply(L,function(x){
  (mk(x$lidar,sprintf("%s — LiDAR LAI",x$site)) | mk(x$s2,sprintf("%s — Sentinel-2 LAI",x$site)) |
   mk(x$blend,sprintf("%s — blended LAI",x$site)) | mk(x$w,sprintf("%s — LiDAR weight w",x$site),c(0,1),"magma","w"))})
fig<-wrap_plots(rows,ncol=1)+plot_annotation(
  title="Cross-site deployment of the Blois-calibrated LiDAR/Sentinel-2 fusion (10 m)",
  subtitle="Same gate (FORMS-H 17.8 m, s=3). Deployability demonstration — the threshold is Blois-calibrated and untested elsewhere; microclimate validation exists only at Blois.")
ggsave(file.path(OUT,"walltowall_crosssite.png"),fig,width=12,height=11,dpi=170,bg="white")
ggsave(file.path(OUT,"walltowall_crosssite.pdf"),fig,width=12,height=11,device=cairo_pdf,bg="white")
cat(sprintf("\nwrote %s/walltowall_crosssite.{png,pdf} + per-site blend tifs\n",OUT))
