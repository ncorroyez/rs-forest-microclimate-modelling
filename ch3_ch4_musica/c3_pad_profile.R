# ==============================================================================
# B4 — mean vertical PAD profile, dense vs open plots, with the effective optical
# depth d_opt = 7 m marked, and the fraction of leaf area BELOW d_opt (the layer
# Sentinel-2 cannot see). Turns the "sub-d_opt canopy buffers too" claim into a
# quantified figure. Reads the LAD_Layer_*.tif stack (Not_Masked), extracts at the
# 53 plots, splits at the chosen near-median regime threshold (LAI 3.86).
#   Rscript c3_pad_profile.R  -> tables/figL_pad_profile.csv + figures/FigL_pad_profile
# ==============================================================================
suppressPackageStartupMessages({library(terra);library(stringr);library(data.table)})
source("Chapter3_config.R")
NM<-"/home/corroyez/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked"
TAB<-"/home/corroyez/Documents/NC_Full/manuscripts/ch3/tables"
DOPT<-CFG_C3$d_opt_m; CROSS<-3.86
df<-as.data.table(readRDS(file.path(CFG_C3$out_dir,"lai_prep","df_plots_real53.rds")))
pts<-vect(as.data.frame(df[,.(x,y)]),geom=c("x","y"),crs="EPSG:32631")
files<-list.files(NM,pattern="^LAD_Layer_[0-9.]+_res_10_m\\.tif$",full.names=TRUE)
h<-as.numeric(str_extract(basename(files),"(?<=LAD_Layer_)[0-9.]+")); files<-files[order(h)]; h<-sort(h)
cat(sprintf("%d LAD layers, heights %.1f-%.1f m\n",length(h),min(h),max(h)))
val<-as.data.frame(terra::extract(rast(files),pts))[,-1,drop=FALSE]  # 53 x nlayers, PAD per layer
val[is.na(val)]<-0
de<-df$LAI_ALS>=CROSS
# per-plot canopy top (highest layer with PAD>0.01) and fraction below d_opt
frac_below<-sapply(1:nrow(val),function(i){p<-as.numeric(val[i,]); tot<-sum(p); if(tot<=0)return(NA)
  top<-max(h[p>0.01],na.rm=TRUE); below<-p[h < (top-DOPT)]; sum(below)/tot})
cat(sprintf("Fraction of LAI BELOW d_opt (invisible to S2):\n  dense mean %.2f | open mean %.2f\n",
    mean(frac_below[de],na.rm=T), mean(frac_below[!de],na.rm=T)))
# mean profile by stratum (PAD per height, normalized as density)
prof<-rbind(
  data.table(height=h, pad=colMeans(val[de,,drop=FALSE]), stratum=sprintf("Dense (LAI>=3.86, n=%d)",sum(de))),
  data.table(height=h, pad=colMeans(val[!de,,drop=FALSE]), stratum=sprintf("Open (LAI<3.86, n=%d)",sum(!de))))
fwrite(prof,file.path(TAB,"figL_pad_profile.csv"))
# top of dense canopy for d_opt band
top_dense<-mean(sapply(which(de),function(i){p<-as.numeric(val[i,]);max(h[p>0.01],na.rm=TRUE)}))
cat(sprintf("mean dense canopy top = %.1f m -> d_opt sees %.1f-%.1f m\n",top_dense,top_dense-DOPT,top_dense))
fwrite(data.table(dense_top=top_dense,dopt=DOPT,frac_below_dense=mean(frac_below[de],na.rm=T),
  frac_below_open=mean(frac_below[!de],na.rm=T)),file.path(TAB,"figL_pad_summary.csv"))
cat("wrote figL_pad_profile.csv + figL_pad_summary.csv\n")
