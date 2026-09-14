# ==============================================================================
# QUICK: variable importance on the MEAN PROFILE of each NEW FPC cluster
# (archetype-centroid sensitivity), 4 centroids x 8 tags = 32 sims. Legacy v3.2.0
# + ERA5 (apples-to-apples). ISOLATED. Runs while the 400-plot full re-run goes.
#   Rscript c1_fpc_archetype_sensitivity.R
# Out: out_files/Chapter1/nc_archetype_fpc/ + tab_fpc_archetype_sensitivity.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
CFG_C3$musica_cmd <- "/home/corroyez/Documents/musica/musica"; stopifnot(file.exists(CFG_C3$musica_cmd))

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_fpc_v1.rds"))
ladc <- grep("LAD_Layer_", names(samp), value=TRUE)
# centroid per archetype: mean traits + mean LAD profile
cent <- samp[, c(.(LAI=mean(LAI), Hmax=mean(Hmax), fCover=mean(fCover), VCI=mean(VCI), x=0, y=0),
                 lapply(.SD, mean, na.rm=TRUE)), by=P, .SDcols=ladc][order(P)]
SDt  <- samp[, .(sLAI=sd(LAI/2,na.rm=T), sFC=sd(fCover,na.rm=T)/0.1, sHM=sd(Hmax,na.rm=T)/5, sVCI=sd(VCI,na.rm=T)), by=P]

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day"); Z_FIX<-1.0; SHIFT<-2L
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds))
metrics_one <- function(path){
  if(!file.exists(path)||file.size(path)<1000) return(NA_real_)
  nc<-try(nc_open(path),silent=TRUE); if(inherits(nc,"try-error")) return(NA_real_); on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var))) return(NA_real_)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
  vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh; Z<-Z_FIX
  if(Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  tv<-t0+dhours(th)-lubridate::hours(SHIFT); Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
  dd<-data.table(date=as.Date(floor_date(tv,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
  merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }

NCDIR<-"out_files/Chapter1/nc_archetype_fpc"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
DLAI<-2; DHMAX<-5; DFC<-0.1; LAI_FLOOR<-0.1; HMAX_FLOOR<-3
run_get <- function(prow, tag, lai, hmax, fcov, ladf=make_lad_real){
  nc<-file.path(NCDIR, sprintf("%s_%s.nc", prow$P, tag))
  sc<-list(lai_fn=function(p) lai, hmax_fn=function(p) hmax, fcover_fn=function(p) fcov, lad_fn=ladf, phenology_fn=NULL)
  run_musica_one(prow, sc, nc, CFG_C3$forcing_file, CFG_C3$musica_cmd)
  if(!file.exists(nc)||file.size(nc)<1000) return(NA_real_); metrics_one(nc) }

rows <- rbindlist(lapply(seq_len(nrow(cent)), function(i){
  prow <- as.data.frame(cent[i])
  base<-run_get(prow,"base",prow$LAI,prow$Hmax,prow$fCover)
  unif<-run_get(prow,"unifLAD",prow$LAI,prow$Hmax,prow$fCover,ladf=make_lad_uniform)
  laip<-run_get(prow,"LAIp",prow$LAI+DLAI,prow$Hmax,prow$fCover)
  laim<-run_get(prow,"LAIm",max(prow$LAI-DLAI,LAI_FLOOR),prow$Hmax,prow$fCover)
  hmp <-run_get(prow,"Hmaxp",prow$LAI,prow$Hmax+DHMAX,prow$fCover)
  hmm <-run_get(prow,"Hmaxm",prow$LAI,max(prow$Hmax-DHMAX,HMAX_FLOOR),prow$fCover)
  fcp <-run_get(prow,"fCovp",prow$LAI,prow$Hmax,min(prow$fCover+DFC,1))
  fcm <-run_get(prow,"fCovm",prow$LAI,prow$Hmax,max(prow$fCover-DFC,0.5))
  laim_v<-max(prow$LAI-DLAI,LAI_FLOOR); hmm_v<-max(prow$Hmax-DHMAX,HMAX_FLOOR); fcm_v<-max(prow$fCover-DFC,0.5); fcp_v<-min(prow$fCover+DFC,1)
  data.table(P=prow$P, LAI=prow$LAI, Hmax=prow$Hmax, fCover=prow$fCover, VCI=prow$VCI, base=base,
    natLAI=(laip-laim)/((DLAI+(prow$LAI-laim_v))/2)/2,    # ~per +1 one-sided LAI (central)
    natLAD=base-unif,
    natfC=(fcp-fcm)/(( (fcp_v-prow$fCover)+(prow$fCover-fcm_v) )/DFC)/2,
    natHM=(hmp-hmm)/((DHMAX+(prow$Hmax-hmm_v))/5)/2)
}))
rows <- merge(rows, SDt, by="P")
rows[, `:=`(LAI_perSD=round(natLAI*sLAI,3), fCover_perSD=round(natfC*sFC,3),
            Hmax_perSD=round(natHM*sHM,3), LAD_perSD=round(natLAD/pmax(1-VCI,0.05)*sVCI,3))]
out <- rows[, .(P, LAI=round(LAI,2), base=round(base,2),
                LAI_perSD, fCover_perSD, Hmax_perSD, LAD_perSD,
                natLAI=round(natLAI,3), natLAD=round(natLAD,3))][order(P)]
cat("\n=== FPC archetype-CENTROID sensitivity (mean profile per cluster) — ΔTmax ===\n"); print(out)
fwrite(out, "out_files/Chapter1/tables/tab_fpc_archetype_sensitivity.csv")
cat("\nper-SD = native slope at the mean profile x within-cluster SD. Compare LAI_perSD vs LAD_perSD in P4 (co-lead?).\n")
