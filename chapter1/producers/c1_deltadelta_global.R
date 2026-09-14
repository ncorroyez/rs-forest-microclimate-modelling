# ==============================================================================
# Réu 2026-06-26 #12 — Delta/delta ΔTmax / Δtrait (PAS de GLM), échelle globale
# ET par cluster. Per-unit = (ΔTmax(+SD) − ΔTmax(−SD)) / (2·SD_native), en
# °C par unité native du trait. Les 4 clusters viennent de tab_oat_sensitivity ;
# le point global = 7 sims supplémentaires au centroïde global (LAD = moyen).
#   Rscript c1_deltadelta_global.R
# Out: tab_deltadelta.csv + FigStation_deltadelta.png
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(purrr)
  library(parallel); library(ggplot2); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
CFG_C3$musica_cmd <- "/home/corroyez/Documents/musica/musica"
PALt <- c(LAI="#1B7837", Hmax="#9E9E9E", fCover="#7FBC41"); TR <- c("LAI","Hmax","fCover")

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
LADM <- make_lad_mean_factory(as.data.frame(samp))    # global mean LAD shape
g <- samp[, .(LAI=mean(LAI), LAI_sd=sd(LAI), Hmax=mean(Hmax), Hmax_sd=sd(Hmax),
              fCover=mean(fCover), fCover_sd=sd(fCover), x=mean(x), y=mean(y))]

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds)); Z_FIX<-1.0; SHIFT<-2L
metrics_one <- function(path){ nc<-try(nc_open(path),silent=TRUE); if(inherits(nc,"try-error"))return(NA_real_); on.exit(nc_close(nc))
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
  vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh; Z<-Z_FIX
  if(Z<=zl[1]){ilo<-1L;ihi<-1L;w<-0}else if(Z>=zl[length(zl)]){ilo<-length(zl);ihi<-ilo;w<-0}else{ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec<-t0+dhours(th)-lubridate::hours(SHIFT); Tc<-((1-w)*Tk[ilo,]+w*Tk[ihi,])-273.15
  dd<-data.table(date=as.Date(floor_date(tvec,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
  merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }

NCDIR<-"out_files/Chapter1/nc_oat_global"; dir.create(NCDIR,recursive=TRUE,showWarnings=FALSE)
pr <- data.frame(LAI=g$LAI, Hmax=g$Hmax, fCover=g$fCover, x=g$x, y=g$y)
runG <- function(tag,lai,hm,fc){ nc<-file.path(NCDIR,sprintf("All_%s.nc",tag))
  if(!(file.exists(nc)&&file.size(nc)>1000)){
    sc<-list(lai_fn=function(p)lai,hmax_fn=function(p)hm,fcover_fn=function(p)fc,lad_fn=LADM,phenology_fn=NULL)
    tryCatch(run_musica_one(pr,sc,nc,CFG_C3$forcing_file,CFG_C3$musica_cmd),error=function(e)cat("ERR",tag,e$message,"\n"))}
  if(file.exists(nc)&&file.size(nc)>1000) metrics_one(nc) else NA_real_ }
J<-list(c("base",g$LAI,g$Hmax,g$fCover),
        c("LAI_p",g$LAI+g$LAI_sd,g$Hmax,g$fCover), c("LAI_m",max(g$LAI-g$LAI_sd,0.1),g$Hmax,g$fCover),
        c("Hmax_p",g$LAI,g$Hmax+g$Hmax_sd,g$fCover), c("Hmax_m",g$LAI,max(g$Hmax-g$Hmax_sd,3),g$fCover),
        c("fCover_p",g$LAI,g$Hmax,min(g$fCover+g$fCover_sd,1)), c("fCover_m",g$LAI,g$Hmax,max(g$fCover-g$fCover_sd,0.5)))
gv<-rbindlist(mclapply(J,function(j)data.table(tag=j[1],dTmax=runG(j[1],as.numeric(j[2]),as.numeric(j[3]),as.numeric(j[4]))),mc.cores=4))
gd<-function(t) (gv[tag==paste0(t,"_p"),dTmax]-gv[tag==paste0(t,"_m"),dTmax])/(2*g[[paste0(t,"_sd")]])
glob <- data.table(P="All", trait=TR, per_unit=sapply(TR,gd), base=gv[tag=="base",dTmax])

# combine with per-cluster (#11/#13 table)
cl <- fread("out_files/Chapter1/tables/tab_oat_sensitivity.csv")[, .(P=as.character(P), trait, per_unit, base)]
D <- rbind(glob, cl); D[, P := factor(P, levels=c("All","P1","P2","P3","P4"))][, trait:=factor(trait,levels=TR)]
fwrite(D[order(trait,P)], "out_files/Chapter1/tables/tab_deltadelta.csv")
cat("=== delta/delta ΔTmax par unité native (global + clusters) ===\n")
print(dcast(D, P~trait, value.var="per_unit")[order(P)])

UNI <- c(LAI="par +1 LAI", Hmax="par +1 m Hmax", fCover="par +0,1 fCover")
D[, pu := ifelse(trait=="fCover", per_unit*0.1, per_unit)]   # fCover affiché par +0.1
p <- ggplot(D, aes(P, pu, colour=trait)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_segment(aes(xend=P, yend=0), linewidth=0.9) +
  geom_point(aes(shape=P=="All"), size=3.4) +
  scale_shape_manual(values=c(`FALSE`=16, `TRUE`=18), guide="none") +
  facet_wrap(~trait, scales="free_y", nrow=1,
     labeller=labeller(trait=function(x) paste0(x, "\n(ΔTmax ", UNI[x], ")"))) +
  scale_colour_manual(values=PALt, guide="none") +
  labs(x=NULL, y="ΔTmax / Δtrait (°C par unité native)",
       title="Delta/delta : sensibilité de ΔTmax par unité de trait — global ET par cluster",
       subtitle="Sans GLM. Négatif = refroidit. LAI = levier dominant et décroissant ; Hmax ≈ 0 ; fCover surtout en P1 ouvert") +
  theme_bw(base_size=11) + theme(strip.text=element_text(face="bold",size=8.5),
       plot.title=element_text(size=10.5,face="bold"), plot.subtitle=element_text(size=8,colour="grey35"))
ggsave("out_files/Chapter1/figures/FigStation_deltadelta.png", p, width=11, height=4.4, dpi=200, bg="white")
cat("DONE -> tab_deltadelta.csv + FigStation_deltadelta.png\n")
