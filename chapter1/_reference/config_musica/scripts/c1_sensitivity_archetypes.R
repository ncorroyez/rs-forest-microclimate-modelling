# ==============================================================================
# Sensibilité par unité "P1-P4 ONLY" : sur les 4 PROFILS MOYENS d'archétype
# (moyennes LAI/Hmax/fCover + profil LAD moyen par cluster), complément du readout
# cLHS (400 plots). Mêmes perturbations (±1 LAI, ±0.1 fCover, ±5 m Hmax + LAD
# réel-vs-uniforme), pour v3.2.0 (legacy) ET v3.2.3 iter (yoyo, vraie BLH).
# 4 archétypes × 8 simus × 2 versions = 64 simus. Resumable.
#   Rscript c1_sensitivity_archetypes.R
# Out: tab_sensitivity_archetypes.csv + FigCmp_sensitivity_archetypes_v320_v323iter.png
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(ggplot2)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

# ---- 4 profils moyens d'archétype ----
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp[, P := relabel_cluster(Cluster)]
ladc <- grep("^LAD_Layer", names(samp), value=TRUE)
arch <- samp[, c(.(LAI=mean(LAI), Hmax=mean(Hmax), fCover=mean(fCover)),
                 lapply(.SD, mean)), .SDcols=ladc, by=P][order(P)]

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds)); Z<-1; SH<-2L
dT <- function(p){ if(!file.exists(p)||file.size(p)<1e5) return(NA_real_)
  nc<-try(nc_open(p),silent=TRUE); if(inherits(nc,"try-error"))return(NA_real_); on.exit(nc_close(nc))
  if(!all(c("Tair_z","relative_height","veget_height_top")%in%names(nc$var)))return(NA_real_)
  tu<-ncatt_get(nc,"time","units")$value; t0<-as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th<-ncvar_get(nc,"time"); Tk<-ncvar_get(nc,"Tair_z"); rh<-ncvar_get(nc,"relative_height")
  vh<-stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE); zl<-rh*vh
  if(Z<=zl[1]){il<-1L;ih<-1L;w<-0}else if(Z>=zl[length(zl)]){il<-length(zl);ih<-il;w<-0}else{il<-max(which(zl<=Z));ih<-il+1L;w<-(Z-zl[il])/(zl[ih]-zl[il])}
  tv<-t0+dhours(th)-lubridate::hours(SH); Tc<-((1-w)*Tk[il,]+w*Tk[ih,])-273.15
  dd<-data.table(date=as.Date(floor_date(tv,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
  merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }

mk_phen <- function(lai1) function(p) {
  ph <- as.data.frame(calc_phenology(list.year=2020:2022, nleafage=1, budburst_date=115,
        leaf_age_max_in=0.56, relative_age_firstmax=0.10, relative_age_lastmax=0.75, LAI_max_per_cohort=lai1))
  d <- ph[ph$year==2020 & ph$Julian_day==365,,drop=FALSE]; d$Julian_day<-366; rbind(ph,d) }

CFG_V <- list(
  "v3.2.0"    = list(bin="/home/corroyez/Documents/musica/musica", forc=CFG_C3$forcing_file, abl=NULL, phen=FALSE),
  "v3.2.3 iter"= list(bin=normalizePath("in_files/model-3.2.3/musica"), forc="in_files/musica_in_Blois_pblh.nc",
                      abl=list("abl_flag"='"iter"'), phen=TRUE))
DLAI<-2; DHMAX<-5; DFC<-0.1; LAI_FLOOR<-0.1; HMAX_FLOOR<-3
nrm <- function(x,s) if(is.finite(x)&&is.finite(s)&&s>1e-9) x/s else NA_real_

run_one <- function(ver, prow, tag, lai, hmax, fcov, ladf=make_lad_real) {
  cfg <- CFG_V[[ver]]; d <- sprintf("out_files/Chapter3/nc_archetypes/%s", gsub("[^A-Za-z0-9]","",ver))
  dir.create(d, recursive=TRUE, showWarnings=FALSE); nc <- file.path(d, sprintf("%s_%s.nc", prow$P, tag))
  sc <- list(lai_fn=function(p)lai, hmax_fn=function(p)hmax, fcover_fn=function(p)fcov, lad_fn=ladf,
             phenology_fn=if(cfg$phen) mk_phen(lai) else NULL)
  run_musica_one(prow, sc, nc, cfg$forc, cfg$bin, extra_setup=if(is.null(cfg$abl)) list() else cfg$abl)
  if(!file.exists(nc)||file.size(nc)<1e5) return(NA_real_); dT(nc)
}

res <- rbindlist(lapply(names(CFG_V), function(ver) rbindlist(lapply(seq_len(nrow(arch)), function(i){
  prow <- as.data.frame(arch[i])
  base <- run_one(ver, prow, "base", prow$LAI, prow$Hmax, prow$fCover)
  unif <- run_one(ver, prow, "unifLAD", prow$LAI, prow$Hmax, prow$fCover, make_lad_uniform)
  laim_v<-max(prow$LAI-DLAI,LAI_FLOOR); hmm_v<-max(prow$Hmax-DHMAX,HMAX_FLOOR); fcp_v<-min(prow$fCover+DFC,1); fcm_v<-max(prow$fCover-DFC,0.5)
  laip<-run_one(ver,prow,"LAIp",prow$LAI+DLAI,prow$Hmax,prow$fCover); laim<-run_one(ver,prow,"LAIm",laim_v,prow$Hmax,prow$fCover)
  hmp<-run_one(ver,prow,"Hmaxp",prow$LAI,prow$Hmax+DHMAX,prow$fCover); hmm<-run_one(ver,prow,"Hmaxm",prow$LAI,hmm_v,prow$fCover)
  fcp<-run_one(ver,prow,"fCovp",prow$LAI,prow$Hmax,fcp_v); fcm<-run_one(ver,prow,"fCovm",prow$LAI,prow$Hmax,fcm_v)
  data.table(version=ver, P=prow$P, base=base,
    LAI =nrm(laip-base, DLAI/2),  LAI_rem =nrm(laim-base,(prow$LAI-laim_v)/2),
    Hmax=nrm(hmp-base, DHMAX/5),  Hmax_rem=nrm(hmm-base,(prow$Hmax-hmm_v)),
    fCov=nrm(fcp-base,(fcp_v-prow$fCover)/DFC), fCov_rem=nrm(fcm-base,(prow$fCover-fcm_v)/DFC),
    LAD =base-unif)
}))))
fwrite(res, "out_files/Chapter3/tables/tab_sensitivity_archetypes.csv")
cat("=== sensibilité P1-P4 only (effet par unité, ΔTmax °C) ===\n"); print(res)

# magnitude per-unit symétrique pour la figure
M <- res[, .(LAI=abs((LAI-LAI_rem)/2), fCover=abs((fCov-fCov_rem)/2)/0.1*0.1, Hmax=abs((Hmax-Hmax_rem)/2),
             LAD=abs(LAD)), by=.(version,P)]
L <- melt(M, id.vars=c("version","P"), variable.name="trait", value.name="val")
L[, trait:=factor(trait, levels=c("LAI","fCover","LAD","Hmax"))]
L[, version:=factor(version, levels=names(CFG_V))]
PAL <- c(LAI="#1B7837",fCover="#7FBC41",LAD="#762A83",Hmax="#BDBDBD")
p <- ggplot(L, aes(P,val,fill=trait)) + geom_col(position=position_dodge(0.8),width=0.74) +
  geom_text(aes(label=sprintf("%.2f",val)),position=position_dodge(0.8),vjust=-0.3,size=2.3,colour="grey30") +
  facet_wrap(~version,nrow=1) + scale_fill_manual(values=PAL,name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL, y="|effet par unité| sur ΔTmax (°C)",
       title="Sensibilité « P1-P4 only » (4 profils moyens d'archétype) : v3.2.0 vs v3.2.3 iter",
       subtitle="Effet par unité au point de fonctionnement MOYEN de chaque archétype (LAI/unit, fCover/0.1, Hmax/5m, LAD réel-vs-uniforme). Complément du readout cLHS (400 plots).") +
  theme_bw(base_size=12) + theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
    legend.position="bottom",strip.text=element_text(face="bold"),plot.title=element_text(size=10.5,face="bold"),
    plot.subtitle=element_text(size=8,colour="grey35"))
out <- "out_files/Chapter3/figures/FigCmp_sensitivity_archetypes_v320_v323iter.png"
ggsave(out, p, width=9.5, height=5, dpi=200, bg="white")
file.copy(out, "Chapitre1/comparaison_versions/figures/", overwrite=TRUE)
cat("DONE ->", out, "\n")
