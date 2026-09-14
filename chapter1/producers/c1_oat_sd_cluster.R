# ==============================================================================
# Réu 2026-06-26 (nouveau scénario) — Sensibilité OAT par ±SD propre au cluster.
# Pour chaque archétype P1-P4 : baseline = centroïde cLHS (moyenne LAI/Hmax/fCover
# + profil LAD moyen du cluster). On fait varier UN trait à la fois de ± la SD
# intra-cluster de ce trait (les autres figés au centroïde), et on lit ΔTmax.
#   => barres d'erreur (#11) + sensibilité du modèle à la structure (#13)
#   => "la variation n'est pas la même par cluster" : SD propre à chaque cluster.
# Binaire validé v3.2.0, ΔTmax fixe 1 m / −2 h / vs macro (HOBO-comparable).
#   Rscript c1_oat_sd_cluster.R
# Out: tab_oat_sd_cluster.csv (+ tab_cluster_meansd.csv) + nc_oat_sd_cluster/
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr); library(purrr)
  library(parallel); library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")
CFG_C3$musica_cmd <- "/home/corroyez/Documents/musica/musica"   # validated legacy v3.2.0

# ---- per-cluster mean & SD (LAI/Hmax/fCover/VCI), keyed on RAW Cluster code ---
# Cluster codes 1-4 are remapped to P1(open)..P4(dense) by relabel_cluster
# (ascending LAI): 3→P1, 2→P2, 4→P3, 1→P4. LAD shape = cluster-mean profile
# rescaled to each (hmax,lai) via make_lad_cluster_type_factory.
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
LADF <- make_lad_cluster_type_factory(as.data.frame(samp))   # keys on raw Cluster
ms <- samp[, .(LAI_m=mean(LAI), LAI_sd=sd(LAI), Hmax_m=mean(Hmax), Hmax_sd=sd(Hmax),
               fCover_m=mean(fCover), fCover_sd=sd(fCover), VCI_m=mean(VCI), VCI_sd=sd(VCI),
               x=mean(x), y=mean(y)), by=Cluster]
ms[, P := relabel_cluster(Cluster)]; setorder(ms, P)
fwrite(ms[, .(P,Cluster,LAI_m,LAI_sd,Hmax_m,Hmax_sd,fCover_m,fCover_sd,VCI_m,VCI_sd)],
       "out_files/Chapter1/tables/tab_cluster_meansd.csv")
cat("=== centroïde ± SD par cluster ===\n")
print(ms[, .(P, LAI=sprintf("%.2f±%.2f",LAI_m,LAI_sd), Hmax=sprintf("%.1f±%.1f",Hmax_m,Hmax_sd),
             fCover=sprintf("%.3f±%.3f",fCover_m,fCover_sd))])

# ---- ΔTmax extractor (identique c1_sensitivity_perplot_chunk) -----------------
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by="day")
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds)); Z_FIX <- 1.0; SHIFT <- 2L
metrics_one <- function(path) {
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(NA_real_); on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu),tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"),na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]) {ilo<-1L;ihi<-1L;w<-0} else if (Z>=zl[length(zl)]) {ilo<-length(zl);ihi<-ilo;w<-0} else
    {ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec <- t0+dhours(th)-lubridate::hours(SHIFT); Tc <- ((1-w)*Tk[ilo,]+w*Tk[ihi,])-273.15
  dd <- data.table(date=as.Date(floor_date(tvec,"hour")),Tc=Tc)[date%in%ds,.(Tmax=max(Tc)),by=date]
  merge(dd,dm,by="date")[,mean(Tmax-Tmax_macro,na.rm=TRUE)] }

NCDIR <- "out_files/Chapter1/nc_oat_sd_cluster"; dir.create(NCDIR, recursive=TRUE, showWarnings=FALSE)
LAI_FLOOR <- 0.1; HMAX_FLOOR <- 3

# centroid plot_row per cluster : raw Cluster (for LADF) + plain trait columns
crow <- function(p) {
  m <- ms[P==p]
  data.frame(Cluster=m$Cluster, LAI=m$LAI_m, Hmax=m$Hmax_m, fCover=m$fCover_m,
             VCI=m$VCI_m, x=m$x, y=m$y)
}

# ---- build job list : per cluster, baseline + ±SD on each of 3 traits --------
jobs <- list()
for (p in ms$P) {
  m <- ms[P==p]
  add <- function(tag,lai,hm,fc) jobs[[length(jobs)+1]] <<- list(P=p,tag=tag,lai=lai,hm=hm,fc=fc)
  add("base",     m$LAI_m,                       m$Hmax_m,                      m$fCover_m)
  add("LAI_p",    m$LAI_m + m$LAI_sd,            m$Hmax_m,                      m$fCover_m)
  add("LAI_m",    max(m$LAI_m - m$LAI_sd, LAI_FLOOR), m$Hmax_m,                 m$fCover_m)
  add("Hmax_p",   m$LAI_m,                       m$Hmax_m + m$Hmax_sd,          m$fCover_m)
  add("Hmax_m",   m$LAI_m,                       max(m$Hmax_m - m$Hmax_sd, HMAX_FLOOR), m$fCover_m)
  add("fCover_p", m$LAI_m,                       m$Hmax_m,                      min(m$fCover_m + m$fCover_sd, 1))
  add("fCover_m", m$LAI_m,                       m$Hmax_m,                      max(m$fCover_m - m$fCover_sd, 0.5))
}
cat(sprintf("\n%d sims (%d clusters × 7)\n", length(jobs), nrow(ms)))

run_one <- function(jb) {
  nc <- file.path(NCDIR, sprintf("%s_%s.nc", jb$P, jb$tag))
  pr <- crow(jb$P)
  if (!(file.exists(nc) && file.size(nc) > 1000)) {
    sc <- list(lai_fn=function(p) jb$lai, hmax_fn=function(p) jb$hm, fcover_fn=function(p) jb$fc,
               lad_fn=LADF, phenology_fn=NULL)   # cluster-mean LAD profile, rescaled to (hmax,lai)
    tryCatch(run_musica_one(pr, sc, nc, CFG_C3$forcing_file, CFG_C3$musica_cmd),
             error=function(e) cat(sprintf("ERR %s_%s: %s\n", jb$P, jb$tag, e$message)))
  }
  dT <- if (file.exists(nc) && file.size(nc) > 1000) metrics_one(nc) else NA_real_
  data.table(P=jb$P, tag=jb$tag, LAI=jb$lai, Hmax=jb$hm, fCover=jb$fc, dTmax=dT)
}
res <- rbindlist(mclapply(jobs, run_one, mc.cores=4, mc.preschedule=FALSE))
fwrite(res, "out_files/Chapter1/tables/tab_oat_sd_cluster.csv")
cat("\n=== ΔTmax par sim ===\n"); print(res[order(P,tag)])
cat(sprintf("\nnc produits : %d/%d\n", sum(file.exists(file.path(NCDIR, paste0(res$P,"_",res$tag,".nc")))), nrow(res)))
cat("DONE -> tab_oat_sd_cluster.csv + tab_cluster_meansd.csv\n")
