# ==============================================================================
# Two methods-justification tests (supervisor 2026-06-18), designed cheap.
#
# TEST 1 — MEDIAN vs MEAN baseline. The Shapley baseline is parameterized by 3
#   scalars (cluster-central LAI/Hmax/fCover) + a uniform LAD (= LAI/depth, no
#   profile median). So mean-vs-median enters ONLY through those scalars; if
#   mean≈median the 15 baseline-containing coalitions are materially identical
#   and φ is unchanged (analytic argument). We CONFIRM by running the baseline
#   coalition (uniform LAD + central values) under mean vs median (legacy binary,
#   internal comparison -> binary choice irrelevant). Watch fCover (floored/capped).
#
# TEST 2 — cLHS distribution vs single MEAN archetype (Jensen). ΔTmax saturates in
#   LAI, so f(mean structure) ≠ mean_i f(structure_i). We compare:
#     f(mean archetype)  = mean LAI/Hmax/fCover + cluster-mean LAD shape, run on
#                          the ATTRIBUTION harness (Blois binary CFG_C3) to match
#     mean_i f(real_i)   = mean over cLHS plots of the cached `full` (1111) column
#                          in metrics_parts (also Blois binary).
#   A non-trivial gap (largest in open P1) justifies the cLHS-distribution as MAIN.
#   Out: out_files/Chapter1/tables/tab_methodo_tests.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))                 # CFG (legacy binary), run_musica_one, make_lad_*, era5/macro
source("Chapter3_config.R")                    # CFG_C3 (Blois binary = attribution harness)
TABDIR <- "out_files/Chapter1/tables"
NCDIR  <- "out_files/Chapter1/nc_methodo"; dir.create(NCDIR, recursive = TRUE, showWarnings = FALSE)
ds <- CFG$date_seq

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
ctr  <- samp[, .(LAI_mean=mean(LAI), Hmax_mean=mean(Hmax), fCov_mean=mean(fCover),
                 LAI_med=median(LAI), Hmax_med=median(Hmax), fCov_med=median(fCover)), by=Cluster][order(Cluster)]
LAD_CLUSTER <- make_lad_cluster_type_factory(as.data.frame(samp))

# ---- metric extractor (ported from c1_metrics_chunk.R) ----------------------
esat_hpa <- function(Tc) 6.108*exp(17.27*Tc/(Tc+237.3)); Z_FIX<-1.0; SHIFT<-2L; P_HPA<-1013; HOTQ<-0.90
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
hot  <- dm[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm=TRUE), date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds)); MET <- c("Tmax_all","slope_all","VPD_all")
metrics_one <- function(path) {
  out <- setNames(rep(NA_real_,3), MET)
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(out); on.exit(nc_close(nc))
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z"); wmr <- ncvar_get(nc,"wair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]) {ilo<-1L;ihi<-1L;w<-0} else if (Z>=zl[length(zl)]) {ilo<-length(zl);ihi<-ilo;w<-0} else {
    ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec <- t0+dhours(th)-lubridate::hours(SHIFT); Tc <- ((1-w)*Tk[ilo,]+w*Tk[ihi,])-273.15
  wv <- (1-w)*wmr[ilo,]+w*wmr[ihi,]; vpd <- pmax(esat_hpa(Tc)-(wv/(1+wv))*P_HPA,0)/10
  dd <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc, vpd=vpd)[date %in% ds, .(Tmax=max(Tc), VPDmax=max(vpd)), by=date]
  m1 <- merge(dd, dm, by="date"); out["Tmax_all"] <- mean(m1$Tmax-m1$Tmax_macro, na.rm=TRUE); out["VPD_all"] <- mean(dd$VPDmax, na.rm=TRUE)
  mic <- data.table(time=floor_date(t0+dhours(th),"hour"), Tmic=Tc); mm <- merge(mic[as.Date(time) %in% ds], era5, by="time")
  if (nrow(mm)>10) out["slope_all"] <- as.numeric(coef(lm(Tmic~Tair_era5, mm))[2]); out
}
run_get <- function(tag, cid, lai, hmax, fcov, shape, binary) {
  nc <- file.path(NCDIR, sprintf("%s_C%d.nc", tag, cid)); prow <- data.frame(Cluster=cid, x=0, y=0)
  sc <- list(lai_fn=function(p) lai, hmax_fn=function(p) hmax, fcover_fn=function(p) fcov,
             lad_fn=if (shape=="uniform") make_lad_uniform else LAD_CLUSTER, phenology_fn=NULL)
  run_musica_one(prow, sc, nc, CFG$forcing_file, binary); if (!file.exists(nc)) return(NULL); as.list(metrics_one(nc))
}

# ---- TEST 1: uniform-LAD baseline, MEAN vs MEDIAN central values (legacy) ----
t1 <- rbindlist(lapply(seq_len(nrow(ctr)), function(i) {
  cid <- ctr$Cluster[i]; b <- ctr[i]
  mn <- run_get("base_mean", cid, b$LAI_mean, b$Hmax_mean, b$fCov_mean, "uniform", CFG$musica_cmd)
  md <- run_get("base_med",  cid, b$LAI_med,  b$Hmax_med,  b$fCov_med,  "uniform", CFG$musica_cmd)
  data.table(Cluster=cid, Tmax_mean=mn$Tmax_all, Tmax_med=md$Tmax_all,
             VPD_mean=mn$VPD_all, VPD_med=md$VPD_all, slope_mean=mn$slope_all, slope_med=md$slope_all)
}))
t1[, `:=`(dTmax=Tmax_med-Tmax_mean, dVPD=VPD_med-VPD_mean, dslope=slope_med-slope_mean)]
fwrite(t1, file.path(TABDIR, "tab_methodo_test1_median_vs_mean.csv"))
cat("\n=== TEST 1: baseline MEDIAN - MEAN (uniform LAD, per cluster) ===\n")
print(t1[, .(Cluster, Tmax_mean=round(Tmax_mean,3), Tmax_med=round(Tmax_med,3), dTmax=round(dTmax,3), dVPD=round(dVPD,3), dslope=round(dslope,4))])

# ---- TEST 2: f(mean archetype) vs mean_i f(real_i), Blois harness ------------
parts <- rbindlist(lapply(list.files(file.path(TABDIR,"metrics_parts"),"part_.*csv$",full.names=TRUE), fread))
parts[, Cluster := as.integer(as.character(Cluster))]
meanf <- parts[metric=="Tmax_all", .(meanf_Tmax=mean(full, na.rm=TRUE), n=.N), by=Cluster][order(Cluster)]
mv    <- parts[metric=="VPD_all", .(meanf_VPD=mean(full, na.rm=TRUE)), by=Cluster]
ms    <- parts[metric=="slope_all", .(meanf_slope=mean(full, na.rm=TRUE)), by=Cluster]
t2 <- rbindlist(lapply(seq_len(nrow(ctr)), function(i) {
  cid <- ctr$Cluster[i]; b <- ctr[i]
  fm <- run_get("arch_mean", cid, b$LAI_mean, b$Hmax_mean, b$fCov_mean, "cluster", CFG_C3$musica_cmd)
  data.table(Cluster=as.integer(as.character(cid)), fmean_Tmax=fm$Tmax_all, fmean_VPD=fm$VPD_all, fmean_slope=fm$slope_all)
}))
t2 <- Reduce(function(a,b) merge(a,b,by="Cluster"), list(t2, meanf, mv, ms))
t2[, `:=`(gap_Tmax=fmean_Tmax-meanf_Tmax, gap_VPD=fmean_VPD-meanf_VPD, gap_slope=fmean_slope-meanf_slope)]
fwrite(t2, file.path(TABDIR, "tab_methodo_test2_clhs_vs_archetype.csv"))
cat("\n=== TEST 2: f(mean archetype) vs mean_i f(real_i) [Jensen], per cluster ===\n")
print(t2[, .(Cluster, n, fmean_Tmax=round(fmean_Tmax,3), meanf_Tmax=round(meanf_Tmax,3), gap_Tmax=round(gap_Tmax,3))])
cat("\nDONE -> tab_methodo_test1_median_vs_mean.csv + tab_methodo_test2_clhs_vs_archetype.csv\n")
