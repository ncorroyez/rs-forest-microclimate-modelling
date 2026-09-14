# ==============================================================================
# Per-unit sensitivity for ALL SIX metric×period combos (ΔTmax, micro/macro slope,
# ΔVPDmax × all-days / hot-10%) from the EXISTING per-unit nc — NO new sims, pure
# re-extraction. Same per-unit perturbation design & step normalisation as
# c1_sensitivity_perplot_chunk*.R, applied to each metric.
#   Rscript c1_metrics6_perplot.R <chunk> <K> <ver>     ver = v320 | iter
# Out: out_files/Chapter1/tables/metrics6_<ver>/part_<chunk>.csv (long: per metric)
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table); library(dplyr) })
args <- commandArgs(trailingOnly = TRUE); chunk <- as.integer(args[1]); K <- as.integer(args[2]); VER <- args[3]
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))                                  # CFG, build_era5_hourly, extract_macro_daily
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
Z_FIX <- 1.0; SHIFT <- 2L; P_HPA <- 1013; HOTQ <- 0.90
ds <- CFG$date_seq
MET <- c("Tmax_all","Tmax_hot","slope_all","slope_hot","VPD_all","VPD_hot")
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
hot  <- dm[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm=TRUE), date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))

metrics_one <- function(path) {
  out <- setNames(rep(NA_real_, 6), MET)
  if (!file.exists(path) || file.size(path) < 1e5) return(out)
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(out); on.exit(nc_close(nc))
  if (!all(c("Tair_z","wair_z","relative_height","veget_height_top") %in% names(nc$var))) return(out)
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z"); wmr <- ncvar_get(nc,"wair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]){ilo<-1L;ihi<-1L;w<-0}else if(Z>=zl[length(zl)]){ilo<-length(zl);ihi<-ilo;w<-0}else{ilo<-max(which(zl<=Z));ihi<-ilo+1L;w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo])}
  tvec <- t0 + dhours(th) - lubridate::hours(SHIFT)
  Tc <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  wv <- (1-w)*wmr[ilo,] + w*wmr[ihi,]
  vpd <- pmax(esat_hpa(Tc) - (wv/(1+wv))*P_HPA, 0)/10
  dd <- data.table(date=as.Date(floor_date(tvec,"hour")), Tc=Tc, vpd=vpd)[date %in% ds, .(Tmax=max(Tc), VPDmax=max(vpd)), by=date]
  m1 <- merge(dd, dm, by="date")
  out["Tmax_all"]<-mean(m1$Tmax-m1$Tmax_macro,na.rm=TRUE); out["Tmax_hot"]<-mean(m1[date%in%hot,Tmax-Tmax_macro],na.rm=TRUE)
  out["VPD_all"] <-mean(dd$VPDmax,na.rm=TRUE);             out["VPD_hot"] <-mean(dd[date%in%hot,VPDmax],na.rm=TRUE)
  mic <- data.table(time=floor_date(t0+dhours(th),"hour"), Tmic=Tc); mm <- merge(mic[as.Date(time)%in%ds], era5, by="time")
  if (nrow(mm)>10) out["slope_all"]<-as.numeric(coef(lm(Tmic~Tair_era5,mm))[2])
  mh <- mm[as.Date(time)%in%hot]; if (nrow(mh)>10) out["slope_hot"]<-as.numeric(coef(lm(Tmic~Tair_era5,mh))[2])
  out
}

# nc path resolver per version
ncp <- function(pid, tag) {
  if (VER == "iter") return(file.path("out_files/Chapter1/nc_sensitivity_perplot_iter", sprintf("%s_%s.nc", pid, tag)))
  # v320: base/unif from the Shapley lattice (1111/1110), perturbations from nc_sensitivity_perplot
  if (tag == "base")    return(file.path("out_files/Chapter1/nc_shapley2x", sprintf("%s_1111.nc", pid)))
  if (tag == "unifLAD") return(file.path("out_files/Chapter1/nc_shapley2x", sprintf("%s_1110.nc", pid)))
  file.path("out_files/Chapter1/nc_sensitivity_perplot", sprintf("%s_%s.nc", pid, tag))
}

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI)&is.finite(Hmax)&is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N,min(.N,100))], by=Cluster]; sub[, pid:=sprintf("S%04d",.I)]
my <- which(((seq_len(nrow(sub))-1) %% K) == (chunk-1))
cat(sprintf("[%s] chunk %d/%d : %d plots\n", VER, chunk, K, length(my)))
DLAI<-2; DHMAX<-5; DFC<-0.1; LAI_FLOOR<-0.1; HMAX_FLOOR<-3
nrm <- function(x,s) if (is.finite(x)&&is.finite(s)&&s>1e-9) x/s else NA_real_

pdir <- sprintf("out_files/Chapter1/tables/metrics6_%s", VER); dir.create(pdir, recursive=TRUE, showWarnings=FALSE)
pf <- file.path(pdir, sprintf("part_%02d.csv", chunk)); if (file.exists(pf)) file.remove(pf)
nok <- 0L
for (j in seq_along(my)) {
  prow <- as.data.frame(sub[my[j]]); pid <- prow$pid
  r <- tryCatch({
    M <- sapply(c("base","unifLAD","LAIp","LAIm","Hmaxp","Hmaxm","fCovp","fCovm"),
                function(tg) metrics_one(ncp(pid, tg)))   # 6 metrics x 8 tags
    laim_v<-max(prow$LAI-DLAI,LAI_FLOOR); hmm_v<-max(prow$Hmax-DHMAX,HMAX_FLOOR)
    fcp_v<-min(prow$fCover+DFC,1); fcm_v<-max(prow$fCover-DFC,0.5)
    stp_lai_add<-DLAI/2; stp_lai_rem<-(prow$LAI-laim_v)/2
    stp_hm_add<-DHMAX;   stp_hm_rem<-(prow$Hmax-hmm_v)
    stp_fc_add<-(fcp_v-prow$fCover)/DFC; stp_fc_rem<-(prow$fCover-fcm_v)/DFC
    rbindlist(lapply(MET, function(mt) data.table(
      pid=pid, Cluster=prow$Cluster, ver=VER, metric=mt,
      LAI=prow$LAI, Hmax=prow$Hmax, fCover=prow$fCover, base=M[mt,"base"],
      LAI_add =nrm(M[mt,"LAIp"] -M[mt,"base"], stp_lai_add),
      LAI_rem =nrm(M[mt,"LAIm"] -M[mt,"base"], stp_lai_rem),
      Hmax_add=nrm(M[mt,"Hmaxp"]-M[mt,"base"], stp_hm_add/5),
      Hmax_rem=nrm(M[mt,"Hmaxm"]-M[mt,"base"], stp_hm_rem),
      fCov_add=nrm(M[mt,"fCovp"]-M[mt,"base"], stp_fc_add),
      fCov_rem=nrm(M[mt,"fCovm"]-M[mt,"base"], stp_fc_rem),
      dT_LAD  =M[mt,"base"]-M[mt,"unifLAD"])))
  }, error=function(e) NULL)
  if (!is.null(r)) { fwrite(r, pf, append=file.exists(pf)); nok<-nok+1L }
  if (j %% 10 == 0) cat(sprintf("  [%s] chunk %d: %d/%d\n", VER, chunk, j, length(my)))
}
cat(sprintf("[%s] chunk %d DONE (%d/%d)\n", VER, chunk, nok, length(my)))
