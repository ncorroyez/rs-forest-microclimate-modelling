# ==============================================================================
# Re-extraction of the native20 per-plot sensitivity metrics under the TIME-MATCHED
# convention (convention B), replacing the previous independent-daily-maxima one.
#
# Two changes vs c1_metrics6_perplot_frblo.R, everything else identical:
#   1. NO -2 h clock shift (the model output and the forcing share one clock; the
#      shift was an ad-hoc cross-clock alignment).
#   2. ΔTmax / ΔVPDmax are TIME-MATCHED: for each day we locate the hour of the
#      MACROCLIMATIC daily maximum and read the sub-canopy value at that same hour
#      (Bouwen 2025: "for each day, maximum air temperature at forcing height ...
#      was identified [and] we extracted the corresponding simulated vertical
#      profiles"), instead of taking each side's own daily maximum.
# The micro-macro slope is unchanged (it is already an hourly paired regression).
#
# NO new MuSICA runs: pure re-extraction from the existing 3200 NetCDFs.
# Writes to a NEW directory so the previous metrics remain intact for comparison.
#   Rscript c1_metrics6_perplot_convB.R [ncores]
# Out: out_files/Chapter1/tables/metrics6_native20_convB/part_all.csv
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table); library(parallel) })
args <- commandArgs(trailingOnly = TRUE)
NC <- if (length(args) >= 1) as.integer(args[1]) else 6L
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
Z_FIX <- 1.0; P_HPA <- 1013; HOTQ <- 0.90
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
MET <- c("Tmax_all","Tmax_hot","slope_all","slope_hot","VPD_all","VPD_hot")
FORC <- "in_files/FR-Blo_2021_v2.nc"      # same forcing as the native20 headline runs

# ---- macro, hourly, on the model clock: daily max AND the hour it occurs -------
ncf <- nc_open(FORC); Tmac <- as.numeric(ncvar_get(ncf, "Tair")) - 273.15; nc_close(ncf)
tmac <- force_utc_nc(FORC, "time")
MAC  <- data.table(th = floor_date(tmac, "hour"), Tmac = Tmac)[, date := as.Date(th)][date %in% ds]
MAC  <- MAC[, .(Tmac = mean(Tmac, na.rm = TRUE)), by = .(th, date)]
MACd <- MAC[, .(Tmax_macro = max(Tmac, na.rm = TRUE), th_max = th[which.max(Tmac)]), by = date]
hot  <- MACd[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm = TRUE), date]
cat(sprintf("macro: %d days, %d hot days\n", nrow(MACd), length(hot)))

metrics_one <- function(path) {
  out <- setNames(rep(NA_real_, 6), MET)
  if (!file.exists(path) || file.size(path) < 1e5) return(out)
  nc <- try(nc_open(path), silent=TRUE); if (inherits(nc,"try-error")) return(out); on.exit(nc_close(nc))
  if (!all(c("Tair_z","wair_z","relative_height","veget_height_top") %in% names(nc$var))) return(out)
  tu <- ncatt_get(nc,"time","units")$value; t0 <- as.POSIXct(sub("hours since ","",tu), tz="UTC")
  th <- ncvar_get(nc,"time"); Tk <- ncvar_get(nc,"Tair_z"); wmr <- ncvar_get(nc,"wair_z")
  rh <- ncvar_get(nc,"relative_height"); vh <- stats::median(ncvar_get(nc,"veget_height_top"), na.rm=TRUE)
  zl <- rh*vh; Z <- Z_FIX
  if (Z<=zl[1]){ilo<-1L;ihi<-1L;w<-0} else if(Z>=zl[length(zl)]){ilo<-length(zl);ihi<-ilo;w<-0} else {
    ilo<-max(which(zl<=Z)); ihi<-ilo+1L; w<-(Z-zl[ilo])/(zl[ihi]-zl[ilo]) }
  Tc  <- ((1-w)*Tk[ilo,] + w*Tk[ihi,]) - 273.15
  wv  <- (1-w)*wmr[ilo,] + w*wmr[ihi,]
  vpd <- pmax(esat_hpa(Tc) - (wv/(1+wv))*P_HPA, 0)/10
  M <- data.table(th = floor_date(t0 + dhours(th), "hour"), Tc = Tc, vpd = vpd)   # NO shift
  M[, date := as.Date(th)]
  # --- time-matched offsets: sub-canopy value AT the hour of the macro daily max
  B <- merge(M[date %in% ds], MACd, by = "date")[th == th_max]
  if (nrow(B) > 10) {
    out["Tmax_all"] <- mean(B$Tc  - B$Tmax_macro, na.rm=TRUE)
    out["VPD_all"]  <- mean(B$vpd, na.rm=TRUE)
    bh <- B[date %in% hot]
    if (nrow(bh) > 3) { out["Tmax_hot"] <- mean(bh$Tc - bh$Tmax_macro, na.rm=TRUE)
                        out["VPD_hot"]  <- mean(bh$vpd, na.rm=TRUE) }
  }
  # --- slope: hourly paired regression micro ~ macro (unchanged in principle)
  mm <- merge(M[date %in% ds], MAC, by = c("th","date"))
  if (nrow(mm) > 10) out["slope_all"] <- as.numeric(coef(lm(Tc ~ Tmac, mm))[2])
  mh <- mm[date %in% hot]; if (nrow(mh) > 10) out["slope_hot"] <- as.numeric(coef(lm(Tc ~ Tmac, mh))[2])
  out
}

D <- "out_files/Chapter1/nc_sensitivity_perplot_native20"
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N,100))], by = Cluster]; sub[, pid := sprintf("S%04d", .I)]
# per-cluster SDs: identical recipe to c1_sensitivity_perplot_chunk_native20.R
SDc <- sub[, .(sdLAI=sd(LAI,na.rm=TRUE), sdHmax=sd(Hmax,na.rm=TRUE), sdFcov=sd(fCover,na.rm=TRUE)), by=Cluster]
setkey(SDc, Cluster)
LAI_FLOOR <- 0.1; HMAX_FLOOR <- 3
per_sd <- function(sim, base, actual_step, SD) {
  if (!is.finite(sim)||!is.finite(base)||!is.finite(actual_step)||actual_step<1e-6||SD<1e-9) return(0)
  (sim - base) * (SD / actual_step)
}
cat(sprintf("plots: %d | cores: %d\n", nrow(sub), NC))

one_plot <- function(i) {
  p <- sub[i]; sdr <- SDc[.(p$Cluster)]
  dL <- sdr$sdLAI; dH <- sdr$sdHmax; dF <- sdr$sdFcov
  g <- function(tg) metrics_one(file.path(D, sprintf("%s_%s.nc", p$pid, tg)))
  M <- sapply(c("base","unifLAD","LAIp","LAIm","Hmaxp","Hmaxm","fCovp","fCovm"), g)
  laip_v <- p$LAI+dL; laim_v <- max(p$LAI-dL, LAI_FLOOR)
  hmp_v  <- p$Hmax+dH; hmm_v  <- max(p$Hmax-dH, HMAX_FLOOR)
  fcp_v  <- min(p$fCover+dF,1); fcm_v <- max(p$fCover-dF,0.5)
  rbindlist(lapply(MET, function(mt) data.table(
    pid=p$pid, Cluster=p$Cluster, ver="native20_convB", metric=mt,
    LAI=p$LAI, Hmax=p$Hmax, fCover=p$fCover, sdLAI=dL, sdHmax=dH, sdFcov=dF,
    base=M[mt,"base"],
    LAI_add =per_sd(M[mt,"LAIp"],  M[mt,"base"], laip_v-p$LAI,    dL),
    LAI_rem =per_sd(M[mt,"LAIm"],  M[mt,"base"], p$LAI-laim_v,    dL),
    Hmax_add=per_sd(M[mt,"Hmaxp"], M[mt,"base"], hmp_v-p$Hmax,    dH),
    Hmax_rem=per_sd(M[mt,"Hmaxm"], M[mt,"base"], p$Hmax-hmm_v,    dH),
    fCov_add=per_sd(M[mt,"fCovp"], M[mt,"base"], fcp_v-p$fCover,  dF),
    fCov_rem=per_sd(M[mt,"fCovm"], M[mt,"base"], p$fCover-fcm_v,  dF),
    dT_LAD  =M[mt,"base"]-M[mt,"unifLAD"])))
}
res <- rbindlist(mclapply(seq_len(nrow(sub)), function(i)
         tryCatch(one_plot(i), error=function(e) NULL), mc.cores = NC), fill=TRUE)
outdir <- "out_files/Chapter1/tables/metrics6_native20_convB"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
fwrite(res, file.path(outdir, "part_all.csv"))
cat(sprintf("DONE: %d rows, %d plots -> %s\n", nrow(res), uniqueN(res$pid), outdir))
