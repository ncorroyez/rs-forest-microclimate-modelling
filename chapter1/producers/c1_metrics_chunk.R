# ==============================================================================
# Chapter 1 — per-point cLHS Shapley on SIX metric×period combinations, from the
# CACHED MuSICA nc (per-cluster baseline run, floor05_v2). NO MuSICA re-run: we
# re-extract from out_files/Chapter1/nc_shapley2x/<pid>_<bits>.nc.
# Metrics: ΔTmax (1 m interp, −2 h shift) | micro/macro slope (1 m interp, hourly,
#          no shift) | ΔVPDmax (1 m, reconstructed). All three at the SAME 1 m height.
#          Each over ALL period and HOT days
#          (top 10% macro Tmax). Per-cluster baseline → report mean|φ|.
# CHUNKED multi-process:  Rscript c1_metrics_chunk.R <chunk> <K>
# Writes out_files/Chapter1/tables/metrics_parts/part_<chunk>.csv (long: one row
# per plot × metric × period, with the four φ).
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE); chunk <- as.integer(args[1]); K <- as.integer(args[2])
src <- list.files("R", pattern = "\\.R$", full.names = TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))                              # CFG, build_era5_hourly, extract_macro_daily
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
Z_FIX <- 1.0; SHIFT <- 2L; P_HPA <- 1013; HOTQ <- 0.90
ds <- CFG$date_seq; ncdir <- "out_files/Chapter1/nc_shapley2x"

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N, 100))], by = Cluster]; sub[, pid := sprintf("S%04d", .I)]
my <- which(((seq_len(nrow(sub)) - 1) %% K) == (chunk - 1))
cat(sprintf("chunk %d/%d : %d plots\n", chunk, K, length(my)))

dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
hot  <- dm[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm = TRUE), date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))   # time, Tair_era5
if (chunk == 1) cat(sprintf("hot days: %d (>= %.1f C macro Tmax)\n", length(hot), quantile(dm$Tmax_macro, HOTQ, na.rm = TRUE)))

Fv <- c("LAI","Hmax","fCover","LAD"); n <- 4
coal <- as.matrix(expand.grid(LAI=0:1, Hmax=0:1, fCover=0:1, LAD=0:1))
wt_sh <- function(s) factorial(s) * factorial(n - s - 1) / factorial(n)
MET <- c("Tmax_all","Tmax_hot","slope_all","slope_hot","VPD_all","VPD_hot")

# --- 6 metrics from one cached nc (single open) -------------------------------
metrics_one <- function(path) {
  out <- setNames(rep(NA_real_, 6), MET)
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(out)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z"); wmr <- ncvar_get(nc, "wair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh; Z <- Z_FIX
  if (Z <= zl[1]) { ilo<-1L; ihi<-1L; w<-0 } else if (Z >= zl[length(zl)]) { ilo<-length(zl); ihi<-ilo; w<-0 } else {
    ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  # ΔTmax / ΔVPDmax : 1 m interp, −2 h shift, daily max, vs macro
  tvec <- t0 + dhours(th) - lubridate::hours(SHIFT)
  Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  wv <- (1 - w) * wmr[ilo, ] + w * wmr[ihi, ]
  vpd <- pmax(esat_hpa(Tc) - (wv / (1 + wv)) * P_HPA, 0) / 10
  dd <- data.table(date = as.Date(floor_date(tvec, "hour")), Tc = Tc, vpd = vpd)[
    date %in% ds, .(Tmax = max(Tc), VPDmax = max(vpd)), by = date]
  m1 <- merge(dd, dm, by = "date")
  out["Tmax_all"] <- mean(m1$Tmax - m1$Tmax_macro, na.rm = TRUE)
  out["Tmax_hot"] <- mean(m1[date %in% hot, Tmax - Tmax_macro], na.rm = TRUE)
  out["VPD_all"]  <- mean(dd$VPDmax, na.rm = TRUE)
  out["VPD_hot"]  <- mean(dd[date %in% hot, VPDmax], na.rm = TRUE)
  # micro/macro slope : SAME 1 m interp height (Tc), hourly regression vs forcing, no shift
  mic <- data.table(time = floor_date(t0 + dhours(th), "hour"), Tmic = Tc)
  mm <- merge(mic[as.Date(time) %in% ds], era5, by = "time")
  if (nrow(mm) > 10) out["slope_all"] <- as.numeric(coef(lm(Tmic ~ Tair_era5, mm))[2])
  mh <- mm[as.Date(time) %in% hot]
  if (nrow(mh) > 10) out["slope_hot"] <- as.numeric(coef(lm(Tmic ~ Tair_era5, mh))[2])
  out
}

# --- exact Shapley of a 16-value vector (named by bits) -----------------------
shap <- function(vals) {
  phi <- setNames(numeric(n), Fv)
  for (v in Fv) { others <- setdiff(Fv, v); acc <- 0
    for (m in 0:length(others)) for (S in (if (m==0) list(character(0)) else combn(others, m, simplify=FALSE))) {
      b0 <- setNames(integer(n), Fv); b0[S] <- 1L; b1 <- b0; b1[v] <- 1L
      f1 <- vals[paste(b1[Fv], collapse="")]; f0 <- vals[paste(b0[Fv], collapse="")]
      if (is.finite(f1) && is.finite(f0)) acc <- acc + wt_sh(m) * (f1 - f0) }
    phi[v] <- acc }
  phi
}

one_plot <- function(i) {
  prow <- as.data.frame(sub[i]); bits_all <- apply(coal, 1, paste, collapse="")
  V <- matrix(NA_real_, 16, 6, dimnames = list(bits_all, MET))
  for (k in 1:16) V[bits_all[k], ] <- metrics_one(file.path(ncdir, sprintf("%s_%s.nc", prow$pid, bits_all[k])))
  rbindlist(lapply(MET, function(mt) {
    phi <- shap(V[, mt])
    data.table(pid = prow$pid, Cluster = prow$Cluster, metric = mt,
               full = V["1111", mt], base = V["0000", mt],
               LAI = phi["LAI"], Hmax = phi["Hmax"], fCover = phi["fCover"], LAD = phi["LAD"])
  }))
}

pdir <- "out_files/Chapter1/tables/metrics_parts"; dir.create(pdir, recursive=TRUE, showWarnings=FALSE)
pf <- file.path(pdir, sprintf("part_%02d.csv", chunk)); if (file.exists(pf)) file.remove(pf)
nok <- 0L
for (j in seq_along(my)) {
  r <- tryCatch(one_plot(my[j]), error=function(e) NULL)
  if (!is.null(r)) { fwrite(r, pf, append = file.exists(pf)); nok <- nok + 1L }
  if (j %% 5 == 0) cat(sprintf("  chunk %d: %d/%d\n", chunk, j, length(my)))
}
cat(sprintf("chunk %d DONE (%d/%d plots)\n", chunk, nok, length(my)))
