# ==============================================================================
# Per-cluster LOCAL SENSITIVITY of the three article microclimate metrics
# (ΔTmax, VPDmax, micro/macro slope; ALL / HOT days) to each structural variable.
#
# Design (validated with advisor, 2026-06-18):
#  * Operating point = each cluster's MEAN structure (LAI, Hmax, fCover) + that
#    cluster's MEAN real LAD shape (make_lad_cluster_type_factory).
#  * For X in {LAI, Hmax, fCover}: FINITE DIFFERENCE across X's FULL-SAMPLE range
#    (p5 -> p95), the OTHER two held at the cluster mean, LAD shape = cluster mean.
#    Sensitivity = metric(X_high) - metric(X_low), in native units. Because the
#    response saturates, a point derivative would be misleading; the full-range
#    finite difference evaluated at the cluster operating point captures the
#    cluster-dependent potency (e.g. LAI potent in open P1, saturated in dense P4).
#  * Hmax perturbation STRETCHES the profile (homothety, .lad_rescale preserves
#    shape and integral=LAI), it does NOT truncate.
#  * LAD term = metric(real cluster shape) - metric(uniform), the effect of the
#    cluster's vertical structuring (VCI_bar -> 0); also reported per unit VCI.
#  * fCover enters only via clumping_factor (LAD untouched), so its sweep is clean.
#
# Legacy binary (CFG$musica_cmd, md5 5307...), parametric phenology (NULL), the
# article forcing and 1 m metric extraction (ported verbatim from
# c1_metrics_chunk.R). 8 runs/cluster x 4 = 32 cached MuSICA runs.
#   Rscript c1_sensitivity_percluster.R
# Out: out_files/Chapter1/tables/tab_sensitivity_percluster.csv (+ nc cache)
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)            # callmusica, calc_phenology
})
src <- list.files("R", pattern = "\\.R$", full.names = TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))                      # CFG, run_musica_one, make_lad_*, era5/macro helpers

NCDIR <- "out_files/Chapter1/nc_sensitivity"; dir.create(NCDIR, recursive = TRUE, showWarnings = FALSE)
TABDIR <- "out_files/Chapter1/tables"; dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
ds <- CFG$date_seq

# ---- sample, cluster means, full-sample range (p5-p95) ----------------------
samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
cl_mean <- samp[, .(LAI = mean(LAI), Hmax = mean(Hmax), fCover = mean(fCover),
                    VCI = mean(VCI)), by = Cluster][order(Cluster)]
# full-sample range (p5-p95) used to NORMALIZE every sensitivity to "°C per full range"
rng <- list(
  LAI    = as.numeric(quantile(samp$LAI,    c(.05, .95), na.rm = TRUE)),
  Hmax   = as.numeric(quantile(samp$Hmax,   c(.05, .95), na.rm = TRUE)),
  fCover = as.numeric(quantile(samp$fCover, c(.05, .95), na.rm = TRUE)),
  VCI    = as.numeric(quantile(samp$VCI,    c(.05, .95), na.rm = TRUE)))
span <- sapply(rng, function(z) diff(z))
# within-cluster p10-p90 = the LOCAL sweep window at each cluster's operating point
loc <- samp[, .(LAI10 = quantile(LAI,.10), LAI90 = quantile(LAI,.90),
                Hmax10 = quantile(Hmax,.10), Hmax90 = quantile(Hmax,.90),
                fCov10 = quantile(fCover,.10), fCov90 = quantile(fCover,.90)), by = Cluster][order(Cluster)]
cat("=== full-sample range (p5,p95), stored LAI convention ===\n"); print(rng)
cat("=== full-sample span ===\n"); print(round(span,3))
cat("=== cluster means ===\n"); print(cl_mean)
cat("=== within-cluster local window (p10,p90) ===\n"); print(loc)

LAD_CLUSTER <- make_lad_cluster_type_factory(as.data.frame(samp))   # picks shape by plot_row$Cluster

# ---- metric extractor: ported verbatim from c1_metrics_chunk.R --------------
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
Z_FIX <- 1.0; SHIFT <- 2L; P_HPA <- 1013; HOTQ <- 0.90
dm   <- as.data.table(extract_macro_daily(CFG$forcing_file, ds))
hot  <- dm[Tmax_macro >= quantile(Tmax_macro, HOTQ, na.rm = TRUE), date]
era5 <- as.data.table(build_era5_hourly(CFG$forcing_file, ds))
MET  <- c("Tmax_all","Tmax_hot","slope_all","slope_hot","VPD_all","VPD_hot")
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
  mic <- data.table(time = floor_date(t0 + dhours(th), "hour"), Tmic = Tc)
  mm <- merge(mic[as.Date(time) %in% ds], era5, by = "time")
  if (nrow(mm) > 10) out["slope_all"] <- as.numeric(coef(lm(Tmic ~ Tair_era5, mm))[2])
  mh <- mm[as.Date(time) %in% hot]
  if (nrow(mh) > 10) out["slope_hot"] <- as.numeric(coef(lm(Tmic ~ Tair_era5, mh))[2])
  out
}

# ---- scenario factory + runner ----------------------------------------------
mk_scen <- function(lai, hmax, fcover, shape = "cluster") {
  lad_fn <- if (shape == "uniform") make_lad_uniform else LAD_CLUSTER
  list(lai_fn = function(pr) lai, hmax_fn = function(pr) hmax,
       fcover_fn = function(pr) fcover, lad_fn = lad_fn, phenology_fn = NULL)
}
run_point <- function(cid, tag, lai, hmax, fcover, shape = "cluster") {
  nc <- file.path(NCDIR, sprintf("C%d_%s.nc", cid, tag))
  prow <- data.frame(Cluster = cid, x = 0, y = 0)
  run_musica_one(prow, mk_scen(lai, hmax, fcover, shape), nc, CFG$forcing_file, CFG$musica_cmd)
  if (!file.exists(nc)) return(NULL)
  as.list(metrics_one(nc))
}

# ---- run grid: per cluster, baseline + uniform + 6 endpoints ----------------
rows <- list()
for (i in seq_len(nrow(cl_mean))) {
  cid <- cl_mean$Cluster[i]; b <- cl_mean[i]; lc <- loc[Cluster == cid]
  cat(sprintf("\n--- cluster %d (mean LAI=%.2f Hmax=%.1f fCover=%.2f VCI=%.2f) ---\n",
              cid, b$LAI, b$Hmax, b$fCover, b$VCI))
  pts <- list(
    baseline = list(b$LAI, b$Hmax, b$fCover, "cluster"),
    uniform  = list(b$LAI, b$Hmax, b$fCover, "uniform"),
    # GLOBAL endpoints -> intrinsic potency (cluster-invariant complement)
    LAI_lo   = list(rng$LAI[1],  b$Hmax, b$fCover, "cluster"),
    LAI_hi   = list(rng$LAI[2],  b$Hmax, b$fCover, "cluster"),
    Hmax_lo  = list(b$LAI, rng$Hmax[1],  b$fCover, "cluster"),
    Hmax_hi  = list(b$LAI, rng$Hmax[2],  b$fCover, "cluster"),
    fCov_lo  = list(b$LAI, b$Hmax, rng$fCover[1], "cluster"),
    fCov_hi  = list(b$LAI, b$Hmax, rng$fCover[2], "cluster"),
    # LOCAL within-cluster (p10,p90) endpoints -> operating-point sensitivity (saturation)
    LAIloc_lo  = list(lc$LAI10,  b$Hmax, b$fCover, "cluster"),
    LAIloc_hi  = list(lc$LAI90,  b$Hmax, b$fCover, "cluster"),
    Hmaxloc_lo = list(b$LAI, lc$Hmax10, b$fCover, "cluster"),
    Hmaxloc_hi = list(b$LAI, lc$Hmax90, b$fCover, "cluster"),
    fCovloc_lo = list(b$LAI, b$Hmax, lc$fCov10, "cluster"),
    fCovloc_hi = list(b$LAI, b$Hmax, lc$fCov90, "cluster"))
  for (tag in names(pts)) {
    p <- pts[[tag]]
    m <- run_point(cid, tag, p[[1]], p[[2]], p[[3]], p[[4]])
    if (is.null(m)) { cat(sprintf("  %-9s FAILED\n", tag)); next }
    rows[[length(rows)+1L]] <- data.table(Cluster = cid, tag = tag,
      lai = p[[1]], hmax = p[[2]], fcover = p[[3]], shape = p[[4]],
      VCI_bar = b$VCI, as.data.table(m))
    cat(sprintf("  %-9s Tmax_all=%.2f VPD_all=%.2f slope_all=%.3f\n", tag, m$Tmax_all, m$VPD_all, m$slope_all))
  }
}
R <- rbindlist(rows)
fwrite(R, file.path(TABDIR, "tab_sensitivity_points.csv"))

# ---- sensitivities -----------------------------------------------------------
# (1) LOCAL, operating-point sensitivity = HEADLINE. Local slope = Δmetric /
#     Δvariable across the within-cluster p10-p90 window, then × full-sample span
#     -> "metric response per full range of the variable" (comparable, shows
#     saturation: LAI steep in open clusters, flat in dense). LAD = local slope
#     wrt VCI [(real-uniform)/VCI_bar] × full-sample VCI span.
# (2) GLOBAL "intrinsic potency" = COMPLEMENT, cluster-invariant (proves LAI is
#     the most potent lever everywhere; only its within-cluster room-to-move is
#     local, which is why per-cluster Shapley reorders).
diff_met <- function(d, hi, lo) setNames(unlist(d[tag==hi, ..MET]) - unlist(d[tag==lo, ..MET]), MET)

mk_sens <- function(scope) rbindlist(lapply(unique(R$Cluster), function(cid) {
  d <- R[Cluster == cid]; vci <- d$VCI_bar[1]
  lad_met <- diff_met(d, "baseline", "uniform")
  if (scope == "local") {
    sl <- function(var, hi, lo, col) diff_met(d, hi, lo) / (d[tag==hi, get(col)][1] - d[tag==lo, get(col)][1]) * span[[var]]
    rbindlist(list(
      data.table(Cluster=cid, scope=scope, variable="LAI",    as.data.table(as.list(sl("LAI","LAIloc_hi","LAIloc_lo","lai")))),
      data.table(Cluster=cid, scope=scope, variable="Hmax",   as.data.table(as.list(sl("Hmax","Hmaxloc_hi","Hmaxloc_lo","hmax")))),
      data.table(Cluster=cid, scope=scope, variable="fCover", as.data.table(as.list(sl("fCover","fCovloc_hi","fCovloc_lo","fcover")))),
      data.table(Cluster=cid, scope=scope, variable="LAD",    as.data.table(as.list(lad_met / vci * span[["VCI"]])))
    ))
  } else {
    rbindlist(list(
      data.table(Cluster=cid, scope=scope, variable="LAI",    as.data.table(as.list(diff_met(d,"LAI_hi","LAI_lo")))),
      data.table(Cluster=cid, scope=scope, variable="Hmax",   as.data.table(as.list(diff_met(d,"Hmax_hi","Hmax_lo")))),
      data.table(Cluster=cid, scope=scope, variable="fCover", as.data.table(as.list(diff_met(d,"fCov_hi","fCov_lo")))),
      data.table(Cluster=cid, scope=scope, variable="LAD",    as.data.table(as.list(lad_met)))
    ))
  }
}))
sens <- rbind(mk_sens("local"), mk_sens("global"))
fwrite(sens, file.path(TABDIR, "tab_sensitivity_percluster.csv"))
cat("\n=== LOCAL ΔTmax_all sensitivity (°C per full range, operating point) — HEADLINE ===\n")
print(dcast(sens[scope=="local"], Cluster ~ variable, value.var="Tmax_all")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)])
cat("\n=== GLOBAL intrinsic potency (°C, full-range diff) — complement ===\n")
print(dcast(sens[scope=="global"], Cluster ~ variable, value.var="Tmax_all")[, lapply(.SD, function(x) if(is.numeric(x)) round(x,2) else x)])
cat("\nDONE -> tab_sensitivity_points.csv + tab_sensitivity_percluster.csv\n")
