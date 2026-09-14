# ==============================================================================
# c1_radius_clhs_extract.R
#
# Per-plot buffering metrics for the cLHS radius sweep, under the canonical
# convention (R/dtmax_convention.R): dTmax time-matched to the hour of the
# macroclimatic daily maximum, no clock shift. The cLHS design carries no
# sub-canopy observation, so no observation clock offset applies and the
# forcing clock (MREF) is the only reference.
#
# The sweep asks how far the clipping radius moves the simulated buffering of a
# plot. Each radius is compared to the 20 m reference plot by plot, so the
# paired difference removes the between-plot variance and isolates the radius
# effect. Differences are then read archetype by archetype, since Chapter 1
# reads everything that way.
#
# Pure re-extraction from existing NetCDFs; no MuSICA run, no LiDAR.
#
# Inputs : out_files/radius_test_clhs/<r>m/musica_out_clhs_<nnn>.nc
#          out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc  (forcing)
#          out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds (archetype)
# Outputs: out_files/Chapter1/tables/radius_clhs_perplot.csv
#          out_files/Chapter1/tables/radius_clhs_by_archetype.csv
#   Rscript c1_radius_clhs_extract.R
# ==============================================================================

suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(parallel) })
src <- list.files("R", "\\.R$", full.names = TRUE)
src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source))
source("pipeline/00_config.R")

BASE <- "out_files/radius_test_clhs"
REF  <- 20                     # reference radius of the reported model
ds   <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
Z    <- 1
FORC <- "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"

MREF <- macro_ref(FORC, ds)
macH <- {
  nc <- nc_open(FORC); tu <- ncatt_get(nc, "time", "units")$value
  th <- ncvar_get(nc, "time"); t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"),
                  Tmac = as.numeric(ncvar_get(nc, "Tair")) - 273.15)
  nc_close(nc)
  d[as.Date(time) %in% ds][, .(Tmac = mean(Tmac, na.rm = TRUE)), by = time]
}

# ---- archetype of each cLHS plot ---------------------------------------------
# The rds carries raw k-means codes; relabel_cluster() maps them to the display
# labels P1..P4 ordered by ascending LAI (R/cluster_relabel.R).
CL <- readRDS("out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds")
ARCH <- data.table(id_plot = sprintf("clhs_%03d", seq_len(nrow(CL))),
                   P = as.character(relabel_cluster(CL$Cluster)))

# ---- one plot ----------------------------------------------------------------
# A NetCDF left truncated by an interrupted run reads as an error, not as NULL,
# and one bad file would otherwise abort the whole radius. Score defensively and
# let the plot drop out as NA.
score <- function(f) {
  m <- tryCatch(micro_hourly_at(f, Z), error = function(e) NULL)
  if (is.null(m)) return(list(dt = NA_real_, sl = NA_real_))
  mm <- merge(m, macH, by = "time")
  list(dt = delta_tmax_mean(m, MREF, ds),
       sl = if (nrow(mm) > 50) as.numeric(coef(lm(Tmic ~ Tmac, mm))[2]) else NA_real_)
}

dirs <- list.dirs(BASE, recursive = FALSE)
dirs <- dirs[grepl("m$", basename(dirs))]
# A radius left partly computed cannot be read: the sweep writes plots in index
# order, so a partial directory is a biased subset of archetypes, not a sample.
MIN_PLOTS <- 350
nfil  <- vapply(dirs, function(d) length(list.files(d, "\\.nc$")), integer(1))
drop  <- basename(dirs)[nfil < MIN_PLOTS]
if (length(drop))
  cat(sprintf("skipping partial radii (< %d plots): %s\n", MIN_PLOTS, paste(drop, collapse = ", ")))
dirs  <- dirs[nfil >= MIN_PLOTS]
radii <- sort(as.numeric(sub("m$", "", basename(dirs))))
cat(sprintf("radii used: %s\n", paste0(radii, " m", collapse = ", ")))

PP <- rbindlist(lapply(radii, function(r) {
  d  <- file.path(BASE, paste0(r, "m"))
  fs <- list.files(d, "\\.nc$", full.names = TRUE)
  cat(sprintf("  %5s m : %d files\n", r, length(fs)))
  rbindlist(mclapply(fs, function(f) {
    id <- sub("^musica_out_", "", sub("\\.nc$", "", basename(f)))
    s  <- tryCatch(score(f), error = function(e) list(dt = NA_real_, sl = NA_real_))
    data.table(radius = r, id_plot = id, dt = s$dt, sl = s$sl)
  }, mc.cores = max(1, detectCores() - 2)))
}))
PP <- merge(PP, ARCH, by = "id_plot", all.x = TRUE)
fwrite(PP, "out_files/Chapter1/tables/radius_clhs_perplot.csv")

# ---- paired difference against the 20 m reference ----------------------------
if (!REF %in% radii) stop("reference radius ", REF, " m is missing")
R0 <- PP[radius == REF, .(id_plot, dt0 = dt, sl0 = sl)]
D  <- merge(PP[radius != REF], R0, by = "id_plot")
D[, `:=`(d_dt = dt - dt0, d_sl = sl - sl0)]

pair <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(list(n = length(x), mean = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
  tt <- t.test(x)
  list(n = length(x), mean = mean(x), lo = tt$conf.int[1], hi = tt$conf.int[2], p = tt$p.value)
}

AGG <- rbindlist(lapply(sort(unique(D$radius)), function(r) {
  rbindlist(lapply(c("ALL", sort(unique(na.omit(D$P)))), function(g) {
    d <- if (g == "ALL") D[radius == r] else D[radius == r & P == g]
    a <- pair(d$d_dt); b <- pair(d$d_sl)
    data.table(radius = r, group = g, n = a$n,
               d_dtmax = round(a$mean, 3), dt_lo = round(a$lo, 3), dt_hi = round(a$hi, 3),
               dt_p = signif(a$p, 3),
               d_slope = round(b$mean, 4), sl_lo = round(b$lo, 4), sl_hi = round(b$hi, 4),
               sl_p = signif(b$p, 3))
  }))
}))
fwrite(AGG, "out_files/Chapter1/tables/radius_clhs_by_archetype.csv")

cat(sprintf("\n=== paired difference against the %g m reference (dTmax, degC) ===\n", REF))
print(AGG[, .(radius, group, n, d_dtmax, dt_lo, dt_hi, dt_p)])
cat("\n=== paired difference, coupling slope ===\n")
print(AGG[, .(radius, group, n, d_slope, sl_lo, sl_hi, sl_p)])
cat("\nwrote out_files/Chapter1/tables/radius_clhs_perplot.csv\n")
cat("wrote out_files/Chapter1/tables/radius_clhs_by_archetype.csv\n")
