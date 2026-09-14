# ==============================================================================
# Per-plot LAD effect (real vertical profile vs uniform), the categorical
# counterpart to the +/- unit sensitivities (LAI/Hmax/fCover) of the per-plot
# analysis. LAD is not a continuous unit, so its marginal effect is the contrast
# real-vs-uniform at the plot's OWN real LAI/Hmax/fCover:
#     dT_LAD = ΔTmax(1111, real LAD) − ΔTmax(1110, uniform LAD)
# Both coalitions already exist in nc_shapley2x (Shapley lattice) -> NO new sims,
# just re-extraction with the same fixed-1m / −2h metric. pid mapping mirrors
# c3_shapley_chunk exactly.
#   Rscript c1_sensitivity_perplot_lad.R
# Out: out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42)
sub <- samp[, .SD[sample(.N, min(.N, 100))], by = Cluster]
sub[, pid := sprintf("S%04d", .I)]

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds))
Z_FIX <- 1.0; SHIFT <- 2L
metrics_one <- function(path) {
  if (!file.exists(path) || file.size(path) < 1000) return(NA_real_)
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NA_real_)
  on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh; Z <- Z_FIX
  if (Z <= zl[1]) { ilo <- 1L; ihi <- 1L; w <- 0 }
  else if (Z >= zl[length(zl)]) { ilo <- length(zl); ihi <- ilo; w <- 0 }
  else { ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  tvec <- t0 + dhours(th) - lubridate::hours(SHIFT); Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  dd <- data.table(date = as.Date(floor_date(tvec, "hour")), Tc = Tc)[date %in% ds, .(Tmax = max(Tc)), by = date]
  m1 <- merge(dd, dm, by = "date"); mean(m1$Tmax - m1$Tmax_macro, na.rm = TRUE)
}
ND <- "out_files/Chapter1/nc_shapley2x"
res <- rbindlist(lapply(seq_len(nrow(sub)), function(i) {
  pid <- sub$pid[i]
  real <- metrics_one(file.path(ND, sprintf("%s_1111.nc", pid)))   # real LAD
  unif <- metrics_one(file.path(ND, sprintf("%s_1110.nc", pid)))   # uniform LAD (real LAI/Hmax/fCover)
  data.table(pid = pid, Cluster = sub$Cluster[i], dT_real = real, dT_unif = unif, dT_LAD = real - unif)
}))
fwrite(res, "out_files/Chapter1/tables/tab_sensitivity_perplot_lad.csv")
res[, P := relabel_cluster(Cluster)]
cat("=== per-plot LAD effect (real − uniform LAD), ΔTmax °C ===\n")
print(res[, .(median=round(median(dT_LAD,na.rm=TRUE),3), q25=round(quantile(dT_LAD,.25,na.rm=TRUE),3),
              q75=round(quantile(dT_LAD,.75,na.rm=TRUE),3), n=sum(is.finite(dT_LAD))), by=P][order(P)])
cat("DONE -> tab_sensitivity_perplot_lad.csv\n")
