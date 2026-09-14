# ==============================================================================
# SMOKE TEST: does MuSICA v3.2.3 (in_files/model-3.2.3/musica) run through the
# existing run_musica_one machinery and produce a format-compatible NetCDF?
# Runs ONE plot (S0001) real full canopy, compares ΔTmax to the legacy v3.2.0
# coalition 1111 nc for the same plot. NO batch until this passes.
#   Rscript c1_smoke_v323.R
# ==============================================================================
suppressPackageStartupMessages({
  library(ncdf4); library(lubridate); library(data.table); library(dplyr)
  library(rmusica); library(musica.tools)
})
src <- list.files("R", pattern="\\.R$", full.names=TRUE); src <- src[!grepl("/(h1_|lovb_)", src)]
invisible(lapply(src, source)); source("Chapter3_config.R")

MUSICA_NEW <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
MUSICA_LEG <- "/home/corroyez/Documents/musica/musica"
cat("new binary:", MUSICA_NEW, "\nlegacy     :", MUSICA_LEG, "\n\n")

samp <- as.data.table(readRDS("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"))
samp <- samp[is.finite(LAI) & is.finite(Hmax) & is.finite(fCover)]
set.seed(42); sub <- samp[, .SD[sample(.N, min(.N, 100))], by = Cluster]; sub[, pid := sprintf("S%04d", .I)]
prow <- as.data.frame(sub[1])   # S0001 (P-? first row)
cat(sprintf("plot %s: LAI=%.2f Hmax=%.1f fCover=%.2f\n", prow$pid, prow$LAI, prow$Hmax, prow$fCover))

ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
dm <- as.data.table(extract_macro_daily(CFG_C3$forcing_file, ds))
Z_FIX <- 1.0; SHIFT <- 2L
metrics_one <- function(path) {
  if (!file.exists(path) || file.size(path) < 1000) return(NA_real_)
  nc <- try(nc_open(path), silent = TRUE); if (inherits(nc, "try-error")) return(NA_real_)
  on.exit(nc_close(nc))
  vars <- names(nc$var)
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% vars)) {
    cat("  !! missing expected vars; nc has:", paste(head(vars,20),collapse=", "), "\n"); return(NA_real_)
  }
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

OUT <- "out_files/Chapter1/nc_smoke_v323"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
nc_new <- file.path(OUT, sprintf("%s_1111_v323.nc", prow$pid))
if (file.exists(nc_new)) file.remove(nc_new)
sc <- list(lai_fn = function(p) prow$LAI, hmax_fn = function(p) prow$Hmax,
           fcover_fn = function(p) prow$fCover, lad_fn = make_lad_real, phenology_fn = NULL)

cat("\n--- running v3.2.3 (one plot) ---\n")
t0 <- proc.time()[3]
run_musica_one(prow, sc, nc_new, CFG_C3$forcing_file, MUSICA_NEW)
cat(sprintf("elapsed: %.0f s | nc exists: %s | size: %s\n", proc.time()[3]-t0,
            file.exists(nc_new), if (file.exists(nc_new)) file.size(nc_new) else NA))

dT_new <- metrics_one(nc_new)
dT_leg <- metrics_one(file.path("out_files/Chapter1/nc_shapley2x", sprintf("%s_1111.nc", prow$pid)))
cat(sprintf("\n=== RESULT ===\nΔTmax v3.2.3 = %+.3f °C\nΔTmax legacy = %+.3f °C\ndiff (new-legacy) = %+.3f °C\n",
            dT_new, dT_leg, dT_new - dT_leg))
if (is.finite(dT_new)) cat("\nSMOKE PASS: v3.2.3 runs and produces a format-compatible nc.\n") else
  cat("\nSMOKE FAIL: see messages above (run error or incompatible output format).\n")
