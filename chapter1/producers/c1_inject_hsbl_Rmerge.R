# ==============================================================================
# Inject h_sbl (real MERRA-2 PBLH) into the true CHS41 R-MERGED station forcing
# (Tair/RH/precip = CHS41 station, secondaries = ERA5). Produces an iter-capable
# hybrid forcing distinct from the ERA5-only tool output.
#
# Backstory. `c1_inject_hsbl_station.R` sourced from
# `prep_site_forcing/MuSICA_in_CHS41-Blois_2021-station-era5.nc`, which the
# prep_site_forcing tool wrote as pure ERA5 despite the CHS41 filename (verified
# 2026-08-24: mean Tair ≡ ERA5-v2 to within float precision after time-alignment).
# `c1_build_station_forcing_Rmerge.R` bypasses that tool bug and writes the true
# graft to `out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc`. This script
# is the injector for that Rmerge file.
#
#   Rscript c1_inject_hsbl_Rmerge.R
# In : out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc
#      out_files/Chapter1/merra2_pblh.csv (datetime, pblh)
# Out: out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(data.table) })
src <- "out_files/MuSICA_in_CHS41-Blois_2021-station_Rmerge.nc"
dst <- "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc"

stopifnot(file.exists(src))
nc <- nc_open(src)
tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
fdt <- format(t0 + th * 3600, "%Y-%m-%dT%H:%M", tz = "UTC")     # HH:30 stamps

pblh <- fread("out_files/Chapter1/merra2_pblh.csv"); setkey(pblh, datetime)
h <- pblh[.(fdt), pblh]
nmiss <- sum(!is.finite(h))
cat(sprintf("forcing steps: %d | matched PBLH: %d | missing: %d\n", length(fdt), sum(is.finite(h)), nmiss))
if (nmiss > 0) {
  idx <- which(is.finite(h))
  h <- approx(idx, h[idx], xout = seq_along(h), method = "constant", rule = 2)$y
  cat("  gaps filled by nearest-neighbour\n")
}
cat(sprintf("h_sbl range: %.0f-%.0f m | median %.0f | daytime(>800m) share %.2f\n",
            min(h), max(h), median(h), mean(h > 800)))

dims <- nc$dim
newvars <- list()
for (v in nc$var) newvars[[v$name]] <- ncvar_def(v$name, v$units, v$dim,
  missval = if (is.null(v$missval)) NA else v$missval, prec = if (v$prec == "double") "double" else "float")
newvars[["h_sbl"]] <- ncvar_def("h_sbl", "m", list(dims[["x"]], dims[["y"]], dims[["time"]]),
                                missval = -9999, prec = "double")
out <- nc_create(dst, newvars)
for (v in nc$var) ncvar_put(out, v$name, ncvar_get(nc, v$name))
ncvar_put(out, "h_sbl", array(h, dim = c(1, 1, length(h))))
nc_close(out); nc_close(nc)
cat("WROTE", dst, "\n")

# ---- sanity: Tair should equal Rmerge (grafted CHS41) exactly ---------------
n1 <- nc_open(src); Ta1 <- as.numeric(ncvar_get(n1, "Tair")); nc_close(n1)
n2 <- nc_open(dst); Ta2 <- as.numeric(ncvar_get(n2, "Tair")); nc_close(n2)
cat(sprintf("sanity Tair: src mean=%.5f dst mean=%.5f | sum|diff|=%.4f (should be ~0)\n",
            mean(Ta1), mean(Ta2), sum(abs(Ta1 - Ta2))))
