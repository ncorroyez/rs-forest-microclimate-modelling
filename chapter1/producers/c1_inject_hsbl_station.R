# ==============================================================================
# Inject h_sbl (real MERRA-2 PBLH) into the TOOL station forcing → iter-capable.
# Same recipe as build_forcing_pblh.R, applied to the station forcing so v3.2.3
# yoyo/iter can run with the station meteorology. Nothing overwritten.
#   Rscript c1_inject_hsbl_station.R
# In : prep_site_forcing/MuSICA_in_CHS41-Blois_2021-station-era5.nc
#      out_files/Chapter1/merra2_pblh.csv  (datetime, pblh)
# Out: out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl.nc
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(data.table) })
src <- "prep_site_forcing/MuSICA_in_CHS41-Blois_2021-station-era5.nc"
dst <- "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl.nc"

nc <- nc_open(src)
tu <- ncatt_get(nc, "time", "units")$value; th <- ncvar_get(nc, "time")
t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
fdt <- format(t0 + th * 3600, "%Y-%m-%dT%H:%M", tz = "UTC")     # HH:30 stamps

pblh <- fread("out_files/Chapter1/merra2_pblh.csv"); setkey(pblh, datetime)
h <- pblh[.(fdt), pblh]
nmiss <- sum(!is.finite(h))
cat(sprintf("forcing steps: %d | matched PBLH: %d | missing: %d\n", length(fdt), sum(is.finite(h)), nmiss))
if (nmiss > 0) {                                               # LOCF/NOCB nearest (as build_forcing_pblh.R)
  idx <- which(is.finite(h))
  h <- approx(idx, h[idx], xout = seq_along(h), method = "constant", rule = 2)$y
  cat("  gaps filled by nearest-neighbour\n")
}
cat(sprintf("h_sbl range: %.0f-%.0f m | median %.0f | daytime(>800m) share %.2f\n",
            min(h), max(h), median(h), mean(h > 800)))

# clone all vars + add h_sbl(x,y,time), matching musica_in_Blois_pblh.nc layout
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
