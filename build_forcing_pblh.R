# ==============================================================================
# Build a forcing file with REAL h_sbl = MERRA-2 PBLH (from get_merra2_pblh.sh),
# aligned to musica_in_Blois.nc's hourly half-hour-centred time axis. Replaces the
# fabricated constant-1000 m placeholder. MERRA-2 tavg1 stamps (HH:30) match the
# forcing stamps exactly, so the join is by datetime string.
#   Rscript build_forcing_pblh.R
# Out: in_files/musica_in_Blois_pblh.nc
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(data.table) })

src <- "in_files/musica_in_Blois.nc"; dst <- "in_files/musica_in_Blois_pblh.nc"
nc <- nc_open(src)
tu <- ncatt_get(nc, "time", "units")$value
th <- ncvar_get(nc, "time")
t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
fdt <- format(t0 + th * 3600, "%Y-%m-%dT%H:%M", tz = "UTC")   # e.g. 2021-01-01T00:30

pblh <- fread("out_files/Chapter3/merra2_pblh.csv")            # datetime, pblh
setkey(pblh, datetime)
h <- pblh[.(fdt), pblh]                                        # aligned to forcing steps
nmiss <- sum(!is.finite(h))
cat(sprintf("forcing steps: %d | matched PBLH: %d | missing: %d\n", length(fdt), sum(is.finite(h)), nmiss))
if (nmiss > 0) {                                              # fill gaps: LOCF then NOCB
  idx <- which(is.finite(h))
  h <- approx(idx, h[idx], xout = seq_along(h), method = "constant", rule = 2)$y
  cat("  gaps filled by nearest-neighbour\n")
}
cat(sprintf("h_sbl (PBLH) range: %.0f - %.0f m | median %.0f | daytime(>800m) share %.2f\n",
            min(h), max(h), median(h), mean(h > 800)))

# rebuild forcing = all original vars + h_sbl(x,y,time)
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
cat("wrote", dst, "\n")
