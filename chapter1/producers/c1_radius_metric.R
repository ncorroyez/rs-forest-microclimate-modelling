# ==============================================================================
# Per-radius fit of MuSICA Tair @1m vs HOBO @1m (hourly): R2 / RMSE / MAE.
# Reads the radius nc (circle or square), interpolates Tair at 1 m (−2 h shift),
# merges with the 53 HOBO loggers, computes per-plot metrics, aggregates per radius.
#   Rscript c1_radius_metric.R <circle|square>
# Out: out_files/Chapter1/tables/radius_metrics_<geom>_{perplot,summary}.csv
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(lubridate); library(data.table) })
geom <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(geom)) geom <- "circle"
base <- if (geom == "square") "out_files/radius_test_square" else "out_files/radius_test_circle"
ds <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day"); Z <- 1.0; SHIFT <- 2L
ids_rm <- c("41_13","41_14","41_20","41_41","41_50","41_51","41_53")

h <- fread("in_files/Blois_data_temperature.csv")[position_sensor == "a"]
h[, time := floor_date(as.POSIXct(datetime, tz = "UTC"), "hour")]
h <- h[as.Date(time) %in% ds, .(Tair_obs = mean(t_hobo, na.rm = TRUE)), by = .(id_plot, time)]

t1m <- function(f) {
  nc <- try(nc_open(f), silent = TRUE); if (inherits(nc, "try-error")) return(NULL); on.exit(nc_close(nc))
  tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height"); vh <- median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh
  if (Z <= zl[1]) { ilo<-1L; ihi<-1L; w<-0 } else if (Z >= zl[length(zl)]) { ilo<-length(zl); ihi<-ilo; w<-0 } else {
    ilo <- max(which(zl <= Z)); ihi <- ilo + 1L; w <- (Z - zl[ilo]) / (zl[ihi] - zl[ilo]) }
  Tc <- ((1 - w) * Tk[ilo, ] + w * Tk[ihi, ]) - 273.15
  data.table(time = floor_date(t0 + dhours(th) - hours(SHIFT), "hour"), Tair_sim = Tc)[
    as.Date(time) %in% ds, .(Tair_sim = mean(Tair_sim, na.rm = TRUE)), by = time]
}

dirs <- list.dirs(base, recursive = FALSE); dirs <- dirs[grepl("m$", basename(dirs))]
res <- rbindlist(lapply(dirs, function(d) {
  r <- as.numeric(sub("m", "", basename(d)))
  rbindlist(lapply(list.files(d, "\\.nc$", full.names = TRUE), function(f) {
    pid <- sub(".*pt_(.+)_radius.*", "\\1", basename(f)); if (pid %in% ids_rm) return(NULL)
    s <- t1m(f); if (is.null(s) || nrow(s) < 24) return(NULL)
    m <- merge(s, h[id_plot == pid], by = "time"); if (nrow(m) < 24) return(NULL)
    data.table(radius = r, id_plot = pid, n = nrow(m),
               r2 = cor(m$Tair_sim, m$Tair_obs, use = "complete.obs")^2,
               rmse = sqrt(mean((m$Tair_sim - m$Tair_obs)^2, na.rm = TRUE)),
               mae = mean(abs(m$Tair_sim - m$Tair_obs), na.rm = TRUE)) }))
}))
agg <- res[, .(n_plots = .N, R2 = round(mean(r2, na.rm=TRUE),3), RMSE = round(mean(rmse, na.rm=TRUE),3),
               MAE = round(mean(mae, na.rm=TRUE),3)), by = radius][order(radius)]
cat(sprintf("=== %s : per-radius (Tair 1m, hourly, vs HOBO) ===\n", toupper(geom))); print(agg)
fwrite(res, sprintf("out_files/Chapter1/tables/radius_metrics_%s_perplot.csv", geom))
fwrite(agg, sprintf("out_files/Chapter1/tables/radius_metrics_%s_summary.csv", geom))
