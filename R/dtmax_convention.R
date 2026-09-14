# ==============================================================================
# CANONICAL ΔTmax convention (single source of truth) — Chapter 1.
#
# WHY THIS FILE EXISTS. ΔTmax used to be computed independently in several scripts,
# which had silently diverged: the simulated side applied a -2 h clock shift and the
# observed side applied none, and each side took its OWN daily maximum. Both are
# fixed here, once.
#
# CONVENTION (time-matched, "B"), following Bouwen (2025), who computes ΔTmax with
# the same model: for each day the hour of the MACROCLIMATIC daily maximum is
# located, and the sub-canopy value is read AT THAT SAME HOUR. Taking each side's
# own daily maximum instead lets the micro maximum float to any hour (a sun fleck,
# sensor noise), which understates buffering in the observations.
#
# NO clock shift is applied anywhere: model output, forcing and loggers are each
# used on their own native clock, and micro/macro are paired by timestamp.
#
# Usage:
#   src <- list.files("R","\\.R$",full.names=TRUE); invisible(lapply(src, source))
#   MREF <- macro_ref(forcing_nc, dates)              # date, Tmax_macro, t_max
#   d    <- delta_tmax(micro_hourly, MREF)            # per-date Delta
#   val  <- mean(d$Delta, na.rm = TRUE)
# ==============================================================================
suppressPackageStartupMessages({ library(ncdf4); library(data.table); library(lubridate) })

#' Macroclimatic daily reference: daily maximum AND the hour at which it occurs.
#' The forcing may be sub-hourly; timestamps are floored to the hour and the max is
#' taken over all records of the day (identical to the production recipe).
#' @return data.table(date, Tmax_macro, t_max)
macro_ref <- function(forcing_nc, dates, varname = "Tair") {
  nc <- nc_open(forcing_nc)
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub(".*since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time")
  Tm <- as.numeric(ncvar_get(nc, varname)) - 273.15
  nc_close(nc)
  d <- data.table(time = floor_date(t0 + th * 3600, "hour"), Tmac = Tm)
  d <- d[as.Date(time) %in% dates]
  d[, .(Tmax_macro = max(Tmac, na.rm = TRUE),
        t_max      = time[which.max(Tmac)]), by = .(date = as.Date(time))]
}

#' Time-matched ΔTmax, per date.
#' @param micro data.table(time, Tmic) on the same clock as the forcing, hourly.
#' @param mref  output of macro_ref().
#' @return data.table(date, Tmic, Tmax_macro, Delta)
delta_tmax <- function(micro, mref) {
  m <- as.data.table(micro)[, .(time = floor_date(time, "hour"), Tmic)]
  m <- m[, .(Tmic = mean(Tmic, na.rm = TRUE)), by = time][, date := as.Date(time)]
  j <- merge(m, mref, by = "date")[time == t_max]
  j[, .(date, Tmic, Tmax_macro, Delta = Tmic - Tmax_macro)]
}

#' Mean time-matched ΔTmax over a set of dates (NA if too few matched days).
delta_tmax_mean <- function(micro, mref, dates = NULL, min_days = 30) {
  d <- delta_tmax(micro, mref)
  if (!is.null(dates)) d <- d[date %in% dates]
  if (nrow(d) < min_days) return(NA_real_)
  mean(d$Delta, na.rm = TRUE)
}

#' Sub-canopy hourly series at a fixed height, interpolated from a MuSICA output.
#' No clock shift: the returned time is the model's own clock.
#' @return data.table(time, Tmic) or NULL
micro_hourly_at <- function(nc_path, Z = 1.0) {
  if (!file.exists(nc_path) || file.size(nc_path) < 1e5) return(NULL)
  nc <- try(nc_open(nc_path), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL); on.exit(nc_close(nc))
  if (!all(c("Tair_z","relative_height","veget_height_top") %in% names(nc$var))) return(NULL)
  tu <- ncatt_get(nc, "time", "units")$value
  t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
  th <- ncvar_get(nc, "time"); Tk <- ncvar_get(nc, "Tair_z")
  rh <- ncvar_get(nc, "relative_height")
  vh <- stats::median(ncvar_get(nc, "veget_height_top"), na.rm = TRUE)
  zl <- rh * vh
  if (Z <= zl[1]) { il <- 1L; ih <- 1L; w <- 0 }
  else if (Z >= zl[length(zl)]) { il <- length(zl); ih <- il; w <- 0 }
  else { il <- max(which(zl <= Z)); ih <- il + 1L; w <- (Z - zl[il]) / (zl[ih] - zl[il]) }
  data.table(time = floor_date(t0 + dhours(th), "hour"),
             Tmic = ((1 - w) * Tk[il, ] + w * Tk[ih, ]) - 273.15)
}

# --- CLOCK ALIGNMENT (established 2026-07-29, confirmed by construction) ------
# The station forcing was built on SOLAR time, which for France is UTC+1 (GMT+1),
# whereas the HOBO loggers record in UTC. The forcing timestamps therefore run
# exactly ONE HOUR AHEAD of the logger timestamps.
#
# This is confirmed independently by the data: hourly cross-correlation of the two
# series over JJAS peaks at a -1 h shift of the forcing (r = 0.983 against 0.973
# unshifted), and the mean diurnal maximum sits at 15 h in the forcing against 14 h
# in the loggers. Construction and measurement agree.
#
# Consequence for the time-matched convention. A MuSICA output inherits the clock of
# its forcing, so simulations are sampled at the macro-max time as-is (both are
# UTC+1). OBSERVATIONS are in UTC and must be sampled one hour earlier to hit the
# same physical instant. Ignoring this inflates the warm bias by about 0.5 C
# (+1.51 instead of +0.98 C) because the sub-canopy is then read an hour past the
# macroclimatic peak.
OBS_CLOCK_OFFSET_H <- -1L

#' Macro reference shifted onto the LOGGER clock, for scoring observations.
macro_ref_obs <- function(...) {
  m <- macro_ref(...)
  m[, t_max := t_max + lubridate::hours(OBS_CLOCK_OFFSET_H)][]
}
