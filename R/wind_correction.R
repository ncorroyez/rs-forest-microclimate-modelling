# ==============================================================================
# Wind correction for the MuSICA forcing (Chap 1).
#
# ERA5 provides the 10 m open-field wind, but MuSICA needs the wind just above each
# canopy (h + 2 m). Under a neutral surface layer a logarithmic profile maps one to
# the other:
#     U(h+2m) = U(10m) * [ln(2 + h - d) - ln(z0)] / [ln(10) - ln(z0_ERA)]
# with zero-plane displacement d = 0.7 h, canopy roughness z0 = 0.1 h, and z0_ERA the
# ERA5 surface roughness of the Blois grid cell (forecast_surface_roughness ~ 0.44 m).
#
# The correction is applied OFFLINE, per plot, by scaling the two wind components
# (Wind_E, Wind_N) of the forcing by the factor computed at the plot's *base* canopy
# height. In the sensitivity analysis the factor is held FIXED at the base height, so
# perturbing H_max does not also change the wind (no height<->wind entanglement).
# ==============================================================================
suppressPackageStartupMessages(library(ncdf4))

# EXTERNAL INPUT, not derivable from anything in this repository (noted 2026-07-31).
# Retrieved once from the Copernicus Climate Data Store: ERA5-Land / ERA5 single levels,
# variable `forecast_surface_roughness`, grid cell 47.50 N 1.25 E (Blois), hourly,
# 1 June to 30 September 2021, then averaged. The series is near-constant over the
# window (0.438 to 0.443 m), which is why a single scalar is used. No fsr file is kept
# in `in_files/`, so this literal cannot be recomputed here; treat it as a documented
# constant. Its influence is bounded rather than assumed: c1_wind_profile_correction.R
# sweeps z0,ERA over 0.01 (grassland) to 1 m (closed forest) and Fig. A2 reports the
# resulting spread in the correction factor.
Z0_ERA_BLOIS <- 0.44   # ERA5 fsr, Blois grid cell (47.50N,1.25E), Jun-Sep 2021 mean

#' Neutral log-profile wind factor U(h+2m) / U(10m)
#'
#' Maps the ERA5 open-field 10 m wind onto the wind just above a canopy of height h,
#' under a neutral surface layer, with zero-plane displacement d = 0.7h and canopy
#' roughness z0 = 0.1h. Over the 13 to 33 m canopies of this study the factor runs
#' about 0.41 to 0.48, i.e. the raw ERA5 wind is roughly 2.1 to 2.4 times too strong.
#'
#' @param h Canopy height in m (scalar or vector).
#' @param z0ERA ERA5 surface roughness of the grid cell, m. Defaults to the documented
#'   Blois constant `Z0_ERA_BLOIS`; `c1_wind_profile_correction.R` sweeps it to bound
#'   its influence.
#' @return Multiplicative factor to apply to the 10 m wind, same length as `h`.
wind_factor <- function(h, z0ERA = Z0_ERA_BLOIS) {
  d <- 0.7 * h; z0 <- 0.1 * h
  (log(2 + h - d) - log(z0)) / (log(10) - log(z0ERA))
}

#' Forcing file whose wind is scaled to a plot's canopy height
#'
#' Copies the base forcing and multiplies both wind components by `wind_factor(hmax)`.
#' Results are cached by the factor rounded to 2 dp, so about 40 files cover all 400
#' plots. The returned path is PROJECT-RELATIVE on purpose: MuSICA's `setup_musica`
#' stages the forcing with `ln -s ../<forcing>`, which an absolute path would break.
#'
#' NOT thread-safe. Two plots whose factors round to the same 2 dp target the same file,
#' so calling this concurrently (for instance inside `mclapply`) can have one worker
#' copying a file another is reading. Build every forcing serially before any fan-out;
#' `scripts/c1_blh_insensitivity.R` shows the pattern.
#'
#' @param hmax Plot canopy height in m; the factor is evaluated here and held fixed, so
#'   perturbing Hmax elsewhere does not also change the wind.
#' @param forc_base Path to the base forcing NetCDF.
#' @param cache_dir Directory for the scaled copies.
#' @param z0ERA ERA5 surface roughness, m; passed through to `wind_factor()`.
#' @return Project-relative path to the wind-corrected forcing NetCDF.
windcorr_forcing <- function(hmax, forc_base = "out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc",
                             cache_dir = "out_files/windcorr_forc", z0ERA = Z0_ERA_BLOIS) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  f  <- wind_factor(hmax, z0ERA)
  ff <- file.path(cache_dir, sprintf("forc_f%.2f.nc", round(f, 2)))
  if (!file.exists(ff) || file.size(ff) < 1e5) {
    file.copy(forc_base, ff, overwrite = TRUE)
    nc <- nc_open(ff, write = TRUE)
    for (v in c("Wind_E", "Wind_N")) ncvar_put(nc, v, ncvar_get(nc, v) * f)
    nc_close(nc)
  }
  ff
}
