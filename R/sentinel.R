# ==============================================================================
# Chapter 1 — Sentinel-2 + FORMS-H annex (H3)
# Chapter 3 — Sentinel-2 LAI time series extraction and smoothing
# ==============================================================================

#' Build the S2 + FORMS-H naive scenario by extracting optical proxy inputs.
#'
#' Replaces LiDAR LAI with Sentinel-2 LAI and LiDAR Hmax with FORMS-H canopy
#' height. Uses a uniform LAD profile. Missing values are imputed by the
#' sample mean.
#'
#' @param df_sample  cLHS sample dataframe with x, y columns.
#' @param in_dir     Input directory containing the S2 and FORMS-H rasters.
#' @param agg_factor Aggregation factor (default 2).
#' @return List with elements $scenario (scenario list) and $df_sample (augmented).
make_s2_formsh_scenario <- function(df_sample, in_dir, agg_factor = 2) {
  r_s2     <- rast(file.path(in_dir, "s2lai_summer_atbd_res_10_m.tif"))
  r_formsh <- rast(file.path(in_dir, "FORMS-H_Blois.tif"))
  r_s2_20  <- terra::aggregate(r_s2, fact = agg_factor, fun = "mean", na.rm = TRUE)
  r_fh_20  <- resample(r_formsh, r_s2_20, method = "bilinear")
  pts <- vect(as.matrix(df_sample[, c("x", "y")]), crs = crs(r_s2_20))
  df_sample$S2_LAI  <- as.numeric(terra::extract(r_s2_20, pts)[, 2])
  df_sample$FORMS_H <- as.numeric(terra::extract(r_fh_20,  pts)[, 2])
  na_count <- sum(is.na(df_sample$S2_LAI) | is.na(df_sample$FORMS_H))
  if (na_count > 0) {
    df_sample$S2_LAI[is.na(df_sample$S2_LAI)]   <- mean(df_sample$S2_LAI,   na.rm = TRUE)
    df_sample$FORMS_H[is.na(df_sample$FORMS_H)] <- mean(df_sample$FORMS_H, na.rm = TRUE)
  }
  scenario <- list(
    name      = "ANNEX_S2_FORMSH",
    lai_fn    = function(plot_row) as.numeric(plot_row$S2_LAI),
    hmax_fn   = function(plot_row) as.numeric(plot_row$FORMS_H),
    fcover_fn = function(plot_row) 1.0,
    lad_fn    = make_lad_uniform
  )
  list(scenario = scenario, df_sample = df_sample)
}

# ==============================================================================
# Chapter 3 — S2 LAI time series extraction and smoothing
#
# DATA REQUIREMENT — multi-date S2 LAI stack:
#   Expected at:  in_files/s2lai_ts_blois.tif
#   Format:       SpatRaster, one layer per S2 acquisition date.
#                 Layer names must be parseable as ISO dates ("2021-03-15" etc.)
#                 or Julian day strings ("2021-074").
#   Units:        LAI in m²/m², same scale as s2lai_summer_atbd_res_10_m.tif.
#   Resolution:   10 m UTM 31N, same extent as the single-date file.
#   How to produce:
#     → Run NC_Full/02_CODES/Sentinel_2/4_create_s2_metric_time_series.R
#     → Apply Whittaker smoother (see NC_Full/02_CODES/Sentinel_2/5_smooth_s2_ts.R)
#     → Export to in_files/ as s2lai_ts_blois.tif
#   Until available, load_s2_ts() will stop with an informative error.
# ==============================================================================

#' Load a multi-date Sentinel-2 LAI raster stack from disk.
#'
#' @title Load S2 LAI time series stack
#' @description Reads the multi-layer SpatRaster expected at `s2lai_ts_blois.tif`
#'   inside `in_dir`. Layer names must be ISO dates or "YYYY-DOY" strings.
#'   Returns the raster with names coerced to `Date` objects stored as attributes.
#'
#' @param in_dir Directory containing `s2lai_ts_blois.tif` (default: `in_files/`).
#' @param filename Filename of the stack (default: `"s2lai_ts_blois.tif"`).
#' @return SpatRaster with layer names as ISO date strings and attribute `$dates`
#'   containing the corresponding `Date` vector.
load_s2_ts <- function(in_dir, filename = "s2lai_ts_blois.tif") {
  path <- file.path(in_dir, filename)
  if (!file.exists(path)) {
    stop(sprintf(
      paste0(
        "S2 LAI time series stack not found: %s\n",
        "  Produce it from NC_Full/02_CODES/Sentinel_2/ and copy to in_files/.\n",
        "  See the DATA REQUIREMENT comment in R/sentinel.R for the expected format."
      ),
      path
    ))
  }
  r <- terra::rast(path)
  raw_names <- names(r)
  # Try ISO format first ("2021-03-15"), then YYYY-DOY ("2021-074")
  dates <- suppressWarnings(as.Date(raw_names))
  if (any(is.na(dates))) {
    dates <- suppressWarnings(as.Date(raw_names, format = "%Y-%j"))
  }
  if (any(is.na(dates))) {
    stop(sprintf(
      "Cannot parse S2 stack layer names as dates (tried ISO and YYYY-DOY). Got: %s ...",
      paste(head(raw_names, 3), collapse = ", ")
    ))
  }
  names(r) <- as.character(dates)
  attr(r, "dates") <- dates
  r
}

#' Extract Sentinel-2 LAI time series for all sample plots.
#'
#' @title Extract per-plot S2 LAI time series
#' @description For each plot in `df_sample`, extracts the LAI value at each
#'   acquisition date from the multi-layer raster stack. Returns a long-format
#'   data.frame suitable for passing to `smooth_s2_ts()`. Cloud-masked pixels
#'   (NA) are retained so the smoother can handle gaps.
#'
#' @param df_sample cLHS sample dataframe with `x` and `y` columns (UTM 31N).
#' @param s2_stack  SpatRaster from `load_s2_ts()` with ISO-date layer names.
#' @return Long data.frame with columns `plot_id, date, doy, lai`.
#'   `plot_id` is `"X{round(x)}_Y{round(y)}"`, consistent with `plot_id_from_row()`.
#' @seealso [load_s2_ts()], [smooth_s2_ts()], [plot_id_from_row()]
extract_s2_ts_per_plot <- function(df_sample, s2_stack) {
  dates <- as.Date(names(s2_stack))
  pts   <- terra::vect(as.matrix(df_sample[, c("x", "y")]),
                       crs = terra::crs(s2_stack))
  vals  <- as.data.frame(terra::extract(s2_stack, pts))
  vals  <- vals[, -1, drop = FALSE]   # drop terra's ID column
  colnames(vals) <- as.character(dates)

  plot_ids <- sprintf("X%d_Y%d",
                      round(as.numeric(df_sample$x)),
                      round(as.numeric(df_sample$y)))

  long <- do.call(rbind, lapply(seq_along(dates), function(i) {
    data.frame(
      plot_id = plot_ids,
      date    = dates[i],
      doy     = as.integer(format(dates[i], "%j")),
      lai     = vals[[i]],
      stringsAsFactors = FALSE
    )
  }))
  long[order(long$plot_id, long$date), ]
}

#' Smooth S2 LAI time series to daily resolution using a cyclic GAM spline.
#'
#' @title Smooth per-plot S2 LAI time series
#' @description Fits a GAM with a cyclic cubic regression spline (period 365 days)
#'   per plot, interpolating sparse S2 acquisitions to a complete 1–365 DOY grid.
#'   NA observations (cloud contamination) are dropped before fitting.
#'   Plots with fewer than `min_obs` valid observations receive a warning and are
#'   excluded from the output.
#'
#' @param df_ts   Long data.frame from `extract_s2_ts_per_plot()` with columns
#'   `plot_id, doy, lai`. NAs allowed.
#' @param k       Number of basis functions for the cyclic spline (default 10).
#'   Reduce if fewer acquisitions are available.
#' @param min_obs Minimum number of non-NA observations required to fit (default 4).
#' @return Named list: `plot_id → data.frame(doy, lai)` with 365 rows each.
#'   Plots that fail to fit are omitted (warning issued).
#' @seealso [extract_s2_ts_per_plot()], [make_phenology_from_s2()]
smooth_s2_ts <- function(df_ts, k = 10, min_obs = 4) {
  ids <- unique(df_ts$plot_id)
  out <- lapply(ids, function(pid) {
    sub <- df_ts[df_ts$plot_id == pid & !is.na(df_ts$lai), ]
    if (nrow(sub) < min_obs) {
      warning(sprintf(
        "Plot '%s': only %d valid S2 obs (need >= %d) — excluded.",
        pid, nrow(sub), min_obs
      ))
      return(NULL)
    }
    fit <- tryCatch(
      mgcv::gam(lai ~ s(doy, bs = "cc", k = min(k, nrow(sub) - 1L)),
                data = sub, knots = list(doy = c(1, 365))),
      error = function(e) {
        warning(sprintf("Plot '%s': GAM failed (%s) — excluded.", pid, e$message))
        NULL
      }
    )
    if (is.null(fit)) return(NULL)
    data.frame(doy = 1:365,
               lai = pmax(as.numeric(predict(fit, data.frame(doy = 1:365))), 0))
  })
  names(out) <- ids
  Filter(Negate(is.null), out)
}

