# ==============================================================================
# Chapter 3 — S2-driven phenology constructors
# ==============================================================================

#' Build a MuSICA-compatible phenology data.frame from a smoothed S2 LAI series.
#'
#' @title Build phenology from Sentinel-2 LAI time series
#' @description Produces the phenology data.frame expected by `rmusica::callmusica()`
#'   with the following design choices:
#'
#'   - `Leaf_area_1yr` is set to the smoothed S2 LAI time series (actual observed
#'     temporal variation, m²/m²). This is the only quantity that differs from the
#'     Chapter 3 Baseline scenario, allowing a clean H2 test.
#'   - Leaf age is derived via `calc_phenology()` with the same parameters as the
#'     Baseline scenario, but with `budburst_date` aligned to the S2 green-up DOY
#'     (first DOY at which LAI exceeds 5% of peak). This keeps leaf-age physiology
#'     consistent between scenarios and prevents confounding.
#'   - The output schema (columns, nleafage=1 layout) matches `calc_phenology()`
#'     exactly, since the function is called internally and only `Leaf_area_1yr`
#'     is replaced.
#'
#' @param s2_ts      Data.frame with columns `doy` (integer 1–365) and `lai`
#'   (smoothed LAI in m²/m²). May have gaps; linear interpolation fills them.
#' @param list_year  Integer vector of years to replicate (e.g., `c(2021, 2022)`).
#' @param leaf_age_max_in       Passed to `calc_phenology` (default 0.56).
#' @param relative_age_firstmax Passed to `calc_phenology` (default 0.10).
#' @param relative_age_lastmax  Passed to `calc_phenology` (default 0.75).
#' @return Data.frame with the same column schema as `calc_phenology(nleafage=1, ...)`,
#'   with `Leaf_area_1yr` replaced by S2 LAI values.
#' @references
#'   Dufrêne & Breda (1995) for the leaf-age formulation in `calc_phenology()`.
make_phenology_from_s2 <- function(s2_ts, list_year,
                                    leaf_age_max_in        = 0.56,
                                    relative_age_firstmax  = 0.10,
                                    relative_age_lastmax   = 0.75) {
  stopifnot(all(c("doy", "lai") %in% names(s2_ts)),
            length(list_year) >= 1)

  all_doy  <- 1:365
  lai_full <- approx(s2_ts$doy, s2_ts$lai, xout = all_doy, rule = 2)$y
  lai_full <- pmax(lai_full, 0)

  lai_max <- max(lai_full, na.rm = TRUE)

  # Derive budburst DOY from S2 green-up (first day > 5% of peak LAI)
  if (lai_max < 1e-6) {
    budburst <- 115L
  } else {
    above_thresh <- which(lai_full > 0.05 * lai_max)
    budburst     <- if (length(above_thresh) == 0L) 115L else above_thresh[1L]
  }

  # Build leaf-age curve with same physiology as Baseline, aligned to S2 timing.
  # Inherits the exact column schema of calc_phenology() (preserves nleafage=1 layout).
  pheno <- calc_phenology(
    list.year             = list_year,
    nleafage              = 1L,
    budburst_date         = budburst,
    leaf_age_max_in       = leaf_age_max_in,
    relative_age_firstmax = relative_age_firstmax,
    relative_age_lastmax  = relative_age_lastmax,
    LAI_max_per_cohort    = max(lai_max, 1e-6)
  )

  # Replace Leaf_area_1yr with actual S2 LAI values (same for all years)
  for (yr in list_year) {
    idx <- pheno$year == yr
    if (sum(idx) != 365L) {
      warning(sprintf("Year %d: expected 365 rows in calc_phenology output, got %d.",
                      yr, sum(idx)))
    }
    pheno$Leaf_area_1yr[idx] <- lai_full[seq_len(sum(idx))]
  }

  pheno
}

#' Factory: build a plot-level phenology_fn from precomputed per-plot S2 series.
#'
#' @title Phenology function factory from per-plot S2 LAI series
#' @description Given a named list of smoothed per-plot S2 LAI series, returns a
#'   closure with signature `function(plot_row) → phenology data.frame` for use
#'   in Chapter 3 scenario lists (the `phenology_fn` slot). If a plot has no
#'   matching time series entry, the function returns NULL and `run_musica_one()`
#'   falls back to the parametric `calc_phenology()`.
#'
#' @param ts_by_plot  Named list: `plot_id → data.frame(doy, lai)`.
#'   Names must match `plot_id_from_row()` output: `"X{round(x)}_Y{round(y)}"`.
#' @param list_year   Integer vector of simulation years (e.g., `c(2021, 2022)`).
#' @param leaf_age_max_in       Leaf-age parameter (default 0.56, matches Baseline).
#' @param relative_age_firstmax Leaf-age parameter (default 0.10, matches Baseline).
#' @param relative_age_lastmax  Leaf-age parameter (default 0.75, matches Baseline).
#' @return Closure `f(plot_row) → data.frame | NULL`.
#' @seealso [make_phenology_from_s2()], [plot_id_from_row()]
make_phenology_fn_factory <- function(ts_by_plot, list_year,
                                       leaf_age_max_in        = 0.56,
                                       relative_age_firstmax  = 0.10,
                                       relative_age_lastmax   = 0.75) {
  force(ts_by_plot); force(list_year)
  force(leaf_age_max_in); force(relative_age_firstmax); force(relative_age_lastmax)
  function(plot_row) {
    pid <- plot_id_from_row(plot_row)
    if (!pid %in% names(ts_by_plot)) {
      return(NULL)
    }
    make_phenology_from_s2(ts_by_plot[[pid]], list_year,
                            leaf_age_max_in        = leaf_age_max_in,
                            relative_age_firstmax  = relative_age_firstmax,
                            relative_age_lastmax   = relative_age_lastmax)
  }
}

#' Derive a consistent plot identifier from a single-row plot dataframe.
#'
#' @title Derive plot ID from x/y coordinates
#' @description Builds the canonical plot ID string `"X{round(x)}_Y{round(y)}"`.
#'   Used as dictionary key for per-plot S2 time series and for matching
#'   `ts_by_plot` to sample rows. Consistent with the `sim_id` x/y suffix
#'   convention in `run_musica_scenario()`.
#'
#' @param plot_row Single-row dataframe with numeric `x` and `y` columns (UTM 31N).
#' @return Character scalar, e.g., `"X552341_Y4789012"`.
plot_id_from_row <- function(plot_row) {
  sprintf("X%d_Y%d", round(as.numeric(plot_row$x)), round(as.numeric(plot_row$y)))
}
