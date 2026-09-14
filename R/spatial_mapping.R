# ==============================================================================
# Chapter 3 — Spatial mapping of microclimatic attenuation
#
# Produces domain-wide maps of predicted ΔTmax by applying a GAM emulator
# (fitted on HOBO validation results) to full-domain LiDAR rasters.
#
# Workflow:
#   1. fit_microclimate_emulator() — fit GAM on HOBO validation daily data
#   2. predict_dtmax_domain()      — apply GAM to full raster stack
#   3. map_correction_impact()     — difference map (corrected − reference)
# ==============================================================================

#' Fit a microclimate GAM emulator from HOBO validation results.
#'
#' @title Fit microclimate GAM emulator
#' @description Fits `ΔTmax_sim ~ s(LAI) + s(Hmax) + s(fCover) + s(Tmax_macro)`
#'   on the daily simulation × HOBO-plot data from one scenario. The emulator
#'   can then be applied to a full-domain raster stack to map predicted ΔTmax
#'   for any day without running MuSICA pixel-by-pixel.
#'
#'   Use the best-performing scenario's `val_out$daily` as training data.
#'   Structural covariates (LAI, Hmax, fCover) vary by plot; `Tmax_macro`
#'   is the ERA5 daily maximum for the target day.
#'
#' @param df_daily    Daily output from `validate_scenarios_at_hobos()$daily`,
#'   filtered to one scenario. Must have columns:
#'   `Delta_sim, LAI, Hmax, fCover, Tmax_macro`.
#'   `LAI` can be any of the six correction columns — use the one matching
#'   the scenario's `lai_fn`.
#' @param lai_col     Name of the LAI column to use as predictor (default `"LAI"`).
#' @param k           Basis dimension for each smooth (default 5; reduce if n is small).
#' @return mgcv::gam object.
#' @seealso [predict_dtmax_domain()], [map_correction_impact()]
fit_microclimate_emulator <- function(df_daily, lai_col = "LAI", k = 5) {
  stopifnot(all(c("Delta_sim", "Hmax", "fCover", "Tmax_macro") %in% names(df_daily)),
            lai_col %in% names(df_daily))
  df_daily$LAI_fit <- df_daily[[lai_col]]
  k_eff <- min(k, floor(length(unique(df_daily$LAI_fit)) / 2))
  mgcv::gam(
    Delta_sim ~ s(LAI_fit, k = k_eff) +
                s(Hmax,       k = k_eff) +
                s(fCover,     k = k_eff) +
                s(Tmax_macro, k = k_eff),
    data   = df_daily,
    method = "REML"
  )
}

#' Apply the microclimate emulator to a full-domain raster stack.
#'
#' @title Predict ΔTmax domain raster
#' @description For a specified target ERA5 day (`tmax_day`), predicts ΔTmax
#'   at every pixel using the fitted GAM emulator and LiDAR covariate rasters.
#'   Pixels with NA in any covariate layer are returned as NA.
#'
#' @param gam_model   mgcv::gam from `fit_microclimate_emulator()`.
#' @param r_lai       SpatRaster — LAI (single layer, target correction).
#' @param r_hmax      SpatRaster — Hmax (p95, m).
#' @param r_fcover    SpatRaster — fCover.
#' @param tmax_day    Numeric scalar — ERA5 Tmax for the target day (°C).
#' @param agg_factor  Aggregation factor applied before prediction (default 2).
#' @return SpatRaster of predicted ΔTmax (°C), same extent/resolution as inputs.
predict_dtmax_domain <- function(gam_model, r_lai, r_hmax, r_fcover,
                                  tmax_day, agg_factor = 2L) {
  if (agg_factor > 1) {
    r_lai    <- terra::aggregate(r_lai,    fact = agg_factor, fun = "mean", na.rm = TRUE)
    r_hmax   <- terra::aggregate(r_hmax,   fact = agg_factor, fun = "mean", na.rm = TRUE)
    r_fcover <- terra::aggregate(r_fcover, fact = agg_factor, fun = "mean", na.rm = TRUE)
  }
  # Align geometries
  r_hmax   <- terra::resample(r_hmax,   r_lai, method = "bilinear")
  r_fcover <- terra::resample(r_fcover, r_lai, method = "bilinear")

  r_stack <- c(r_lai, r_hmax, r_fcover)
  names(r_stack) <- c("LAI_fit", "Hmax", "fCover")

  df_domain <- as.data.frame(r_stack, cells = TRUE, na.rm = FALSE)
  valid     <- complete.cases(df_domain[, c("LAI_fit", "Hmax", "fCover")])
  df_domain$Tmax_macro <- tmax_day

  pred_vec <- rep(NA_real_, nrow(df_domain))
  pred_vec[valid] <- as.numeric(
    stats::predict(gam_model, newdata = df_domain[valid, ], type = "response")
  )

  r_out <- r_lai[[1]]
  terra::values(r_out) <- NA_real_
  r_out[df_domain$cell] <- pred_vec
  names(r_out) <- sprintf("Delta_Tmax_pred_%s", format(tmax_day, nsmall = 1))
  r_out
}

#' Compute a correction-impact difference map.
#'
#' @title Map LAI correction impact on predicted ΔTmax
#' @description Computes `ΔTmax(LAI_corrected) − ΔTmax(LAI_reference)` per pixel,
#'   showing where and by how much the LAI correction changes the predicted
#'   microclimate.  Positive values = warmer with correction; negative = cooler.
#'
#' @param gam_model   mgcv::gam from `fit_microclimate_emulator()`.
#' @param r_lai_ref   SpatRaster — reference LAI (e.g., LAI_S2_ATBD).
#' @param r_lai_corr  SpatRaster — corrected LAI (e.g., LAI_S2_DOPT).
#' @param r_hmax      SpatRaster — Hmax.
#' @param r_fcover    SpatRaster — fCover.
#' @param tmax_day    ERA5 Tmax for the target day (°C).
#' @param agg_factor  Aggregation factor (default 2).
#' @return SpatRaster with the signed difference in predicted ΔTmax.
map_correction_impact <- function(gam_model, r_lai_ref, r_lai_corr,
                                   r_hmax, r_fcover, tmax_day, agg_factor = 2L) {
  r_ref  <- predict_dtmax_domain(gam_model, r_lai_ref,  r_hmax, r_fcover,
                                  tmax_day, agg_factor)
  r_corr <- predict_dtmax_domain(gam_model, r_lai_corr, r_hmax, r_fcover,
                                  tmax_day, agg_factor)
  r_diff <- r_corr - r_ref
  names(r_diff) <- "Delta_Tmax_correction_impact"
  r_diff
}

#' Load full-domain LiDAR covariate rasters for Blois.
#'
#' @title Load domain rasters for spatial mapping
#' @description Loads the LAI, Hmax (p95), fCover rasters covering the full
#'   Blois domain from in_files/.
#'
#' @param in_dir  Path to in_files/.
#' @return Named list: `$lai_als, $hmax, $fcover`.
load_domain_rasters <- function(in_dir) {
  list(
    lai_als = terra::rast(file.path(in_dir, "lai_z1_res_10_m.tif")),
    hmax    = terra::rast(file.path(in_dir, "hmax_p95_res_10_m.tif")),
    fcover  = terra::rast(file.path(in_dir, "fCover_res_10_m.tif"))
  )
}

#' Identify heatwave days from ERA5 macroclimate data.
#'
#' @title Select heatwave reference days
#' @description Returns the `n` hottest days in `df_macro` as candidate reference
#'   dates for spatial mapping. Heatwave days stress-test the correction's impact.
#'
#' @param df_macro  Daily macroclimate tibble from `extract_macro_daily()`.
#' @param n         Number of days to return (default 3).
#' @return Date vector of the `n` hottest days (by `Tmax_macro`).
pick_heatwave_days <- function(df_macro, n = 3L) {
  df_macro %>%
    dplyr::arrange(dplyr::desc(Tmax_macro)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::pull(date)
}

#' Plot correction-impact map with ggplot2.
#'
#' @title Plot ΔTmax correction impact
#' @description Converts a correction-impact SpatRaster to a ggplot tile map
#'   with a diverging palette centred at 0.
#'
#' @param r_diff     SpatRaster from `map_correction_impact()`.
#' @param title      Plot title.
#' @param date_label Label for the target date (used in subtitle).
#' @return ggplot object.
plot_correction_impact <- function(r_diff, title = "LAI correction impact on ΔTmax",
                                    date_label = "") {
  df_map <- as.data.frame(r_diff, xy = TRUE)
  names(df_map)[3] <- "delta"
  lim <- max(abs(df_map$delta), na.rm = TRUE)
  ggplot(df_map, aes(x = x, y = y, fill = delta)) +
    geom_raster() +
    scale_fill_gradient2(low = "#31688e", mid = "white", high = "#d8576b",
                          midpoint = 0, limits = c(-lim, lim), na.value = "grey90",
                          name = "ΔΔTmax\n(°C)") +
    coord_equal() +
    labs(title    = title,
         subtitle = sprintf("Corrected − Reference | %s", date_label),
         x = "Easting (m)", y = "Northing (m)") +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(size = 7))
}
