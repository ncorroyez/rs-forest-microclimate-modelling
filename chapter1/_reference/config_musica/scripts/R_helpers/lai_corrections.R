# ==============================================================================
# Chapter 3 — LAI source loading, correction, and per-plot extraction
#
# Six LAI sources are prepared for use in MuSICA scenarios:
#   LAI_ALS          : LiDAR full-depth (reference)
#   LAI_ALS_DOPT     : LiDAR at d_opt layers from canopy top (Pareto criterion)
#   LAI_S2_ATBD      : S2 PROSAIL inversion, ATBD parameterisation
#   LAI_S2_DOPT      : S2 PROSAIL inversion, d_opt-optimised parameterisation
#   LAI_S2_RESCALED  : LAI_S2_ATBD rescaled so plot peak matches LAI_ALS
#   LAI_RF           : Random-forest correction of LAI_S2_ATBD using LiDAR covariates
#
# Functions are grouped as:
#   1. Raster loaders (static maps and time series assembly)
#   2. Correction computation (rescaled, RF)
#   3. Per-plot extraction (produces data.frame with all six LAI columns)
# ==============================================================================

# ---- 1. Raster loaders -------------------------------------------------------

#' Assemble a multi-date S2 LAI raster stack from individual date files.
#'
#' @title Assemble S2 LAI time series from NC_Full individual-date files
#' @description Scans `res_dir` for files matching the pattern
#'   `s2lai_YYYY-MM-DD_{distrib}_res_10_m.tif`, filters to `date_range`,
#'   stacks them, and sets layer names to ISO date strings.
#'
#' @param res_dir    Directory containing per-date S2 LAI rasters
#'   (e.g., `NC_Full/03_RESULTS/Blois/Metrics/Deciduous_Only/`).
#' @param distrib    PROSAIL parameterisation tag (e.g., `"atbd"`,
#'   `"optim_Blois"`, `"atbd_optim_common"`).
#' @param date_range Length-2 Date vector `c(start, end)`. Default covers
#'   the full 2021 growing season (Apr → Oct).
#' @return SpatRaster with one layer per available date, names as ISO dates.
#'   Returns NULL with a warning if no files are found.
assemble_s2_ts_stack <- function(res_dir, distrib,
                                  date_range = c(as.Date("2021-04-01"),
                                                 as.Date("2021-10-31"))) {
  pattern <- sprintf("^s2lai_\\d{4}-\\d{2}-\\d{2}_%s_res_10_m\\.tif$", distrib)
  files   <- list.files(res_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    warning(sprintf("No S2 LAI files found for distrib='%s' in %s", distrib, res_dir))
    return(NULL)
  }
  dates   <- as.Date(stringr::str_extract(basename(files), "\\d{4}-\\d{2}-\\d{2}"))
  keep    <- !is.na(dates) & dates >= date_range[1] & dates <= date_range[2]
  files   <- files[keep]; dates <- dates[keep]
  if (length(files) == 0) { warning("No files in requested date range."); return(NULL) }
  ord     <- order(dates)
  r       <- terra::rast(files[ord])
  names(r)<- as.character(dates[ord])
  r
}

#' Load a static (single-date) LAI raster.
#'
#' @title Load static LAI raster
#' @description Convenience wrapper around `terra::rast()` with an informative
#'   error if the file is missing. Used for LAI_ALS, LAI_ALS_DOPT, LAI_S2_DOPT.
#'
#' @param path    Full path to the .tif file.
#' @param label   Human-readable label for error messages (default: basename).
#' @return SpatRaster (single layer).
load_static_lai <- function(path, label = basename(path)) {
  if (!file.exists(path))
    stop(sprintf("LAI raster not found (%s): %s", label, path))
  terra::rast(path)
}

# ---- 2. Correction computation -----------------------------------------------

#' Compute per-pixel rescaled S2 LAI: S2 × (LAI_ALS / LAI_S2_summer).
#'
#' @title Rescale S2 LAI raster to LiDAR amplitude
#' @description Pixel-wise multiplicative rescaling so the S2 LAI map matches
#'   the LiDAR LAI amplitude on average. Both rasters must have the same extent
#'   and resolution; if not, `r_s2` is resampled to `r_als` before computing.
#'   The scale factor is clamped to [0.1, 10] to avoid outliers driving extreme
#'   corrections in sparse-canopy pixels.
#'
#' @param r_als   SpatRaster — LiDAR LAI (single layer, LAI_ALS).
#' @param r_s2    SpatRaster — S2 LAI summer snapshot (single layer, LAI_S2_ATBD).
#' @return SpatRaster of rescaled S2 LAI.
compute_lai_s2_rescaled <- function(r_als, r_s2) {
  if (!terra::compareGeom(r_als, r_s2, stopOnError = FALSE)) {
    r_s2 <- terra::resample(r_s2, r_als, method = "bilinear")
  }
  scale <- r_als / r_s2
  scale <- terra::clamp(scale, 0.1, 10)
  r_out <- r_s2 * scale
  names(r_out) <- "LAI_S2_RESCALED"
  r_out
}

#' Apply rescaling factor to a multi-date S2 LAI time series stack.
#'
#' @title Rescale S2 LAI time series to LiDAR amplitude
#' @description Applies the pixel-wise scale factor `r_als / r_s2_summer` to
#'   every layer in the time series stack. Useful for injecting the correct
#'   seasonal amplitude while retaining the S2 phenological shape.
#'
#' @param r_als       SpatRaster — LiDAR LAI (reference amplitude).
#' @param r_s2_summer SpatRaster — S2 LAI summer snapshot (same parameterisation
#'   as `r_s2_ts`, used to derive the scale factor).
#' @param r_s2_ts     Multi-layer SpatRaster — S2 LAI time series.
#' @return Multi-layer SpatRaster of rescaled time series.
rescale_s2_ts <- function(r_als, r_s2_summer, r_s2_ts) {
  if (!terra::compareGeom(r_als, r_s2_summer, stopOnError = FALSE))
    r_s2_summer <- terra::resample(r_s2_summer, r_als, method = "bilinear")
  if (!terra::compareGeom(r_als, r_s2_ts, stopOnError = FALSE))
    r_s2_ts     <- terra::resample(r_s2_ts, r_als, method = "bilinear")
  scale <- terra::clamp(r_als / r_s2_summer, 0.1, 10)
  r_s2_ts * scale
}

#' Train a Random Forest to correct S2 LAI bias using LiDAR structural covariates.
#'
#' @title Train RF LAI bias correction model
#' @description Trains an RF model that predicts LiDAR LAI from S2 LAI plus
#'   structural covariates (Hmax, fCover, VCI, LCV). Designed for a small
#'   sample (HOBO plot locations) so `ntree = 500` and `mtry = 2` are used.
#'   The model can then be applied to any plot with the same covariates.
#'
#' @param df_plots  Data.frame with one row per plot, columns:
#'   `LAI_ALS, LAI_S2_ATBD, Hmax, fCover, VCI, LCV`.
#' @param seed      Random seed (default 42).
#' @return Named list with elements `$model` (randomForest object) and
#'   `$importance` (variable importance data.frame).
train_lai_rf <- function(df_plots, seed = 42, target = "LAI_ALS",
                         features = c("LAI_S2_ATBD", "Hmax", "fCover", "VCI", "LCV")) {
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' required. Install with install.packages('randomForest').")
  required_cols <- c(target, features)
  missing <- setdiff(required_cols, names(df_plots))
  if (length(missing) > 0)
    stop("Missing columns for RF training: ", paste(missing, collapse = ", "))
  df_clean <- stats::na.omit(df_plots[, required_cols])
  if (nrow(df_clean) < 10)
    stop(sprintf("Too few complete observations for RF (%d < 10).", nrow(df_clean)))
  set.seed(seed)
  rf <- randomForest::randomForest(
    x         = df_clean[, features, drop = FALSE],
    y         = df_clean[[target]],
    ntree     = 500,
    mtry      = 2L,
    importance= TRUE
  )
  imp_df <- as.data.frame(randomForest::importance(rf))
  imp_df$variable <- rownames(imp_df)
  cat(sprintf("  RF (target=%s) trained on %d plots | OOB R² = %.3f | OOB RMSE = %.3f\n",
              target, nrow(df_clean),
              1 - rf$mse[rf$ntree] / var(df_clean[[target]]),
              sqrt(rf$mse[rf$ntree])))
  list(model = rf, importance = imp_df)
}

#' Apply RF correction to a per-plot LAI data.frame (single date or time series).
#'
#' @title Apply RF LAI correction per plot
#' @description Predicts RF-corrected LAI for each row of `df_lai`, using
#'   `LAI_S2` as the input S2 column and fixed structural covariates.
#'
#' @param df_lai    Data.frame with columns `LAI_S2` (single date or series row),
#'   `Hmax, fCover, VCI, LCV`.
#' @param rf_model  randomForest object from `train_lai_rf()$model`.
#' @param s2_col    Name of the S2 LAI column to use as RF input (default `"LAI_S2_ATBD"`).
#' @return Numeric vector of RF-predicted LAI values, same length as `nrow(df_lai)`.
predict_lai_rf <- function(df_lai, rf_model, s2_col = "LAI_S2_ATBD") {
  df_feat <- df_lai[, c(s2_col, "Hmax", "fCover", "VCI", "LCV"), drop = FALSE]
  names(df_feat)[1] <- "LAI_S2_ATBD"
  as.numeric(stats::predict(rf_model, newdata = df_feat))
}

# ---- 3. Per-plot extraction and LAI table ------------------------------------

#' Extract all six LAI correction values at plot locations.
#'
#' @title Build per-plot LAI correction table
#' @description Extracts LAI values from each raster source at the supplied
#'   plot coordinates and assembles a data.frame with all six LAI variants.
#'   RF correction is trained inline on the extracted values.
#'   Missing values are imputed by plot-mean before RF training.
#'
#' @param df_plots    Data.frame with `x, y` columns (UTM 31N) and `Hmax, fCover`.
#' @param paths       Named list of raster paths (see `make_lai_paths()`).
#' @param in_dir      Path to `in_files/` for covariates not in `paths`.
#' @param agg_factor  Aggregation factor applied to 10 m rasters (default 2 = 20 m).
#' @return `df_plots` augmented with columns:
#'   `LAI_ALS, LAI_ALS_DOPT, LAI_S2_ATBD, LAI_S2_DOPT, LAI_S2_RESCALED, LAI_RF`.
build_lai_correction_table <- function(df_plots, paths, in_dir, agg_factor = 2L) {
  # Helper: aggregate then extract
  extr <- function(r, pts) {
    if (agg_factor > 1)
      r <- terra::aggregate(r, fact = agg_factor, fun = "mean", na.rm = TRUE)
    as.numeric(terra::extract(r, pts)[, 2])
  }
  pts <- terra::vect(as.matrix(df_plots[, c("x", "y")]),
                     crs = "EPSG:32631")

  cat("  Extracting LAI_ALS...\n")
  df_plots$LAI_ALS      <- extr(load_static_lai(paths$lai_als,      "LAI_ALS"),      pts)

  cat("  Extracting LAI_ALS_DOPT...\n")
  df_plots$LAI_ALS_DOPT <- extr(load_static_lai(paths$lai_als_dopt, "LAI_ALS_DOPT"), pts)

  cat("  Extracting LAI_S2_ATBD...\n")
  df_plots$LAI_S2_ATBD  <- extr(load_static_lai(paths$lai_s2_atbd,  "LAI_S2_ATBD"),  pts)

  cat("  Extracting LAI_S2_DOPT...\n")
  df_plots$LAI_S2_DOPT  <- extr(load_static_lai(paths$lai_s2_dopt,  "LAI_S2_DOPT"),  pts)

  # Rescaled: S2_ATBD × (ALS / S2_ATBD) — computed pixel-wise then extracted
  cat("  Computing LAI_S2_RESCALED...\n")
  r_als    <- load_static_lai(paths$lai_als, "LAI_ALS")
  r_s2_ref <- load_static_lai(paths$lai_s2_atbd, "LAI_S2_ATBD")
  r_rescaled <- compute_lai_s2_rescaled(r_als, r_s2_ref)
  if (agg_factor > 1)
    r_rescaled <- terra::aggregate(r_rescaled, fact = agg_factor, fun = "mean", na.rm = TRUE)
  df_plots$LAI_S2_RESCALED <- as.numeric(terra::extract(r_rescaled, pts)[, 2])

  # Covariates for RF (VCI, LCV — extracted from in_files)
  cat("  Extracting RF covariates (VCI, LCV)...\n")
  r_vci <- terra::rast(file.path(in_dir, "vci_res_10_m.tif"))
  r_lcv <- terra::rast(file.path(in_dir, "lcv_res_10_m.tif"))
  if (agg_factor > 1) {
    r_vci <- terra::aggregate(r_vci, fact = agg_factor, fun = "mean", na.rm = TRUE)
    r_lcv <- terra::aggregate(r_lcv, fact = agg_factor, fun = "mean", na.rm = TRUE)
  }
  df_plots$VCI <- as.numeric(terra::extract(r_vci, pts)[, 2])
  df_plots$LCV <- as.numeric(terra::extract(r_lcv, pts)[, 2])

  # Impute NAs (mean) before RF
  for (col in c("LAI_ALS", "LAI_S2_ATBD", "Hmax", "fCover", "VCI", "LCV")) {
    df_plots[[col]][is.na(df_plots[[col]])] <- mean(df_plots[[col]], na.rm = TRUE)
  }

  # Train RF and predict (target = full LiDAR canopy LAI)
  cat("  Training RF correction model (target=LAI_ALS)...\n")
  rf_out <- train_lai_rf(df_plots, target = "LAI_ALS")
  df_plots$LAI_RF <- predict_lai_rf(df_plots, rf_out$model)
  attr(df_plots, "rf_model")     <- rf_out$model
  attr(df_plots, "rf_importance")<- rf_out$importance

  # Second RF targeting the d_opt sensor-consistency LAI (Ch2). ML analogue of
  # the linear LAI_S2_RESCALED_DOPT anchor: tests whether a random forest can
  # extract the d_opt-consistency signal from S2 + structure better than a
  # single-factor rescaling. Trained on plots with non-NA LAI_ALS_DOPT.
  cat("  Training RF correction model (target=LAI_ALS_DOPT)...\n")
  rf_dopt_out <- train_lai_rf(df_plots, target = "LAI_ALS_DOPT")
  df_plots$LAI_RF_DOPT <- predict_lai_rf(df_plots, rf_dopt_out$model)
  attr(df_plots, "rf_model_dopt")      <- rf_dopt_out$model
  attr(df_plots, "rf_importance_dopt") <- rf_dopt_out$importance

  df_plots
}

#' Build the named list of LAI raster paths from a Chapter 3 config list.
#'
#' @title Assemble LAI raster path list
#' @description Collects the four static LAI raster paths expected by
#'   `build_lai_correction_table()`. Pass the returned list as `paths`.
#'
#' @param CFG_C3  Chapter 3 configuration list (see `Chapter3_main.R`).
#' @return Named list: `lai_als, lai_als_dopt, lai_s2_atbd, lai_s2_dopt`.
make_lai_paths <- function(CFG_C3) {
  atbd_path <- if (!is.null(CFG_C3$nc_s2_atbd)) CFG_C3$nc_s2_atbd else
    file.path(CFG_C3$in_dir, "s2lai_summer_atbd_res_10_m.tif")
  list(
    lai_als      = file.path(CFG_C3$in_dir, "lai_z1_res_10_m.tif"),
    lai_als_dopt = file.path(CFG_C3$nc_lai_als_dopt,
                             sprintf("LAI_ALS_dopt_%s.tif", CFG_C3$dopt_variant)),
    lai_s2_atbd  = atbd_path,
    lai_s2_dopt  = file.path(CFG_C3$nc_s2_dopt,
                             sprintf("s2lai_summer_opt_%s_res_10_m.tif", CFG_C3$dopt_variant))
  )
}

# ---- 4. LOO-CV for RF evaluation --------------------------------------------

#' Leave-one-out cross-validation for the RF LAI correction.
#'
#' @title LOO-CV for RF LAI correction
#' @description For each plot, trains the RF on all other plots and predicts
#'   on the held-out plot. This gives an honest evaluation of RF performance
#'   without data leakage: the same plots cannot be used for both training and
#'   validation of the MuSICA simulations.
#'
#'   Returns a data.frame with per-plot predictions from three methods:
#'   (1) raw S2, (2) LOO-RF corrected, (3) LiDAR reference (target).
#'   Summary metrics (R², RMSE, Bias) are stored as an attribute `$metrics`.
#'
#' @param df_plots  Data.frame with columns `LAI_ALS, LAI_S2_ATBD, Hmax,
#'   fCover, VCI, LCV` (output of `build_lai_correction_table()`).
#' @param seed      Random seed (default 42).
#' @return Data.frame with columns `plot_id, LAI_ALS, LAI_S2_ATBD,
#'   LAI_RF_loo` and attribute `$metrics` (data.frame with R², RMSE, Bias
#'   for both `raw_s2` and `rf_loo`).
loo_cv_rf <- function(df_plots, seed = 42) {
  required <- c("LAI_ALS", "LAI_S2_ATBD", "Hmax", "fCover", "VCI", "LCV")
  df_clean <- stats::na.omit(df_plots[, c("plot_id", required)])
  n        <- nrow(df_clean)
  if (n < 5) stop(sprintf("Too few complete plots for LOO-CV (%d < 5).", n))

  features <- c("LAI_S2_ATBD", "Hmax", "fCover", "VCI", "LCV")
  set.seed(seed)
  preds <- numeric(n)

  cat(sprintf("  LOO-CV RF: %d folds...\n", n))
  for (i in seq_len(n)) {
    train <- df_clean[-i, features, drop = FALSE]
    y_tr  <- df_clean$LAI_ALS[-i]
    test  <- df_clean[i,  features, drop = FALSE]
    rf_i  <- randomForest::randomForest(x = train, y = y_tr, ntree = 500,
                                         mtry = 2L, importance = FALSE)
    preds[i] <- as.numeric(stats::predict(rf_i, newdata = test))
  }

  df_loo <- data.frame(
    plot_id     = df_clean$plot_id,
    LAI_ALS     = df_clean$LAI_ALS,
    LAI_S2_ATBD = df_clean$LAI_S2_ATBD,
    LAI_RF_loo  = pmax(preds, 0)
  )

  calc_metrics <- function(obs, pred) {
    data.frame(
      r2   = round(cor(obs, pred, use = "complete.obs")^2, 3),
      rmse = round(sqrt(mean((obs - pred)^2, na.rm = TRUE)), 3),
      bias = round(mean(pred - obs, na.rm = TRUE), 3)
    )
  }
  metrics <- rbind(
    cbind(method = "raw_s2", calc_metrics(df_loo$LAI_ALS, df_loo$LAI_S2_ATBD)),
    cbind(method = "rf_loo", calc_metrics(df_loo$LAI_ALS, df_loo$LAI_RF_loo))
  )
  attr(df_loo, "metrics") <- metrics
  df_loo
}

# ---- 5. Sensitivity expectation (Ch1 → Ch3 bridge) --------------------------

#' Estimate expected RMSE gain from each LAI correction via thermal sensitivity.
#'
#' @title Expected RMSE gain from LAI correction
#' @description Uses the linear relationship `∂(ΔTmax)/∂(LAI)` estimated from
#'   the HOBO validation results (slope of simulated ΔTmax vs LAI across plots)
#'   to set a theoretical upper bound on RMSE improvement from each correction.
#'   This connects Chapter 1 sensitivity analysis to Chapter 3 validation.
#'
#'   Formula: `expected_gain ≈ |sensitivity| × mean(|LAI_corrected − LAI_S2_ATBD|)`
#'
#' @param df_plots    Data.frame with LAI correction columns.
#' @param sensitivity ΔTmax per LAI unit (°C / m² m⁻²). Estimated from Chapter 1
#'   LOO scenarios or provided manually (e.g., from GAMM marginal effects).
#'   If NULL, estimated empirically from a simple linear fit on `df_plots`.
#' @return Data.frame with columns `correction, mean_lai_delta, expected_rmse_gain`.
compute_sensitivity_expectation <- function(df_plots, sensitivity = NULL) {
  lai_cols <- c("LAI_ALS", "LAI_ALS_DOPT", "LAI_S2_DOPT",
                "LAI_S2_RESCALED", "LAI_RF")
  lai_cols <- intersect(lai_cols, names(df_plots))

  if (is.null(sensitivity)) {
    # Empirical estimate from HOBO plots: slope of LAI_ALS vs LAI_S2_ATBD
    # as a proxy for ΔΔTmax per ΔLAI (rough approximation)
    fit <- stats::lm(LAI_ALS ~ LAI_S2_ATBD, data = df_plots)
    sensitivity <- NA_real_
    message("sensitivity not provided — set manually from Chapter 1 GAMM marginal effects.")
  }

  do.call(rbind, lapply(lai_cols, function(col) {
    delta <- df_plots[[col]] - df_plots$LAI_S2_ATBD
    data.frame(
      correction         = col,
      mean_lai_delta     = round(mean(delta, na.rm = TRUE), 3),
      sd_lai_delta       = round(sd(delta, na.rm = TRUE), 3),
      expected_rmse_gain = round(abs(sensitivity) * mean(abs(delta), na.rm = TRUE), 3)
    )
  }))
}

# ---- 6. Build per-plot S2 time series tables for all LAI correction types ---

#' Structure-conditioned S2 phenology emulator (RF, target = LAI_S2).
#'
#' @title RF phenology emulator conditioned on LiDAR structure + d_opt
#' @description Pools all (plot x date) S2 ATBD observations and fits a random
#'   forest predicting the daily S2 LAI from the day-of-year and static LiDAR
#'   structure predictors (d_opt LAI, Hmax, fCover, VCI, LCV). Re-predicting on
#'   `doy = 1:365` per plot yields a smooth, gap-free phenology at S2 magnitude
#'   that borrows strength across structurally similar plots and tests whether
#'   the Ch2 d_opt consistency depth explains the S2 signal. The target (S2) is
#'   independent of the LiDAR predictors, so — unlike using d_opt to predict
#'   LAI_ALS — there is no leakage. Daily variation comes from `doy` (the only
#'   non-static feature); the actual per-day S2 observation is intentionally not
#'   used, so the output is a denoised structural trajectory rather than the raw
#'   (noisy) series.
#'
#' @param stk_atbd     SpatRaster multi-date stack of S2 ATBD LAI (names = dates).
#' @param df_plots     Per-plot table with `plot_id` and the feature columns.
#' @param pts          SpatVector of plot locations (UTM 31N), aligned to df_plots.
#' @param feature_cols Static predictors (default d_opt + structure metrics).
#' @param seed         RF seed (default 42).
#' @return Named list `plot_id -> data.frame(doy, lai)`, same format as
#'   `smooth_s2_ts()`. Plots with any NA feature are dropped.
build_rf_pheno_ts <- function(stk_atbd, df_plots, pts,
                              feature_cols = c("LAI_ALS_DOPT", "Hmax",
                                               "fCover", "VCI", "LCV"),
                              seed = 42) {
  if (is.null(stk_atbd)) return(NULL)
  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' required.")
  dates <- as.Date(names(stk_atbd))
  vals  <- as.data.frame(terra::extract(stk_atbd, pts))[, -1, drop = FALSE]
  long <- do.call(rbind, lapply(seq_along(dates), function(i) {
    data.frame(plot_id = df_plots$plot_id,
               doy = as.integer(format(dates[i], "%j")),
               lai = vals[[i]], stringsAsFactors = FALSE)
  }))
  long <- merge(long, df_plots[, c("plot_id", feature_cols)], by = "plot_id")
  train <- stats::na.omit(long[, c("lai", "doy", feature_cols)])
  set.seed(seed)
  rf <- randomForest::randomForest(
    x = train[, c("doy", feature_cols), drop = FALSE], y = train$lai,
    ntree = 500, mtry = 2L)
  cat(sprintf("  RF-pheno (target=LAI_S2, +d_opt/structure) on %d plot-date rows | OOB R²=%.3f\n",
              nrow(train), 1 - rf$mse[rf$ntree] / var(train$lai)))
  ids <- unique(df_plots$plot_id)
  out <- lapply(ids, function(pid) {
    pr <- df_plots[df_plots$plot_id == pid, , drop = FALSE]
    if (any(is.na(pr[1, feature_cols]))) return(NULL)
    nd <- data.frame(doy = 1:365)
    for (fc in feature_cols) nd[[fc]] <- pr[[fc]][1]
    data.frame(doy = 1:365,
               lai = pmax(as.numeric(stats::predict(rf, nd)), 0))
  })
  names(out) <- ids
  Filter(Negate(is.null), out)
}

#' Train S2->full-LAI correction RFs on the cLHS LiDAR sample (off-HOBO).
#'
#' @title cLHS-trained LAI correction models (structure, +d_opt)
#' @description Trains two random forests on the independent cLHS LiDAR sample
#'   (closed-canopy plots), target = full LiDAR LAI, for out-of-sample transfer
#'   to the HOBO plots: one on structure only, one with the Ch2 d_opt feature.
#'   Returns NULL if the cLHS sample is unavailable.
#' @param clhs_path  Path to clhs_sample.rds (x,y,LAI,Hmax,fCover,VCI).
#' @param r_s2_path  Summer S2 ATBD LAI raster (training S2 feature).
#' @param r_dopt_path LAI_ALS_dopt raster (training d_opt feature).
#' @return list(struct, dopt) of randomForest models, or NULL.
train_clhs_rf <- function(clhs_path, r_s2_path, r_dopt_path) {
  if (!file.exists(clhs_path) || !requireNamespace("randomForest", quietly = TRUE))
    return(NULL)
  cl <- as.data.frame(readRDS(clhs_path))
  pts <- terra::vect(cbind(cl$x, cl$y), crs = "EPSG:32631")
  cl$LAI_S2_ATBD  <- as.numeric(terra::extract(terra::rast(r_s2_path),  pts)[, 2])
  cl$LAI_ALS_DOPT <- as.numeric(terra::extract(terra::rast(r_dopt_path), pts)[, 2])
  cl$LAI_ALS <- cl$LAI
  base <- c("LAI_S2_ATBD", "Hmax", "fCover", "VCI")
  d1 <- stats::na.omit(cl[, c("LAI_ALS", base)])
  d2 <- stats::na.omit(cl[, c("LAI_ALS", base, "LAI_ALS_DOPT")])
  set.seed(42); m_struct <- randomForest::randomForest(x = d1[, base], y = d1$LAI_ALS, ntree = 800)
  set.seed(42); m_dopt   <- randomForest::randomForest(x = d2[, c(base, "LAI_ALS_DOPT")], y = d2$LAI_ALS, ntree = 800)
  cat(sprintf("  cLHS RF trained: struct n=%d | +d_opt n=%d\n", nrow(d1), nrow(d2)))
  list(struct = m_struct, dopt = m_dopt)
}

#' Strictly-consistent d_opt pair: S2_DOPT phenological shape at LAI_ALS_dopt magnitude.
#'
#' @title d_opt-consistent dynamic LAI series
#' @description Takes the LAI_S2_DOPT multi-date series (the d_opt-optimised S2
#'   product), and per plot rescales it so its summer peak equals LAI_ALS_dopt —
#'   making both the LiDAR anchor and the S2 dynamics live on the same d_opt
#'   scale. Robust to series/snapshot parameterisation mismatch because the
#'   rescaling targets the series' OWN summer peak, not a snapshot column.
#'
#' @param stk_dopt SpatRaster multi-date LAI_S2_DOPT stack (names = dates).
#' @param df_plots Per-plot table with plot_id and LAI_ALS_DOPT.
#' @param pts      SpatVector of plot locations.
#' @return Named list plot_id -> data.frame(doy, lai), smooth_s2_ts format.
build_dopt_pure_ts <- function(stk_dopt, df_plots, pts) {
  if (is.null(stk_dopt)) return(NULL)
  dates <- as.Date(names(stk_dopt))
  summer <- format(dates, "%m") %in% c("06", "07", "08", "09")
  vals  <- as.data.frame(terra::extract(stk_dopt, pts))[, -1, drop = FALSE]
  long <- do.call(rbind, lapply(seq_along(dates), function(i) {
    data.frame(plot_id = df_plots$plot_id, date = dates[i],
               doy = as.integer(format(dates[i], "%j")),
               raw = vals[[i]], stringsAsFactors = FALSE)
  }))
  # per-plot summer peak of the DOPT series -> scale so peak == LAI_ALS_dopt
  peak <- tapply(long$raw[long$date %in% dates[summer]],
                 long$plot_id[long$date %in% dates[summer]], max, na.rm = TRUE)
  anchor <- setNames(df_plots$LAI_ALS_DOPT, df_plots$plot_id)
  long$lai <- pmax(long$raw * (anchor[long$plot_id] / peak[long$plot_id]), 0)
  long <- long[is.finite(long$lai) & !is.na(long$lai), ]
  long <- long[order(long$plot_id, long$date), ]
  smooth_s2_ts(long, k = min(10L, floor(length(dates) * 0.8)))
}

#' Build per-plot S2 time series tables for all LAI correction types.
#'
#' @title Build per-plot S2 LAI time series for all corrections
#' @description Assembles multi-date stacks, extracts values at plot locations,
#'   applies rescaling and RF correction per date, and returns a named list of
#'   `smooth_s2_ts()`-compatible long data.frames (one per correction).
#'
#' @param df_plots     Data.frame with `x, y, plot_id` and RF/rescale covariates
#'   (output of `build_lai_correction_table()`).
#' @param CFG_C3       Chapter 3 configuration list.
#' @return Named list: `atbd, dopt, rescaled, rf` — each a named list of
#'   per-plot smoothed time series `data.frame(doy, lai)` (output of `smooth_s2_ts()`).
build_all_s2_ts <- function(df_plots, CFG_C3) {
  nc_res_blois <- CFG_C3$nc_results_blois

  # Stacks per parameterisation
  cat("  Assembling S2 time series stacks...\n")
  stk_atbd <- assemble_s2_ts_stack(nc_res_blois, "atbd",
                                    date_range = CFG_C3$ts_date_range)
  stk_dopt <- assemble_s2_ts_stack(nc_res_blois,
                                    CFG_C3$s2_dopt_distrib,
                                    date_range = CFG_C3$ts_date_range)

  pts <- terra::vect(as.matrix(df_plots[, c("x", "y")]),
                     crs = "EPSG:32631")

  # Helper: extract long table then smooth
  extract_and_smooth <- function(stk, label, s2_col = NULL,
                                  rescale_num = NULL, rescale_den = "LAI_S2_ATBD",
                                  rf_model = NULL, rf_extra_cols = NULL) {
    if (is.null(stk)) {
      warning(sprintf("Stack '%s' is NULL — skipping.", label))
      return(NULL)
    }
    cat(sprintf("  Extracting %s time series...\n", label))
    dates <- as.Date(names(stk))
    vals  <- as.data.frame(terra::extract(stk, pts))[, -1, drop = FALSE]
    colnames(vals) <- as.character(dates)

    long <- do.call(rbind, lapply(seq_along(dates), function(i) {
      lai_vals <- vals[[i]]

      if (!is.null(rescale_num)) {
        # rescaled: apply per-plot scale factor (LiDAR anchor / S2 summer ref).
        # rescale_num = "LAI_ALS"      → anchor to full LiDAR canopy LAI.
        # rescale_num = "LAI_ALS_DOPT" → anchor to the d_opt sensor-consistency
        #   LAI (Ch2): propagates the S2-LiDAR consistency at the LiDAR date T
        #   along the S2 phenological shape.
        scale_fac <- df_plots[[rescale_num]] / df_plots[[rescale_den]]
        scale_fac <- pmin(pmax(scale_fac, 0.1), 10)
        lai_vals  <- lai_vals * scale_fac
      }
      if (!is.null(rf_model)) {
        df_feat <- data.frame(
          LAI_S2_ATBD = lai_vals,
          Hmax        = df_plots$Hmax,
          fCover      = df_plots$fCover,
          VCI         = df_plots$VCI,
          LCV         = df_plots$LCV
        )
        # extra static (non-daily) features, e.g. d_opt LAI as a per-plot offset
        for (col in rf_extra_cols) df_feat[[col]] <- df_plots[[col]]
        valid    <- complete.cases(df_feat)
        lai_pred <- rep(NA_real_, nrow(df_feat))
        if (any(valid))
          lai_pred[valid] <- as.numeric(
            stats::predict(rf_model, newdata = df_feat[valid, , drop = FALSE])
          )
        lai_vals <- lai_pred
      }
      data.frame(
        plot_id = df_plots$plot_id,
        date    = dates[i],
        doy     = as.integer(format(dates[i], "%j")),
        lai     = pmax(lai_vals, 0),
        stringsAsFactors = FALSE
      )
    }))
    long <- long[order(long$plot_id, long$date), ]
    cat(sprintf("  Smoothing %s time series...\n", label))
    smooth_s2_ts(long, k = min(10L, floor(length(dates) * 0.8)))
  }

  rf_model      <- attr(df_plots, "rf_model")
  rf_model_dopt <- attr(df_plots, "rf_model_dopt")
  # cLHS-trained transfer models (off-HOBO). Uses the same S2 ATBD snapshot
  # raster as the HOBO LAI_S2_ATBD feature, and the common d_opt raster.
  clhs_rf <- train_clhs_rf(
    clhs_path  = "out_files/Sensitivity_Analysis/clhs_sample.rds",
    r_s2_path  = CFG_C3$nc_s2_atbd,
    r_dopt_path = file.path(CFG_C3$nc_lai_als_dopt,
                            sprintf("LAI_ALS_dopt_%s.tif", CFG_C3$dopt_variant)))
  list(
    atbd          = extract_and_smooth(stk_atbd, "S2_ATBD"),
    dopt          = extract_and_smooth(stk_dopt, "S2_DOPT"),
    rescaled      = extract_and_smooth(stk_atbd, "S2_RESCALED",
                                       rescale_num = "LAI_ALS"),
    # Sensor-consistency temporal correction: same S2 (ATBD) phenological shape
    # as `rescaled`, but amplitude anchored to LAI_ALS_dopt instead of the full
    # LiDAR LAI. Single-variable contrast vs `rescaled` = the anchor (Ch2 d_opt
    # consistency at date T vs full canopy LAI). Den = LAI_S2_ATBD keeps the
    # snapshot/series parameterisation consistent.
    rescaled_dopt = extract_and_smooth(stk_atbd, "S2_RESCALED_DOPT",
                                       rescale_num = "LAI_ALS_DOPT"),
    rf            = if (!is.null(rf_model))
                 extract_and_smooth(stk_atbd, "S2_RF", rf_model = rf_model)
               else NULL,
    # ML sensor-consistency temporal correction: same S2 (ATBD) series fed
    # through the RF trained on LAI_ALS_DOPT. Contrast vs `rescaled_dopt` = ML
    # mapping vs single-factor linear rescaling, both anchored on d_opt (Ch2).
    rf_dopt       = if (!is.null(rf_model_dopt))
                 extract_and_smooth(stk_atbd, "S2_RF_DOPT", rf_model = rf_model_dopt)
               else NULL,
    # Structure-conditioned S2 phenology emulator (target = LAI_S2, features =
    # DOY + d_opt + structure). Denoised daily series at S2 magnitude; tests
    # whether d_opt/structure conditioning beats raw-smoothed S2 (DYN_S2_ATBD).
    rf_pheno        = build_rf_pheno_ts(stk_atbd, df_plots, pts),
    # Ablation: same emulator WITHOUT d_opt (DOY + structure only). Contrast vs
    # rf_pheno isolates the marginal value of the Ch2 d_opt feature.
    rf_pheno_nodopt = build_rf_pheno_ts(stk_atbd, df_plots, pts,
                                        feature_cols = c("Hmax", "fCover", "VCI", "LCV")),
    # Strictly-consistent d_opt pair: LAI_S2_DOPT phenological SHAPE rescaled
    # per-plot to the LAI_ALS_dopt summer magnitude (both sensors on d_opt scale).
    rescaled_dopt_pure = build_dopt_pure_ts(stk_dopt, df_plots, pts),
    # cLHS-trained correction applied dynamically to HOBO (out-of-sample): tests
    # whether the off-HOBO + d_opt model improves MuSICA vs HOBO-trained DYN_RF.
    rf_clhs      = if (!is.null(clhs_rf))
                 extract_and_smooth(stk_atbd, "S2_RF_CLHS", rf_model = clhs_rf$struct)
               else NULL,
    rf_clhs_dopt = if (!is.null(clhs_rf))
                 extract_and_smooth(stk_atbd, "S2_RF_CLHS_DOPT", rf_model = clhs_rf$dopt,
                                    rf_extra_cols = "LAI_ALS_DOPT")
               else NULL
  )
}

# ---- 7. Load preprocessed LAI outputs (Step 1 → Step 2 handoff) --------------

#' Load all LAI correction outputs produced by Chapter3_01_lai_corrections.R.
#'
#' @title Load Chapter 3 LAI preparation outputs
#' @description Single entry-point for `Chapter3_main.R` to consume all LAI
#'   correction artefacts written by the preprocessing script. Stops with an
#'   informative message if required files are missing.
#'
#' @param cfg  Chapter 3 configuration list (i.e. `CFG_C3`). Must contain
#'   `$lai_prep_dir`.
#' @return Named list:
#'   \describe{
#'     \item{`df_plots`}{Data.frame with 6 LAI columns + HOBO covariates.}
#'     \item{`ts_by_plot`}{Named list `atbd/dopt/rescaled/rf` of per-plot
#'       smoothed S2 time series `data.frame(doy, lai)`.}
#'     \item{`loo_cv`}{LOO-CV data.frame (or `NULL` if not produced).}
#'   }
load_lai_prep <- function(cfg) {
  prep_dir <- cfg$lai_prep_dir
  required <- c("df_plots_lai.rds", "ts_by_plot.rds")
  missing  <- required[!file.exists(file.path(prep_dir, required))]
  if (length(missing) > 0) {
    stop(
      "LAI preparation outputs not found in ", prep_dir, ":\n",
      paste(" ", missing, collapse = "\n"), "\n",
      "Run Chapter2bis_lai_corrections.R first."
    )
  }
  loo_path <- file.path(prep_dir, "rf_loo_cv.rds")
  list(
    df_plots   = readRDS(file.path(prep_dir, "df_plots_lai.rds")),
    ts_by_plot = readRDS(file.path(prep_dir, "ts_by_plot.rds")),
    loo_cv     = if (file.exists(loo_path)) readRDS(loo_path) else NULL
  )
}
