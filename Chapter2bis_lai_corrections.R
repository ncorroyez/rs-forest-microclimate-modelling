# ==============================================================================
# Chapter 3 — Step 1: LAI correction preprocessing
#
# Standalone script. Run ONCE before Chapter3_main.R.
# Produces all LAI-corrected inputs for the MuSICA scenarios and diagnostic
# figures comparing the six correction methods.
#
# Inputs (read-only):
#   in_files/               — LiDAR rasters, covariates
#   NC_Full/03_RESULTS/...  — per-date S2 LAI rasters (4 parameterisations)
#   NC_Full/revision/...    — LAI_ALS_dopt, LAI_S2_dopt rasters
#
# Outputs → out_files/Chapter3/lai_prep/:
#   df_plots_lai.rds        — per-plot data.frame, 6 LAI columns + RF model
#   ts_by_plot.rds          — smoothed S2 time series (atbd/dopt/rescaled/rf)
#   rf_loo_cv.rds           — LOO-CV results with metrics
#   sensitivity_table.csv   — theoretical RMSE-gain upper bounds per correction
#
# Diagnostic figures → out_files/Chapter3/figs/:
#   c3_lai_distribution.png  — density plots of 6 LAI sources at HOBO plots
#   c3_rf_loo_cv.png         — LOO-CV scatter (RF vs raw S2 vs LiDAR)
#   c3_phenology_sample.png  — seasonal LAI curves for 4 corrections × 4 plots
#   c3_lai_correction_map.png — spatial maps of LAI correction delta rasters
# ==============================================================================

# ---- 0. Setup ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(mgcv)
  library(stringr); library(purrr); library(scales)
  library(rmusica)
  library(musica.tools)
})

invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
source("Chapter3_config.R")

FLAGS_PREP <- list(
  RELOAD_LAI   = TRUE,   # recompute LAI extraction even if cache exists
  RELOAD_TS    = TRUE,   # recompute S2 time series smoothing
  RUN_LOO_CV   = TRUE,   # LOO-CV for RF (~2 min for 43 plots)
  SAVE_FIGURES = TRUE
)

for (d in c(CFG_C3$out_dir, CFG_C3$lai_prep_dir,
            CFG_C3$figs_dir, CFG_C3$tables_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("=================================================================\n")
cat("  Chapter 3 — Step 1: LAI corrections\n")
cat("=================================================================\n\n")

# ---- 1. Load HOBO plot locations + LiDAR base covariates --------------------

cat("-- 1. Loading LiDAR rasters\n")
r_lidar <- load_lidar_rasters(CFG_C3$in_dir, agg_factor = CFG_C3$agg_factor)
df_plots <- build_hobo_inputs(CFG_C3$hobo_geojson, r_lidar$stack, CFG_C3$ids_to_remove)
df_plots$plot_id <- sprintf("X%d_Y%d", round(df_plots$x), round(df_plots$y))
cat(sprintf("   %d HOBO plots after exclusions.\n\n", nrow(df_plots)))

# ---- 2. Build the LAI correction table (6 sources) --------------------------

lai_rds <- file.path(CFG_C3$lai_prep_dir, "df_plots_lai.rds")

if (!FLAGS_PREP$RELOAD_LAI && file.exists(lai_rds)) {
  df_plots <- readRDS(lai_rds)
  cat(sprintf("-- 2. LAI correction table loaded from cache (%d plots)\n\n",
              nrow(df_plots)))
} else {
  cat("-- 2. Building LAI correction table\n")
  cat("   Sources: ALS | ALS_DOPT | S2_ATBD | S2_DOPT | S2_RESCALED | RF\n")
  df_plots <- build_lai_correction_table(
    df_plots   = df_plots,
    paths      = make_lai_paths(CFG_C3),
    in_dir     = CFG_C3$in_dir,
    agg_factor = CFG_C3$agg_factor
  )
  saveRDS(df_plots, lai_rds)
  cat(sprintf("\n   Saved → %s\n\n", lai_rds))
}

# Summary table
lai_cols <- c("LAI_ALS", "LAI_ALS_DOPT", "LAI_S2_ATBD",
              "LAI_S2_DOPT", "LAI_S2_RESCALED", "LAI_RF")
lai_cols <- intersect(lai_cols, names(df_plots))
df_summary <- do.call(rbind, lapply(lai_cols, function(col) {
  x <- df_plots[[col]]
  data.frame(source = col,
             mean = round(mean(x, na.rm=T), 2), sd = round(sd(x, na.rm=T), 2),
             min  = round(min(x, na.rm=T), 2),  max = round(max(x, na.rm=T), 2))
}))
cat("   LAI statistics across HOBO plots:\n")
print(df_summary, row.names = FALSE)
cat("\n")

# ---- 3. LOO-CV for the RF correction ----------------------------------------
# Run before MuSICA to expose RF performance without data leakage.

loo_rds <- file.path(CFG_C3$lai_prep_dir, "rf_loo_cv.rds")

if (FLAGS_PREP$RUN_LOO_CV && "LAI_RF" %in% names(df_plots)) {
  cat("-- 3. LOO-CV for RF correction\n")
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    cat("   [SKIP] 'randomForest' not installed\n\n")
    df_loo <- NULL
  } else {
    df_loo <- loo_cv_rf(df_plots)
    saveRDS(df_loo, loo_rds)
    cat("\n   LOO-CV metrics:\n")
    print(attr(df_loo, "metrics"))
    cat(sprintf("   Saved → %s\n\n", loo_rds))
  }
} else {
  df_loo <- if (file.exists(loo_rds)) readRDS(loo_rds) else NULL
  if (!is.null(df_loo)) cat("-- 3. LOO-CV loaded from cache\n\n")
}

# ---- 4. Sensitivity expectation (Ch2 → Ch3 bridge) -------------------------
# Theoretical upper bound on RMSE gain from each correction.
# Set CFG_C3$ch1_sensitivity from Chapter 1 GAMM marginal effects.

cat("-- 4. Sensitivity expectation table\n")
df_sens <- compute_sensitivity_expectation(df_plots, sensitivity = CFG_C3$ch1_sensitivity)
cat("   Mean ΔLAI and expected RMSE gain per correction:\n")
print(df_sens, row.names = FALSE)
write.csv(df_sens,
          file.path(CFG_C3$tables_dir, "c3_sensitivity_expectation.csv"),
          row.names = FALSE)
cat(sprintf("   Saved → %s\n\n",
            file.path(CFG_C3$tables_dir, "c3_sensitivity_expectation.csv")))

# ---- 5. S2 LAI time series (all 4 correction variants) ----------------------

ts_rds <- file.path(CFG_C3$lai_prep_dir, "ts_by_plot.rds")

if (!FLAGS_PREP$RELOAD_TS && file.exists(ts_rds)) {
  ts_list <- readRDS(ts_rds)
  cat(sprintf("-- 5. S2 time series loaded from cache (%s)\n\n",
              paste(names(ts_list),
                    sapply(ts_list, function(x) if (is.null(x)) "NULL" else length(x)),
                    sep = "=", collapse = ", ")))
} else {
  cat("-- 5. Building S2 LAI time series (4 corrections)\n")
  cat("   ATBD | DOPT (optim_Blois) | RESCALED | RF\n")
  ts_list <- build_all_s2_ts(df_plots, CFG_C3)
  saveRDS(ts_list, ts_rds)
  n_per <- sapply(ts_list, function(x) if (is.null(x)) 0L else length(x))
  cat(sprintf("\n   Plots smoothed: %s\n",
              paste(names(n_per), n_per, sep = "=", collapse = " | ")))
  cat(sprintf("   Saved → %s\n\n", ts_rds))
}

# ---- 6. Diagnostic figures ---------------------------------------------------

cat("-- 6. Generating diagnostic figures\n")

# 6a. LAI distribution by correction method
df_lai_long <- tidyr::pivot_longer(df_plots[, lai_cols, drop = FALSE],
                                    everything(), names_to = "source", values_to = "LAI")
p_dist <- ggplot(df_lai_long, aes(x = LAI, fill = source, colour = source)) +
  geom_density(alpha = 0.28, linewidth = 0.7) +
  geom_rug(alpha = 0.4, linewidth = 0.4) +
  scale_fill_viridis_d(option = "D") +
  scale_colour_viridis_d(option = "D") +
  labs(title    = "Chapter 3 — LAI correction comparison at HOBO plots",
       subtitle = sprintf("n = %d plots, Blois 2021", nrow(df_plots)),
       x = expression(LAI ~ (m^2 ~ m^{-2})), y = "Density",
       fill = "Source", colour = "Source") +
  theme_bw(base_size = 11)

# 6b. LOO-CV scatter
if (!is.null(df_loo)) {
  metrics_loo <- attr(df_loo, "metrics")
  lab_s2 <- with(metrics_loo[metrics_loo$method == "raw_s2", ],
                  sprintf("S2 ATBD (R²=%.2f, RMSE=%.2f)", r2, rmse))
  lab_rf <- with(metrics_loo[metrics_loo$method == "rf_loo", ],
                  sprintf("RF LOO-CV (R²=%.2f, RMSE=%.2f)", r2, rmse))
  p_loo <- ggplot(
    tidyr::pivot_longer(df_loo, c(LAI_S2_ATBD, LAI_RF_loo), names_to = "method", values_to = "pred"),
    aes(x = LAI_ALS, y = pred, colour = method)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.65, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, alpha = 0.15) +
    scale_colour_manual(values = c("LAI_S2_ATBD" = "#31688e", "LAI_RF_loo" = "#d8576b"),
                        labels = c("LAI_S2_ATBD" = lab_s2, "LAI_RF_loo" = lab_rf)) +
    labs(title    = "Chapter 3 — RF correction: LOO-CV vs raw S2",
         subtitle = "Each plot predicted by a model NOT trained on it",
         x = expression(LAI[LiDAR] ~ (m^2 ~ m^{-2})),
         y = expression(LAI[predicted] ~ (m^2 ~ m^{-2})),
         colour = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
} else {
  p_loo <- ggplot() + labs(title = "LOO-CV not available") + theme_bw()
}

# 6c. Phenology comparison: S2 corrections vs parametric, 4 sample plots
n_show   <- min(4L, nrow(df_plots))
ids_show <- df_plots$plot_id[seq_len(n_show)]
pal_pheno <- c("Parametric_ALS" = "black", "S2_ATBD" = "#31688e",
               "S2_DOPT" = "#35b779", "S2_RESCALED" = "#fde725", "S2_RF" = "#d8576b")

df_pheno_long <- do.call(rbind, lapply(ids_show, function(pid) {
  pr <- df_plots[df_plots$plot_id == pid, ]
  rows <- lapply(c("atbd", "dopt", "rescaled", "rf"), function(src) {
    if (!is.null(ts_list[[src]]) && pid %in% names(ts_list[[src]])) {
      ts <- ts_list[[src]][[pid]]
      data.frame(doy = ts$doy, lai = ts$lai,
                 source = paste0("S2_", toupper(src)), plot_id = pid)
    }
  })
  if (nrow(pr) > 0) {
    ph <- calc_phenology(list.year = CFG_C3$list_year, nleafage = 1L,
                          budburst_date = 115L, leaf_age_max_in = 0.56,
                          relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
                          LAI_max_per_cohort = pr$LAI_ALS[1])
    rows$param <- data.frame(
      doy = ph$Julian_day[ph$year == CFG_C3$list_year[1]],
      lai = ph$Leaf_area_1yr[ph$year == CFG_C3$list_year[1]],
      source = "Parametric_ALS", plot_id = pid)
  }
  do.call(rbind, Filter(Negate(is.null), rows))
}))

p_pheno <- ggplot(df_pheno_long, aes(x = doy, y = lai, colour = source)) +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~plot_id, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = pal_pheno, na.value = "grey60") +
  labs(title    = "Chapter 3 — Phenology curves per LAI correction",
       subtitle = sprintf("Sample of %d HOBO plots", n_show),
       x = "Day of year", y = expression(LAI ~ (m^2 ~ m^{-2})), colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# 6d. Spatial correction maps: delta LAI rasters (S2_DOPT − S2_ATBD, ALS − S2_ATBD)
r_als  <- terra::aggregate(terra::rast(file.path(CFG_C3$in_dir, "lai_z1_res_10_m.tif")),
                            fact = CFG_C3$agg_factor, fun = "mean", na.rm = TRUE)
r_s2   <- terra::aggregate(terra::rast(file.path(CFG_C3$in_dir,
                                                   "s2lai_summer_atbd_res_10_m.tif")),
                            fact = CFG_C3$agg_factor, fun = "mean", na.rm = TRUE)
r_s2   <- terra::resample(r_s2, r_als, method = "bilinear")
r_dopt <- terra::aggregate(
  terra::rast(file.path(CFG_C3$nc_s2_dopt,
                         sprintf("s2lai_summer_opt_%s_res_10_m.tif", CFG_C3$dopt_variant))),
  fact = CFG_C3$agg_factor, fun = "mean", na.rm = TRUE)
r_dopt <- terra::resample(r_dopt, r_als, method = "bilinear")

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_delta_map <- function(r_a, r_b, title, lim = NULL) {
  df <- as.data.frame(r_a - r_b, xy = TRUE); names(df)[3] <- "delta"
  lim <- lim %||% max(abs(df$delta), na.rm = TRUE)
  ggplot(df, aes(x, y, fill = delta)) + geom_raster() +
    scale_fill_gradient2(low="#31688e", mid="white", high="#d8576b", midpoint=0,
                          limits=c(-lim,lim), na.value="grey90", name="ΔLAI\n(m²/m²)") +
    coord_equal() + labs(title=title, x="Easting (m)", y="Northing (m)") +
    theme_bw(base_size=9) + theme(axis.text=element_text(size=6))
}

lim_shared <- max(
  max(abs(terra::values(r_als - r_s2)), na.rm=TRUE),
  max(abs(terra::values(r_dopt - r_s2)), na.rm=TRUE)
)
p_map_als  <- make_delta_map(r_als,  r_s2, "LAI_ALS − LAI_S2_ATBD",  lim_shared)
p_map_dopt <- make_delta_map(r_dopt, r_s2, "LAI_S2_DOPT − LAI_S2_ATBD", lim_shared)

p_corr_maps <- (p_map_als | p_map_dopt) +
  patchwork::plot_annotation(
    title    = "Chapter 3 — LAI correction spatial patterns",
    subtitle = "Positive = correction gives more LAI than raw S2 | Blois, summer 2021"
  )

# 6e. Paired scatter: all correction vs ALS reference
df_scatter <- tidyr::pivot_longer(
  df_plots[, c("LAI_ALS", intersect(lai_cols[-1], names(df_plots))), drop=FALSE],
  -LAI_ALS, names_to = "method", values_to = "LAI_pred"
)
p_scatter_lai <- ggplot(df_scatter, aes(x = LAI_ALS, y = LAI_pred, colour = method)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(alpha=0.6, size=2) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.7) +
  facet_wrap(~method, ncol=3) +
  scale_colour_viridis_d(option="D") +
  labs(title    = "Chapter 3 — All LAI corrections vs LiDAR reference",
       subtitle = sprintf("n = %d HOBO plots | dashed = 1:1 line", nrow(df_plots)),
       x = expression(LAI[ALS] ~ (m^2 ~ m^{-2})),
       y = expression(LAI[corrected] ~ (m^2 ~ m^{-2})),
       colour = NULL) +
  theme_bw(base_size=10) + theme(legend.position="none")

# Save all figures
if (FLAGS_PREP$SAVE_FIGURES) {
  save_plot(p_dist,        file.path(CFG_C3$figs_dir, "c3_lai_distribution.png"),      9, 5)
  save_plot(p_loo,         file.path(CFG_C3$figs_dir, "c3_rf_loo_cv.png"),             7, 5)
  save_plot(p_pheno,       file.path(CFG_C3$figs_dir, "c3_phenology_sample.png"),     10, 7)
  save_plot(p_corr_maps,   file.path(CFG_C3$figs_dir, "c3_lai_correction_map.png"),   14, 6, dpi=300)
  save_plot(p_scatter_lai, file.path(CFG_C3$figs_dir, "c3_lai_scatter_all.png"),      12, 8)
  cat("   Figures saved.\n")
}

cat("\n=================================================================\n")
cat("  Step 1 complete. Outputs in:", CFG_C3$lai_prep_dir, "\n")
cat("  Next: run Chapter3_main.R\n")
cat("=================================================================\n")
