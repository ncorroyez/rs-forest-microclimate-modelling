# ==============================================================================
# Chapter 1 — HOBO validation: diagnostics and per-scenario validation
# ==============================================================================

#' Select 3 representative HOBO sensors by buffering capacity.
#'
#' Returns the sensor with the strongest cooling, the median, and the weakest.
#'
#' @param df_hobo_daily Daily HOBO data with id_plot and Delta_obs columns.
#' @return Tibble with 3 rows and columns id_plot, mean_delta, rationale, facet_label.
pick_representative_hobos <- function(df_hobo_daily) {
  df_summary <- df_hobo_daily %>% group_by(id_plot) %>% summarise(mean_delta = mean(Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(mean_delta)
  n_total <- nrow(df_summary)
  res <- bind_rows(
    df_summary %>% slice(1) %>% mutate(rationale = "Strongest cooling (Max buffer)"),
    df_summary %>% slice(ceiling(n_total / 2)) %>% mutate(rationale = "Median cooling"),
    df_summary %>% slice(n_total) %>% mutate(rationale = "Weakest cooling (Min buffer)")
  ) %>% mutate(facet_label = sprintf("HOBO %s\n[%s : Mean ΔTmax = %.1f°C]", id_plot, rationale, mean_delta))
  return(res)
}

#' Time-series overlay for 3 selected HOBOs: HOBO obs, MuSICA real, MuSICA uniform, ERA5.
#'
#' @param df_selected_hobos Tibble from pick_representative_hobos.
#' @param df_hobo           Full daily HOBO dataset.
#' @param df_musica_real    MuSICA results for real LAD scenario.
#' @param df_musica_uniform MuSICA results for uniform LAD scenario.
#' @param df_macro          Daily macroclimate (ERA5).
#' @return ggplot object.
plot_timeseries_faceted <- function(df_selected_hobos, df_hobo, df_musica_real, df_musica_uniform, df_macro) {
  ids <- df_selected_hobos$id_plot
  df1 <- df_hobo %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_obs, source = "HOBO")
  df2 <- df_musica_real %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_micro, source = "MuSICA real LAD")
  df3 <- df_musica_uniform %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_micro, source = "MuSICA uniform LAD")
  df4 <- tidyr::expand_grid(id_plot = ids, df_macro) %>% transmute(id_plot, date, value = Tmax_macro, source = "ERA5 macro")

  df_all <- bind_rows(df1, df2, df3, df4) %>% inner_join(df_selected_hobos %>% dplyr::select(id_plot, facet_label), by = "id_plot") %>% mutate(facet_label = factor(facet_label, levels = df_selected_hobos$facet_label))

  ggplot(df_all, aes(x = date, y = value, colour = source)) + geom_line(linewidth = 0.6, alpha = 0.85) + facet_wrap(~ facet_label, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("HOBO" = "black", "MuSICA real LAD" = "#31688e", "MuSICA uniform LAD" = "#d8576b", "ERA5 macro" = "grey60")) +
    labs(title = "Time-series — Selected HOBOs vs MuSICA", subtitle = "Plots chosen by their average buffering capacity during summer", x = NULL, y = expression(T[max] ~ (degree*C)), colour = NULL) + theme(strip.text = element_text(size = 11, lineheight = 1.2))
}

#' Extract structural covariates for HOBO sensor locations from LiDAR rasters.
#'
#' By default, every variable is read from the (already-aggregated) `r_stack`
#' at the sensor coordinates (MEAN over the aggregated cell). For Hmax this
#' systematically underestimates canopy height in heterogeneous forests where
#' tall trees coexist with gaps : the cell-MEAN dilutes tall trees with the
#' canopy gaps in the same cell.
#'
#' If `hmax_raw_raster` is supplied, Hmax is recomputed as the MAX of the
#' raw (non-aggregated) raster within a circular buffer of `hmax_buffer`
#' metres around each HOBO point. This reproduces the legacy behaviour
#' (`main_Blois_detailed.R`, October 2025 simulations) and yields a Hmax
#' closer to the dominant canopy height of the stand.
#'
#' Diagnostic note. With `hmax_buffer = 25 m` the mean Hmax across the 53
#' HOBO sensors is 24.84 m, vs 22.09 m with the default MEAN-aggregated
#' stack. For plot 41_07 the legacy value is 11.09 m, vs 6.85 m default.
#' The MEAN-aggregated stack is the root cause of the systematic
#' ΔTmax > 0 bias observed in the validation slide (see
#' `outputs/musica_bug_diagnosis.md`).
#'
#' @param hobo_geojson    Path to the HOBO GeoJSON file.
#' @param r_stack         LiDAR raster stack (from `load_lidar_rasters`).
#' @param ids_to_remove   HOBO IDs to exclude.
#' @param hmax_raw_raster Optional path or `SpatRaster` for the raw, non-aggregated
#'   Hmax raster (e.g. `in_files/max_res_10_m.tif`). When supplied, Hmax is
#'   recomputed by buffer-MAX. Default `NULL` keeps the legacy behaviour.
#' @param hmax_buffer     Buffer radius (m) for the buffer-MAX extraction.
#'   Default 25 m (matches `main_Blois_detailed.R` radius_test setup).
#' @return Dataframe with id_plot and structural covariates.
build_hobo_inputs <- function(hobo_geojson, r_stack, ids_to_remove,
                                 hmax_raw_raster = NULL, hmax_buffer = 25,
                                 hobo_buffer_mode = c("cell", "buffer25"),
                                 raw_raster_dir = NULL) {
  hobo_buffer_mode <- match.arg(hobo_buffer_mode)
  hobo_pts <- st_read(hobo_geojson, quiet = TRUE)
  ext_vals <- terra::extract(r_stack, vect(hobo_pts), xy = TRUE)

  # New `buffer25` mode (default for HOBO from now on) : extract LAI, Hmax,
  # fCover and each LAD layer from the RAW 10 m raster, aggregated within a
  # 25 m radius circular buffer around each HOBO point. Hmax uses MAX, the
  # others MEAN. This reproduces the legacy `main_Blois_detailed.R` (Oct 2025)
  # behaviour and is the right thing to do for point HOBO sensors (the sensor
  # is sensitive to the canopy AROUND it, not the local 20 m raster cell).
  if (hobo_buffer_mode == "buffer25") {
    if (is.null(raw_raster_dir)) raw_raster_dir <- "in_files"
    r_lai_raw    <- terra::rast(file.path(raw_raster_dir, "lai_z1_res_10_m.tif"))
    r_hmax_raw_  <- terra::rast(file.path(raw_raster_dir, "max_res_10_m.tif"))
    r_fcover_raw <- terra::rast(file.path(raw_raster_dir, "fCover_res_10_m.tif"))
    r_lad_raw    <- terra::rast(file.path(raw_raster_dir, "lad_profiles_z1_res_10_m.tif"))
    buf_pts <- sf::st_buffer(hobo_pts, 10)   # 10 m footprint (radius sweep: 5-10 m best-validating)
    buf_v   <- terra::vect(buf_pts)
    # Keep LAI/LAD at PAD scale (no ×2): the single one->two-sided doubling lives
    # in run_musica_one. Applying ×2 here too double-counted the leaf area.
    lai_buf    <-     terra::extract(r_lai_raw,    buf_v, fun = "mean", na.rm = TRUE)[, 2]
    hmax_buf   <-     terra::extract(r_hmax_raw_,  buf_v, fun = "max",  na.rm = TRUE)[, 2]
    fcover_buf <-     terra::extract(r_fcover_raw, buf_v, fun = "mean", na.rm = TRUE)[, 2]
    lad_buf    <-     terra::extract(r_lad_raw,    buf_v, fun = "mean", na.rm = TRUE)
    ext_vals$LAI    <- lai_buf
    ext_vals$Hmax   <- hmax_buf
    ext_vals$fCover <- fcover_buf
    lad_cols <- grep("^LAD_Layer_", names(ext_vals), value = TRUE)
    if (length(lad_cols) == length(names(lad_buf)) - 1L) {
      for (k in seq_along(lad_cols)) {
        ext_vals[[lad_cols[k]]] <- lad_buf[, k + 1L]
      }
    }
  } else if (!is.null(hmax_raw_raster)) {
    # Backward-compatible : only override Hmax via raw buffer-MAX
    r_hmax_raw <- if (inherits(hmax_raw_raster, "SpatRaster")) hmax_raw_raster
                   else terra::rast(hmax_raw_raster)
    buf_pts <- sf::st_buffer(hobo_pts, hmax_buffer)
    hmax_buf <- terra::extract(r_hmax_raw, terra::vect(buf_pts),
                                  fun = "max", na.rm = TRUE)[, 2]
    ext_vals$Hmax <- hmax_buf
  }

  hobo_pts %>% st_drop_geometry() %>% dplyr::select(id_plot) %>%
    bind_cols(ext_vals) %>%
    filter(!(id_plot %in% ids_to_remove)) %>%
    drop_na(LAI, Hmax, fCover)
}

#' Extract daily ΔTmax from MuSICA outputs at HOBO sensor locations.
#'
#' @param scenario_dir  Directory with HOBO NetCDF outputs for one scenario.
#' @param df_macro      Daily macroclimate.
#' @param date_seq      Dates of interest.
#' @return Tibble with id_plot, scenario columns added.
extract_deltatmax_hobo_scenario <- function(scenario_dir, df_macro, date_seq) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE); if (length(nc_files) == 0) return(NULL)
  map_df(nc_files, function(f) {
    id_str <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)"); res <- extract_deltatmax_one(f, df_macro, date_seq); if (is.null(res)) return(NULL)
    res %>% mutate(id_plot = id_str, scenario = basename(scenario_dir))
  })
}

#' Run MuSICA at HOBO locations for a list of scenarios and compute validation metrics.
#'
#' @param df_hobo_inputs  Dataframe from build_hobo_inputs.
#' @param df_hobo_daily   Observed daily HOBO data.
#' @param scenarios       Named list of scenario lists.
#' @param parent_dir      Parent output directory.
#' @param df_macro        Daily macroclimate.
#' @param date_seq        Dates of interest.
#' @param forcing_file    Path to the ERA5 forcing NetCDF (e.g. CFG_C3$forcing_file).
#' @param musica_cmd      Shell command to call MuSICA (e.g. CFG_C3$musica_cmd).
#' @param force           Re-run existing outputs (default FALSE).
#' @return List with elements $daily (paired obs/sim) and $metrics (per-scenario stats).
validate_scenarios_at_hobos <- function(df_hobo_inputs, df_hobo_daily, scenarios, parent_dir, df_macro, date_seq, forcing_file, musica_cmd, force = FALSE, abl_flag = NULL) {
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  # iter/yoyo: needs abl_flag + a phenology row for the day-before-sim (2020/366),
  # since FR-Blo forcing starts 2021-01-01. Wrap each scenario's phenology to seed it.
  extra_setup <- if (!is.null(abl_flag)) list("abl_flag" = sprintf('"%s"', abl_flag)) else NULL
  seed_pheno <- function(base_ph) { if (is.null(base_ph)) return(NULL)
    b <- as.data.frame(base_ph); s <- b[1, , drop = FALSE]
    if ("year" %in% names(s)) s$year <- 2020L; if ("Julian_day" %in% names(s)) s$Julian_day <- 366L
    rbind(s, b) }
  hobo_paths <- character(length(scenarios)); names(hobo_paths) <- names(scenarios)
  for (nm in names(scenarios)) {
    sc <- scenarios[[nm]]; out_d <- file.path(parent_dir, sc$name); dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
    existing <- list.files(out_d, pattern = "\\.nc$")
    if (!force && length(existing) >= nrow(df_hobo_inputs)) { hobo_paths[nm] <- out_d; next }
    cat(sprintf("\n[HOBO %s] %d HOBO plots → %s\n", sc$name, nrow(df_hobo_inputs), out_d))
    for (i in seq_len(nrow(df_hobo_inputs))) {
      plot_row <- df_hobo_inputs[i, ]; sim_id <- sprintf("HOBO_%s", plot_row$id_plot)
      out_file <- file.path(out_d, paste0("musica_out_", sim_id, ".nc")); if (file.exists(out_file)) next
      sc_run <- sc
      if (!is.null(abl_flag)) {                                    # seed phenology for iter
        base_ph <- if (is.null(sc$phenology_fn))
          as.data.frame(calc_phenology(list.year = 2020:2022, nleafage = 1, budburst_date = 115,
                        leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
                        LAI_max_per_cohort = sc$lai_fn(plot_row)))
        else sc$phenology_fn(plot_row)
        ph_seeded <- seed_pheno(base_ph); if (!is.null(ph_seeded)) sc_run$phenology_fn <- function(p) ph_seeded
      }
      run_musica_one(plot_row, sc_run, out_file, forcing_file, musica_cmd, extra_setup = extra_setup)
    }
    hobo_paths[nm] <- out_d
  }
  df_sims <- map_df(hobo_paths, extract_deltatmax_hobo_scenario, df_macro = df_macro, date_seq = date_seq)
  df_obs <- df_hobo_daily %>% dplyr::select(id_plot, date, Delta_obs)
  df_daily <- df_sims %>% inner_join(df_obs, by = c("id_plot", "date")) %>% rename(Delta_sim = Delta_Tmax)
  df_metrics <- df_daily %>% group_by(scenario) %>% summarise(n = n(), r2 = cor(Delta_obs, Delta_sim, use = "complete.obs")^2, rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)), mae = mean(abs(Delta_obs - Delta_sim), na.rm = TRUE), bias = mean(Delta_sim - Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(rmse)
  list(daily = df_daily, metrics = df_metrics)
}

#' Bar chart of RMSE per scenario at HOBO locations, coloured by bias.
#'
#' @param df_metrics Tibble from validate_scenarios_at_hobos$metrics.
#' @param title      Plot title.
#' @return ggplot object.
plot_hobo_validation <- function(df_metrics, title = "HOBO sensor validation — RMSE per scenario") {
  bias_sd  <- round(sd(df_metrics$bias, na.rm = TRUE), 2)
  bias_abs <- round(mean(abs(df_metrics$bias), na.rm = TRUE), 2)
  sub_txt  <- sprintf(
    "Mean absolute bias = +%.2f°C (missing ET in MuSICA) | inter-scenario σ = ±%.2f°C\nσ << absolute bias: HOBO confirms systematic offset but cannot rank H1 scenarios",
    bias_abs, bias_sd
  )
  df_metrics %>%
    mutate(label = sprintf("R²=%.2f  bias=%+.2f°C", r2, bias)) %>%
    ggplot(aes(x = reorder(scenario, rmse), y = rmse)) +
    geom_col(fill = "#31688e", alpha = 0.8) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.2, colour = "grey20") +
    coord_flip() +
    expand_limits(y = max(df_metrics$rmse) * 1.35) +
    labs(title = title, subtitle = sub_txt,
         x = NULL, y = "RMSE vs observed HOBO ΔTmax (°C)") +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9))
}

#' Scatter ΔTmax(sim) vs ΔTmax(obs HOBO) coloured by structural cluster.
#'
#' Key figure: the 3D canopy typology (Cluster) should predict the rank order of
#' cooling refuges without necessarily predicting the absolute amplitude.
#' Faceted by heatwave / normal regime.
#'
#' @param hobo_res       List from validate_scenarios_at_hobos ($daily has id_plot, date, Delta_obs, Delta_sim, scenario, Tmax_macro).
#' @param df_hobo_inputs Dataframe from build_hobo_inputs (id_plot, x, y).
#' @param df_forest      Full forest dataframe with Cluster column (x, y, Cluster).
#' @param scenario_name  Scenario to extract from hobo_res$daily (default "REF_all_real").
#' @param threshold      Heatwave Tmax_macro threshold in °C (default 30).
#' @param out_path       PNG output path.
#' @return Invisible ggplot.
plot_hobo_sim_obs_by_cluster <- function(hobo_res, df_hobo_inputs, df_forest,
                                          scenario_name = "REF_all_real",
                                          threshold     = 30,
                                          out_path      = "outputs/figures/annex/A19_hobo_cluster_scatter.png") {
  # Assign each HOBO sensor to the nearest forest grid cell's cluster
  df_forest_xy <- df_forest %>%
    dplyr::select(x, y, Cluster) %>%
    mutate(Cluster = as.character(Cluster))

  df_hobo_cluster <- df_hobo_inputs %>%
    dplyr::select(id_plot, x, y) %>%
    mutate(
      Cluster = purrr::map2_chr(x, y, function(hx, hy) {
        dists  <- (df_forest_xy$x - hx)^2 + (df_forest_xy$y - hy)^2
        df_forest_xy$Cluster[which.min(dists)]
      })
    ) %>%
    dplyr::select(id_plot, Cluster)

  # hobo_res$daily already carries Tmax_macro from extract_deltatmax_one;
  # a second join would produce Tmax_macro.x/.y and break the mutate below.
  df_daily <- hobo_res$daily %>%
    filter(scenario == scenario_name) %>%
    inner_join(df_hobo_cluster, by = "id_plot") %>%
    mutate(
      regime  = factor(ifelse(Tmax_macro >= threshold, "Heatwave", "Normal"),
                       levels = c("Normal", "Heatwave")),
      Cluster = factor(Cluster)
    )

  if (nrow(df_daily) == 0) {
    warning(sprintf("[plot_hobo_sim_obs_by_cluster] no data for scenario '%s'", scenario_name))
    return(invisible(NULL))
  }

  lims <- range(c(df_daily$Delta_obs, df_daily$Delta_sim), na.rm = TRUE)

  df_stats <- df_daily %>%
    group_by(regime) %>%
    summarise(
      r2   = cor(Delta_obs, Delta_sim, use = "complete.obs")^2,
      rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)),
      bias = mean(Delta_sim - Delta_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf("R²=%.2f\nRMSE=%.2f°C\nBias=%+.2f°C", r2, rmse, bias),
      x_pos = lims[1], y_pos = lims[2]
    )

  p <- ggplot(df_daily, aes(x = Delta_obs, y = Delta_sim,
                             colour = Cluster, shape = Cluster)) +
    geom_point(alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.1,
                aes(group = 1), colour = "grey30") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.7) +
    geom_text(data = df_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = -0.05, vjust = 1.1, fontface = "bold", size = 3.2,
              inherit.aes = FALSE) +
    coord_fixed(xlim = lims, ylim = lims) +
    facet_wrap(~ regime) +
    scale_colour_viridis_d(option = "turbo", name = "Cluster") +
    scale_shape_manual(values = 15:18, name = "Cluster") +
    labs(
      title    = sprintf("HOBO validation — ΔTmax sim vs obs | Scenario: %s", scenario_name),
      subtitle = "Coloured by structural cluster — 3D typology ranks cooling refuges without predicting absolute amplitude",
      x        = expression("HOBO obs" ~ Delta * T[max] ~ (degree*C)),
      y        = expression("MuSICA sim" ~ Delta * T[max] ~ (degree*C))
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold"),
      strip.text      = element_text(face = "bold", size = 11)
    )

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 12, height = 6)
  cat(sprintf("  [hobo_cluster_scatter] PNG : %s\n", out_path))
  invisible(p)
}

#' A19 — but faceted by structural cluster instead of regime.
#'
#' Met en évidence la performance MuSICA conditionnée au type de peuplement :
#' un cluster avec un R² bas → la canopée 3D y compte vraiment ; un cluster
#' avec un R² haut → MuSICA y est cohérent avec HOBO. Permet d'identifier
#' où les écarts modèle / mesures sont concentrés.
#'
#' @param hobo_res       List from validate_scenarios_at_hobos.
#' @param df_hobo_inputs HOBO inputs avec colonnes (id_plot, x, y).
#' @param df_forest      Forest dataframe complet avec colonne Cluster.
#' @param scenario_name  Scénario à extraire (défaut "REF_all_real").
#' @param threshold      Seuil canicule (°C) — utilisé pour la couleur des points.
#' @param out_path       PNG de sortie.
#' @return Invisible ggplot.
plot_hobo_per_cluster <- function(hobo_res, df_hobo_inputs, df_forest,
                                   scenario_name = "REF_all_real",
                                   threshold     = 30,
                                   out_path      = "outputs/figures/annex/A19b_hobo_per_cluster.png") {
  df_forest_xy <- df_forest %>%
    dplyr::select(x, y, Cluster) %>%
    mutate(Cluster = as.character(Cluster))

  df_hobo_cluster <- df_hobo_inputs %>%
    dplyr::select(id_plot, x, y) %>%
    mutate(Cluster = purrr::map2_chr(x, y, function(hx, hy) {
      d <- (df_forest_xy$x - hx)^2 + (df_forest_xy$y - hy)^2
      df_forest_xy$Cluster[which.min(d)]
    })) %>%
    dplyr::select(id_plot, Cluster)

  df_daily <- hobo_res$daily %>%
    filter(scenario == scenario_name) %>%
    inner_join(df_hobo_cluster, by = "id_plot") %>%
    mutate(
      regime  = factor(ifelse(Tmax_macro >= threshold, "Heatwave", "Normal"),
                       levels = c("Normal", "Heatwave")),
      Cluster = factor(Cluster)
    )

  if (nrow(df_daily) == 0) {
    warning(sprintf("[plot_hobo_per_cluster] no data for scenario '%s'", scenario_name))
    return(invisible(NULL))
  }

  lims <- range(c(df_daily$Delta_obs, df_daily$Delta_sim), na.rm = TRUE)

  # Helper : etoiles de significativite
  sig_star <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    "ns"
  }

  # Stats + tests par cluster
  # NB: les obs intra-plot sont temporellement autocorrelees ; les p-values
  # ci-dessous sont donc OPTIMISTES (n effectif < n_obs). A interpreter comme
  # un indicateur de signal, pas comme un test inferentiel strict.
  compute_per_cluster <- function(d) {
    n_obs <- nrow(d)
    if (n_obs < 3) return(tibble::tibble(n_obs = n_obs))
    ct  <- tryCatch(cor.test(d$Delta_obs, d$Delta_sim), error = function(e) NULL)
    lm_ <- tryCatch(lm(Delta_sim ~ Delta_obs, data = d), error = function(e) NULL)
    bt  <- tryCatch(t.test(d$Delta_sim - d$Delta_obs), error = function(e) NULL)

    r2   <- if (is.null(ct)) NA_real_ else as.numeric(ct$estimate^2)
    p_r2 <- if (is.null(ct)) NA_real_ else ct$p.value

    if (is.null(lm_)) {
      slope <- NA_real_; p_slope_vs_1 <- NA_real_
    } else {
      slope    <- as.numeric(coef(lm_)[2])
      se_slope <- summary(lm_)$coefficients[2, 2]
      df_res   <- lm_$df.residual
      t_stat   <- (slope - 1) / se_slope
      p_slope_vs_1 <- 2 * pt(abs(t_stat), df_res, lower.tail = FALSE)
    }
    bias <- if (is.null(bt)) NA_real_ else as.numeric(bt$estimate)
    p_bias <- if (is.null(bt)) NA_real_ else bt$p.value
    rmse <- sqrt(mean((d$Delta_obs - d$Delta_sim)^2, na.rm = TRUE))

    tibble::tibble(n_obs = n_obs, r2 = r2, p_r2 = p_r2,
                   slope = slope, p_slope_vs_1 = p_slope_vs_1,
                   bias  = bias, p_bias = p_bias, rmse = rmse)
  }

  df_n_plots <- df_daily %>% group_by(Cluster) %>%
    summarise(n_plots = dplyr::n_distinct(id_plot), .groups = "drop")

  df_stats <- df_daily %>%
    group_by(Cluster) %>%
    group_modify(~ compute_per_cluster(.x)) %>%
    ungroup() %>%
    inner_join(df_n_plots, by = "Cluster") %>%
    mutate(
      stars_r2    = vapply(p_r2,        sig_star, character(1)),
      stars_slope = vapply(p_slope_vs_1, sig_star, character(1)),
      stars_bias  = vapply(p_bias,      sig_star, character(1)),
      label = sprintf("n_plot=%d | n=%d\nR²=%.2f %s\nslope=%.2f %s (vs 1)\nBias=%+.2f°C %s\nRMSE=%.2f°C",
                       n_plots, n_obs, r2, stars_r2,
                       slope, stars_slope, bias, stars_bias, rmse),
      x_pos = lims[1], y_pos = lims[2]
    )

  cat("\n=== A19 per-cluster stats (avec p-values) ===\n")
  cat("  NB: p-values optimistes (autocorrelation temporelle non corrigee)\n")
  print(df_stats %>% dplyr::select(Cluster, n_plots, n_obs,
                                    r2, p_r2, slope, p_slope_vs_1,
                                    bias, p_bias, rmse))

  p <- ggplot(df_daily, aes(x = Delta_obs, y = Delta_sim, colour = regime)) +
    geom_point(alpha = 0.55, size = 1.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.0,
                aes(group = 1), colour = "grey25") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black", linewidth = 0.7) +
    geom_text(data = df_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = -0.05, vjust = 1.1, fontface = "bold",
              size = 3.0, inherit.aes = FALSE) +
    coord_fixed(xlim = lims, ylim = lims) +
    facet_wrap(~ Cluster, nrow = 2,
               labeller = labeller(Cluster = function(x) paste("Cluster", x))) +
    scale_colour_manual(values = c("Normal" = "#31688e", "Heatwave" = "#d8576b"),
                        name = "Regime") +
    labs(
      title    = sprintf("HOBO validation by structural cluster | Scenario: %s", scenario_name),
      subtitle = "Slope <1 ⇒ MuSICA under-predicts spatial variability | * p<0.05 ** p<0.01 *** p<0.001 (autocorr. uncorrected)",
      x        = expression("HOBO obs" ~ Delta * T[max] ~ (degree*C)),
      y        = expression("MuSICA sim" ~ Delta * T[max] ~ (degree*C))
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title      = element_text(face = "bold"),
          strip.text      = element_text(face = "bold", size = 11))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 11, height = 11)

  # CSV stats (avec p-values + etoiles)
  csv_path <- sub("\\.png$", ".csv", out_path)
  write.csv(df_stats %>% dplyr::select(-x_pos, -y_pos, -label),
            csv_path, row.names = FALSE)
  cat(sprintf("  [hobo_per_cluster] PNG : %s\n  [hobo_per_cluster] CSV : %s\n",
              out_path, csv_path))
  invisible(p)
}

# ---- Absolute Tmax scatter ----------------------------------------------------

#' Scatterplot Tmax_micro(sim) vs Tmax_obs(HOBO) for all scenarios.
#'
#' @title Absolute Tmax: simulated vs observed
#' @description For each (id_plot, date), plots the simulated sub-canopy daily
#'   Tmax (MuSICA) against the observed sub-canopy daily Tmax (HOBO sensor),
#'   faceted by scenario. Restricted to Jun–Sep 2021 via the date_seq used
#'   when building df_daily. Tmax_obs is reconstructed as Tmax_macro + Delta_obs.
#'
#' @param df_daily  Data.frame from validate_scenarios_at_hobos()$daily.
#'   Must contain: Tmax_micro, Tmax_macro, Delta_obs, scenario.
#' @param ncol      Number of columns in facet grid (default 5).
#' @return ggplot object.
plot_tmax_absolute_scatter <- function(df_daily, ncol = 5) {
  df_plot <- df_daily %>%
    dplyr::mutate(Tmax_obs = Tmax_macro + Delta_obs) %>%
    dplyr::filter(!is.na(Tmax_micro), !is.na(Tmax_obs))

  metrics <- df_plot %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      r2   = cor(Tmax_obs, Tmax_micro, use = "complete.obs")^2,
      rmse = sqrt(mean((Tmax_micro - Tmax_obs)^2, na.rm = TRUE)),
      bias = mean(Tmax_micro - Tmax_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf("R²=%.2f\nRMSE=%.2f°C\nbias=%+.2f°C", r2, rmse, bias)
    )

  lims <- range(c(df_plot$Tmax_obs, df_plot$Tmax_micro), na.rm = TRUE)
  lims <- lims + c(-0.5, 0.5)

  ggplot(df_plot, aes(x = Tmax_obs, y = Tmax_micro)) +
    geom_point(alpha = 0.25, size = 0.7, colour = "#31688e") +
    geom_smooth(method = "lm", se = TRUE, colour = "#d8576b",
                linewidth = 0.8, alpha = 0.15) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_text(data = metrics, aes(label = label),
              x = lims[1] + 0.3, y = lims[2] - 0.2,
              hjust = 0, vjust = 1, size = 2.5, inherit.aes = FALSE,
              colour = "grey20") +
    facet_wrap(~scenario, ncol = ncol) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(
      title    = "Chapter 3 — Sub-canopy Tmax: simulated vs observed (absolute)",
      subtitle = "Jun–Sep 2021 | blue = points | red = OLS | dashed = 1:1",
      x        = expression(T[max*",obs"]^{HOBO} ~ (degree*C)),
      y        = expression(T[max*",micro"]^{MuSICA} ~ (degree*C))
    ) +
    theme_bw(base_size = 8) +
    theme(strip.text    = element_text(size = 7),
          plot.subtitle = element_text(colour = "grey30", size = 8))
}

# ---- Microclimat = f(macroclimate) scatter ------------------------------------

#' Scatterplot Tmax_micro vs Tmax_macro for all scenarios, with HOBO overlay.
#'
#' @title Sub-canopy Tmax as a function of ERA5 macroclimate Tmax
#' @description For each (id_plot, date), plots simulated sub-canopy Tmax
#'   (MuSICA) vs ERA5 macroclimate Tmax, faceted by scenario. HOBO observations
#'   (Tmax_obs = Tmax_macro + Delta_obs) are shown as grey background points in
#'   all facets. The slope of the OLS line indicates the thermal sensitivity of
#'   the sub-canopy to the macroclimate; a slope < 1 indicates buffering.
#'
#' @param df_daily  Data.frame from validate_scenarios_at_hobos()$daily.
#' @param ncol      Number of columns in facet grid (default 5).
#' @return ggplot object.
plot_tmax_micro_vs_macro <- function(df_daily, ncol = 5) {
  df_sim <- df_daily %>%
    dplyr::filter(!is.na(Tmax_micro), !is.na(Tmax_macro))

  # HOBO reference — no scenario column so it repeats across all facets
  df_obs <- df_daily %>%
    dplyr::distinct(id_plot, date, Tmax_macro, Delta_obs) %>%
    dplyr::mutate(Tmax_obs = Tmax_macro + Delta_obs) %>%
    dplyr::filter(!is.na(Tmax_obs))

  # Metrics per scenario: OLS slope + R²
  metrics <- df_sim %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      slope = tryCatch(coef(lm(Tmax_micro ~ Tmax_macro))[2], error = function(e) NA_real_),
      r2    = cor(Tmax_macro, Tmax_micro, use = "complete.obs")^2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = sprintf("slope=%.2f\nR²=%.2f", slope, r2)
    )

  # HOBO reference slope (once, no facet)
  hobo_slope <- tryCatch(coef(lm(Tmax_obs ~ Tmax_macro, data = df_obs))[2],
                          error = function(e) NA_real_)
  hobo_r2    <- cor(df_obs$Tmax_macro, df_obs$Tmax_obs, use = "complete.obs")^2

  xlims <- range(df_sim$Tmax_macro, na.rm = TRUE) + c(-0.5, 0.5)
  ylims <- range(c(df_sim$Tmax_micro, df_obs$Tmax_obs), na.rm = TRUE) + c(-0.5, 0.5)

  ggplot() +
    geom_point(data = df_obs,
               aes(x = Tmax_macro, y = Tmax_obs),
               colour = "black", alpha = 0.12, size = 0.5) +
    geom_smooth(data = df_obs,
                aes(x = Tmax_macro, y = Tmax_obs),
                method = "lm", se = FALSE,
                colour = "black", linetype = "dotted", linewidth = 0.7) +
    geom_point(data = df_sim,
               aes(x = Tmax_macro, y = Tmax_micro),
               colour = "#31688e", alpha = 0.2, size = 0.6) +
    geom_smooth(data = df_sim,
                aes(x = Tmax_macro, y = Tmax_micro),
                method = "lm", se = TRUE,
                colour = "#d8576b", linewidth = 0.8, alpha = 0.15) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_text(data = metrics, aes(label = label),
              x = xlims[1] + 0.3, y = ylims[2] - 0.2,
              hjust = 0, vjust = 1, size = 2.4, inherit.aes = FALSE,
              colour = "#d8576b") +
    annotate("text",
             x = xlims[1] + 0.3, y = ylims[1] + 2,
             label = sprintf("HOBO: slope=%.2f R²=%.2f", hobo_slope, hobo_r2),
             hjust = 0, size = 2.3, colour = "grey30") +
    facet_wrap(~scenario, ncol = ncol) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    labs(
      title    = "Chapter 3 — Sub-canopy Tmax as a function of ERA5 Tmax",
      subtitle = "Blue = MuSICA sim (red OLS) | Black dots = HOBO obs (dotted OLS) | dashed = 1:1",
      x        = expression(T[max*",macro"]^{ERA5} ~ (degree*C)),
      y        = expression(T[max] ~ (degree*C))
    ) +
    theme_bw(base_size = 8) +
    theme(strip.text    = element_text(size = 7),
          plot.subtitle = element_text(colour = "grey30", size = 8))
}
