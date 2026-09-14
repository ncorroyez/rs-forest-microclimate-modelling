# ==============================================================================
# Chapter 3 — Heatwave temporal shift analysis
#
# Tests whether dynamic scenarios (DYN_*) better reproduce the observed lag
# between sub-canopy Tmax and ERA5 Tmax on heatwave days. A positive shift
# (sub-canopy Tmax peaks later than ERA5 Tmax) is physically expected: the
# canopy absorbs and re-emits radiation, delaying the under-canopy thermal
# signal by 1–3 h. Static scenarios (fixed LAI d'été) cannot modulate this
# delay seasonally; dynamic scenarios can.
#
# Requires rmusica / musica.tools to be loaded (for get_variable()).
# ==============================================================================

# ---- 1. Per-file extraction ---------------------------------------------------

#' Extract the hour of daily Tmax from a single MuSICA NetCDF output.
#'
#' @title Hour of sub-canopy Tmax from MuSICA output
#' @description Reads Tair interpolated at a fixed absolute height
#'   (default `CFG$tair_target_height` = 1 m, via `get_tair_at_z()`),
#'   and for each date returns the hour at which Tmax occurs.
#'   Plot-invariant et aligné sur la hauteur des sondes HOBO.
#'   Time is used as-is (Europe/Paris), matching the ERA5 forcing convention.
#'
#' @param nc_file   Path to a MuSICA NetCDF output.
#' @param date_seq  Date vector to filter on.
#' @param z_target  Hauteur absolue (m) ; defaut CFG$tair_target_height (1.0).
#' @return Data.frame with columns `date`, `Tmax_micro` (°C), `hour_of_tmax`.
extract_tmax_timing_one <- function(nc_file, date_seq,
                                     z_target = if (exists("CFG") && !is.null(CFG$tair_target_height)) CFG$tair_target_height else 1.0) {
  nc  <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  res <- get_tair_at_z(nc, z_target = z_target)
  nc_close(nc)
  if (is.null(res) || nrow(res) == 0) return(NULL)

  # NB: NO timezone correction here.
  # The ERA5 forcing file is timestamped in Europe/Paris local time (despite
  # force_utc_nc labelling it UTC), and MuSICA NC inherits the same convention
  # (get_variable() returns Europe/Paris). Subtracting 2h would systematically
  # shift MuSICA -2h relative to ERA5, reversing the sign of the physical lag.
  res %>%
    mutate(
      date = as.Date(time),
      hour = as.integer(format(time, "%H"))
    ) %>%
    filter(date %in% date_seq) %>%
    group_by(date) %>%
    filter(Tair_sim == max(Tair_sim, na.rm = TRUE)) %>%
    slice(1) %>%
    summarise(Tmax_micro = Tair_sim, hour_of_tmax = hour, .groups = "drop")
}

#' Extract the hour of daily Tmax from HOBO sub-canopy temperature loggers.
#'
#' @title Hour of HOBO Tmax
#' @description Reads the raw hourly HOBO CSV without pre-aggregating to daily,
#'   and for each (id_plot, date) returns the hour at which maximum temperature
#'   was observed.
#'
#' @param csv_file       Path to the HOBO temperature CSV.
#' @param date_seq       Date vector to filter on.
#' @param ids_to_remove  HOBO plot IDs to exclude.
#' @return Data.frame with columns `id_plot`, `date`, `Tmax_obs` (°C),
#'   `hour_of_tmax` (integer, UTC).
read_hobo_tmax_timing <- function(csv_file, date_seq, ids_to_remove) {
  read.csv(csv_file) %>%
    mutate(
      datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      date     = as.Date(datetime),
      hour     = as.integer(format(datetime, "%H"))
    ) %>%
    filter(position_sensor == "a",
           date %in% date_seq,
           !(id_plot %in% ids_to_remove)) %>%
    group_by(id_plot, date) %>%
    filter(t_hobo == max(t_hobo, na.rm = TRUE)) %>%
    slice(1) %>%
    summarise(Tmax_obs = t_hobo, hour_of_tmax = hour, .groups = "drop")
}

#' Extract the hour of daily Tmax from ERA5 macroclimate forcing.
#'
#' @title Hour of ERA5 Tmax
#' @description Reads the `Tair` variable from the ERA5 NetCDF and returns,
#'   for each date, the hour (UTC) at which the daily maximum occurs.
#'
#' @param forcing_file  Path to the ERA5 NetCDF forcing file.
#' @param date_seq      Date vector to filter on.
#' @return Data.frame with columns `date`, `Tmax_macro` (°C),
#'   `hour_of_tmax_macro` (integer, UTC).
extract_macro_tmax_timing <- function(forcing_file, date_seq) {
  nc         <- nc_open(forcing_file)
  macro_tair <- ncvar_get(nc, "Tair")
  nc_close(nc)
  macro_time <- force_utc_nc(forcing_file, "time")

  data.frame(time = macro_time, Tair = macro_tair - 273.15) %>%
    mutate(date = as.Date(time), hour = as.integer(format(time, "%H"))) %>%
    filter(date %in% date_seq) %>%
    group_by(date) %>%
    filter(Tair == max(Tair, na.rm = TRUE)) %>%
    slice(1) %>%
    summarise(Tmax_macro = Tair, hour_of_tmax_macro = hour, .groups = "drop")
}

# ---- 2. Orchestration --------------------------------------------------------

#' Compute heatwave temporal shifts for all scenarios and HOBO observations.
#'
#' @title Compute Tmax temporal shifts on heatwave days
#' @description For all heatwave days (ERA5 Tmax ≥ `heatwave_threshold`),
#'   computes the lag (in hours) between the hour of sub-canopy Tmax and the
#'   hour of ERA5 Tmax, for both observed HOBO sensors and each MuSICA scenario.
#'
#'   A positive shift indicates that sub-canopy temperatures peak later than
#'   the macroclimate — the expected signature of canopy thermal buffering.
#'   Comparing this lag across HOBO observations and scenarios reveals whether
#'   dynamic parameterisation better reproduces the observed temporal structure.
#'
#' @param scenario_names     Character vector of scenario names (e.g. from
#'   `names(scenarios_c3)`).
#' @param scenario_nc_root   Parent directory of scenario NC outputs
#'   (sub-directories are named by scenario: `<root>/<scenario_name>/`).
#' @param csv_file           Path to the HOBO temperature CSV.
#' @param forcing_file       Path to the ERA5 NetCDF forcing file.
#' @param date_seq           Date vector (full season, used to filter ERA5 timing).
#' @param ids_to_remove      HOBO plot IDs to exclude.
#' @param heatwave_threshold ERA5 Tmax threshold (°C) for heatwave days (default 30).
#' @return Named list:
#'   \describe{
#'     \item{`hobo_shift`}{Data.frame: id_plot, date, Tmax_obs, hour_of_tmax,
#'       hour_of_tmax_macro, shift (hours).}
#'     \item{`sim_shift`}{Data.frame: id_plot, date, scenario, Tmax_micro,
#'       hour_of_tmax, hour_of_tmax_macro, Tmax_macro, shift (hours).}
#'     \item{`hw_dates`}{Date vector of heatwave days.}
#'     \item{`df_macro_timing`}{ERA5 timing data.frame for all heatwave days.}
#'   }
#'   Returns NULL if no heatwave days are found.
compute_tmax_shift_scenarios <- function(scenario_names, scenario_nc_root,
                                          csv_file, forcing_file,
                                          date_seq, ids_to_remove,
                                          heatwave_threshold = 30) {
  # 1. ERA5 timing — full season, then filter to heatwave days
  cat(sprintf("  Computing ERA5 Tmax timing (threshold ≥ %g°C)...\n",
              heatwave_threshold))
  df_macro_timing <- extract_macro_tmax_timing(forcing_file, date_seq)
  hw_dates <- df_macro_timing$date[df_macro_timing$Tmax_macro >= heatwave_threshold]

  if (length(hw_dates) == 0) {
    warning(sprintf(
      "No heatwave days found (ERA5 Tmax ≥ %g°C) in the supplied date_seq.",
      heatwave_threshold
    ))
    return(NULL)
  }
  cat(sprintf("  %d heatwave days identified.\n", length(hw_dates)))
  df_macro_hw <- df_macro_timing[df_macro_timing$date %in% hw_dates, ]

  # 2. HOBO timing on heatwave days
  cat("  Extracting HOBO Tmax timing...\n")
  df_hobo_timing <- read_hobo_tmax_timing(csv_file, hw_dates, ids_to_remove)
  df_hobo_shift  <- df_hobo_timing %>%
    inner_join(df_macro_hw[, c("date", "hour_of_tmax_macro")], by = "date") %>%
    mutate(shift = hour_of_tmax - hour_of_tmax_macro)

  # 3. MuSICA timing per scenario on heatwave days
  cat(sprintf("  Extracting MuSICA Tmax timing (%d scenarios)...\n",
              length(scenario_names)))
  df_sim_shift <- do.call(rbind, lapply(scenario_names, function(sc) {
    nc_dir   <- file.path(scenario_nc_root, sc)
    nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
    if (length(nc_files) == 0) {
      message(sprintf("  [timing] No NC files found for scenario '%s' in %s", sc, nc_dir))
      return(NULL)
    }
    do.call(rbind, lapply(nc_files, function(f) {
      id_str <- sub(".*HOBO_([^.]+)\\.nc$", "\\1", basename(f))
      res    <- extract_tmax_timing_one(f, hw_dates)
      if (is.null(res)) return(NULL)
      res$id_plot  <- id_str
      res$scenario <- sc
      res
    }))
  }))

  if (nrow(df_sim_shift) == 0) {
    warning("No MuSICA timing data extracted — NC files may be missing or unreadable.")
    return(NULL)
  }

  df_sim_shift <- df_sim_shift %>%
    inner_join(df_macro_hw[, c("date", "hour_of_tmax_macro", "Tmax_macro")],
               by = "date") %>%
    mutate(shift = hour_of_tmax - hour_of_tmax_macro)

  list(
    hobo_shift      = df_hobo_shift,
    sim_shift       = df_sim_shift,
    hw_dates        = hw_dates,
    df_macro_timing = df_macro_hw
  )
}

# ---- 3. Visualisation --------------------------------------------------------

#' Boxplot of heatwave Tmax temporal shifts: HOBO vs MuSICA scenarios.
#'
#' @title Plot heatwave Tmax temporal shift (F6)
#' @description For each source (HOBO obs + one box per scenario), shows the
#'   distribution of hourly lags (shift = hour_sub-canopy_Tmax − hour_ERA5_Tmax)
#'   on heatwave days. A positive shift = sub-canopy Tmax peaks after ERA5 Tmax
#'   (thermal buffering lag). HOBO observations set the empirical target; DYN
#'   scenarios are expected to outperform STATIC scenarios in reproducing this lag.
#'
#' @param shift_res          List from `compute_tmax_shift_scenarios()`.
#' @param pal_sc             Named character vector of colours per scenario (same
#'   as `pal_sc` built in `Chapter3_main.R`). NULL uses built-in defaults.
#' @param heatwave_threshold Threshold used (for subtitle label), default 30.
#' @return ggplot object.
plot_tmax_shift <- function(shift_res, pal_sc = NULL, heatwave_threshold = 30) {
  df_obs <- shift_res$hobo_shift %>%
    dplyr::transmute(id_plot, date, shift, source = "HOBO obs")

  df_sim <- shift_res$sim_shift %>%
    dplyr::transmute(id_plot, date, shift, source = scenario)

  sc_order <- c("HOBO obs", sort(unique(df_sim$source)))
  df_all   <- dplyr::bind_rows(df_obs, df_sim) %>%
    dplyr::mutate(source = factor(source, levels = sc_order))

  # Default palette — extend with HOBO black
  default_pal <- c(
    "HOBO obs"       = "black",
    "STATIC_ALS"     = "#31688e",
    "STATIC_S2_ATBD" = "#35b779",
    "STATIC_S2_DOPT" = "#21908c",
    "DYN_S2_ATBD"    = "#fde725",
    "DYN_S2_DOPT"    = "#d8576b"
  )
  pal_use <- if (!is.null(pal_sc)) c("HOBO obs" = "black", pal_sc) else default_pal

  n_plots <- length(unique(df_all$id_plot[df_all$source != "HOBO obs"]))
  n_days  <- length(unique(df_all$date))

  ggplot2::ggplot(df_all, ggplot2::aes(x = source, y = shift, fill = source)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_boxplot(alpha = 0.72, outlier.size = 0.6, outlier.alpha = 0.4,
                          width = 0.65) +
    ggplot2::scale_fill_manual(values = pal_use, na.value = "grey70") +
    ggplot2::scale_y_continuous(
      breaks = seq(-6, 8, 2),
      labels = function(x) sprintf("%+dh", as.integer(x))
    ) +
    ggplot2::labs(
      title    = "Temporal lag of sub-canopy Tmax vs ERA5 — heatwave days",
      subtitle = sprintf(
        "ERA5 Tmax ≥ %g°C | %d heatwave days | %d HOBO plots",
        heatwave_threshold, n_days, n_plots
      ),
      x    = NULL,
      y    = "Sub-canopy Tmax − ERA5 Tmax lag (h)",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x     = ggplot2::element_text(angle = 20, hjust = 1, size = 10)
    )
}

#' Summary table of timing shifts per source on heatwave days.
#'
#' @title Timing shift summary statistics
#' @description Returns a data.frame with median, mean, sd and IQR of the
#'   hourly Tmax lag for HOBO observations and each MuSICA scenario.
#'   Designed for direct reporting in the manuscript (Table S in annexe or
#'   in-text values for section 3.4/3.5).
#'
#' @param shift_res  List from `compute_tmax_shift_scenarios()`.
#' @return Data.frame with columns `source, n, median_h, mean_h, sd_h, iqr_h`.
summarise_tmax_shift <- function(shift_res) {
  df_obs <- shift_res$hobo_shift %>%
    dplyr::transmute(source = "HOBO_obs", shift)
  df_sim <- shift_res$sim_shift %>%
    dplyr::transmute(source = scenario,   shift)

  dplyr::bind_rows(df_obs, df_sim) %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(
      n        = dplyr::n(),
      median_h = round(median(shift, na.rm = TRUE), 1),
      mean_h   = round(mean(shift, na.rm = TRUE), 1),
      sd_h     = round(sd(shift, na.rm = TRUE), 1),
      iqr_h    = round(IQR(shift, na.rm = TRUE), 1),
      .groups  = "drop"
    ) %>%
    dplyr::arrange(source)
}
