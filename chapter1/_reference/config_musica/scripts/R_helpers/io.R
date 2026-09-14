# ==============================================================================
# Chapter 1 — I/O helpers
# ==============================================================================

#' Force UTC timezone on a NetCDF time variable.
#'
#' @param ncfile   Path to the NetCDF file.
#' @param varname  Name of the time variable (default "time").
#' @return POSIXct vector in UTC.
force_utc_nc <- function(ncfile, varname = "time") {
  if (!file.exists(ncfile)) stop(paste("File not found:", ncfile))
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  time_val  <- ncvar_get(nc, varname)
  units_str <- ncatt_get(nc, varname, "units")$value
  parts     <- strsplit(units_str, " since ")[[1]]
  unit_type <- parts[1]
  origin_str<- gsub("\\s*\\(.*\\)| UTC", "", parts[2])
  origin_dt <- as.POSIXct(origin_str, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  mult <- switch(trimws(unit_type),
                 "hours" = 3600, "minutes" = 60,
                 "seconds" = 1,  "days"    = 86400, 3600)
  with_tz(origin_dt + (time_val * mult), "UTC")
}

#' Load and aggregate the LiDAR rasters.
#'
#' @param in_dir     Directory containing the .tif files.
#' @param agg_factor Aggregation factor (default 2 = 10 m → 20 m).
#' @return List with elements $metrics, $lad, $stack.
load_lidar_rasters <- function(in_dir, agg_factor = 2) {
  r_lai     <- rast(file.path(in_dir, "lai_z1_res_10_m.tif"))
  r_vci     <- rast(file.path(in_dir, "vci_res_10_m.tif"))
  r_hmax    <- rast(file.path(in_dir, "max_res_10_m.tif"))
  r_fcover  <- rast(file.path(in_dir, "fCover_res_10_m.tif"))
  r_lad_st  <- rast(file.path(in_dir, "lad_profiles_z1_res_10_m.tif"))

  # PAD → LAI conversion. The `lai_z1_res_10_m.tif` and `lad_profiles_z1_res_10_m.tif`
  # rasters were generated with lidR's `LAD()` which returns Plant Area Density
  # (PAD). For a spherical leaf angle distribution PAD ≈ LAI/2, so multiplying
  # by 2 recovers LAI. The legacy `main_Blois_detailed.R` pipeline applied this
  # factor when computing `plantareaindex = 2 * sum(allometry$density)`;
  # forgetting it in the refactored pipeline halved the leaf area fed to
  # MuSICA, breaking the buffering signal.
  r_lai    <- 2 * r_lai
  r_lad_st <- 2 * r_lad_st

  # LAI / VCI / fCover : integrated or fractional quantities → MEAN aggregation
  r_lvf_agg <- terra::aggregate(c(r_lai, r_vci, r_fcover),
                                   fact = agg_factor, fun = "mean", na.rm = TRUE)
  # Hmax : canopy ceiling (height of the tallest tree in the cell) → MAX
  # aggregation. MEAN dilutes the dominant trees with canopy gaps in the same
  # cell and systematically under-estimates the canopy top in heterogeneous
  # stands. MAX preserves the local canopy ceiling.
  r_hmax_agg <- terra::aggregate(r_hmax,
                                    fact = agg_factor, fun = "max", na.rm = TRUE)
  r_metrics <- c(r_lvf_agg[[1]], r_lvf_agg[[2]], r_hmax_agg, r_lvf_agg[[3]])
  names(r_metrics) <- c("LAI", "VCI", "Hmax", "fCover")

  r_lad_agg <- terra::aggregate(r_lad_st, fact = agg_factor, fun = "mean", na.rm = TRUE)

  list(metrics = r_metrics, lad = r_lad_agg, stack = c(r_metrics, r_lad_agg))
}

#' Extract daily macroclimate drivers from ERA5 forcing.
#'
#' Computes Tmax, Wind, Rad, VPD, Rain aggregated to daily resolution.
#'
#' @param forcing_file Path to the ERA5 NetCDF forcing file.
#' @param date_seq     Date sequence to filter on.
#' @return Daily tibble with columns: date, Tmax_macro, Wind_mean, Rad_mean, VPD_mean, Rain_sum.
extract_macro_daily <- function(forcing_file, date_seq) {
  nc <- nc_open(forcing_file)

  get_var <- function(v) {
    val <- try(ncvar_get(nc, v), silent = TRUE)
    if (inherits(val, "try-error")) return(rep(NA, nc$dim$time$len))
    return(val)
  }

  macro_tair   <- get_var("Tair")
  macro_qair   <- get_var("Qair")
  macro_psurf  <- get_var("PSurf")
  macro_wind_n <- get_var("Wind_N")
  macro_wind_e <- get_var("Wind_E")
  macro_swdown <- get_var("SWdown")
  macro_rain   <- get_var("Rainf")

  nc_close(nc)
  macro_time <- force_utc_nc(forcing_file, "time")

  t_celsius  <- macro_tair - 273.15
  wind_abs   <- sqrt(macro_wind_n^2 + macro_wind_e^2)

  e_s <- 611.2 * exp((17.67 * t_celsius) / (t_celsius + 243.5))
  e_a <- (macro_qair * macro_psurf) / (0.622 + 0.378 * macro_qair)
  vpd <- pmax(e_s - e_a, 0)

  df_hourly <- data.frame(
    time       = macro_time,
    Tair_macro = t_celsius,
    Wind_macro = wind_abs,
    Rad_macro  = macro_swdown,
    VPD_macro  = vpd,
    Rain_macro = macro_rain
  )

  df_hourly %>%
    filter(as.Date(time) %in% date_seq) %>%
    mutate(date = as.Date(time)) %>%
    group_by(date) %>%
    summarise(
      Tmax_macro = max(Tair_macro, na.rm = TRUE),
      Wind_mean  = mean(Wind_macro, na.rm = TRUE),
      Rad_mean   = mean(Rad_macro, na.rm = TRUE),
      VPD_mean   = mean(VPD_macro, na.rm = TRUE),
      Rain_sum   = sum(Rain_macro * 3600, na.rm = TRUE),
      .groups    = "drop"
    )
}

#' Save a ggplot to disk and return the path invisibly.
#'
#' @param p      ggplot object.
#' @param path   Output file path (PNG).
#' @param width  Width in inches (default 10).
#' @param height Height in inches (default 7).
#' @param dpi    Resolution (default 150).
#' @return Invisible path character.
save_plot <- function(p, path, width = 10, height = 7, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = dpi)

  review_dir <- getOption("fig_review_dir", NULL)
  if (!is.null(review_dir) && grepl("^outputs/", path)) {
    rel  <- sub("^outputs/", "", path)
    dest <- file.path(review_dir, rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    file.copy(path, dest, overwrite = TRUE)
  }

  invisible(path)
}

#' Create the standard outputs/ sub-directory tree.
#'
#' @param root Root directory (default "outputs").
#' @return Invisible root path.
setup_outputs_dirs <- function(root = "outputs") {
  for (sd in c("clusters", "h2", "h1", "fpca", "gamm", "hobo", "s2_annex", "audit",
               "figures", "figures/annex")) {
    dir.create(file.path(root, sd), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(root)
}

#' Read HOBO sub-canopy temperatures and compute daily Delta_obs vs ERA5.
#'
#' @param csv_file      Path to the HOBO temperature CSV.
#' @param date_seq      Date sequence for filtering.
#' @param df_macro      Daily macroclimate dataframe (from extract_macro_daily).
#' @param ids_to_remove Character vector of HOBO plot IDs to exclude.
#' @return Daily tibble with columns: id_plot, date, Tmax_obs, Delta_obs.
read_hobo_daily <- function(csv_file, date_seq, df_macro, ids_to_remove) {
  df_raw <- read.csv(csv_file) %>%
    mutate(
      datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      date     = as.Date(datetime)
    ) %>%
    filter(position_sensor == "a", date %in% date_seq,
           !(id_plot %in% ids_to_remove)) %>%
    group_by(id_plot, date) %>%
    summarise(Tmax_obs = max(t_hobo, na.rm = TRUE), .groups = "drop")

  n_before <- nrow(df_raw)
  df_out <- df_raw %>%
    inner_join(df_macro, by = "date") %>%
    mutate(Delta_obs = Tmax_obs - Tmax_macro)
  n_dropped <- n_before - nrow(df_out)
  if (n_dropped > 0)
    message(sprintf("[read_hobo_daily] %d rows dropped: dates present in HOBO CSV but missing in ERA5 forcing",
                    n_dropped))
  df_out
}
