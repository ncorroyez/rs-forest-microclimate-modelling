# ==============================================================================
# Annex Analysis: Naive Integration of Sentinel-2 and FORMS-H in MuSICA
# Description: Evaluates the limitations of using standard 2D optical data (S2) 
#              and a canopy height model (FORMS-H) to simulate microclimates, 
#              comparing outputs to in-situ HOBO observations.
# ==============================================================================

# ==============================================================================
# 0. CONFIGURATION & LIBRARIES
# ==============================================================================

# ----------------------------- Clear environment ------------------------------
# rm(list=ls(all=TRUE))
# gc()

# ----------------------------- Working Directory ------------------------------
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  cat("Working directory set to:", getwd(), "\n")
}

library(musica.tools)
library(rmusica)
library(terra)
library(sf)
library(tidyverse)
library(lubridate)
library(Metrics)
library(patchwork)
library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(hexbin)

# --- File paths ---
in_dir  <- "in_files"
out_dir <- "out_files/Annexe_S2_FORMS-H"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

macro_nc_file <- file.path(in_dir, "musica_in_Blois.nc")
date_seq      <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
hobo_file     <- file.path(in_dir, "Blois_data_temperature.csv")

# ==============================================================================
# 1. DATA LOADING & EXTRACTION AT HOBO PLOTS
# ==============================================================================

r_s2_lai         <- rast(file.path(in_dir, "s2lai_summer_atbd_res_10_m.tif"))
r_formsh         <- rast(file.path(in_dir, "FORMS-H_Blois.tif"))
r_formsh_aligned <- resample(r_formsh, r_s2_lai, method = "bilinear")

# Aggregate Sentinel-2 LAI from 10m to 20m
# fact = 2 means 2x2 pixels (10m -> 20m). fun = "mean" averages the LAI.
cat("Aggregating S2 LAI to 20m...\n")
r_s2_lai_20m <- terra::aggregate(r_s2_lai, 
                                 fact = 2, fun = "mean", na.rm = TRUE)

# 2. Resample FORMS-H to exactly match the new 20m S2 grid
cat("Aligning FORMS-H to the new 20m grid...\n")
r_formsh_20m <- resample(r_formsh, r_s2_lai_20m, method = "bilinear")

# 3. Stack the required rasters
r_optical_stack <- c(r_s2_lai_20m, r_formsh_20m)
names(r_optical_stack) <- c("LAI_S2", "Hmax_FORMSH")

# Stack the required rasters (fCover is omitted as we assume it is 1.0)
# r_optical_stack <- c(r_s2_lai, r_formsh_aligned)
# names(r_optical_stack) <- c("LAI_S2", "Hmax_FORMSH")

# Load HOBO plots and extract values
hobo_pts <- st_read(file.path(in_dir, "data_Blois_utm31n.geojson"), 
                    quiet = TRUE)
df_hobo_extract <- terra::extract(r_optical_stack, 
                                  vect(hobo_pts), xy = TRUE)
ids_to_remove <- c("41_13", "41_14", "41_20", "41_41", 
                   "41_50", "41_51", "41_53")

df_naive_inputs <- hobo_pts %>%
  st_drop_geometry() %>%
  dplyr::select(id_plot) %>%
  bind_cols(df_hobo_extract) %>%
  filter(!(id_plot %in% ids_to_remove)) %>%
  drop_na(LAI_S2, Hmax_FORMSH) 

# ==============================================================================
# 2. NAIVE LAD PROFILE GENERATION
# ==============================================================================

#' Generate a naive uniform LAD profile based on total LAI and Hmax
generate_naive_lad <- function(total_lai, hmax) {
  max_layer <- ceiling(hmax)
  if (max_layer <= 0 || is.na(max_layer)) return(data.frame(height = 1, 
                                                            density = 0))
  
  canopy_base  <- 2
  canopy_depth <- max(max_layer - canopy_base, 1)
  lad_value    <- total_lai / canopy_depth
  
  heights   <- 1:max_layer
  densities <- ifelse(heights >= canopy_base, lad_value, 0)
  
  data.frame(height = heights, density = densities)
}

# ==============================================================================
# 3. MuSICA SIMULATIONS WITH S2 + FORMS-H
# ==============================================================================

# musica.cmd  <- "bash -i -c musica"
# total_plots <- nrow(df_naive_inputs)
# cat(sprintf("Starting MuSICA batch for %d plots using Naive inputs...\n", 
#             total_plots))
# 
# for (i in seq_len(total_plots)) {
#   current_plot <- df_naive_inputs[i, ]
#   sim_id       <- sprintf("HOBO_NAIVE_%s", current_plot$id_plot)
#   history_file <- file.path(out_dir, 
#                             "S2_20m",
#                             paste0("musica_out_", sim_id, ".nc"))
#   
#   if (file.exists(history_file)) next
#   
#   cat(sprintf("Processing %d/%d: %s\n", i, total_plots, sim_id))
#   
#   canopy_height_top <- current_plot$Hmax_FORMSH
#   plantareaindex    <- current_plot$LAI_S2
#   clumping          <- 1 
#   
#   allometry <- generate_naive_lad(plantareaindex, canopy_height_top)
#   
#   phenology <- calc_phenology(
#     list.year             = c(2021, 2022),
#     nleafage              = 1,
#     budburst_date         = 115,
#     leaf_age_max_in       = 0.56,
#     relative_age_firstmax = 0.10,
#     relative_age_lastmax  = 0.75,
#     LAI_max_per_cohort    = plantareaindex
#   )
#   
#   tryCatch({
#     callmusica(
#       musica.param = list(
#         "setupctl" = list(
#           "clumping_factor"  = clumping,
#           "forcing_filename" = macro_nc_file,
#           "history_filename" = history_file,
#           "forcing_height"   = canopy_height_top + 2
#         )
#       ),
#       leaf.param = list(
#         "musica_veg1" = list(
#           "phenology"        = phenology,
#           "allometry"        = allometry,
#           "leafphenologyctl" = list("lai_max_per_cohort" = plantareaindex),
#           "leafallometryctl" = list(
#             "canopy_height_top"    = canopy_height_top,
#             "canopy_height_bottom" = 1
#           ),
#           "leafmusicactl" = list("canopy_height_top" = canopy_height_top)
#         )
#       ),
#       musica.cmd = musica.cmd,
#       keep.tmp   = FALSE,
#       out.netcdf = TRUE
#     )
#   }, error = function(e) {
#     cat(sprintf("ERROR on %s: %s\n", sim_id, e$message))
#   })
# }

# ==============================================================================
# 4. POST-PROCESSING & VALIDATION
# ==============================================================================

# cat("\nExtracting Macroclimate (ERA5)...\n")
# # Helper function for correct UTC time mapping
# force_utc_nc <- function(ncfile, varname = "time") {
#   nc        <- nc_open(ncfile)
#   on.exit(nc_close(nc))
#   time_val  <- ncvar_get(nc, varname)
#   units_str <- ncatt_get(nc, varname, "units")$value
#   parts     <- strsplit(units_str, " since ")[[1]]
#   origin_dt <- as.POSIXct(gsub("\\s*\\(.*\\)| UTC", "", parts[2]), tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
#   mult      <- switch(trimws(parts[1]), "hours" = 3600, "minutes" = 60, "seconds" = 1, "days" = 86400, 3600)
#   with_tz(origin_dt + (time_val * mult), "UTC")
# }
# 
# nc_force   <- nc_open(macro_nc_file)
# macro_tair <- ncvar_get(nc_force, "Tair")
# nc_close(nc_force)
# 
# df_macro <- data.frame(time = force_utc_nc(macro_nc_file), Tair_macro = macro_tair - 273.15) %>%
#   filter(as.Date(time) %in% date_seq) %>%
#   mutate(date = as.Date(time)) %>%
#   group_by(date) %>%
#   summarise(Tmax_macro = max(Tair_macro, na.rm = TRUE), .groups = "drop")
# 
# cat("Processing In-Situ HOBO observations...\n")
# df_hobo_temp <- read.csv(hobo_file) %>%
#   mutate(
#     datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
#     date     = as.Date(datetime)
#   ) %>%
#   filter(position_sensor == "a", date %in% date_seq) %>%
#   group_by(id_plot, date) %>%
#   summarise(Tmax_obs = max(t_hobo, na.rm = TRUE), .groups = "drop") %>%
#   inner_join(df_macro, by = "date") %>%
#   mutate(Delta_obs = Tmax_obs - Tmax_macro)
# 
# cat("Extracting MuSICA Naive Simulations...\n")
# nc_files_naive <- list.files(out_dir, pattern = "\\.nc$", full.names = TRUE)
# 
# df_musica_naive <- map_df(nc_files_naive, function(f) {
#   filename <- basename(f)
#   id_str   <- str_extract(filename, "(?<=HOBO_NAIVE_).*(?=\\.nc)")
#   
#   nc <- try(nc_open(f), silent = TRUE)
#   if (inherits(nc, "try-error")) return(NULL)
#   raw_data <- try(get_variable(nc, "Tair_z"), silent = TRUE)
#   nc_close(nc)
#   if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
#   
#   raw_data %>%
#     filter(nair == 1) %>%
#     mutate(Tair_sim = Tair_z - 273.15) %>%
#     filter(as.Date(time) %in% date_seq) %>%
#     mutate(time = time - hours(2), date = as.Date(time)) %>%
#     group_by(date) %>%
#     summarise(Tmax_naive = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
#     mutate(id_plot = id_str)
# })
# 
# # Merge and calculate naive Delta Tmax
# df_eval_naive <- df_musica_naive %>%
#   inner_join(df_macro, by = "date") %>%
#   mutate(Delta_naive = Tmax_naive - Tmax_macro) %>%
#   # Drop Tmax_macro from df_hobo_temp before joining to prevent duplication
#   inner_join(df_hobo_temp %>% dplyr::select(-Tmax_macro), by = c("id_plot", "date"))
# 
# # ==============================================================================
# # 4.B EXTENDED VALIDATION & VISUALIZATIONS
# # ==============================================================================
# 
# # --- 1. Calculate Evaluation Metrics ---
# 
# # Helper function to compute linear regression slope
# calc_slope <- function(y, x) {
#   model <- lm(y ~ x)
#   return(coef(model)[2])
# }
# 
# # A. Metrics for Absolute Tmicro
# r2_tmicro    <- cor(df_eval_naive$Tmax_obs, df_eval_naive$Tmax_naive, use = "complete.obs")^2
# rmse_tmicro  <- rmse(df_eval_naive$Tmax_obs, df_eval_naive$Tmax_naive)
# bias_tmicro  <- mean(df_eval_naive$Tmax_naive - df_eval_naive$Tmax_obs, na.rm = TRUE)
# 
# # B. Metrics for Delta Tmax (Cooling Effect)
# r2_delta    <- cor(df_eval_naive$Delta_obs, df_eval_naive$Delta_naive, use = "complete.obs")^2
# rmse_delta  <- rmse(df_eval_naive$Delta_obs, df_eval_naive$Delta_naive)
# bias_delta  <- mean(df_eval_naive$Delta_naive - df_eval_naive$Delta_obs, na.rm = TRUE)
# 
# # C. Metrics for Thermal Buffering (Slope of Tmicro ~ Tmacro)
# slope_obs   <- calc_slope(df_eval_naive$Tmax_obs, df_eval_naive$Tmax_macro)
# slope_naive <- calc_slope(df_eval_naive$Tmax_naive, df_eval_naive$Tmax_macro)
# 
# log_slope_obs   <- log(slope_obs)
# log_slope_naive <- log(slope_naive)
# 
# # Print summary to console
# cat("--------------------------------------------------\n")
# cat("--- EXTENDED PERFORMANCE METRICS (NAIVE MODEL) ---\n")
# cat("--------------------------------------------------\n")
# cat(sprintf("T_MICRO   | R²: %.2f | RMSE: %.2f | Bias: %.2f\n", r2_tmicro, rmse_tmicro, bias_tmicro))
# cat(sprintf("DELTA_MAX | R²: %.2f | RMSE: %.2f | Bias: %.2f\n", r2_delta, rmse_delta, bias_delta))
# cat(sprintf("SLOPES    | Obs: %.2f (log: %.2f) | Naive: %.2f (log: %.2f)\n", 
#             slope_obs, log_slope_obs, slope_naive, log_slope_naive))
# cat("--------------------------------------------------\n")
# 
# # --- 2. Build the Plots ---
# 
# # Common theme for all panels
# theme_eval <- theme_bw(base_size = 12) +
#   theme(plot.title = element_text(face = "bold", size = 12),
#         plot.subtitle = element_text(size = 10, color = "grey30"))
# 
# # Panel A: Absolute Tmicro
# p_tmicro <- ggplot(df_eval_naive, aes(x = Tmax_obs, y = Tmax_naive)) +
#   geom_point(alpha = 0.3, color = "#31688e") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
#   geom_smooth(method = "lm", color = "red", se = FALSE) +
#   labs(
#     title = "A. Absolute Sub-Canopy Temperature",
#     subtitle = sprintf("R² = %.2f | RMSE = %.2f °C\nBias = %.2f °C", 
#                        r2_tmicro, rmse_tmicro, bias_tmicro),
#     x = "Observed Tmicro (°C)",
#     y = "Simulated Tmicro (°C)"
#   ) +
#   theme_eval
# 
# # Panel B: Delta Tmax (Cooling Effect)
# p_delta <- ggplot(df_eval_naive, aes(x = Delta_obs, y = Delta_naive)) +
#   geom_point(alpha = 0.3, color = "#d8576b") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
#   geom_smooth(method = "lm", color = "red", se = FALSE) +
#   labs(
#     title = "B. Canopy Cooling Effect (\U0394Tmax)",
#     subtitle = sprintf("R² = %.2f | RMSE = %.2f °C\nBias = %.2f °C", 
#                        r2_delta, rmse_delta, bias_delta),
#     x = "Observed \U0394Tmax (°C)",
#     y = "Simulated \U0394Tmax (°C)"
#   ) +
#   theme_eval
# 
# # Prepare data for Panel C
# df_plot_slope <- df_eval_naive %>%
#   dplyr::select(date, id_plot, Tmax_macro, Tmax_obs, Tmax_naive) %>%
#   pivot_longer(
#     cols = c(Tmax_obs, Tmax_naive),
#     names_to = "Data_Source",
#     values_to = "Tmax_micro"
#   ) %>%
#   mutate(
#     Data_Source = recode_factor(Data_Source,
#                                 "Tmax_obs"   = "Observed (HOBO)",
#                                 "Tmax_naive" = "Simulated (Naive S2)")
#   )
# 
# # Panel C: Thermal Buffering (Tmicro ~ Tmacro)
# p_buffering <- ggplot(df_plot_slope, aes(x = Tmax_macro, y = Tmax_micro, color = Data_Source)) +
#   geom_point(alpha = 0.2, size = 1.5) +
#   # 1:1 Line representing no forest buffering
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
#   geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
#   scale_color_manual(values = c("Observed (HOBO)" = "#31688e", 
#                                 "Simulated (Naive S2)" = "#d8576b")) +
#   labs(
#     title = "C. Thermal Buffering (Tmicro ~ Tmacro)",
#     subtitle = sprintf("Slope Obs = %.2f (log: %.2f) | Slope Sim = %.2f (log: %.2f)", 
#                        slope_obs, log_slope_obs, slope_naive, log_slope_naive),
#     x = "Macroclimate Tmax (ERA5, °C)",
#     y = "Sub-Canopy Tmax (°C)",
#     color = NULL
#   ) +
#   theme_eval +
#   theme(legend.position = c(0.2, 0.85),
#         legend.background = element_rect(fill = "white", color = "black"))
# 
# # --- 3. Combine Plots with Patchwork ---
# 
# combined_plot <- (p_tmicro | p_delta) / p_buffering +
#   plot_annotation(
#     title = "Naive Optical Model vs. Reality: The Missing 3D Structure",
#     theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
#   )
# 
# print(combined_plot)





























# ==============================================================================
# 4. POST-PROCESSING & VALIDATION (HOURLY SLOPES + DAILY MAX)
# ==============================================================================

# --- 1. Extract Hourly Macroclimate (ERA5) using rmusica ---
cat("\nExtracting Hourly Macroclimate (ERA5)...\n")

nc_force <- nc_open(macro_nc_file)
raw_macro <- get_variable(nc_force, "Tair")
nc_close(nc_force)

df_macro_hourly <- raw_macro %>%
  mutate(
    # Round to nearest hour to align with HOBO
    datetime = round_date(time, "hour"), 
    Tair = Tair - 273.15
  ) %>%
  filter(as.Date(datetime) %in% date_seq) %>%
  # Aggregate in case the forcing is at 30-min intervals
  group_by(datetime) %>%
  summarise(Tair_macro = mean(Tair, na.rm = TRUE), .groups = "drop")

# --- 2. Extract Hourly In-Situ HOBO observations ---
cat("Processing Hourly In-Situ HOBO observations...\n")
df_hobo_hourly <- read.csv(hobo_file) %>%
  mutate(
    # Ensure HOBO timestamps are strictly rounded to the hour as well
    datetime = round_date(as.POSIXct(datetime, 
                                     format = "%Y-%m-%d %H:%M:%S", 
                                     tz = "UTC"), 
                          "hour")
  ) %>%
  filter(position_sensor == "a", as.Date(datetime) %in% date_seq) %>%
  dplyr::select(id_plot, datetime, t_hobo) %>%
  rename(Tair_obs = t_hobo) %>%
  group_by(id_plot, datetime) %>%
  summarise(Tair_obs = mean(Tair_obs, na.rm = TRUE), .groups = "drop")

# --- 3. Function to Extract Hourly MuSICA outputs with Cache ---
extract_musica_hourly <- function(nc_dir, id_regex, rds_name) {
  rds_file <- file.path(nc_dir, rds_name)
  
  if (file.exists(rds_file)) {
    cat(sprintf("-> Cached data found. Loading %s...\n", rds_name))
    return(readRDS(rds_file))
  }
  
  cat(sprintf("-> Extracting Hourly NetCDF from %s...\n", nc_dir))
  nc_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) stop(paste("No .nc files found in", nc_dir))
  
  df_extracted <- map_df(nc_files, function(f) {
    id_str <- str_extract(basename(f), id_regex)
    nc <- try(nc_open(f), silent = TRUE)
    if (inherits(nc, "try-error")) return(NULL)
    raw_data <- try(get_variable(nc, "Tair_z"), silent = TRUE)
    nc_close(nc)
    if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
    
    raw_data %>%
      filter(nair == 1) %>%
      mutate(
        Tair_sim = Tair_z - 273.15,
        # Adjust UTC and round simulated time to the nearest hour
        datetime = round_date(time - hours(2), "hour") 
      ) %>%
      filter(as.Date(datetime) %in% date_seq) %>%
      group_by(datetime) %>%
      summarise(Tair_sim = mean(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
      mutate(id_plot = id_str)
  })
  
  saveRDS(df_extracted, rds_file)
  return(df_extracted)
}

# --- 4. Load & Merge HOURLY Data ---
dir_s2    <- "out_files/Annexe_S2_FORMS-H/S2_20m"
dir_lidar <- "out_files/Annexe_S2_FORMS-H/LiDAR_20m"

df_sim_s2_hr    <- extract_musica_hourly(dir_s2, 
                                         "(?<=HOBO_NAIVE_).*(?=\\.nc)", "extracted_s2_hourly_aligned.rds") %>% rename(Tair_S2 = Tair_sim)
df_sim_lidar_hr <- extract_musica_hourly(dir_lidar, 
                                         "(?<=HOBO_).*(?=\\.nc)", "extracted_lidar_hourly_aligned.rds") %>% rename(Tair_LiDAR = Tair_sim)

df_hourly_eval <- df_hobo_hourly %>%
  inner_join(df_macro_hourly, by = "datetime") %>%
  inner_join(df_sim_s2_hr, by = c("id_plot", "datetime")) %>%
  inner_join(df_sim_lidar_hr, by = c("id_plot", "datetime"))

cat(sprintf("Time alignment successful! Total overlapping hourly records: %d\n",
            nrow(df_hourly_eval)))

# --- 5. Aggregate to DAILY Data (for Delta Tmax and Absolute Tmax plots) ---
df_daily_eval <- df_hourly_eval %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(id_plot, date) %>%
  summarise(
    Tmax_obs   = max(Tair_obs, na.rm = TRUE),
    Tmax_macro = max(Tair_macro, na.rm = TRUE),
    Tmax_S2    = max(Tair_S2, na.rm = TRUE),
    Tmax_LiDAR = max(Tair_LiDAR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Delta_obs   = Tmax_obs - Tmax_macro,
    Delta_S2    = Tmax_S2 - Tmax_macro,
    Delta_LiDAR = Tmax_LiDAR - Tmax_macro
  )

# ==============================================================================
# 4.B METRICS CALCULATION
# ==============================================================================

# A) Daily Metrics (R2 & RMSE on Tmax and DeltaTmax)
r2_t_s2   <- cor(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_S2, 
                 use = "complete.obs")^2
rmse_t_s2 <- rmse(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_S2)
r2_d_s2   <- cor(df_daily_eval$Delta_obs, df_daily_eval$Delta_S2, 
                 use = "complete.obs")^2
rmse_d_s2 <- rmse(df_daily_eval$Delta_obs, df_daily_eval$Delta_S2)

r2_t_ld   <- cor(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_LiDAR, 
                 use = "complete.obs")^2
rmse_t_ld <- rmse(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_LiDAR)
r2_d_ld   <- cor(df_daily_eval$Delta_obs, df_daily_eval$Delta_LiDAR, 
                 use = "complete.obs")^2
rmse_d_ld <- rmse(df_daily_eval$Delta_obs, df_daily_eval$Delta_LiDAR)

# B) Hourly Metrics (Thermal Buffering Slopes)
calc_slope <- function(y, x) { coef(lm(y ~ x))[2] }

slope_obs <- calc_slope(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_macro)
slope_s2  <- calc_slope(df_hourly_eval$Tair_S2, df_hourly_eval$Tair_macro)
slope_ld  <- calc_slope(df_hourly_eval$Tair_LiDAR, df_hourly_eval$Tair_macro)


# ==============================================================================
# 4.C VISUALIZATION: 4-PANEL COMPARISON
# ==============================================================================

# --- 1. Additional Metrics for Hourly Tmicro (Panel A) ---
r2_hr_s2   <- cor(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_S2, use = "complete.obs")^2
rmse_hr_s2 <- rmse(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_S2)
r2_hr_ld   <- cor(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_LiDAR, use = "complete.obs")^2
rmse_hr_ld <- rmse(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_LiDAR)

# --- 2. Aesthetics ---
theme_eval <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 16),
        legend.position = "bottom", legend.title = element_blank())

colors_model <- c("Naive S2 + FORMS-H" = "#d8576b", 
                  "Detailed LiDAR 3D" = "#35b779", 
                  "Observed (HOBO)" = "#31688e")

# --- 3. Prepare Data for Panel A (Hourly Tmicro Scatter) ---
df_plot_hourly_scatter <- df_hourly_eval %>%
  dplyr::select(id_plot, datetime, Tair_obs, Tair_S2, Tair_LiDAR) %>%
  pivot_longer(cols = c(Tair_S2, Tair_LiDAR), names_to = "Model", values_to = "Tair_sim") %>%
  mutate(Model = recode_factor(Model, 
                               "Tair_S2" = "Naive S2 + FORMS-H", 
                               "Tair_LiDAR" = "Detailed LiDAR 3D"))

# Panel A: Absolute Tmicro (Hourly)
p_tmicro_hr <- ggplot(df_plot_hourly_scatter, aes(x = Tair_obs, y = Tair_sim, color = Model)) +
  geom_point(alpha = 0.05, size = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = colors_model) +
  labs(
    title = "A. Hourly Sub-Canopy Temperature (Tmicro)",
    subtitle = sprintf("LiDAR: R²=%.2f, RMSE=%.2f\nS2: R²=%.2f, RMSE=%.2f", 
                       r2_hr_ld, rmse_hr_ld, r2_hr_s2, rmse_hr_s2),
    x = "Observed Hourly Tmicro (°C)", y = "Simulated Hourly Tmicro (°C)"
  ) + theme_eval + theme(legend.position = "none")


# --- 4. Prepare Data for Panels B & C (Daily Max Data) ---
df_plot_daily <- df_daily_eval %>%
  dplyr::select(id_plot, date, Tmax_obs, Delta_obs, Tmax_S2, Delta_S2, Tmax_LiDAR, Delta_LiDAR) %>%
  pivot_longer(cols = c(Tmax_S2, Tmax_LiDAR), names_to = "Model", values_to = "Tmax_sim") %>%
  mutate(
    Delta_sim = ifelse(Model == "Tmax_S2", Delta_S2, Delta_LiDAR),
    Model     = recode_factor(Model, 
                              "Tmax_S2" = "Naive S2 + FORMS-H", 
                              "Tmax_LiDAR" = "Detailed LiDAR 3D")
  )

# Panel B: Absolute Tmicro_max (Daily Max)
p_tmicro_max <- ggplot(df_plot_daily, aes(x = Tmax_obs, y = Tmax_sim, color = Model)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = colors_model) +
  labs(
    title = "B. Daily Max Sub-Canopy Temp (Tmicro_max)",
    subtitle = sprintf("LiDAR: R²=%.2f, RMSE=%.2f\nS2: R²=%.2f, RMSE=%.2f", 
                       r2_t_ld, rmse_t_ld, r2_t_s2, rmse_t_s2),
    x = "Observed Tmax (°C)", y = "Simulated Tmax (°C)"
  ) + theme_eval + theme(legend.position = "none")

# Panel C: Delta Tmax (Cooling Effect)
p_delta <- ggplot(df_plot_daily, aes(x = Delta_obs, y = Delta_sim, color = Model)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = colors_model) +
  labs(
    title = "C. Canopy Cooling Effect (\U0394Tmax)",
    subtitle = sprintf("LiDAR: R²=%.2f, RMSE=%.2f\nS2: R²=%.2f, RMSE=%.2f", 
                       r2_d_ld, rmse_d_ld, r2_d_s2, rmse_d_s2),
    x = "Observed \U0394Tmax (°C)", y = "Simulated \U0394Tmax (°C)"
  ) + theme_eval + theme(legend.position = "none")


# --- 5. Prepare Data for Panel D (log(slope) buffering) ---
df_slopes <- data.frame(
  Model = factor(c("Observed (HOBO)", "Detailed LiDAR 3D", "Naive S2 + FORMS-H"),
                 levels = c("Observed (HOBO)", "Detailed LiDAR 3D", "Naive S2 + FORMS-H")),
  Log_Slope = c(log(slope_obs), log(slope_ld), log(slope_s2))
)

# Panel D: log(slope) Barplot
p_slope <- ggplot(df_slopes, aes(x = Model, y = Log_Slope, fill = Model)) +
  geom_col(color = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", Log_Slope)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = colors_model) +
  labs(
    title = "D. Thermal Buffering Capacity",
    subtitle = "log(slope) of Tmicro ~ Tmacro (Hourly)",
    x = NULL, y = "log(Slope)"
  ) + theme_eval + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 15, hjust = 1))

# --- 6. Combine Plots with Patchwork (2x2 Grid) ---
combined_plot <- (p_tmicro_hr | p_tmicro_max) / (p_delta | p_slope) +
  plot_annotation(
    title = "Impact of 3D Canopy Structure on Microclimate Simulations",
    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))
  ) + plot_layout(guides = "collect") & theme(legend.position = "none")

print(combined_plot)







































# ==============================================================================
# 4.C VISUALIZATION: 7-PANEL HEXBIN COMPARISON (S2 vs LiDAR)
# ==============================================================================

# --- 1. Additional Metrics for Hourly Tmicro ---
r2_hr_s2   <- cor(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_S2, use = "complete.obs")^2
rmse_hr_s2 <- rmse(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_S2)
bias_hr_s2 <- mean(df_hourly_eval$Tair_S2 - df_hourly_eval$Tair_obs, na.rm = TRUE)

r2_hr_ld   <- cor(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_LiDAR, use = "complete.obs")^2
rmse_hr_ld <- rmse(df_hourly_eval$Tair_obs, df_hourly_eval$Tair_LiDAR)
bias_hr_ld <- mean(df_hourly_eval$Tair_LiDAR - df_hourly_eval$Tair_obs, na.rm = TRUE)

# --- 1.B Compute Performance Metrics for Daily Max (Tmicro_max) ---
r2_t_s2   <- cor(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_S2, use = "complete.obs")^2
rmse_t_s2 <- rmse(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_S2)
bias_t_s2 <- mean(df_daily_eval$Tmax_S2 - df_daily_eval$Tmax_obs, na.rm = TRUE)

r2_t_ld   <- cor(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_LiDAR, use = "complete.obs")^2
rmse_t_ld <- rmse(df_daily_eval$Tmax_obs, df_daily_eval$Tmax_LiDAR)
bias_t_ld <- mean(df_daily_eval$Tmax_LiDAR - df_daily_eval$Tmax_obs, na.rm = TRUE)

# --- 1.C Compute Performance Metrics for Cooling Effect (Delta Tmax) ---
r2_d_s2   <- cor(df_daily_eval$Delta_obs, df_daily_eval$Delta_S2, use = "complete.obs")^2
rmse_d_s2 <- rmse(df_daily_eval$Delta_obs, df_daily_eval$Delta_S2)
bias_d_s2 <- mean(df_daily_eval$Delta_S2 - df_daily_eval$Delta_obs, na.rm = TRUE)

r2_d_ld   <- cor(df_daily_eval$Delta_obs, df_daily_eval$Delta_LiDAR, use = "complete.obs")^2
rmse_d_ld <- rmse(df_daily_eval$Delta_obs, df_daily_eval$Delta_LiDAR)
bias_d_ld <- mean(df_daily_eval$Delta_LiDAR - df_daily_eval$Delta_obs, na.rm = TRUE)

# --- 2. Aesthetics ---
theme_hex <- theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right",
        legend.key.width = unit(0.3, "cm"),
        plot.margin = margin(5, 5, 5, 5))

lims_hr    <- c(0, 40)
lims_max   <- c(15, 42)
lims_delta <- c(-6, 12)

# Function to add performance metrics inside the plot
add_metrics <- function(r2, rmse, bias) {
  annotate("text", x = -Inf, y = Inf, 
           label = sprintf("R² = %.2f\nRMSE = %.2f\nBias = %.2f", r2, rmse, bias),
           hjust = -0.1, vjust = 1.1, fontface = "bold", size = 5)
}
add_slope <- function(slope_val) {
  annotate("text", x = Inf, y = -Inf, 
           label = sprintf("log(Slope) = %.3f ", log(slope_val)),
           hjust = 1.1, vjust = -0.5, fontface = "italic", size = 4.5, color = "black")
}

# --- 3. Panel A: Hourly Tmicro ---
p_a1 <- ggplot(df_hourly_eval, aes(x = Tair_obs, y = Tair_S2)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "magma", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#d8576b", linewidth = 1.2) +
  add_metrics(r2_hr_s2, rmse_hr_s2, bias_hr_s2) +
  add_slope(slope_s2) +
  coord_fixed(ratio = 1, xlim = lims_hr, ylim = lims_hr) +
  labs(title = "A1. Hourly Tmicro (MuSICA-S2)", x = NULL, y = "Simulated") + theme_hex

p_a2 <- ggplot(df_hourly_eval, aes(x = Tair_obs, y = Tair_LiDAR)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "mako", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#35b779", linewidth = 1.2) +
  add_metrics(r2_hr_ld, rmse_hr_ld, bias_hr_ld) +
  add_slope(slope_ld) +
  coord_fixed(ratio = 1, xlim = lims_hr, ylim = lims_hr) +
  labs(title = "A2. Hourly Tmicro (MuSICA-LiDAR)", x = NULL, y = NULL) + theme_hex


# --- 4. Panel B: Daily Tmicro_max ---
p_b1 <- ggplot(df_daily_eval, aes(x = Tmax_obs, y = Tmax_S2)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(option = "magma", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#d8576b", linewidth = 1.2) +
  add_metrics(r2_t_s2, rmse_t_s2, bias_t_s2) +
  coord_fixed(ratio = 1, xlim = lims_max, ylim = lims_max) +
  labs(title = "B1. Daily Tmax (MuSICA-S2)", x = NULL, y = "Simulated") + theme_hex

p_b2 <- ggplot(df_daily_eval, aes(x = Tmax_obs, y = Tmax_LiDAR)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(option = "mako", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#35b779", linewidth = 1.2) +
  add_metrics(r2_t_ld, rmse_t_ld, bias_t_ld) + 
  coord_fixed(ratio = 1, xlim = lims_max, ylim = lims_max) +
  labs(title = "B2. Daily Tmax (MuSICA-LiDAR)", x = NULL, y = NULL) + theme_hex


# --- 5. Panel C: Canopy Cooling Effect (Delta Tmax) ---
p_c1 <- ggplot(df_daily_eval, aes(x = Delta_obs, y = Delta_S2)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(option = "magma", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#d8576b", linewidth = 1.2) +
  add_metrics(r2_d_s2, rmse_d_s2, bias_d_s2) + 
  coord_fixed(ratio = 1, xlim = lims_delta, ylim = lims_delta) +
  labs(title = "C1. ΔTmax (MuSICA-S2)", x = "HOBOs ΔTmax (°C)", y = "Simulated") + theme_hex

p_c2 <- ggplot(df_daily_eval, aes(x = Delta_obs, y = Delta_LiDAR)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(option = "mako", trans = "log10") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#35b779", linewidth = 1.2) +
  add_metrics(r2_d_ld, rmse_d_ld, bias_d_ld) + 
  coord_fixed(ratio = 1, xlim = lims_delta, ylim = lims_delta) +
  labs(title = "C2. ΔTmax (MuSICA-LiDAR)", x = "HOBOs ΔTmax (°C)", y = NULL) + theme_hex


# --- 6. Panel D: Thermal Buffering (log slope) ---
df_slopes <- data.frame(
  Model = factor(c("HOBOs", "LiDAR", "S2 + FORMS-H"),
                 levels = c("HOBOs", "LiDAR", "S2 + FORMS-H")),
  Log_Slope = c(log(slope_obs), log(slope_ld), log(slope_s2))
)

p_slope <- ggplot(df_slopes, aes(x = Model, y = Log_Slope, fill = Model)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Log_Slope)), 
            vjust = ifelse(df_slopes$Log_Slope < 0, 1.2, -0.5), size = 4, fontface = "bold") +
  scale_fill_manual(values = colors_model) +
  labs(
    title = "D. Thermal Buffering Capacity",
    subtitle = "log(slope) of Hourly Tmicro ~ Tmacro",
    x = NULL, y = "log(Slope)"
  ) + theme_bw(base_size = 16) + 
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 16))


# --- 7. Combine Plots with Patchwork ---
combined_hex_plot <- (p_a1 | p_a2) / (p_b1 | p_b2) / (p_c1 | p_c2) +
  plot_annotation(
    title = "Impact of 3D Canopy Structure on Microclimate Simulations",
    subtitle = paste("Observed buffering log(slope) =", sprintf("%.3f", log(slope_obs))), # <--- Ajout de l'obs ici
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic", color = "grey30"))
  )

print(combined_hex_plot)