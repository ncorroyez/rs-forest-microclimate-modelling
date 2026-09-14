# ==============================================================================
# SCRIPT: Sensitivity Analysis - LiDAR Radius & MuSICA Execution
# OUTPUT: /out_files/radius_test2/{radius}m/musica_out_Blois_pt_{id}_radius_{r}.nc
# ==============================================================================

rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(musica.tools)
library(rmusica)
library(lidR)
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)
library(purrr)
library(tidyverse)

# --- 1. CONFIGURATION ---

ctg_path <- "/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm"
geojson_path <- "in_files/data_Blois_utm31n.geojson"
forcing_file <- "./in_files/musica_in_Blois.nc" 
musica.cmd <- "bash -i -c musica"
base_out_dir <- "./out_files/radius_test_square"
hobo_file     <- "in_files/Blois_data_temperature.csv"
lidar_file    <- "in_files/allometry_profiles.csv"
metrics_file  <- "in_files/metrics_results_25.csv"
date_seq      <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
ids_to_remove <- c("41_13", "41_14", "41_20", 
                   "41_41", "41_50", "41_51", "41_53")

# Load catalog and geometry
ctg <- readLAScatalog(ctg_path)
plots_sf <- st_read(geojson_path)

# Extract plot coordinates
df_plots <- plots_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(x = X, y = Y) %>%
  mutate(id_plot = plots_sf$id_plot)

# Define test radii
radii <- c(5, 10, 15, 20, 25, 30, 50) 
all_metrics_list <- list()

# --- 2. PROCESSING & MODEL LOOP ---

for (radius in radii) {
  message(sprintf("\n========== RADIUS: %s m ==========", radius))
  
  # Create radius-specific output directory
  radius_dir <- file.path(base_out_dir, paste0(radius, "m"))
  if (!dir.exists(radius_dir)) dir.create(radius_dir, recursive = TRUE)
  
  for (i in seq_len(nrow(df_plots))) {
    current_id <- df_plots$id_plot[i]
    message(sprintf("Processing Plot: %s (Radius: %sm)", current_id, radius))
    
    # Define output file path
    history_file <- file.path(radius_dir, paste0("musica_out_Blois_pt_", current_id, "_radius_", radius, ".nc"))
    
    # Skip simulation if file already exists
    if (file.exists(history_file)) {
      message(" -> File exists, skipping MuSICA run.")
      next 
    }
    
    # A. LiDAR Extraction & Normalization
    # pt <- clip_circle(ctg, df_plots$x[i], df_plots$y[i], radius)
    x_center <- df_plots$x[i]
    y_center <- df_plots$y[i]
    
    # Clip a square: x +/- radius, y +/- radius
    # For a radius of 10, this creates a 20x20 bounding box
    pt <- clip_rectangle(
      ctg, 
      xleft   = x_center - radius, 
      ybottom = y_center - radius, 
      xright  = x_center + radius, 
      ytop    = y_center + radius
    )
    
    if (is.null(pt) || npoints(pt) < 10) {
      message(" -> Not enough LiDAR points, skipping.")
      next
    }
    
    # Use rasterize_* functions for lidR >= 4.0.0
    dtm <- rasterize_terrain(pt, res = 1, algorithm = tin())
    pt_norm <- normalize_height(pt, dtm)
    
    # B. Structural Parameters
    # 1. LAD & PAI
    lad_data <- LAD(pt_norm$Z, z0 = 1) 
    allometry <- data.frame(height = lad_data$z, density = lad_data$lad)
    pai_val <- round(2 * sum(allometry$density, na.rm = TRUE), 2)
    
    # 2. Canopy Heights
    h_max <- max(pt_norm$Z, na.rm = TRUE)
    
    # 3. Clumping (fCover at 30% Hmax)
    chm <- rasterize_canopy(pt_norm, res = 0.5, pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1)))
    threshold_30 <- 0.30 * h_max
    
    valid_pixels <- sum(!is.na(values(chm)))
    
    # Prevent division by zero
    if (valid_pixels > 0) {
      clumping_30h <- round(sum(values(chm) > threshold_30, na.rm = TRUE) / valid_pixels, 3)
    } else {
      clumping_30h <- 0
    }
    
    # C. MuSICA Setup
    pheno <- calc_phenology(list.year = c(2021, 2022),
                            nleafage = 1,
                            budburst_date = 115,
                            leaf_age_max_in = 0.56,
                            relative_age_firstmax = 0.10,
                            relative_age_lastmax = 0.75,
                            LAI_max_per_cohort = pai_val)
    
    # D. CALL MUSICA
    tryCatch({
      musica_out <- callmusica(
        musica.param = list(
          "setupctl" = list(
            "clumping_factor" = clumping_30h,
            "forcing_filename" = forcing_file,
            "history_filename" = history_file,
            "forcing_height" = h_max + 2 
          )
        ),
        leaf.param = list(
          "musica_veg1" = list(
            "phenology" = pheno,
            "allometry" = allometry,
            "leafphenologyctl" = list("lai_max_per_cohort" = pai_val),
            "leafallometryctl" = list("canopy_height_top" = h_max, "canopy_height_bottom" = 2),
            "leafmusicactl" = list("canopy_height_top" = h_max)
          )
        ),
        musica.cmd = musica.cmd,
        keep.tmp = FALSE,
        out.netcdf = TRUE
      )
      
      # Append metrics only if execution succeeds
      all_metrics_list[[length(all_metrics_list) + 1]] <- data.frame(
        id_plot = current_id,
        radius = radius,
        pai = pai_val,
        clumping = clumping_30h,
        h_max = round(h_max, 2),
        file_path = history_file
      )
      
    }, error = function(e) {
      message(sprintf(" -> Error during MuSICA execution for plot %s: %s", current_id, e$message))
    })
  }
}

# --- 3. FINAL SUMMARY ---

if (length(all_metrics_list) > 0) {
  df_summary <- bind_rows(all_metrics_list)
  summary_file <- file.path(base_out_dir, "simulation_summary.csv")
  
  # Append to existing summary or create new
  if (file.exists(summary_file)) {
    write.table(df_summary, summary_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  } else {
    write.csv(df_summary, summary_file, row.names = FALSE)
  }
  message("\nSimulations finished. Summary appended/saved in ", base_out_dir)
} else {
  message("\nNo new simulations were successfully completed.")
}











# Process HOBO data once outside the loop for maximum efficiency
df_hobo_clean <- read.csv(hobo_file) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  # Filter by sensor position, date range, and remove specific problematic IDs
  filter(
    position_sensor == "a", 
    as.Date(datetime) %in% date_seq,
    !(id_plot %in% ids_to_remove)
  ) %>%
  # Hourly aggregation to match MuSICA frequency
  group_by(id_plot, time = floor_date(datetime, "hour")) %>%
  summarise(Tair_obs = mean(t_hobo, na.rm = TRUE), .groups = "drop")

# --- 2. Comparison Loop over Radii and IDs ---
message("Starting NetCDF extraction and comparison...")

all_metrics <- map_df(radii, function(r) {
  
  # 1. List all .nc files for the current radius directory
  current_dir <- file.path(base_out_dir, paste0(r, "m"))
  files <- list.files(current_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning(sprintf("No NetCDF files found in: %s", current_dir))
    return(NULL)
  }
  
  message(sprintf(" -> Processing Radius: %sm (%d files)", r, length(files)))
  
  # 2. Process each file (each id_plot) within the radius
  map_df(files, function(f) {
    
    # Extract ID from filename (e.g., extracts '41_01' from 'musica_out_Blois_pt_41_01_radius_5.nc')
    pid <- str_extract(basename(f), "(?<=pt_).+(?=_radius)")
    
    # Skip if ID is in the exclusion list or extraction failed
    if (is.na(pid) || pid %in% ids_to_remove) return(NULL)
    
    # --- Extraction with get_variable ---
    nc <- try(nc_open(f))
    if (inherits(nc, "try-error")) return(NULL)
    
    raw_data <- try(get_variable(nc, "Tair_z"))
    nc_close(nc)
    
    if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
    
    # --- Processing MuSICA Data ---
    df_sim <- raw_data %>%
      filter(nair == 1) %>%                             # Bottom atmospheric layer
      mutate(Tair_sim = Tair_z - 273.15) %>%            # Convert Kelvin to Celsius
      filter(as.Date(time) %in% date_seq) %>%           # Filter to target date sequence
      mutate(time = time - hours(2)) %>%                # Apply -2h offset (e.g., UTC matching)
      mutate(time = floor_date(time, "hour")) %>%       # Hourly floor
      mutate(time = force_tz(time, "UTC")) %>%          # Force UTC timezone
      group_by(time) %>%
      summarise(Tair_sim = mean(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        radius = r,
        id_plot = pid
      )
    
    # --- Merge with HOBO and Metrics Calculation ---
    combined <- df_sim %>%
      inner_join(df_hobo_clean, by = c("time", "id_plot"))
    
    # Skip if no overlapping data is found
    if (nrow(combined) == 0) return(NULL)
    
    metrics <- combined %>%
      summarise(
        id_plot  = pid,
        radius   = r,
        n_points = n(),
        rmse     = round(sqrt(mean((Tair_sim - Tair_obs)^2, na.rm = TRUE)), 2),
        mae      = round(mean(abs(Tair_sim - Tair_obs), na.rm = TRUE), 2),
        bias     = round(mean(Tair_sim - Tair_obs, na.rm = TRUE), 2),
        r2       = round(cor(Tair_sim, Tair_obs, use = "complete.obs")^2, 2)
      )
    
    return(metrics)
  })
})

# --- 3. Save detailed metrics ---
# write.csv(all_metrics, file.path(base_out_dir, "detailed_plot_metrics.csv"), row.names = FALSE)
# message("\nDetailed metrics saved to 'detailed_plot_metrics.csv'")

# --- 4. Summary of metrics per radius ---
radius_comparison <- all_metrics %>%
  group_by(radius) %>%
  summarise(
    n_plots      = n(),                              # Number of plots successfully compared
    avg_rmse     = round(mean(rmse, na.rm = TRUE), 3),
    avg_mae      = round(mean(mae, na.rm = TRUE), 3),
    avg_bias     = round(mean(bias, na.rm = TRUE), 3),
    avg_r2       = round(mean(r2, na.rm = TRUE), 3),
    sd_rmse      = round(sd(rmse, na.rm = TRUE), 3)  # Variability between plots for this radius
  ) %>%
  arrange(avg_rmse) # Sort by best overall performance (lowest RMSE)

print(radius_comparison)
# write.csv(radius_comparison, file.path(base_out_dir, "summary_radius_comparison.csv"), row.names = FALSE)
# message("Summary metrics saved to 'summary_radius_comparison.csv'")















# ==============================================================================
# SCRIPT: Compare Delta Tmax (Canopy Cooling) between MuSICA and HOBO
# ==============================================================================

# --- 1. Configuration & Pre-processing ---

# Add the path to your forcing file (Macroclimate)
nc_force_file <- "./in_files/musica_in_Blois.nc" 

message("Loading and processing Macroclimate (Forcing) data...")
nc_force <- try(nc_open(nc_force_file))
if (inherits(nc_force, "try-error")) stop("Cannot open forcing file.")

# Extract forcing temperature (check if your variable is named 'Tair' or 'Tair_f')
raw_macro <- get_variable(nc_force, "Tair") 
nc_close(nc_force)

# Process Macroclimate to Daily Tmax
df_macro <- raw_macro %>%
  mutate(Tair_macro = Tair - 273.15) %>% 
  filter(as.Date(time) %in% date_seq) %>%
  mutate(time = time - hours(2)) %>%
  mutate(date = as.Date(force_tz(time, "UTC"))) %>%
  group_by(date) %>%
  summarise(Tmax_macro = max(Tair_macro, na.rm = TRUE), .groups = "drop")

message("Loading and processing HOBO data (Observed Microclimate Tmax)...")
df_hobo_clean <- read.csv(hobo_file) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(
    position_sensor == "a", 
    as.Date(datetime) %in% date_seq,
    !(id_plot %in% ids_to_remove)
  ) %>%
  # Aggregate to daily maximum instead of hourly mean
  mutate(date = as.Date(datetime)) %>%
  group_by(id_plot, date) %>%
  summarise(Tmax_obs = max(t_hobo, na.rm = TRUE), .groups = "drop")


# --- 2. Comparison Loop over Radii and IDs ---
message("Starting NetCDF extraction and Delta calculation...")

all_metrics <- map_df(radii, function(r) {
  
  current_dir <- file.path(base_out_dir, paste0(r, "m"))
  files <- list.files(current_dir, pattern = "\\.nc$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning(sprintf("No NetCDF files found in: %s", current_dir))
    return(NULL)
  }
  
  message(sprintf(" -> Processing Radius: %sm (%d files)", r, length(files)))
  
  map_df(files, function(f) {
    
    pid <- str_extract(basename(f), "(?<=pt_).+(?=_radius)")
    if (is.na(pid) || pid %in% ids_to_remove) return(NULL)
    
    # --- Extract MuSICA ---
    nc <- try(nc_open(f), silent = TRUE)
    if (inherits(nc, "try-error")) return(NULL)
    
    raw_data <- try(get_variable(nc, "Tair_z"), silent = TRUE)
    nc_close(nc)
    
    if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
    
    # --- Process MuSICA to Daily Tmax ---
    df_sim <- raw_data %>%
      filter(nair == 1) %>%                             
      mutate(Tair_sim = Tair_z - 273.15) %>%            
      filter(as.Date(time) %in% date_seq) %>%           
      mutate(time = time - hours(2)) %>%                
      mutate(date = as.Date(force_tz(time, "UTC"))) %>%
      group_by(date) %>%
      summarise(Tmax_sim = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        radius = r,
        id_plot = pid
      )
    
    # --- Merge all 3 sources and calculate Deltas ---
    combined <- df_sim %>%
      inner_join(df_macro, by = "date") %>%
      inner_join(df_hobo_clean, by = c("date", "id_plot")) %>%
      mutate(
        # Calculate offsets (Negative = canopy is cooler than macroclimate)
        Delta_sim = Tmax_sim - Tmax_macro,
        Delta_obs = Tmax_obs - Tmax_macro
      )
    
    if (nrow(combined) == 0) return(NULL)
    
    # --- Calculate Metrics on Deltas ---
    metrics <- combined %>%
      summarise(
        id_plot  = pid,
        radius   = r,
        n_days   = n(), # Number of valid days compared
        rmse     = round(sqrt(mean((Delta_sim - Delta_obs)^2, na.rm = TRUE)), 3),
        mae      = round(mean(abs(Delta_sim - Delta_obs), na.rm = TRUE), 3),
        bias     = round(mean(Delta_sim - Delta_obs, na.rm = TRUE), 3),
        r2       = round(cor(Delta_sim, Delta_obs, use = "complete.obs")^2, 3)
      )
    
    return(metrics)
  })
})

# --- 3. Summary of Delta metrics per radius ---
radius_comparison <- all_metrics %>%
  group_by(radius) %>%
  summarise(
    n_plots      = n(),                              
    avg_rmse     = round(mean(rmse, na.rm = TRUE), 3),
    avg_mae      = round(mean(mae, na.rm = TRUE), 3),
    avg_bias     = round(mean(bias, na.rm = TRUE), 3),
    avg_r2       = round(mean(r2, na.rm = TRUE), 3),
    sd_rmse      = round(sd(rmse, na.rm = TRUE), 3)  
  ) %>%
  arrange(avg_rmse) 

print(radius_comparison)