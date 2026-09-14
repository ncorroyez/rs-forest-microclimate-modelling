rm(list=ls(all=TRUE))
gc()

if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# --- Libraries ---
library(musica.tools) # Assuming this handles get_variable
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(terra)
library(sf)
library(lubridate)
library(purrr)
library(broom)
library(tibble)

source("functions.R") # Assuming get_equilibrium_stats is here

# --- 1. Preparation & Data Loading ---

# Dates
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")

# MuSICA parameters
target_radius <- 25

# --- Load Plot Coordinates & Metrics (LAD, VCI, etc) ---
# (Keeping your original logic here)
LAI_stack <- terra::rast("in_files/ladstack_classic_no_na.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)

df_coords <- json %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(coord_x_utm31n = X, coord_y_utm31n = Y)

lad_values <- terra::extract(LAI_stack, df_coords)
df_lad <- as.data.frame(lad_values)
df_lad <- cbind(df_lad, df_coords)
df_lad$ID <- json$id_plot

# Filter valid plots
df_lad <- df_lad %>%
  filter(rowSums(!is.na(dplyr::select(., -ID, -coord_x_utm31n, -coord_y_utm31n))) > 0)

n_field <- df_lad$ID
df_lad$ID = NULL

# CHS41 IGN
ign <- read.table("MetHor2021.txt", sep = ";", header = TRUE, 
                  stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1")

ign_Tair_all <- data.frame()
for (date_study in date_seq) {
  date_study <- as.Date(date_study)
  
  ign_Tair_filtered <- ign %>%
    filter(Date == format(date_study, "%d/%m/%Y")) %>%
    filter(Code.placette == "CHS 41") %>%
    dplyr::select(Code.placette, Date, Heure..TU., Température.instantanée...C.) %>%
    rename(Tair_ign = Température.instantanée...C.) %>%
    mutate(time = as.POSIXct(paste(Date, sprintf("%04d", Heure..TU.)), 
                             format = "%d/%m/%Y %H%M", tz = "UTC")) %>%
    dplyr::select(-Date, -Heure..TU., -Code.placette) %>%
    dplyr::select(time, Tair_ign)
  
  ign_Tair_all <- bind_rows(ign_Tair_all, ign_Tair_filtered)
}

# Metrics from CSV (VCI, PAI, Max)
# Re-loading metrics to ensure we have plot characteristics
metrics_csv <- read.csv("out_files/radius_test/metrics_results.csv", header = TRUE)
metrics_df <- metrics_csv %>%
  filter(radius == target_radius) %>%
  dplyr::select(id_plot, pai, vci, max, mean, lcv) # Select relevant metrics

# --- Load Observed Data (HOBO) ---
site <- "Blois"
csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
temperature_csv <- read.csv(csv_file)

temperature_df <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

temperature_all_days <- data.frame()
for (date_study in date_seq) {
  temp_filtered <- temperature_df %>%
    filter(as.Date(datetime) == as.Date(date_study),
           position_sensor == "a")
  temperature_all_days <- bind_rows(temperature_all_days, temp_filtered)
}
available_plots <- unique(temperature_all_days$id_plot)
# stop()
# ids_to_remove <- c("41_13", "41_14", "41_20", "41_34", 
#                    "41_41", "41_50", "41_51", "41_53",
#                    "41_17", "41_18", "41_19", "41_27",
#                    "41_30", "41_39", "41_47", "41_49",
#                    "41_55")
# available_plots <- available_plots[available_plots != ids_to_remove]

# Filter metrics to available plots
metrics_df <- metrics_df %>% filter(id_plot %in% available_plots)

# --- Load ERA5 (Macroclimate) ---
era5 <- nc_open("in_files/musica_in_Blois.nc")
era5_var <- get_variable(era5, "Tair")
nc_close(era5)

era5_Tair_all <- era5_var %>%
  filter(as.Date(time) %in% date_seq) %>%
  mutate(Tair = round(Tair - 273.15, 1)) %>%
  mutate(time = floor_date(time, unit = "hour"))

# --- Load SAFRAN (Macroclimate) ---
safran <- nc_open("in_files/musica_in_Safran_Blois_2021.nc")
safran_var <- get_variable(safran, "Tair")
nc_close(safran)

safran_Tair_all <- safran_var %>%
  filter(as.Date(time) %in% date_seq) %>%
  mutate(Tair = round(Tair - 273.15, 1)) %>%
  mutate(time = floor_date(time, unit = "hour"))

# --- 2. Load Simulated Data ---

# Define the three directories/scenarios
scenarios <- list(
  "ATBD_Prorata"                 = "out_files/TS/smooth/corr/atbd",
  "3_sites_Opt_Prorata"           = "out_files/TS/smooth/corr/atbd_optim_common",
  "3_sites_Opt_Brown_Prorata"     = "out_files/TS/smooth/corr/atbd_optim_common_brownmodif",
  "Blois_Only_Prorata"          = "out_files/TS/smooth/corr/optim_Blois",
  # "Optim_Blois_b_corr"        = "out_files/TS/corr/optim_Blois_brownmodif",
  
  # "ATBD_corrD"                = "out_files/TS/corr_d/atbd",
  # "3_sites_Opt_corrD"          = "out_files/TS/corr_d/atbd_optim_common",
  # "3_sites_Opt_Brown_corrD"    = "out_files/TS/corr_d/atbd_optim_common_brownmodif",
  # "Blois_Only_corrD"         = "out_files/TS/corr_d/optim_Blois",
  # "Optim_Blois_b_corrD"       = "out_files/TS/corr_d/optim_Blois_brownmodif",
  
  "ATBD_none"                  = "out_files/TS/smooth/gap/atbd",
  "3_sites_Opt_none"            = "out_files/TS/smooth/gap/atbd_optim_common",
  "3_sites_Opt_Brown_none"      = "out_files/TS/smooth/gap/atbd_optim_common_brownmodif",
  "Blois_Only_none"           = "out_files/TS/smooth/gap/optim_Blois",
  # "Optim_Blois_b_gap"         = "out_files/TS/gap/optim_Blois_brownmodif",
  
  "ATBD_RF"              = "out_files/TS/smooth/RF/LAI_ALS/atbd",
  "3_sites_Opt_RF"        = "out_files/TS/smooth/RF/LAI_ALS/atbd_optim_common",
  "3_sites_Opt_Brown_RF"  = "out_files/TS/smooth/RF/LAI_ALS/atbd_optim_common_brownmodif",
  "Blois_Only_RF"       = "out_files/TS/smooth/RF/LAI_ALS/optim_Blois",
  # "Optim_Blois_b_corr_RF"     = "out_files/TS/RF/LAI_ALS/optim_Blois_brownmodif",
  
  # "ATBD_corrD_RF"             = "out_files/TS/RF/LAI_ALS_dopt/atbd",
  # "ATBD_Optim_corrD_RF"       = "out_files/TS/RF/LAI_ALS_dopt/atbd_optim_common",
  # "ATBD_Optim_Brown_corrD_RF" = "out_files/TS/RF/LAI_ALS_dopt/atbd_optim_common_brownmodif",
  # "Optim_Blois_corrD_RF"      = "out_files/TS/RF/LAI_ALS_dopt/optim_Blois",
  # "Optim_Blois_b_corrD_RF"    = "out_files/TS/RF/LAI_ALS_dopt/optim_Blois_brownmodif",
  
  "Lidar_LAI_4"               = "out_files/TS/lidar/lai_4",
  "Lidar_Constant"            = "out_files/TS/lidar/lidar_constant"
)

df_Tair_all <- data.frame()
for (i in 1:nrow(df_lad)) {
  id_plot <- n_field[i]
  
  # Skip if plot not in observations
  if (!id_plot %in% available_plots) next
  
  if(i %% 5 == 0) print(paste("Processing plot:", id_plot, "(", i, "/", length(n_field), ")"))
  
  # Loop through the 3 scenarios for this plot
  for (sim_name in names(scenarios)) {
    
    base_path <- scenarios[[sim_name]]
    
    # Construct filename (Assuming standard naming convention inside the folders)
    # Check if 'remaining' meant subfolders or just the file. 
    # Usually: musica_out_Blois_pt_[ID]_radius_[R].nc
    nc_file <- file.path(base_path, 
                         paste0("musica_out_ERA5_Blois_2021_pt_", id_plot, ".nc"))
    
    if (file.exists(nc_file)) {
      out <- tryCatch(nc_open(nc_file), error = function(e) NULL)
      
      if (!is.null(out)) {
        # Extract Tair_z
        # Note: get_variable usually returns a dataframe with time, value, etc.
        # We need to filter for the specific dates to save memory/time
        
        # Extracting broadly then filtering
        raw_data <- get_variable(out, "Tair_z") 
        
        df_Tair <- raw_data %>%
          filter(as.Date(time) %in% date_seq, nair == 1) %>%
          mutate(
            Tair_z = Tair_z - 273.15,
            id_plot = id_plot,
            sim_type = sim_name # Add the scenario identifier
          ) %>%
          dplyr::select(time, Tair_z, id_plot, sim_type)
        
        df_Tair_all <- bind_rows(df_Tair_all, df_Tair)
        nc_close(out)
      }
    } else {
      warning(paste("File not found:", nc_file))
    }
  }
}

# Clean simulated data
df_Tair_all <- df_Tair_all %>%
  mutate(time = floor_date(time, unit = "hour")) %>%
  arrange(sim_type, id_plot, time)

# --- 3. Merging Data for Analysis ---

# 3.1 Aggregate simulated to hourly (if not already)
df_sim_hourly <- df_Tair_all %>%
  group_by(sim_type, id_plot, time = floor_date(time, "hour")) %>%
  summarise(Tair_sim = mean(Tair_z, na.rm = TRUE), .groups = "drop")

# 3.2 Aggregate observed to hourly
df_obs_hourly <- temperature_all_days %>%
  group_by(id_plot, time = floor_date(datetime, "hour")) %>%
  summarise(Tair_obs = mean(t_hobo, na.rm = TRUE), .groups = "drop")

# Add HOBO as scenario
df_obs_as_sim <- df_obs_hourly %>%
  rename(Tair_sim = Tair_obs) %>%
  mutate(sim_type = "HOBOs")
df_sim_hourly <- bind_rows(df_sim_hourly, df_obs_as_sim)

# 3.3 Merge Sim, Obs, and ERA5
# First merge Sim and Obs
combined_all <- merge(df_sim_hourly, df_obs_hourly, by = c("id_plot", "time"))

# Then merge with ERA5 & SAFRAN & IGN
combined_all <- merge(combined_all, era5_Tair_all, by = "time") %>% 
  rename(Tair_era5 = Tair)
combined_all <- merge(combined_all, safran_Tair_all, by = "time") %>% 
  rename(Tair_safran = Tair)
combined_all <- merge(combined_all, ign_Tair_all, by = "time")

# Finally merge with Metrics
combined_all <- merge(combined_all, metrics_df, by = "id_plot")

# Months and simulations
sim_levels_order <- c(
  "HOBOs",
  "Lidar_Constant",
  "Lidar_LAI_4",
  
  "ATBD_none",
  "ATBD_Prorata",
  "ATBD_RF",
  
  "3_sites_Opt_none",
  "3_sites_Opt_Prorata",
  "3_sites_Opt_RF",
  
  "3_sites_Opt_Brown_none",
  "3_sites_Opt_Brown_Prorata",
  "3_sites_Opt_Brown_RF",
  
  "Blois_Only_none",
  "Blois_Only_Prorata",
  "Blois_Only_RF"
)

combined_all <- combined_all %>%
  mutate(
    Month = factor(format(time, "%B"), 
                   levels = c("June", "July", "August", "September")),
    sim_type = factor(sim_type, levels = intersect(sim_levels_order, 
                                                   unique(sim_type)))
  ) %>%
  drop_na()



# target_dates <- as.Date(c(
#   # "2021-06-09",
#   "2021-06-14"
#   # "2021-07-19", "2021-07-22",
#   # "2021-08-26",
#   # "2021-09-02", "2021-09-07", "2021-09-22", "2021-09-30"
# ))
# combined_all <- combined_all %>%
#   filter(as.Date(time) %in% target_dates)

# --- 4. Statistical Analysis (Per Scenario) ---

# Calculate global stats per simulation type
global_stats <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R = cor(Tair_sim, Tair_obs),
    RMSE = sqrt(mean((Tair_sim - Tair_obs)^2)),
    Bias = mean(Tair_sim - Tair_obs),
    MAE = mean(abs(Tair_sim - Tair_obs)),
    Slope = coef(lm(Tair_sim ~ Tair_obs))[2],
    .groups = "drop"
  )

print("Global Statistics Sim per Model:")
print(global_stats)

global_stats_era5 <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R = cor(Tair_sim, Tair_era5),
    RMSE = sqrt(mean((Tair_sim - Tair_era5)^2)),
    Bias = mean(Tair_sim - Tair_era5),
    MAE = mean(abs(Tair_sim - Tair_era5)),
    Slope = coef(lm(Tair_sim ~ Tair_era5))[2],
    .groups = "drop"
  )

print("Global Statistics ERA5 per Model:")
print(global_stats)

global_stats_ign <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R = cor(Tair_sim, Tair_ign),
    RMSE = sqrt(mean((Tair_sim - Tair_ign)^2)),
    Bias = mean(Tair_sim - Tair_ign),
    MAE = mean(abs(Tair_sim - Tair_ign)),
    Slope = coef(lm(Tair_sim ~ Tair_ign))[2],
    .groups = "drop"
  )

print("Global Statistics IGN per Model:")
print(global_stats_ign)

global_stats_safran <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R = cor(Tair_sim, Tair_safran),
    RMSE = sqrt(mean((Tair_sim - Tair_safran)^2)),
    Bias = mean(Tair_sim - Tair_safran),
    MAE = mean(abs(Tair_sim - Tair_safran)),
    Slope = coef(lm(Tair_sim ~ Tair_safran))[2],
    .groups = "drop"
  )

print("Global Statistics SAFRAN per Model:")
print(global_stats_safran)

# --- 5. Slope Analysis (Gril et al. 2023 approach) ---
# Calculate the slope (T_micro ~ T_macro) for every plot AND every simulation type

slope_analysis <- combined_all %>%
  group_by(sim_type, id_plot) %>%
  nest() %>%
  mutate(
    # Model: Sim ~ ERA5 (Microclimate buffering relative to Macro)
    model = map(data, ~lm(Tair_sim ~ Tair_era5, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Tair_era5") %>%
  dplyr::select(sim_type, id_plot, estimate) %>%
  rename(slope = estimate) %>%
  mutate(log_slope = log(abs(slope))) %>%
  ungroup()

# Merge metrics back into slope analysis
slope_analysis <- merge(slope_analysis, metrics_df, by = "id_plot")

# --- 6. Plotting ---

# 
# Plot 1: Observed vs Simulated (Hexbin) faceted by Scenario
stats_labels <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R     = cor(Tair_sim, Tair_obs, use = "complete.obs"),
    RMSE  = sqrt(mean((Tair_sim - Tair_obs)^2, na.rm = TRUE)),
    Bias  = mean(Tair_sim - Tair_obs, na.rm = TRUE),
    Slope = coef(lm(Tair_sim ~ Tair_obs))[2],
    .groups = "drop"
  ) %>%
  mutate(
    # Création du texte formaté
    label_text = sprintf("R: %.2f\nRMSE: %.2f\nBias: %.2f\nSlope: %.2f", 
                         R, RMSE, Bias, Slope)
  )
p_scatter <- ggplot(combined_all, aes(x = Tair_obs, y = Tair_sim)) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "blue", linewidth = 0.5, se = FALSE) +
  facet_wrap(~sim_type, ncol = 3) +
  geom_text(data = stats_labels, aes(label = label_text),
            x = -Inf, y = Inf,
            hjust = -0.1,
            vjust = 1.1,
            size = 2,
            color = "black",
            fontface = "bold",
            inherit.aes = FALSE,
            check_overlap = TRUE) +
  
  labs(
    title = "Simulated vs Observed Temperature",
    x = "Observed Temperature (HOBO) [°C]",
    y = "Simulated Temperature (MuSICA) [°C]",
    fill = "Count"
  ) +
  theme_bw() +
  coord_fixed()

print(p_scatter)
ggsave("Plots/Silvilaser/1sim_vs_obs_comparison.png", 
       p_scatter, width = 8, height = 8, dpi = 300)

# ---------------
combined_with_slope <- combined_all %>%
  left_join(slope_analysis, by = c("sim_type", "id_plot")) %>%
  # Create the categorical variable for plotting
  mutate(
    buffering_cat = case_when(
      log_slope < 0 ~ "Buffering (Slope < 1)",
      log_slope > 0 ~ "Amplification (Slope > 1)",
      TRUE ~ "Neutral"
    )
  )

p_slope_split <- ggplot(combined_with_slope, aes(x = Tair_obs, y = Tair_sim)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "magma") +
  
  # Reference line (1:1)
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  
  # Trend line
  geom_smooth(method = "lm", color = "white", linewidth = 0.5, se = FALSE) +
  
  # FACET GRID: Rows = Sim Type, Cols = Buffering Regime
  facet_grid(sim_type ~ buffering_cat) +
  
  labs(
    title = "Simulated vs Observed Temperature (By Buffering Regime)",
    # subtitle = "Separated by Log Slope of Sim vs ERA5",
    x = "Observed Temperature (HOBO) [°C]",
    y = "Simulated Temperature (MuSICA) [°C]",
    fill = "Count"
  ) +
  theme_bw() +
  coord_fixed() +
  theme(legend.position = "right")

print(p_slope_split)

split_stats <- combined_with_slope %>%
  group_by(sim_type, buffering_cat) %>%
  summarise(
    R = cor(Tair_sim, Tair_obs),
    RMSE = sqrt(mean((Tair_sim - Tair_obs)^2)),
    Bias = mean(Tair_sim - Tair_obs),
    MAE = mean(abs(Tair_sim - Tair_obs)),
    Slope = coef(lm(Tair_sim ~ Tair_obs))[2],
    .groups = "drop"
  )

print("Global Statistics per Model:")
print(split_stats)

# 
# Plot 2: Relationship between VCI and Buffering (log_slope)
# Comparing the 3 models against each other regarding structural effect

p_vci_slope <- ggplot(slope_analysis, aes(x = vci, y = log_slope, color = sim_type, fill = sim_type)) +
  geom_point(alpha = 0.3) +
  # Add trend lines
  geom_smooth(method = "lm", formula = y ~ stats::poly(x, 2), se = F, alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  # scale_color_brewer(palette = "Dark2") +
  # scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Effect of Vertical Complexity (VCI) on Microclimate Buffering",
    subtitle = "Comparison of 3 MuSICA Configurations",
    x = "Vertical Complexity Index (VCI)",
    y = "log(slope) [Buffering Capacity]",
    color = "Simulation",
    fill = "Simulation"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

print(p_vci_slope)
# ggsave("Plots/Silvilaser/vci_vs_slope_comparison.png", p_vci_slope, width = 8, height = 6)


# 
# Plot 3: Boxplot of RMSE per plot to see stability
rmse_per_plot <- combined_all %>%
  group_by(sim_type, id_plot) %>%
  summarise(RMSE = sqrt(mean((Tair_sim - Tair_obs)^2)), .groups = "drop")

p_rmse <- ggplot(rmse_per_plot, aes(x = sim_type, y = RMSE, fill = sim_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  # scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of RMSE per Plot",
    x = "Simulation Configuration",
    y = "RMSE (°C)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(p_rmse)
# ggsave("Plots/Silvilaser/rmse_boxplot.png", p_rmse, width = 6, height = 6)

# --- 7. Create Composite Figure for Silvilaser ---

# Let's create a composite of the Scatter (Obs vs Sim) and the Slope Analysis
# We grab the scatter plot and the PAI/VCI/MAX vs Slope plots for the best model vs others

# Create PAI vs Slope plot
p_pai_slope <- ggplot(slope_analysis, aes(x = pai, y = log_slope, color = sim_type)) +
  geom_smooth(method = "lm", formula = y ~ stats::poly(x, 2), se = FALSE) +
  geom_point(alpha = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_color_brewer(palette = "Dark2") +
  labs(x = "PAI", y = "log(slope)") +
  theme_bw() +
  theme(legend.position = "none")

# Create Max Height vs Slope plot
p_max_slope <- ggplot(slope_analysis, aes(x = max, y = log_slope, color = sim_type)) +
  geom_smooth(method = "lm", formula = y ~ stats::poly(x, 2), se = FALSE) +
  geom_point(alpha = 0.2, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  # scale_color_brewer(palette = "Dark2") +
  labs(x = "Max Height (m)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

# Update VCI plot for grid
p_vci_clean <- p_vci_slope + 
  theme(legend.position = "bottom") + 
  labs(title = NULL, subtitle = NULL, y = NULL)

# Combine structural drivers plots
drivers_row <- grid.arrange(
  p_pai_slope, 
  p_max_slope, 
  p_vci_clean, 
  nrow = 1, 
  widths = c(1, 1, 1.5) # Give more space to VCI as it has the legend
)

# Final arrangement
final_grid <- grid.arrange(p_scatter, drivers_row, nrow = 2, heights = c(1, 0.8))

# ggsave("Plots/Silvilaser/final_comparison_figure.png", final_grid, width = 12, height = 10, dpi = 300)

# --- 4. Statistical Analysis per Month ---

# Global stats (All months combined)
global_stats <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R = round(cor(Tair_sim, Tair_obs), 2),
    RMSE = round(sqrt(mean((Tair_sim - Tair_obs)^2)), 2),
    Bias = round(mean(Tair_sim - Tair_obs), 2),
    .groups = "drop"
  )
print("Global Statistics per Model:")
print(global_stats)

# Monthly stats (New table)
monthly_stats <- combined_all %>%
  group_by(sim_type, Month) %>%
  summarise(
    R = round(cor(Tair_sim, Tair_obs), 2),
    RMSE = round(sqrt(mean((Tair_sim - Tair_obs)^2)), 2),
    Bias = round(mean(Tair_sim - Tair_obs), 2),
    .groups = "drop"
  )

print("Monthly Statistics per Model:")
print(monthly_stats)

# --- 5. Slope Analysis per Month ---
slope_analysis_monthly <- combined_all %>%
  group_by(sim_type, id_plot, Month) %>% # Added Month to grouping
  nest() %>%
  mutate(
    # Model: Sim ~ ERA5
    model = map(data, ~lm(Tair_sim ~ Tair_era5, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Tair_era5") %>%
  dplyr::select(sim_type, id_plot, Month, estimate) %>%
  rename(slope = estimate) %>%
  mutate(log_slope = log(abs(slope))) %>%
  ungroup()

my_colors <- c(
  # --- GROUPE DE RÉFÉRENCE (Noir / Gris Foncé) ---
  "HOBOs"                      = "black",
  "Lidar_Constant"             = "#333333",
  "Lidar_LAI_4"                = "#333333", 
  
  # --- FAMILLE 1 : ATBD (Tons Rouges/Oranges) ---
  "ATBD_none"                  = "#FCAE91", # Clair (Saumon)
  "ATBD_Prorata"               = "#FB6A4A", # Moyen
  "ATBD_RF"                    = "#CB181D", # Foncé (Rouge vif)
  
  # --- FAMILLE 2 : 3_sites_Opt (Tons Bleus) ---
  "3_sites_Opt_none"           = "#BDD7E7", # Clair
  "3_sites_Opt_Prorata"        = "#6BAED6", # Moyen
  "3_sites_Opt_RF"             = "#2171B5", # Foncé
  
  # --- FAMILLE 3 : 3_sites_Opt_Brown (Tons Violets) ---
  "3_sites_Opt_Brown_none"     = "#CBC9E2", # Clair
  "3_sites_Opt_Brown_Prorata"  = "#9E9AC8", # Moyen
  "3_sites_Opt_Brown_RF"       = "#6A51A3", # Foncé
  
  # --- FAMILLE 4 : Blois_Only (Tons Verts) ---
  "Blois_Only_none"            = "#BAE4B3", # Clair
  "Blois_Only_Prorata"         = "#74C476", # Moyen
  "Blois_Only_RF"              = "#238B45"  # Foncé
)

# Plot Drivers vs Buffering (Split by Month) 
# This helps verify if PAI impact on buffering is consistent or fades in September
p_drivers_monthly <- ggplot(slope_analysis_monthly, aes(x = sim_type, y = log_slope, fill = sim_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Evolution of Buffering Capacity (log slope) by Month",
    y = "log(slope) [ <0 = Buffering ]",
    x = NULL
  ) +
  facet_wrap(~Month, ncol = 2) +
  theme_bw()

print(p_drivers_monthly)

p_drivers_monthly2 <- ggplot(slope_analysis_monthly, aes(x = Month, y = log_slope, fill = sim_type)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = my_colors) +
  labs(
    title = "Evolution of Buffering Capacity (log slope) by sim_type",
    y = "log(slope) [ <0 = Buffering ]",
    x = NULL
  ) +
  facet_wrap(~sim_type, ncol = 3) +
  theme_bw()

print(p_drivers_monthly2)
# stop()
# --- 6. Plotting (Updated for Monthly Split) ---

# Plot 1: Observed vs Simulated (Hexbin) faceted by Scenario AND Month
p_scatter_monthly <- ggplot(combined_all, aes(x = Tair_obs, y = Tair_sim)) +
  geom_hex(bins = 60) + # Reduced bins slightly for smaller panels
  scale_fill_viridis_c(option = "magma") +
  
  # 1:1 Reference Line
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  
  # Trend line
  geom_smooth(method = "lm", color = "white", linewidth = 0.5, se = FALSE) +
  
  # THE CHANGE IS HERE: Grid with Sim on rows, Month on columns
  facet_grid(sim_type ~ Month) +
  
  labs(
    title = "Simulated vs Observed Temperature (Monthly Evolution)",
    subtitle = "Comparing performance across Summer 2021",
    x = "Observed Temperature (HOBO) [°C]",
    y = "Simulated Temperature (MuSICA) [°C]",
    fill = "Count"
  ) +
  theme_bw() +
  coord_fixed() + 
  theme(strip.text = element_text(size = 10, face = "bold"))

print(p_scatter_monthly)

# --- 7. Diurnal Cycle Analysis ---

# Extract Hour
combined_all <- combined_all %>%
  mutate(Hour = hour(time))

# Calculate stats per Hour and Sim Type
diurnal_stats <- combined_all %>%
  group_by(sim_type, Hour) %>%
  summarise(
    Bias = mean(Tair_sim - Tair_obs),
    RMSE = sqrt(mean((Tair_sim - Tair_obs)^2)),
    .groups = "drop"
  )

diurnal_10 <- diurnal_stats %>%
  filter(Hour == 10)
diurnal_14 <- diurnal_stats %>%
  filter(Hour == 14)

# Plot Bias per Hour
p_diurnal <- ggplot(diurnal_stats, aes(x = Hour, y = Bias, color = sim_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_fill_manual(values = my_colors) +
  scale_x_continuous(breaks = seq(0, 23, 2)) +
  labs(
    title = "Diurnal Cycle of Model Bias",
    x = "Hour of Day (UTC)",
    y = "Bias (Sim - Obs) [°C]"
  ) +
  # xlim(c(8, 12)) +
  # ylim(c(3, 4.2)) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_diurnal)
# 

# --- 8. Bias vs Structure Analysis ---

# Calculate mean Bias per plot
bias_per_plot <- combined_all %>%
  group_by(sim_type, id_plot) %>%
  summarise(
    Mean_Bias = mean(Tair_sim - Tair_obs),
    PAI = first(pai), # PAI is constant per plot
    VCI = first(vci),
    .groups = "drop"
  )

# Plot Bias vs PAI
p_bias_pai <- ggplot(bias_per_plot, aes(x = PAI, y = Mean_Bias, color = sim_type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + # Add trend line
  scale_color_manual(values = my_colors) +
  labs(
    title = "Systematic Bias vs Canopy Density (PAI)",
    subtitle = "Ideally, trend lines should be flat (no structural bias)",
    x = "Plant Area Index (PAI)",
    y = "Mean Bias (°C)"
  ) +
  theme_bw()

print(p_bias_pai)

# --- 9. Extreme Events Analysis ---

# Define "Hot Days" as days where Obs Temp > 90th percentile
threshold_hot <- quantile(combined_all$Tair_obs, 0.90)

hot_stats <- combined_all %>%
  filter(Tair_obs > threshold_hot) %>%
  group_by(sim_type) %>%
  summarise(
    RMSE_Hot = sqrt(mean((Tair_sim - Tair_obs)^2)),
    Bias_Hot = mean(Tair_sim - Tair_obs),
    count = n()
  )

print("Performance during Extreme Heat (> 90th percentile):")
print(hot_stats)
















# --- 1. Classify Regimes (Buffering vs Amplification) ---
# We use the previously calculated slope analysis
slope_classified <- slope_analysis_monthly %>%
  mutate(
    Regime = case_when(
      log_slope < 0 ~ "Buffering (Slope < 1)",
      log_slope > 0 ~ "Amplification (Slope > 1)",
      TRUE ~ "Neutral"
    )
  ) %>%
  dplyr::select(sim_type, id_plot, Month, Regime)

# --- 2. Merge with Hourly Data ---
# Join regime info to the main dataset
combined_regime <- combined_all %>%
  left_join(slope_classified, by = c("sim_type", "id_plot", "Month")) %>%
  filter(!is.na(Regime)) %>% 
  # EXCLUDE HOBO FROM SCENARIOS
  filter(sim_type != "HOBO_Ref")

# --- 3. Calculate Statistics for Labels ---
# We compute stats per group to display them on the plot
stats_labels <- combined_regime %>%
  group_by(sim_type, Regime) %>%
  summarise(
    R     = cor(Tair_sim, Tair_obs, use = "complete.obs"),
    RMSE  = sqrt(mean((Tair_sim - Tair_obs)^2, na.rm = TRUE)),
    Bias  = mean(Tair_sim - Tair_obs, na.rm = TRUE),
    Slope = coef(lm(Tair_sim ~ Tair_obs))[2], # Regression slope Sim ~ Obs
    .groups = "drop"
  ) %>%
  mutate(
    # Create a single string for the label
    label_text = sprintf("R: %.2f\nRMSE: %.2f\nBias: %.2f\nSlope: %.2f", 
                         R, RMSE, Bias, Slope)
  )

# --- 4. Generate the Plot ---
p_scatter_metrics <- ggplot(combined_regime, aes(x = Tair_obs, y = Tair_sim)) +
  
  # A. Hexbin density
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "magma", name = "Count") +
  
  # B. Reference Line (1:1) - Red dashed
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.6) +
  
  # C. Trend Line (Model behavior) - Cyan
  geom_smooth(method = "lm", color = "cyan", linewidth = 0.5, se = FALSE) +
  
  # D. Faceting: Rows = Simulation, Columns = Regime
  facet_grid(sim_type ~ Regime) +
  
  # E. Add Metrics Text
  # We use the pre-calculated stats_labels dataframe
  geom_text(data = stats_labels, aes(label = label_text),
            x = -Inf, y = Inf, # Position at top-left
            hjust = -0.1,      # Slight padding from left edge
            vjust = 1.1,       # Slight padding from top edge
            size = 3, 
            color = "black",   # White text stands out on dark hexbins; use black if background is light
            fontface = "bold") +
  
  # F. Aesthetics
  labs(
    title = "MuSICA Performance by Microclimatic Regime",
    subtitle = "Comparison of Scenarios under Buffering (left) vs. Amplification (right) conditions",
    x = "Observed Temperature (HOBO) [°C]",
    y = "Simulated Temperature (MuSICA) [°C]"
  ) +
  theme_bw() +
  # coord_fixed() + # Ensures square aspect ratio (important for 1:1 comparison)
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"), # Horizontal scenario labels
    legend.position = "bottom"
  )

# Print
print(p_scatter_metrics)

# Save (Recommended for large grids)
ggsave("Plots/TS/scatter_regime_metrics.png", p_scatter_metrics, 
       width = 10, height = 10, dpi = 100)







# ERA5
# 1. Calculate the Slope (Buffering Capacity) for each scenario
# This quantifies the "decoupling" from the macroclimate.
stats_buffering <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R     = cor(Tair_sim, Tair_era5, use = "complete.obs"),
    RMSE  = sqrt(mean((Tair_sim - Tair_era5)^2, na.rm = TRUE)),
    Bias  = mean(Tair_sim - Tair_era5, na.rm = TRUE),
    Slope = coef(lm(Tair_sim ~ Tair_era5))[2],
    .groups = "drop"
  ) %>%
  mutate(
    label_text = sprintf("R: %.2f\nRMSE: %.2f\nBias: %.2f\nSlope: %.2f", 
                         R, RMSE, Bias, Slope)
  )

# 2. Create the Scatterplot
p_micro_macro <- ggplot(combined_all, aes(x = Tair_era5, y = Tair_sim)) +
  
  # A. Hexbin for density (handling millions of points)
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "cividis", name = "Count") +
  
  # B. Reference Line (1:1) -> Represents the Macroclimate (No buffering)
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  
  # C. Trend Line (Actual Buffering)
  geom_smooth(method = "lm", color = "cyan", linewidth = 0.6, se = FALSE) +
  
  # D. Facet by Scenario (including HOBO_Ref if previously added)
  facet_wrap(~sim_type, ncol = 5) +
  
  # E. Add Slope Statistics text
  geom_text(data = stats_buffering, aes(label = label_text),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            size = 3.5, fontface = "bold", color = "black") +
  
  # F. Aesthetics
  labs(
    title = "Microclimate vs Macroclimate (Buffering Analysis)",
    x = "Macroclimate Temperature (ERA5) [°C]",
    y = "Microclimate Temperature (Sim or Obs) [°C]"
  ) +
  theme_bw() +
  coord_fixed() # Square aspect ratio is crucial for comparing slopes

print(p_micro_macro)

# IGN
# 1. Calculate the Slope (Buffering Capacity) for each scenario
# This quantifies the "decoupling" from the macroclimate.
stats_buffering <- combined_all %>%
  group_by(sim_type) %>%
  summarise(
    R     = cor(Tair_sim, Tair_ign, use = "complete.obs"),
    RMSE  = sqrt(mean((Tair_sim - Tair_ign)^2, na.rm = TRUE)),
    Bias  = mean(Tair_sim - Tair_ign, na.rm = TRUE),
    Slope = coef(lm(Tair_sim ~ Tair_ign))[2],
    .groups = "drop"
  ) %>%
  mutate(
    label_text = sprintf("R: %.2f\nRMSE: %.2f\nBias: %.2f\nSlope: %.2f", 
                         R, RMSE, Bias, Slope)
  )

# 2. Create the Scatterplot
p_micro_macro <- ggplot(combined_all, aes(x = Tair_ign, y = Tair_sim)) +
  
  # A. Hexbin for density (handling millions of points)
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "cividis", name = "Count") +
  
  # B. Reference Line (1:1) -> Represents the Macroclimate (No buffering)
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
  
  # C. Trend Line (Actual Buffering)
  geom_smooth(method = "lm", color = "cyan", linewidth = 0.6, se = FALSE) +
  
  # D. Facet by Scenario (including HOBO_Ref if previously added)
  facet_wrap(~sim_type, ncol = 5) +
  
  # E. Add Slope Statistics text
  geom_text(data = stats_buffering, aes(label = label_text),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            size = 3.5, fontface = "bold", color = "black") +
  
  # F. Aesthetics
  labs(
    title = "Microclimate vs Macroclimate (Buffering Analysis)",
    x = "Macroclimate Temperature (ERA5) [°C]",
    y = "Microclimate Temperature (Sim or Obs) [°C]"
  ) +
  theme_bw() +
  coord_fixed() # Square aspect ratio is crucial for comparing slopes

print(p_micro_macro)







slope_analysis_all <- combined_all %>%
  # We group by Simulation, Plot, and Month to capture seasonal changes in regime
  group_by(sim_type, id_plot, Month) %>%
  nest() %>%
  mutate(
    # Model: Micro (Tair_sim) ~ Macro (Tair_era5)
    model = map(data, ~lm(Tair_sim ~ Tair_era5, data = .)),
    tidied = map(model, tidy)
  ) %>%
  unnest(tidied) %>%
  filter(term == "Tair_era5") %>%
  dplyr::select(sim_type, id_plot, Month, estimate) %>%
  rename(slope = estimate) %>%
  mutate(
    log_slope = log(abs(slope)),
    # Create the classification Label
    Regime = case_when(
      log_slope < 0 ~ "Buffering (Slope < 1)",
      log_slope > 0 ~ "Amplification (Slope > 1)",
      TRUE ~ "Neutral"
    )
  ) %>%
  ungroup()

# --- 2. Merge Regime back to Hourly Data ---
combined_regime_macro <- combined_all %>%
  inner_join(slope_classified %>% 
               dplyr::select(sim_type, id_plot, Month, Regime), 
             by = c("sim_type", "id_plot", "Month"))

# --- 3. Calculate Stats for Display (Slope & R2 per panel) ---
stats_macro_split <- combined_regime_macro %>%
  group_by(sim_type, Regime) %>%
  summarise(
    Slope_Global = coef(lm(Tair_sim ~ Tair_ign))[2],
    R2 = cor(Tair_sim, Tair_ign)^2,
    .groups = "drop"
  ) %>%
  mutate(
    label_text = sprintf("Slope: %.2f\nR²: %.2f", Slope_Global, R2)
  )

# --- 4. Plot Micro vs Macro (Split by Regime) ---
p_micro_macro_split <- ggplot(combined_regime_macro, aes(x = Tair_ign, y = Tair_sim)) +
  
  # A. Hexbin Density
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "cividis", name = "Count") +
  
  # B. Reference Line 1:1 (Macroclimate)
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.6) +
  
  # C. Trend Line (Microclimate behavior)
  geom_smooth(method = "lm", color = "cyan", linewidth = 0.6, se = FALSE) +
  
  # D. Facet Grid: Rows = Scenario, Cols = Regime
  facet_grid(sim_type ~ Regime) +
  
  # E. Statistics Labels
  geom_text(data = stats_macro_split, aes(label = label_text),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            size = 3, fontface = "bold", color = "black") +
  
  # F. Aesthetics
  labs(
    title = "Microclimate vs. Macroclimate: Buffering vs. Amplification",
    subtitle = "Left: Plots/Months where forest acts as a buffer. Right: Where it amplifies extremes.\nRed Line = Open Field (ERA5). Cyan Line = Forest Trend.",
    x = "Macroclimate Temperature (ERA5) [°C]",
    y = "Microclimate Temperature (Sim or Obs) [°C]"
  ) +
  theme_bw() +
  # coord_fixed() +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold", size = 8),
    legend.position = "bottom"
  )

print(p_micro_macro_split)

ggsave("Plots/TS/p_micro_macro.png", p_micro_macro_split, 
       width = 10, height = 10, dpi = 100)
