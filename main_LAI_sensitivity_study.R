rm(list=ls(all=TRUE))
gc()

# --- 1. Setup & Libraries ---
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(musica.tools) # Required for get_variable
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(broom)
library(viridis) # For color gradients

# --- 2. Configuration & Constants ---
# Define the date range for the study
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")

LAI_stack <- terra::rast("in_files/ladstack_classic_no_na.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)
metrics <- read.csv("in_files/metrics_results_25.csv") %>%
  dplyr::select(!pai) 

df_coords <- json %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(coord_x_utm31n = X, coord_y_utm31n = Y)

# lad_values <- terra::extract(LAI_stack, df_coords)
# df_lad <- as.data.frame(lad_values)
# df_lad <- cbind(df_lad, df_coords)
# df_lad$ID <- json$id_plot

df_lad <- read.csv("in_files/allometry_profiles.csv")
n_field <- df_lad$id_plot

# Filter valid plots
# df_lad <- df_lad %>%
#   filter(rowSums(!is.na(dplyr::select(., -ID, -coord_x_utm31n, -coord_y_utm31n))) > 0)
# 
# n_field <- df_lad$ID
# df_lad$ID = NULL

# Define Sensitivity Scenarios (LAI 1 to 8 + Constant)
scenarios <- list()
for (i in 1:8) {
  scenarios[[paste0("LAI_", i)]] <- file.path("out_files/TS/lidar", paste0("lai_", i))
}
scenarios[["Lidar_Constant"]] <- "out_files/TS/lidar/lidar_constant"

# Create a specific factor order for plotting (Low LAI -> High LAI -> Constant)
sim_order <- c(paste0("LAI_", 1:8), "Lidar_Constant", "HOBO_Ref")

# --- 3. Load Observational Data (HOBOs) ---
print("Loading HOBO data...")
site <- "Blois"
csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
temperature_csv <- read.csv(csv_file)

# Format dates and filter
temperature_df <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a") %>%
  filter(as.Date(datetime) %in% date_seq)

# Aggregate HOBOs to Hourly
df_obs_hourly <- temperature_df %>%
  group_by(id_plot, time = floor_date(datetime, "hour")) %>%
  summarise(Tair_obs = mean(t_hobo, na.rm = TRUE), .groups = "drop")

# List of valid plots based on observations
available_plots <- unique(df_obs_hourly$id_plot)

# --- 4. Load Macroclimate Data (IGN) ---
# Used as the reference for calculating slopes (Micro vs Macro)
print("Loading IGN Macroclimate data...")
ign <- read.table("MetHor2021.txt", sep = ";", header = TRUE, 
                  stringsAsFactors = FALSE, fileEncoding = "ISO-8859-1")

ign_Tair_all <- data.frame()
for (date_study in date_seq) {
  # Format date for IGN file matching
  ign_filtered <- ign %>%
    filter(Date == format(as.Date(date_study), "%d/%m/%Y")) %>%
    filter(Code.placette == "CHS 41") %>%
    rename(Tair_ign = Température.instantanée...C.) %>%
    mutate(time = as.POSIXct(paste(Date, sprintf("%04d", Heure..TU.)), 
                             format = "%d/%m/%Y %H%M", tz = "UTC")) %>%
    dplyr::select(time, Tair_ign)
  
  ign_Tair_all <- bind_rows(ign_Tair_all, ign_filtered)
}

# --- 5. Load Simulated Data (MuSICA) ---
print("Loading Simulation data (this may take time)...")
df_sim_raw <- data.frame()

# Iterate over observed plots
for (id_plot in available_plots) {
  
  # Progress indicator
  if(which(available_plots == id_plot) %% 5 == 0) {
    print(paste("Processing plot:", id_plot))
  }
  
  # Iterate over scenarios
  for (sim_name in names(scenarios)) {
    
    # Construct file path
    base_path <- scenarios[[sim_name]]
    nc_file <- file.path(base_path, paste0("musica_out_ERA5_Blois_2021_pt_", id_plot, ".nc"))
    
    if (file.exists(nc_file)) {
      # Safely open NetCDF
      out <- tryCatch(nc_open(nc_file), error = function(e) NULL)
      
      if (!is.null(out)) {
        # Extract Temperature Variable
        # Note: Depending on musica.tools version, check if 'get_variable' works or use 'ncvar_get'
        raw_data <- tryCatch(get_variable(out, "Tair_z"), error = function(e) NULL)
        
        if(!is.null(raw_data)) {
          # Filter and Process
          df_chunk <- raw_data %>%
            filter(as.Date(time) %in% date_seq, nair == 1) %>% # nair=1 usually lowest layer
            mutate(
              Tair_z = Tair_z - 273.15, # Kelvin to Celsius
              id_plot = id_plot,
              sim_type = sim_name
            ) %>%
            dplyr::select(time, Tair_z, id_plot, sim_type)
          
          df_sim_raw <- bind_rows(df_sim_raw, df_chunk)
        }
        nc_close(out)
      }
    }
  }
}

# Aggregate Simulations to Hourly
df_sim_hourly <- df_sim_raw %>%
  mutate(time = floor_date(time, unit = "hour")) %>%
  group_by(sim_type, id_plot, time) %>%
  summarise(Tair_sim = mean(Tair_z, na.rm = TRUE), .groups = "drop")

# --- 6. Data Merging ---
print("Merging datasets...")

# Merge Sim and Obs
combined_sensi <- merge(df_sim_hourly, df_obs_hourly, by = c("id_plot", "time"))

# Merge with Macro (IGN)
combined_sensi <- merge(combined_sensi, ign_Tair_all, by = "time")

# Set Factor Levels for clean plotting order
combined_sensi$sim_type <- factor(combined_sensi$sim_type, levels = sim_order)

# ------------------------------------------------------------------------------
# --- ANALYSIS 1: Slope Scatterplot (Sim vs Obs) ---
# ------------------------------------------------------------------------------
print("Generating Slope Analysis...")

# 1. Calculate Observed Slope (Reference X-axis)
# Linear Model: HOBO ~ IGN
obs_slopes <- combined_sensi %>%
  dplyr::select(id_plot, time, Tair_obs, Tair_ign) %>%
  distinct() %>%
  group_by(id_plot) %>%
  nest() %>%
  mutate(model = map(data, ~lm(Tair_obs ~ Tair_ign, data = .)),
         slope_obs = map_dbl(model, ~coef(.)[2])) %>%
  dplyr::select(id_plot, slope_obs)

# 2. Calculate Simulated Slope (Y-axis) per Scenario
# Linear Model: Sim ~ IGN
sim_slopes <- combined_sensi %>%
  group_by(sim_type, id_plot) %>%
  nest() %>%
  mutate(model = map(data, ~lm(Tair_sim ~ Tair_ign, data = .)),
         slope_sim = map_dbl(model, ~coef(.)[2])) %>%
  dplyr::select(sim_type, id_plot, slope_sim)

# 3. Merge
slope_data <- inner_join(sim_slopes, obs_slopes, by = "id_plot")

# 4. Plot
p_slope <- ggplot(slope_data, aes(x = slope_obs, y = slope_sim, color = sim_type)) +
  # 1:1 Reference Line
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  
  # Points
  geom_point(alpha = 0.6, size = 2) +
  
  # Trend lines per scenario
  geom_smooth(method = "lm", se = FALSE, size = 0.8, alpha = 0.5) +
  
  # Colors: Plasma is good for ordered sequences (LAI 1->8)
  scale_color_viridis_d(option = "plasma", name = "Scenario") +
  
  labs(
    title = "Buffering Capacity Sensitivity Analysis",
    subtitle = "Relationship between Observed and Simulated Coupling with Macroclimate",
    x = "Observed Slope (HOBO vs IGN)",
    y = "Simulated Slope (MuSICA vs IGN)"
  ) +
  theme_bw() +
  coord_fixed(ratio = 1, xlim = c(0.4, 1.2), ylim = c(0.4, 1.2))

print(p_slope)
# ggsave("Plots/Sensitivity/1_Slope_Comparison_LAI.png", p_slope, width = 8, height = 7)

# ------------------------------------------------------------------------------
# --- ANALYSIS 2: Hourly Time Series (Mean + SD Ribbon) ---
# ------------------------------------------------------------------------------
print("Generating Hourly Time Series...")

# Calculate Stats for Simulations
ts_sim_stats <- combined_sensi %>%
  group_by(time, sim_type) %>%
  summarise(
    mean_val = mean(Tair_sim, na.rm = TRUE),
    sd_val = sd(Tair_sim, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Stats for HOBO (Reference)
ts_obs_stats <- combined_sensi %>%
  group_by(time) %>%
  summarise(
    mean_val = mean(Tair_obs, na.rm = TRUE),
    sd_val = sd(Tair_obs, na.rm = TRUE),
    sim_type = "HOBO_Ref",
    .groups = "drop"
  )

# Combine
ts_all <- bind_rows(ts_sim_stats, ts_obs_stats)
ts_all$sim_type <- factor(ts_all$sim_type, levels = sim_order)

# Define colors manually to ensure HOBO is Black
# Get n colors from plasma for the scenarios
n_scenarios <- length(levels(ts_all$sim_type)) - 1
pal <- viridis::plasma(n_scenarios)
color_map <- c(pal, "black") # Append black for HOBO
names(color_map) <- levels(ts_all$sim_type)

# Plot
p_ts_hourly <- ggplot(ts_all, aes(x = time, y = mean_val, color = sim_type, fill = sim_type)) +
  # Ribbon for SD (Spatial Variability)
  geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val), 
              alpha = 0.1, color = NA) +
  
  # Mean Line
  geom_line(size = 0.6) +
  
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  
  labs(
    title = "Hourly Temperature Evolution (Mean ± SD)",
    subtitle = "Comparison of LAI Scenarios vs HOBO Observations (Black)",
    x = NULL,
    y = "Temperature (°C)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_ts_hourly)
# ggsave("Plots/Sensitivity/2_TimeSeries_Hourly.png", p_ts_hourly, width = 12, height = 6)

# ------------------------------------------------------------------------------
# --- ANALYSIS 3: Daily Metrics (Mean & Max) ---
# ------------------------------------------------------------------------------
print("Generating Daily Metrics Analysis...")

# Calculate Daily Stats (Sim)
daily_sim <- combined_sensi %>%
  mutate(date = as.Date(time)) %>%
  group_by(date, sim_type) %>%
  summarise(
    Daily_Mean = mean(Tair_sim, na.rm = TRUE),
    Daily_Max = max(Tair_sim, na.rm = TRUE),
    .groups = "drop"
  )

# Calculate Daily Stats (Obs)
daily_obs <- combined_sensi %>%
  mutate(date = as.Date(time)) %>%
  group_by(date) %>%
  summarise(
    Daily_Mean = mean(Tair_obs, na.rm = TRUE),
    Daily_Max = max(Tair_obs, na.rm = TRUE),
    sim_type = "HOBO_Ref",
    .groups = "drop"
  )

# Combine & Reshape for Faceting
daily_all <- bind_rows(daily_sim, daily_obs) %>%
  pivot_longer(cols = c(Daily_Mean, Daily_Max), 
               names_to = "Metric", 
               values_to = "Temperature")

daily_all$sim_type <- factor(daily_all$sim_type, levels = sim_order)

# Plot
p_daily <- ggplot(daily_all, aes(x = date, y = Temperature, color = sim_type)) +
  geom_line(size = 0.8, alpha = 0.8) +
  
  scale_color_manual(values = color_map) +
  
  facet_wrap(~Metric, ncol = 1, scales = "free_y") +
  
  labs(
    title = "Daily Aggregated Metrics",
    subtitle = "Sensitivity of Daily Mean and Maximum Temperature to LAI",
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_daily)
# ggsave("Plots/Sensitivity/3_TimeSeries_Daily.png", p_daily, width = 10, height = 8)

print("Done.")

# ------------------------------------------------------------------------------
# --- ANALYSIS 4: LAD Profile Analysis ---
# ------------------------------------------------------------------------------
print("Analyzing LAD Profiles...")

# --- A. Data Preparation (Reshape LAD) ---
# Assuming columns in df_lad (except ID/coords) are height layers ordered from ground up
# We need to reshape from Wide to Long format

# Identify LAD columns (exclude ID and coordinates if present)
lad_cols <- names(df_lad)[!names(df_lad) %in% c("ID",
                                                "coord_x_utm31n",
                                                "coord_y_utm31n", 
                                                "id_plot",
                                                "radius")]

# Add ID back if it was removed in your snippet (using n_field)
df_lad_clean <- df_lad
df_lad_clean$id_plot <- n_field

# Pivot to Long format for plotting/analysis
lad_long <- df_lad_clean %>%
  dplyr::select(id_plot, all_of(lad_cols)) %>%
  pivot_longer(cols = -id_plot, names_to = "layer_name", values_to = "LAD") %>%
  group_by(id_plot) %>%
  mutate(
    # Assuming columns are ordered by height. 
    # If layers are 1m thick, layer_index ~ height in meters.
    height = row_number() + 2
  ) %>%
  ungroup()

# --- B. Calculate Profile Metrics ---
# We characterize each profile by a few key numbers to find trends
profile_metrics <- lad_long %>%
  group_by(id_plot) %>%
  summarise(
    Total_PAI = sum(LAD, na.rm = TRUE),          # Total density
    Max_LAD = max(LAD, na.rm = TRUE),            # Peak density
    Height_Max_LAD = height[which.max(LAD)],     # Height of peak density
    # Center of Gravity (weighted mean height)
    Center_Gravity = sum(height * LAD, na.rm = TRUE) / sum(LAD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(metrics, by = "id_plot")

# --- C. Visualize Profiles ---
# Plot 1: All vertical profiles
p_profiles <- ggplot(lad_long, aes(x = LAD, y = height, group = id_plot, color = id_plot)) +
  geom_path(alpha = 0.5) +
  labs(
    title = "Vertical Leaf Area Density (LAD) Profiles",
    subtitle = "Vegetation structure for all study plots",
    x = "LAD (m²/m³)",
    y = "Height Layer (index)"
  ) +
  theme_bw() +
  theme(legend.position = "none") # Hide legend as there are too many plots

print(p_profiles)
# ggsave("Plots/Sensitivity/4_LAD_Profiles_Raw.png", p_profiles, width = 6, height = 8)

# --- D. Link Structure to Sensitivity ---
# We merge the profile metrics with the Slope results from Analysis 1
# 'slope_data' comes from the previous block (Analysis 1)

# Filter slope_data to keep only one scenario (e.g., the best performing or the standard LAI)
# OR we can look at the "Range" of sensitivity across LAI scenarios
slope_sensitivity_metrics <- slope_data %>%
  filter(sim_type == "LAI_4") %>% # Example: Picking LAI_4 as representative
  left_join(profile_metrics, by = "id_plot")

# Plot 2: How does Structure affect Buffering (Slope)?
# We use Center of Gravity to see if "Top-heavy" vs "Bottom-heavy" profiles behave differently
p_structure_effect <- ggplot(slope_sensitivity_metrics, aes(x = Center_Gravity, y = slope_sim)) +
  geom_point(aes(size = Total_PAI, color = Total_PAI), alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  scale_color_viridis_c(option = "viridis", name = "Total PAI") +
  labs(
    title = "Impact of Vertical Structure on Buffering",
    subtitle = "Are top-heavy canopies (High Center of Gravity) better buffers?",
    x = "Center of Gravity (Height Index)",
    y = "Simulated Slope (Micro vs Macro)",
    size = "Total PAI"
  ) +
  theme_bw()

print(p_structure_effect)
# ggsave("Plots/Sensitivity/5_Structure_vs_Slope.png", p_structure_effect, width = 8, height = 6)

# --- E. Advanced Trend: Sensitivity Range vs Structure ---
# Does the LAI parameter choice matter MORE for dense plots or sparse plots?

# Calculate the variability (Standard Deviation) of the slope across all LAI 1-8 scenarios per plot
slope_variability <- slope_data %>%
  filter(sim_type != "HOBO_Ref" & sim_type != "Lidar_Constant") %>% # Only LAI scenarios
  group_by(id_plot) %>%
  summarise(
    Slope_SD = sd(slope_sim),      # How much does the slope change when LAI changes?
    Slope_Range = max(slope_sim) - min(slope_sim),
    .groups = "drop"
  ) %>%
  left_join(profile_metrics, by = "id_plot")

p_sensitivity_structure <- ggplot(slope_variability, aes(x = Total_PAI, y = Slope_Range)) +
  geom_point(color = "darkred", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  labs(
    title = "Model Sensitivity vs Canopy Density",
    subtitle = "Do denser plots have higher uncertainty when changing LAI parameter?",
    x = "Total PAI (from LAD profile)",
    y = "Range of Slope Variation (Max - Min across LAI 1-8)"
  ) +
  theme_bw()

print(p_sensitivity_structure)
# ggsave("Plots/Sensitivity/6_Sensitivity_Range_vs_PAI.png", p_sensitivity_structure, width = 8, height = 6)

# ------------------------------------------------------------------------------
# --- ANALYSIS 5: Sensitivity Range vs Structure ---
# ------------------------------------------------------------------------------
print("Calculating Sensitivity Range across LAI scenarios...")

# 1. Calculate the Variability of the Slope (The "Range")
# We look at how much the buffering capacity changes when we tweak LAI from 1 to 8.
slope_variability <- slope_data %>%
  # Filter to keep ONLY the LAI sensitivity scenarios (exclude Constant & HOBO)
  filter(grepl("LAI_", sim_type)) %>% 
  group_by(id_plot) %>%
  summarise(
    # The Gap between the strongest and weakest buffering predicted by the model
    Slope_Range = max(slope_sim) - min(slope_sim),
    
    # Standard Deviation (another measure of stability)
    Slope_SD = sd(slope_sim),
    
    # Keep the mean slope just for context
    Slope_Mean = mean(slope_sim),
    .groups = "drop"
  )

# 2. Merge with Structural Metrics (PAI, Center of Gravity)
sensitivity_vs_structure <- slope_variability %>%
  left_join(profile_metrics, by = "id_plot")

# ------------------------------------------------------------------------------
# --- Plot A: Sensitivity vs Canopy Density (PAI) ---
# Question: Is the model more "nervous" (sensitive) in dense or sparse forests?
# ------------------------------------------------------------------------------
p_sens_pai <- ggplot(sensitivity_vs_structure, aes(x = Total_PAI, y = Slope_Range)) +
  # Points colored by their vertical center of gravity
  geom_point(aes(color = Center_Gravity), size = 3, alpha = 0.8) +
  
  # Trend line
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  
  scale_color_viridis_c(option = "magma", name = "Center of Gravity\n(Height Index)") +
  
  labs(
    title = "Model Sensitivity vs. Canopy Density (PAI)",
    subtitle = "Y-axis: Difference in buffering between LAI_1 and LAI_8.\nHigh Y = The choice of LAI parameter is CRITICAL for these plots.",
    x = "Total PAI (Leaf Area Density Sum)",
    y = "Sensitivity Range (Max Slope - Min Slope)"
  ) +
  theme_bw()

print(p_sens_pai)
# ggsave("Plots/Sensitivity/5_Sensitivity_vs_PAI.png", p_sens_pai, width = 8, height = 6)

# ------------------------------------------------------------------------------
# --- Plot B: Sensitivity vs Vertical Arrangement (Center of Gravity) ---
# Question: Is the model more sensitive when biomass is high up vs low down?
# ------------------------------------------------------------------------------
p_sens_cog <- ggplot(sensitivity_vs_structure, aes(x = Center_Gravity, y = Slope_Range)) +
  geom_point(aes(size = Total_PAI), color = "#2171B5", alpha = 0.7) +
  
  geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se = FALSE) +
  
  labs(
    title = "Model Sensitivity vs. Vertical Structure",
    subtitle = "Does the height of the vegetation affect model stability?",
    x = "Center of Gravity (Height Index)",
    y = "Sensitivity Range (Max Slope - Min Slope)",
    size = "Total PAI"
  ) +
  theme_bw()

print(p_sens_cog)
# ggsave("Plots/Sensitivity/6_Sensitivity_vs_COG.png", p_sens_cog, width = 8, height = 6)

# ------------------------------------------------------------------------------
# --- Plot C: Mean Buffering & Uncertainty vs Structure ---
# ------------------------------------------------------------------------------

p_mean_error <- ggplot(sensitivity_vs_structure, aes(x = Total_PAI, y = Slope_Mean)) +
  geom_errorbar(aes(ymin = Slope_Mean - Slope_SD, 
                    ymax = Slope_Mean + Slope_SD), 
                width = 0.1, color = "grey60") +
  geom_point(aes(color = Center_Gravity), size = 3) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed", size = 0.8) +
  scale_color_viridis_c(option = "viridis", name = "Gravity Center\n(Height)") +
  labs(
    title = "Mean Buffering Capacity and Incertainty vs PAI",
    subtitle = "Points = Mean of LAI scnearios. Bars = Standard deviation.",
    x = "Plant Area Index (PAI)",
    y = "Mean Slope Micro/Macro (± SD)"
  ) +
  theme_bw()

print(p_mean_error)
# ggsave("Plots/Sensitivity/5b_Mean_SD_vs_PAI.png", p_mean_error, width = 8, height = 6)

# ------------------------------------------------------------------------------
# --- STATS: Correlations ---
# ------------------------------------------------------------------------------
# Quick check to see which structural factor drives sensitivity the most
cor_pai <- cor(sensitivity_vs_structure$Slope_Range, sensitivity_vs_structure$Total_PAI, use = "complete.obs")
cor_cog <- cor(sensitivity_vs_structure$Slope_Range, sensitivity_vs_structure$Center_Gravity, use = "complete.obs")

print(paste("Correlation (Sensitivity ~ PAI):", round(cor_pai, 3)))
print(paste("Correlation (Sensitivity ~ CoG):", round(cor_cog, 3)))

# ------------------------------------------------------------------------------
# --- ANALYSIS 6: Sensitivity of Delta Tmax (Cooling Capacity) ---
# ------------------------------------------------------------------------------
print("Calculating Sensitivity for Delta Tmax...")

# --- A. Data Preparation ---

# 1. Calculate Daily Max per Plot/Scenario
daily_delta <- combined_sensi %>%
  mutate(date = as.Date(time)) %>%
  group_by(sim_type, id_plot, date) %>%
  summarise(
    Tmax_sim = max(Tair_sim, na.rm = TRUE),
    Tmax_ign = max(Tair_ign, na.rm = TRUE), # Macro reference
    .groups = "drop"
  ) %>%
  mutate(
    # Delta Tmax = Micro - Macro
    # This represents the "Cooling Effect" (usually negative)
    Delta_Tmax = Tmax_sim - Tmax_ign
  )

# 2. Aggregate to Plot Level (Mean Cooling over the period)
# We get one "Average Cooling Value" per plot for each scenario
plot_delta_stats <- daily_delta %>%
  group_by(sim_type, id_plot) %>%
  summarise(
    Avg_Delta_Tmax = mean(Delta_Tmax, na.rm = TRUE),
    .groups = "drop"
  )

# 3. Calculate Sensitivity across LAI 1-8 Scenarios
delta_variability <- plot_delta_stats %>%
  filter(grepl("LAI_", sim_type)) %>% # Only keep LAI scenarios
  group_by(id_plot) %>%
  summarise(
    # How much does the cooling prediction change between LAI 1 and 8?
    Metric_Range = max(Avg_Delta_Tmax) - min(Avg_Delta_Tmax),
    
    # Standard deviation of the cooling prediction
    Metric_SD = sd(Avg_Delta_Tmax),
    
    # Average cooling predicted across all scenarios
    Metric_Mean = mean(Avg_Delta_Tmax),
    .groups = "drop"
  ) %>%
  # Merge with Structure Metrics (PAI, Center of Gravity)
  left_join(profile_metrics, by = "id_plot")

# --- B. Plotting ---

# Plot A: Sensitivity vs PAI
# ------------------------------------------------------------------------------
p_delta_pai <- ggplot(delta_variability, aes(x = Total_PAI, y = Metric_Range)) +
  geom_point(aes(color = Center_Gravity), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  scale_color_viridis_c(option = "magma", name = "Gravity Center\n(Height)") +
  labs(
    title = "Sensitivity of Delta Tmax vs. Canopy Density (PAI)",
    subtitle = "Y-axis: Range of variation in predicted Cooling (Delta Tmax) across LAI scenarios.\nHigh Y = High Uncertainty on the cooling effect.",
    x = "Total PAI",
    y = "Range of Delta Tmax (°C)"
  ) +
  theme_bw()

print(p_delta_pai)
# ggsave("Plots/Sensitivity/6a_DeltaTmax_Sens_vs_PAI.png", p_delta_pai, width = 8, height = 6)

# Plot B: Sensitivity vs Center of Gravity
# ------------------------------------------------------------------------------
p_delta_cog <- ggplot(delta_variability, aes(x = Center_Gravity, y = Metric_Range)) +
  geom_point(aes(size = Total_PAI), color = "#2171B5", alpha = 0.7) +
  geom_smooth(method = "lm", color = "darkred", linetype = "dashed", se = FALSE) +
  labs(
    title = "Sensitivity of Delta Tmax vs. Vertical Structure",
    subtitle = "Is the cooling prediction more uncertain for bottom-heavy or top-heavy plots?",
    x = "Center of Gravity (Height Index)",
    y = "Range of Delta Tmax (°C)",
    size = "Total PAI"
  ) +
  theme_bw()

print(p_delta_cog)
# ggsave("Plots/Sensitivity/6b_DeltaTmax_Sens_vs_COG.png", p_delta_cog, width = 8, height = 6)

# Plot C: Mean Cooling & Uncertainty vs PAI
# ------------------------------------------------------------------------------
p_delta_mean_error <- ggplot(delta_variability, aes(x = Total_PAI, y = Metric_Mean)) +
  
  # Error bars (Uncertainty)
  geom_errorbar(aes(ymin = Metric_Mean - Metric_SD, 
                    ymax = Metric_Mean + Metric_SD), 
                width = 0.1, color = "grey60") +
  
  # Points (Mean Performance)
  geom_point(aes(color = Center_Gravity), size = 3) +
  
  # Trend line
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed", size = 0.8) +
  
  scale_color_viridis_c(option = "viridis", name = "Gravity Center\n(Height)") +
  
  labs(
    title = "Mean Cooling Capacity (Delta Tmax) and Uncertainty vs PAI",
    subtitle = "Points = Mean Delta Tmax. Bars = Standard Deviation across LAI scenarios.\nLower values (more negative) = Stronger Cooling.",
    x = "Plant Area Index (PAI)",
    y = "Mean Delta Tmax (Micro - Macro) [°C]"
  ) +
  theme_bw()

print(p_delta_mean_error)
# ggsave("Plots/Sensitivity/6c_Mean_DeltaTmax_vs_PAI.png", p_delta_mean_error, width = 8, height = 6)

# --- Stats Check ---
cor_pai_delta <- cor(delta_variability$Metric_Range, delta_variability$Total_PAI, use = "complete.obs")
cor_cog_delta <- cor(delta_variability$Metric_Range, delta_variability$Center_Gravity, use = "complete.obs")

print(paste("Correlation (DeltaT Sensitivity ~ PAI):", round(cor_pai_delta, 3)))
print(paste("Correlation (DeltaT Sensitivity ~ CoG):", round(cor_cog_delta, 3)))

# ------------------------------------------------------------------------------
# --- ANALYSIS 7: Single Plot Time Series (Impact of LAI) ---
# ------------------------------------------------------------------------------

# --- 1. Plot Selection Strategy ---

# Option A: Select the "Median" plot (Representative of the dataset)
# We find the plot with the Total PAI closest to the median PAI of all plots.
target_val <- median(profile_metrics$Total_PAI, na.rm = TRUE)
chosen_plot_info <- profile_metrics %>%
  mutate(diff = abs(Total_PAI - target_val)) %>%
  arrange(diff) %>%
  slice(1)

target_id <- chosen_plot_info$id_plot
print(paste("Selected Plot (Median PAI Profile):", target_id))
print(paste("Plot PAI:", round(chosen_plot_info$Total_PAI, 2)))

# Option B: Force a specific plot ID (Uncomment to use)
# target_id <- "41_13" 

# --- 2. Data Preparation for the Selected Plot ---

single_plot_data <- combined_sensi %>%
  filter(id_plot == target_id) %>%
  # Keep only LAI scenarios and the HOBO reference
  filter(grepl("LAI_|HOBO", sim_type)) 

# Ensure factor order for clean legend (LAI 1 -> 8 -> HOBO)
level_order <- c(paste0("LAI_", 1:8), "HOBO_Ref")
single_plot_data$sim_type <- factor(single_plot_data$sim_type, 
                                    levels = level_order)

# Create custom color palette
# LAI scenarios in Plasma gradient (Orange -> Purple), HOBO in solid Black
n_lai <- 8
lai_colors <- viridis::plasma(n_lai, end = 0.85) # end=0.85 avoids very light yellows
names(lai_colors) <- paste0("LAI_", 1:8)
custom_colors <- c(lai_colors, "HOBO_Ref" = "black")

zoom_start <- min(single_plot_data$time)
zoom_end   <- max(single_plot_data$time)

# --- 3. Visualization: Heatwave Zoom (Absolute Temperatures) ---

# 1. Prepare Simulation Data (LAI Scenarios)
data_zoom_sim <- single_plot_data %>%
  filter(time >= zoom_start & time <= zoom_end) %>%
  filter(grepl("LAI_", sim_type)) # Ensure we only have LAI scenarios here

# 2. Prepare Observation Data (HOBO)
# We extract unique observation values to draw a single black line
data_zoom_obs <- data_zoom_sim %>%
  dplyr::select(time, Tair_obs) %>%
  distinct()

data_zoom <- single_plot_data %>%
  filter(time >= zoom_start & time <= zoom_end)

# 3. Plot
p_timeseries_zoom <- ggplot() +
  
  # A. Simulation Lines (Colored by LAI)
  geom_line(data = data_zoom_sim, 
            aes(x = time, y = Tair_sim, color = sim_type), 
            linewidth = 0.7, alpha = 0.8) +
  
  # B. HOBO Observation Line (Black & Dashed)
  geom_line(data = data_zoom_obs, 
            aes(x = time, y = Tair_obs), 
            color = "black", linewidth = 1.2, linetype = "dashed") +
  
  # C. Colors & Labels
  scale_color_manual(values = lai_colors, name = "Scenario") +
  
  labs(
    title = paste("Impact of LAI on Temporal Dynamics (Plot:", target_id, ")"),
    subtitle = "Zoom on a June week. Dashed Black = HOBO (Observed). Colors = LAI Scenarios.",
    x = NULL,
    y = "Temperature (°C)"
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p_timeseries_zoom)
# ggsave(paste0("Plots/Sensitivity/7_TimeSeries_Zoom_", target_id, ".png"), 
#        p_timeseries_zoom, width = 10, height = 6)

# --- 4. Visualization: Instantaneous Bias (Sim - Obs) ---

# Since 'data_zoom' comes from combined_sensi, it should have both columns:
# Tair_sim (Simulation) and Tair_obs (Observation)

data_zoom_diff <- data_zoom %>%
  # Calculate bias directly (Sim - Obs)
  mutate(Bias = Tair_sim - Tair_obs) %>%
  # Select only what is needed for plotting
  dplyr::select(time, sim_type, Bias)

# Ensure factor order for the legend
data_zoom_diff$sim_type <- factor(data_zoom_diff$sim_type, levels = paste0("LAI_", 1:8))

# Plot
p_bias_zoom <- ggplot(data_zoom_diff, aes(x = time, y = Bias, color = sim_type)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) + # Reference line (Perfect fit)
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = lai_colors, name = "Scenario") +
  labs(
    title = "Instantaneous Bias (Simulation - Observation)",
    subtitle = "Bias > 0: Model overestimates T°. Bias < 0: Model underestimates T°.",
    x = NULL,
    y = "Bias (°C)"
  ) +
  theme_bw()

print(p_bias_zoom)
# ggsave(paste0("Plots/Sensitivity/7b_Bias_Zoom_", target_id, ".png"), 
#        p_bias_zoom, width = 10, height = 6)
