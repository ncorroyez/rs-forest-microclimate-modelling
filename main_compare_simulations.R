rm(list=ls(all=TRUE))
gc()

# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(musica.tools)
library(ncdf4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(purrr)
library(viridis)
library(scales)
library(ggridges)

# Create output directory for plots
plot_dir <- "./Plots/Sensitivity"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ==============================================================================
# 2. CONFIGURATION & CONSTANTS
# ==============================================================================

# Root directory of simulation results
results_dir <- './out_files/Sensitivity_Analysis'

# Simulation Parameters
target_heights <- c(10, 20, 30, 40)
profile_shapes <- c("uniform", 
                    "top_heavy", 
                    "bottom_heavy",
                    "top_bottom_heavy", 
                    "symmetric",
                    "top_peak_uniform",
                    "extreme_bottom",
                    "low_symmetric")
clumping_vals <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1) 
lai_values     <- seq(1, 8, by = 1)
forcing_type   <- "ERA5"

# Study Period (Summer 2021)
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

#' Generate synthetic vertical LAD profile on a fixed grid
#' 
#' @param height_target Target canopy height (e.g., 10, 20, 30, 40)
#' @param type String: "uniform", "top_heavy", "bottom_heavy", "symmetric"
#' @return A dataframe with height (z) and normalized density
get_synthetic_profile <- function(height_target, type) {
  
  # Create vector from 2.5 to (Height - 0.5) with step of 1
  # Example for 10m: 2.5, 3.5, ..., 9.5
  z_levels <- seq(2.5, height_target - 0.5, by = 1)
  n <- length(z_levels)
  
  # Normalized height for Beta distribution calculation (0 to 1)
  x_seq <- seq(0.01, 0.99, length.out = n) 
  
  if (type == "uniform") {
    dens <- rep(1, n)
  } else if (type == "top_heavy") {
    # Skewed towards top
    dens <- dbeta(x_seq, shape1 = 6, shape2 = 2.5)
  } else if (type == "bottom_heavy") {
    # Skewed towards bottom
    dens <- dbeta(x_seq, shape1 = 2.5, shape2 = 6)
  } else if (type == "symmetric") {
    # Gaussian-like
    dens <- dbeta(x_seq, shape1 = 5, shape2 = 5)
  } else if (type == "top_bottom_heavy") {
    # Bimodal: Sum of Top Heavy and Bottom Heavy distributions
    dens <- dbeta(x_seq, shape1 = 4, shape2 = 20) + dbeta(x_seq, shape1 = 20, shape2 = 5)
  } else if (type == "top_peak_uniform") {
    # Description: Peak at the very top, roughly uniform below.
    # Logic: A strong top-skewed Beta added to a constant baseline (uniform).
    # shape1=15, shape2=2 pushes the peak to the very top.
    # + 0.5 adds the "uniform" base layer.
    dens <- dbeta(x_seq, shape1 = 15, shape2 = 2) + 0.5
  } else if (type == "extreme_bottom") {
    # Description: Weak all along except huge increase at the bottom.
    # Logic: A very strong bottom-skewed Beta with a tiny baseline.
    # shape1=1.2 (pushes peak near 0), shape2=15 (decays very fast).
    # + 0.05 represents the "weak" vegetation along the trunk/branches.
    dens <- dbeta(x_seq, shape1 = 1.2, shape2 = 15) + 0.05
  } else if (type == "low_symmetric") {
    # "Ventru vers le bas": Gaussian-like but peaked at ~35% height instead of 50%
    # Seen in plots 41_57, 41_08
    dens <- dbeta(x_seq, shape1 = 3, shape2 = 5)
    
  } else {
    stop(paste("Unknown profile type:", type))
  }
  
  # Normalize to sum = 1 (Density shape only)
  if (sum(dens) > 0) {
    dens <- dens / sum(dens)
  }
  
  return(data.frame(height = z_levels, density = dens))
}

# Function to calculate Center of Gravity (CoG)
calc_cog <- function(shape_name) {
  # We use 30m as the standard reference height
  prof <- get_synthetic_profile(30, shape_name) 
  return(sum(prof$height * prof$density) / sum(prof$density))
}

# ==============================================================================
# 4. DATA LOADING: MACROCLIMATE (ERA5 REFERENCE)
# ==============================================================================
cat("--> Loading ERA5 Macroclimate data (musica_in_Blois.nc)...\n")

era5_file <- "in_files/musica_in_Blois.nc"
if (!file.exists(era5_file)) stop(paste("File not found:", era5_file))

# 1. Extraction
era5 <- nc_open(era5_file)
era5_var <- get_variable(era5, "Tair")
nc_close(era5)

# 2. Cleaning & Timezone Forcing
df_macro <- era5_var %>%
  rename(Tair_macro = Tair) %>%
  mutate(
    Tair_macro = Tair_macro - 273.15, # Kelvin -> Celsius
    time = as.POSIXct(time, tz = "UTC")
  ) %>%
  filter(as.Date(time) %in% date_seq) %>%
  mutate(time = floor_date(time, "hour")) %>%
  # CRITICAL: Force UTC to match Simulation data later
  mutate(time = force_tz(time, "UTC")) %>% 
  group_by(time) %>%
  summarise(Tair_macro = mean(Tair_macro, na.rm = TRUE), .groups = "drop")

cat(paste("Macroclimate (ERA5) loaded:", nrow(df_macro), "hours.\n"))

# ==============================================================================
# 5. DATA LOADING: MICROCLIMATE (SYNTHETIC SIMULATIONS)
# ==============================================================================
cat("--> Loading Synthetic Simulation data (All Parameters)...\n")

df_sim_all <- data.frame()

for (h in target_heights) {
  for (shape in profile_shapes) {
    for (clump in clumping_vals) {
      for (lai in lai_values) {
        
        # Path must match the simulation output structure
        nc_path <- file.path(results_dir, 
                             paste0("height_", h, "m"),
                             shape,
                             paste0("clump_", clump),
                             paste0("lai_", lai),
                             paste0("musica_out_", 
                                    forcing_type, "_synthetic.nc"))
        
        if (file.exists(nc_path)) {
          out <- tryCatch(nc_open(nc_path), error = function(e) NULL)
          
          if (!is.null(out)) {
            # Extract Tair_z at ground level (nair = 1)
            raw_data <- tryCatch(get_variable(out, "Tair_z"), 
                                 error = function(e) NULL)
            
            if (!is.null(raw_data)) {
              df_chunk <- raw_data %>%
                filter(as.Date(time) %in% date_seq, nair == 1) %>% 
                mutate(
                  Tair_micro = Tair_z - 273.15,
                  Height_Scenario = h,
                  Shape_Scenario = shape,
                  Clump_Scenario = clump, # Capture Clumping
                  LAI_Scenario = lai,
                  # Unique ID for grouping
                  Run_ID = paste(h, shape, clump, lai, sep="_")
                ) %>%
                dplyr::select(time, 
                              Tair_micro, 
                              Height_Scenario, 
                              Shape_Scenario, 
                              Clump_Scenario, 
                              LAI_Scenario, 
                              Run_ID)
              
              df_sim_all <- bind_rows(df_sim_all, df_chunk)
            }
            nc_close(out)
          }
        }
      }
    }
  }
  cat(paste("Finished processing height:", h, "m\n"))
}

# Aggregate to Hourly with Time Shift & Timezone Fix
df_sim_hourly <- df_sim_all %>%
  # 1. Apply time shift (-2h)
  mutate(time = time - hours(2)) %>%
  # 2. Floor to nearest hour
  mutate(time = floor_date(time, unit = "hour")) %>%
  # 3. CRITICAL: Force UTC to match Macroclimate
  mutate(time = force_tz(time, "UTC")) %>%
  group_by(Run_ID, 
           Height_Scenario, 
           Shape_Scenario, 
           Clump_Scenario,
           LAI_Scenario, 
           time
  ) %>%
  summarise(Tair_micro = mean(Tair_micro, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# 6. MERGE DATASETS
# ==============================================================================
cat("--> Merging datasets...\n")

df_combined <- inner_join(df_macro, df_sim_hourly, by = "time") %>%
  mutate(
    # Cooling Effect: Negative values = Micro is cooler than Macro
    Delta_T = Tair_micro - Tair_macro
  )

# Define Focus Dataset (30m Standard) for most comparisons
df_focus_30m <- df_combined %>% filter(Height_Scenario == 10)

cat("Data Ready. Starting Analysis...\n")

# ==============================================================================
# SECTION A: CONTEXT & PHYSICS (VALIDATION)
# ==============================================================================

# PLOT A0: Visualizing the Vertical Profiles (LAD)
# ------------------------------------------------------------------------------
# Context: What do the simulated trees actually look like?
cat("Plot A0: Vertical Profiles...\n")

df_profiles_viz <- data.frame()
for (s in profile_shapes) {
  # Generate for 30m example
  tmp <- get_synthetic_profile(30, s)
  tmp$Shape_Scenario <- s
  df_profiles_viz <- rbind(df_profiles_viz, tmp)
}

p_a0 <- ggplot(df_profiles_viz, aes(x = density, y = height, color = Shape_Scenario)) +
  geom_path(size = 2, alpha = 0.8) +
  scale_color_viridis_d(option = "turbo", name = "Profile Shape") +
  labs(
    title = "Vertical Leaf Area Density (LAD) Profiles",
    subtitle = "Normalized shape distribution for a 30m canopy.",
    x = "Relative Leaf Density",
    y = "Height (m)"
  ) +
  theme_bw()

print(p_a0)
ggsave(file.path(plot_dir, "A0_Vertical_Profiles.png"), p_a0, width = 6, height = 6)

# 1. Generate Data for All Heights and Shapes
# ------------------------------------------------------------------------------
target_heights <- c(10, 20, 30, 40)
df_profiles_all <- data.frame()

for (h in target_heights) {
  for (s in profile_shapes) {
    # Generate profile for specific height h
    tmp <- get_synthetic_profile(h, s)
    tmp$Shape_Scenario <- s
    tmp$Height_Scenario <- h
    
    # Create a nice label for faceting
    tmp$Height_Label <- paste0(h, "m Canopy")
    
    df_profiles_all <- rbind(df_profiles_all, tmp)
  }
}

# Ensure factors are ordered correctly for plotting
df_profiles_all$Height_Label <- factor(df_profiles_all$Height_Label, 
                                       levels = paste0(target_heights, "m Canopy"))

# ------------------------------------------------------------------------------
# VARIATION 1: WRAP BY HEIGHT (Compare Shapes)
# ------------------------------------------------------------------------------
# Goal: For a 10m tree, how do the shapes differ?
p_a0_by_height <- ggplot(df_profiles_all, aes(x = density, y = height, color = Shape_Scenario)) +
  geom_path(size = 1.5, alpha = 0.8) +
  
  # Facet by Height, allowing Y-axis to scale to the tree size
  facet_wrap(~Height_Label, scales = "free_y", ncol = 4) +
  
  scale_color_viridis_d(option = "turbo", name = "Profile Shape") +
  labs(
    title = "Vertical Profiles: Comparison by Canopy Height",
    subtitle = "Normalized Leaf Area Density (LAD) distribution.",
    x = "Relative Leaf Density",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_a0_by_height)
ggsave(file.path(plot_dir, "A0_1_Profiles_Wrap_Height.png"), p_a0_by_height, width = 12, height = 6)


# ------------------------------------------------------------------------------
# VARIATION 2: WRAP BY PROFILE SHAPE (Compare Heights)
# ------------------------------------------------------------------------------
# Goal: How does a 'Top Heavy' profile look at 10m vs 40m?
p_a0_by_shape <- ggplot(df_profiles_all, aes(x = density, y = height, color = as.factor(Height_Scenario))) +
  geom_path(size = 1.5, alpha = 0.8) +
  
  # Facet by Shape
  facet_wrap(~Shape_Scenario, ncol = 4) +
  
  # Use Magma palette to distinguish heights (Yellow = Tall, Dark = Short)
  scale_color_viridis_d(option = "magma", end = 0.9, name = "Height (m)") +
  
  labs(
    title = "Vertical Profiles: Comparison by Shape Architecture",
    subtitle = "Evolution of the profile geometry across different target heights.",
    x = "Relative Leaf Density",
    y = "Height (m)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_a0_by_shape)
ggsave(file.path(plot_dir, "A0_2_Profiles_Wrap_Shape.png"), p_a0_by_shape, width = 12, height = 6)

# PLOT A1: Intrinsic Buffering Capacity (Micro vs Macro)
# ------------------------------------------------------------------------------
# Physics: Does the model dampen temperature extremes?
cat("Plot A1: Buffering Capacity...\n")

df_buffer <- df_focus_30m %>% filter(LAI_Scenario == 5)

p_a1 <- ggplot(df_buffer, aes(x = Tair_macro, y = Tair_micro, color = Shape_Scenario)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  scale_color_viridis_d(option = "turbo", name = "Shape") +
  labs(
    title = "Intrinsic Buffering Capacity (Input vs Output)",
    subtitle = "Height 30m, LAI 5. Slope < 1 indicates active buffering.",
    x = "ERA5 Forcing T° (°C)", y = "Simulated Understory T° (°C)"
  ) +
  coord_fixed(xlim=c(10,35), ylim=c(10,35)) + theme_bw()

print(p_a1)
ggsave(file.path(plot_dir, "A1_Buffering_Capacity.png"), p_a1, width = 8, height = 8)


# ==============================================================================
# SECTION B: STRUCTURAL DRIVERS (SHAPE IMPACT)
# ==============================================================================

# PLOT B1: Center of Gravity Impact
# ------------------------------------------------------------------------------
# Hypothesis: Higher biomass interception (Higher CoG) = Better Cooling.
cat("Plot B1: Center of Gravity...\n")

df_cog_vals <- data.frame(Shape_Scenario = profile_shapes, Center_Gravity = sapply(profile_shapes, calc_cog))

metrics_cog <- df_focus_30m %>%
  filter(LAI_Scenario == 8) %>%
  mutate(date = as.Date(time)) %>%
  group_by(Shape_Scenario, date) %>%
  summarise(Daily_Delta = min(Delta_T), .groups="drop") %>%
  group_by(Shape_Scenario) %>%
  summarise(Mean_Max_Cooling = mean(Daily_Delta), .groups="drop") %>%
  left_join(df_cog_vals, by="Shape_Scenario")

p_b1 <- ggplot(metrics_cog, aes(x = Center_Gravity, y = Mean_Max_Cooling)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey", linetype = "dashed") +
  geom_point(aes(color = Shape_Scenario), size = 6) +
  geom_text(aes(label = Shape_Scenario), vjust = -1.5, fontface = "bold") +
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "Impact of Biomass Height on Cooling",
    # subtitle = "Does intercepting light higher up (High CoG) cool the ground better? (LAI=5)",
    x = "Center of Gravity Height (m)", y = "Mean Max Cooling (°C)"
  ) +
  theme_bw() + theme(legend.position = "none")

print(p_b1)
ggsave(file.path(plot_dir, "B1_CoG_Impact.png"), p_b1, width = 8, height = 6)


# PLOT B2: Structural Anomaly
# ------------------------------------------------------------------------------
# Comparison: How much warmer/cooler is a shape compared to the 'Uniform' baseline?
cat("Plot B2: Structural Anomaly...\n")

stats_abs <- df_focus_30m %>%
  mutate(date = as.Date(time)) %>%
  group_by(Shape_Scenario, LAI_Scenario, date) %>%
  summarise(Daily_Max = max(Tair_micro), .groups="drop") %>%
  group_by(Shape_Scenario, LAI_Scenario) %>%
  summarise(Avg_Max = mean(Daily_Max), .groups="drop")

df_anomaly <- stats_abs %>%
  pivot_wider(names_from = Shape_Scenario, values_from = Avg_Max) %>%
  mutate(
    Diff_Top = top_heavy - uniform,
    Diff_Bottom = bottom_heavy - uniform,
    Diff_Sym = symmetric - uniform
  ) %>%
  select(LAI_Scenario, Diff_Top, Diff_Bottom, Diff_Sym) %>%
  pivot_longer(cols = starts_with("Diff"), names_to = "Comp", values_to = "Val")

p_b2 <- ggplot(df_anomaly, aes(x = LAI_Scenario, y = Val, fill = Comp)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#E69F00", "#999999", "#56B4E9")) +
  labs(
    title = "Structural Anomaly (vs Uniform Profile)",
    subtitle = "Positive = Warmer than Uniform. Negative = Cooler.",
    x = "LAI", y = "Temperature Anomaly (°C)"
  ) + theme_bw()

print(p_b2)
ggsave(file.path(plot_dir, "B2_Structural_Anomaly.png"), p_b2, width = 8, height = 6)


# ==============================================================================
# SECTION C: DENSITY DRIVERS (LAI IMPACT)
# ==============================================================================

# PLOT C1: Cooling vs LAI Curve
# ------------------------------------------------------------------------------
# Efficiency: Diminishing returns of adding leaves.
cat("Plot C1: Cooling Efficiency...\n")

stats_cool <- df_focus_30m %>%
  mutate(date = as.Date(time)) %>%
  group_by(Shape_Scenario, LAI_Scenario, date) %>%
  summarise(D_Tmax = max(Tair_micro) - max(Tair_macro), .groups="drop") %>%
  group_by(Shape_Scenario, LAI_Scenario) %>%
  summarise(Mean_D_Tmax = mean(D_Tmax), .groups="drop")

p_c1 <- ggplot(stats_cool, aes(x = LAI_Scenario, y = Mean_D_Tmax, color = Shape_Scenario)) +
  geom_hline(yintercept = 0) +
  geom_line(size=1.2) + geom_point(size=3) +
  scale_color_viridis_d(option="turbo") +
  scale_x_continuous(breaks = 1:8) +
  labs(title="Cooling Efficiency vs LAI", subtitle="Reduction of Daily Heat Peak (Micro - ERA5)", 
       x="LAI", y="Mean Delta Tmax (°C)") + theme_bw()

print(p_c1)
ggsave(file.path(plot_dir, "C1_Cooling_Efficiency.png"), p_c1, width = 8, height = 6)


# PLOT C2: Log Slope (Decoupling)
# ------------------------------------------------------------------------------
# Metric: How independent is the microclimate from the macroclimate?
cat("Plot C2: Log Slope...\n")

stats_slope <- df_focus_30m %>%
  group_by(Shape_Scenario, LAI_Scenario) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(Tair_micro ~ Tair_macro, data=.)),
    Slope = map_dbl(model, ~coef(.)[2]),
    Log_Slope = log(Slope)
  ) %>% select(-data, -model)

p_c2 <- ggplot(stats_slope, aes(x=LAI_Scenario, y=Log_Slope, color=Shape_Scenario)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_line(size=1.2) + geom_point(size=3) +
  scale_color_viridis_d(option="turbo") +
  scale_x_continuous(breaks = 1:8) +
  labs(title="Thermal Decoupling (Log Slope)", subtitle="More negative = Stronger Isolation/Buffering",
       x="LAI", y="Log(Slope)") + theme_bw()

print(p_c2)
ggsave(file.path(plot_dir, "C2_LogSlope.png"), p_c2, width = 8, height = 6)


# PLOT C3: Sensitivity Envelope (The Ribbon)
# ------------------------------------------------------------------------------
# Uncertainty: How much does LAI choice matter for a given shape?
cat("Plot C3: Sensitivity Envelope...\n")

# Zoom on a heatwave for clarity
zoom_env <- df_focus_30m %>% filter(time >= as.Date("2021-06-14") & time <= as.Date("2021-06-18"))

stats_env <- zoom_env %>%
  group_by(Shape_Scenario, time) %>%
  summarise(T_Mean=mean(Tair_micro), T_Min=min(Tair_micro), T_Max=max(Tair_micro), 
            T_Macro=mean(Tair_macro), .groups="drop")

p_c3 <- ggplot(stats_env, aes(x=time)) +
  geom_line(aes(y=T_Macro), linetype="dashed", color="black") +
  geom_ribbon(aes(ymin=T_Min, ymax=T_Max, fill=Shape_Scenario), alpha=0.4) +
  geom_line(aes(y=T_Mean, color=Shape_Scenario), size=1) +
  facet_wrap(~Shape_Scenario, ncol=1) +
  scale_fill_viridis_d(option="turbo") + scale_color_viridis_d(option="turbo") +
  labs(title="Sensitivity Envelope (LAI 1-8)", 
       subtitle="Ribbon width = Temperature uncertainty due to LAI choice.",
       y="Temperature (°C)") + theme_bw()

print(p_c3)
ggsave(file.path(plot_dir, "C3_Sensitivity_Envelope.png"), p_c3, width = 10, height = 10)

# ==============================================================================
# SECTION C: DENSITY DRIVERS (LAI IMPACT) - MULTI-HEIGHT ANALYSIS
# ==============================================================================

# 1. PREPARE DATA FOR C1 (Cooling Efficiency)
# ------------------------------------------------------------------------------
# We calculate the daily maximum cooling (Micro Max - Macro Max)
# for ALL heights, shapes, and LAI values.

cat("--> Processing Data for Plot C1 (Cooling Efficiency)...\n")

stats_cool_all <- df_combined %>%
  mutate(date = as.Date(time)) %>%
  # Calculate Daily Max for Micro and Macro
  group_by(Height_Scenario, Shape_Scenario, LAI_Scenario, date) %>%
  summarise(
    Tmax_micro = max(Tair_micro), 
    Tmax_macro = max(Tair_macro),
    .groups = "drop"
  ) %>%
  # Calculate Delta and Average over the period
  mutate(Delta_Tmax = Tmax_micro - Tmax_macro) %>%
  group_by(Height_Scenario, Shape_Scenario, LAI_Scenario) %>%
  summarise(
    Mean_Delta_Tmax = mean(Delta_Tmax), 
    .groups = "drop"
  ) %>%
  # Create nice labels for plotting
  mutate(
    Height_Label = paste0(Height_Scenario, "m Canopy"),
    Shape_Label  = paste0("Shape: ", Shape_Scenario)
  )

# PLOT C1 - VARIANT A: FACET BY SHAPE
# ------------------------------------------------------------------------------
# Question: For a specific Shape (e.g., Top Heavy), how does Tree Height affect cooling?
cat("Plot C1_A: Cooling Efficiency (Facet by Shape)...\n")

p_c1_a <- ggplot(stats_cool_all, aes(x = LAI_Scenario, y = Mean_Delta_Tmax, color = as.factor(Height_Scenario))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  
  # Facet by Shape
  facet_wrap(~Shape_Scenario, ncol = 3) +
  
  # Color by Height (Magma is good for sequential magnitude)
  scale_color_viridis_d(option = "magma", end = 0.9, direction = -1, name = "Tree Height (m)") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  
  labs(
    title = "C1-A: Cooling Efficiency vs LAI (grouped by Shape)",
    subtitle = "Comparing how tree height influences cooling for a given foliage distribution.",
    x = "Leaf Area Index (LAI)", 
    y = "Mean Delta Tmax (°C) [Micro - Macro]"
  ) + 
  theme_bw() +
  theme(legend.position = "bottom")

print(p_c1_a)
ggsave(file.path(plot_dir, "C1_A_Cooling_Facet_Shape.png"), p_c1_a, 
       width = 10, height = 6, dpi = 300)


# PLOT C1 - VARIANT B: FACET BY HEIGHT
# ------------------------------------------------------------------------------
# Question: For a specific Height (e.g., 20m), which Shape cools best?
cat("Plot C1_B: Cooling Efficiency (Facet by Height)...\n")

p_c1_b <- ggplot(stats_cool_all, aes(x = LAI_Scenario, y = Mean_Delta_Tmax, color = Shape_Scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  
  # Facet by Height
  facet_wrap(~Height_Label, ncol = 2) +
  
  # Color by Shape (Turbo is good for categorical distinction)
  scale_color_viridis_d(option = "turbo", name = "Profile Shape") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  
  labs(
    title = "C1-B: Cooling Efficiency vs LAI (grouped by Height)",
    subtitle = "Comparing performance of foliage distributions within specific canopy heights.",
    x = "Leaf Area Index (LAI)", 
    y = "Mean Delta Tmax (°C) [Micro - Macro]"
  ) + 
  theme_bw() +
  theme(legend.position = "bottom")

print(p_c1_b)
ggsave(file.path(plot_dir, "C1_B_Cooling_Facet_Height.png"), p_c1_b, 
       width = 7, height = 7, dpi = 300)


# 2. PREPARE DATA FOR C2 (Thermal Decoupling / Log Slope)
# ------------------------------------------------------------------------------
# We calculate the slope of the linear regression (Micro ~ Macro)
# Lower slope = Stronger decoupling/buffering.

cat("--> Processing Data for Plot C2 (Thermal Decoupling)...\n")

stats_slope_all <- df_combined %>%
  group_by(Height_Scenario, Shape_Scenario, LAI_Scenario) %>%
  nest() %>%
  mutate(
    # Fit Linear Model: Micro ~ Macro
    model = map(data, ~lm(Tair_micro ~ Tair_macro, data = .)),
    # Extract Slope coefficient
    Slope = map_dbl(model, ~coef(.)[2]),
    # Log transform for visualization (optional, but requested in original)
    Log_Slope = log(Slope)
  ) %>% 
  select(-data, -model) %>%
  # Create nice labels
  mutate(
    Height_Label = paste0(Height_Scenario, "m Canopy")
  )

# PLOT C2 - VARIANT A: FACET BY SHAPE
# ------------------------------------------------------------------------------
cat("Plot C2_A: Decoupling (Facet by Shape)...\n")

p_c2_a <- ggplot(stats_slope_all, aes(x = LAI_Scenario, y = Log_Slope, color = as.factor(Height_Scenario))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Slope 1 = No buffering
  geom_line(size = 1) + 
  geom_point(size = 2) +
  
  facet_wrap(~Shape_Scenario, ncol = 3) +
  
  scale_color_viridis_d(option = "magma", end = 0.9, direction = -1, name = "Tree Height (m)") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  
  labs(
    title = "C2-A: Thermal Decoupling (Log_Slope) vs LAI (grouped by Shape)",
    subtitle = "Log_Slope < 0 indicates buffering. Lower values = Stronger decoupling.",
    x = "Leaf Area Index (LAI)", 
    y = "Regression Log_Slope (Micro vs Macro)"
  ) + 
  theme_bw() +
  theme(legend.position = "bottom")

print(p_c2_a)
ggsave(file.path(plot_dir, "C2_A_Decoupling_Facet_Shape.png"), p_c2_a, 
       width = 10, height = 6, dpi = 300)


# PLOT C2 - VARIANT B: FACET BY HEIGHT
# ------------------------------------------------------------------------------
cat("Plot C2_B: Decoupling (Facet by Height)...\n")

p_c2_b <- ggplot(stats_slope_all, aes(x = LAI_Scenario, y = Log_Slope, color = Shape_Scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_line(size = 1) + 
  geom_point(size = 2) +
  
  facet_wrap(~Height_Label, ncol = 2) +
  
  scale_color_viridis_d(option = "turbo", name = "Profile Shape") +
  scale_x_continuous(breaks = seq(1, 8, 1)) +
  
  labs(
    title = "C2-B: Thermal Decoupling vs LAI (grouped by Height)",
    subtitle = "Comparing buffering capacity of different shapes within specific heights.",
    x = "Leaf Area Index (LAI)", 
    y = "Regression Log_Slope (Micro vs Macro)"
  ) + 
  theme_bw() +
  theme(legend.position = "bottom")

print(p_c2_b)
ggsave(file.path(plot_dir, "C2_B_Decoupling_Facet_Height.png"), p_c2_b, 
       width = 7, height = 7, dpi = 300)

cat("Analysis Complete. All 4 plots saved in:", plot_dir, "\n")


# ==============================================================================
# SECTION D: TEMPORAL DYNAMICS
# ==============================================================================

# PLOT D1: Hourly Power Heatmap
# ------------------------------------------------------------------------------
# Timing: When is the LAI most effective?
cat("Plot D1: Hourly Power...\n")

df_hour <- df_focus_30m %>%
  mutate(Hour = hour(time)) %>%
  group_by(Shape_Scenario, Hour, LAI_Scenario) %>%
  summarise(Mean_T = mean(Tair_micro), .groups="drop") %>%
  group_by(Shape_Scenario, Hour) %>%
  summarise(Power = Mean_T[LAI_Scenario==8] - Mean_T[LAI_Scenario==1], .groups="drop")

p_d1 <- ggplot(df_hour, aes(x=Hour, y=Shape_Scenario, fill=Power)) +
  geom_tile() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, name="LAI Impact") +
  scale_x_continuous(breaks=seq(0,23,2)) +
  labs(title="Hourly Sensitivity Heatmap", 
       subtitle="Temperature diff: Dense (LAI 8) - Sparse (LAI 1). Blue = Cooling.",
       x="Hour (UTC)") + theme_minimal()

print(p_d1)
ggsave(file.path(plot_dir, "D1_Hourly_Power.png"), p_d1, width = 8, height = 5)


# PLOT D2: Seasonal Day/Night Dynamics
# ------------------------------------------------------------------------------
# Seasonality: Evolution over summer months.
cat("Plot D2: Seasonal Dynamics...\n")

df_seas <- df_focus_30m %>%
  mutate(Month = month(time, label=TRUE), Hour = hour(time),
         Period = ifelse(Hour>=7 & Hour<=21, "Day", "Night")) %>%
  group_by(Shape_Scenario, LAI_Scenario, Month) %>%
  summarise(
    Delta_Day = mean(Delta_T[Period=="Day"], na.rm=TRUE),
    Delta_Night = mean(Delta_T[Period=="Night"], na.rm=TRUE), .groups="drop"
  ) %>%
  pivot_longer(cols=c(Delta_Day, Delta_Night), names_to="Per", values_to="Val")

p_d2 <- ggplot(df_seas, aes(x=LAI_Scenario, y=Val, color=Shape_Scenario)) +
  geom_hline(yintercept=0, linetype="dashed") + geom_line() +
  facet_grid(Per ~ Month, scales="free_y") +
  scale_color_viridis_d(option="turbo") +
  labs(title="Seasonal Day/Night Dynamics", y="Mean Delta T (°C)") + theme_bw()

print(p_d2)
ggsave(file.path(plot_dir, "D2_Seasonal_Dynamics.png"), p_d2, width = 12, height = 8)


# PLOT D3: Time Series Zoom
# ------------------------------------------------------------------------------
# Raw Data: Visual check on a heatwave event.
cat("Plot D3: Time Series Zoom...\n")

target_zoom_data <- df_focus_30m %>%
  filter(Shape_Scenario == "top_heavy") %>% # Focus on one shape to keep it readable
  filter(time >= as.Date("2021-06-12") & time <= as.Date("2021-06-19"))

# Define colors for LAI
lai_colors <- viridis::plasma(8, end = 0.9)
names(lai_colors) <- as.character(1:8)

p_e3 <- ggplot() +
  geom_line(data = target_zoom_data, 
            aes(x = time, y = Tair_micro, color = as.factor(LAI_Scenario), group = LAI_Scenario),
            size = 0.8, alpha = 0.8) +
  geom_line(data = target_zoom_data %>% select(time, Tair_macro) %>% distinct(),
            aes(x = time, y = Tair_macro), 
            color = "black", size = 1.2, linetype = "dashed") +
  scale_color_manual(values = lai_colors, name = "LAI") +
  labs(
    title = "Heatwave Dynamics Zoom (Top Heavy, 30m)",
    subtitle = "Black Dashed = ERA5 Macroclimate. Colored = Microclimate (LAI 1-8).",
    x = NULL, y = "Temperature (°C)"
  ) +
  theme_bw()

print(p_e3)
ggsave(file.path(plot_dir, "D3_Zoom_TimeSeries.png"), p_e3, width = 10, height = 6)

# SECTION E: FOCUS ANALYSIS - DISTRIBUTION DYNAMICS
# (Case Study: Impact of LAI on a specific Profile & Height)
# ------------------------------------------------------------------------------
target_h_focus <- 30
target_shape_focus <- "top_heavy" # Change to "uniform", "bottom_heavy", etc.

cat(paste("--> Generating Focus Analysis for:", target_shape_focus, "@", target_h_focus, "m...\n"))

# Data Preparation
# ------------------------------------------------------------------------------
# Filter the global dataset (df_combined)
df_focus_lai <- df_combined %>%
  filter(Height_Scenario == target_h_focus, Shape_Scenario == target_shape_focus) %>%
  mutate(
    # Create Month Label
    Month_Label = lubridate::month(time, label = TRUE, abbr = FALSE, locale = "C"),
    # Convert LAI to factor for discrete Y-axis plotting
    LAI_Factor = as.factor(LAI_Scenario)
  )

#  PLOT E1: RIDGELINE PLOT (Temperature Distributions)
# ------------------------------------------------------------------------------
# Shows how the temperature curve shifts and squeezes as LAI increases.

p_ridge_lai <- ggplot(df_focus_lai, aes(x = Tair_micro, y = LAI_Factor, fill = stat(x))) +
  
  # Create density ridges with median lines
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01, 
                               quantile_lines = TRUE, quantiles = 2) +
  
  # Color by Temperature
  scale_fill_viridis_c(name = "Temp. (°C)", option = "magma") +
  
  # Facet by Month
  facet_wrap(~Month_Label, ncol = 2) +
  
  labs(
    title = paste("Impact of LAI on Temperature Distribution (", target_shape_focus, ", ", target_h_focus, "m)", sep=""),
    subtitle = "Monthly evolution of thermal density. Black line indicates median.\nNote the shift to the left (cooling) as LAI increases (going up Y-axis).",
    x = "Microclimate Temperature (°C)",
    y = "Leaf Area Index (LAI)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_ridge_lai)
ggsave(file.path(plot_dir, paste0("E1_Focus_", target_shape_focus, "_Ridgeline.png")), 
       p_ridge_lai, width = 10, height = 8)

# ==============================================================================
# PLOT E1 (ALT): RIDGELINE PLOT (DAILY TMAX DISTRIBUTION)
# ==============================================================================
# Goal: Show how the distribution of daily PEAK heat shifts as LAI increases.

# 1. Prepare Daily Tmax Data
# ------------------------------------------------------------------------------
df_focus_daily_tmax <- df_focus_lai %>%
  mutate(date = as.Date(time)) %>%
  group_by(Month_Label, LAI_Factor, date) %>%
  summarise(
    Daily_Tmax = max(Tair_micro), # Extract the daily maximum
    .groups = "drop"
  )

# 2. Plotting
# ------------------------------------------------------------------------------
p_ridge_tmax <- ggplot(df_focus_daily_tmax, aes(x = Daily_Tmax, y = LAI_Factor, fill = stat(x))) +
  
  # Create density ridges
  geom_density_ridges_gradient(
    scale = 2.5, 
    rel_min_height = 0.01, 
    quantile_lines = TRUE, 
    quantiles = 2 # Shows the median
  ) +
  
  # Color by Temperature (Magma is good for heat)
  scale_fill_viridis_c(name = "Tmax (°C)", option = "magma") +
  
  # Facet by Month
  facet_wrap(~Month_Label, ncol = 2) +
  
  labs(
    title = paste("Impact of LAI on Daily Maximum Temperatures (", target_shape_focus, ", ", target_h_focus, "m)", sep=""),
    subtitle = "Distribution of Daily Tmax. Shifting left means reducing heatwave intensity.",
    x = "Daily Maximum Temperature (°C)",
    y = "Leaf Area Index (LAI)"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_ridge_tmax)
ggsave(file.path(plot_dir, paste0("E1b_Focus_", target_shape_focus, "_Ridgeline_Tmax.png")), p_ridge_tmax, width = 10, height = 8)


# PLOT E2: EFFICIENCY HEATMAP (Mean Delta T)
# ------------------------------------------------------------------------------
# Summary Matrix: How much cooling do we get per LAI unit per Month?

stats_matrix <- df_focus_lai %>%
  group_by(Month_Label, LAI_Factor) %>%
  summarise(
    Mean_Cooling = mean(Delta_T), # Delta T = Micro - Macro
    .groups = "drop"
  )

p_matrix_lai <- ggplot(stats_matrix, aes(x = Month_Label, y = LAI_Factor, fill = Mean_Cooling)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Mean_Cooling, 2)), color = "white", size = 4) +
  
  # Darker colors = Stronger Cooling (Negative values)
  scale_fill_viridis_c(option = "viridis", direction = -1, name = "Delta T (°C)") +
  
  labs(
    title = paste("Thermal Efficiency Matrix (", target_shape_focus, ")", sep=""),
    subtitle = "Mean Cooling (Micro - ERA5) by Month and LAI.\nDarker/Negative values indicate higher efficiency.",
    x = NULL,
    y = "Leaf Area Index (LAI)"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p_matrix_lai)
ggsave(file.path(plot_dir, paste0("E2_Focus_", target_shape_focus, "_Heatmap.png")), 
       p_matrix_lai, width = 8, height = 6)

cat("=== All Analyses Complete ===\n")















# ==============================================================================
# SECTION F: STATISTICAL MODELLING (Variable Importance)
# ==============================================================================
library(lme4)      # For Linear Mixed Models
library(lmerTest)  # For p-values in LMM
library(sjPlot)    # For plotting estimates (Optional but recommended)
library(broom.mixed) # For tidying model outputs

cat("--> Starting Statistical Analysis...\n")

# 1. DATA PREPARATION FOR MODELLING
# ------------------------------------------------------------------------------
# We focus on "Daily Maximum Cooling" (Delta Tmax) because it's the most
# relevant metric for heat mitigation. Hourly data is too noisy/autocorrelated.

df_stats <- df_combined %>%
  mutate(date = as.Date(time)) %>%
  group_by(Height_Scenario, 
           Shape_Scenario, 
           Clump_Scenario,
           LAI_Scenario, 
           date) %>%
  summarise(
    # Cooling = Micro - Macro. We want negative values (cooling)
    Delta_Tmax = max(Tair_micro) - max(Tair_macro),
    T_Macro_Max = max(Tair_macro), # Covariate: It cools more when it's hotter outside
    .groups = "drop"
  ) %>%
  # Standardize continuous inputs (Z-score) to compare coefficients
  mutate(
    LAI_Scaled = scale(LAI_Scenario),
    Height_Scaled = scale(Height_Scenario),
    Clump_Scaled = scale(Clump_Scenario),
    T_Macro_Scaled = scale(T_Macro_Max)
  )

# Set Reference Level for Shape (Uniform is usually the baseline)
df_stats$Shape_Scenario <- factor(df_stats$Shape_Scenario)
df_stats$Shape_Scenario <- relevel(df_stats$Shape_Scenario, ref = "uniform")

# 2. LINEAR MIXED MODEL (LMM)
# ------------------------------------------------------------------------------
# Formula: Cooling ~ Variables + (1|Date)
# We use (1|date) as a random intercept because all scenarios share the same weather
# on a given day. This isolates the structural effect from the weather effect.

cat("Fitting LMM Model...\n")

model_lmm <- lmer(Delta_Tmax ~ 
                    poly(LAI_Scaled, 2) + 
                    Height_Scaled + 
                    Clump_Scaled +
                    Shape_Scenario + 
                    T_Macro_Scaled + # Control for ambient temperature
                    (1 | date),      # Random effect for day-to-day weather variability
                  data = df_stats)

# 3. EXTRACT RESULTS & COEFFICIENTS
# ------------------------------------------------------------------------------
summary(model_lmm)

# Calculate ANOVA table (Type II Wald F tests with Kenward-Roger df)
# This gives the "Importance" (F-value) of each variable
anova_res <- anova(model_lmm)
print(anova_res)

# 4. VISUALIZATION OF EFFECTS (Forest Plot)
# ------------------------------------------------------------------------------
# We extract the fixed effects estimates
effects <- tidy(model_lmm, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = recode(term, 
                  "LAI_Scaled" = "LAI (1 SD)",
                  "Height_Scaled" = "Height (1 SD)",
                  "Clump_Scaled" = "Clumping (1 SD)",
                  "T_Macro_Scaled" = "T_Macro (1 SD)")
  )

p_stat_1 <- ggplot(effects, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, size = 1) +
  geom_point(size = 3, color = "darkblue") +
  labs(
    title = "Statistical Impact Drivers (LMM Coefficients)",
    subtitle = "Effect on Delta Tmax (Negative = More Cooling). Scaled inputs.",
    x = "Coefficient Estimate (Change in °C per SD)",
    y = "Variable"
  ) +
  theme_bw()

print(p_stat_1)
# ggsave(file.path(plot_dir, "F1_Statistical_Effects.png"), p_stat_1, width = 8, height = 6)


# 5. INTERACTION ANALYSIS (LAI x CLUMPING)
# ------------------------------------------------------------------------------
# Does Clumping matter more when LAI is high? (Physics says yes)

cat("Fitting Interaction Model...\n")
model_int <- lmer(Delta_Tmax ~ LAI_Scaled * Clump_Scaled + (1 | date), data = df_stats)

# Visualizing the interaction
# We create a prediction grid to visualize the result
grid_int <- expand.grid(
  LAI_Scenario = c(1, 2, 3, 4, 5, 6, 7, 8),
  Clump_Scenario = seq(0.5, 1, 0.1)
) %>%
  mutate(
    LAI_Scaled = (LAI_Scenario - mean(df_stats$LAI_Scenario)) / sd(df_stats$LAI_Scenario),
    Clump_Scaled = (Clump_Scenario - mean(df_stats$Clump_Scenario)) / sd(df_stats$Clump_Scenario)
  )

grid_int$Predicted_Delta <- predict(model_int, newdata = grid_int, re.form = NA)

p_stat_2 <- ggplot(grid_int, aes(x = Clump_Scenario, y = Predicted_Delta, color = as.factor(LAI_Scenario))) +
  geom_line(size = 1.2) +
  scale_color_viridis_d(option = "plasma", name = "LAI Level") +
  labs(
    title = "Interaction Effect: LAI x Clumping",
    subtitle = "Does aggregation (Clumping) hurt cooling efficiency more at high LAI?",
    x = "Clumping Factor (1 = Random, 0.5 = Clumped)",
    y = "Predicted Delta Tmax (°C)"
  ) +
  theme_bw()

print(p_stat_2)
# ggsave(file.path(plot_dir, "F2_Interaction_LAIxClump.png"), p_stat_2, width = 8, height = 6)

cat("Statistical Analysis Complete.\n")
































# ==============================================================================
# SECTION G: VISUALISATION SYNTHÉTIQUE (PHYSIQUE)
# ==============================================================================

# 1. Création d'un jeu de données résumé (Moyenne sur l'été)
# On s'intéresse au "Delta Tmax" moyen (Refroidissement maximal quotidien)
df_viz <- df_combined %>%
  mutate(date = as.Date(time)) %>%
  group_by(Height_Scenario, Shape_Scenario, Clump_Scenario, LAI_Scenario, date) %>%
  summarise(
    Daily_Delta = max(Tair_micro) - max(Tair_macro), # Delta Tmax du jour
    .groups = "drop"
  ) %>%
  # On moyenne sur toute la saison pour avoir une valeur par scénario
  group_by(Height_Scenario, Shape_Scenario, Clump_Scenario, LAI_Scenario) %>%
  summarise(
    Mean_Cooling = mean(Daily_Delta),
    SD_Cooling = sd(Daily_Delta),
    .groups = "drop"
  )

# Pour les plots, on ordonne les formes de la plus "performante" à la moins performante
# (Basé sur tes résultats précédents ou une logique visuelle)
shape_order <- c("top_peak_uniform", "top_heavy", "uniform", "symmetric", 
                 "top_bottom_heavy", "low_symmetric", "bottom_heavy", "extreme_bottom")

df_viz$Shape_Scenario <- factor(df_viz$Shape_Scenario, levels = shape_order)

# --- FIGURE 1: HIERARCHIE DES FORMES (Architecture vs LAI) ---
# On moyenne les hauteurs et clumping pour voir l'effet "pur" de la forme
p_viz1 <- ggplot(df_viz, aes(x = LAI_Scenario, y = Mean_Cooling, color = Shape_Scenario)) +
  # Ligne de tendance lissée
  geom_smooth(se = FALSE, size = 1.5, span = 0.5) + 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_viridis_d(option = "turbo", name = "Architecture (Profil)") +
  scale_x_continuous(breaks = 1:8) +
  
  labs(
    title = "Performance Thermique par Architecture",
    subtitle = "Comparaison de l'efficacité de refroidissement en fonction de la densité de feuillage (LAI).",
    x = "Leaf Area Index (LAI)",
    y = "Refroidissement Moyen (Micro - Macro) [°C]",
    caption = "Moyenne sur toutes les hauteurs et clumpings."
  ) +
  theme_bw() +
  theme(legend.position = "right")

print(p_viz1)
# ggsave(file.path(plot_dir, "G1_Synthese_Formes.png"), p_viz1, width = 10, height = 6)

# --- FIGURE 2: L'EFFET CLUMPING (Facet par Forme) ---
# On fixe une hauteur standard (ex: 20m) pour ne pas surcharger le graph
df_viz_20m <- df_viz %>% filter(Height_Scenario == 20)

p_viz2 <- ggplot(df_viz_20m, aes(x = LAI_Scenario, y = Mean_Cooling, 
                                 group = Clump_Scenario, color = as.factor(Clump_Scenario))) +
  
  geom_line(size = 1, alpha = 0.8) +
  
  # Facette par Forme pour voir qui est sensible au clumping
  facet_wrap(~Shape_Scenario, ncol = 4) +
  
  # Couleurs : Jaune (Clump=1, Continu) -> Violet (Clump=0.5, Troué)
  scale_color_viridis_d(option = "plasma", direction = -1, name = "Clumping Factor\n(1=Continu, 0.5=Troué)") +
  scale_x_continuous(breaks = c(2, 4, 6, 8)) +
  
  labs(
    title = "Sensibilité à la Continuité du Feuillage (Clumping)",
    subtitle = "Impact de l'agrégation sur le refroidissement pour un arbre de 20m.",
    x = "LAI",
    y = "Refroidissement Moyen (°C)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_viz2)
# ggsave(file.path(plot_dir, "G2_Synthese_Clumping.png"), p_viz2, width = 12, height = 7)

# --- FIGURE 3: HEATMAP LAI vs CLUMPING ---
# On agrège toutes les formes et hauteurs pour voir la tendance globale du modèle
df_heatmap <- df_viz %>%
  group_by(LAI_Scenario, Clump_Scenario) %>%
  summarise(Global_Cooling = mean(Mean_Cooling), .groups = "drop")

p_viz3 <- ggplot(df_heatmap, aes(x = as.factor(LAI_Scenario), y = as.factor(Clump_Scenario), fill = Global_Cooling)) +
  geom_tile(color = "white", size = 0.5) +
  
  # Ajout du texte pour lire les valeurs exactes
  geom_text(aes(label = round(Global_Cooling, 1)), color = "white", size = 3.5, fontface = "bold") +
  
  scale_fill_viridis_c(option = "magma", direction = 1, name = "Delta T (°C)") +
  
  labs(
    title = "Matrice de Performance Globale (LAI vs Clumping)",
    subtitle = "Moyenne de tous les profils et hauteurs. Plus foncé = Plus frais.",
    x = "Leaf Area Index (LAI)",
    y = "Clumping Factor"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p_viz3)
# ggsave(file.path(plot_dir, "G3_Heatmap_Global.png"), p_viz3, width = 8, height = 6)