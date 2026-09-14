rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}
# Analyse des données LIDAR
library(musica.tools)
library(rmusica)
library(ncdf4)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(cowplot)
library(lidR)
library(raster)
library(sf)
library(gstat)
library(terra)
library(lubridate)
library(lmom)
library(doParallel)
library(foreach)
library(purrr)
library(broom)
source("functions.R")

# Preparation
# date <- "27/06/2021"
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-23"), by = "day")
# "41_34" misses between 2021-09-24 and 2021-09-30

# LAD stack
LAI_stack <- terra::rast("in_files/ladstack_classic_no_na.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson")


ids_to_remove <- c("41_13", "41_14", "41_20", "41_34", 
                   "41_41", "41_50", "41_51", "41_53",
                   "41_17", "41_18", "41_19", "41_27",
                   "41_30", "41_39", "41_47", "41_49",
                   "41_55")
json <- json %>%
  filter(!id_plot %in% ids_to_remove) %>%
  arrange(id_plot)


df_coords <- json %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(coord_x_utm31n = X, coord_y_utm31n = Y)

height_levels <- seq(2.5, 2.5 + (nlyr(LAI_stack) - 1), by = 1)
lad_values <- terra::extract(LAI_stack, df_coords)
df_lad <- as.data.frame(lad_values)
df_lad <- cbind(df_lad, df_coords)
df_lad$ID <- json$id_plot
df_lad <- df_lad %>%
  filter(rowSums(!is.na(dplyr::select(., -ID, -coord_x_utm31n, -coord_y_utm31n))) > 0)

n_field <- df_lad$ID
df_lad$ID = NULL

# Observed
site <- "Blois"
csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
temperature_csv <- read.csv(csv_file)

# Convert the datetime column to POSIXct
temperature_df <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", 
                               tz = "UTC")) %>%
  arrange(as.numeric(sub("41_", "", id_plot)))

temperature_all_days <- data.frame()
for (date_study in date_seq) {
  temp_filtered <- temperature_df %>%
    filter(as.Date(datetime) == as.Date(date_study),
           position_sensor == "a")
  
  temperature_all_days <- bind_rows(temperature_all_days, temp_filtered)
}
available_plots <- unique(temperature_all_days$id_plot)

# MuSICA parameters
radius <- met_radius <- 25
# clumping <- met_clumping <- 0.9

# Metrics
# metrics_df <- read.csv("out_files/radius_test/metrics_results.csv", 
#                        header = TRUE)
metrics_csv <- read.csv("out_files/radius_test/metrics_results.csv", 
                        header = TRUE)
metrics_df <- metrics_csv %>%
  filter(id_plot %in% available_plots) %>%
  filter(radius == met_radius)
# filter(clumping == met_clumping)

# Simulated
df_Tair_all <- data.frame()
out.all <- list()

# date_study <- "2021-06-20"
# date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-06-30"), by = "day")

for (i in 1:nrow(df_lad)) {
  print(i)
  id_plot <- n_field[i]
  
  # Open the NetCDF file for this field
  nc_file <- paste0("out_files/radius_test/", radius,
  "m/musica_out_Blois_pt_", id_plot,
  "_radius_", radius, ".nc")
  # nc_file <- paste0("out_files/TS/lidar_constant",
  #                   "/musica_out_ERA5_Blois_2021_pt_", id_plot, ".nc")
  out <- nc_open(nc_file)
  out.all[[paste0("pt_", id_plot)]] <- out
  
  raw_data <- get_variable(out, "Tair_z")
  df_Tair <- raw_data %>%
    filter(as.Date(time) %in% date_seq, nair == 1) %>%
    mutate(
      Tair_z = Tair_z - 273.15,
      id_plot = id_plot
    )
  
  # for (date_study in date_seq) {
  #   # Extract the variable
  #   df_Tair <- get_variable(out, "Tair_z") %>%
  #     filter(as.Date(time) == date_study, nair == 1) %>%
  #     mutate(Tair_z = Tair_z - 273.15,
  #            id_plot = id_plot)
  
  # ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
  df_Tair_all <- bind_rows(df_Tair_all, df_Tair)
  # }
}
df_Tair_all <- df_Tair_all %>%
  mutate(time = floor_date(time, unit = "hour")) %>%
  arrange(time) %>%
  filter(id_plot %in% available_plots)

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
    rename(Tair = Température.instantanée...C.) %>%
    mutate(time = as.POSIXct(paste(Date, sprintf("%04d", Heure..TU.)), 
                             format = "%d/%m/%Y %H%M", tz = "UTC")) %>%
    dplyr::select(-Date, -Heure..TU., -Code.placette) %>%
    dplyr::select(time, Tair)
  
  ign_Tair_all <- bind_rows(ign_Tair_all, ign_Tair_filtered)
}

# ERA5
era5 <- nc_open("in_files/musica_in_Blois.nc")
era5_var <- get_variable(era5, "Tair")

era5_Tair_all <- data.frame()
for (date_study in date_seq) {
  date_study <- as.Date(date_study)
  
  era5_Tair_filtered <- era5_var %>%
    filter(format(time, "%Y-%m-%d") == format(date_study, "%Y-%m-%d")) %>%
    mutate(Tair = round(Tair - 273.15, 1)) %>%
    mutate(time = floor_date(time, unit = "hour")) # %>%
  # dplyr::select(-z)
  
  era5_Tair_all <- bind_rows(era5_Tair_all, era5_Tair_filtered)
}

# ------------------------------ Initial Plots ---------------------------------
# CHS41 vs ERA5
combined_data <- merge(era5_Tair_all, ign_Tair_all, by = "time", suffixes = c("_era5", "_ign"))
r_value <- cor(combined_data$Tair_era5, combined_data$Tair_ign)
bias <- mean(combined_data$Tair_ign - combined_data$Tair_era5)
nrmse <- sqrt(mean((combined_data$Tair_ign - combined_data$Tair_era5)^2)) / IQR(combined_data$Tair_ign)

# Fit a linear model
model <- lm(Tair_ign ~ Tair_era5, data = combined_data)
r2 <- summary(model)$r.squared
mae <- mean(abs(combined_data$Tair_ign - combined_data$Tair_era5))
slope <- coef(model)[2]
intercept <- coef(model)[1]

p <- ggplot(combined_data, aes(x = Tair_ign, y = Tair_era5)) +
  # geom_point(size = 2, alpha = 0.8) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "ERA5 (°C)",
       # title = "Air Temperature of ERA5 vs Local weather station",
       # subtitle = paste("R =", round(r_value, 2),
       #                  "| R² =", round(r2, 2),
       #                  "| NRMSE =", round(nrmse, 2),
       #                  "| Bias =", round(bias, 2),
       #                  "| Slope =", round(slope, 2),
       #                  "| Intercept =", round(intercept, 2))
  ) +
  annotate(
    "text",
    x = 1, y = 40,  # Top-left corner
    label = paste(
      "R² =", round(r2, 2),
      "\nNRMSE =", round(nrmse, 2),
      "\nMAE =", round(mae, 2), "°C"
    ),
    hjust = 0,  # Nudge left
    vjust = 1.5,    # Nudge down
    size = 5,       # Adjust text size
    color = "black",
    fontface = "bold"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme(
    legend.position = "right",
    # legend.text = element_text(size = 12),
    # axis.title = element_text(size = 14),
    # plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  theme_bw(base_size = 16) +
  coord_fixed(ratio = 1)
print(p)
ggsave("Plots/Silvilaser/air_era5.png", p, width = 6, height = 6, dpi = 300)

# Observed vs Simulated
# Merge observed and simulated data by time and id_plot
combined_obs_sim <- merge(df_Tair_all, temperature_all_days, by.x = c("time", "id_plot"), 
                          by.y = c("datetime", "id_plot"), suffixes = c("_sim", "_obs"))

# Calculate the correlation coefficient (R), Bias, and NRMSE
r_value_obs_sim <- cor(combined_obs_sim$Tair_z, combined_obs_sim$t_hobo)
bias_obs_sim <- mean(combined_obs_sim$t_hobo - combined_obs_sim$Tair_z)
nrmse_obs_sim <- sqrt(mean((combined_obs_sim$t_hobo - combined_obs_sim$Tair_z)^2)) / IQR(combined_obs_sim$t_hobo)

# Fit a linear model
model_obs_sim <- lm(t_hobo ~ Tair_z, data = combined_obs_sim)
r2_obs_sim <- summary(model_obs_sim)$r.squared
slope_obs_sim <- coef(model_obs_sim)[2]
intercept_obs_sim <- coef(model_obs_sim)[1]

# Plot observed vs simulated
ggplot(combined_obs_sim, aes(x = Tair_z, y = t_hobo)) +
  geom_point(aes(color = id_plot), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Simulated (°C)",
       y = "Observed (°C)",
       title = "Observed vs Simulated Tair",
       subtitle = paste("R =", round(r_value_obs_sim, 2),
                        "| R² =", round(r2_obs_sim, 2),
                        "| NRMSE =", round(nrmse_obs_sim, 2),
                        "| Bias =", round(bias_obs_sim, 2),
                        "| Slope =", round(slope_obs_sim, 2),
                        "| Intercept =", round(intercept_obs_sim, 2))) +
  xlim(c(min(combined_obs_sim$Tair_z, combined_obs_sim$t_hobo),
         max(combined_obs_sim$Tair_z, combined_obs_sim$t_hobo))) + 
  ylim(c(min(combined_obs_sim$Tair_z, combined_obs_sim$t_hobo),
         max(combined_obs_sim$Tair_z, combined_obs_sim$t_hobo))) + 
  theme_bw()


# Observed vs CHS 41
# Aggregate temperature_all_days to get mean t_hobo per hour
temperature_hourly <- temperature_all_days %>%
  group_by(time = floor_date(datetime, "hour")) %>%  # Round datetime to the hour
  summarise(t_hobo_mean = mean(t_hobo, na.rm = TRUE)) %>%
  ungroup()

# Merge with IGN data
combined_temp_ign <- merge(temperature_hourly, ign_Tair_all, by = "time")

# Compute statistical metrics
r_value_obs <- cor(combined_temp_ign$Tair, combined_temp_ign$t_hobo_mean)
bias_obs <- mean(combined_temp_ign$t_hobo_mean - combined_temp_ign$Tair)
nrmse_obs <- sqrt(mean((combined_temp_ign$t_hobo_mean - combined_temp_ign$Tair)^2)) / IQR(combined_temp_ign$Tair)

# Fit a linear model
model_obs <- lm(t_hobo_mean ~ Tair, data = combined_temp_ign)
r2_obs <- summary(model_obs)$r.squared
slope_obs <- coef(model_obs)[2]
intercept_obs <- coef(model_obs)[1]

# Plot observed (HOBO) vs CHS41 (IGN)
ggplot(combined_temp_ign, aes(x = Tair, y = t_hobo_mean)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "HOBO Mean Temperature (°C)",
       title = "Air Temperature of HOBO vs Local weather station",
       subtitle = paste("R =", round(r_value_obs, 2),
                        "| R² =", round(r2_obs, 2),
                        "| NRMSE =", round(nrmse_obs, 2),
                        "| Bias =", round(bias_obs, 2),
                        "| Slope =", round(slope_obs, 2),
                        "| Intercept =", round(intercept_obs, 2))) +
  xlim(c(min(combined_temp_ign$Tair, combined_temp_ign$t_hobo_mean),
         max(combined_temp_ign$Tair, combined_temp_ign$t_hobo_mean))) + 
  ylim(c(min(combined_temp_ign$Tair, combined_temp_ign$t_hobo_mean),
         max(combined_temp_ign$Tair, combined_temp_ign$t_hobo_mean))) + 
  theme_bw()


# Simulated vs ERA5
df_Tair_hourly <- df_Tair_all %>%
  group_by(time = floor_date(time, "hour")) %>%  # Round time to the hour
  summarise(Tair_sim_mean = mean(Tair_z, na.rm = TRUE)) %>%
  ungroup()

# Merge with ERA5 data
combined_sim <- merge(df_Tair_hourly, era5_Tair_all, by = "time")

# Compute statistical metrics
r_value_sim <- cor(combined_sim$Tair_sim_mean, combined_sim$Tair)
bias_sim <- mean(combined_sim$Tair - combined_sim$Tair_sim_mean)
nrmse_sim <- sqrt(mean((combined_sim$Tair - combined_sim$Tair_sim_mean)^2)) / IQR(combined_sim$Tair)

# Fit a linear model
model_sim <- lm(Tair ~ Tair_sim_mean, data = combined_sim)
r2_sim <- summary(model_sim)$r.squared
slope_sim <- coef(model_sim)[2]
intercept_sim <- coef(model_sim)[1]

# Plot simulated vs ERA5
ggplot(combined_sim, aes(x = Tair_sim_mean, y = Tair)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Simulated (°C)",
       y = "ERA5 (°C)",
       title = "Simulated vs ERA5 Tair",
       subtitle = paste("R =", round(r_value_sim, 2),
                        "| R² =", round(r2_sim, 2),
                        "| NRMSE =", round(nrmse_sim, 2),
                        "| Bias =", round(bias_sim, 2),
                        "| Slope =", round(slope_sim, 2),
                        "| Intercept =", round(intercept_sim, 2))) +
  xlim(c(min(combined_sim$Tair_sim_mean, combined_sim$Tair),
         max(combined_sim$Tair_sim_mean, combined_sim$Tair))) + 
  ylim(c(min(combined_sim$Tair_sim_mean, combined_sim$Tair),
         max(combined_sim$Tair_sim_mean, combined_sim$Tair))) + 
  theme_bw()
# ----------------------- log(slope) Gril et al. 2023 --------------------------
# Aggregate df_Tair_all to get mean Tair_sim per hour
# df_Tair_hourly <- df_Tair_all %>%
#   group_by(time = floor_date(time, "hour")) %>%  # Round time to the hour
#   summarise(Tair_sim_mean = mean(Tair_z, na.rm = TRUE)) %>%
#   ungroup()

# ----------------------------- ERA5 vs CHS 41 ---------------------------------
combined_clim <- merge(era5_Tair_all, ign_Tair_all, by = "time", suffixes = c("_era5", "_ign"))

# Compute statistical metrics
r_value <- cor(combined_clim$Tair_era5, combined_clim$Tair_ign)
bias <- mean(combined_clim$Tair_era5 - combined_clim$Tair_ign)
nrmse <- sqrt(mean((combined_clim$Tair_ign - combined_clim$Tair_era5)^2)) / IQR(combined_clim$Tair_ign)

# Fit a linear model
model <- lm(Tair_era5 ~ Tair_ign, data = combined_clim)
r2 <- summary(model)$r.squared
slope <- coef(model)[2]
intercept <- coef(model)[1]

ggplot(combined_clim, aes(x = Tair_ign, y = Tair_era5)) +
  geom_point(size = 2, alpha = 0.8, color = "#128d84") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "ERA5 (°C)",
       title = "Air Temperature of ERA5 vs Local weather station from June to September 2021",
       subtitle = paste("R =", round(r_value, 2),
                        "| R² =", round(r2, 2),
                        "| NRMSE =", round(nrmse, 2),
                        "| Bias =", round(bias, 2),
                        "| Slope =", round(slope, 2),
                        "| Intercept =", round(intercept, 2))) +
  xlim(c(min(combined_data$Tair_ign, combined_data$Tair_era5),
         max(combined_data$Tair_ign, combined_data$Tair_era5))) + 
  ylim(c(min(combined_data$Tair_ign, combined_data$Tair_era5),
         max(combined_data$Tair_ign, combined_data$Tair_era5))) + 
  theme_bw(base_size = 16)

# ---------------------------- Simulated vs ERA5 -------------------------------
df_Tair_hourly <- df_Tair_all %>%
  mutate(time = floor_date(time, "hour"))

# Merge with ERA5 data
combined_sim1 <- merge(df_Tair_hourly, era5_Tair_all, by = "time")
combined_sim <- merge(combined_sim1, metrics_df, by = "id_plot")

# Compute statistical metrics
r_value_sim <- cor(combined_sim$Tair_z, combined_sim$Tair)
bias_sim <- mean(combined_sim$Tair_z - combined_sim$Tair)
nrmse_sim <- sqrt(mean((combined_sim$Tair - combined_sim$Tair_z)^2)) / IQR(combined_sim$Tair)
mae_sim <- mean(abs(combined_sim$Tair_z - combined_sim$Tair))
r_squared_sim <- summary(lm(Tair_z ~ Tair, 
                            data = combined_sim))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_sim <- combined_sim %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(Tair_z ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(Tair_z ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = Tair_z - Tair
  ) %>%
  ungroup()

# Split the data into two groups based on slope
combined_sim_neg <- combined_sim %>% filter(log_slope < 0)
r_value_sim_neg <- cor(combined_sim_neg$Tair_z, combined_sim_neg$Tair)
bias_sim_neg <- mean(combined_sim_neg$Tair_z - combined_sim_neg$Tair)
nrmse_sim_neg <- sqrt(mean((combined_sim_neg$Tair - combined_sim_neg$Tair_z)^2)) / IQR(combined_sim_neg$Tair)
mae_sim_neg <- mean(abs(combined_sim_neg$Tair_z - combined_sim_neg$Tair))
r_squared_sim_neg <- summary(lm(Tair_z ~ Tair, 
                                data = combined_sim_neg))$r.squared

combined_sim_pos <- combined_sim %>% filter(log_slope > 0)
r_value_sim_pos <- cor(combined_sim_pos$Tair_z, combined_sim_pos$Tair)
bias_sim_pos <- mean(combined_sim_pos$Tair_z - combined_sim_pos$Tair)
nrmse_sim_pos <- sqrt(mean((combined_sim_pos$Tair - combined_sim_pos$Tair_z)^2)) / IQR(combined_sim_pos$Tair)
mae_sim_pos <- mean(abs(combined_sim_pos$Tair_z - combined_sim_pos$Tair))
r_squared_sim_pos <- summary(lm(Tair_z ~ Tair, 
                                data = combined_sim_pos))$r.squared

# Plot for slope < 0
p1 <- ggplot(combined_sim_neg, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly macroclimate measures (°C)",
       # x = "ERA5 (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "(a)",
       # title = "Predicted by MuSICA vs ERA5 (Slope < 0)",
       # subtitle = paste("R =", round(r_value_sim_neg, 2),
       #                  "| NRMSE =", round(nrmse_sim_neg, 2),
       #                  "| MAE =", round(mae_sim_neg, 2)
       # ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = 1, y = 40,  # Top-left corner
    label = paste(
      "R² =", round(r_squared_sim_neg, 2),
      "\nNRMSE =", round(nrmse_sim_neg, 2),
      "\nMAE =", round(mae_sim_neg, 2)
    ),
    hjust = 0,  # Nudge left
    vjust = 1.5,    # Nudge down
    size = 4,       # Adjust text size
    color = "black",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 39, y = 0,  # Bottom-right corner
    label = paste("N =", n_distinct(combined_sim_neg$id_plot), "plots",
                  "\nPAI =", round(mean(combined_sim_neg$pai), 2),
                  "\nVCI =", round(mean(combined_sim_neg$vci), 2),
                  "\nHmax =", round(mean(combined_sim_neg$max), 2)
    ),
    hjust = 1,  # Nudge right
    vjust = -0.5, # Nudge up
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#d9e2d1"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed(ratio = 1)

# Plot for slope > 0
p2 <- ggplot(combined_sim_pos, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly macroclimate measures (°C)",
       # x = "ERA5 (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "(b)",
       # title = "Predicted by MuSICA vs ERA5 (Slope > 0)",
       # subtitle = paste("R =", round(r_value_sim_pos, 2),
       #                  "| NRMSE =", round(nrmse_sim_pos, 2),
       #                  "| MAE =", round(mae_sim_neg, 2)
       # ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = 1, y = 40,  # Top-left corner
    label = paste(
      "R² =", round(r_squared_sim_pos, 2),
      "\nNRMSE =", round(nrmse_sim_pos, 2),
      "\nMAE =", round(mae_sim_pos, 2)
    ),
    hjust = 0,  # Nudge left
    vjust = 1.5,    # Nudge down
    size = 4,       # Adjust text size
    color = "black",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 39, y = 0,  # Bottom-right corner
    label = paste("N =", n_distinct(combined_sim_pos$id_plot), "plots",
                  "\nPAI =", round(mean(combined_sim_pos$pai), 2),
                  "\nVCI =", round(mean(combined_sim_pos$vci), 2),
                  "\nHmax =", round(mean(combined_sim_pos$max), 2)
    ),
    hjust = 1,  # Nudge right
    vjust = -0.5, # Nudge up
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  scale_color_gradient2(
    low = "#f6efcb", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed(ratio = 1)

grid.arrange(p1, p2, ncol = 2)

# Final plot
p12 <- ggplot(combined_sim, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  
  # geom_smooth(method = "lm", se = FALSE, color = "black") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "ERA5 (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "Predicted by MuSICA vs ERA5",
       subtitle = paste("R =", round(r_value_sim, 2),
                        "| NRMSE =", round(nrmse_sim, 2),
                        "| Bias =", round(bias_sim, 2)
       ),
       color = "Log Slope"
  ) +
  # annotate(
  #   "text",
  #   x = Inf, y = -Inf,
  #   label = paste("N =", n_distinct(combined_sim$id_plot), "plots"),
  #   hjust = 1.1, vjust = -1.1,
  #   color = "black",
  #   size = 5
  # ) +
  scale_color_gradient2(
    low = "#128d84", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)
print(p12)

# Equilibrium
equilibrium_sim <- get_equilibrium_stats(combined_sim)
print(equilibrium_sim)

print(mean(combined_sim_pos$vci))
print(mean(combined_sim_neg$vci))

print(mean(combined_sim_pos$lcv))
print(mean(combined_sim_neg$lcv))

print(mean(combined_sim_pos$mean))
print(mean(combined_sim_neg$mean))

print(mean(combined_sim_pos$max))
print(mean(combined_sim_neg$max))

print(mean(combined_sim_pos$pai))
print(mean(combined_sim_neg$pai))

print(mean(combined_sim_pos$gf))
print(mean(combined_sim_neg$gf))

result <- combined_sim_neg %>%
  group_by(vci) %>%
  filter(n_distinct(id_plot) > 1) %>%
  dplyr::select(id_plot, vci) %>%
  distinct() %>%
  slice(1:2)
print(result, n = 20)

df_filtered <- combined_sim_neg %>%
  filter(id_plot %in% c("41_12", "41_23"))
mean_offset <- df_filtered %>%
  group_by(id_plot) %>%
  summarise(mean_offset = mean(offset, na.rm = TRUE))
print(mean_offset)

# --------------------------- Simulated vs CHS 41 ------------------------------
df_Tair_hourly <- df_Tair_all %>%
  mutate(time = floor_date(time, "hour"))

# Merge with ERA5 data
combined_sim_chs <- merge(df_Tair_hourly, ign_Tair_all, by = "time")
combined_sim_chs41 <- merge(combined_sim_chs, metrics_df, by = "id_plot")

# Compute statistical metrics
r_value_sim_chs41 <- cor(combined_sim_chs41$Tair_z, combined_sim_chs41$Tair)
bias_sim_chs41 <- mean(combined_sim_chs41$Tair_z - combined_sim_chs41$Tair)
nrmse_sim_chs41 <- sqrt(mean((combined_sim_chs41$Tair - combined_sim_chs41$Tair_z)^2)) / IQR(combined_sim_chs41$Tair)

# Fit a linear model and calculate log(slope) for each point
combined_sim_chs41 <- combined_sim_chs41 %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(Tair_z ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(Tair_z ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = Tair_z - Tair
  ) %>%
  ungroup()

# Split the data into two groups based on slope
combined_sim_chs41_neg <- combined_sim_chs41 %>% filter(log_slope < 0)
r_value_sim_chs41_neg <- cor(combined_sim_chs41_neg$Tair_z, combined_sim_chs41_neg$Tair)
bias_sim_chs41_neg <- mean(combined_sim_chs41_neg$Tair_z - combined_sim_chs41_neg$Tair)
nrmse_sim_chs41_neg <- sqrt(mean((combined_sim_chs41_neg$Tair - combined_sim_chs41_neg$Tair_z)^2)) / IQR(combined_sim_chs41_neg$Tair)

combined_sim_chs41_pos <- combined_sim_chs41 %>% filter(log_slope > 0)
r_value_sim_chs41_pos <- cor(combined_sim_chs41_pos$Tair_z, combined_sim_chs41_pos$Tair)
bias_sim_chs41_pos <- mean(combined_sim_chs41_pos$Tair_z - combined_sim_chs41_pos$Tair)
nrmse_sim_chs41_pos <- sqrt(mean((combined_sim_chs41_pos$Tair - combined_sim_chs41_pos$Tair_z)^2)) / IQR(combined_sim_chs41_pos$Tair)

# Plot for slope < 0
plot_neg <- ggplot(combined_sim_chs41_neg, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "Predicted hourly by MuSICA (°C)",
       title = "Predicted by MuSICA vs Local weather station (Slope < 0)",
       subtitle = paste("R =", round(r_value_sim_chs41_neg, 2),
                        "| NRMSE =", round(nrmse_sim_chs41_neg, 2),
                        "| Bias =", round(bias_sim_chs41_neg, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_sim_chs41_neg$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#d9e2d1"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Plot for slope > 0
plot_pos <- ggplot(combined_sim_chs41_pos, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "CHS 41 (°C)",
       y = "Predicted hourly by MuSICA (°C)",
       title = "Predicted by MuSICA vs CHS 41 (Slope > 0)",
       subtitle = paste("R =", round(r_value_sim_chs41_pos, 2),
                        "| NRMSE =", round(nrmse_sim_chs41_pos, 2),
                        "| Bias =", round(bias_sim_chs41_pos, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_sim_chs41_pos$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#f6efcb", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

grid.arrange(plot_neg, plot_pos, ncol = 2)

# Final plot
ggplot(combined_sim_chs41, aes(x = Tair, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  
  # geom_smooth(method = "lm", se = FALSE, color = "black") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "CHS 41 (°C)",
       y = "Predicted hourly by MuSICA (°C)",
       title = "Predicted by MuSICA vs CHS 41",
       subtitle = paste("R =", round(r_value_sim_chs41, 2),
                        "| NRMSE =", round(nrmse_sim_chs41, 2),
                        "| Bias =", round(bias_sim_chs41, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_sim_chs41$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Equilibrium
equilibrium_sim_chs41 <- get_equilibrium_stats(combined_sim_chs41)
print(equilibrium_sim_chs41)

print(mean(combined_sim_chs41_pos$vci))
print(mean(combined_sim_chs41_neg$vci))

print(mean(combined_sim_chs41_pos$lcv))
print(mean(combined_sim_chs41_neg$lcv))

print(mean(combined_sim_chs41_pos$mean))
print(mean(combined_sim_chs41_neg$mean))

print(mean(combined_sim_chs41_pos$max))
print(mean(combined_sim_chs41_neg$max))

print(mean(combined_sim_chs41_pos$pai))
print(mean(combined_sim_chs41_neg$pai))

print(mean(combined_sim_chs41_pos$gf))
print(mean(combined_sim_chs41_neg$gf))

# --------------------------- Observed vs CHS 41 -------------------------------
t_hobo_hourly <- temperature_all_days %>%
  mutate(time = floor_date(datetime, "hour"))
combined_obs1 <- merge(t_hobo_hourly, ign_Tair_all, by = "time")
combined_obs <- merge(combined_obs1, metrics_df, by = "id_plot")

# Compute statistical metrics
r_value_obs <- cor(combined_obs$t_hobo, combined_obs$Tair)
bias_obs <- mean(combined_obs$t_hobo - combined_obs$Tair)
nrmse_obs <- sqrt(mean((combined_obs$Tair - combined_obs$t_hobo)^2)) / IQR(combined_obs$Tair)

# Fit a linear model and calculate log(slope) for each point
combined_obs <- combined_obs %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(t_hobo ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(t_hobo ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = t_hobo - Tair
  ) %>%
  ungroup()

# Split the data into two groups based on slope
combined_obs_neg <- combined_obs %>% filter(log_slope < 0)
r_value_obs_neg <- cor(combined_obs_neg$t_hobo, combined_obs_neg$Tair)
bias_obs_neg <- mean(combined_obs_neg$t_hobo - combined_obs_neg$Tair)
nrmse_obs_neg <- sqrt(mean((combined_obs_neg$Tair - combined_obs_neg$t_hobo)^2)) / IQR(combined_obs_neg$Tair)

combined_obs_pos <- combined_obs %>% filter(log_slope > 0)
r_value_obs_pos <- cor(combined_obs_pos$t_hobo, combined_obs_pos$Tair)
bias_obs_pos <- mean(combined_obs_pos$t_hobo - combined_obs_pos$Tair)
nrmse_obs_pos <- sqrt(mean((combined_obs_pos$Tair - combined_obs_pos$t_hobo)^2)) / IQR(combined_obs_pos$Tair)

# Plot for slope < 0
plot_neg <- ggplot(combined_obs_neg, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "CHS 41 (°C)",
       y = "Hourly HOBO Sensors Measures (°C)",
       title = "HOBO Sensors Measures vs CHS 41 (Slope < 0)",
       subtitle = paste("R =", round(r_value_obs_neg, 2),
                        "| NRMSE =", round(nrmse_obs_neg, 2),
                        "| Bias =", round(bias_obs_neg, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_neg$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#d9e2d1"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Plot for slope > 0
plot_pos <- ggplot(combined_obs_pos, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "CHS 41 (°C)",
       y = "Hourly HOBO Sensors Measures (°C)",
       title = "HOBO Sensors Measures vs CHS 41 (Slope > 0)",
       subtitle = paste("R =", round(r_value_obs_pos, 2),
                        "| NRMSE =", round(nrmse_obs_pos, 2),
                        "| Bias =", round(bias_obs_pos, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_pos$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#f6efcb", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

grid.arrange(plot_neg, plot_pos, ncol = 2)

# Final plot
ggplot(combined_obs, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "CHS 41 (°C)",
       y = "Hourly HOBO Sensors Measures (°C)",
       title = "HOBO Sensors Measures vs CHS 41",
       subtitle = paste("R =", round(r_value_obs, 2),
                        "| NRMSE =", round(nrmse_obs, 2),
                        "| Bias =", round(bias_obs, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Equilibrium
equilibrium_obs <- get_equilibrium_stats(combined_obs)
print(equilibrium_obs)

print(mean(combined_obs_pos$vci))
print(mean(combined_obs_neg$vci))

print(mean(combined_obs_pos$lcv))
print(mean(combined_obs_neg$lcv))

print(mean(combined_obs_pos$mean))
print(mean(combined_obs_neg$mean))

print(mean(combined_obs_pos$max))
print(mean(combined_obs_neg$max))

print(mean(combined_obs_pos$pai))
print(mean(combined_obs_neg$pai))

print(mean(combined_obs_pos$gf))
print(mean(combined_obs_neg$gf))

# ---------------------------- Observed vs ERA5 --------------------------------
combined_obs_era1 <- merge(t_hobo_hourly, era5_Tair_all, by = "time")
combined_obs_era5 <- merge(combined_obs_era1, metrics_df, by = "id_plot")

# Compute statistical metrics
r_value_obs_era5 <- cor(combined_obs_era5$t_hobo, combined_obs_era5$Tair)
bias_obs_era5 <- mean(combined_obs_era5$t_hobo - combined_obs_era5$Tair)
nrmse_obs_era5 <- sqrt(mean((combined_obs_era5$Tair - combined_obs_era5$t_hobo)^2)) / IQR(combined_obs_era5$Tair)
mae_obs_era5 <- mean(abs(combined_obs_era5$Tair - combined_obs_era5$t_hobo))
r_squared_obs_era5 <- summary(lm(t_hobo ~ Tair, 
                                 data = combined_obs_era5))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_obs_era5 <- combined_obs_era5 %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(t_hobo ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(t_hobo ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = t_hobo - Tair
  ) %>%
  ungroup()

# Split the data into two groups based on slope
combined_obs_era5_neg <- combined_obs_era5 %>% filter(log_slope < 0)
r_value_obs_era5_neg <- cor(combined_obs_era5_neg$t_hobo, combined_obs_era5_neg$Tair)
bias_obs_era5_neg <- mean(combined_obs_era5_neg$t_hobo - combined_obs_era5_neg$Tair)
nrmse_obs_era5_neg <- sqrt(mean((combined_obs_era5_neg$Tair - combined_obs_era5_neg$t_hobo)^2)) / IQR(combined_obs_era5_neg$Tair)
mae_obs_era5_neg <- mean(abs(combined_obs_era5_neg$Tair - combined_obs_era5_neg$t_hobo))
r_squared_obs_era5_neg <- summary(lm(t_hobo ~ Tair, 
                                     data = combined_obs_era5_neg))$r.squared

combined_obs_era5_pos <- combined_obs_era5 %>% filter(log_slope > 0)
r_value_obs_era5_pos <- cor(combined_obs_era5_pos$t_hobo, combined_obs_era5_pos$Tair)
bias_obs_era5_pos <- mean(combined_obs_era5_pos$t_hobo - combined_obs_era5_pos$Tair)
nrmse_obs_era5_pos <- sqrt(mean((combined_obs_era5_pos$Tair - combined_obs_era5_pos$t_hobo)^2)) / IQR(combined_obs_era5_pos$Tair)
mae_obs_era5_pos <- mean(abs(combined_obs_era5_pos$Tair - combined_obs_era5_pos$t_hobo))
r_squared_obs_era5_pos <- summary(lm(t_hobo ~ Tair, 
                                     data = combined_obs_era5_pos))$r.squared

# Plot for slope < 0
plot_neg <- ggplot(combined_obs_era5_neg, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "ERA5 (°C)",
       y = "Hourly HOBO Sensors Measures (°C)",
       title = "HOBO Sensors Measures vs ERA5 (Slope < 0)",
       subtitle = paste("R =", round(r_value_obs_era5_neg, 2),
                        "| NRMSE =", round(nrmse_obs_era5_neg, 2),
                        "| MAE =", round(mae_obs_era5_neg, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_era5_neg$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#d9e2d1"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Plot for slope > 0
plot_pos <- ggplot(combined_obs_era5_pos, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "ERA5 (°C)",
       y = "Hourly HOBO Sensors Measures (°C)",
       title = "HOBO Sensors Measures vs ERA5 (Slope > 0)",
       subtitle = paste("R =", round(r_value_obs_era5_pos, 2),
                        "| NRMSE =", round(nrmse_obs_era5_pos, 2),
                        "| MAE =", round(mae_obs_era5_pos, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_era5_pos$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#f6efcb", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

grid.arrange(plot_neg, plot_pos, ncol = 2)

# Final plot
p3 <- ggplot(combined_obs_era5, aes(x = Tair, y = t_hobo, color = log_slope)) +
  geom_point(size = 2, alpha = 0.5) +
  
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly macroclimate measures (°C)",
       # x = "ERA5 (°C)",
       y = "Hourly HOBO sensors measures (°C)",
       title = "(c)",
       # title = "HOBO Sensors Measures vs ERA5",
       # subtitle = paste("R =", round(r_value_obs_era5, 2),
       #                  "| NRMSE =", round(nrmse_obs_era5, 2),
       #                  "| MAE =", round(mae_obs_era5, 2)
       # ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = 1, y = 40,  # Top-left corner
    label = paste(
      "R² =", round(r_squared_obs_era5, 2),
      "\nNRMSE =", round(nrmse_obs_era5, 2),
      "\nMAE =", round(mae_obs_era5, 2)
    ),
    hjust = 0,  # Nudge left
    vjust = 1.5,    # Nudge down
    size = 4,       # Adjust text size
    color = "black",
    fontface = "bold"
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "right",
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed(ratio = 1)
print(p3)

# Equilibrium
equilibrium_obs_era5 <- get_equilibrium_stats(combined_obs_era5)
print(equilibrium_obs_era5)

print(mean(combined_obs_era5_pos$vci))
print(mean(combined_obs_era5_neg$vci))

print(mean(combined_obs_era5_pos$lcv))
print(mean(combined_obs_era5_neg$lcv))

print(mean(combined_obs_era5_pos$mean))
print(mean(combined_obs_era5_neg$mean))

print(mean(combined_obs_era5_pos$max))
print(mean(combined_obs_era5_neg$max))

print(mean(combined_obs_era5_pos$pai))
print(mean(combined_obs_era5_neg$pai))

print(mean(combined_obs_era5_pos$gf))
print(mean(combined_obs_era5_neg$gf))

# -------------------------- Simulated vs Observed -----------------------------
df_Tair_hourly <- df_Tair_hourly %>%
  dplyr::select(id_plot, time, Tair_z)
t_hobo_hourly <- t_hobo_hourly %>%
  dplyr::select(id_plot, time, t_hobo)

combined_obs_sim1 <- merge(df_Tair_hourly, t_hobo_hourly, by = c("id_plot", "time"))
combined_obs_sim <- merge(combined_obs_sim1, metrics_df, by = "id_plot") %>%
  mutate(hour = hour(time))

# combined_obs_sim <- combined_obs_sim %>%
#   filter(hour %in% c(7:21))

r_value_obs_sim <- cor(combined_obs_sim$t_hobo, combined_obs_sim$Tair_z)
bias_obs_sim <- mean(combined_obs_sim$Tair_z - combined_obs_sim$t_hobo)
nrmse_obs_sim <- sqrt(mean((combined_obs_sim$t_hobo - combined_obs_sim$Tair_z)^2)) / IQR(combined_obs_sim$t_hobo)
mae_obs_sim <- mean(abs(combined_obs_sim$Tair_z - combined_obs_sim$t_hobo))
r_squared_obs_sim <- summary(lm(Tair_z ~ t_hobo, 
                                data = combined_obs_sim))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_obs_sim <- combined_obs_sim %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(Tair_z ~ t_hobo, data = cur_data()))[2],
    equilibrium = coef(lm(Tair_z ~ t_hobo, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = Tair_z - t_hobo
  ) %>%
  ungroup()

# combined_obs_sim <- combined_obs_sim %>%
#   filter(hour %in% c(0, 1, 2, 3, 4, 5, 6, 22, 23))

# Split the data into two groups based on slope
combined_obs_sim_neg <- combined_obs_sim %>% filter(log_slope < 0)
r_value_obs_sim_neg <- cor(combined_obs_sim_neg$t_hobo, combined_obs_sim_neg$Tair_z)
bias_obs_sim_neg <- mean(combined_obs_sim_neg$Tair_z - combined_obs_sim_neg$t_hobo)
nrmse_obs_sim_neg <- sqrt(mean((combined_obs_sim_neg$t_hobo - combined_obs_sim_neg$Tair_z)^2)) / IQR(combined_obs_sim_neg$t_hobo)
mae_obs_sim_neg <- mean(abs(combined_obs_sim_neg$Tair_z - combined_obs_sim_neg$t_hobo))
r_squared_obs_sim_neg <- summary(lm(Tair_z ~ t_hobo, 
                                    data = combined_obs_sim_neg))$r.squared

combined_obs_sim_pos <- combined_obs_sim %>% filter(log_slope > 0)
r_value_obs_sim_pos <- cor(combined_obs_sim_pos$t_hobo, combined_obs_sim_pos$Tair_z)
bias_obs_sim_pos <- mean(combined_obs_sim_pos$Tair_z - combined_obs_sim_pos$t_hobo)
nrmse_obs_sim_pos <- sqrt(mean((combined_obs_sim_pos$t_hobo - combined_obs_sim_pos$Tair_z)^2)) / IQR(combined_obs_sim_pos$t_hobo)
mae_obs_sim_pos <- mean(abs(combined_obs_sim_pos$Tair_z - combined_obs_sim_pos$t_hobo))
r_squared_obs_sim_pos <- summary(lm(Tair_z ~ t_hobo, 
                                    data = combined_obs_sim_pos))$r.squared

# Create the plot for slope < 0
plot_neg <- ggplot(combined_obs_sim_neg, aes(x = t_hobo, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly HOBO sensors measures (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "MuSICA Predictions vs HOBO Sensors (Slope < 0)",
       subtitle = paste("R =", round(r_value_obs_sim_neg, 2),
                        "| NRMSE =", round(nrmse_obs_sim_neg, 2),
                        "| Bias =", round(bias_obs_sim_neg, 2),
                        "| MAE =", round(mae_obs_sim_neg, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_sim_neg$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#d9e2d1"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Create the plot for slope > 0
plot_pos <- ggplot(combined_obs_sim_pos, aes(x = t_hobo, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly HOBO sensors measures (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "MuSICA Predictions vs HOBO Sensors (Slope > 0)",
       subtitle = paste("R =", round(r_value_obs_sim_pos, 2),
                        "| NRMSE =", round(nrmse_obs_sim_pos, 2),
                        "| Bias =", round(bias_obs_sim_pos, 2),
                        "| MAE =", round(mae_obs_sim_pos, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_sim_pos$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#f6efcb", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

grid.arrange(plot_neg, plot_pos, ncol = 2)

# Final plot
ggplot(combined_obs_sim, aes(x = t_hobo, y = Tair_z, color = log_slope)) +
  geom_point(size = 2, alpha = 0.8) +
  
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly HOBO sensors measures (°C)",
       y = "Hourly MuSICA predictions (°C)",
       title = "MuSICA Predictions vs HOBO Sensors",
       subtitle = paste("R =", round(r_value_obs_sim, 2),
                        "| NRMSE =", round(nrmse_obs_sim, 2),
                        "| MAE =", round(mae_obs_sim, 2)
       ),
       color = "Log Slope"
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = paste("N =", n_distinct(combined_obs_sim$id_plot), "plots"),
    hjust = 1.1, vjust = -1.1,
    color = "black",
    size = 5
  ) +
  scale_color_gradient2(
    low = "#128d84", high = "#e9c624"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm")
  ) +
  coord_fixed(ratio = 1)

# Equilibrium
equilibrium_obs_sim <- get_equilibrium_stats(combined_obs_sim)
print(equilibrium_obs_sim)

print(mean(combined_obs_sim_pos$vci))
print(mean(combined_obs_sim_neg$vci))

print(mean(combined_obs_sim_pos$lcv))
print(mean(combined_obs_sim_neg$lcv))

print(mean(combined_obs_sim_pos$mean))
print(mean(combined_obs_sim_neg$mean))

print(mean(combined_obs_sim_pos$max))
print(mean(combined_obs_sim_neg$max))

print(mean(combined_obs_sim_pos$pai))
print(mean(combined_obs_sim_neg$pai))

print(mean(combined_obs_sim_pos$gf))
print(mean(combined_obs_sim_neg$gf))

# hist(combined_obs_sim_pos$lcv)
# hist(combined_obs_sim_neg$lcv)
# plot(combined_obs_sim$mean, combined_obs_sim$lcv)
# 
# print(mean(combined_obs_sim_pos$offset))
# print(mean(combined_obs_sim_neg$offset))
# 
# mean_offset_per_plot <- combined_obs_sim %>%
#   group_by(id_plot) %>%
#   summarise(mean_offset = mean(offset, na.rm = TRUE))
# 
# print(mean_offset_per_plot)

p4 <- ggplot(combined_obs_sim, aes(x = t_hobo, y = Tair_z)) +
  # geom_point(size = 2, alpha = 0.8, color = "#128d84") +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Hourly field sensors measures (°C)",
       y = "Hourly MuSICA predictions (°C)"
       # title = "(c)"
       # title = "MuSICA Predictions vs HOBO Sensors",
       # subtitle = paste("R =", round(r_value_obs_sim, 2),
       #                  "| NRMSE =", round(nrmse_obs_sim, 2),
       #                  "| MAE =", round(mae_obs_sim, 2))
  ) +
  annotate(
    "text",
    x = 1, y = 40,  # Top-left corner
    label = paste(
      "R² =", round(r_squared_obs_sim, 2),
      "\nNRMSE =", round(nrmse_obs_sim, 2),
      "\nMAE =", round(mae_obs_sim, 2), "°C"
    ),
    hjust = 0,  # Nudge left
    vjust = 1.5,    # Nudge down
    size = 5,       # Adjust text size
    color = "black",
    fontface = "bold"
  ) +
  xlim(c(0, 40)) +
  ylim(c(0, 40)) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "right",
    # legend.text = element_text(size = 12),
    # axis.title = element_text(size = 14),
    # plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed(ratio = 1)
print(p4)
ggsave("Plots/Silvilaser/musica_hobo.png", p4, width = 6, height = 6, dpi = 300)

# -------------------------- Figure Silvilaser2025 -----------------------------
# title <- textGrob(paste("(a): MuSICA predictions as a function of macroclimate measures for buffering plots,",
#                   "\n(b): MuSICA predictions as a function of macroclimate measures for amplifying plots,", 
#                   "\n(c): HOBO sensors measures as a function of macroclimate measures,", 
#                   "\n(d): MuSICA predictions as a function of HOBO sensors measures",
#                   "\nReferences: (1): 10.1016/j.rse.2023.113820"),
#                   gp = gpar(fontface = "bold", fontsize = 14))
# # Add ref git rmusica ?
# 
# g <- arrangeGrob(p1, p2, p3, p4, nrow = 2, ncol = 2,
#                  bottom = title)
# ggsave(file.path("Plots/Silvilaser/aabig_figure.png"), g, 
#        width = 10, height = 10, dpi = 300)

# For MuSICA (combined_sim)
mean_sim_vci <- combined_sim %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA)")

mean_sim_max <- combined_sim %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA)")

mean_sim_pai <- combined_sim %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA)")

# For HOBO (combined_obs_era5)
mean_obs_era5_vci <- combined_obs_era5 %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model")

mean_obs_era5_max <- combined_obs_era5 %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model")

mean_obs_era5_pai <- combined_obs_era5 %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model")

# Combine vci data from MuSICA and HOBO
mean_vci <- bind_rows(
  mean_sim_vci %>% mutate(microclimate = "Physics-based model (MuSICA)"),
  mean_obs_era5_vci %>% mutate(microclimate = "Statistical model")
)

# Combine max data from MuSICA and HOBO
mean_max <- bind_rows(
  mean_sim_max %>% mutate(microclimate = "Physics-based model (MuSICA)"),
  mean_obs_era5_max %>% mutate(microclimate = "Statistical model")
)

# Combine pai data from MuSICA and HOBO
mean_pai <- bind_rows(
  mean_sim_pai %>% mutate(microclimate = "Physics-based model (MuSICA)"),
  mean_obs_era5_pai %>% mutate(microclimate = "Statistical model")
)

plot_vci <- ggplot(filter(mean_vci, microclimate == "Statistical model"),
                   aes(x = vci, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Statistical model" = "#128d84")) +
  labs(x = "Vertical Complexity Index",
       # y = NULL,
       y = "    ",
       color = NULL) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = diff(range(mean_vci$vci)) / diff(range(mean_vci$mean_log_slope)))

plot_max <- ggplot(filter(mean_max, microclimate == "Statistical model"),
                   aes(x = max, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Statistical model" = "#128d84")) +
  labs(x = "Maximum Height (m)",
       y = "log(slope)",
       color = NULL) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = diff(range(mean_max$max)) / diff(range(mean_max$mean_log_slope)))

plot_pai <- ggplot(filter(mean_pai, microclimate == "Statistical model"),
                   aes(x = pai, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Statistical model" = "#128d84")) +
  labs(x = "Plant Area Index",
       # y = NULL,
       y = "    ",
       color = NULL) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = diff(range(mean_pai$pai)) / diff(range(mean_pai$mean_log_slope)))

# Arrange the plots
plog <- grid.arrange(plot_max, plot_vci, plot_pai, nrow = 1, widths = c(1, 1, 1))
ggsave("Plots/Silvilaser/comp1.png", plog, width = 10, height = 4, dpi = 300)

# 3. Create the plots for each variable (vci, max, pai)
plot_vci <- ggplot(mean_vci, aes(x = vci, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Physics-based model (MuSICA)" = "#E69F00", "Statistical model" = "#128d84")) +
  labs(
    # title = " ",
    x = "Vertical Complexity Index",
    # y = "log(slope)",
    y = "    ",
    color = NULL
  ) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2)) +
  coord_fixed(ratio = diff(range(mean_vci$vci)) / diff(range(mean_vci$mean_log_slope)))
# print(plot_vci)

plot_max <- ggplot(mean_max, aes(x = max, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Physics-based model (MuSICA)" = "#E69F00", "Statistical model" = "#128d84")) +
  labs(
    # title = " ",
    x = "Maximum Height (m)",
    y = "log(slope)",
    # y = NULL,
    color = NULL
  ) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2)) +
  coord_fixed(ratio = diff(range(mean_max$max)) / diff(range(mean_max$mean_log_slope)))
# coord_fixed(ratio = 51)

plot_pai <- ggplot(mean_pai, aes(x = pai, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Physics-based model (MuSICA)" = "#E69F00", "Statistical model" = "#128d84")) +
  labs(
    # title = " ",
    x = "Plant Area Index",
    y = "    ",
    # y = NULL,
    color = NULL
  ) +
  ylim(c(-0.4, 0.4)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2)) +
  coord_fixed(ratio = diff(range(mean_pai$pai)) / diff(range(mean_pai$mean_log_slope)))

# Arrange the updated plots in a grid
plog <- grid.arrange(plot_max, plot_vci, plot_pai, nrow = 1, widths = c(1, 1, 1))
ggsave("Plots/Silvilaser/compfull.png", plog, width = 10, height = 4, dpi = 300)
# plog <- plot_grid(plot_max, plot_vci, plot_pai, nrow = 1, align = "hv", rel_widths = c(1, 1, 1))

png("Plots/Silvilaser/aabig_figure.png", width = 18, height = 12, res = 300)
grid.draw(plog)
dev.off()

stop()

title <- richtext_grob(
  "<b>(a):</b> HOBO sensors measures as a function of macroclimate measures,<br><br>
   <b>(b):</b> MuSICA predictions as a function of HOBO sensors measuress,<br><br>
   <b>&#40;c&#41;:</b> MuSICA predictions as a function of macroclimate measures,<br><br>
   <b>(d):</b> Effect of LiDAR-derived variables related to canopy structure on the log of the slope between microclimate and macroclimate temperatures,<br><br>
   <i>References: (1): 10.1016/j.rse.2023.113820</i>",
  gp = gpar(fontsize = 12),
  hjust = 0,  # Left-align text
  x = unit(0.05, "npc")  # Avoid text clipping
)
# Add ref git rmusica ?

g <- arrangeGrob(p3, p4, p12, p2, nrow = 2, ncol = 2,
                 bottom = title)
g <- arrangeGrob(plot_vci, p3, p4, nrow = 1, ncol = 3,
                 bottom = title)
png("Plots/Silvilaser/aabig_figure.png", width = 10, height = 10, units = "in", res = 300)
grid.draw(g)
dev.off()

# ------------------------ Temperatures time-series ----------------------------

# Merge
# df_merged <- merge(df_Tair_hourly, t_hobo_hourly, by = c("id_plot", "time"), all = TRUE)
# df_merged <- df_merged %>%
#   group_by(time) %>%
#   summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
#   dplyr::select(-id_plot)
# df_merged <- merge(df_merged, era5_Tair_all, by = "time", all = TRUE)
# df_merged <- merge(df_merged, ign_Tair_all, by = "time", all = TRUE) 
# colnames(df_merged) <- c("time", "Tair_z", "T_hobo", "Tair_era5", "Tair_ign")

# Define the selected id_plot values
id_plot_selected <- c("41_07", "41_18", "41_19", "41_27", 
                      "41_30", "41_39", "41_47", "41_49")

# Split the datasets
df_with_plots <- df_Tair_hourly %>% filter(id_plot %in% id_plot_selected)
df_without_plots <- df_Tair_hourly %>% filter(!id_plot %in% id_plot_selected)
t_hobo_with_plots <- t_hobo_hourly %>% filter(id_plot %in% id_plot_selected)
t_hobo_without_plots <- t_hobo_hourly %>% filter(!id_plot %in% id_plot_selected)
df_merged_with <- merge_and_process(df_with_plots, t_hobo_with_plots)
df_merged_without <- merge_and_process(df_without_plots, t_hobo_without_plots)

# Calculate max and min for both datasets
extremes_with <- calculate_daily_extremes(df_merged_with)
extremes_without <- calculate_daily_extremes(df_merged_without)
plot_extremes(extremes_with$max, extremes_with$min, "With Selected id_plots")
plot_extremes(extremes_without$max, extremes_without$min, "Without Selected id_plots")

# Analyze hourly data for both cases
start_date <- as.POSIXct("2021-06-10 00:00:00")
end_date <- as.POSIXct("2021-06-17 23:00:00")
analyze_hourly_week(df_merged_with, start_date, end_date)
analyze_hourly_week(df_merged_without, start_date, end_date)