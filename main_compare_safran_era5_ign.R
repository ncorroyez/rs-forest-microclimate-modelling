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

read_musica_plot <- function(id_plot, source_type, date_seq) {
  
  dir_name <- if(source_type == "ERA5") "ERA5" else "Safran"
  file_pattern <- paste0("musica_out_", dir_name, "_Blois_2021_pt_", id_plot, ".nc")
  nc_file <- file.path("out_files/musica", dir_name, file_pattern)
  # nc_file <- file.path("out_files", file_pattern)
  
  if (!file.exists(nc_file)) {
    warning(paste("Fichier introuvable:", nc_file))
    return(NULL)
  }
  
  nc <- nc_open(nc_file)
  on.exit(nc_close(nc))
  raw_data <- get_variable(nc, "Tair_z")
  
  processed_data <- raw_data %>%
    filter(nair == 1) %>%
    mutate(
      # time = with_tz(time, tzone = "UTC"),
      time = lubridate::ymd_hms(time),
      date_only = as.Date(time)
    ) %>% 
    filter(date_only %in% date_seq) %>%
    mutate(
      Tair_z = Tair_z - 273.15,
      id_plot = id_plot,
      source = source_type
    ) %>%
    dplyr::select(-date_only)
  
  return(processed_data)
}

# Preparation
# date <- "27/06/2021"
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
# "41_34" misses between 2021-09-24 and 2021-09-30

# LAD stack
LAI_stack <- terra::rast("in_files/ladstack_classic_no_na.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson")

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
available_plots <- available_plots[available_plots != "41_34"]
# stop()
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
df_MuSICA_ERA5_Tair_all <- df_MuSICA_Safran_Tair_all <- data.frame()
df_MuSICA_ERA5_old_Tair_all <- data.frame()
# out.all <- list()

# date_study <- "2021-06-20"
# date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-06-30"), by = "day")

# ERA5 old
# for (i in 1:nrow(df_lad)) {
#   print(i)
#   id_plot <- n_field[i]
#   
#   # Open the NetCDF file for this field
#   nc_file <- paste0("out_files/radius_test/", radius,
#                     "m/musica_out_Blois_pt_", id_plot,
#                     "_radius_", radius, ".nc")
#   out <- nc_open(nc_file)
#   # out.all[[paste0("pt_", id_plot)]] <- out
#   
#   for (date_study in date_seq) {
#     # Extract the variable
#     df_Tair <- get_variable(out, "Tair_z") %>%
#       filter(as.Date(time) == date_study, nair == 1) %>%
#       mutate(Tair_z = Tair_z - 273.15,
#              id_plot = id_plot)
#     
#     # ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
#     df_MuSICA_ERA5_old_Tair_all <- bind_rows(df_MuSICA_ERA5_old_Tair_all, df_Tair)
#   }
# }
# df_MuSICA_ERA5_old_Tair_all <- df_MuSICA_ERA5_old_Tair_all %>%
#   mutate(time = floor_date(time, unit = "hour")) %>%
#   arrange(time) %>%
#   filter(id_plot %in% available_plots)

# ERA5
# for (i in 1:nrow(df_lad)) {
#   print(i)
#   id_plot <- n_field[i]
#   
#   # Open the NetCDF file for this field
#   nc_file <- paste0("out_files/musica/ERA5",
#                     "/musica_out_ERA5_Blois_2021_pt_", id_plot, ".nc")
#   out <- nc_open(nc_file)
#   # out.all[[paste0("pt_", id_plot)]] <- out
#   
#   for (date_study in date_seq) {
#     # Extract the variable
#     df_Tair <- get_variable(out, "Tair_z") %>%
#       filter(as.Date(time) == date_study, nair == 1) %>%
#       mutate(Tair_z = Tair_z - 273.15,
#              id_plot = id_plot)
#     
#     # ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
#     df_MuSICA_ERA5_Tair_all <- bind_rows(df_MuSICA_ERA5_Tair_all, df_Tair)
#   }
# }
# df_MuSICA_ERA5_Tair_all <- df_MuSICA_ERA5_Tair_all %>%
#   mutate(time = floor_date(time, unit = "hour")) %>%
#   arrange(time) %>%
#   filter(id_plot %in% available_plots)
df_MuSICA_ERA5_Tair_all <- map_dfr(n_field, read_musica_plot, 
                                   source_type = "ERA5", 
                                   date_seq = date_seq) %>%
  mutate(time = floor_date(time, unit = "hour")) %>%
  arrange(time) %>%
  filter(id_plot %in% available_plots)
print(attr(df_MuSICA_ERA5_Tair_all$time, "tzone"))

# SAFRAN
# for (i in 1:nrow(df_lad)) {
#   print(i)
#   id_plot <- n_field[i]
#   
#   # Open the NetCDF file for this field
#   nc_file <- paste0("out_files/musica/Safran",
#                     "/musica_out_Safran_Blois_2021_pt_", id_plot, ".nc")
#   out <- nc_open(nc_file)
#   # out.all[[paste0("pt_", id_plot)]] <- out
#   
#   for (date_study in date_seq) {
#     # Extract the variable
#     df_Tair <- get_variable(out, "Tair_z") %>%
#       filter(as.Date(time) == date_study, nair == 1) %>%
#       mutate(Tair_z = Tair_z - 273.15,
#              id_plot = id_plot)
#     
#     # ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
#     df_MuSICA_Safran_Tair_all <- bind_rows(df_MuSICA_Safran_Tair_all, df_Tair)
#   }
# }
# df_MuSICA_Safran_Tair_all <- df_MuSICA_Safran_Tair_all %>%
#   mutate(time = floor_date(time, unit = "hour")) %>%
#   arrange(time) %>%
#   filter(id_plot %in% available_plots)
df_MuSICA_Safran_Tair_all <- map_dfr(n_field, read_musica_plot, 
                                     source_type = "Safran", 
                                     date_seq = date_seq) %>%
  mutate(time = floor_date(time, unit = "hour")) %>%
  arrange(time) %>%
  filter(id_plot %in% available_plots)
print(attr(df_MuSICA_Safran_Tair_all$time, "tzone"))

stop()
# CHS41
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
era5_var <- era5_var %>%
  mutate(time = lubridate::force_tz(time, tzone = "UTC"))

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

# SAFRAN
safran <- nc_open("in_files/musica_in_Safran_Blois_2021.nc")
safran_var <- get_variable(safran, "Tair")
safran_var <- safran_var %>%
  mutate(time = lubridate::force_tz(time, tzone = "UTC"))

safran_Tair_all <- data.frame()
for (date_study in date_seq) {
  date_study <- as.Date(date_study)
  
  safran_Tair_filtered <- safran_var %>%
    filter(format(time, "%Y-%m-%d") == format(date_study, "%Y-%m-%d")) %>%
    mutate(Tair = round(Tair - 273.15, 1)) %>%
    mutate(time = floor_date(time, unit = "hour")) # %>%
  # dplyr::select(-z)
  
  safran_Tair_all <- bind_rows(safran_Tair_all, safran_Tair_filtered)
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
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "ERA5 (°C)",
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
    legend.position = "right"
  ) +
  theme_bw(base_size = 16) +
  coord_fixed(ratio = 1)
print(p)

# CHS41 vs SAFRAN
combined_data <- merge(safran_Tair_all, ign_Tair_all, by = "time", suffixes = c("_safran", "_ign"))
r_value <- cor(combined_data$Tair_safran, combined_data$Tair_safran)
bias <- mean(combined_data$Tair_ign - combined_data$Tair_safran)
nrmse <- sqrt(mean((combined_data$Tair_ign - combined_data$Tair_safran)^2)) / IQR(combined_data$Tair_ign)

# Fit a linear model
model <- lm(Tair_ign ~ Tair_safran, data = combined_data)
r2 <- summary(model)$r.squared
mae <- mean(abs(combined_data$Tair_ign - combined_data$Tair_safran))
slope <- coef(model)[2]
intercept <- coef(model)[1]

p <- ggplot(combined_data, aes(x = Tair_ign, y = Tair_safran)) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Local weather station (°C)",
       y = "SAFRAN (°C)",
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
    legend.position = "right"
  ) +
  theme_bw(base_size = 16) +
  coord_fixed(ratio = 1)
print(p)

# ERA5 vs SAFRAN
combined_data <- merge(era5_Tair_all, safran_Tair_all, by = "time", suffixes = c("_era5", "_safran"))
r_value <- cor(combined_data$Tair_era5, combined_data$Tair_safran)
bias <- mean(combined_data$Tair_safran - combined_data$Tair_era5)
nrmse <- sqrt(mean((combined_data$Tair_safran - combined_data$Tair_era5)^2)) / IQR(combined_data$Tair_safran)

# Fit a linear model
model <- lm(Tair_safran ~ Tair_era5, data = combined_data)
r2 <- summary(model)$r.squared
mae <- mean(abs(combined_data$Tair_safran - combined_data$Tair_era5))
slope <- coef(model)[2]
intercept <- coef(model)[1]

p <- ggplot(combined_data, aes(x = Tair_safran, y = Tair_era5)) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c(option = "viridis") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SAFRAN (°C)",
       y = "ERA5 (°C)",
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
    legend.position = "right"
  ) +
  theme_bw(base_size = 16) +
  coord_fixed(ratio = 1)
print(p)

# ------------------------------- MuSICA Plots ---------------------------------
df_MuSICA_Safran_Tair_hourly <- df_MuSICA_Safran_Tair_all %>%
  mutate(time = floor_date(time, "hour"))
df_MuSICA_ERA5_Tair_hourly <- df_MuSICA_ERA5_Tair_all %>%
  mutate(time = floor_date(time, "hour"))
# df_MuSICA_ERA5_old_Tair_hourly <- df_MuSICA_ERA5_old_Tair_all %>%
#   mutate(time = floor_date(time, "hour"))
t_hobo_hourly <- temperature_all_days %>%
  mutate(time = floor_date(datetime, "hour"))

# Microclimate: MuSICA forced ERA5, macroclimate: ERA5
combined_ERA5_hobo <- merge(merge(df_MuSICA_ERA5_Tair_hourly, era5_Tair_all,
                                  by = c("time")), metrics_df, 
                            by = "id_plot") %>% mutate(hour = hour(time))

# combined_ERA5_hobo <- combined_ERA5_hobo %>%
#   filter(hour %in% c(7:21))

r_value_obs_sim <- cor(combined_ERA5_hobo$Tair, combined_ERA5_hobo$Tair_z)
bias_obs_sim <- mean(combined_ERA5_hobo$Tair_z - combined_ERA5_hobo$Tair)
nrmse_obs_sim <- sqrt(mean((combined_ERA5_hobo$Tair - combined_ERA5_hobo$Tair_z)^2)) / IQR(combined_ERA5_hobo$Tair)
mae_obs_sim <- mean(abs(combined_ERA5_hobo$Tair_z - combined_ERA5_hobo$Tair))
r_squared_obs_sim <- summary(lm(Tair_z ~ Tair, 
                                data = combined_ERA5_hobo))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_ERA5_hobo <- combined_ERA5_hobo %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(Tair_z ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(Tair_z ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = Tair_z - Tair
  ) %>%
  ungroup()

p1 <- ggplot(combined_ERA5_hobo, aes(x = Tair, y = Tair_z)) +
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
print(p1)

# Microclimate: MuSICA forced SAFRAN, macroclimate: SAFRAN
combined_Safran_hobo <- merge(merge(df_MuSICA_Safran_Tair_hourly, safran_Tair_all,
                                    by = c("time")), metrics_df, 
                              by = "id_plot") %>% mutate(hour = hour(time))

# combined_Safran_hobo <- combined_Safran_hobo %>%
#   filter(hour %in% c(7:21))

r_value_obs_sim <- cor(combined_Safran_hobo$Tair, combined_Safran_hobo$Tair_z)
bias_obs_sim <- mean(combined_Safran_hobo$Tair_z - combined_Safran_hobo$Tair)
nrmse_obs_sim <- sqrt(mean((combined_Safran_hobo$Tair - combined_Safran_hobo$Tair_z)^2)) / IQR(combined_Safran_hobo$Tair)
mae_obs_sim <- mean(abs(combined_Safran_hobo$Tair_z - combined_Safran_hobo$Tair))
r_squared_obs_sim <- summary(lm(Tair_z ~ Tair, 
                                data = combined_Safran_hobo))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_Safran_hobo <- combined_Safran_hobo %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(Tair_z ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(Tair_z ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = Tair_z - Tair
  ) %>%
  ungroup()

p2 <- ggplot(combined_Safran_hobo, aes(x = Tair, y = Tair_z)) +
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
print(p2)

# Microclimate: HOBO, macroclimate: ERA5
combined_obs_ERA5_hobo <- merge(merge(t_hobo_hourly, era5_Tair_all,
                                    by =  "time"), metrics_df, 
                              by = "id_plot") %>% mutate(hour = hour(time))

# combined_obs_ERA5_hobo <- combined_obs_ERA5_hobo %>%
#   filter(hour %in% c(7:21))

r_value_obs_sim <- cor(combined_obs_ERA5_hobo$t_hobo, combined_obs_ERA5_hobo$Tair)
bias_obs_sim <- mean(combined_obs_ERA5_hobo$Tair - combined_obs_ERA5_hobo$t_hobo)
nrmse_obs_sim <- sqrt(mean((combined_obs_ERA5_hobo$t_hobo - combined_obs_ERA5_hobo$Tair)^2)) / IQR(combined_Safran_hobo$t_hobo)
mae_obs_sim <- mean(abs(combined_obs_ERA5_hobo$Tair - combined_obs_ERA5_hobo$t_hobo))
r_squared_obs_sim <- summary(lm(t_hobo ~ Tair, 
                                data = combined_obs_ERA5_hobo))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_obs_ERA5_hobo <- combined_obs_ERA5_hobo %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(t_hobo ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(t_hobo ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = t_hobo - Tair
  ) %>%
  ungroup()

# Microclimate: HOBO, macroclimate: SAFRAN
combined_obs_Safran_hobo <- merge(merge(t_hobo_hourly, safran_Tair_all,
                                        by =  "time"), metrics_df, 
                                  by = "id_plot") %>% mutate(hour = hour(time))

# combined_obs_Safran_hobo <- combined_obs_Safran_hobo %>%
#   filter(hour %in% c(7:21))

r_value_obs_sim <- cor(combined_obs_Safran_hobo$t_hobo, combined_obs_Safran_hobo$Tair)
bias_obs_sim <- mean(combined_obs_Safran_hobo$Tair - combined_obs_Safran_hobo$t_hobo)
nrmse_obs_sim <- sqrt(mean((combined_obs_Safran_hobo$t_hobo - combined_obs_Safran_hobo$Tair)^2)) / IQR(combined_Safran_hobo$t_hobo)
mae_obs_sim <- mean(abs(combined_obs_Safran_hobo$Tair - combined_obs_Safran_hobo$t_hobo))
r_squared_obs_sim <- summary(lm(t_hobo ~ Tair, 
                                data = combined_obs_Safran_hobo))$r.squared

# Fit a linear model and calculate log(slope) for each point
combined_obs_Safran_hobo <- combined_obs_Safran_hobo %>%
  group_by(id_plot) %>%
  mutate(
    slope = coef(lm(t_hobo ~ Tair, data = cur_data()))[2],
    equilibrium = coef(lm(t_hobo ~ Tair, data = cur_data()))[1] / (1 - slope),
    log_slope = log(abs(slope)),
    offset = t_hobo - Tair
  ) %>%
  ungroup()

# For MuSICA ERA5
mean_musica_era5_vci <- combined_ERA5_hobo %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, ERA5)")

mean_musica_era5_max <- combined_ERA5_hobo %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, ERA5)")

mean_musica_era5_pai <- combined_ERA5_hobo %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, ERA5)")

# For MuSICA SAFRAN
mean_musica_safran_vci <- combined_Safran_hobo %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)")

mean_musica_safran_max <- combined_Safran_hobo %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)")

mean_musica_safran_pai <- combined_Safran_hobo %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)")

# For HOBO (ERA5)
mean_obs_era5_vci <- combined_obs_ERA5_hobo %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (ERA5)")

mean_obs_era5_max <- combined_obs_ERA5_hobo %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (ERA5)")

mean_obs_era5_pai <- combined_obs_ERA5_hobo %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (ERA5)")

# For HOBO (SAFRAN)
mean_obs_safran_vci <- combined_obs_Safran_hobo %>%
  group_by(id_plot, vci) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (SAFRAN)")

mean_obs_safran_max <- combined_obs_Safran_hobo %>%
  group_by(id_plot, max) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (SAFRAN)")

mean_obs_safran_pai <- combined_obs_Safran_hobo %>%
  group_by(id_plot, pai) %>%
  summarise(mean_log_slope = mean(log_slope, na.rm = TRUE), .groups = "drop") %>%
  dplyr::select(-id_plot) %>%
  mutate(microclimate = "Statistical model (SAFRAN)")

# Combine vci data from MuSICA and HOBO
mean_vci <- bind_rows(
  mean_musica_era5_vci   %>% mutate(microclimate = "Physics-based model (MuSICA, ERA5)"),
  mean_musica_safran_vci %>% mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)"),
  mean_obs_era5_vci      %>% mutate(microclimate = "Statistical model (ERA5)"),
  mean_obs_safran_vci    %>% mutate(microclimate = "Statistical model (SAFRAN)")
)

# Combine max data from MuSICA and HOBO
mean_max <- bind_rows(
  mean_musica_era5_max   %>% mutate(microclimate = "Physics-based model (MuSICA, ERA5)"),
  mean_musica_safran_max %>% mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)"),
  mean_obs_era5_max      %>% mutate(microclimate = "Statistical model (ERA5)"),
  mean_obs_safran_max    %>% mutate(microclimate = "Statistical model (SAFRAN)")
)

# Combine pai data from MuSICA and HOBO
mean_pai <- bind_rows(
  mean_musica_era5_pai   %>% mutate(microclimate = "Physics-based model (MuSICA, ERA5)"),
  mean_musica_safran_pai %>% mutate(microclimate = "Physics-based model (MuSICA, SAFRAN)"),
  mean_obs_era5_pai      %>% mutate(microclimate = "Statistical model (ERA5)"),
  mean_obs_safran_pai    %>% mutate(microclimate = "Statistical model (SAFRAN)")
)


# 3. Create the plots for each variable (vci, max, pai)
plot_vci <- ggplot(mean_vci, aes(x = vci, y = mean_log_slope, color = microclimate)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 1.2) +
  scale_color_manual(values = c("Physics-based model (MuSICA, ERA5)" = "#E69F00",
                                "Physics-based model (MuSICA, SAFRAN)" = "red",
                                "Statistical model (ERA5)" = "#128d84",
                                "Statistical model (SAFRAN)" = "blue")) +
  labs(
    # title = " ",
    x = "Vertical Complexity Index",
    # y = "log(slope)",
    y = "    ",
    color = NULL
  ) +
  # ylim(c(-0.4, 0.4)) +
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
  scale_color_manual(values = c("Physics-based model (MuSICA, ERA5)" = "#E69F00",
                                "Physics-based model (MuSICA, SAFRAN)" = "red",
                                "Statistical model (ERA5)" = "#128d84",
                                "Statistical model (SAFRAN)" = "blue")) +
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
  scale_color_manual(values = c("Physics-based model (MuSICA, ERA5)" = "#E69F00",
                                "Physics-based model (MuSICA, SAFRAN)" = "red",
                                "Statistical model (ERA5)" = "#128d84",
                                "Statistical model (SAFRAN)" = "blue")) +
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
print(plog)

# Microclimate: y: MuSICA (ERA5), x: