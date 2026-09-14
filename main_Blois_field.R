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
library(lidR)
library(raster)
library(sf)
library(gstat)
library(terra)
library(lubridate)

# Initialize an empty list to store the results for each file
results <- list()
allometry_list <- list()
musica_list <- list()
out.all <- list()
microclimate_height <- 1
inter_crown_clumping <- 0.7 # Loop w/ =/= CCI ?
musica.cmd <- "bash -i -c musica"
sites <- c("Blois")

ctg <- readLAScatalog("/media/corroyez/My Passport/01_DATA/Blois/LiDAR/2-las_utm")
LAI_stack <- terra::rast("in_files/ladstack_classic.tif")
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
# radius_list <- c(5, 10, 20)
radius_list <- c(20)
# df_lad[is.na(df_lad)] <- 0.0

for (radius in radius_list){
# df_lad <- df_lad[, -c(ncol(df_lad)-1, ncol(df_lad))]
# for (i in 1:nrow(df_lad)) {
  for (i in seq_len(length(n_field))){
  # lad_values <- as.numeric(df_lad[i, ])
  # allometry <- data.frame(
  #   height = height_levels,
  #   density = lad_values
  # )
  # allometry <- allometry %>%
  #   filter(!is.na(height) & !is.na(density))
  # 
  # canopy_height_top <- min(40, max(allometry$height[allometry$density > 0] + 0.5, na.rm = TRUE))
  
  pt <- clip_circle(ctg,
                    xcenter = df_coords[i,]$coord_x_utm31n,
                    ycenter = df_coords[i,]$coord_y_utm31n,
                    radius = radius
  )

  dtm <- grid_terrain(pt, 1, kriging(k = 10L))
  pt <- normalize_height(pt, dtm)
  chm <- grid_canopy(pt, res = 1, dsmtin())

  lad <- LAD(pt@data$Z)
  allometry <- data.frame(
    height = lad$z,
    density = lad$lad
  )

  canopy_height_top <- max(cellStats(chm, "max"), allometry$height + 1)
  
  PAI <- 2 * sum(allometry$density, na.rm = TRUE)
  
  phenology <- calc_phenology(list.year = c(2021, 2022),
                              nleafage = 3,
                              budburst_date = 115,
                              leaf_age_max_in = 0.56,
                              relative_age_firstmax = 0.10,
                              relative_age_lastmax = 0.75,
                              LAI_max_per_cohort = PAI)
  
  musica_out <- callmusica_vegetation_structure(
    crown_clumping_factor = inter_crown_clumping,
    microclimate_height = microclimate_height,
    leaf_param = list("musica_veg1" =
                        list(phenology = phenology,
                             allometry = allometry,
                             PAI = PAI,
                             canopy_height_top = canopy_height_top,
                             canopy_height_bottom = 2
                        )),
    run.name = paste0("pt_", n_field[i], "_radius_", radius),
    # run.name = paste0("pt_", n_field[i]),
    musica.cmd = musica.cmd,
    save.output = TRUE,
    keep.tmp = FALSE)
  
  # File copy
    dir.create(paste0("out_files/radius_test/", radius, "m"),
               recursive = TRUE, showWarnings = FALSE)
    file_path <- paste0("out_files/musica_out_Blois_pt_", n_field[i],
                        "_radius_", radius, ".nc")
    new_file_path <- paste0("out_files/radius_test/", radius,
                            "m/musica_out_Blois_pt_", n_field[i],
                            "_radius_", radius, ".nc")
    if (file.copy(file_path, new_file_path)) {
      file.remove(file_path)
    }
  }
}

# Tair_z == 1 (0-1.3m,field = 1m), Tsoil == 4 (68.9-102.7mm, field = 80mm)

# --------------------------------- Air ----------------------------------------
# for (n_field in seq_len(dim(df_coords)[1])){
df_Tair_all <- data.frame()
radius <- 20
# Loop over each field
for (i in 1:nrow(df_lad)) {
  id_plot <- n_field[i]
  
  # Open the NetCDF file for this field
  # nc_file <- paste0("out_files/musica_out_Blois_pt_", id_plot, ".nc")
  nc_file <- paste0("out_files/radius_test/", radius, "m/",
                    "musica_out_Blois_pt_", id_plot,
                    "_radius_", radius, ".nc")
  out <- nc_open(nc_file)
  out.all[[paste0("pt_", id_plot)]] <- out
  
  # Extract the variable
  df_Tair <- get_variable(out, "Tair_z") %>%
    filter(as.Date(time) == "2021-06-17", nair == 1) %>%
    mutate(Tair_z = Tair_z - 273.15,
           id_plot = id_plot)
  
  ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
  df_Tair_all <- bind_rows(df_Tair_all, df_Tair)
}
# shinymusica(out.all)
df_Tair_all <- df_Tair_all %>%
  mutate(time = floor_date(time, unit = "hour"))

temperature_all <- data.frame()
# for (site in sites) {
# Read the CSV file for the site
site <- "Blois"
csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
temperature_csv <- read.csv(csv_file)

# Convert the datetime column to POSIXct
temperature_df <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", 
                               tz = "UTC"))

# Filter for the desired date and sensor position (e.g., "a")
temperature_df <- temperature_df %>%
  filter(as.Date(datetime) == "2021-06-17",
         id_plot %in% df_Tair_all$id_plot,
         position_sensor == "a"
  )

missing_in_temperature <- setdiff(df_Tair_all$id_plot, temperature_df$id_plot)
print(missing_in_temperature)

# Find common IDs
common_ids <- intersect(df_Tair_all$id_plot, temperature_df$id_plot)

# Filter both data frames to keep only the common IDs
df_Tair_all <- df_Tair_all %>% filter(id_plot %in% common_ids)
temperature_df <- temperature_df %>% filter(id_plot %in% common_ids)

# Filter and rename
df_Tair_all <- df_Tair_all %>%
  dplyr::select(id_plot, datetime = time, Tair_z)
temperature_df <- temperature_df %>%
  dplyr::select(id_plot, datetime, Tair_z = t_hobo)

# Analysis
df_joined <- left_join(df_Tair_all, temperature_df,
                       by = c("id_plot", "datetime"),
                       suffix = c("_sim", "_obs"))
df_joined <- df_joined %>%
  filter(id_plot != "41_28")

# Model
model <- lm(Tair_z_obs ~ Tair_z_sim, data = df_joined)
model_summary <- summary(model)

# Metrics
slope <- coef(model)[2]
intercept <- coef(model)[1]
r_value <- cor(df_joined$Tair_z_sim, df_joined$Tair_z_obs, use = "complete.obs")
r2 <- model_summary$r.squared
rmse <- sqrt(mean((df_joined$Tair_z_obs - df_joined$Tair_z_sim)^2, 
                  na.rm = TRUE))
obs_iqr <- IQR(df_joined$Tair_z_obs, na.rm = TRUE)
nrmse <- rmse / obs_iqr
bias <- mean(df_joined$Tair_z_sim - df_joined$Tair_z_obs, na.rm = TRUE)

# Plot
ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
  geom_point(aes(color = id_plot), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Observed Tair_z (°C)",
       y = "Simulated Tair_z (°C)",
       title = "Scatterplot: Observed vs. Simulated Tair_z",
       subtitle = paste("R =", round(r_value, 2),
                        "| R² =", round(r2, 2),
                        "| NRMSE =", round(nrmse, 2),
                        "| Bias =", round(bias, 2),
                        "| Slope =", round(slope, 2),
                        "| Intercept =", round(intercept, 2))) +
  # xlim(17, 21) +
  # ylim(17, 21) +
  theme_bw()

# Append to the overall CSV data frame
# temperature_all <- bind_rows(temperature_all, temperature_df)
# }

# id_plot facet
ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ id_plot, ncol = 4) +  # Adjust as needed
  labs(x = "Observed Tair_z (°C)",
       y = "Simulated Tair_z (°C)",
       title = "Scatterplot: Observed vs. Simulated Tair_z by Field") +
  theme_bw()

# datetime facet
ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ datetime, scales = "free_x", ncol = 4) +  # Facet by full datetime
  labs(x = "Observed Tair_z (°C)",
       y = "Simulated Tair_z (°C)",
       title = "Scatterplot: Observed vs. Simulated Tair_z by Time") +
  theme_bw()

# Color hours
df_joined <- df_joined %>%
  mutate(hour = hour(datetime))
ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim, color = factor(hour))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ id_plot, ncol = 4) +  # One scatterplot per `id_plot`
  scale_color_viridis_d() +  # Better color scheme for hours
  labs(x = "Observed Tair_z (°C)",
       y = "Simulated Tair_z (°C)",
       title = "Scatterplot: Observed vs. Simulated Tair_z",
       subtitle = "Colored by Hour of the Day",
       color = "Hour") +
  theme_bw()

# Loop over each unique `id_plot`
plot_list <- list()
for (plot_id in unique(df_joined$id_plot)) {
  p <- ggplot(df_joined %>% filter(id_plot == plot_id), 
              aes(x = Tair_z_obs, y = Tair_z_sim, color = factor(hour))) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_color_viridis_d() +  
    labs(x = "Observed Tair_z (°C)",
         y = "Simulated Tair_z (°C)",
         title = paste("Scatterplot for id_plot:", plot_id),
         subtitle = "Colored by Hour of the Day",
         color = "Hour") +
    theme_bw()
  
  # Store the plot in the list
  plot_list[[plot_id]] <- p
}

# Print all plots
for (p in plot_list) {
  print(p)
}

# --------------------------------- Soil ---------------------------------------
# for (n_field in seq_len(dim(df_coords)[1])){
df_T_soil_all <- data.frame()
# Loop over each field
for (i in 1:nrow(df_lad)) {
  id_plot <- n_field[i]
  
  # Open the NetCDF file for this field
  nc_file <- paste0("out_files/musica_out_Blois_pt_", id_plot, ".nc")
  out <- nc_open(nc_file)
  out.all[[paste0("pt_", id_plot)]] <- out
  
  # Extract the variable
  df_T_soil <- get_variable(out, "T_soil") %>%
    filter(as.Date(time) == "2021-06-17", nsoil == 4) %>%
    mutate(T_soil = T_soil - 273.15,
           id_plot = id_plot)
  
  ggplot_variable(filter(df_T_soil, nsoil == 4), out.type = "standard")
  df_T_soil_all <- bind_rows(df_T_soil_all, df_T_soil)
}
df_T_soil_all <- df_T_soil_all %>%
  mutate(time = floor_date(time, unit = "hour"))

soil_temperature_all <- data.frame()
# for (site in sites) {
# Read the CSV file for the site
site <- "Blois"
csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
temperature_csv <- read.csv(csv_file)

# Convert the datetime column to POSIXct
soil_temperature_df <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", 
                               tz = "UTC"))

# Filter for the desired date and sensor position (e.g., "a")
soil_temperature_df <- soil_temperature_df %>%
  filter(as.Date(datetime) == "2021-06-17",
         id_plot %in% df_T_soil_all$id_plot,
         position_sensor == "s"
  )

# Floor hours to get even
# soil_temperature_df <- soil_temperature_df %>%
#   mutate(datetime = datetime - hours(hour(datetime) %% 2))

missing_in_temperature <- setdiff(df_T_soil_all$id_plot, soil_temperature_df$id_plot)
print(missing_in_temperature)

# Find common IDs
common_ids <- intersect(df_T_soil_all$id_plot, soil_temperature_df$id_plot)

# Filter both data frames to keep only the common IDs
df_T_soil_all <- df_T_soil_all %>% filter(id_plot %in% common_ids)
soil_temperature_df <- soil_temperature_df %>% filter(id_plot %in% common_ids)

# Filter and rename
df_T_soil_all <- df_T_soil_all %>%
  dplyr::select(id_plot, datetime = time, T_soil)
soil_temperature_df <- soil_temperature_df %>%
  dplyr::select(id_plot, datetime, T_soil = t_hobo)

# Analysis
df_soil_joined <- left_join(df_T_soil_all, soil_temperature_df,
                            by = c("id_plot", "datetime"),
                            suffix = c("_sim", "_obs"))
df_soil_joined <- df_soil_joined %>%
  filter(!is.na(T_soil_obs))
# df_soil_joined <- df_soil_joined %>%
#   filter(id_plot != "41_28")

# Model
model <- lm(T_soil_obs ~ T_soil_sim, data = df_soil_joined)
model_summary <- summary(model)

# Metrics
slope <- coef(model)[2]
intercept <- coef(model)[1]
r_value <- cor(df_soil_joined$T_soil_sim, df_soil_joined$T_soil_obs, use = "complete.obs")
r2 <- model_summary$r.squared
rmse <- sqrt(mean((df_soil_joined$T_soil_obs - df_soil_joined$T_soil_sim)^2, 
                  na.rm = TRUE))
obs_iqr <- IQR(df_soil_joined$T_soil_obs, na.rm = TRUE)
nrmse <- rmse / obs_iqr
bias <- mean(df_soil_joined$T_soil_sim - df_soil_joined$T_soil_obs, na.rm = TRUE)

# Plot
ggplot(df_soil_joined, aes(x = T_soil_sim, y = T_soil_obs)) +
  geom_point(aes(color = id_plot), size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Simulated T_soil (°C)",
       y = "Observed T_soil (°C)",
       title = "Simulated vs. Observed T_soil",
       subtitle = paste("R =", round(r_value, 2),
                        "| R² =", round(r2, 2),
                        "| NRMSE =", round(nrmse, 2),
                        "| Bias =", round(bias, 2),
                        "| Slope =", round(slope, 2),
                        "| Intercept =", round(intercept, 2))) +
  # xlim(18, 32) +
  # ylim(18, 32) +
  theme_bw()

# Append to the overall CSV data frame
# temperature_all <- bind_rows(temperature_all, temperature_df)
# }

