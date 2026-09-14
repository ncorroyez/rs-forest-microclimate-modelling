rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

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
library(lmom)
library(doParallel)
library(foreach)
library(purrr)
library(exactextractr)

# force_utc_nc <- function(ncfile, varname = "time") {
#   if (!file.exists(ncfile)) stop(paste("Fichier introuvable:", ncfile))
# 
#   nc <- nc_open(ncfile)
#   on.exit(nc_close(nc)) # Fermeture automatique garantie
# 
#   # Récupération des données brutes
#   time_val <- ncvar_get(nc, varname)
#   units_str <- ncatt_get(nc, varname, "units")$value
# 
#   # Parsing robuste de la chaîne "hours since YYYY-MM-DD..."
#   parts <- strsplit(units_str, " since ")[[1]]
#   if (length(parts) < 2) stop("Format de date NetCDF non reconnu")
# 
#   unit_type <- parts[1] # ex: "hours"
#   origin_str <- parts[2] # ex: "2021-01-01 00:00:00"
# 
#   # Nettoyage : retirer d'éventuels suffixes comme "(GMT+00)" ou "UTC" à la fin
#   origin_str <- gsub("\\s*\\(.*\\)", "", origin_str)
#   origin_str <- gsub(" UTC", "", origin_str)
# 
#   # Conversion de l'origine en POSIXct UTC strict
#   origin_dt <- as.POSIXct(origin_str, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
# 
#   # Facteur de conversion en secondes
#   mult <- switch(trimws(unit_type),
#                  "hours" = 3600,
#                  "minutes" = 60,
#                  "seconds" = 1,
#                  "days" = 86400,
#                  3600) # Défaut
# 
#   # Calcul final et forçage UTC
#   final_time <- origin_dt + (time_val * mult)
#   return(with_tz(final_time, "UTC"))
# }

# Initialize an empty list to store the results for each file
results <- list()
allometry_list <- list()
musica_list <- list()
out.all <- list()
metrics_list <- list()
microclimate_height <- 1
# inter_crown_clumping <- 0.7 # Loop w/ =/= CCI ?
# musica.cmd <- "bash -i -c musica"
musica.cmd <- "bash -i -c musica > musica.log"
sites <- c("Blois")
# forcing_types <- c("ERA5", "Safran")
forcing_types <- c("ERA5")
results_dir <- './03_RESULTS'

ctg <- readLAScatalog("/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/2-las_utm")
LAI_stack <- terra::rast("in_files/ladstack_classic_no_na.tif")
clumping <- terra::rast("in_files/fCover_res_10_m.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson")
# ids_to_remove <- c("41_13", "41_14", "41_20", "41_34", 
#                    "41_41", "41_50", "41_51", "41_53",
#                    "41_17", "41_18", "41_19", "41_27",
#                    "41_30", "41_39", "41_47", "41_49",
#                    "41_55")
# json <- json %>%
#   filter(!id_plot %in% ids_to_remove) %>%
#   arrange(id_plot)

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

df_lad$clumping_value <- terra::extract(
  clumping,
  terra::vect(df_lad, geom = c("coord_x_utm31n", "coord_y_utm31n"), crs = crs(clumping))
)[, 2]
# df_lad$clumping_value <- 0.8

n_field <- df_lad$ID
# stop()
# n_field <- n_field[[1]]
df_lad$ID = NULL

# lai <- terra::rast("in_files/lidarlai_res_10_m.tif")
# lai_values <- terra::extract(lai, df_coords)

# distribs <- c(
#               'atbd',
#               'atbd_optim_common',
#               'atbd_optim_common_brownmodif',
#               'optim_Blois'
#               # 'optim_Blois_brownmodif'
#               )
distribs <- c('atbd') # atbd atbd_optim_common

for (distrib in distribs){
  
  path_s2_corr <- file.path("../NC_Full/03_RESULTS",
                            "Blois/Metrics/Not_Masked/Smooth_TS",
                            "smooth",
                            # "RF/LAI_ALS",
                            "gap",
                            paste0("LAI_Gapfilled_10d_", distrib, ".tif"))
  
  if (!file.exists(path_s2_corr)) stop("Raster S2 corrigé introuvable")
  r_s2_stack <- terra::rast(path_s2_corr)
  dates_s2 <- terra::time(r_s2_stack)
  if (is.null(dates_s2)) stop("Le raster S2 n'a pas de dates associées")
  
  # distrib <- 'lai_4'
  # distrib <- 'lidar_constant'
  # radius_list <- c(5, 10, 15, 20, 25, 50)
  radius_list <- c(25)
  
  for (radius in radius_list) {
    cat("Processing radius:", radius, "\n")
    for (forcing_type in forcing_types) {
      for (i in seq_len(length(n_field))) {
        
        pt <- clip_circle(ctg,
                          xcenter = df_coords[i,]$coord_x_utm31n,
                          ycenter = df_coords[i,]$coord_y_utm31n,
                          radius = radius
        )
        
        # Metrics
        # lmoms <- cloud_metrics(pt, ~as.list(lmom::samlmu(Z, nmom=3, ratios=F)))
        # lcv <- lmoms[[2]] / lmoms[[1]]
        # lskew <- lmoms[[3]] / lmoms[[2]]
        # vci <- cloud_metrics(pt, ~VCI(Z, zmax = max(Z)))
        # gf <- npoints(filter_ground(pt)) / npoints(pt)
        # 
        # # Normalize height and compute CHM
        dtm <- grid_terrain(pt, 1, tin(extrapolate = kriging()))
        pt <- normalize_height(pt, dtm)
        chm <- grid_canopy(pt, res = 1, pitfree(thresholds = c(0, 10, 20),
                                                max_edge = c(0, 1)
                                                # subcircle = 0.5
        ))
        mean <- cellStats(chm, "mean")
        
        # LAD and allometry
        lad <- LAD(pt@data$Z)
        lad_profile <- lad$lad
        sum_lad <- sum(lad_profile, na.rm = TRUE)
        
        if (sum_lad > 0) {
          lad_normalized <- lad_profile / sum_lad
        } else {
          lad_normalized <- lad_profile
        }
        
        allometry <- data.frame(
          height = lad$z,
          density = lad_normalized 
        )
        
        # MuSICA parameters
        canopy_height_top <- max(cellStats(chm, "max"), allometry$height + 1)
        # canopy_height_top <- max(allometry$height)
        plantareaindex <- round(2 * sum(allometry$density, na.rm = TRUE), 2)
        pai_real <- round(2 * sum_lad, 2)
        clumping <- df_lad$clumping_value[i]
        
        if (clumping < 0.5){
          if (mean < 2 || pai_real < 1){
            clumping <- 1
          } 
          else {
            clumping <- 0.5
          }
        }
        
        
        pt_vect <- vect(df_coords[i, c("coord_x_utm31n", "coord_y_utm31n")], 
                        geom = NULL, crs = crs(r_s2_stack))
        pt_buffer <- buffer(pt_vect, width = radius)
        pt_sf <- st_as_sf(pt_buffer)
        vals_s2_df <- exact_extract(r_s2_stack, pt_sf, fun = 'mean', progress = FALSE)
        vals_s2 <- as.numeric(vals_s2_df)
        
        # Création d'un vecteur de dates journalier pour l'année de simulation (2021)
        # MuSICA attend généralement une phénologie journalière (365 jours)
        sim_days <- seq(as.Date("2021-01-01"), as.Date("2021-12-31"), by = "day")
        
        # Interpolation linéaire des valeurs S2 (10 jours) vers Journalier
        # rule = 2 permet de prolonger la valeur (plateau) si la simu dépasse les dates S2 (janvier/décembre)
        # daily_lai <- approx(x = dates_s2, 
        #                     y = vals_s2, 
        #                     xout = sim_days, 
        #                     rule = 2)$y
        # dpai_adjusted <- 2 * daily_lai[168] # 2021-06-17
        # daily_lai <- 1
        
        # stop()
        # Store metrics in a list
        # metrics_list[[paste0("plot_", n_field[i], "_radius_", 
        #                      radius, "_clumping_", clumping)]] <- data.frame(
        #                        id_plot = n_field[i],
        #                        radius = radius,
        #                        clumping = round(clumping, 2),
        #                        lcv = round(lcv, 2),
        #                        lskew = round(lskew, 2),
        #                        vci = round(vci, 2),
        #                        gf = round(gf, 2),
        #                        mean = round(mean, 2),
        #                        max = round(canopy_height_top, 2),
        #                        pai = round(plantareaindex, 2)
        #                      )
        
        phenology <- calc_phenology(list.year = c(2021, 2022),
                                    nleafage = 1, # ---> 1
                                    budburst_date = 115,
                                    leaf_age_max_in = 0.56,
                                    relative_age_firstmax = 0.10,
                                    relative_age_lastmax = 0.75,
                                    LAI_max_per_cohort = plantareaindex)
        
        lais <- seq(1, 8, by = 1)
        for (daily_lai in lais) {
          phenology$Leaf_area_1yr[which(phenology$year == 2021)] <- 2 * daily_lai
          plantareaindex_max_s2 <- 2 * daily_lai
          
          
          # phenology$Leaf_area_1yr[which(phenology$year == 2021)] <- 2 * daily_lai
          # # plantareaindex_max_s2 <- dpai_adjusted
          # plantareaindex_max_s2 <- 2 * daily_lai
          # plantareaindex_max_s2 <- plantareaindex
          
          # library(ggplot2)
          # library(dplyr)
          # 
          # # Filter data for 2021 to check the injected profile
          # pheno_2021 <- phenology %>%
          #   filter(year == 2021) # %>%
          #   # filter(Julian_day >= 151 & Julian_day <= 274)
          # 
          # # Plot
          # p <- ggplot(pheno_2021, aes(x = Julian_day, y = Leaf_area_1yr)) +
          #   geom_line(color = "#2c7bb6", size = 1.2) + # Blue line
          #   geom_point(size = 0.8, alpha = 0.5) +      # Points to see daily steps
          #   labs(title = "MuSICA Input Phenology (June-Sept 2021)",
          #        subtitle = "Sentinel-2 Corrected TS (atbd)",
          #        # subtitle = "Sentinel-2 Corrected TS (optim common)",
          #        # subtitle = "Same LiDAR LAI over the period",
          #        x = "Julian Day (DOY)",
          #        y = "Leaf Area Index (m²/m²)") +
          #   theme_minimal() +
          #   # ylim(c(3, 7)) +
          #   theme(plot.title = element_text(face = "bold"))
          # plot(p)
          # 
          # stop()
          # musica_out <- callmusica(
          #   crown_clumping_factor = clumping,
          #   microclimate_height = microclimate_height,
          #   leaf_param = list("musica_veg1" =
          #                       list(phenology = phenology,
          #                            allometry = allometry,
          #                            PAI = plantareaindex,
          #                            canopy_height_top = canopy_height_top,
          #                            canopy_height_bottom = 2
          #                       )),
          #   run.name = paste0("pt_", n_field[i], "_radius_", radius),
          #   musica.cmd = musica.cmd,
          #   forcing_height = canopy_height_top + 2,
          #   save.output = TRUE,
          #   keep.tmp = FALSE)
          
          # musica_out <- callmusica(
          #   musica.param = list(
          #     "setupctl" = list(
          #       "clumping_factor" = clumping,
          #       "forcing_filename" = "./in_files/musica_in_Safran_Blois_2021.nc",
          #       "forcing_height" = canopy_height_top + 2
          #     )
          #   ),
          #   leaf.param = list("musica_veg1" =
          #                       list(phenology = phenology,
          #                            allometry = allometry,
          #                            PAI = PAI,
          #                            canopy_height_top = canopy_height_top,
          #                            canopy_height_bottom = 2
          #                            )),
          #   run.name = paste0("pt_", n_field[i], "_radius_", radius),
          #   musica.cmd = musica.cmd,
          #   out.df = TRUE,
          #   keep.tmp = FALSE
          # )
          
          # 
          # # File copy
          # dir.create(paste0("out_files/radius_test/", radius, "m"),
          #            recursive = TRUE, showWarnings = FALSE)
          # file_path <- paste0("out_files/musica_out_Blois_pt_", n_field[i],
          #                     "_radius_", radius, ".nc")
          # new_file_path <- paste0("out_files/radius_test/", radius,
          #                         "m/musica_out_Blois_pt_", n_field[i],
          #                         "_radius_", radius, ".nc")
          # if (file.copy(file_path, new_file_path)) {
          #   file.remove(file_path)
          # }
          
          forcing_file <- if (forcing_type == "Safran") {
            "./in_files/musica_in_Safran_Blois_2021.nc"
          } else {
            "./in_files/musica_in_Blois.nc"
          }
          
          # tryCatch({
          #   forcing_time_utc <- force_utc_nc(forcing_file, varname = "time")
          #   
          #   cat("--------------------------------------------------\n")
          #   cat("--> Vérification Forçage :", basename(forcing_file), "\n")
          #   cat("    Unités détectées     :", ncatt_get(nc_open(forcing_file), 
          #                                               "time", "units")$value, "\n")
          #   cat("    Premier pas de temps :", format(forcing_time_utc[1], 
          #                                            "%Y-%m-%d %H:%M:%S %Z"), "\n")
          #   cat("    Dernier pas de temps :", format(tail(forcing_time_utc, 1),
          #                                            "%Y-%m-%d %H:%M:%S %Z"), "\n")
          #   cat("--------------------------------------------------\n")
          #   tmp_nc <- nc_open(forcing_file); nc_close(tmp_nc) 
          #   
          # }, error = function(e) {
          #   cat("Erreur lors de la lecture du temps NetCDF :", e$message, "\n")
          #   stop("Arrêt critique : Problème de temps dans le fichier forçage.")
          # })
          # stop()
          history_file <- paste0("./out_files/TS/",
                                 "smooth",
                                 # "/RF/LAI_ALS/",
                                 # "/corr/",
                                 # "/gap/",
                                 "/lidar/",
                                 # distrib,
                                 paste0("lai_", daily_lai),
                                 "/musica_out_", forcing_type,
                                 "_Blois_2021_pt_", n_field[i], ".nc")
          
          if (!dir.exists(dirname(history_file))) {
            dir.create(dirname(history_file), showWarnings = FALSE, recursive = TRUE)
          }
          
          musica_out <- callmusica(
            musica.param = list(
              "setupctl" = list(
                "clumping_factor" = clumping,
                "forcing_filename" = forcing_file,
                "history_filename" = history_file,
                "forcing_height" = canopy_height_top + 2
              )
            ),
            leaf.param = 
              list( # species list
                "musica_veg1" =
                  list( # phenology, allometry and namelist control
                    "phenology" = phenology,
                    "allometry" = allometry,
                    "leafphenologyctl" = 
                      list(
                        "lai_max_per_cohort" = plantareaindex_max_s2
                      ),
                    "leafallometryctl" =
                      list(
                        "canopy_height_top" = canopy_height_top,
                        "canopy_height_bottom" = 2
                      ),
                    "leafmusicactl" =
                      list(
                        "canopy_height_top" = canopy_height_top
                      )
                  )
              ),
            # run.name = "run1",
            musica.cmd = musica.cmd,
            keep.tmp = T,
            out.netcdf = T
            # out.var = c("Tair_z","wair_z", "T_soil", "w_soil", "h_canopy"),
            # out.subset = list("zair" = c(1,3,10),
            #                   "zsoil" = c(0.1, 0.5, 1),
            #                   "nspecies" = "all",
            #                   "nleafage" = "all",
            #                   "time" = "all",
            #                   "Tair_z" = list("zair" = c(1,3)))
          )
          # stop()
        } # lais loops
      }
    }
  }
}
# metrics_df <- do.call(rbind, metrics_list)
# 
# file_path <- "out_files/radius_test/metrics_results.csv"
# if (!file.exists(file_path)) {
#   write.csv(metrics_df, file_path, row.names = FALSE)
# } else {
#   write.table(metrics_df, file_path, sep = ",", row.names = FALSE, 
#               col.names = FALSE, append = TRUE)
# }
# stop()
# metrics_df <- read.csv("out_files/radius_test/metrics_results.csv", 
#                        header = TRUE)
# stop()

# --------------------------------- Air ----------------------------------------
# for (n_field in seq_len(dim(df_coords)[1])){
# df_Tair_all <- data.frame()
# radius <- 5
# clumping <- 0.7
# date_study <- "2021-06-20"
# # Loop over each field
# for (i in 1:nrow(df_lad)) {
#   id_plot <- n_field[i]
#   
#   # Open the NetCDF file for this field
#   # nc_file <- paste0("out_files/musica_out_Blois_pt_", id_plot, ".nc")
#   nc_file <- paste0("out_files/radius_test/", radius, "m/clumping_", clumping,
#                     "/musica_out_Blois_pt_", id_plot,
#                     "_radius_", radius, "_clumping_", clumping, ".nc")
#   out <- nc_open(nc_file)
#   out.all[[paste0("pt_", id_plot)]] <- out
#   
#   # Extract the variable
#   df_Tair <- get_variable(out, "Tair_z") %>%
#     filter(as.Date(time) == date_study, nair == 1) %>%
#     mutate(Tair_z = Tair_z - 273.15,
#            id_plot = id_plot)
#   
#   ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
#   df_Tair_all <- bind_rows(df_Tair_all, df_Tair)
# }
# # shinymusica(out.all)
# df_Tair_all <- df_Tair_all %>%
#   mutate(time = floor_date(time, unit = "hour"))
# 
# temperature_all <- data.frame()
# # for (site in sites) {
# # Read the CSV file for the site
# site <- "Blois"
# csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
# temperature_csv <- read.csv(csv_file)
# 
# # Convert the datetime column to POSIXct
# temperature_df <- temperature_csv %>%
#   mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", 
#                                tz = "UTC"))
# 
# # Filter for the desired date and sensor position (e.g., "a")
# temperature_df <- temperature_df %>%
#   filter(as.Date(datetime) == as.Date(date_study),
#          id_plot %in% df_Tair_all$id_plot,
#          position_sensor == "a"
#   )
# 
# missing_in_temperature <- setdiff(df_Tair_all$id_plot, temperature_df$id_plot)
# print(missing_in_temperature)
# 
# # Find common IDs
# common_ids <- intersect(df_Tair_all$id_plot, temperature_df$id_plot)
# 
# # Filter both data frames to keep only the common IDs
# df_Tair_all <- df_Tair_all %>% filter(id_plot %in% common_ids)
# temperature_df <- temperature_df %>% filter(id_plot %in% common_ids)
# 
# # Filter and rename
# df_Tair_all <- df_Tair_all %>%
#   dplyr::select(id_plot, datetime = time, Tair_z)
# temperature_df <- temperature_df %>%
#   dplyr::select(id_plot, datetime, Tair_z = t_hobo)
# 
# # Analysis
# remove_id_plot <- c("41_28", "41_48", "41_36", "41_40", "41_42", "41_56", "41_59")
# df_joined <- left_join(df_Tair_all, temperature_df,
#                        by = c("id_plot", "datetime"),
#                        suffix = c("_sim", "_obs"))
# df_joined <- df_joined %>%
#   filter(!id_plot %in% remove_id_plot)
# df_joined <- df_joined %>% drop_na()
# 
# # Model
# model <- lm(Tair_z_obs ~ Tair_z_sim, data = df_joined)
# model_summary <- summary(model)
# 
# # Metrics
# slope <- coef(model)[2]
# intercept <- coef(model)[1]
# r_value <- cor(df_joined$Tair_z_sim, df_joined$Tair_z_obs, use = "complete.obs")
# r2 <- model_summary$r.squared
# rmse <- sqrt(mean((df_joined$Tair_z_obs - df_joined$Tair_z_sim)^2, 
#                   na.rm = TRUE))
# obs_iqr <- IQR(df_joined$Tair_z_obs, na.rm = TRUE)
# nrmse <- rmse / obs_iqr
# bias <- mean(df_joined$Tair_z_sim - df_joined$Tair_z_obs, na.rm = TRUE)
# 
# # Plot
# ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
#   geom_point(aes(color = id_plot), size = 2, alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   labs(x = "Observed Tair_z (°C)",
#        y = "Simulated Tair_z (°C)",
#        title = paste("Observed vs Simulated Tair_z in", site, "for date", date),
#        subtitle = paste("R =", round(r_value, 2),
#                         "| R² =", round(r2, 2),
#                         "| NRMSE =", round(nrmse, 2),
#                         "| Bias =", round(bias, 2),
#                         "| Slope =", round(slope, 2),
#                         "| Intercept =", round(intercept, 2))) +
#   # xlim(17, 21) +
#   # ylim(17, 21) +
#   theme_bw()
# 
# # Append to the overall CSV data frame
# # temperature_all <- bind_rows(temperature_all, temperature_df)
# # }
# 
# # ---------------------------id_plot facet -------------------------------------
# ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
#   geom_point(size = 2, alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   facet_wrap(~ id_plot, ncol = 4) +  # Adjust as needed
#   labs(x = "Observed Tair_z (°C)",
#        y = "Simulated Tair_z (°C)",
#        title = "Scatterplot: Observed vs. Simulated Tair_z by Field") +
#   theme_bw()
# 
# # -------------------------- datetime facet ------------------------------------
# ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim)) +
#   geom_point(size = 2, alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   facet_wrap(~ datetime, scales = "free_x", ncol = 4) +  # Facet by full datetime
#   labs(x = "Observed Tair_z (°C)",
#        y = "Simulated Tair_z (°C)",
#        title = "Scatterplot: Observed vs. Simulated Tair_z by Time") +
#   theme_bw()
# 
# # --------------------------- Color hours --------------------------------------
# df_joined <- df_joined %>%
#   mutate(hour = hour(datetime))
# ggplot(df_joined, aes(x = Tair_z_obs, y = Tair_z_sim, color = factor(hour))) +
#   geom_point(size = 2, alpha = 0.8) +
#   geom_smooth(method = "lm", se = FALSE, color = "black") +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#   facet_wrap(~ id_plot, ncol = 4) +  # One scatterplot per `id_plot`
#   scale_color_viridis_d() +  # Better color scheme for hours
#   labs(x = "Observed Tair_z (°C)",
#        y = "Simulated Tair_z (°C)",
#        title = "Scatterplot: Observed vs. Simulated Tair_z",
#        subtitle = "Colored by Hour of the Day",
#        color = "Hour") +
#   theme_bw()
# 
# # --------------------------- id_plot loop -------------------------------------
# plot_list <- list()
# for (plot_id in unique(df_joined$id_plot)) {
#   p <- ggplot(df_joined %>% filter(id_plot == plot_id), 
#               aes(x = Tair_z_obs, y = Tair_z_sim, color = factor(hour))) +
#     geom_point(size = 2, alpha = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, color = "black") +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
#     scale_color_viridis_d() +  
#     labs(x = "Observed Tair_z (°C)",
#          y = "Simulated Tair_z (°C)",
#          title = paste("Scatterplot for id_plot:", plot_id),
#          subtitle = "Colored by Hour of the Day",
#          color = "Hour") +
#     theme_bw()
#   
#   # Store the plot in the list
#   plot_list[[plot_id]] <- p
# }
# 
# # Print all plots
# for (p in plot_list) {
#   print(p)
# }
# 
# # ------------------------------- Monthly --------------------------------------
# # Parameters
# radius <- 5
# clumping <- 0.7
# remove_id_plot <- c("41_28", "41_48", "41_36", "41_40", "41_42", "41_56", "41_59")
# # Define the sequence of dates in June 2021
# dates_june <- seq.Date(as.Date("2021-06-01"), as.Date("2021-06-30"), by = "day")
# 
# # Initialize an empty data frame to store simulated temperature data
# df_Tair_all <- data.frame()
# 
# # Loop over each day in June
# for (current_date in dates_june) {
#   current_date <- as.Date(current_date, origin = "1970-01-01")
#   cat("Processing date:", as.character(current_date), "\n")
#   
#   # Loop over each field (assuming df_lad and n_field are defined)
#   for (i in 1:nrow(df_lad)) {
#     id_plot <- n_field[i]
#     
#     # Construct the path to the NetCDF file for this field
#     nc_file <- paste0("out_files/radius_test/", radius, "m/clumping_", clumping,
#                       "/musica_out_Blois_pt_", id_plot,
#                       "_radius_", radius, "_clumping_", clumping, ".nc")
#     
#     # Open the NetCDF file and extract the "Tair_z" variable
#     nc <- nc_open(nc_file)
#     # (Optional) store the nc object if needed:
#     # out.all[[paste0("pt_", id_plot)]] <- nc
#     
#     df_Tair <- get_variable(nc, "Tair_z") %>%
#       filter(as.Date(time) == current_date, nair == 1) %>%  # filter for the current date and sensor condition
#       mutate(Tair_z = Tair_z - 273.15,       # convert from Kelvin to Celsius
#              id_plot = id_plot,
#              date = current_date)            # add the current date as a new column
#     
#     df_Tair_all <- bind_rows(df_Tair_all, df_Tair)
#   }
# }
# 
# # Round timestamps to the nearest hour and extract the hour of day
# df_Tair_all <- df_Tair_all %>%
#   mutate(datetime_hour = floor_date(time, unit = "hour"),
#          hour = hour(time)) %>%
#   filter(!id_plot %in% remove_id_plot)
# 
# # Aggregate simulated temperatures: mean per hour per day
# df_simulated <- df_Tair_all %>%
#   group_by(date, hour) %>%
#   summarise(Tair_z_mean = mean(Tair_z, na.rm = TRUE), .groups = "drop")
# 
# # --- Load Observed Temperature Data ---
# csv_file <- file.path("in_files", paste0(site, "_data_temperature.csv"))
# temperature_csv <- read.csv(csv_file)
# 
# temperature_df <- temperature_csv %>%
#   mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
#   filter(as.Date(datetime) %in% dates_june, position_sensor == "a") %>%
#   # filter(as.Date(datetime) %in% "2021-06-17", position_sensor == "a") %>%
#   mutate(date = as.Date(datetime), hour = hour(datetime))
# 
# # Aggregate observed temperatures: mean per hour over ALL plots
# df_observed_mean <- temperature_df %>%
#   group_by(hour) %>%
#   summarise(Tair_z_obs_mean = mean(t_hobo, na.rm = TRUE), .groups = "drop")
# 
# # --- Plot: Simulated Curves + Observed Mean Curve ---
# ggplot() +
#   geom_line(data = df_simulated, aes(x = hour, y = Tair_z_mean, color = factor(date)), size = 1, alpha = 0.7) +
#   geom_line(data = df_observed_mean, aes(x = hour, y = Tair_z_obs_mean), color = "black", size = 1.2, linetype = "solid") +
#   labs(x = "Hour of Day", y = "Temperature (°C)",
#        title = "Simulated vs Observed Temperature (June)",
#        subtitle = "Black curve = mean observed temperature across all plots",
#        color = "Date") +
#   theme_minimal()
# 
# # Compute absolute differences between simulated and observed temperatures
# df_comparison <- df_simulated %>%
#   left_join(df_observed_mean, by = "hour") %>%
#   mutate(abs_diff = abs(Tair_z_mean - Tair_z_obs_mean))
# 
# # Sum differences per simulated date
# df_total_diff <- df_comparison %>%
#   group_by(date) %>%
#   summarise(total_diff = sum(abs_diff, na.rm = TRUE), .groups = "drop")
# 
# # Find the simulated date with the smallest total difference
# best_match_date <- df_total_diff %>%
#   filter(total_diff == min(total_diff)) %>%
#   pull(date)
# 
# best_match_date <- as.Date(as.numeric(best_match_date), origin = "1970-01-01")
# print(paste("The simulated date closest to observed temperatures is:", best_match_date))
# 
# # ----------------------- LiDAR Metrics Insights ------------------------------
# radius <- 5
# clumping <- 0.7
# metrics_df_filtered <- metrics_df %>%
#   filter(radius == 5, clumping == 0.5)
# 
# # Compute Per-Plot Simulation Errors
# 
# error_metrics <- df_joined %>%
#   group_by(id_plot) %>%
#   summarise(
#     mean_bias = mean(Tair_z_sim - Tair_z_obs, na.rm = TRUE),
#     nrmse = sqrt(mean((Tair_z_sim - Tair_z_obs)^2, na.rm = TRUE)) / IQR(Tair_z_obs, na.rm = TRUE),
#     r = cor(Tair_z_sim, Tair_z_obs, use = "complete.obs"),
#     n = n(),
#     .groups = "drop"
#   )
# 
# # Merge Error Metrics with Other Metrics
# error_with_metrics <- error_metrics %>%
#   inner_join(metrics_df_filtered, by = "id_plot")
# 
# # Visualize Relationships Between Simulation Errors and Metrics
# metrics <- c("lcv", "lskew", "vci", "gf", "mean")
# 
# # Function to create scatter plots for RMSE and Bias
# plot_relationships <- function(metric) {
#   
#   p_r <- ggplot(error_with_metrics, aes_string(x = metric, y = "r", color = "id_plot")) +
#     geom_point(size = 3, alpha = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, color = "red") +
#     labs(x = metric,
#          y = "Pearson Correlation",
#          title = paste("Relationship between", metric, "and Pearson Correlation"),
#          subtitle = "Each point represents one plot") +
#     theme_bw() +
#     scale_color_discrete(name = "Plot ID")
#   
#   p_nrmse <- ggplot(error_with_metrics, aes_string(x = metric, y = "nrmse", color = "id_plot")) +
#     geom_point(size = 3, alpha = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, color = "red") +
#     labs(x = metric,
#          y = "NRMSE (°C)",
#          title = paste("Relationship between", metric, "and NRMSE"),
#          subtitle = "Each point represents one plot") +
#     theme_bw() +
#     scale_color_discrete(name = "Plot ID")
#   
#   p_bias <- ggplot(error_with_metrics, aes_string(x = metric, y = "mean_bias", color = "id_plot")) +
#     geom_point(size = 3, alpha = 0.8) +
#     geom_smooth(method = "lm", se = FALSE, color = "red") +
#     labs(x = metric,
#          y = "Mean Bias (°C)",
#          title = paste("Relationship between", metric, "and Simulation Bias"),
#          subtitle = "Each point represents one plot") +
#     theme_bw() +
#     scale_color_discrete(name = "Plot ID")
#   
#   list(p_r = p_r, p_nrmse = p_nrmse, p_bias = p_bias)
# }
# 
# # Generate plots for all metrics
# plot_list <- map(metrics, plot_relationships)
# walk(plot_list, ~ print(.x$p_r))
# walk(plot_list, ~ print(.x$p_nrmse))
# walk(plot_list, ~ print(.x$p_bias))
# 
# # Visualizing Simulation Temperature Curves Colored by a Metric 
# df_sim_with_metrics <- df_Tair_all %>%
#   inner_join(metrics_df_filtered, by = "id_plot")
# 
# # For example, plot simulation curves for a given date (or for all data) with color indicating the 'mean' metric:
# ggplot(df_sim_with_metrics, aes(x = datetime, y = Tair_z, group = id_plot, color = mean)) +
#   geom_line(alpha = 0.8) +
#   scale_color_gradient(low = "blue", high = "red") +
#   labs(x = "Datetime",
#        y = "Simulated Temperature (°C)",
#        title = "Simulation Temperature Curves Colored by Mean Metric",
#        color = "Mean") +
#   theme_minimal()
