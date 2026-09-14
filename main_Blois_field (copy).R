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

# Initialize an empty list to store the results for each file
results <- list()
allometry_list <- list()
musica_list <- list()
out.all <- list()
microclimate_height <- 1
inter_crown_clumping <- 0.7 # Loop w/ =/= CCI ?
musica.cmd <- "bash -i -c musica"

LAI_stack <- terra::rast("in_files/ladstack.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson")

df_coords <- json %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(coord_x_utm31n = X, coord_y_utm31n = Y)

height_levels <- seq(2.5, 2.5 + (nlyr(LAI_stack) - 1), by = 1)
lad_values <- terra::extract(LAI_stack, df_coords)
df_lad <- as.data.frame(lad_values)
df_lad$ID <- json$id_plot
df_lad <- df_lad %>%
  filter(rowSums(!is.na(dplyr::select(., -ID))) > 0)

n_field <- df_lad$ID
df_lad$ID = NULL
df_lad[is.na(df_lad)] <- 0.0

for (i in 1:nrow(df_lad)) {
  lad_values <- as.numeric(df_lad[i, ])
  allometry <- data.frame(
    height = height_levels,
    density = lad_values
  )
  canopy_height_top <- min(40, max(allometry$height[allometry$density > 0] + 0.5, na.rm = TRUE))

  PAI <- 2 * sum(allometry$density, na.rm = TRUE)
  
  phenology <- calc_phenology(list.year = c(2021),
                              nleafage = 3,
                              budburst_date = 115,
                              leaf_age_max_in = 0.56,
                              relative_age_firstmax = 0.10,
                              relative_age_lastmax = 0.75,
                              LAI_max_per_cohort = PAI/3)
  
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
    run.name = paste0("pt_", n_field[i]),
    musica.cmd = musica.cmd,
    save.output = TRUE,
    keep.tmp = FALSE)
}

# Tair == 1 (0-1.3m,field = 1m), Tsoil == 4 (68.9-102.7mm, field = 80mm)

out.all <- list()
# for (n_field in seq_len(dim(df_coords)[1])){
for (i in n_field) {
  out <- nc_open(paste0("out_files/musica_out_Blois_pt_", i, ".nc"))
  out.all[[paste0("pt_", i)]] <- out
  df_Tair <- get_variable(out, "Tair_z")
  df_Tair$Tair_z <- df_Tair$Tair_z - 273.15
  names(attributes(df_Tair))
  attr(df_Tair, "relative_height")
  ggplot_variable(filter(df_Tair, nair == 1), out.type = "standard")
  a
}




shinymusica(out.all)

# Tools
t1 <- nc_open("out_files/musica_out_Blois_pt1.nc")
t2 <- nc_open("out_files/musica_out_Blois_pt2.nc")
t22 <- nc_open("out_files/musica_out_Blois_pt22.nc")
t47 <- nc_open("out_files/musica_out_Blois_pt47.nc")

df_Tair <- get_variable(t47, "Tair_z")
names(attributes(df_Tair))
attr(df_Tair, "relative_height")
ggplot_variable(filter(df_Tair, nair == 5), out.type = "standard")

filter(df_Tair, nair == 5) %>%
  select(-nair) %>%
  ggplot_variable(out.type = "daily_heatmap")

ggplot_variable(df_Tair, out.type = "heatmap")

# pt1 <- clip_rectangle(ctg,
#                       xleft = df_coords[1,]$coord_x_l93 - 5,
#                       ybottom = df_coords[1,]$coord_y_l93 - 5,
#                       xright = df_coords[1,]$coord_x_l93 + 5,
#                       ytop = df_coords[1,]$coord_y_l93 + 5
# )
# dtm <- grid_terrain(pt1, 1, kriging(k = 10L))
# pt1 <- normalize_height(pt1, dtm)
# chm <- grid_canopy(pt1, res = 1, dsmtin())
#
# lad <- LAD(pt1@data$Z)
# allometry <- data.frame(
#   height = lad$z,
#   density = lad$lad
# )
#
# canopy_top <- cellStats(chm, "max")
# plant_area_index <- 2 * sum(allometry$density)
#
# musica_out <- callmusica_vegetation_structure(
#   allometry = allometry,
#   clumping_factor = inter_crown_clumping,
#   PAI = plant_area_index,
#   canopy_height_top = canopy_top,
#   microclimate_height = microclimate_height,
#   species.nml = "in_files/musica_veg1.nml",
#   run.name = "pt1",
#   musica.cmd = musica.cmd,
#   save.output = TRUE)

# out <- nc_open("out_files/musica_out_Blois.nc")
# out.all <- list("pt1" = out
# )
#
# shinymusica(out.all)

# Loop through each file
# for (file in files) {
# 
#   run_name <- tools::file_path_sans_ext(basename(file))
#   las <- readLAS(file, filter = "-keep_first")
#   if (is.null(las)) next
# 
#   dtm <- grid_terrain(las, 1, kriging(k = 10L))
#   las_norm <- normalize_height(las, dtm)
#   chm <- grid_canopy(las_norm, res = 1, dsmtin())
#   print(chm)
# 
#   lad <- LAD(las_norm@data$Z)
#   allometry <- data.frame(
#     height = lad$z,
#     density = lad$lad
#   )
# 
#   canopy_top <- cellStats(chm, "max")
#   plant_area_index <- 2 * sum(allometry$density)
# 
#   results[[run_name]] <- list(
#     canopy_top = canopy_top,
#     plant_area_index = plant_area_index,
#     inter_crown_clumping = inter_crown_clumping,
#     microclimate_height = microclimate_height,
#     allometry = allometry
#   )
# 
#   # Loop MuSICA
#   musica_out <- callmusica_vegetation_structure(
#     allometry = allometry,
#     clumping_factor = inter_crown_clumping,
#     PAI = plant_area_index,
#     canopy_height_top = canopy_top,
#     microclimate_height = microclimate_height,
#     species.nml = "in_files/musica_veg_FR-Bil_pinus_pinaster.nml",
#     run.name = run_name,
#     musica.cmd = musica.cmd,
#     save.output = TRUE)
# 
#   musica_list[[run_name]] <- musica_out
#   allometry$run <- run_name
#   allometry_list[[run_name]] <- allometry
# }
# 
# combined_allometry <- do.call(rbind, allometry_list)
# ggplot(combined_allometry) +
#   geom_path(aes(x = density, y = height, color = run))
# 
# out_conv1 <- nc_open("out_files/musica_out_FR-Bil_2022_CONV1_15m.nc")
# out_conv2 <- nc_open("out_files/musica_out_FR-Bil_2022_CONV2_15m.nc")
# out_faible1 <- nc_open("out_files/musica_out_FR-Bil_2022_FAIBLE1_15m.nc")
# out_faible2 <- nc_open("out_files/musica_out_FR-Bil_2022_FAIBLE2_15m.nc")
# out_vieux1 <- nc_open("out_files/musica_out_FR-Bil_2022_VIEUX1_15m.nc")
# out_vieux2 <- nc_open("out_files/musica_out_FR-Bil_2022_VIEUX2_15m.nc")
# 
# out.all <- list("out_conv1" = out_conv1,
#                 "out_conv2" = out_conv2,
#                 "out_faible1" = out_faible1,
#                 "out_faible2" = out_faible2,
#                 "out_vieux1" = out_vieux1,
#                 "out_vieux2" = out_vieux2
# )
# 
# shinymusica(out.all)
# 
# 
# 
# conv1_allometry <- subset(combined_allometry, run == "VIEUX2_15m")
# conv1_allometry$run <- NULL
# musica_out1 <- callmusica_vegetation_structure(
#   allometry = conv1_allometry,
#   clumping_factor = inter_crown_clumping,
#   PAI = plant_area_index,
#   canopy_height_top = canopy_top,
#   microclimate_height = microclimate_height,
#   species.nml = "in_files/musica_veg_FR-Bil_pinus_pinaster.nml",
#   run.name = "VIEUX2_15m",
#   musica.cmd = musica.cmd,
#   save.output = TRUE)
# 
# out <- nc_open("/home/corroyez/Documents/visite_bdx_nov24/Bilos/out_files/musica_out_FR-Bil_2022_VIEUX2_15m.nc")
# df_Tair <- get_variable(out, "Tair_z")
# shinymusica(list('run' = out))







# las <- readLAS("../Lidar_Louchats/MAST/CONV1_15m.las", filter = "-keep_first")
# dtm <- grid_terrain(las, 1, kriging(k = 10L))
# las_norm <- normalize_height(las, dtm)
# chm <- grid_canopy(las_norm, res = 1, dsmtin())
# lad <- LAD(las_norm@data$Z)
# 
# microclimate_height <- 1.5
# allometry <- data.frame(
#   height = lad$z,
#   density = lad$lad
# )
# allometry$run <- "CONV1"
# # canopy_top <- cellStats(chm, "max") # meters above ground
# canopy_top <- 20
# inter_crown_clumping <- 0.8 # unitless (between 0 and 1)
# plant_area_index <- 2*sum(allometry$density) # m2/m2
# 
# musica.cmd <- "bash -i -c musica"
# musica_out <- callmusica_vegetation_structure(
#   allometry = allometry,
#   clumping_factor = inter_crown_clumping,
#   PAI = plant_area_index,
#   canopy_height_top = canopy_top,
#   microclimate_height = microclimate_height,
#   species.nml = "in_files/musica_veg_FR-Bil_pinus_pinaster.nml",
#   run.name = "CONV1",
#   musica.cmd = musica.cmd,
#   save.output = TRUE)
