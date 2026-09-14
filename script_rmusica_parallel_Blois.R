rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

library(terra)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)

# Data preparation
ctg <- readLAScatalog("/media/corroyez/My Passport/01_DATA/Blois/LiDAR/2-las_utm")
LAI_stack <- terra::rast("in_files/ladstack_classic.tif")
json <- st_read("in_files/data_Blois_utm31n.geojson")

height_levels <- seq(2.5, 2.5 + (nlyr(LAI_stack) - 1), by = 1)
LAI_df <- as.data.frame(LAI_stack, xy = TRUE, na.rm = FALSE) 
valid_pixels <- rowSums(is.na(LAI_df[ , -(1:2)])) < ncol(LAI_stack) 
LAI_df <- LAI_df[valid_pixels, ]
allometry_df <- LAI_df %>%
  pivot_longer(cols = -c(x, y), names_to = "layer", values_to = "density") %>%
  mutate(height = rep(height_levels, each = nrow(LAI_df))) %>%
  select(x, y, height, density)

registerDoParallel(cl <- makeCluster(3))
fullMuSICA <- foreach(PAI = 6:8, .combine = "rbind") %dopar% {
  ######################################@
  ## Define call function
  library(musica.tools)
  library(rmusica)
  library(ncdf4)
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(terra)
  ######################################@
  ## MAIN PROGRAM
  ######################################@
  musica.cmd <- "bash -i -c musica"

  ## targeted microclimate height
  microclimate_height <- 1 # meters above ground

  ## compute microclimate for canopy structure data from site_loc #1
  ## note: allometry$height should have at least one value > canopy_top


  # Setup parameters --------------------------------------------------------
  # allometry1 <- data.frame(
  #   height = seq(2, 22, by = 2),
  #   density = c(0.1, 0.1, 0.1, 0.1, 0.2, # 2-10
  #               0.3, 0.5, 0.5, 0.1, 0.0, # 12-20
  #               0.0) # 22
  # )
  # list of all non-NA pixels, with PAD and coordinates
  # one script for field points (and compare with field T°)
  # one script for map
  allometry1 <- data.frame(
    height = seq(2.5, 39.5, by = 1),
    density = ladstack$LAD
  )
  PAI <- sum(allometry1, na.rm = TRUE)
  inter_crown_clumping1 <- 0.8 # unitless (between 0 and 1)

  phenology1 <- calc_phenology(list.year = c(2021),
                               nleafage = 3,
                               budburst_date = 103,
                               leaf_age_max_in = 0.62,
                               relative_age_firstmax = 0.127,
                               relative_age_lastmax = 0.69,
                               LAI_max_per_cohort = PAI/3)

  musica_out <- callmusica_vegetation_structure(
    crown_clumping_factor = inter_crown_clumping1,
    microclimate_height = microclimate_height,
    leaf_param = list("in_files/musica_veg1.nml" =
                        list(phenology = phenology1,
                             allometry = allometry1,
                             PAI = PAI)),
    run.name = paste0(PAI),
    musica.cmd = musica.cmd,
    save.output = FALSE,
    keep.tmp = FALSE)
  musica_out$PAI <- PAI
  musica_out
}
