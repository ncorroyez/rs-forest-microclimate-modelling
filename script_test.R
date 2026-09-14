rm(list=ls(all=TRUE)) # Clear the global environment (remove all objects)
gc() # Trigger the garbage collector to free up memory
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  getwd()
}

######################################
## Define call function
library(musica.tools)
library(rmusica)
library(ncdf4)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
######################################@
## MAIN PROGRAM
######################################@
musica.cmd <- "bash -i -c musica"

## targeted microclimate height
# no more used ???
microclimate_height <- 1.5 # meters above ground

## compute microclimate for canopy structure data from site_loc #1
## note: allometry$height should have at least one value > canopy_top


# Setup parameters --------------------------------------------------------
allometry1 <- data.frame(
  height = seq(2, 22, by = 2),
  density = c(0.1, 0.1, 0.1, 0.1, 0.2, # 2-10
              0.3, 0.5, 0.5, 0.1, 0.0, # 12-20
              0.0) # 22
)
inter_crown_clumping1 <- 0.8 # unitless (between 0 and 1)
plant_area_index1 <- 7.0 # m2/m2


# plot allometry ----------------------------------------------------------

tmp1 <- allometry1
tmp1$run <- "run1"

# setup phenology ---------------------------------------------------------

phenology1 <- calc_phenology(list.year = c(2021, 2022),
                             nleafage = 3,
                             budburst_date = 110,
                             leaf_age_max_in = 2.8,
                             relative_age_firstmax = 0.05,
                             relative_age_lastmax = 0.83,
                             LAI_max_per_cohort = plant_area_index1/3)

phenology1 <- calc_phenology(list.year = c(2021, 2022),
                             nleafage = 3, # ---> 1
                             budburst_date = 115,
                             leaf_age_max_in = 0.56,
                             relative_age_firstmax = 0.10,
                             relative_age_lastmax = 0.75,
                             LAI_max_per_cohort = plant_area_index1/3)
# Run 1 -------------------------------------------------------------------

musica_out <- callmusica(
  musica.param = list(
    "setupctl" = list(
      "clumping_factor" = inter_crown_clumping1,
      "forcing_height" = max(allometry1$height) + 2
    )
  ),
  leaf.param = 
    list(
      "musica_veg1" =
        list(
          "phenology" = phenology1,
          "allometry" = allometry1,
          "leafphenologyctl" = 
            list(
              "lai_max_per_cohort" = plant_area_index1/3
            ),
          "leafallometryctl" = 
            list(
              "canopy_height_top" = max(allometry1$height),
              "canopy_height_bottom" = min(allometry1$height) # =2
            ),
          "leafmusicactl" = 
            list(
              "canopy_height_top" = max(allometry1$height)
            )
        )
    ),
  run.name = "run1",
  musica.cmd = musica.cmd,
  keep.tmp = TRUE,
  out.netcdf = F,
  # out.df = F,
  out.var = c("Tair_z","wair_z", "T_soil", "w_soil", "h_canopy"),
  out.subset = list("zair" = c(1,3,10),
                    "zsoil" = c(0.1, 0.5, 1),
                    "nspecies" = "all",
                    "nleafage" = "all",
                    "time" = "all",
                    "Tair_z" = list("zair" = c(1,3)))
)
