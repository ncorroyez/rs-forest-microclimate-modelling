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


allometry2 <- data.frame(
  height = seq(2, 22, by = 2),
  density = c(0.1, 0.3, 0.4, 0.4, 0.3, # 2-10
              0.3, 0.2, 0.1, 0.0, 0.0, # 12-20
              0.0) # 22
)
inter_crown_clumping2 <- 0.3 # unitless (between 0 and 1)
plant_area_index2 <- 7.0 # m2/m2


# plot allometry ----------------------------------------------------------

tmp1 <- allometry1
tmp1$run <- "run1"
tmp2 <- allometry2
tmp2$run <- "run2"
allometry <- rbind(tmp1, tmp2)
rm(tmp1, tmp2)

ggplot(allometry) +
  geom_path(aes(x = density, y = height, color = run))

# setup phenology ---------------------------------------------------------

phenology1 <- calc_phenology(list.year = c(2022),
                             nleafage = 3,
                             budburst_date = 110,
                             leaf_age_max_in = 2.8,
                             relative_age_firstmax = 0.05,
                             relative_age_lastmax = 0.83,
                             LAI_max_per_cohort = plant_area_index1/3)

phenology2 <- calc_phenology(list.year = c(2022),
                             nleafage = 3,
                             budburst_date = 120,
                             leaf_age_max_in = 2.8,
                             relative_age_firstmax = 0.12,
                             relative_age_lastmax = 0.90,
                             LAI_max_per_cohort = plant_area_index2/3)


# plot phenology ----------------------------------------------------------


tmp1 <-
  phenology1 %>%
  format_phenology() %>%
  mutate(run = "run1")
tmp2 <-
  phenology2 %>%
  format_phenology() %>%
  mutate(run = "run2")
phenology <-
  rbind(tmp1, tmp2)
phenology %>%
  ggplot() +
  geom_line(aes(x = relative_age, y = LAI, color = run))

phenology %>%
  ggplot() +
  geom_line(aes(x = Julian_day, y = LAI, color = run)) +
  facet_wrap(~leaf_age, ncol = 1)


# Run 1 -------------------------------------------------------------------


musica_out1 <- callmusica_vegetation_structure(
  crown_clumping_factor = inter_crown_clumping1,
  microclimate_height = microclimate_height,
  leaf_param = list("musica_veg_FR-Bil_pinus_pinaster" =
                      list(phenology = phenology1,
                           allometry = allometry1,
                           PAI = plant_area_index1)),
  run.name = "run1",
  musica.cmd = musica.cmd,
  save.output = TRUE)


# run 2 -------------------------------------------------------------------
## compute microclimate for canopy structure data from site_loc #2
## note: allometry$height should have at least one value > canopy_top
## site_loc <- 2
##

musica_out2 <- callmusica_vegetation_structure(
  crown_clumping_factor = inter_crown_clumping2,
  microclimate_height = microclimate_height,
  leaf_param = list("musica_veg_FR-Bil_pinus_pinaster" =
                      list(phenology = phenology2,
                           allometry = allometry2,
                           PAI = plant_area_index2)),
  run.name = "run2",
  musica.cmd = musica.cmd,
  save.output = TRUE)



musica_out1$run <- "run1"
musica_out2$run <- "run2"

df.out <- rbind(musica_out1, musica_out2)

ggplot(df.out) +
  geom_line(aes(x = time, y = Tair_1.5), color = "#1b9e77") +
  geom_line(aes(x = time, y = T_soil_0.06), color = "#d95f02") +
  facet_wrap(~run) +
  ylab("Temperature (K)")

df.out %>%
  pivot_wider(id_cols = "time", names_from = "run",
              values_from = c("Tair_1.5", "wair_1.5", "wind_1.5",
                              "T_soil_0.06", "w_soil_0.06") ) %>%
  mutate(diffTair = Tair_1.5_run1 - Tair_1.5_run2) %>%
  ggplot() +
  geom_line(aes(x = time, y = diffTair)) +
  ylab("Air Temperature differences")
