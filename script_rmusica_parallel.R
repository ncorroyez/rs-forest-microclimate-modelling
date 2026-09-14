library(foreach)
library(doParallel)

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

  phenology1 <- calc_phenology(list.year = c(2022),
                               nleafage = 3,
                               budburst_date = 110,
                               leaf_age_max_in = 2.8,
                               relative_age_firstmax = 0.05,
                               relative_age_lastmax = 0.83,
                               LAI_max_per_cohort = PAI/3)

  musica_out <- callmusica_vegetation_structure(
    crown_clumping_factor = inter_crown_clumping1,
    microclimate_height = microclimate_height,
    leaf_param = list("musica_veg_FR-Bil_pinus_pinaster" =
                        list(phenology = phenology1,
                             allometry = allometry1,
                             PAI = PAI)),
    run.name = paste0(PAI),
    musica.cmd = musica.cmd,
    save.output = TRUE,
    keep.tmp = FALSE)
  musica_out$PAI <- PAI
  musica_out
}
