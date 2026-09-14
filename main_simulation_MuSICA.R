rm(list=ls(all=TRUE)) # Clear the global environment
gc() # Trigger garbage collector

# --- 1. Setup and Libraries ---

if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
cat("Working Directory:", getwd(), "\n")

library(musica.tools)
library(rmusica)
library(ncdf4)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# --- 2. Helper Functions ---

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
    dens <- dbeta(x_seq, shape1 = 7, shape2 = 3)
  } else if (type == "bottom_heavy") {
    # Skewed towards bottom
    dens <- dbeta(x_seq, shape1 = 2, shape2 = 8)
  } else if (type == "symmetric") {
    # Gaussian-like
    dens <- dbeta(x_seq, shape1 = 5, shape2 = 5)
  } else if (type == "top_bottom_heavy") {
    # Bimodal: Sum of Top Heavy and Bottom Heavy distributions
    dens <- dbeta(x_seq, shape1 = 4, shape2 = 20) + dbeta(x_seq, shape1 = 20, shape2 = 4) + 2
  } else if (type == "top_peak_uniform") {
    # Description: Peak at the very top, roughly uniform below.
    # Logic: A strong top-skewed Beta added to a constant baseline (uniform).
    # shape1=15, shape2=2 pushes the peak to the very top.
    # + 0.5 adds the "uniform" base layer.
    dens <- dbeta(x_seq, shape1 = 8, shape2 = 2) + 1
  } else if (type == "extreme_bottom") {
    # Description: Weak all along except huge increase at the bottom.
    # Logic: A very strong bottom-skewed Beta with a tiny baseline.
    # shape1=1.2 (pushes peak near 0), shape2=15 (decays very fast).
    # + 0.05 represents the "weak" vegetation along the trunk/branches.
    dens <- dbeta(x_seq, shape1 = 1.1, shape2 = 15)
  } else if (type == "low_symmetric") {
    # "Ventru vers le bas": Gaussian-like but peaked at ~35% height instead of 50%
    # Seen in plots 41_57, 41_08
    dens <- dbeta(x_seq, shape1 = 3, shape2 = 5)
    
  } else if (type == "high_symmetric") {
    # "Ventru vers le haut": Gaussian-like but peaked at ~65% height instead of 50%
    dens <- dbeta(x_seq, shape1 = 7, shape2 = 5)
    
  } else {
    stop(paste("Unknown profile type:", type))
  }
  
  # Normalize to sum = 1 (Density shape only)
  if (sum(dens) > 0) {
    dens <- dens / sum(dens)
  }
  
  return(data.frame(height = z_levels, density = dens))
}

# --- 3. Configuration & Parameters ---

# Output directory
results_dir <- './out_files/Sensitivity_Analysis'
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# Simulation command
musica.cmd <- "bash -i -c musica > musica.log"

# Forcing file configuration
# NOTE: Ensure this path points to a valid generic forcing file (e.g. at 2m or above canopy)
forcing_types <- c("ERA5") 
forcing_file_base <- "./in_files/musica_in_Blois.nc" # Or Safran path if preferred

# Sensitivity Parameters
target_heights <- c(10, 20, 30, 40)
profile_shapes <- c(
                    "uniform",
                    # "top_heavy", 
                    # "bottom_heavy",
                    # "top_bottom_heavy", 
                    "symmetric",
                    # "top_peak_uniform",
                    # "extreme_bottom",
                    "low_symmetric"
                    # "high_symmetric"
                    )
# profile_shapes <- c("top_bottom_heavy")
lai_values      <- seq(1, 8, by = 1)
clumping_vals   <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)

# --- 4. Initialize Base Phenology ---

# We generate a generic phenology object once.
# The 'Leaf_area_1yr' (Magnitude) will be overwritten inside the loop.
phenology_base <- calc_phenology(
  list.year = c(2021, 2022),
  nleafage = 1,
  budburst_date = 115,            # DOY 115 (approx late April)
  leaf_age_max_in = 0.56,         # ~200 days lifespan
  relative_age_firstmax = 0.10,
  relative_age_lastmax = 0.75,
  LAI_max_per_cohort = 1          # Placeholder, updated in loop
)

# --- 5. Main Execution Loop ---

cat("Starting Sensitivity Analysis Loops...\n")

for (forcing_type in forcing_types) {
  
  # Select forcing file (modify logic if filenames differ by type)
  forcing_file <- if (forcing_type == "Safran") {
    "./in_files/musica_in_Safran_Blois_2021.nc"
  } else {
    "./in_files/musica_in_Blois.nc"
  }
  
  # Loop 1: Canopy Height
  for (h_canopy in target_heights) {
    cat(sprintf("> Height: %d m\n", h_canopy))
    
    # Loop 2: Profile Shape
    for (p_shape in profile_shapes) {
      cat(sprintf("  >> Shape: %s\n", p_shape))
      
      # Generate the synthetic allometry (Z and Density)
      # z = 2.5, 3.5 ... h_canopy-0.5
      allometry <- get_synthetic_profile(height_target = h_canopy, type = p_shape)
      
      for (clumping_val in clumping_vals) {
        cat(sprintf("    >>> Clumping: %.2f\n", clumping_val))
        
        # Loop 3: LAI Magnitude
        for (daily_lai in lai_values) {
          cat(sprintf("      >>>> LAI: %.2f\n", daily_lai))
          
          # --- Prepare Parameters ---
          
          # Calculate PAI (MuSICA often assumes PAI = 2 * LAI for random distribution approximation)
          plantareaindex_forced <- 2 * daily_lai
          
          # Reset phenology object to avoid pointer issues
          phenology_run <- phenology_base
          
          # Inject the specific PAI magnitude into the 2021 timeline
          phenology_run$Leaf_area_1yr[which(phenology_run$year == 2021)] <- plantareaindex_forced
          
          # Define output filename structure
          # Structure: .../height_XM/shape/lai_X/...
          history_file <- file.path(
            results_dir,
            paste0("height_", h_canopy, "m"),
            p_shape,
            paste0("clump_", clumping_val),
            paste0("lai_", daily_lai),
            paste0("musica_out_", forcing_type, "_synthetic.nc")
          )
          
          # Ensure directory exists
          if (!dir.exists(dirname(history_file))) {
            dir.create(dirname(history_file), recursive = TRUE, showWarnings = FALSE)
          }
          
          # --- Run MuSICA ---
          
          tryCatch({
            musica_out <- callmusica(
              musica.param = list(
                "setupctl" = list(
                  "clumping_factor"  = clumping_val,
                  "forcing_filename" = forcing_file,
                  "history_filename" = history_file,
                  "forcing_height"   = h_canopy + 2
                )
              ),
              leaf.param = list(
                "musica_veg1" = list(
                  "phenology" = phenology_run,
                  "allometry" = allometry,
                  "leafphenologyctl" = list(
                    "lai_max_per_cohort" = plantareaindex_forced
                  ),
                  "leafallometryctl" = list(
                    "canopy_height_top"    = h_canopy,
                    "canopy_height_bottom" = 2
                  ),
                  "leafmusicactl" = list(
                    "canopy_height_top" = h_canopy
                  )
                )
              ),
              musica.cmd = musica.cmd,
              keep.tmp = FALSE,
              out.netcdf = TRUE
            )
          }, error = function(e) {
            cat(sprintf("Error running MuSICA for H=%d, Shape=%s, LAI=%d: %s\n", 
                        h_canopy, p_shape, clumping_val, daily_lai, e$message))
          })
          
          # Clean memory after every run
          gc()
        } # End LAI loop
      } # End Clumping loop
    } # End Shape loop
  } # End Height loop
} # End Forcing loop

cat("Simulation complete.\n")

# ==============================================================================
# 6. VISUALIZATION OF GENERATED VERTICAL PROFILES
# ==============================================================================

# cat("Generating profile visualization plot...\n")
# 
# # 1. Compile data from all scenarios
# df_plot_all <- data.frame()
# 
# for (h in target_heights) {
#   for (shape in profile_shapes) {
#     
#     # Retrieve profile (normalized density)
#     tmp_profile <- get_synthetic_profile(height_target = h, type = shape)
#     
#     # Add metadata for plotting
#     tmp_profile$Target_Height <- paste0(h, " m")
#     tmp_profile$Shape <- shape
#     
#     # Example conversion for visualization (arbitrary LAI scaler)
#     tmp_profile$LAD_Example <- tmp_profile$density * 5
#     
#     df_plot_all <- rbind(df_plot_all, tmp_profile)
#   }
# }
# 
# # Order factors for logical display
# df_plot_all$Target_Height <- factor(df_plot_all$Target_Height,
#                                     levels = paste0(target_heights, " m"))
# 
# # 2. Create the plot
# p <- ggplot(df_plot_all, aes(x = density, y = height, color = Shape)) +
#   # Lines connecting height steps
#   geom_path(size = 1, alpha = 0.8) +
#   # Points to highlight discrete layers
#   geom_point(size = 1.5) +
#   # Facet by canopy height (free Y scales)
#   facet_wrap(~Target_Height, scales = "free_y", ncol = 5) +
#   # Aesthetics
#   theme_bw() +
#   labs(
#     title = "Synthetic Vertical Profiles Tested (Normalized Density)",
#     subtitle = "Comparison of shapes (Uniform, Top/Bottom Heavy, Symmetric) by canopy height",
#     x = "Probability Density (sum = 1)",
#     y = "Height (m)",
#     color = "Profile Type"
#   ) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(face = "bold", size = 14),
#     strip.text = element_text(face = "bold", size = 12)
#   )
# 
# # 3. Print and Save
# print(p)
# ggsave(filename = file.path(results_dir, "Synthetic_Profiles_Visualization.png"),
#        width = 12, height = 8, dpi = 300)
# 
# # 4. Create the plot
# p2 <- ggplot(df_plot_all, aes(x = density, y = height, color = Target_Height)) +
#   # Lines connecting height steps
#   geom_path(size = 1, alpha = 0.8) +
#   # Points to highlight discrete layers
#   geom_point(size = 1.5) +
#   # Facet by canopy height (free Y scales)
#   facet_wrap(~Shape, scales = "free_y", ncol = 5) +
#   # Aesthetics
#   theme_bw() +
#   labs(
#     title = "Synthetic Vertical Profiles Tested (Normalized Density)",
#     subtitle = "Comparison of shapes (Uniform, Top/Bottom Heavy, Symmetric) by canopy height",
#     x = "Probability Density (sum = 1)",
#     y = "Height (m)",
#     color = "Profile Type"
#   ) +
#   theme(
#     legend.position = "bottom",
#     plot.title = element_text(face = "bold", size = 14),
#     strip.text = element_text(face = "bold", size = 12)
#   )
# 
# # 5. Print and Save
# print(p2)
# ggsave(filename = file.path(results_dir, "Synthetic_Profiles_Visualization2.png"),
#        width = 12, height = 8, dpi = 300)
