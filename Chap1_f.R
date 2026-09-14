# ==============================================================================
# Forest Structure and Microclimate: A Coupled LiDAR and Physical Modeling Approach
# Author: Nathan Corroyez
# Description: Sensitivity analysis of sub-canopy temperatures (ΔTmax) to vertical
#              Leaf Area Density (LAD) profiles using MuSICA mechanistic model,
#              LiDAR structural metrics, and GAMM statistical modeling.
# ==============================================================================


# ==============================================================================
# 0. CONFIGURATION & LIBRARIES
# ==============================================================================

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(musica.tools)
library(rmusica)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(broom.mixed)
library(ggeffects)
library(patchwork)
library(viridis)
library(mgcv)
library(MASS)
library(terra)
library(fda)
library(corrplot)
library(sf)
library(clhs)
library(vip)
library(Metrics)
library(tidyterra)

# --- Graphical theme ---
theme_set(
  theme_bw() +
    theme(
      text             = element_text(size = 12),
      plot.title       = element_text(face = "bold", size = 14),
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey95"),
      strip.text       = element_text(face = "bold")
    )
)

set.seed(42)

# --- File paths (adjust to local environment) ---
in_dir  <- "in_files"
out_dir <- "out_files/Sensitivity_Analysis"

# --- Simulation parameters ---
macro_nc_file <- file.path(in_dir, "musica_in_Blois.nc")
metrics_file  <- file.path(in_dir, "metrics_results_25.csv")
date_seq      <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")
ids_to_remove <- c("41_13", "41_14", "41_20", "41_41", "41_50", "41_51", "41_53")


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Compute uniformity metrics for a sampling distribution
#'
#' @param data         A data frame containing the variable.
#' @param variable_name  Name of the column to evaluate.
#' @param n_bins       Number of bins for the histogram (default 15).
#' @return A data frame with KS distance and Shannon evenness index.
calc_uniformity <- function(data, variable_name, n_bins = 15) {
  x <- na.omit(data[[variable_name]])
  if (length(x) == 0) return(NULL)
  
  # Kolmogorov-Smirnov distance vs theoretical uniform distribution
  ks_res <- suppressWarnings(ks.test(x, "punif", min(x), max(x)))
  d_stat <- unname(ks_res$statistic)
  
  # Shannon evenness index
  breaks  <- seq(min(x), max(x), length.out = n_bins + 1)
  counts  <- hist(x, breaks = breaks, plot = FALSE)$counts
  p       <- counts[counts > 0] / sum(counts)
  H       <- -sum(p * log(p))
  evenness <- H / log(n_bins)   # J = H / H_max; 0 = skewed, 1 = perfectly uniform
  
  data.frame(
    Variable        = variable_name,
    KS_Distance     = round(d_stat, 3),   # closer to 0 is better
    Shannon_Evenness = round(evenness, 3)  # closer to 1 is better
  )
}

#' Normalize and prepare LAD profiles for representativeness comparison
#'
#' @param df_input  Data frame with Archetype, Hmax, and LAD_Layer_* columns.
#' @return Long-format data frame with relative and absolute height columns.
prepare_profiles <- function(df_input) {
  df_input %>%
    filter(!is.na(Archetype)) %>%
    mutate(
      plot_id      = row_number(),
      Main_Cluster = str_extract(Archetype, "^[^_]+")
    ) %>%
    dplyr::select(plot_id, Main_Cluster, Hmax, starts_with("LAD_Layer_")) %>%
    pivot_longer(
      cols      = starts_with("LAD_Layer_"),
      names_to  = "Layer",
      values_to = "LAD"
    ) %>%
    mutate(
      abs_height = as.numeric(str_replace(Layer, "LAD_Layer_", "")),
      LAD        = replace_na(LAD, 0)
    ) %>%
    filter(abs_height <= ceiling(Hmax) + 1) %>%
    group_by(plot_id) %>%
    mutate(
      rel_height = pmin(abs_height / Hmax, 1),
      sum_LAD    = sum(LAD, na.rm = TRUE),
      rel_LAD    = if_else(sum_LAD > 0, LAD / sum_LAD, 0)
    ) %>%
    ungroup()
}

#' Force UTC timezone on a NetCDF time variable
#'
#' @param ncfile   Path to the NetCDF file.
#' @param varname  Name of the time variable (default "time").
#' @return POSIXct vector in UTC.
force_utc_nc <- function(ncfile, varname = "time") {
  if (!file.exists(ncfile)) stop(paste("File not found:", ncfile))
  nc       <- nc_open(ncfile)
  on.exit(nc_close(nc))
  time_val  <- ncvar_get(nc, varname)
  units_str <- ncatt_get(nc, varname, "units")$value
  parts     <- strsplit(units_str, " since ")[[1]]
  if (length(parts) < 2) stop("Unrecognized NetCDF date format")
  unit_type  <- parts[1]
  origin_str <- gsub("\\s*\\(.*\\)", "", parts[2])
  origin_str <- gsub(" UTC", "", origin_str)
  origin_dt  <- as.POSIXct(origin_str, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  mult <- switch(trimws(unit_type),
                 "hours"   = 3600,
                 "minutes" = 60,
                 "seconds" = 1,
                 "days"    = 86400,
                 3600
  )
  with_tz(origin_dt + (time_val * mult), "UTC")
}


# ==============================================================================
# SECTION 1. STRUCTURAL CHARACTERIZATION (TYPOLOGY)
# ==============================================================================

# --- 1.1  Load 10-m LiDAR rasters ---
r_lai_10m    <- rast(file.path(in_dir, "lai_z1_res_10_m.tif"))
r_vci_10m    <- rast(file.path(in_dir, "vci_res_10_m.tif"))
r_hmax_10m   <- rast(file.path(in_dir, "max_res_10_m.tif"))
r_fcover_10m <- rast(file.path(in_dir, "fCover_res_10_m.tif"))
r_lad_stack_10m <- rast(file.path(in_dir, "lad_profiles_z1_res_10_m.tif"))

fact <- 2 # 1

# --- 1.2  Aggregate to 20-m resolution (matches S2 optical fetch) ---
r_metrics_20m <- terra::aggregate(
  c(r_lai_10m, r_vci_10m, r_hmax_10m, r_fcover_10m),
  fact = fact, fun = "mean", na.rm = TRUE
)
names(r_metrics_20m) <- c("LAI", "VCI", "Hmax", "fCover")

r_lad_stack_20m <- terra::aggregate(r_lad_stack_10m, 
                                    fact = fact, fun = "mean", na.rm = TRUE)
r_all_20m       <- c(r_metrics_20m, r_lad_stack_20m)

# --- 1.3  Convert to data frame, keep valid forest pixels ---
df_all_pixels <- as.data.frame(r_all_20m, xy = TRUE, na.rm = FALSE) %>%
  filter(if_any(-c(x, y), ~ !is.na(.)))

mat_lad   <- df_all_pixels %>% dplyr::select(starts_with("LAD_Layer_")) %>% as.matrix()
row_sums  <- rowSums(mat_lad, na.rm = TRUE)
valid_idx <- row_sums > 0.5

df_forest       <- df_all_pixels[valid_idx, ]
mat_lad         <- mat_lad[valid_idx, ]
mat_lad[is.na(mat_lad)] <- 0
row_sums        <- row_sums[valid_idx]


# --- 1.4  Functional PCA on normalized LAD profiles (pure shape) ---
# Toggle: TRUE = relative height (Z/Hmax), FALSE = absolute height
normalize_height <- TRUE

# Normalize LAD so area under curve = 1 (controls for LAI)
mat_lad_norm  <- sweep(mat_lad, 1, row_sums, FUN = "/")
lad_col_names <- colnames(mat_lad_norm)
z_breaks      <- as.numeric(gsub("LAD_Layer_", "", lad_col_names))

if (normalize_height) {
  cat("Normalizing height to relative scale (Z/Hmax)...\n")
  n_bins      <- length(z_breaks)
  z_rel_grid  <- seq(0, 1, length.out = n_bins)
  mat_lad_final <- matrix(0, nrow = nrow(mat_lad_norm), ncol = n_bins)
  
  for (i in seq_len(nrow(mat_lad_norm))) {
    pixel_hmax  <- max(df_forest$Hmax[i], 1)
    pixel_z_rel <- z_breaks / pixel_hmax
    interp      <- approx(x = pixel_z_rel, y = mat_lad_norm[i, ], xout = z_rel_grid, rule = 2)
    mat_lad_final[i, ] <- pmax(interp$y, 0)
  }
  
  # Re-normalize after interpolation
  rel_row_sums <- rowSums(mat_lad_final, na.rm = TRUE)
  rel_row_sums[rel_row_sums == 0] <- 1
  mat_lad_final <- sweep(mat_lad_final, 1, rel_row_sums, FUN = "/")
  arg_vals      <- z_rel_grid
  
} else {
  cat("Using absolute height...\n")
  mat_lad_final <- mat_lad_norm
  arg_vals      <- z_breaks
}

# Functional data smoothing with B-spline basis
# Define the range of basis functions to test
nbasis_range <- 4:15
gcv_scores   <- numeric(length(nbasis_range))

cat("Calculating Mean GCV for nbasis optimization...\n")
for (i in seq_along(nbasis_range)) {
  test_basis <- create.bspline.basis(rangeval = c(min(arg_vals), max(arg_vals)), nbasis = nbasis_range[i])
  smooth_obj <- smooth.basis(argvals = arg_vals, y = t(mat_lad_final), fdParobj = test_basis)
  gcv_scores[i] <- mean(smooth_obj$gcv)
}

# Mathematical elbow: maximum perpendicular distance to the chord for GCV
p1_gcv <- c(nbasis_range[1], gcv_scores[1])
n_range <- length(nbasis_range)
p2_gcv <- c(nbasis_range[n_range], gcv_scores[n_range])

distances_gcv <- sapply(1:n_range, function(i) {
  p0 <- c(nbasis_range[i], gcv_scores[i])
  abs((p2_gcv[2] - p1_gcv[2]) * p0[1] - (p2_gcv[1] - p1_gcv[1]) * p0[2] + p2_gcv[1] * p1_gcv[2] - p2_gcv[2] * p1_gcv[1]) /
    sqrt((p2_gcv[2] - p1_gcv[2])^2 + (p2_gcv[1] - p1_gcv[1])^2)
})

optimal_nbasis <- nbasis_range[which.max(distances_gcv)]
cat(sprintf("Mathematical Elbow detected at nbasis = %d\n", optimal_nbasis))

# Plot GCV Elbow
plot(nbasis_range, gcv_scores, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of basis functions (nbasis)",
     ylab = "Mean GCV Score",
     main = "Elbow Method for nbasis Selection")
abline(v = optimal_nbasis, col = "red", linetype = "dashed", lwd = 2)

# Functional data smoothing with the OPTIMAL B-spline basis
basis    <- create.bspline.basis(rangeval = c(min(arg_vals), max(arg_vals)), nbasis = optimal_nbasis)
fd_obj   <- Data2fd(argvals = arg_vals, y = t(mat_lad_final), basisobj = basis)

# Extract top 3 functional principal components
fpca_res <- pca.fd(fd_obj, nharm = 3)
cat("Variance of pure shape captured by FPC 1-3:\n")
print(round(fpca_res$varprop * 100, 1))

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot.pca.fd(fpca_res, harm = 1:3, expand = 0.5, lwd = 2)
par(mfrow = c(1, 1))

df_forest$FPC1 <- fpca_res$scores[, 1]
df_forest$FPC2 <- fpca_res$scores[, 2]
df_forest$FPC3 <- fpca_res$scores[, 3]

# ==============================================================================
# --- VISUAL CHECK: EMPIRICAL MEAN PROFILES FOR ALL FPCs (2x3 GRID) ---
# ==============================================================================

# Helper function to generate Relative and Absolute plots for a given FPC
plot_fpc_empirical <- function(fpc_scores, fpc_name, title_desc, col_low, col_high, hmax_vec) {
  
  # 1. Identify extreme pixels (lowest 5% and highest 5%)
  q_low  <- quantile(fpc_scores, 0.05, na.rm = TRUE)
  q_high <- quantile(fpc_scores, 0.95, na.rm = TRUE)
  idx_low  <- which(fpc_scores <= q_low)
  idx_high <- which(fpc_scores >= q_high)
  
  # 2. Extract Relative Profiles (Normalized shape)
  prof_rel_low  <- colMeans(mat_lad_final[idx_low, ], na.rm = TRUE)
  prof_rel_mean <- colMeans(mat_lad_final, na.rm = TRUE)
  prof_rel_high <- colMeans(mat_lad_final[idx_high, ], na.rm = TRUE)
  
  df_rel <- data.frame(
    Height = rep(arg_vals, 3), 
    LAD    = c(prof_rel_low, prof_rel_mean, prof_rel_high),
    Group  = factor(rep(c("Lowest 5%", "Average", "Highest 5%"), each = length(arg_vals)),
                    levels = c("Lowest 5%", "Average", "Highest 5%"))
  )
  
  # 3. Extract Absolute Profiles (True physical volume)
  prof_abs_low  <- colMeans(mat_lad[idx_low, ], na.rm = TRUE)
  prof_abs_mean <- colMeans(mat_lad, na.rm = TRUE)
  prof_abs_high <- colMeans(mat_lad[idx_high, ], na.rm = TRUE)
  
  # Calculate mean Hmax for each group to truncate trailing zeros
  mean_hmax_low  <- mean(hmax_vec[idx_low], na.rm = TRUE)
  mean_hmax_mean <- mean(hmax_vec, na.rm = TRUE)
  mean_hmax_high <- mean(hmax_vec[idx_high], na.rm = TRUE)
  
  # Build individual data frames and filter out the empty air above the canopy
  df_abs_low  <- data.frame(Height = z_breaks, LAD = prof_abs_low, Group = "Lowest 5%")
  df_abs_low  <- df_abs_low[df_abs_low$Height <= mean_hmax_low, ]
  
  df_abs_mean <- data.frame(Height = z_breaks, LAD = prof_abs_mean, Group = "Average")
  df_abs_mean <- df_abs_mean[df_abs_mean$Height <= mean_hmax_mean, ]
  
  df_abs_high <- data.frame(Height = z_breaks, LAD = prof_abs_high, Group = "Highest 5%")
  df_abs_high <- df_abs_high[df_abs_high$Height <= mean_hmax_high, ]
  
  # Combine truncated profiles
  df_abs <- rbind(df_abs_low, df_abs_mean, df_abs_high)
  df_abs$Group <- factor(df_abs$Group, levels = c("Lowest 5%", "Average", "Highest 5%"))
  
  # 4. Plot Relative (Top Row)
  p_rel <- ggplot(df_rel, aes(x = LAD, y = Height, color = Group, linetype = Group)) +
    geom_path(linewidth = 1.2) +
    scale_color_manual(values = c(col_low, "grey50", col_high)) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    labs(title = paste(fpc_name, "-", title_desc),
         subtitle = "Relative Space (Pure Shape)",
         x = "Relative LAD (Normalized)", y = "Relative Height (Z/Hmax)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, face = "italic"))
  
  # 5. Plot Absolute (Bottom Row)
  p_abs <- ggplot(df_abs, aes(x = LAD, y = Height, color = Group, linetype = Group)) +
    geom_path(linewidth = 1.2) +
    scale_color_manual(values = c(col_low, "grey50", col_high)) +
    scale_linetype_manual(values = c("solid", "dashed", "solid")) +
    labs(title = NULL, 
         subtitle = "Absolute Space (Physical Volume)",
         x = expression("Absolute LAD (" * m^2 / m^3 * ")"), y = "Height (m)") +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(hjust = 0.5, face = "italic"))
  
  return(list(rel = p_rel, abs = p_abs))
}

# --- Generate the lists of plots for FPC1, FPC2, FPC3 ---
# Il faut maintenant passer la variable df_forest$Hmax à la fonction
plots_fpc1 <- plot_fpc_empirical(df_forest$FPC1, "FPC1", "Top vs Bottom Heavy", "#fca50a", "#31688e", df_forest$Hmax)
plots_fpc2 <- plot_fpc_empirical(df_forest$FPC2, "FPC2", "Vertical Dispersion", "#fca50a", "#31688e", df_forest$Hmax)
plots_fpc3 <- plot_fpc_empirical(df_forest$FPC3, "FPC3", "Micro-Complexity", "#fca50a", "#31688e", df_forest$Hmax)

# --- Combine everything into a 2x3 Grid using Patchwork ---
library(patchwork)
combined_fpcs <- (plots_fpc1$rel | plots_fpc2$rel | plots_fpc3$rel) /
  (plots_fpc1$abs | plots_fpc2$abs | plots_fpc3$abs) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title = element_blank())

combined_fpcs <- combined_fpcs + plot_annotation(
  title = "Empirical Validation of Functional Principal Components (FPCs)",
  subtitle = "Averaged real LAD profiles for the extreme 5% of scores (Top: Relative / Bottom: Absolute)",
  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                plot.subtitle = element_text(size = 12, hjust = 0.5))
)

print(combined_fpcs)

# Save the plot
# ggsave("FPC_Validation_2x3.png", plot = combined_fpcs, width = 12, height = 9, dpi = 300)

# --- 1.5  K-Means clustering: elbow method ---
df_cluster_input <- df_forest %>%
  dplyr::select(LAI, Hmax, FPC1, FPC2, FPC3)

is_clean  <- complete.cases(df_cluster_input)
df_clean  <- df_cluster_input[is_clean, ]
df_scaled <- scale(df_clean)

# Correlation matrix of structural variables
cor_matrix <- cor(df_clean, method = "pearson")
corrplot(
  cor_matrix,
  method = "color", type = "upper", addCoef.col = "black",
  tl.col = "black", tl.srt = 45, diag = FALSE,
  col   = colorRampPalette(c("#4477AA", "#FFFFFF", "#BB4444"))(200),
  title = "Correlation Matrix of Structural Variables",
  mar   = c(0, 0, 2, 0)
)

# Compute total within-cluster SS for k = 1:10
cat("Calculating WSS for the Elbow Method on all valid pixels...\n")
max_k <- 10
wss   <- numeric(max_k)
for (k in 1:max_k) {
  km     <- kmeans(df_scaled, centers = k, iter.max = 100, nstart = 25)
  wss[k] <- km$tot.withinss
}

# Mathematical elbow: maximum perpendicular distance to the chord
p1 <- c(1, wss[1])
p2 <- c(max_k, wss[max_k])
distances <- sapply(1:max_k, function(i) {
  p0 <- c(i, wss[i])
  abs((p2[2] - p1[2]) * p0[1] - (p2[1] - p1[1]) * p0[2] + p2[1] * p1[2] - p2[2] * p1[1]) /
    sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
})
optimal_k <- which.max(distances)
cat(sprintf("Mathematical Elbow detected at k = %d\n", optimal_k))

plot(1:max_k, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares",
     main = "Elbow Method")
abline(v = optimal_k, col = "red", linetype = "dashed", lwd = 2)

# --- 1.6  Final K-Means & spatial mapping ---
final_km <- kmeans(df_scaled, centers = optimal_k, iter.max = 100, nstart = 25)

df_forest$Cluster            <- NA
df_forest$Cluster[is_clean]  <- as.factor(final_km$cluster)

r_landscape_clusters         <- rast(r_all_20m[[1]])
values(r_landscape_clusters) <- NA

df_mapped <- df_forest %>% filter(!is.na(Cluster))
r_landscape_clusters[cellFromXY(r_landscape_clusters, df_mapped[, c("x", "y")])] <-
  as.numeric(df_mapped$Cluster)

plot(
  r_landscape_clusters,
  col  = viridis::turbo(optimal_k),
  main = sprintf("Forest Archetypes (k = %d)", optimal_k),
  type = "classes"
)

# Pixel count per cluster
cluster_counts <- df_forest %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarise(
    Pixel_Count = n(),
    Percentage  = round((n() / sum(!is.na(df_forest$Cluster))) * 100, 1)
  ) %>%
  arrange(desc(Pixel_Count))
print(cluster_counts)


# --- 1.7  Ecological profiling: structural boxplots per cluster ---
df_plot <- df_forest %>%
  filter(!is.na(Cluster)) %>%
  mutate(Cluster = as.factor(Cluster)) %>%
  dplyr::select(Cluster, LAI, Hmax, fCover, FPC1, FPC2, FPC3) %>%
  pivot_longer(cols = -Cluster, names_to = "Metric", values_to = "Value")

p_profiles <- ggplot(df_plot, aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.alpha = 0.05) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title = "Structural Identity of the Forest Archetypes",
    x     = "K-Means Cluster",
    y     = "Metric Value"
  ) +
  theme(legend.position = "none")
print(p_profiles)

# Assign ecological archetype labels
archetype_names <- setNames(as.character(1:5), as.character(1:5))  # update with real names
df_forest <- df_forest %>%
  mutate(Archetype = recode_factor(as.character(Cluster), !!!archetype_names))

cat("\nLandscape composition by Archetype:\n")
print(table(df_forest$Archetype))


# --- 1.8  Intra-cluster correlation matrices ---
vars_to_correlate <- c("LAI", "Hmax", "fCover", "FPC1", "FPC2", "FPC3")

par(mfrow = c(2, 2))
for (i in 1:4) {
  df_sub_clean <- df_forest %>%
    filter(Cluster == i) %>%
    dplyr::select(all_of(vars_to_correlate)) %>%
    na.omit()
  
  cor_mat <- cor(df_sub_clean, method = "pearson")
  corrplot(
    cor_mat,
    method = "color", type = "upper", addCoef.col = "black",
    number.cex = 0.7, tl.col = "black", tl.srt = 45, diag = FALSE,
    col   = colorRampPalette(c("#4477AA", "#FFFFFF", "#BB4444"))(200),
    title = paste("Cluster", i, "Correlations"),
    mar   = c(0, 0, 2, 0)
  )
}
par(mfrow = c(1, 1))


# --- 1.9  Spatial mapping of Functional Principal Components ---
r_template         <- rast(r_metrics_20m[[1]]); values(r_template) <- NA
r_fpc1 <- r_fpc2 <- r_fpc3 <- r_template

cell_idx        <- cellFromXY(r_template, df_forest[, c("x", "y")])
r_fpc1[cell_idx] <- df_forest$FPC1
r_fpc2[cell_idx] <- df_forest$FPC2
r_fpc3[cell_idx] <- df_forest$FPC3

r_fpcs       <- c(r_fpc1, r_fpc2, r_fpc3)
names(r_fpcs) <- c(
  "FPC1",
  "FPC2",
  "FPC3"
)
fpc_names <- c("FPC1", "FPC2", "FPC3")
# Divergent palette: blue = negative, white = 0, red = positive
div_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
par(mfrow = c(1, 3), mar = c(2, 2, 4, 4))
for(i in 1:3) {
  r_current <- r_fpcs[[i]]
  val_max <- max(abs(minmax(r_current)), na.rm = TRUE)
  plot(
    r_current,
    col   = div_pal,
    range = c(-val_max, val_max),
    axes  = FALSE, 
    box   = FALSE,
    main  = fpc_names[i],
    plg   = list(title = "Score", cex = 0.8)
  )
}
par(mfrow = c(1, 1))

# ==============================================================================
# SECTION 2. SAMPLING DESIGN FOR MuSICA
# ==============================================================================

# ------------------------------------------------------------------------------
# Option A — Ecological stratification (coarse grid)
# ------------------------------------------------------------------------------

df_stratified_A <- df_forest %>%
  filter(!is.na(Archetype)) %>%
  mutate(
    LAI_class  = cut(LAI,  breaks = c(0, 2, 4, 6, Inf),
                     labels = c("Low", "Med", "High", "VeryHigh"), include.lowest = TRUE),
    Hmax_class = cut(Hmax, breaks = c(0, 10, 20, 30, Inf),
                     labels = c("Regen", "Pole", "Mature", "Giant"), include.lowest = TRUE),
    Stratum_ID = paste(Archetype, LAI_class, Hmax_class, sep = "_")
  )

# Keep strata with >= 20 pixels
valid_strata_A <- df_stratified_A %>%
  count(Stratum_ID) %>%
  filter(n >= 20) %>%
  pull(Stratum_ID)

df_valid_strata_A <- df_stratified_A %>% filter(Stratum_ID %in% valid_strata_A)
cat(sprintf(
  "Number of valid ecological strata identified: %d out of %d possible combinations\n",
  length(valid_strata_A), 5 * 4 * 4
))

# Balanced random sampling (N = 10 per stratum)
df_musica_sample_A <- df_valid_strata_A %>%
  group_by(Stratum_ID) %>%
  sample_n(size = min(n(), 10)) %>%
  ungroup()

cat(sprintf("Total LAD profiles extracted (Option A): %d\n", nrow(df_musica_sample_A)))

# Export data frame
df_export_A <- df_musica_sample_A %>%
  dplyr::select(x, y, Archetype, LAI, Hmax, fCover, VCI,
                FPC1, FPC2, FPC3, LAI_class, Hmax_class, Stratum_ID,
                starts_with("LAD_Layer_"))

cat("\nSampling Matrix (Archetype × LAI):\n");  print(table(df_export_A$Archetype, df_export_A$LAI_class))
cat("\nSampling Matrix (Archetype × Height):\n"); print(table(df_export_A$Archetype, df_export_A$Hmax_class))

# Histograms of sampled structural space
df_hist_A <- df_musica_sample_A %>%
  dplyr::select(LAI, Hmax, fCover, FPC1, FPC2, FPC3) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

p_hist_A <- ggplot(df_hist_A, aes(x = Value, fill = Metric)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free", ncol = 3) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
  labs(
    title    = "Structural Distribution — Option A",
    # subtitle = "Validating balanced stratification of the physical space",
    x = "Metric Value", y = "Frequency (Number of Pixels)"
  ) +
  theme(legend.position = "none")
print(p_hist_A)

# Uniformity metrics for Option A
metrics_A <- bind_rows(lapply(
  c("Hmax", "LAI", "FPC1", "FPC2", "FPC3"),
  calc_uniformity, data = df_musica_sample_A
))
print(metrics_A)


# ------------------------------------------------------------------------------
# Option B — Fine-grained uniform stratification
# ------------------------------------------------------------------------------

df_stratified_B <- df_forest %>%
  filter(!is.na(Archetype)) %>%
  mutate(
    LAI_bin  = cut(LAI,  breaks = seq(0, ceiling(max(LAI,  na.rm = TRUE)), by = 1), include.lowest = TRUE),
    Hmax_bin = cut(Hmax, breaks = seq(0, ceiling(max(Hmax, na.rm = TRUE)), by = 5), include.lowest = TRUE),
    Stratum_ID = paste(Archetype, LAI_bin, Hmax_bin, sep = "_")
  ) %>%
  filter(!is.na(Stratum_ID))

valid_strata_B <- df_stratified_B %>%
  count(Stratum_ID) %>%
  filter(n >= 3) %>%
  pull(Stratum_ID)

cat(sprintf("Number of valid fine ecological strata: %d\n", length(valid_strata_B)))

df_musica_sample_B <- df_stratified_B %>%
  filter(Stratum_ID %in% valid_strata_B) %>%
  group_by(Stratum_ID) %>%
  sample_n(size = min(n(), 2)) %>%
  ungroup()

cat(sprintf("Total LAD profiles extracted (Option B): %d\n", nrow(df_musica_sample_B)))

df_export_B <- df_musica_sample_B %>%
  dplyr::select(x, y, Archetype, LAI, Hmax, fCover, VCI,
                FPC1, FPC2, FPC3, Stratum_ID, starts_with("LAD_Layer_"))

# Spatial map
plot(r_landscape_clusters, col = viridis::turbo(5),
     main = "Spatial Distribution — Option B (Uniform Grid)",
     type = "classes", legend = FALSE, axes = FALSE, box = FALSE)
points(df_export_B$x, df_export_B$y, pch = 21, bg = "white", col = "black", cex = 0.8, lwd = 1.2)
legend("topleft", legend = as.character(1:5), fill = viridis::turbo(5), bty = "n", cex = 0.8, inset = 0.02)

# Histograms
df_hist_B <- df_musica_sample_B %>%
  dplyr::select(LAI, Hmax, fCover, FPC1, FPC2, FPC3) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

p_hist_B <- ggplot(df_hist_B, aes(x = Value, fill = Metric)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free", ncol = 3) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
  labs(
    title    = "Structural Distribution — Option B (Uniform Grid)",
    # subtitle = "Flattening distributions via fine-grained multivariate stratification",
    x = "Metric Value", y = "Frequency (Number of Pixels)"
  ) +
  theme(legend.position = "none")
print(p_hist_B)

# Absolute LAD profiles per archetype
df_binned_B <- df_export_B %>%
  mutate(plot_id = row_number(), Main_Cluster = str_extract(Stratum_ID, "^[^_]+")) %>%
  dplyr::select(plot_id, Main_Cluster, Hmax, starts_with("LAD_Layer_")) %>%
  pivot_longer(cols = starts_with("LAD_Layer_"), names_to = "Layer", values_to = "LAD") %>%
  mutate(
    abs_height = as.numeric(str_replace(Layer, "LAD_Layer_", "")),
    LAD        = replace_na(LAD, 0)
  ) %>%
  filter(abs_height <= ceiling(Hmax) + 1) %>%
  mutate(height_bin = round(abs_height)) %>%
  group_by(Main_Cluster, height_bin) %>%
  summarise(
    mean_LAD_bin = mean(LAD,                    na.rm = TRUE),
    q25          = quantile(LAD, 0.25, na.rm = TRUE),
    q75          = quantile(LAD, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

p_profiles_B <- ggplot(df_binned_B,
                       aes(y = height_bin, color = Main_Cluster, fill = Main_Cluster)) +
  geom_ribbon(aes(xmin = q25, xmax = q75), alpha = 0.3, color = NA) +
  geom_path(aes(x = mean_LAD_bin), linewidth = 1.2) +
  facet_wrap(~ Main_Cluster, scales = "free_x") +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  labs(
    title    = "Absolute Vertical Foliage Profiles (Summarized)",
    subtitle = "Solid line: Mean LAD | Shaded area: IQR (25th–75th percentile)",
    x = expression("Absolute LAD (m"^2 * "/m"^3 * ")"),
    y = "Height Above Ground (m)"
  ) +
  theme(legend.position = "none")
# print(p_profiles_B)

# Uniformity metrics for Option B
metrics_B <- bind_rows(lapply(
  c("Hmax", "LAI", "FPC1", "FPC2", "FPC3"),
  calc_uniformity, data = df_musica_sample_B
))
print(metrics_B)


# ------------------------------------------------------------------------------
# Option C — Stratified conditioned Latin Hypercube Sampling (cLHS)
# ------------------------------------------------------------------------------

n_per_archetype <- 100

df_musica_sample_C <- df_forest %>%
  filter(!is.na(Archetype)) %>%
  split(.$Archetype) %>%
  map_dfr(function(df_subset) {
    df_lhs_vars <- df_subset %>% dplyr::select(
                                               x, y, 
                                               LAI, Hmax,
                                               FPC1, FPC2, FPC3
                                               )
    n_target    <- min(n_per_archetype, nrow(df_subset))
    lhs_res     <- clhs(df_lhs_vars, size = n_target, iter = 10000,
                        simple = FALSE, progress = FALSE)
    df_subset[lhs_res$index_samples, ]
  })

cat("\nSampling Matrix — cLHS (plots per Archetype):\n")
print(table(df_musica_sample_C$Archetype))

df_export_C <- df_musica_sample_C %>%
  dplyr::select(x, y, Archetype, LAI, Hmax, fCover, VCI,
                FPC1, FPC2, FPC3, starts_with("LAD_Layer_"))

# Histograms
df_hist_C <- df_musica_sample_C %>%
  dplyr::select(LAI, Hmax, fCover, FPC1, FPC2, FPC3) %>%
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

p_hist_C <- ggplot(df_hist_C, aes(x = Value, fill = Metric)) +
  geom_histogram(bins = 30, color = "black", alpha = 0.8) +
  facet_wrap(~ Metric, scales = "free", ncol = 3) +
  scale_fill_viridis_d(option = "mako", begin = 0.3, end = 0.8) +
  labs(
    title    = "Structural Distribution — Option C (cLHS)",
    # subtitle = "Optimized uniform sampling across multidimensional space",
    x = "Metric Value", y = "Frequency (Number of Pixels)"
  ) +
  theme(legend.position = "none")
print(p_hist_C)

# Uniformity metrics for Option C
metrics_C <- bind_rows(lapply(
  c("Hmax", "LAI", "FPC1", "FPC2", "FPC3"),
  calc_uniformity, data = df_musica_sample_C
))
print(metrics_C)

# ── Bivariate feature space: cLHS samples vs HOBO sensors ──
hobo_pts          <- st_read(file.path(in_dir, "data_Blois_utm31n.geojson"))
extracted_values  <- terra::extract(r_metrics_20m, vect(hobo_pts))
hobo_clusters    <- terra::extract(r_landscape_clusters, vect(hobo_pts))

df_hobo <- hobo_pts %>%
  bind_cols(extracted_values) %>%
  mutate(
    Archetype = hobo_clusters[, 2],
    LAI       = as.numeric(LAI), 
    Hmax      = as.numeric(Hmax)
  ) %>%
  filter(!(id_plot %in% ids_to_remove)) %>%
  drop_na()

p_scatter_C <- ggplot() +
  # Background points (cLHS samples) colored by Archetype
  geom_point(data = df_musica_sample_C,
             aes(x = LAI, y = Hmax, color = as.factor(Archetype)),
             alpha = 0.5, size = 2) +
  
  # Foreground points (HOBOs) filled by Archetype
  geom_point(data = df_hobo,
             aes(x = LAI, y = Hmax, fill = as.factor(Archetype)),
             color = "black", shape = 24, size = 3.5, stroke = 1) +
  
  # Synchronize color and fill scales to merge the legends
  scale_color_viridis_d(option = "turbo", name = "Archetype") +
  scale_fill_viridis_d(option = "turbo", name = "Archetype") +
  
  labs(
    title = "Bivariate Feature Space: cLHS Samples vs In-Situ HOBOs",
    x     = "Leaf Area Index (LAI)", 
    y     = "Maximum Height (Hmax, m)"
  ) +
  theme(legend.position = "right")

print(p_scatter_C)

# ── cLHS representativeness: sample vs full population profiles ──
# (Using Option C as the active sample for downstream analyses)
df_export <- df_export_C
df_musica_sample <- df_musica_sample_C

hobo_cells <- terra::extract(r_metrics_20m[["LAI"]], 
                             vect(hobo_pts), cells = TRUE)$cell
df_hobo_full <- df_forest[hobo_cells, ]

df_prof_sample <- prepare_profiles(df_musica_sample)
df_prof_pop    <- prepare_profiles(df_forest)
df_prof_hobo <- prepare_profiles(df_hobo_full)

shared_theme <- theme_bw(base_size = 11) + theme(legend.position = "none")
color_scale  <- scale_color_viridis_d(option = "turbo")
fill_scale   <- scale_fill_viridis_d(option = "turbo")

# ------------------------------------------------------------------------------
# LIGNE 1 : PROFILS RELATIFS (Forme pure)
# ------------------------------------------------------------------------------
p_norm_sample <- ggplot(df_prof_sample, aes(x = rel_height, y = rel_LAD, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, 1), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "A) Relative (cLHS Sample)", x = "Relative Height", y = "Relative LAD")

p_norm_pop <- ggplot(df_prof_pop, aes(x = rel_height, y = rel_LAD, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, 1), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "B) Relative (Full Population)", x = "Relative Height", y = "Relative LAD")

p_norm_hobo <- ggplot(df_prof_hobo, aes(x = rel_height, y = rel_LAD, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, 1), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "C) Relative (In-Situ HOBOs)", x = "Relative Height", y = "Relative LAD")


# ------------------------------------------------------------------------------
# LIGNE 2 : PROFILS ABSOLUS (Volume physique)
# ------------------------------------------------------------------------------
max_h <- max(df_prof_pop$abs_height, na.rm = TRUE)

p_abs_sample <- ggplot(df_prof_sample, aes(x = LAD, y = abs_height, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, max_h), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "D) Absolute (cLHS Sample)", x = "Height (m)", y = expression("LAD (" * m^2 / m^3 * ")"))

p_abs_pop <- ggplot(df_prof_pop, aes(x = LAD, y = abs_height, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, max_h), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "E) Absolute (Full Population)", x = "Height (m)", y = expression("LAD (" * m^2 / m^3 * ")"))

p_abs_hobo <- ggplot(df_prof_hobo, aes(x = LAD, y = abs_height, color = Main_Cluster, fill = Main_Cluster)) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 7), alpha = 0.2, linewidth = 1) +
  coord_flip(xlim = c(0, max_h), ylim = c(0, NA)) +
  color_scale + fill_scale + shared_theme +
  labs(title = "F) Absolute (In-Situ HOBOs)", x = "Height (m)", y = expression("LAD (" * m^2 / m^3 * ")"))


# ------------------------------------------------------------------------------
# ASSEMBLAGE FINAL (2x3) AVEC PATCHWORK
# ------------------------------------------------------------------------------
combined_validation <- (p_norm_sample | p_norm_pop | p_norm_hobo) / 
  (p_abs_sample  | p_abs_pop  | p_abs_hobo ) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "Validation of Spatial Sampling vs Entire Forest vs HOBO Network",
    subtitle = "Comparing canopy architectural shapes (Relative) and physical volumes (Absolute)",
    theme    = theme(legend.position = "right", plot.title = element_text(face = "bold"))
  )

print(combined_validation)

# ==============================================================================
# --- SPATIAL MAPPING OF SAMPLES (Option C vs HOBOs) ---
# ==============================================================================

# 1. Préparation de la carte spatiale
# On utilise r_landscape_clusters comme fond de carte
p_map_samples <- ggplot() +
  # 1. Fond de carte des archétypes (Raster converti en facteur)
  geom_spatraster(data = as.factor(r_landscape_clusters)) +
  
  # 2. Points cLHS (Option C) - Modifiés pour être BLANCS et PLUS GROS
  geom_point(data = df_musica_sample_C, 
             aes(x = x, y = y), 
             color = "black", fill = "white", shape = 21, size = 2, stroke = 0.5, alpha = 0.9) +
  
  # 3. Capteurs HOBO (Validation terrain) - Triangles
  geom_sf(data = df_hobo, 
          aes(fill = as.factor(Archetype)), 
          color = "black", shape = 24, size = 3.5, stroke = 1) +
  
  # 4. Esthétique et échelles
  scale_fill_viridis_d(option = "turbo", name = "Archetype", 
                       na.value = "transparent", na.translate = FALSE) +
  
  labs(
    title = "Spatial Distribution of Sampling Design",
    subtitle = "White dots: cLHS training samples | Triangles: HOBO sensors",
    x = "Easting (UTM 31N)", y = "Northing (UTM 31N)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold"),
    panel.grid      = element_blank()
  )

print(p_map_samples)

# ==============================================================================
# --- COMPARISON OF ALL OPTIONS (A, B, C) ---
# ==============================================================================

# Si tu veux comparer visuellement la "couverture" spatiale des 3 options :
# df_all_samples <- bind_rows(
#   df_musica_sample_A %>% mutate(Method = "A. Ecological Stratification"),
#   df_musica_sample_B %>% mutate(Method = "B. Uniform Grid"),
#   df_musica_sample_C %>% mutate(Method = "C. cLHS (Optimized)")
# )
# 
# p_map_comparison <- ggplot() +
#   geom_spatraster(data = r_landscape_clusters, alpha = 0.4) +
#   geom_point(data = df_all_samples, aes(x = x, y = y), size = 0.5) +
#   facet_wrap(~ Method) +
#   scale_fill_viridis_d(option = "turbo", guide = "none") +
#   labs(title = "Spatial Coverage Comparison between Sampling Methods") +
#   theme_void() +
#   theme(strip.text = element_text(face = "bold"))
# 
# print(p_map_comparison)

# ==============================================================================
# SECTION 3. MuSICA BATCH SIMULATIONS
# ==============================================================================

forcing_file <- file.path(in_dir, "musica_in_Blois.nc")
musica.cmd   <- "bash -i -c musica"
out_nc_dir   <- "out_files/musica_option_b_ERA5_results_0104_20mf"
if (!dir.exists(out_nc_dir)) dir.create(out_nc_dir, recursive = TRUE)

# Get the number of simulations to run
total_plots <- nrow(df_export)
# total_plots <- 1

# Optional: Global skip if the directory already contains all expected files
existing_files <- list.files(out_nc_dir, pattern = "\\.nc$")
if (length(existing_files) >= total_plots) {
  
  cat("All simulations appear to be completed already. Skipping to end.\n")
  
} else {
  
  cat(sprintf("Starting MuSICA batch processing for %d plots...\n", total_plots))
  
  for (i in seq_len(nrow(df_export))) {
    current_plot <- df_export[i, ]
    
    # Unique simulation identifier
    sim_id <- sprintf("Sim_%03d_X%d_Y%d", i,
                      round(current_plot$x), round(current_plot$y))
    
    history_file <- file.path(out_nc_dir, paste0("musica_out_", sim_id, ".nc"))
    
    # Skip if file exists (resume capability)
    if (file.exists(history_file)) {
      cat(sprintf("Skipping %s (already exists).\n", sim_id))
      next 
    }
    
    cat(sprintf("\nProcessing %d/%d: %s\n", i, nrow(df_export), sim_id))
    
    # --- Rest of your MuSICA logic ---
    canopy_height_top <- current_plot$Hmax
    plantareaindex    <- current_plot$LAI
    clumping          <- current_plot$fCover
    
    lad_cols   <- grep("LAD_Layer_", names(current_plot), value = TRUE)
    lad_values <- as.numeric(current_plot[lad_cols])
    lad_values[is.na(lad_values)] <- 0
    z_heights  <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
    
    allometry <- data.frame(height = z_heights, density = lad_values) %>%
      filter(height <= ceiling(canopy_height_top) + 1)
    
    phenology <- calc_phenology(
      list.year             = c(2021, 2022),
      nleafage              = 1,
      budburst_date         = 115,
      leaf_age_max_in       = 0.56,
      relative_age_firstmax = 0.10,
      relative_age_lastmax  = 0.75,
      LAI_max_per_cohort    = plantareaindex
    )
    
    tryCatch({
      callmusica(
        musica.param = list(
          "setupctl" = list(
            "clumping_factor"   = clumping,
            "forcing_filename"  = forcing_file,
            "history_filename"  = history_file,
            "forcing_height"    = canopy_height_top + 2
          )
        ),
        leaf.param = list(
          "musica_veg1" = list(
            "phenology"        = phenology,
            "allometry"        = allometry,
            "leafphenologyctl" = list("lai_max_per_cohort" = plantareaindex),
            "leafallometryctl" = list(
              "canopy_height_top"    = canopy_height_top,
              "canopy_height_bottom" = 1
            ),
            "leafmusicactl" = list("canopy_height_top" = canopy_height_top)
          )
        ),
        musica.cmd = musica.cmd,
        keep.tmp   = FALSE,
        out.netcdf = TRUE
      )
    }, error = function(e) {
      cat(sprintf("ERROR on %s: %s\n", sim_id, e$message))
    })
  }
}

cat("\n=== All MuSICA simulations completed. Outputs in:", out_nc_dir, "===\n")

# ==============================================================================
# SECTION 3.B. MuSICA BATCH SIMULATIONS FOR IN-SITU HOBO PLOTS
# ==============================================================================
# Objective: Simulate Tair specifically at the HOBO locations to build the
# "Triangle of Validation" (GAMM vs MuSICA vs Empirical Data).

cat("\nPreparing structural data for HOBO MuSICA simulations...\n")

# 1. Extract all structural metrics (LAI, fCover, Hmax, LAD profiles) at HOBO locations
hobo_pts <- st_read(file.path(in_dir, "data_Blois_utm31n.geojson"), 
                    quiet = TRUE)
hobo_all_extract <- terra::extract(r_all_20m, vect(hobo_pts), xy = TRUE)

# 2. Attach the HOBO plot IDs and clean the dataset
df_hobo_musica <- hobo_pts %>%
  st_drop_geometry() %>%
  dplyr::select(id_plot) %>%
  bind_cols(hobo_all_extract) %>%
  filter(!(id_plot %in% ids_to_remove))

out_nc_dir_hobo <- "out_files/musica_hobo_ERA5_results20m"
if (!dir.exists(out_nc_dir_hobo)) dir.create(out_nc_dir_hobo, recursive = TRUE)

total_hobo_plots <- nrow(df_hobo_musica)

# Optional: Global skip if the directory already contains all expected files
existing_hobo_files <- list.files(out_nc_dir_hobo, pattern = "\\.nc$")
if (length(existing_hobo_files) >= total_hobo_plots) {
  cat("All HOBO simulations appear to be completed already. Skipping to end.\n")
} else {
  cat(sprintf("Starting MuSICA batch processing for %d HOBO plots...\n", 
              total_hobo_plots))
  
  for (i in seq_len(nrow(df_hobo_musica))) {
    current_plot <- df_hobo_musica[i, ]
    plot_id_str  <- current_plot$id_plot
    
    # Unique simulation identifier using the actual HOBO ID
    sim_id <- sprintf("HOBO_%s", plot_id_str)
    history_file <- file.path(out_nc_dir_hobo, 
                              paste0("musica_out_", sim_id, ".nc"))
    
    # Skip if file exists (resume capability)
    if (file.exists(history_file)) {
      cat(sprintf("Skipping %s (already exists).\n", sim_id))
      next 
    }
    
    cat(sprintf("\nProcessing HOBO %d/%d: %s\n", i, total_hobo_plots, sim_id))
    
    # --- Extract MuSICA input parameters ---
    canopy_height_top <- current_plot$Hmax
    plantareaindex    <- current_plot$LAI
    clumping          <- current_plot$fCover
    
    lad_cols   <- grep("LAD_Layer_", names(current_plot), value = TRUE)
    lad_values <- as.numeric(current_plot[lad_cols])
    lad_values[is.na(lad_values)] <- 0
    z_heights  <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
    
    allometry <- data.frame(height = z_heights, density = lad_values) %>%
      filter(height <= ceiling(canopy_height_top) + 1)
    
    phenology <- calc_phenology(
      list.year             = c(2021, 2022),
      nleafage              = 1,
      budburst_date         = 115,
      leaf_age_max_in       = 0.56,
      relative_age_firstmax = 0.10,
      relative_age_lastmax  = 0.75,
      LAI_max_per_cohort    = plantareaindex
    )
    
    # --- Run MuSICA ---
    tryCatch({
      callmusica(
        musica.param = list(
          "setupctl" = list(
            "clumping_factor"  = clumping,
            "forcing_filename" = forcing_file,
            "history_filename" = history_file,
            "forcing_height"   = canopy_height_top + 2
          )
        ),
        leaf.param = list(
          "musica_veg1" = list(
            "phenology"        = phenology,
            "allometry"        = allometry,
            "leafphenologyctl" = list("lai_max_per_cohort" = plantareaindex),
            "leafallometryctl" = list(
              "canopy_height_top"    = canopy_height_top,
              "canopy_height_bottom" = 1
            ),
            "leafmusicactl" = list("canopy_height_top" = canopy_height_top)
          )
        ),
        musica.cmd = musica.cmd,
        keep.tmp   = FALSE,
        out.netcdf = TRUE
      )
    }, error = function(e) {
      cat(sprintf("ERROR on %s: %s\n", sim_id, e$message))
    })
  }
}

cat("\n=== All HOBO MuSICA simulations completed. Outputs in:", 
    out_nc_dir_hobo, "===\n")

# ==============================================================================
# SECTION 4. POST-PROCESSING: MACROCLIMATE REFERENCE
# ==============================================================================

cat("Extracting ERA5 macroclimate reference...\n")

nc_force    <- nc_open(forcing_file)
macro_tair  <- ncvar_get(nc_force, "Tair")
nc_close(nc_force)

macro_time <- force_utc_nc(forcing_file, varname = "time")

df_macro <- data.frame(time = macro_time, Tair_macro = macro_tair - 273.15) %>%
  filter(as.Date(time) %in% date_seq) %>%
  mutate(date = as.Date(time)) %>%
  group_by(date) %>%
  summarise(Tmax_macro = max(Tair_macro, na.rm = TRUE), .groups = "drop")

cat(sprintf("df_macro: %d summer days extracted.\n", nrow(df_macro)))


# ==============================================================================
# SECTION 5. POST-PROCESSING: DAILY MICROCLIMATE EXTRACTION
# ==============================================================================

cat("Extracting daily Tmax from all MuSICA simulations...\n")

nc_files <- list.files(out_nc_dir, pattern = "\\.nc$", full.names = TRUE)

df_daily_raw <- map_df(nc_files, function(f) {
  filename <- basename(f)
  x_val    <- as.numeric(str_extract(filename, "(?<=_X)\\d+"))
  y_val    <- as.numeric(str_extract(filename, "(?<=_Y)\\d+"))
  
  nc       <- try(nc_open(f), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  raw_data <- try(get_variable(nc, "Tair_z"), silent = TRUE)
  nc_close(nc)
  if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
  
  raw_data %>%
    filter(nair == 1) %>%
    mutate(Tair_sim = Tair_z - 273.15) %>%
    filter(as.Date(time) %in% date_seq) %>%
    mutate(time = time - hours(2), date = as.Date(time)) %>%
    group_by(date) %>%
    summarise(Tmax_micro = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
    inner_join(df_macro, by = "date") %>%
    mutate(Delta_Tmax = Tmax_micro - Tmax_macro, x = x_val, y = y_val) %>%
    dplyr::select(x, y, date, Delta_Tmax)
})

# Join with structural LiDAR descriptors
df_daily_metrics <- df_daily_raw %>%
  inner_join(
    df_export %>% mutate(x = x, y = y), 
    by = c("x", "y")
  )

cat(sprintf("Extraction complete: %d daily observations.\n", 
            nrow(df_daily_metrics)))

# ==============================================================================
# SECTION 6. STATISTICAL MODELING: VARIANCE PARTITIONING WITH GAMMs
# ==============================================================================

# --- 6.1  Create plot_id to correct for pseudo-replication ---
df_daily_metrics <- df_daily_metrics %>%
  mutate(plot_id = as.factor(paste0("X", x, "_Y", y)))

cat(sprintf("Unique plots: %d | Days: %d | Total obs: %d\n",
            length(unique(df_daily_metrics$plot_id)),
            length(unique(df_daily_metrics$date)),
            nrow(df_daily_metrics)))

# --- 6.2  Scale predictors ---
df_stats_hp <- df_daily_metrics %>%
  mutate(
    LAI_sc      = as.numeric(scale(LAI)),
    fCover_sc   = as.numeric(scale(fCover)),
    Hmax_sc     = as.numeric(scale(Hmax)),
    FPC1_sc     = as.numeric(scale(FPC1)),
    FPC2_sc     = as.numeric(scale(FPC2)),
    FPC3_sc     = as.numeric(scale(FPC3)),
    date_factor = as.factor(date),
    plot_id     = plot_id
  )

resp_var <- "Delta_Tmax"


# --- 6.3  Full corrected GAMM (with spatial random effect) ---
# s(date_factor) absorbs daily weather variability
# s(plot_id)     absorbs spatial clustering (corrects pseudo-replication)
form_full_corrected <- Delta_Tmax ~
  s(LAI_sc) +
  s(fCover_sc) +
  s(Hmax_sc) +
  s(FPC1_sc) +
  s(FPC2_sc) +
  s(FPC3_sc) +
  s(date_factor, bs = "re") +  # temporal random effect
  s(plot_id,     bs = "re")    # spatial random effect

gam_full_corrected <- bam(
  form_full_corrected, data = df_stats_hp,
  family = scat(), method = "fREML", discrete = TRUE
)

cat("\n--- Corrected full model summary ---\n")
summary(gam_full_corrected)
r2_full_corrected <- summary(gam_full_corrected)$r.sq


# --- 6.4  Original model without spatial RE (for comparison) ---
preds <- c("s(LAI_sc)", "s(fCover_sc)", "s(Hmax_sc)",
           # "ti(FPC1_sc, Hmax_sc)", 
           "s(FPC1_sc)", "s(FPC2_sc)", "s(FPC3_sc)")

form_full <- as.formula(
  paste(resp_var, "~", paste(preds, collapse = " + "), "+ s(date_factor, bs='re')")
)
gam_full <- bam(form_full, data = df_stats_hp, family = scat(), method = "fREML", discrete = TRUE)
r2_full  <- summary(gam_full)$r.sq

cat(sprintf("\nR² original (no plot_id): %.3f\n", r2_full))
cat(sprintf("R² corrected (+ plot_id): %.3f\n",  r2_full_corrected))


# --- 6.5  Concurvity diagnostics ---
cat("\n--- Global concurvity (worst-case) ---\n")
conc_full <- concurvity(gam_full_corrected, full = TRUE)
print(round(conc_full["worst", ], 3))

# Bar plot of global concurvity
conc_df <- data.frame(
  Smooth      = names(conc_full["worst", ]),
  Concurvity  = as.numeric(conc_full["worst", ])
) %>%
  filter(!grepl("re\\)", Smooth)) %>%
  arrange(desc(Concurvity))

p_conc_global <- ggplot(conc_df, aes(x = reorder(Smooth, Concurvity), y = Concurvity,
                                     fill = Concurvity > 0.8)) +
  geom_col() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red",    linewidth = 1) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "orange", linewidth = 0.8) +
  annotate("text", x = 1, y = 0.82, label = "Critical threshold (0.8)", color = "red",    hjust = 0, size = 3.5) +
  annotate("text", x = 1, y = 0.62, label = "Moderate threshold (0.6)", color = "orange", hjust = 0, size = 3.5) +
  scale_fill_manual(values = c("FALSE" = "#31688e", "TRUE" = "#d73027"), guide = "none") +
  coord_flip() +
  labs(
    title    = "Global Concurvity of GAMM Terms",
    subtitle = "Each bar = proportion of a smooth explained by ALL others combined",
    x = NULL, y = "Worst-case Concurvity"
  )
print(p_conc_global)

# Pairwise concurvity matrix
conc_pair       <- concurvity(gam_full_corrected, full = FALSE)
structural_terms <- grep("^s\\((LAI|fCover|Hmax|FPC)",
                         rownames(conc_pair$worst), value = TRUE)
mat_pair <- conc_pair$worst[structural_terms, structural_terms]

cat("\n--- Pairwise worst-case concurvity (structural terms only) ---\n")
print(round(mat_pair, 3))

corrplot(
  mat_pair,
  method = "color", type = "upper", addCoef.col = "black",
  number.cex = 0.85, tl.col = "black", tl.srt = 45, diag = FALSE,
  col   = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
  title = "Pairwise Worst-case Concurvity — Structural Terms",
  mar   = c(0, 0, 2, 0)
)

cat("\n--- Pairs with concurvity > 0.6 ---\n")
pairs_problem <- which(mat_pair > 0.6 & upper.tri(mat_pair), arr.ind = TRUE)
if (nrow(pairs_problem) == 0) {
  cat("No problematic pairs detected (all < 0.6). ✓\n")
} else {
  for (i in seq_len(nrow(pairs_problem))) {
    r <- pairs_problem[i, 1]; c <- pairs_problem[i, 2]
    cat(sprintf("  %-25s vs %-25s : %.3f %s\n",
                rownames(mat_pair)[r], colnames(mat_pair)[c], mat_pair[r, c],
                ifelse(mat_pair[r, c] > 0.8, "⚠ SEVERE", "→ moderate")))
  }
}

# Optional sensitivity test: model without Hmax (if concurvity with FPCs > 0.6)
form_no_hmax <- Delta_Tmax ~
  s(LAI_sc) +
  s(fCover_sc) +
  # s(Hmax_sc)  <-- removed: collinear with FPC2 which encodes canopy height
  s(FPC1_sc) +
  s(FPC2_sc) +
  s(FPC3_sc) +
  s(date_factor, bs = "re") +
  s(plot_id,     bs = "re")

gam_no_hmax <- bam(form_no_hmax, data = df_stats_hp,
                   family = scat(), method = "fREML", discrete = TRUE)
cat(sprintf("R² without Hmax: %.3f (vs %.3f with Hmax)\n",
            summary(gam_no_hmax)$r.sq, r2_full_corrected))
cat("If R² is near-identical → Hmax is redundant with FPC2; removing it is justified.\n")


# --- 6.6  Justification of scat() vs Gaussian ---
cat("\n--- Residual normality diagnostics (Gaussian baseline) ---\n")

gam_gaussian <- bam(form_full_corrected, data = df_stats_hp,
                    family = gaussian(), method = "fREML", discrete = TRUE)
resid_gauss  <- residuals(gam_gaussian, type = "deviance")

# Shapiro-Wilk test (limited to 5000 obs)
idx_sw  <- sample(length(resid_gauss), min(5000, length(resid_gauss)))
sw_test <- shapiro.test(resid_gauss[idx_sw])
cat(sprintf("Shapiro-Wilk (n=%d): W = %.4f, p = %.2e\n",
            length(idx_sw), sw_test$statistic, sw_test$p.value))

# Excess kurtosis
kurtosis_val <- mean((resid_gauss - mean(resid_gauss))^4) / sd(resid_gauss)^4 - 3
cat(sprintf("Excess kurtosis: %.3f\n", kurtosis_val))
cat(ifelse(abs(kurtosis_val) > 1,
           "→ Heavy tails confirmed: scat() JUSTIFIED ✓\n",
           "→ Near-normal distribution: scat() not necessary\n"))

# Diagnostic plots (4 panels)
par(mfrow = c(2, 2))
qqnorm(resid_gauss, main = "QQ-Plot of Residuals (Gaussian GAM)",
       pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3))
qqline(resid_gauss, col = "red", lwd = 2)

hist(resid_gauss, breaks = 80, freq = FALSE,
     main = "Residual Distribution", xlab = "Residuals", col = "lightblue")
curve(dnorm(x, mean = mean(resid_gauss), sd = sd(resid_gauss)),
      add = TRUE, col = "red", lwd = 2)

plot(fitted(gam_gaussian), resid_gauss,
     pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3),
     xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted(gam_gaussian), resid_gauss), col = "blue", lwd = 2)

aic_gauss <- AIC(gam_gaussian)
aic_scat  <- AIC(gam_full_corrected)
barplot(c(Gaussian = aic_gauss, `scat()` = aic_scat),
        col = c("lightblue", "salmon"),
        main = "AIC Comparison", ylab = "AIC (lower = better)")
par(mfrow = c(1, 1))

cat(sprintf("AIC Gaussian: %.1f\n", aic_gauss))
cat(sprintf("AIC scat():   %.1f\n", aic_scat))
cat(sprintf("ΔAIC = %.1f → %s\n", aic_gauss - aic_scat,
            ifelse(aic_gauss - aic_scat > 10,
                   "scat() strongly preferred ✓",
                   "marginal difference, Gaussian may suffice")))

# Methods section text
cat(sprintf(
  "\nMethods note: Residuals of the Gaussian GAM show excess kurtosis of %.2f
(Shapiro-Wilk p < %s). A scaled t-distribution family (scat()) was therefore
adopted for its robustness to heavy tails (ΔAIC = %.1f in favour of scat()).\n",
  kurtosis_val, format(sw_test$p.value, scientific = TRUE, digits = 2),
  aic_gauss - aic_scat
))


# --- 6.7  Date-only baseline model ---
form_date <- as.formula(paste(resp_var, "~ s(date_factor, bs='re')"))
gam_date  <- bam(form_date, data = df_stats_hp, family = scat(), method = "fREML", discrete = TRUE)
r2_date   <- summary(gam_date)$r.sq

# ==============================================================================
# --- 6.8  Leave-one-out variance decomposition ---
# ==============================================================================
# loo_r2      <- numeric(length(preds))
# names(loo_r2) <- c("LAI", "fCover", "Hmax",
#                    "FPC1 x Hmax",
#                    "FPC1", "FPC2", "FPC3")
# cat(sprintf("Running LOO variance decomposition for %s...\n", resp_var))
# for (i in seq_along(preds)) {
#   remaining  <- preds[-i]
#   form_red   <- as.formula(
#     paste(resp_var, "~", paste(remaining, collapse = " + "), "+ s(date_factor, bs='re')")
#   )
#   mod_red    <- bam(form_red, data = df_stats_hp, family = scat(), method = "fREML", discrete = TRUE)
#   loo_r2[i]  <- r2_full - summary(mod_red)$r.sq
# }
# 
# df_res <- data.frame(
#   Component = c("Date (Weather)", names(loo_r2), "Residual"),
#   R2_abs    = c(r2_date, loo_r2, 1 - r2_full)
# )
# df_res <- df_res %>%
#   mutate(
#     Pct = (R2_abs / sum(R2_abs)) * 100,
#     # On définit l'ordre d'affichage (du haut vers le bas dans le barplot)
#     Component = factor(Component, levels = rev(c(
#       "Date (Weather)", "LAI", "fCover", "Hmax", 
#       # "FPC1 x Hmax", 
#       "FPC1", "FPC2", "FPC3", "Residual"
#     )))
#   )

# --- 6.9  Interaction model: LAI × FPC1 (testing H2) ---
cat("\n--- Interaction Model (LAI × FPC1) ---\n") # Et Hmax ?

form_interact <- as.formula(paste(
  resp_var, "~",
  "s(LAI_sc) + s(FPC1_sc) + ti(LAI_sc, FPC1_sc) +",
  "s(fCover_sc) + s(Hmax_sc) +",
  "ti(FPC1_sc, Hmax_sc) +",
  "ti(FPC2_sc, Hmax_sc) +",
  "ti(FPC3_sc, Hmax_sc) +",
  "s(FPC2_sc) + s(FPC3_sc) +",
  "s(date_factor, bs='re')"
))

gam_interact     <- bam(form_interact, data = df_stats_hp,
                        family = scat(), method = "fREML", discrete = TRUE)
summary_interact <- summary(gam_interact)
print(summary_interact$s.table["ti(LAI_sc,FPC1_sc)", ])

# 3D surface visualizations
# par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
# 
# vis.gam(gam_interact,
#         view      = c("LAI_sc", "FPC1_sc"),
#         plot.type = "persp",
#         theta = 45, phi = 30, ticktype = "detailed",
#         color = "topo",
#         zlab  = "Predicted ΔTmax",
#         xlab  = "LAI (Scaled)", ylab = "FPC1 (Scaled)",
#         main  = "3D Interaction Surface: LAI vs Canopy Shape")
# 
# vis.gam(gam_interact,
#         view        = c("LAI_sc", "FPC1_sc"),
#         plot.type   = "contour",
#         color       = "topo", contour.col = "black",
#         xlab = "LAI (Scaled)", ylab = "FPC1 (Scaled)",
#         main = "Contour Map of Cooling Effect")
# points(df_stats_hp$LAI_sc, df_stats_hp$FPC1_sc, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.2))
# 
# par(mfrow = c(1, 1))

# Predict LAI effect across 3 typical FPC1 profiles (-1.2 SD, Mean, +1.2 SD)
pred_interact_2d <- ggpredict(
  gam_interact, 
  terms = c("LAI_sc [all]", "FPC1_sc [-1.2, 0, 1.2]")
) %>%
  as.data.frame() %>%
  # Back-transform LAI to its original scale for x-axis readability
  mutate(Real_LAI = (x * sd(df_daily_metrics$LAI)) + mean(df_daily_metrics$LAI))

# Create the plot with physically accurate labels
# FPC1 < 0 (-1.2 SD) = Bottom-Heavy (Yellow/Orange)
# FPC1 > 0 (+1.2 SD) = Top-Heavy (Blue)
p_interact_2d <- ggplot(pred_interact_2d, aes(x = Real_LAI, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(
    values = c("#fca50a", "#35b779", "#31688e"),
    labels = c("Bottom-Heavy (-1.2 SD)", "Average Shape (0)", "Top-Heavy (+1.2 SD)")
  ) +
  scale_fill_manual(
    values = c("#fca50a", "#35b779", "#31688e"),
    labels = c("Bottom-Heavy (-1.2 SD)", "Average Shape (0)", "Top-Heavy (+1.2 SD)")
  ) +
  labs(
    title    = "Synergistic Cooling: Interaction of LAI and Canopy Shape",
    subtitle = "Top-heavy canopies optimize cooling by facilitating sub-canopy ventilation",
    x        = "Total Foliage (Leaf Area Index)", 
    y        = "Predicted Cooling Effect (ΔTmax °C)",
    color    = "Canopy Architecture", 
    fill     = "Canopy Architecture"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = c(0.75, 0.8),
    legend.background = element_rect(color = "black", fill = "white"),
    plot.title        = element_text(face = "bold")
  )

print(p_interact_2d)

# 1. Create the detailed table with LAI steps of 0.5
df_lai_table <- pred_interact_2d %>%
  # Round the real LAI to the nearest 0.5 step
  mutate(LAI_step = round(Real_LAI * 2) / 2) %>%
  # Filter to keep a clean standard range
  filter(LAI_step >= 0.5 & LAI_step <= 8.0) %>%
  # Calculate the mean prediction for each rounded step and shape group
  group_by(LAI_step, group) %>%
  summarise(Predicted_DeltaT = mean(predicted), .groups = "drop") %>%
  # Pivot to wide format for columns
  pivot_wider(names_from = group, values_from = Predicted_DeltaT) %>%
  rename(
    `Bottom-Heavy` = `-1.2`,
    `Average`      = `0`,
    `Top-Heavy`    = `1.2`
  ) %>%
  # Calculate exact differences (Adding the Bottom - Average column as requested)
  mutate(
    `Bottom vs Avg` = `Bottom-Heavy` - `Average`,
    `Bottom vs Top` = `Bottom-Heavy` - `Top-Heavy`
  ) %>%
  # Round everything to 2 decimal places for clean display
  mutate(across(where(is.numeric), ~round(., 2)))

# 2. Reshape to long format to feed into ggplot2's geom_tile
df_plot_table <- df_lai_table %>%
  pivot_longer(cols = -LAI_step, names_to = "Column", values_to = "Value") %>%
  mutate(
    # Lock the exact column order from left to right
    Column = factor(Column, levels = c(
      "Bottom-Heavy", "Average", "Top-Heavy", "Bottom vs Avg", "Bottom vs Top"
    )),
    # Reverse LAI factor so low values are at the top (like a standard reading table)
    LAI_factor = factor(LAI_step, levels = sort(unique(LAI_step), decreasing = TRUE))
  )

# 3. Create the graphical table using ggplot2
p_table <- ggplot(df_plot_table, aes(x = Column, y = LAI_factor)) +
  # Add alternating background colors
  geom_tile(aes(fill = grepl("vs", Column)), color = "white", linewidth = 0.8) +
  scale_fill_manual(values = c("TRUE" = "#f8f9fa", "FALSE" = "#ffffff"), guide = "none") +
  
  # Conditional text coloring: black for raw values, blue/red for differences
  geom_text(aes(label = sprintf("%.2f", Value),
                color = ifelse(grepl("vs", Column), 
                               ifelse(Value < 0, "#d8576b", "#31688e"),
                               "black")), 
            size = 4.5, fontface = "bold") +
  scale_color_identity() + 
  
  # Clean up the axes
  scale_x_discrete(position = "top") +
  labs(
    title    = "Predicted \u0394Tmax (\u00B0C) and Aerodynamic Differences by LAI",
    subtitle = "Differences: Blue = Negative (Cooler) | Red = Positive (Warmer)",
    x = NULL, y = "Leaf Area Index (LAI)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x.top = element_text(face = "bold", size = 12, color = "black"),
    axis.text.y     = element_text(face = "bold", size = 12, color = "black"),
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40", margin = margin(b=15)),
    plot.margin     = margin(10, 10, 10, 10)
  )

print(p_table)

# ==============================================================================
# --- 3-PANEL INTERACTION (LAI x FPC1 across Real Canopy Heights) ---
# ==============================================================================

# 1. Retrieve the original mean and SD of canopy height (Hmax)
mean_h <- mean(df_daily_metrics$Hmax, na.rm = TRUE)
sd_h   <- sd(df_daily_metrics$Hmax, na.rm = TRUE)

# 2. Define the scaled steps used for the prediction
sd_steps <- c(-1.2, 0, 1.2)

# 3. Calculate the corresponding REAL heights in meters (rounded to 1 decimal)
real_heights <- round((sd_steps * sd_h) + mean_h, 1)

# 4. Create dynamic, physically meaningful labels for the facets
facet_labels <- c(
  sprintf("Short Canopy (~%.1f m)", real_heights[1]),
  sprintf("Average Height (~%.1f m)", real_heights[2]),
  sprintf("Tall Canopy (~%.1f m)", real_heights[3])
)

cat("Real heights used for facets:\n")
print(facet_labels)

# 5. Predict the 3-way interaction using the updated rigorous model
# (Assumes gam_interact already contains ti(FPC1_sc, Hmax_sc))
pred_interact_3way <- ggpredict(
  gam_interact, 
  terms = c("LAI_sc [all]", "FPC1_sc [-1.2, 0, 1.2]", "Hmax_sc [-1.2, 0, 1.2]")
) %>%
  as.data.frame() %>%
  mutate(
    Real_LAI = (x * sd(df_daily_metrics$LAI, na.rm = TRUE)) + mean(df_daily_metrics$LAI, na.rm = TRUE),
    Height_Group = factor(facet, levels = c("-1.2", "0", "1.2"), labels = facet_labels),
    Shape_Group = factor(group, levels = c("-1.2", "0", "1.2"), 
                         labels = c("Bottom-Heavy", "Average", "Top-Heavy"))
  )

# 6. Build the 3-panel plot
p_interact_3way <- ggplot(pred_interact_3way, aes(x = Real_LAI, y = predicted, color = Shape_Group, fill = Shape_Group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~ Height_Group) +
  scale_color_manual(values = c("#fca50a", "#35b779", "#31688e")) +
  scale_fill_manual(values = c("#fca50a", "#35b779", "#31688e")) +
  labs(
    title    = "Synergistic Cooling: Interaction of LAI and Canopy Shape by Height",
    x        = "Leaf Area Index (LAI)", 
    y        = "Predicted Cooling Effect (\u0394Tmax \u00B0C)",
    color    = "Canopy Architecture:", 
    fill     = "Canopy Architecture:"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position   = "bottom",
    plot.title        = element_text(face = "bold"),
    strip.background  = element_rect(fill = "grey90"),
    strip.text        = element_text(face = "bold", size = 12)
  )

print(p_interact_3way)

df_lai_table_3way <- pred_interact_3way %>%
  mutate(LAI_step = round(Real_LAI * 2) / 2) %>%
  filter(LAI_step >= 0.5 & LAI_step <= 8.0) %>%
  # On groupe directement par nos labels propres
  group_by(LAI_step, Shape_Group, Height_Group) %>%
  summarise(Predicted_DeltaT = mean(predicted), .groups = "drop") %>%
  pivot_wider(names_from = Shape_Group, values_from = Predicted_DeltaT) %>%
  # Plus besoin de rename(), les colonnes s'appellent déjà correctement
  mutate(
    `Bottom vs Avg` = `Bottom-Heavy` - `Average`,
    `Bottom vs Top` = `Bottom-Heavy` - `Top-Heavy`
  ) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

# Reshape pour le tableau graphique
df_plot_table_3way <- df_lai_table_3way %>%
  pivot_longer(
    cols = c(`Bottom-Heavy`, `Average`, `Top-Heavy`, `Bottom vs Avg`, `Bottom vs Top`), 
    names_to = "Column", 
    values_to = "Value"
  ) %>%
  mutate(
    Column = factor(Column, levels = c("Bottom-Heavy", "Average", "Top-Heavy", "Bottom vs Avg", "Bottom vs Top")),
    LAI_factor = factor(LAI_step, levels = sort(unique(LAI_step), decreasing = TRUE))
  )

# Plot du tableau
p_table_3way <- ggplot(df_plot_table_3way, aes(x = Column, y = LAI_factor)) +
  geom_tile(aes(fill = grepl("vs", Column)), color = "white", linewidth = 1.2) + # Lignes blanches plus épaisses
  scale_fill_manual(values = c("TRUE" = "#f1f3f5", "FALSE" = "#ffffff"), guide = "none") +
  
  # Augmenter légèrement la taille du texte (size = 4 ou 4.5)
  geom_text(aes(label = sprintf("%.2f", Value),
                color = ifelse(grepl("vs", Column), 
                               ifelse(Value < 0, "#31688e", "#d8576b"),
                               "black")), 
            size = 4.2, fontface = "bold") +
  
  scale_color_identity() + 
  facet_wrap(~ Height_Group) +
  scale_x_discrete(position = "top") +
  
  # labs(title = "Predicted \u0394Tmax (\u00B0C) and Aerodynamic Differences",
  #      subtitle = "Blue: Cooler | Red: Warmer",
  #      x = NULL, y = "LAI") +
  
  theme_minimal(base_size = 16) + # Augmenter la base_size globale
  theme(
    panel.grid      = element_blank(),
    # Améliorer la lisibilité des labels inclinés
    axis.text.x.top = element_text(face = "bold", size = 11, color = "black", 
                                   angle = 35, vjust = 0, hjust = 0),
    axis.text.y     = element_text(face = "bold", size = 12, color = "black"),
    strip.background= element_rect(fill = "grey20", color = NA), # Bandeau plus sombre pour contraster
    strip.text      = element_text(face = "bold", size = 13, color = "white"),
    plot.margin     = margin(20, 20, 20, 20)
  )

print(p_table_3way)

# ==============================================================================
# --- 6.10 Interaction Models: LAI x FPC2 and LAI x FPC3 ---
# ==============================================================================

# --- A. MODEL 2: Interaction LAI x FPC2 (Vertical Dispersion) ---
cat("\n--- Interaction Model (LAI x FPC2) ---\n")
form_interact_fpc2 <- as.formula(paste(
  resp_var, "~",
  "s(LAI_sc, k=5) + s(FPC2_sc, k=5) + ti(LAI_sc, FPC2_sc, k=5) +",
  "s(fCover_sc, k=5) + s(Hmax_sc, k=5) + s(FPC1_sc, k=5) + s(FPC3_sc, k=5) +",
  "s(date_factor, bs='re')"
))

gam_int_fpc2 <- bam(form_interact_fpc2, data = df_stats_hp, 
                    family = scat(), method = "fREML", discrete = TRUE)

# Print significance of the interaction term
print(summary(gam_int_fpc2)$s.table["ti(LAI_sc,FPC2_sc)", ])

# Predict for FPC2: Unimodal (-1.2 SD) vs Bimodal (+1.2 SD)
pred_fpc2_2d <- ggpredict(gam_int_fpc2, terms = c("LAI_sc [all]", "FPC2_sc [-1.2, 0, 1.2]")) %>%
  as.data.frame() %>%
  mutate(Real_LAI = (x * sd(df_daily_metrics$LAI)) + mean(df_daily_metrics$LAI))

p_int_fpc2 <- ggplot(pred_fpc2_2d, aes(x = Real_LAI, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("#d8576b", "#35b779", "#31688e"),
                     labels = c("Unimodal (-1.2 SD)", "Average Shape (0)", "Bimodal (+1.2 SD)")) +
  scale_fill_manual(values = c("#d8576b", "#35b779", "#31688e"),
                    labels = c("Unimodal (-1.2 SD)", "Average Shape (0)", "Bimodal (+1.2 SD)")) +
  labs(title = "Interaction: LAI and Vertical Dispersion (FPC2)",
       x = "Leaf Area Index", y = "Predicted \u0394Tmax (\u00B0C)",
       color = "FPC2 Shape:", fill = "FPC2 Shape:") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")


# --- B. MODEL 3: Interaction LAI x FPC3 (Micro-Complexity) ---
cat("\n--- Interaction Model (LAI x FPC3) ---\n")
form_interact_fpc3 <- as.formula(paste(
  resp_var, "~",
  "s(LAI_sc, k=5) + s(FPC3_sc, k=5) + ti(LAI_sc, FPC3_sc, k=5) +",
  "s(fCover_sc, k=5) + s(Hmax_sc, k=5) + s(FPC1_sc, k=5) + s(FPC2_sc, k=5) +",
  "s(date_factor, bs='re')"
))

gam_int_fpc3 <- bam(form_interact_fpc3, data = df_stats_hp, 
                    family = scat(), method = "fREML", discrete = TRUE)

# Print significance of the interaction term
print(summary(gam_int_fpc3)$s.table["ti(LAI_sc,FPC3_sc)", ])

# Predict for FPC3: Simple (-1.2 SD) vs Complex (+1.2 SD)
# Note: You might want to use larger SDs here (e.g., 2.0) if FPC3 variations are subtle
pred_fpc3_2d <- ggpredict(gam_int_fpc3, terms = c("LAI_sc [all]", "FPC3_sc [-1.2, 0, 1.2]")) %>%
  as.data.frame() %>%
  mutate(Real_LAI = (x * sd(df_daily_metrics$LAI)) + mean(df_daily_metrics$LAI))

p_int_fpc3 <- ggplot(pred_fpc3_2d, aes(x = Real_LAI, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("#fca50a", "#35b779", "#7a0403"),
                     labels = c("Simple (-1.2 SD)", "Average Shape (0)", "Complex (+1.2 SD)")) +
  scale_fill_manual(values = c("#fca50a", "#35b779", "#7a0403"),
                    labels = c("Simple (-1.2 SD)", "Average Shape (0)", "Complex (+1.2 SD)")) +
  labs(title = "Interaction: LAI and Micro-Complexity (FPC3)",
       x = "Leaf Area Index", y = "Predicted \u0394Tmax (\u00B0C)",
       color = "FPC3 Shape:", fill = "FPC3 Shape:") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")


# --- C. Display the plots side-by-side ---
p_combined_interactions <- p_int_fpc2 + p_int_fpc3 + 
  plot_annotation(title = "Synergistic Effects of Total Foliage and Sub-Canopy Architecture")

print(p_combined_interactions)

# ==============================================================================
# SECTION 7. VISUALIZATIONS: MAIN RESULTS
# ==============================================================================

# --- 7.1  Variance decomposition stacked bar ---
p_decomp <- ggplot(df_res, aes(x = 1, y = Pct, fill = Component)) +
  geom_col(position = "stack", color = "white", linewidth = 0.5) +
  geom_text(
    # Only show labels for components > 2% to avoid overlapping
    data = df_res %>% filter(Pct > 2.0), 
    aes(label = sprintf("%s\n%.1f%%", Component, Pct)),
    position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 3
  ) +
  scale_fill_viridis_d(option = "viridis", direction = -1) + # Palette auto-harmorisée
  coord_flip() +
  theme_void() +
  labs(title = "Variance Decomposition of ΔTmax")
print(p_decomp)


# --- 7.2  Partial effect curves: LAI saturation & FPC1 ---
pred_lai <- ggpredict(gam_full, terms = "LAI_sc [all]", type = "fixed") %>%
  as.data.frame() %>%
  mutate(Real_Value = (x * sd(df_daily_metrics$LAI)) + mean(df_daily_metrics$LAI))

p_sat <- ggplot(pred_lai, aes(x = Real_Value, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#31688e") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#31688e") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "A. LAI Saturation", x = "Leaf Area Index (LAI)", y = "Cooling Effect (°C)")

pred_fpc1 <- ggpredict(gam_full, terms = "FPC1_sc [all]") %>%
  as.data.frame() %>%
  mutate(Real_Value = (x * sd(df_daily_metrics$FPC1)) + mean(df_daily_metrics$FPC1))

p_shape_fpc <- ggplot(pred_fpc1, aes(x = Real_Value, y = predicted)) +
  geom_line(linewidth = 1.2, color = "#fde724") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4, fill = "#fde724") +
  labs(
    title    = "B. Profile Shape Effect (FPC1)",
    subtitle = "Negative = Top-Heavy | Positive = Bottom-Heavy",
    x = "FPC1 Score", y = "Relative Impact on Cooling (°C)"
  )

print(p_sat + p_shape_fpc)

# --- 7.3 Marginal effects of FPC1, FPC2, FPC3 (Unscaled & Height-Classed) ---
cat("Computing marginal effects for Functional Principal Components (Real scale)...\n")

# 1. Calcul des vraies valeurs pour les labels de Hauteur (Hmax)
mean_h <- mean(df_daily_metrics$Hmax, na.rm = TRUE)
sd_h   <- sd(df_daily_metrics$Hmax, na.rm = TRUE)

h_labels <- c(
  sprintf("Short (~%.1f m)", (-1.2 * sd_h) + mean_h),
  sprintf("Average (~%.1f m)", (0 * sd_h) + mean_h),
  sprintf("Tall (~%.1f m)", (1.2 * sd_h) + mean_h)
)

# Extraction des moyennes/SD pour dé-scaler les FPCs
mean_fpc1 <- mean(df_daily_metrics$FPC1, na.rm = TRUE)
sd_fpc1   <- sd(df_daily_metrics$FPC1, na.rm = TRUE)

mean_fpc2 <- mean(df_daily_metrics$FPC2, na.rm = TRUE)
sd_fpc2   <- sd(df_daily_metrics$FPC2, na.rm = TRUE)

mean_fpc3 <- mean(df_daily_metrics$FPC3, na.rm = TRUE)
sd_fpc3   <- sd(df_daily_metrics$FPC3, na.rm = TRUE)


# ==============================================================================
# --- PANNEAU A : FPC1 x Hmax (gam_interact -> les courbes vont se croiser)
# ==============================================================================
pred_fpc1_3h <- ggpredict(gam_interact, terms = c("FPC1_sc [all]", "Hmax_sc [-1.2, 0, 1.2]")) %>%
  as.data.frame() %>%
  mutate(Real_FPC1 = (x * sd_fpc1) + mean_fpc1) # Back-transform FPC1

p_fpc1 <- ggplot(pred_fpc1_3h, aes(x = Real_FPC1, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  scale_fill_manual( values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  labs(
    title    = "A. Interaction of Shape and Height",
    subtitle = "Cooling effect of FPC1 depends on tree height (ti(Hmax, FPC1))",
    x = "FPC1 Score\n\u2190 Bottom-Heavy          Top-Heavy \u2192",
    y = "Marginal Effect on \u0394Tmax (\u00B0C)",
    color = "Tree Height", fill = "Tree Height"
  ) +
  theme_bw(base_size = 16) + theme(legend.position = "right")


# ==============================================================================
# --- PANNEAU B : FPC2 pur décliné sur 3 hauteurs (gam_full -> courbes parallèles)
# ==============================================================================
pred_fpc2_3h <- ggpredict(gam_full, terms = c("FPC2_sc [all]", "Hmax_sc [-1.2, 0, 1.2]")) %>%
  as.data.frame() %>%
  mutate(Real_FPC2 = (x * sd_fpc2) + mean_fpc2)

p_fpc2 <- ggplot(pred_fpc2_3h, aes(x = Real_FPC2, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  scale_fill_manual( values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  labs(
    title    = "B. Effect of Vertical Dispersion (FPC2)",
    subtitle = "Additive effect across height classes",
    x = "FPC2 Score\n\u2190 Unimodal          Bimodal \u2192",
    y = "Effect on \u0394Tmax (\u00B0C)",
    color = "Tree Height", fill = "Tree Height"
  ) +
  theme_bw(base_size = 16) + theme(legend.position = "none")


# ==============================================================================
# --- PANNEAU C : FPC3 pur décliné sur 3 hauteurs (gam_full -> courbes parallèles)
# ==============================================================================
pred_fpc3_3h <- ggpredict(gam_full, terms = c("FPC3_sc [all]", "Hmax_sc [-1.2, 0, 1.2]")) %>%
  as.data.frame() %>%
  mutate(Real_FPC3 = (x * sd_fpc3) + mean_fpc3)

p_fpc3 <- ggplot(pred_fpc3_3h, aes(x = Real_FPC3, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  scale_fill_manual( values = c("#fca50a", "#35b779", "#31688e"), labels = h_labels) +
  labs(
    title    = "C. Effect of Micro-Complexity (FPC3)",
    subtitle = "Additive effect across height classes",
    x = "FPC3 Score\n\u2190 Simple          Complex \u2192",
    y = "Effect on \u0394Tmax (\u00B0C)",
    color = "Tree Height", fill = "Tree Height"
  ) +
  theme_bw(base_size = 16) + theme(legend.position = "none")


# ==============================================================================
# --- ASSEMBLAGE DES PANNEAUX
# ==============================================================================
plot_final <- p_fpc1 / (p_fpc2 + p_fpc3)
print(plot_final)




# 1. We include ALL structural predictors now
all_preds <- c("LAI_sc", "fCover_sc", "Hmax_sc", 
               "FPC1_sc", "FPC2_sc", "FPC3_sc")

# 2. Revised LOO Loop including FPC1
loo_r2_all <- numeric(length(all_preds))
names(loo_r2_all) <- c("LAI", "fCover", "Height",
                       "FPC1", "FPC2", "FPC3")

for (i in seq_along(all_preds)) {
  remaining <- all_preds[-i]
  form_red  <- as.formula(
    paste("Delta_Tmax ~", paste(remaining, collapse = " + "), "+ s(date_factor, bs='re')")
  )
  
  mod_red     <- bam(form_red, data = df_stats_hp, family = scat(), method = "fREML", discrete = TRUE)
  loo_r2_all[i] <- max(0, r2_full - summary(mod_red)$r.sq)
}

# 3. Building the dataframe
df_res_complete <- data.frame(
  Component = c("Weather (Date)", names(loo_r2_all), "Residual"),
  R2_abs    = c(r2_date, loo_r2_all, 1 - r2_full)
) %>%
  mutate(
    Pct = (R2_abs / sum(R2_abs)) * 100,
    Component = factor(Component, levels = rev(c(
      "Weather (Date)", "LAI", "fCover", "Height", 
      "FPC1", "FPC2", "FPC3", 
      "Residual"
    )))
  )

# 4. Plotting with FPC1 included
p_variance_complete <- ggplot(df_res_complete, aes(x = 1, y = Pct, fill = Component)) +
  geom_col(position = "stack", color = "white", linewidth = 0.6) +
  
  # Note: Labels for FPC1/2/3 might overlap if Pct is too small (< 1%)
  # We use check_overlap = TRUE or a filter for labels
  geom_text(
    data = df_res_complete %>% filter(Pct > 1.0),
    aes(label = sprintf("%s\n%.1f%%", Component, Pct)),
    position = position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 3.5
  ) +
  
  scale_fill_manual(values = c(
    "Weather (Date)" = "#440154",
    "LAI"            = "#31688e",
    "fCover"         = "#35b779",
    "Height"         = "#90d743",
    "FPC1"           = "#fde724",
    "FPC2"           = "#fca50a",
    "FPC3"           = "#d8576b",
    "Residual"       = "grey80"
  )) +
  
  coord_flip() +
  labs(
    title    = "Complete Variance Decomposition of ΔTmax",
    # subtitle = "Including all Canopy Structural PCAs",
    x = NULL, y = "Percentage of Total Variance (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid      = element_blank(),
    axis.text.y     = element_blank(),
    legend.position = "none",
    plot.title      = element_text(face = "bold")
  )

print(p_variance_complete)

# ==============================================================================
# SECTION 8. HEATWAVE STRESS-TEST
# ==============================================================================

threshold_temp <- 30
df_stats_hw    <- df_stats_hp %>% left_join(df_macro, by = "date")

hot_days_count <- df_macro %>% filter(Tmax_macro >= threshold_temp) %>% nrow()
cat(sprintf("Days with Tmax ≥ %d°C: %d\n", threshold_temp, hot_days_count))

df_heatwave <- df_stats_hw %>%
  filter(Tmax_macro >= threshold_temp) %>%
  mutate(date_factor = droplevels(date_factor))

cat(sprintf("Heatwave observations: %d\n", nrow(df_heatwave)))

# 1. Full model on heatwave subset (STRICTLY matching the 'All Days' model)
preds_updated <- c(
  "s(LAI_sc, k=5)", "s(fCover_sc, k=5)", "s(Hmax_sc, k=5)", 
  "s(FPC1_sc, k=5)", "s(FPC2_sc, k=5)", "s(FPC3_sc, k=5)"
  # Removed ti(FPC1, Hmax) to match the new baseline
)

form_full_updated <- as.formula(paste(
  resp_var, "~",
  paste(preds_updated, collapse = " + "),
  "+ s(date_factor, bs='re')"
))

gam_hw_full <- bam(form_full_updated, data = df_heatwave,
                   family = scat(), method = "fREML", discrete = TRUE)
r2_hw_full  <- summary(gam_hw_full)$r.sq

gam_hw_date <- bam(form_date, data = df_heatwave,
                   family = scat(), method = "fREML", discrete = TRUE)
r2_hw_date  <- summary(gam_hw_date)$r.sq

# 2. LOO decomposition for heatwave model
loo_r2_hw        <- numeric(length(preds_updated))
# Naming MUST match the 'All Days' exact strings
names(loo_r2_hw) <- c("LAI", "fCover", "Height", "FPC1", "FPC2", "FPC3")

cat("Computing LOO R² for heatwave subset...\n")
for (i in seq_along(preds_updated)) {
  remaining <- preds_updated[-i]
  form_red  <- as.formula(
    paste(resp_var, "~", paste(remaining, collapse = " + "), "+ s(date_factor, bs='re')")
  )
  mod_red       <- bam(form_red, data = df_heatwave, family = scat(), method = "fREML", discrete = TRUE)
  loo_r2_hw[i]  <- max(0, r2_hw_full - summary(mod_red)$r.sq)
}

# 3. Build Heatwave dataframe
df_compare <- data.frame(
  Condition = "Heatwaves (> 30°C)",
  Component = c("Weather (Date)", names(loo_r2_hw), "Residual"),
  R2_abs    = c(r2_hw_date, loo_r2_hw, 1 - r2_hw_full)
) %>%
  mutate(
    Pct       = (R2_abs / sum(R2_abs)) * 100,
    Component = factor(Component, levels = rev(c(
      "Weather (Date)", "LAI", "fCover", "Height",
      "FPC1", "FPC2", "FPC3", "Residual"
    )))
  )

# 4. Merge using the CORRECT 'df_res_complete' from Section 7
df_res_all <- df_res_complete %>%
  mutate(Condition = "All Days")

df_final_plot <- bind_rows(df_res_all, df_compare) %>%
  mutate(Condition = factor(Condition, levels = c("All Days", "Heatwaves (> 30°C)")))

# 5. Comparison plot: All Days vs Heatwaves (LOO Variance)
p_compare_clear <- ggplot(df_final_plot, aes(x = Component, y = Pct, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(
    data = df_final_plot %>% filter(Pct > 1.0), # Hide tiny labels to avoid mess
    aes(label = sprintf("%.1f%%", Pct)),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    color = "black", 
    fontface = "bold", 
    size = 3.5
  ) +
  scale_fill_manual(values = c(
    "All Days"           = "grey70",
    "Heatwaves (> 30°C)" = "#d73027"
  )) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(df_final_plot$Pct) * 1.15)) +
  labs(
    title    = "Microclimate Drivers: All Weather vs Extreme Heat",
    subtitle = "Side-by-side comparison of variance explained by each component",
    x = "Model Component", 
    y = "Percentage of Variance Explained (%)",
    fill = "Climatic Condition"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title        = element_text(face = "bold"),
    legend.position   = "bottom",
    panel.grid.major.y = element_blank() 
  )

print(p_compare_clear)

# ------------------------------------------------------------------------------
# VARIABLE IMPORTANCE BY PERMUTATION (VIP)
# ------------------------------------------------------------------------------

# Wrapper de prédiction robuste
pred_wrapper_vip <- function(object, newdata) {
  if (!"plot_id" %in% colnames(newdata)) {
    newdata$plot_id <- df_stats_hp$plot_id[1] 
  }
  predict(object, newdata = newdata, type = "response", exclude = "s(plot_id)")
}

# Liste exacte des prédicteurs physiques et temporels
features_to_test <- c("LAI_sc", "fCover_sc", "Hmax_sc", 
                      "FPC1_sc", "FPC2_sc", "FPC3_sc", "date_factor")

cat("\nComputing VIP for ALL DAYS (RMSE & R²)...\n")
vi_all_rmse <- vi(gam_full_corrected, method = "permute", train = df_stats_hp,
                  target = resp_var, feature_names = features_to_test,
                  metric = "rmse", smaller_is_better = TRUE,
                  pred_wrapper = pred_wrapper_vip, nsim = 10)

vi_all_r2 <- vi(gam_full_corrected, method = "permute", train = df_stats_hp,
                target = resp_var, feature_names = features_to_test,
                metric = "rsq", smaller_is_better = FALSE,
                pred_wrapper = pred_wrapper_vip, nsim = 10)

cat("Computing VIP for HEATWAVES (RMSE & R²)...\n")
vi_hw_rmse <- vi(gam_hw_full, method = "permute", train = df_heatwave,
                 target = resp_var, feature_names = features_to_test,
                 metric = "rmse", smaller_is_better = TRUE,
                 pred_wrapper = pred_wrapper_vip, nsim = 10)

vi_hw_r2 <- vi(gam_hw_full, method = "permute", train = df_heatwave,
               target = resp_var, feature_names = features_to_test,
               metric = "rsq", smaller_is_better = FALSE,
               pred_wrapper = pred_wrapper_vip, nsim = 10)

# Fonction de nettoyage des noms harmonisée avec le LOO
clean_var_names <- function(df) {
  df %>% mutate(
    Variable = recode(Variable,
                      "date_factor" = "Weather (Date)", 
                      "LAI_sc"      = "LAI", 
                      "fCover_sc"   = "fCover",
                      "Hmax_sc"     = "Height", 
                      "FPC1_sc"     = "FPC1", 
                      "FPC2_sc"     = "FPC2", 
                      "FPC3_sc"     = "FPC3"
    ),
    Variable = factor(Variable, levels = rev(c(
      "Weather (Date)", "LAI", "fCover", "Height", "FPC1", "FPC2", "FPC3"
    )))
  )
}

# Fusion des données RMSE
df_vip_rmse <- bind_rows(
  as.data.frame(vi_all_rmse) %>% mutate(Condition = "All Days"),
  as.data.frame(vi_hw_rmse)  %>% mutate(Condition = "Heatwaves (> 30°C)")
) %>%
  clean_var_names() %>%
  mutate(Condition = factor(Condition, levels = c("All Days", "Heatwaves (> 30°C)")))

# Fusion des données R²
df_vip_r2 <- bind_rows(
  as.data.frame(vi_all_r2) %>% mutate(Condition = "All Days"),
  as.data.frame(vi_hw_r2)  %>% mutate(Condition = "Heatwaves (> 30°C)")
) %>%
  clean_var_names() %>%
  mutate(Condition = factor(Condition, levels = c("All Days", "Heatwaves (> 30°C)")))

# --- 1. Calculate Baseline Metrics for both models ---
# Extract R-squared from summaries
r2_all <- summary(gam_full_corrected)$r.sq
r2_hw  <- summary(gam_hw_full)$r.sq

# Calculate RMSE from residuals
rmse_all <- sqrt(mean(residuals(gam_full_corrected, type = "response")^2))
rmse_hw  <- sqrt(mean(residuals(gam_hw_full, type = "response")^2))

# --- 2. Update Plot A: RMSE VIP ---
p_compare_rmse_vip <- ggplot(df_vip_rmse, aes(x = Variable, y = Importance, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("All Days" = "grey70", "Heatwaves (> 30°C)" = "#d73027")) +
  theme_bw(base_size = 14) +
  labs(
    title    = "A. Shift in Predictor Importance (RMSE)",
    # Dynamically inject baseline RMSE metrics into the subtitle
    subtitle = sprintf("Baseline RMSE \u2192 All Days: %.2f \u00B0C | Heatwaves: %.2f \u00B0C", rmse_all, rmse_hw),
    x = NULL, 
    y = "Importance (Increase in RMSE)"
  ) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# --- 3. Update Plot B: R-squared VIP ---
p_compare_r2_vip <- ggplot(df_vip_r2, aes(x = Variable, y = Importance, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = c("All Days" = "grey70", "Heatwaves (> 30°C)" = "#d73027")) +
  theme_bw(base_size = 14) +
  labs(
    title    = "B. Shift in Predictor Importance (R\u00B2)",
    # Dynamically inject baseline R-squared metrics into the subtitle
    subtitle = sprintf("Baseline R\u00B2 \u2192 All Days: %.2f | Heatwaves: %.2f", r2_all, r2_hw),
    x = NULL, 
    y = "Importance (Drop in R\u00B2)"
  ) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# --- 4. Render the combined plot ---
print(p_compare_rmse_vip | p_compare_r2_vip)


# ==============================================================================
# SECTION 9. EMPIRICAL VALIDATION AGAINST IN-SITU HOBO SENSORS
# ==============================================================================

cat("1. Extracting LiDAR metrics at HOBO locations...\n")

sf_hobo <- st_read("in_files/data_Blois_utm31n.geojson", quiet = TRUE)

r_landscape_all <- c(
  r_metrics_20m[["LAI"]], r_metrics_20m[["Hmax"]], r_metrics_20m[["fCover"]],
  r_fpcs
)
names(r_landscape_all) <- c("LAI", "Hmax", "fCover", "FPC1", "FPC2", "FPC3")

df_hobo_spatial <- terra::extract(r_landscape_all, vect(sf_hobo))
df_hobo_spatial$id_plot <- sf_hobo$id_plot
df_hobo_spatial <- df_hobo_spatial %>% dplyr::select(-ID)

cat("2. Processing HOBO temperature time series...\n")

df_hobo_temp <- read.csv("in_files/Blois_data_temperature.csv") %>%
  mutate(
    datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    date     = as.Date(datetime)
  ) %>%
  filter(position_sensor == "a", date %in% date_seq) %>%
  group_by(id_plot, date) %>%
  summarise(Tmax_obs = max(t_hobo, na.rm = TRUE), .groups = "drop") %>%
  inner_join(df_macro, by = "date") %>%
  mutate(Delta_obs = Tmax_obs - Tmax_macro)

cat("3. Merging spatial and temporal data...\n")

df_hobo_real <- df_hobo_temp %>%
  inner_join(df_hobo_spatial, by = "id_plot") %>%
  drop_na(LAI, 
          # FPC1, 
          Delta_obs)

cat(sprintf("df_hobo_real: %d daily observations | %d unique plots\n",
            nrow(df_hobo_real), length(unique(df_hobo_real$id_plot))))


# --- Scale HOBO predictors using training-data parameters ---
df_hobo_test <- df_hobo_real %>%
  mutate(
    LAI_sc      = (LAI    - mean(df_daily_metrics$LAI,    na.rm = TRUE)) / sd(df_daily_metrics$LAI,    na.rm = TRUE),
    fCover_sc   = (fCover - mean(df_daily_metrics$fCover, na.rm = TRUE)) / sd(df_daily_metrics$fCover, na.rm = TRUE),
    Hmax_sc     = (Hmax   - mean(df_daily_metrics$Hmax,   na.rm = TRUE)) / sd(df_daily_metrics$Hmax,   na.rm = TRUE),
    FPC1_sc     = (FPC1   - mean(df_daily_metrics$FPC1,   na.rm = TRUE)) / sd(df_daily_metrics$FPC1,   na.rm = TRUE),
    FPC2_sc     = (FPC2   - mean(df_daily_metrics$FPC2,   na.rm = TRUE)) / sd(df_daily_metrics$FPC2,   na.rm = TRUE),
    FPC3_sc     = (FPC3   - mean(df_daily_metrics$FPC3,   na.rm = TRUE)) / sd(df_daily_metrics$FPC3,   na.rm = TRUE),
    date_factor = as.factor(date)
  ) %>%
  filter(date_factor %in% levels(df_stats_hp$date_factor))

cat("Predicting temperatures for HOBO sensors using the trained GAM...\n")
df_hobo_test$Delta_pred <- predict(gam_full, newdata = df_hobo_test, type = "response")

obs      <- df_hobo_test$Delta_obs
prd      <- df_hobo_test$Delta_pred
rmse_val <- sqrt(mean((obs - prd)^2, na.rm = TRUE))
r2_val   <- cor(obs, prd, use = "complete.obs")^2
bias_val <- mean(prd - obs, na.rm = TRUE)

cat("--------------------------------------------------\n")
cat("GAM IN-SITU VALIDATION METRICS (HOBO sensors)\n")
cat(sprintf("R²   : %.3f\n", r2_val))
cat(sprintf("RMSE : %.2f °C\n", rmse_val))
cat(sprintf("Bias : %.2f °C\n", bias_val))
cat("--------------------------------------------------\n")

p_validation <- ggplot(df_hobo_test, aes(x = Delta_obs, y = Delta_pred)) +
  geom_point(alpha = 0.5, color = "#287D8EFF") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(
    title    = "Validation of the Spatial GAM against HOBO Sensors",
    subtitle = sprintf("R² = %.2f | RMSE = %.2f °C | Bias = %.2f °C",
                       r2_val, rmse_val, bias_val),
    x = "Observed Cooling Effect (HOBO, °C)",
    y = "Predicted Cooling Effect (GAM, °C)"
  ) +
  coord_fixed(ratio = 1)
print(p_validation)

# ==============================================================================
# SECTION 9.B. THE "TRIANGLE OF VALIDATION": GAMM vs MuSICA vs HOBO
# ==============================================================================

cat("1. Extracting daily Tmax from MuSICA HOBO simulations...\n")
out_nc_dir_hobo <- "out_files/musica_hobo_ERA5_results"
nc_files_hobo <- list.files(out_nc_dir_hobo, pattern = "\\.nc$", full.names = TRUE)

df_musica_hobo <- map_df(nc_files_hobo, function(f) {
  filename <- basename(f)
  # Extract id_plot from "musica_out_HOBO_{id_plot}.nc"
  id_str <- str_extract(filename, "(?<=HOBO_).*(?=\\.nc)")
  
  nc <- try(nc_open(f), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  raw_data <- try(get_variable(nc, "Tair_z"), silent = TRUE)
  nc_close(nc)
  if (inherits(raw_data, "try-error") || is.null(raw_data)) return(NULL)
  
  raw_data %>%
    filter(nair == 1) %>%
    mutate(Tair_sim = Tair_z - 273.15) %>%
    filter(as.Date(time) %in% date_seq) %>%
    mutate(time = time - hours(2), date = as.Date(time)) %>%
    group_by(date) %>%
    summarise(Tmax_musica = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
    mutate(id_plot = id_str)
})

cat("2. Merging MuSICA predictions with Macroclimate and empirical HOBO data...\n")
# Calculate Delta_Tmax for MuSICA
df_musica_hobo <- df_musica_hobo %>%
  inner_join(df_macro, by = "date") %>%
  mutate(Delta_musica = Tmax_musica - Tmax_macro)

# Merge all three sources: Empirical (df_hobo_real), MuSICA (df_musica_hobo), and GAMM predictors
df_triangle <- df_hobo_real %>%
  filter(!(id_plot %in% ids_to_remove)) %>%
  inner_join(df_musica_hobo %>% dplyr::select(id_plot, date, Delta_musica), 
             by = c("id_plot", "date")) %>%
  # Scale predictors for GAMM exactly as training data parameters
  mutate(
    LAI_sc      = (LAI    - mean(df_daily_metrics$LAI,    na.rm = TRUE)) / sd(df_daily_metrics$LAI,    na.rm = TRUE),
    fCover_sc   = (fCover - mean(df_daily_metrics$fCover, na.rm = TRUE)) / sd(df_daily_metrics$fCover, na.rm = TRUE),
    Hmax_sc     = (Hmax   - mean(df_daily_metrics$Hmax,   na.rm = TRUE)) / sd(df_daily_metrics$Hmax,   na.rm = TRUE),
    FPC1_sc     = (FPC1   - mean(df_daily_metrics$FPC1,   na.rm = TRUE)) / sd(df_daily_metrics$FPC1,   na.rm = TRUE),
    FPC2_sc     = (FPC2   - mean(df_daily_metrics$FPC2,   na.rm = TRUE)) / sd(df_daily_metrics$FPC2,   na.rm = TRUE),
    FPC3_sc     = (FPC3   - mean(df_daily_metrics$FPC3,   na.rm = TRUE)) / sd(df_daily_metrics$FPC3,   na.rm = TRUE),
    date_factor = as.factor(date)
  ) %>%
  filter(date_factor %in% levels(df_stats_hp$date_factor))

cat("3. Predicting temperatures using the trained GAMM...\n")
# We use the base gam_full (without spatial random effects) for transferability
df_triangle$Delta_gamm <- predict(gam_full, newdata = df_triangle, type = "response")

# --- Metrics computation ---
# A) GAMM vs MuSICA (Emulator fidelity)
r2_gamm_musica   <- cor(df_triangle$Delta_musica, df_triangle$Delta_gamm, use = "complete.obs")^2
rmse_gamm_musica <- sqrt(mean((df_triangle$Delta_musica - df_triangle$Delta_gamm)^2, na.rm = TRUE))
bias_gamm_musica <- mean(df_triangle$Delta_gamm - df_triangle$Delta_musica, na.rm = TRUE)

# B) MuSICA vs HOBO (Mechanistic model vs Reality)
r2_musica_hobo   <- cor(df_triangle$Delta_obs, df_triangle$Delta_musica, use = "complete.obs")^2
rmse_musica_hobo <- sqrt(mean((df_triangle$Delta_obs - df_triangle$Delta_musica)^2, na.rm = TRUE))
bias_musica_hobo <- mean(df_triangle$Delta_musica - df_triangle$Delta_obs, na.rm = TRUE)

# C) GAMM vs HOBO (Statistical Model vs Reality)
r2_gamm_hobo   <- cor(df_triangle$Delta_obs, df_triangle$Delta_gamm, use = "complete.obs")^2
rmse_gamm_hobo <- sqrt(mean((df_triangle$Delta_obs - df_triangle$Delta_gamm)^2, na.rm = TRUE))
bias_gamm_hobo <- mean(df_triangle$Delta_gamm - df_triangle$Delta_obs, na.rm = TRUE)

cat("\n--- TRIANGLE OF VALIDATION METRICS ---\n")
cat("A. GAMM vs MuSICA (Emulator Fidelity):\n")
cat(sprintf("   R²   : %.3f | RMSE : %.2f °C | Bias : %.2f °C\n", r2_gamm_musica, rmse_gamm_musica, bias_gamm_musica))
cat("B. MuSICA vs In-Situ HOBO (Mechanistic Reality Gap):\n")
cat(sprintf("   R²   : %.3f | RMSE : %.2f °C | Bias : %.2f °C\n", r2_musica_hobo, rmse_musica_hobo, bias_musica_hobo))
cat("C. GAMM vs In-Situ HOBO (Statistical Reality Gap):\n")
cat(sprintf("   R²   : %.3f | RMSE : %.2f °C | Bias : %.2f °C\n", r2_gamm_hobo, rmse_gamm_hobo, bias_gamm_hobo))
cat("--------------------------------------\n")

# --- 3. Plotting the 3 facets of validation ---
# --- Theme for all panels ---
theme_tri <- theme_bw(base_size = 18) + 
  theme(plot.title = element_text(face = "bold", size = 18),
        panel.grid.minor = element_blank())

# --- Revised Stats Function (Positions adapted to scales) ---
add_stats_A <- function(r2, rmse, bias) {
  annotate("text", x = -1.5, y = 6, 
           label = sprintf(" R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", r2, rmse, bias),
           hjust = -0.1, vjust = 1.1, fontface = "bold", size = 8)
}

add_stats_BC <- function(r2, rmse, bias) {
  annotate("text", x = -7.5, y = 10.5, 
           label = sprintf("R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", r2, rmse, bias),
           hjust = 0, vjust = 1, fontface = "bold", size = 8)
}

# PANEL A: GAMM vs MuSICA (Automatic Scaling)
p_gamm_musica <- ggplot(df_triangle, aes(x = Delta_musica, y = Delta_gamm)) +
  geom_point(alpha = 0.2, color = "#35b779") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
  add_stats_A(r2_gamm_musica, rmse_gamm_musica, bias_gamm_musica) +
  labs(title = "A. GAMM vs MuSICA",
       x = "Simulated DeltaTmax (MuSICA, °C)", 
       y = "Predicted DeltaTmax (GAMM, °C)") +
  theme_tri # No coord_fixed or scale limits here = Automatic

# --- FIXED RANGE FOR B & C ---
fixed_lims <- c(-8, 12)

# PANEL B: MuSICA vs HOBO (Mechanistic Gap)
p_musica_hobo <- ggplot(df_triangle, aes(x = Delta_obs, y = Delta_musica)) +
  geom_point(alpha = 0.2, color = "#287D8EFF") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
  scale_x_continuous(limits = fixed_lims) +
  scale_y_continuous(limits = fixed_lims) +
  add_stats_BC(r2_musica_hobo, rmse_musica_hobo, bias_musica_hobo) +
  labs(title = "B. MuSICA vs HOBOs",
       x = "Observed DeltaTmax (HOBOs, °C)", 
       y = "Simulated DeltaTmax (MuSICA, °C)") +
  theme_tri

# PANEL C: GAMM vs HOBO (Statistical Gap)
p_gamm_hobo <- ggplot(df_triangle, aes(x = Delta_obs, y = Delta_gamm)) +
  geom_point(alpha = 0.2, color = "#440154") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
  geom_smooth(method = "lm", color = "#d73027", se = TRUE) +
  scale_x_continuous(limits = fixed_lims) +
  scale_y_continuous(limits = fixed_lims) +
  add_stats_BC(r2_gamm_hobo, rmse_gamm_hobo, bias_gamm_hobo) +
  labs(title = "C. GAMM vs HOBOs",
       x = "Observed DeltaTmax (HOBOs, °C)", 
       y = "Predicted DeltaTmax (GAMM, °C)") +
  theme_tri

# --- Final Layout ---
(p_gamm_musica | (p_musica_hobo / p_gamm_hobo)) + 
  plot_annotation(theme = theme(plot.title = element_text(size = 18, 
                                                          face = "bold", 
                                                          hjust = 0.5)))

# ==============================================================================
# MODEL SUMMARY
# ==============================================================================

cat("\n\n=== MODEL COMPARISON SUMMARY ===\n")
cat(sprintf("%-45s R² = %.3f\n", "1. Original (no plot_id RE)",    r2_full))
cat(sprintf("%-45s R² = %.3f\n", "2. Corrected (+ plot_id RE)",    r2_full_corrected))
cat(sprintf("%-45s R² = %.3f\n", "3. Without Hmax (sensitivity)",  summary(gam_no_hmax)$r.sq))
cat("\nRecommendation: use Model 2 (corrected) as the primary model;\n")
cat("use Model 3 as a sensitivity test if Hmax/FPC concurvity > 0.6.\n")