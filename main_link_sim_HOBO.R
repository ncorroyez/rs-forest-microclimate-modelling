rm(list=ls(all=TRUE))
gc()

# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

# Install missing packages if needed
# install.packages(c("dplyr", "tidyr", "ggplot2", "lubridate", "viridis", "lme4", "ncdf4", "sf", "musica.tools"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(viridis)
library(lme4)
library(ncdf4)
library(musica.tools)
library(tibble)
library(ggrepel)
library(patchwork)

# Create output directories
if (!dir.exists("Plots/Validation")) dir.create("Plots/Validation", recursive = TRUE)
if (!dir.exists("03_RESULTS")) dir.create("03_RESULTS", recursive = TRUE)

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

# Study Period
date_seq <- seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day")

# Load Real LiDAR Profiles
# Assumes columns: id_plot, h_2.5, h_3.5, etc.
df_lad_raw <- read.csv("in_files/allometry_profiles.csv")
valid_plots <- unique(df_lad_raw$id_plot)

# Simulation Scenarios Paths
scenarios <- list()
for (i in 1:8) scenarios[[paste0("LAI_", i)]] <- file.path("out_files/TS/lidar", paste0("lai_", i))
scenarios[["Lidar_Constant"]] <- "out_files/TS/lidar/lidar_constant"

# ==============================================================================
# 3. DATA LOADING
# ==============================================================================

# A. HOBO OBSERVATIONS
# ------------------------------------------------------------------------------
cat("Loading HOBO Observation Data...\n")
site <- "Blois"
temperature_csv <- read.csv(file.path("in_files", paste0(site, "_data_temperature.csv")))

df_obs <- temperature_csv %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a") %>%
  filter(as.Date(datetime) %in% date_seq) %>%
  filter(id_plot %in% valid_plots) %>%
  group_by(id_plot, time = floor_date(datetime, "hour")) %>%
  summarise(Tair_obs = mean(t_hobo, na.rm = TRUE), .groups = "drop")

# B. MUSICA SIMULATIONS
# ------------------------------------------------------------------------------
cat("Loading MuSICA Simulations...\n")

df_sim_list <- list()

for (scen_name in names(scenarios)) {
  scen_path <- scenarios[[scen_name]]
  nc_files <- list.files(scen_path, pattern = "\\.nc$", full.names = TRUE)
  
  for (nc_file in nc_files) {
    fname <- basename(nc_file)
    plot_id <- valid_plots[sapply(valid_plots, function(x) grepl(x, fname))][1]
    
    if (!is.na(plot_id)) {
      nc <- tryCatch(nc_open(nc_file), error = function(e) NULL)
      if (!is.null(nc)) {
        raw_t <- tryCatch(get_variable(nc, "Tair_z"), error = function(e) NULL)
        if (!is.null(raw_t)) {
          df_chunk <- raw_t %>%
            filter(nair == 1) %>% # Ground level
            mutate(
              time = as.POSIXct(time, tz="UTC"),
              Tair_sim = Tair_z - 273.15,
              id_plot = plot_id,
              Scenario = scen_name
            ) %>%
            filter(as.Date(time) %in% date_seq) %>%
            select(time, id_plot, Scenario, Tair_sim)
          
          df_sim_list[[paste(scen_name, plot_id)]] <- df_chunk
        }
        nc_close(nc)
      }
    }
  }
}
df_sim_all <- bind_rows(df_sim_list)

# C. MERGE DATASETS (Time Alignment)
# ------------------------------------------------------------------------------
cat("Merging Datasets...\n")
# Align MuSICA (H:30) to HOBO (H:00)
df_sim_all_fixed <- df_sim_all %>% mutate(time = floor_date(time, unit = "hour"))

df_val <- inner_join(df_obs, df_sim_all_fixed, by = c("id_plot", "time")) %>%
  mutate(
    Error = Tair_sim - Tair_obs,
    Month = lubridate::month(time, label = TRUE)
  )

df_lad_raw <- df_lad_raw %>%
  filter(id_plot %in% unique(df_val$id_plot))

# ==============================================================================
# 4. PROFILE PROCESSING (METRICS & CLASSIFICATION)
# ==============================================================================
cat("--- STEP 4: Advanced Profile Analysis ---\n")

# 1. Pivot Wide to Long
df_lad_long <- df_lad_raw %>%
  select(id_plot, starts_with("h_")) %>%
  pivot_longer(cols = starts_with("h_"), names_to = "height_str", values_to = "lad") %>%
  mutate(z = as.numeric(gsub("h_", "", height_str))) %>%
  filter(!is.na(lad)) %>%
  arrange(id_plot, z)

# 2. Compute Metrics on TRIMMED Profiles
# (We cut the profile above the actual tree height to get correct relative shape)
plot_metrics_full <- df_lad_long %>%
  group_by(id_plot) %>%
  mutate(
    # Detect Canopy Height (last non-zero point)
    H_canopy_raw = max(z[lad > 1e-7], na.rm = TRUE)
  ) %>%
  # CRITICAL: Remove empty space above tree
  filter(z <= H_canopy_raw) %>%
  mutate(
    # Relative Coordinates (0=Ground, 1=Top)
    z_norm = (z - min(z)) / (max(z) - min(z)),
    lad_rel = lad / sum(lad)
  ) %>%
  summarise(
    # Absolute Metric
    H_canopy = unique(H_canopy_raw),
    
    # Relative Peak Position
    Peak_Height_Rel = z_norm[which.max(lad_rel)],
    
    # Vertical Mass Distribution
    Mass_Bot = sum(lad_rel[z_norm <= 1/3]),
    Mass_Mid = sum(lad_rel[z_norm > 1/3 & z_norm <= 2/3]),
    Mass_Top = sum(lad_rel[z_norm > 2/3]),
    
    # Heterogeneity
    LAD_SD = sd(lad_rel),
    
    .groups = "drop"
  )

# 3. Define Classification Function
classify_profiles <- function(metrics, 
                              LAD_SD_thr   = 0.015, 
                              Peak_top_thr = 2/3, 
                              Mass_thr     = 0.20) {
  metrics %>%
    mutate(
      Calculated_Shape = case_when(
        H_canopy < 11 ~ "bottom_heavy",
        # 1. Uniform
        LAD_SD < LAD_SD_thr ~ "uniform",
        # 2. Bimodal (Hollow Middle)
        Mass_Bot > Mass_thr & Mass_Top > Mass_thr & 
          # Mass_Mid < Mass_Bot & Mass_Mid < Mass_Top ~ "top_bottom_heavy",
          ((Mass_Bot + Mass_Top) / 2 ) > Mass_Mid ~ "top_bottom_heavy",
        # 3. Top Heavy
        Peak_Height_Rel > Peak_top_thr ~ "top_heavy",
        # 4. Bottom Heavy
        Peak_Height_Rel < (1 - Peak_top_thr) ~ "bottom_heavy",
        # 5. Intermediate
        TRUE ~ "symmetric"
      )
    )
}

# ==============================================================================
# 5. SENSITIVITY ANALYSIS
# ==============================================================================
cat("--- STEP 5: Sensitivity Analysis ---\n")

threshold_grid <- expand.grid(
  LAD_SD_thr   = c(0.012, 0.015, 0.018),
  Peak_top_thr = c(0.55, 0.60, 0.65),
  Mass_thr     = c(0.16, 0.20, 0.24)
)

sens_results <- vector("list", nrow(threshold_grid))

for (i in seq_len(nrow(threshold_grid))) {
  thr <- threshold_grid[i, ]
  
  # Temporary Classification
  tmp_metrics <- classify_profiles(plot_metrics_full, 
                                   thr$LAD_SD_thr, thr$Peak_top_thr, thr$Mass_thr)
  
  # Merge with Error Data
  tmp_df <- df_val %>%
    filter(Scenario == "Lidar_Constant") %>%
    left_join(tmp_metrics %>% select(id_plot, Calculated_Shape), by = "id_plot") %>%
    group_by(Calculated_Shape) %>%
    summarise(
      Bias = mean(Error, na.rm=TRUE),
      RMSE = sqrt(mean(Error^2, na.rm=TRUE)),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(Run_ID = i)
  
  sens_results[[i]] <- tmp_df
}

sens_df <- bind_rows(sens_results)
# Use sens_df to confirm robust trends if needed

# ==============================================================================
# 6. STATISTICAL MODELING
# ==============================================================================
cat("--- STEP 6: Statistical Modeling ---\n")

# Apply Final Classification (Standard Thresholds)
final_metrics <- classify_profiles(plot_metrics_full)
plot_metrics <- final_metrics # Update reference

# Prepare Data for Modeling
df_height <- df_val %>%
  filter(Scenario == "Lidar_Constant") %>%
  left_join(plot_metrics, by = "id_plot")

# A. LINEAR MODELS (Aggregated by Plot)
df_plot <- df_height %>%
  group_by(id_plot, Calculated_Shape, H_canopy) %>%
  summarise(
    Bias = mean(Error, na.rm=TRUE),
    RMSE = sqrt(mean(Error^2, na.rm=TRUE)),
    .groups = "drop"
  )

mod_shape_plot    <- lm(Bias ~ Calculated_Shape, data = df_plot)
mod_interact_plot <- lm(Bias ~ Calculated_Shape * H_canopy, data = df_plot)

cat("\n--- ANOVA (Interaction Significance) ---\n")
print(anova(mod_shape_plot, mod_interact_plot))

# B. MIXED MODELS (Full Hourly Data)
cat("\n--- Mixed Model (lmer) ---\n")
# Scale height for convergence
df_height$H_scaled <- scale(df_height$H_canopy)

mod_mixed <- lmer(
  Error ~ Calculated_Shape * H_scaled + (1 | id_plot),
  data = df_height,
  REML = TRUE
)
print(summary(mod_mixed))

# ==============================================================================
# 7. GENERATE PLOTS (1 to 8)
# ==============================================================================
cat("--- STEP 7: Visualization ---\n")

# PLOT 1: GLOBAL VALIDATION
# ------------------------------------------------------------------------------
df_val_const <- df_val %>% filter(Scenario == "Lidar_Constant")
metrics_global <- df_val_const %>% summarise(
  RMSE = sqrt(mean(Error^2)), Bias = mean(Error), 
  R2 = cor(Tair_sim, Tair_obs)^2
)

p1 <- ggplot(df_val_const, aes(x=Tair_obs, y=Tair_sim)) +
  geom_bin2d(bins=80) + scale_fill_viridis_c(option="magma", trans="log") +
  geom_abline(color="red", linetype="dashed") +
  annotate("text", x=1, y=35, label=paste0("RMSE: ", round(metrics_global$RMSE,2), "\nBias: ", round(metrics_global$Bias,2)), color="red", size=5, hjust=0) +
  labs(title="Global Validation", x="Observed T°C", y="Simulated T°C") + theme_bw()
ggsave("Plots/Validation/1_Global_Validation.png", p1, width=8, height=7)

# PLOT 2: PHENOLOGY
# ------------------------------------------------------------------------------
best_lai <- df_val %>% filter(Scenario!="Lidar_Constant") %>%
  mutate(LAI=as.numeric(gsub("LAI_", "", Scenario))) %>%
  group_by(id_plot, Month) %>% filter(sqrt(mean(Error^2))==min(sqrt(mean(Error^2)))) %>%
  summarise(Opt=mean(LAI), .groups="drop")

p2 <- ggplot(best_lai, aes(x=Month, y=Opt)) + geom_boxplot(fill="lightblue") + geom_smooth(aes(group=1), se=FALSE, color="red") +
  labs(title="Reconstructed Phenology", y="Optimal LAI") + theme_bw()
ggsave("Plots/Validation/2_Seasonal_Phenology.png", p2, width=8, height=6)

# PLOT 3: STRUCTURE BIAS (BOXPLOT)
# ------------------------------------------------------------------------------
p3 <- ggplot(df_plot, aes(x=Calculated_Shape, y=Bias, fill=Calculated_Shape)) +
  geom_hline(yintercept=0, linetype="dashed") + geom_boxplot(alpha=0.7) +
  scale_fill_viridis_d(option="turbo") +
  labs(title="Model Bias by Shape", y="Bias (°C)") + theme_bw()
ggsave("Plots/Validation/3_Structure_Bias.png", p3, width=8, height=6)

# PLOT 4: REALITY CHECK
# ------------------------------------------------------------------------------
df_reality <- df_val_const %>% left_join(plot_metrics, by="id_plot")
p4 <- ggplot(df_reality, aes(x=Calculated_Shape, y=Tair_obs, fill=Calculated_Shape)) +
  geom_boxplot(alpha=0.7) + facet_wrap(~Month) + scale_fill_viridis_d(option="turbo") +
  labs(title="Reality Check: Observed T°C", y="Observed T°C") + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Plots/Validation/4_Reality_Check.png", p4, width=10, height=8)

# PLOT 5: PROFILE COMPARISON (REAL + THEORETICAL)
# ------------------------------------------------------------------------------
cat("Generating Plot 5 (Profiles)...\n")
get_theo_curve <- function(n, type) {
  x <- seq(0,1,len=n)
  if(type=="uniform") d<-rep(1,n)
  else if(type=="top_heavy") d<-dbeta(x,6,2.5)
  else if(type=="bottom_heavy") d<-dbeta(x,2.5,6)
  else if(type=="symmetric") d<-dbeta(x,5,5)
  else d<-dbeta(x,4,20)+dbeta(x,20,4)
  return(d/sum(d))
}

# All
df_viz <- data.frame()
for(pid in unique(df_lad_long$id_plot)) {
  info <- plot_metrics %>% filter(id_plot==pid)
  if(nrow(info)==0) next
  shape <- info$Calculated_Shape
  h <- info$H_canopy
  real <- df_lad_long %>% filter(id_plot==pid) %>% arrange(z) %>% filter(z<=h)
  # if(nrow(real)<3) next
  
  theo <- get_theo_curve(nrow(real), shape)
  df_viz <- rbind(df_viz, data.frame(
    id_plot=pid, z=real$z, real_lad=real$lad/sum(real$lad), theo_lad=theo,
    Shape=shape, Label=paste0(pid,"\n",shape,"\n(H=",round(h,1),"m)")
  ))
}

p5 <- ggplot(df_viz) +
  # geom_area(aes(x=real_lad, y=z), fill="grey30", alpha=0.3) +
  geom_path(aes(x=real_lad, y=z), color="black", size=0.4) +
  geom_path(aes(x=theo_lad, y=z, color=Shape), size=1, linetype="dashed") +
  facet_wrap(~Label, scales="free") +
  scale_color_viridis_d(option="turbo") +
  labs(title="Structural Classification", x="Rel Density", y="Height (m)") +
  theme_bw() + theme(axis.text.x=element_blank(), strip.text=element_text(size=6))
ggsave("Plots/Validation/5_Profile_Comparison_All.png", p5, 
       width=20, height=15, dpi = 300)

# Per shape
shapes_found <- unique(df_viz$Shape)

for (shp in shapes_found) {
  sub_data <- df_viz %>% filter(Shape == shp)
  
  p_focus <- ggplot(sub_data) +
    geom_path(aes(x = real_lad, y = z), color = "black", size = 0.6) +
    geom_path(aes(x = theo_lad, y = z), color = "red", size = 1.2, linetype = "dashed") +
    
    facet_wrap(~Label, scales = "free") +
    labs(title = paste0("Focus: ", shp, " Profiles"), x = "Rel Density", y = "Height (m)") +
    theme_bw() + theme(axis.text.x = element_blank())
  
  # Dynamic height based on number of plots
  n_plots <- length(unique(sub_data$id_plot))
  h_calc <- max(5, ceiling(n_plots/4) * 3)
  
  ggsave(paste0("Plots/Validation/5_Profile_Focus_", shp, ".png"), 
         p_focus, width = 12, height = h_calc, dpi = 300)
}

# PLOT 6: INTERACTION SCATTER (CRITICAL)
# ------------------------------------------------------------------------------
p6 <- ggplot(df_plot, aes(x=H_canopy, y=Bias, color=Calculated_Shape)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(size=3, alpha=0.8) + geom_smooth(method="lm", se=FALSE, size=1) +
  scale_color_viridis_d(option="turbo") +
  labs(title="Height vs Shape Interaction", x="Canopy Height (m)", y="Mean Bias (°C)") + theme_bw()
ggsave("Plots/Validation/6_Height_Interaction.png", p6, width=9, height=7)

# ==============================================================================
# PLOT 6: INTERACTION SCATTER (WITH LABELS)
# ==============================================================================
cat("Generating Plot 6 with Labels...\n")

p6 <- ggplot(df_plot, aes(x = H_canopy, y = Bias, color = Calculated_Shape)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Trend Lines (Seulement si assez de points, sinon ça peut faire moche pour certaines classes)
  geom_smooth(method = "lm", se = FALSE, size = 1, alpha = 0.5) +
  
  # Points
  geom_point(size = 3, alpha = 0.8) +
  
  # Labels (Smart placement with ggrepel)
  geom_text_repel(
    aes(label = id_plot),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 20, # Increase if some labels are missing
    show.legend = FALSE
  ) +
  
  scale_color_viridis_d(option = "turbo") +
  labs(
    title = "Height vs Shape Interaction", 
    subtitle = "Labels indicate Plot ID",
    x = "Canopy Height (m)", 
    y = "Mean Bias (°C)"
  ) +
  theme_bw()

ggsave("Plots/Validation/6_Height_Interaction_Labels.png", p6, width = 10, height = 8)
print(p6)

# PLOT 7: FACET BY HEIGHT CLASS
# ------------------------------------------------------------------------------
df_plot_class <- df_plot %>% mutate(H_Class = cut(H_canopy, breaks=c(0,15,25,50), labels=c("Short","Medium","Tall")))
p7 <- ggplot(df_plot_class, aes(x=Calculated_Shape, y=Bias, fill=Calculated_Shape)) +
  geom_hline(yintercept=0, linetype="dashed") + geom_boxplot() + facet_wrap(~H_Class) +
  scale_fill_viridis_d(option="turbo") + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("Plots/Validation/7_Bias_by_Height_Class.png", p7, width=10, height=6)

# PLOT 8: CONTINUOUS ANALYSIS (Corrected)
# ------------------------------------------------------------------------------
# Safe join method: Aggregate error first, then join metrics
df_error_agg <- df_val_const %>% group_by(id_plot) %>% summarise(Bias=mean(Error), .groups="drop")
df_plot_cont <- df_error_agg %>% left_join(plot_metrics_full, by="id_plot")

p8 <- ggplot(df_plot_cont, aes(x=Peak_Height_Rel, y=Bias)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(aes(size=H_canopy), alpha=0.6) +
  geom_smooth(method="loess", se=TRUE, color="black") +
  theme_bw() +
  labs(title="Continuous Analysis: Peak vs Bias", x="Relative Peak Height (0-1)", y="Bias (°C)")
ggsave("Plots/Validation/8_Continuous_Analysis.png", p8, width=8, height=6)

cat("=== FULL ANALYSIS COMPLETE & SUCCESSFUL ===\n")







# ==============================================================================
# 9. K-MEANS CLUSTERING ANALYSIS (UNSUPERVISED)
# ==============================================================================
cat("\n--- STEP 9: K-Means Clustering Analysis ---\n")

# A. DATA PREPARATION (DECILES)
# ------------------------------------------------------------------------------
# We calculate vertical deciles for each plot to feed the clustering algorithm.
# This creates a standardized "shape signature" for each tree.

# 1. Prepare normalized profiles
df_profiles_norm <- df_lad_raw %>%
  select(id_plot, starts_with("h_")) %>%
  pivot_longer(cols = starts_with("h_"), names_to = "height_str", values_to = "lad") %>%
  mutate(z = as.numeric(gsub("h_", "", height_str))) %>%
  filter(!is.na(lad)) %>%
  group_by(id_plot) %>%
  mutate(H_canopy = max(z[lad > 1e-7], na.rm = TRUE)) %>% # Real height
  filter(z <= H_canopy) %>% # Keep only the tree
  mutate(
    z_rel = (z - min(z)) / (max(z) - min(z)), # Relative Height (0-1)
    lad_rel = lad / sum(lad) # Relative Density (Sum=1)
  ) %>%
  ungroup()

# 2. Aggregate by Decile (10 slices)
df_deciles <- df_profiles_norm %>%
  mutate(
    Decile = cut(z_rel, breaks = seq(0, 1, 0.1), labels = 1:10, include.lowest = TRUE)
  ) %>%
  group_by(id_plot, Decile) %>%
  summarise(LAD_Sum = sum(lad_rel), .groups = "drop") %>%
  # Fill missing deciles with 0
  complete(id_plot, Decile, fill = list(LAD_Sum = 0)) %>%
  group_by(id_plot) %>%
  mutate(LAD_Prop = LAD_Sum / sum(LAD_Sum)) %>% # Final renormalization
  ungroup()

# 3. Create Matrix for K-Means (Rows=Plots, Cols=Deciles)
df_matrix <- df_deciles %>%
  select(id_plot, Decile, LAD_Prop) %>%
  pivot_wider(names_from = Decile, values_from = LAD_Prop, names_prefix = "D") %>%
  column_to_rownames("id_plot")

# B. ELBOW METHOD (FIND OPTIMAL K)
# ------------------------------------------------------------------------------
cat("Running Elbow Method (k = 2 to 10)...\n")

set.seed(42)
k_values <- 2:10
wss <- sapply(k_values, function(k) {
  kmeans(df_matrix, centers = k, nstart = 25)$tot.withinss
})

# Plot Elbow
df_elbow <- data.frame(k = k_values, wss = wss)

p_elbow <- ggplot(df_elbow, aes(x = k, y = wss)) +
  geom_line(color = "blue", size = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = k_values) +
  labs(
    title = "Elbow Method for Optimal k",
    subtitle = "Look for the 'knee' in the curve (typically k=3, 4, or 5)",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares"
  ) +
  theme_bw()

# ggsave("Plots/Validation/9a_Elbow_Method.png", p_elbow, width = 6, height = 4)
print(p_elbow)

# C. RUN K-MEANS (WITH CHOSEN K)
# ------------------------------------------------------------------------------
# We force k=5 to compare with your manual classes, but you can change this 
# based on the Elbow plot result.
k_final <- 5
cat(paste0("Running final K-Means with k = ", k_final, "...\n"))

set.seed(42)
km_res <- kmeans(df_matrix, centers = k_final, nstart = 50)

# Extract Clusters
df_clusters <- data.frame(
  id_plot = names(km_res$cluster),
  Cluster = as.factor(km_res$cluster)
)

# D. VISUALIZE CLUSTER SHAPES (AVERAGE PROFILE)
# ------------------------------------------------------------------------------
# This tells you WHAT each cluster physically represents.

cluster_centers <- df_deciles %>%
  left_join(df_clusters, by = "id_plot") %>%
  group_by(Cluster, Decile) %>%
  summarise(Mean_LAD = mean(LAD_Prop), .groups = "drop") %>%
  mutate(Decile_Num = as.numeric(Decile))

p_shapes <- ggplot(cluster_centers, aes(x = Mean_LAD, y = Decile_Num, color = Cluster)) +
  geom_path(size = 1.5) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = 1:10, name = "Height Decile (1=Bottom, 10=Top)") +
  scale_x_continuous(name = "Relative Foliage Density") +
  scale_color_viridis_d(option = "turbo") +
  facet_wrap(~Cluster) +
  labs(
    title = "Average Profiles by Cluster (Decile Method)", 
    subtitle = "Identify physical shapes: Top Heavy (mass at 8-10) vs Bottom Heavy (mass at 1-3)"
  ) +
  theme_bw()
print(p_shapes)
# ggsave("Plots/Validation/9b_Cluster_Shapes.png", p_shapes, width = 10, height = 6) 

# E. PLOT 9: INTERACTION WITH CLUSTERS (Height vs Bias)
# ------------------------------------------------------------------------------
cat("Generating Plot 9 (Clusters vs Bias)...\n")

# Prepare data: Join Bias + Cluster + Height
# Need absolute height for X-axis
df_heights_ref <- df_profiles_norm %>% 
  group_by(id_plot) %>% 
  summarise(H_canopy = unique(H_canopy), .groups="drop")

df_plot_k <- df_val %>%
  filter(Scenario == "Lidar_Constant") %>%
  group_by(id_plot) %>%
  summarise(Bias = mean(Error, na.rm=TRUE), .groups="drop") %>%
  left_join(df_clusters, by = "id_plot") %>%
  left_join(df_heights_ref, by = "id_plot")

# Plot
p9 <- ggplot(df_plot_k, aes(x = H_canopy, y = Bias, color = Cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Trend lines per cluster
  geom_smooth(method = "lm", se = FALSE, size = 1, alpha = 0.5) +
  
  # Points with ID
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = id_plot), size = 3, max.overlaps = 20, show.legend = FALSE) +
  
  scale_color_viridis_d(option = "turbo", name = "Decile Cluster") +
  labs(
    title = "Interaction: Height vs Bias (Unsupervised Clustering)", 
    subtitle = "Do the unsupervised clusters explain the bias better than manual classes?",
    x = "Canopy Height (m)", 
    y = "Mean Bias (°C)"
  ) +
  theme_bw()

# ggsave("Plots/Validation/9c_Height_Interaction_Clusters.png", p9, width = 10, height = 8)
print(p9)

cat("=== K-MEANS ANALYSIS COMPLETE ===\n")

























# ==============================================================================
# 9-BIS. K-MEANS ON ABSOLUTE LAD (NO ROW NORMALIZATION)
# ==============================================================================
cat("\n--- STEP 9-BIS: Absolute K-Means (Magnitude + Shape) ---\n")

# A. DATA PREPARATION (ABSOLUTE DECILES)
# ------------------------------------------------------------------------------

# 1. Prepare profiles (Normalizing HEIGHT is still necessary to align deciles)
# But we do NOT normalize LAD sum to 1.
df_profiles_abs <- df_lad_raw %>%
  select(id_plot, starts_with("h_")) %>%
  pivot_longer(cols = starts_with("h_"), names_to = "height_str", values_to = "lad") %>%
  mutate(z = as.numeric(gsub("h_", "", height_str))) %>%
  filter(!is.na(lad)) %>%
  group_by(id_plot) %>%
  mutate(H_canopy = max(z[lad > 0.001], na.rm = TRUE)) %>% 
  filter(z <= H_canopy) %>% 
  mutate(
    z_rel = (z - min(z)) / (max(z) - min(z)) # Keep Relative Height for alignment
    # lad_rel calculation REMOVED: We keep raw 'lad'
  ) %>%
  ungroup()

# 2. Aggregate by Decile (Sum of Absolute LAD)
df_deciles_abs <- df_profiles_abs %>%
  mutate(
    Decile = cut(z_rel, breaks = seq(0, 1, 0.1), labels = 1:10, include.lowest = TRUE)
  ) %>%
  group_by(id_plot, Decile) %>%
  summarise(LAD_Sum = sum(lad), .groups = "drop") %>% # Sum of absolute LAD
  complete(id_plot, Decile, fill = list(LAD_Sum = 0)) %>%
  # CRITICAL CHANGE: NO division by sum(LAD_Sum) here.
  ungroup()

# 3. Create Matrix
df_matrix_abs <- df_deciles_abs %>%
  select(id_plot, Decile, LAD_Sum) %>%
  pivot_wider(names_from = Decile, values_from = LAD_Sum, names_prefix = "D") %>%
  column_to_rownames("id_plot")

# B. ELBOW METHOD (ABSOLUTE)
# ------------------------------------------------------------------------------
cat("Running Elbow Method on Absolute Data...\n")
set.seed(42)
wss <- sapply(2:10, function(k) kmeans(df_matrix_abs, centers = k, nstart = 25)$tot.withinss)

p_elbow_abs <- ggplot(data.frame(k=2:10, wss=wss), aes(x=k, y=wss)) +
  geom_line(color="darkred", size=1) + geom_point(size=3) +
  labs(title="Elbow Method (Absolute LAD)", subtitle="Determined by LAI magnitude + Shape", y="Within-SS") +
  theme_bw()
print(p_elbow_abs)
# ggsave("Plots/Validation/9a_Elbow_Absolute.png", p_elbow_abs, width=6, height=4)

# C. RUN K-MEANS
# ------------------------------------------------------------------------------
k_final <- 5 # Adjust based on new elbow if needed
set.seed(42)
km_res <- kmeans(df_matrix_abs, centers = k_final, nstart = 50)

df_clusters_abs <- data.frame(
  id_plot = names(km_res$cluster),
  Cluster_Abs = as.factor(km_res$cluster) # Rename to differentiate
)

# D. VISUALIZE SHAPES (ABSOLUTE QUANTITY)
# ------------------------------------------------------------------------------
cluster_centers_abs <- df_deciles_abs %>%
  left_join(df_clusters_abs, by = "id_plot") %>%
  group_by(Cluster_Abs, Decile) %>%
  summarise(Mean_LAD = mean(LAD_Sum), .groups = "drop") %>%
  mutate(Decile_Num = as.numeric(Decile))

# Calculate Average LAI per cluster to help interpretation
cluster_lai <- df_deciles_abs %>%
  left_join(df_clusters_abs, by="id_plot") %>%
  group_by(Cluster_Abs) %>%
  summarise(Avg_LAI = round(mean(tapply(LAD_Sum, id_plot, sum)), 2))

p_shapes_abs <- ggplot(cluster_centers_abs, aes(x = Mean_LAD, y = Decile_Num, color = Cluster_Abs)) +
  geom_path(size = 1.5) +
  scale_y_continuous(breaks = 1:10, name = "Height Decile") +
  scale_x_continuous(name = "Absolute LAD Sum (m2/m3)") +
  scale_color_viridis_d(option = "turbo") +
  facet_wrap(~Cluster_Abs) +
  labs(
    title = "Absolute Profiles by Cluster", 
    subtitle = "Clusters are now driven by Density (LAI) AND Shape"
  ) +
  geom_text(data=cluster_lai, aes(x=max(cluster_centers_abs$Mean_LAD)*0.7, y=1, label=paste("LAI~", Avg_LAI)), 
            inherit.aes=FALSE, size=3) +
  theme_bw()

# ggsave("Plots/Validation/9b_Cluster_Shapes_Absolute.png", p_shapes_abs, width = 10, height = 6)
print(p_shapes_abs)

# E. PLOT INTERACTION
# ------------------------------------------------------------------------------
# Join with Height & Bias
df_plot_abs <- df_val %>%
  filter(Scenario == "Lidar_Constant") %>%
  group_by(id_plot) %>%
  summarise(Bias = mean(Error, na.rm=TRUE), .groups="drop") %>%
  left_join(df_clusters_abs, by = "id_plot") %>%
  left_join(df_profiles_abs %>% group_by(id_plot) %>% summarise(H=max(H_canopy)), by="id_plot")

p9_abs <- ggplot(df_plot_abs, aes(x = H, y = Bias, color = Cluster_Abs)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "lm", se = FALSE, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label=id_plot), size=3, show.legend=FALSE) +
  scale_color_viridis_d(option = "turbo", name = "Decile Cluster") +
  labs(
    title = "Interaction: Height vs Bias (Absolute Clustering)", 
    subtitle = "Grouping by Absolute Biomass quantity",
    x = "Canopy Height (m)", y = "Bias (°C)"
  ) + theme_bw()

# ggsave("Plots/Validation/9c_Height_Interaction_Absolute.png", p9_abs, width = 10, height = 8)
print(p9_abs)

# ==============================================================================
# 10. COMPARISON PLOT: RELATIVE VS ABSOLUTE CLUSTERING
# ==============================================================================
cat("\n--- STEP 10: Generating Comparison Plot (Shape vs Density) ---\n")

p_shapes_clean <- p_shapes + 
  labs(subtitle = "A. Shape Driven (Relative Deciles)") +
  theme(legend.position = "bottom")

p_shapes_abs_clean <- p_shapes_abs + 
  labs(subtitle = "B. Density Driven (Absolute Deciles)") +
  theme(legend.position = "bottom")

# Assemblage côte à côte
combined_p_shapes <- (p_shapes_clean | p_shapes_abs_clean) +
  plot_annotation(
    title = "Comparison of Clustering Strategies on Mean Profiles",
    subtitle = "Left: Grouping by vertical shape (Profile). Right: Grouping by foliage quantity (Biomass).",
    tag_levels = 'A'
  )

# Sauvegarde
print(combined_p_shapes)

p9_clean <- p9 + 
  labs(subtitle = "A. Shape Driven (Relative Deciles)") +
  theme(legend.position = "bottom")

p9_abs_clean <- p9_abs + 
  labs(subtitle = "B. Density Driven (Absolute Deciles)") +
  theme(legend.position = "bottom")

# Assemblage côte à côte
combined_plot9 <- (p9_clean | p9_abs_clean) +
  plot_annotation(
    title = "Comparison of Clustering Strategies on Model Bias",
    subtitle = "Left: Grouping by vertical shape (Profile). Right: Grouping by foliage quantity (Biomass).",
    tag_levels = 'A'
  )

# Sauvegarde
print(combined_plot)

cat("=== COMPARISON PLOT SAVED ===\n")