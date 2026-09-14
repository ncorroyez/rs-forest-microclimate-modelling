plot_obs_vs_sim <- function(combined_df, slope_threshold = 0) {
  
  # Select relevant columns
  df_Tair <- df_Tair %>% dplyr::select(id_plot, time, Tair_z)
  df_HOBO <- df_HOBO %>% dplyr::select(id_plot, time, t_hobo)
  
  # Merge datasets
  combined <- merge(df_Tair, df_HOBO, by = c("id_plot", "time"))
  
  # Compute global statistics
  r_value <- cor(combined$t_hobo, combined$Tair_z)
  bias <- mean(combined$Tair_z - combined$t_hobo)
  nrmse <- sqrt(mean((combined$t_hobo - combined$Tair_z)^2)) / IQR(combined$t_hobo)
  mae <- mean(abs(combined$Tair_z - combined$t_hobo))
  
  # Compute slopes per id_plot
  combined <- combined %>%
    group_by(id_plot) %>%
    mutate(
      slope = coef(lm(t_hobo ~ Tair_z, data = cur_data()))["Tair_z"],
      equilibrium = coef(lm(t_hobo ~ Tair_z, data = cur_data()))["(Intercept)"] / (1 - slope),
      log_slope = log(abs(slope))
    ) %>%
    ungroup()
  
  # Split data
  combined_neg <- combined %>% filter(log_slope < slope_threshold)
  combined_pos <- combined %>% filter(log_slope > slope_threshold)
  
  # Helper function for plotting
  plot_data <- function(data, title_suffix) {
    r_val <- cor(data$t_hobo, data$Tair_z)
    bias_val <- mean(data$Tair_z - data$t_hobo)
    nrmse_val <- sqrt(mean((data$t_hobo - data$Tair_z)^2)) / IQR(data$t_hobo)
    mae_val <- mean(abs(data$Tair_z - data$t_hobo))
    
    ggplot(data, aes(x = t_hobo, y = Tair_z, color = log_slope)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      labs(
        x = "Hourly HOBO Sensors Measures (°C)",
        y = "Predicted hourly by MuSICA (°C)",
        title = paste("Hourly Predictions by MuSICA vs HOBO Sensors", title_suffix),
        subtitle = paste("R =", round(r_val, 2),
                         "| NRMSE =", round(nrmse_val, 2),
                         "| Bias =", round(bias_val, 2),
                         "| MAE =", round(mae_val, 2)
        ),
        color = "Log Slope"
      ) +
      annotate(
        "text",
        x = Inf, y = -Inf,
        label = paste("N =", n_distinct(data$id_plot), "plots"),
        hjust = 1.1, vjust = -1.1,
        color = "black",
        size = 5
      ) +
      scale_color_gradient2(low = "#128d84", high = "#e9c624") +
      xlim(c(0, 40)) +
      ylim(c(0, 40)) +
      theme_bw() +
      theme(legend.position = "bottom", legend.key.size = unit(1, "cm")) +
      coord_fixed(ratio = 1)
  }
  
  # Generate plots
  plot_neg <- plot_data(combined_neg, "(Slope < 0)")
  plot_pos <- plot_data(combined_pos, "(Slope > 0)")
  plot_all <- plot_data(combined, "")
  
  # Arrange plots
  grid.arrange(plot_neg, plot_pos, ncol = 2)
  
  return(plot_all)
}

get_equilibrium_stats <- function(data, threshold = 0) {
  library(dplyr)
  
  # Overall equilibrium stats (one value per id_plot)
  overall_stats <- data %>% 
    distinct(id_plot, equilibrium) %>%
    summarise(
      mean_equilibrium   = mean(equilibrium, na.rm = TRUE),
      median_equilibrium = median(equilibrium, na.rm = TRUE),
      sd_equilibrium     = sd(equilibrium, na.rm = TRUE)
    ) %>% 
    mutate(group = "overall")
  
  # Equilibrium stats for plots with log_slope > threshold (positive)
  pos_stats <- data %>%
    filter(log_slope > threshold) %>%
    distinct(id_plot, equilibrium) %>%
    summarise(
      mean_equilibrium   = mean(equilibrium, na.rm = TRUE),
      median_equilibrium = median(equilibrium, na.rm = TRUE),
      sd_equilibrium     = sd(equilibrium, na.rm = TRUE)
    ) %>%
    mutate(group = "pos")
  
  # Equilibrium stats for plots with log_slope < threshold (negative)
  neg_stats <- data %>%
    filter(log_slope < threshold) %>%
    distinct(id_plot, equilibrium) %>%
    summarise(
      mean_equilibrium   = mean(equilibrium, na.rm = TRUE),
      median_equilibrium = median(equilibrium, na.rm = TRUE),
      sd_equilibrium     = sd(equilibrium, na.rm = TRUE)
    ) %>%
    mutate(group = "neg")
  
  # Combine the results
  stats_all <- bind_rows(overall_stats, pos_stats, neg_stats)
  return(stats_all)
}

# Function to merge and process data
merge_and_process <- function(df_Tair, t_hobo) {
  df_merged <- merge(df_Tair, t_hobo, by = c("id_plot", "time"), all = TRUE) %>%
    group_by(time) %>%
    summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
    dplyr::select(-id_plot)
  
  df_merged <- merge(df_merged, era5_Tair_all, by = "time", all = TRUE)
  df_merged <- merge(df_merged, ign_Tair_all, by = "time", all = TRUE)
  colnames(df_merged) <- c("time", "Tair_z", "T_hobo", "Tair_era5", "Tair_ign")
  
  return(df_merged)
}

# Function to calculate daily max and min
calculate_daily_extremes <- function(df_merged, filter_month = "06") {
  df_daily_max <- df_merged %>%
    mutate(date = as.Date(time)) %>%
    filter(format(date, "%m") == filter_month) %>%
    group_by(date) %>%
    summarise(
      MuSICA_max = max(Tair_z, na.rm = TRUE),
      HOBO_max = max(T_hobo, na.rm = TRUE),
      ERA5_max = max(Tair_era5, na.rm = TRUE),
      CHS41_max = max(Tair_ign, na.rm = TRUE),
      .groups = "drop"
    )
  
  df_daily_min <- df_merged %>%
    mutate(date = as.Date(time)) %>%
    filter(format(date, "%m") == filter_month) %>%
    group_by(date) %>%
    summarise(
      MuSICA_min = min(Tair_z, na.rm = TRUE),
      HOBO_min = min(T_hobo, na.rm = TRUE),
      ERA5_min = min(Tair_era5, na.rm = TRUE),
      CHS41_min = min(Tair_ign, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(max = df_daily_max, min = df_daily_min))
}

# Function to plot the extremes (Max and Min)
plot_extremes <- function(df_daily_max, df_daily_min, title) {
  df_max_long <- df_daily_max %>%
    pivot_longer(cols = c(MuSICA_max, HOBO_max, ERA5_max, CHS41_max),
                 names_to = "Source", values_to = "Max_Temperature")
  
  df_min_long <- df_daily_min %>%
    pivot_longer(cols = c(MuSICA_min, HOBO_min, ERA5_min, CHS41_min),
                 names_to = "Source", values_to = "Min_Temperature")
  
  # Plot Max temperatures
  pmax <- ggplot(df_max_long, aes(x = date, y = Max_Temperature, color = Source)) +
    geom_line(linewidth = 1.5) +
    labs(title = paste(title, "Max Temperature"), y = "Temperature (°C)", x = "Date") +
    theme_minimal()
  
  # Plot Min temperatures
  pmin <- ggplot(df_min_long, aes(x = date, y = Min_Temperature, color = Source)) +
    geom_line(linewidth = 1.5) +
    labs(title = paste(title, "Min Temperature"), y = "Temperature (°C)", x = "Date") +
    theme_minimal()
  
  print(pmax)
  print(pmin)
}

# Function to analyze hourly data for a specific week
analyze_hourly_week <- function(df_merged, start_date, end_date) {
  df_week <- df_merged %>%
    rename(MuSICA = Tair_z, 
           HOBO = T_hobo, 
           ERA5 = Tair_era5, 
           CHS41 = Tair_ign) %>%
    filter(time >= start_date & time <= end_date)
  
  df_long_week <- df_week %>%
    pivot_longer(cols = c("MuSICA", "HOBO", "ERA5", "CHS41"), 
                 names_to = "Source", values_to = "Temperature")
  
  ggplot(df_long_week, aes(x = time, y = Temperature, color = Source)) +
    geom_line(linewidth = 1.5) +
    labs(title = "Hourly Temperature from June 10 to June 17, 2021",
         y = "Temperature (°C)", x = "Time") +
    theme_minimal()
}