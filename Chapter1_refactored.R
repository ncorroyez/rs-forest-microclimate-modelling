# ==============================================================================
# Chapter 1 — Sensitivity of forest microclimate to 3D canopy architecture
# Author: Nathan Corroyez
#
# Refactor philosophy
# -------------------
# The previous monolithic script mixed (1) structural typology, (2) MuSICA
# batch simulations, (3) statistical emulation via GAMM, (4) HOBO validation,
# and (5) a Sentinel-2 vs LiDAR comparison. Following supervisor feedback,
# this refactor reorganises the code around three principles:
#
#   1. Hypotheses must be TESTABLE, not just interpretable.
#      → H1 (hierarchy of LiDAR-derived inputs: LAI, Hmax, fCover, LAD)
#        is now answered by a SCENARIO-BASED MuSICA sensitivity analysis.
#      → H2 (role of vertical LAD shape) is answered FIRST by a simple
#        uniform-vs-real LAD contrast.
#      → H3 (Sentinel-2 vs LiDAR) is an OPTIONAL annex module.
#
#   2. Every MuSICA call goes through ONE function, parameterised by a
#      "scenario".
#
#   3. The GAMM is a downstream emulator on the reference (full-real) scenario,
#      integrating BOTH structural metrics and physical macroclimate drivers
#      (Wind, Rad, VPD) to open the "time" black box.
#
# Pipeline (see main() at the bottom):
#   Section 1.  Configuration & libraries
#   Section 2.  I/O helpers (including physical ERA5 macroclimate drivers)
#   Section 3.  Forest dataframe construction
#   Section 4.  cLHS sampling
#   Section 5.  LAD profile constructors
#   Section 6.  Scenario registry
#   Section 7.  MuSICA batch runner
#   Section 8.  Post-processing (NetCDF → daily ΔTmax)
#   Section 9.  H2 first-pass: uniform vs real LAD
#   Section 10. H1 sensitivity: forward / LOO
#   Section 11. Shape Metrics (FPCA & Height of Median LAI)
#   Section 12. GAMM emulator (structure + macroclimate)
#   Section 13. Diagnostics: per-HOBO time-series overlays
#   Section 14. HOBO validation
#   Section 15. Sentinel-2 annex (H3)
#   Section 16. main() orchestrator
# ==============================================================================

# ==============================================================================
# SECTION 1.  CONFIGURATION & LIBRARIES
# ==============================================================================

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

suppressPackageStartupMessages({
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
  library(terra)
  library(fda)
  library(corrplot)
  library(sf)
  library(clhs)
  library(vip)
  library(Metrics)
  library(tidyterra)
})

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

# ---- Paths -------------------------------------------------------------------
CFG <- list(
  in_dir          = "in_files",
  out_dir         = "out_files/Sensitivity_Analysis",
  forcing_file    = "in_files/musica_in_Blois.nc",
  hobo_geojson    = "in_files/data_Blois_utm31n.geojson",
  hobo_temp_csv   = "in_files/Blois_data_temperature.csv",
  musica_cmd      = "bash -i -c musica",
  date_seq        = seq(as.Date("2021-06-01"), as.Date("2021-09-30"), by = "day"),
  ids_to_remove   = c("41_13", "41_14", "41_20", "41_41", "41_50", "41_51", "41_53"),
  agg_factor      = 2,        # 10 m → 20 m aggregation
  n_per_cluster   = 100,      # cLHS sample size per cluster
  k_clusters      = NULL,     # K-means cluster assignment per plot
  heatwave_thr    = 30,       # °C threshold for heatwave subset
  out_h2          = "out_files/H2_uniform_vs_real",
  out_h1_forward  = "out_files/H1_forward",
  out_h1_loo      = "out_files/H1_loo",
  out_h1_factorial= "out_files/H1_factorial",
  out_hobo        = "out_files/musica_hobo_validation",
  out_s2_annex    = "out_files/Annex_S2_FORMSH"
)

# ---- Run flags ---------------------------------------------------------------
FLAGS <- list(
  RUN_CLHS              = TRUE,   # build cLHS sample
  RUN_H2_UNIFORM_VS_REAL= TRUE,   # H2 first pass
  RUN_H1_FORWARD        = TRUE,   # forward inclusion
  RUN_H1_LOO            = TRUE,   # leave-one-out

  RUN_GAMM_EMULATOR     = TRUE,   # Statistical emulator
  USE_MACROCLIMATE_GAMM = TRUE,   # Include physical meteorology in GAMM
  GAMM_SHAPE_VAR        = "CLUSTER",  # Shape predictor: "CLUSTER" (cat.), "FPC", "H_MEDIAN", "NONE"

  RUN_HOBO_VALIDATION   = TRUE,   # MuSICA at HOBO locations
  RUN_S2_ANNEX          = FALSE,  # Sentinel-2 annex (H3) — reporter à plus tard (cf. retours encadrants)
  RUN_H1_FACTORIAL      = FALSE   # heavy, off by default
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# SECTION 2.  I/O HELPERS
# ==============================================================================

#' Force UTC timezone on a NetCDF time variable.
force_utc_nc <- function(ncfile, varname = "time") {
  if (!file.exists(ncfile)) stop(paste("File not found:", ncfile))
  nc <- nc_open(ncfile); on.exit(nc_close(nc))
  time_val  <- ncvar_get(nc, varname)
  units_str <- ncatt_get(nc, varname, "units")$value
  parts     <- strsplit(units_str, " since ")[[1]]
  unit_type <- parts[1]
  origin_str<- gsub("\\s*\\(.*\\)| UTC", "", parts[2])
  origin_dt <- as.POSIXct(origin_str, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  mult <- switch(trimws(unit_type),
                 "hours" = 3600, "minutes" = 60,
                 "seconds" = 1,  "days"    = 86400, 3600)
  with_tz(origin_dt + (time_val * mult), "UTC")
}

#' Load and aggregate the LiDAR rasters.
load_lidar_rasters <- function(in_dir, agg_factor = 2) {
  r_lai     <- rast(file.path(in_dir, "lai_z1_res_10_m.tif"))
  r_vci     <- rast(file.path(in_dir, "vci_res_10_m.tif"))
  r_hmax    <- rast(file.path(in_dir, "max_res_10_m.tif"))
  r_fcover  <- rast(file.path(in_dir, "fCover_res_10_m.tif"))
  r_lad_st  <- rast(file.path(in_dir, "lad_profiles_z1_res_10_m.tif"))
  
  r_metrics <- terra::aggregate(
    c(r_lai, r_vci, r_hmax, r_fcover),
    fact = agg_factor, fun = "mean", na.rm = TRUE
  )
  names(r_metrics) <- c("LAI", "VCI", "Hmax", "fCover")
  
  r_lad_agg <- terra::aggregate(r_lad_st, fact = agg_factor, fun = "mean", na.rm = TRUE)
  
  list(metrics = r_metrics, lad = r_lad_agg, stack = c(r_metrics, r_lad_agg))
}

#' Extract daily macroclimate drivers from ERA5 forcing.
extract_macro_daily <- function(forcing_file, date_seq) {
  nc <- nc_open(forcing_file)
  
  get_var <- function(v) {
    val <- try(ncvar_get(nc, v), silent = TRUE)
    if (inherits(val, "try-error")) return(rep(NA, nc$dim$time$len))
    return(val)
  }
  
  macro_tair   <- get_var("Tair")
  macro_qair   <- get_var("Qair")
  macro_psurf  <- get_var("PSurf")
  macro_wind_n <- get_var("Wind_N")
  macro_wind_e <- get_var("Wind_E")
  macro_swdown <- get_var("SWdown")
  macro_rain   <- get_var("Rainf")
  
  nc_close(nc)
  macro_time <- force_utc_nc(forcing_file, "time")
  
  t_celsius  <- macro_tair - 273.15
  wind_abs   <- sqrt(macro_wind_n^2 + macro_wind_e^2)
  
  e_s <- 611.2 * exp((17.67 * t_celsius) / (t_celsius + 243.5))
  e_a <- (macro_qair * macro_psurf) / (0.622 + 0.378 * macro_qair)
  vpd <- pmax(e_s - e_a, 0)
  
  df_hourly <- data.frame(
    time       = macro_time, 
    Tair_macro = t_celsius,
    Wind_macro = wind_abs,
    Rad_macro  = macro_swdown,
    VPD_macro  = vpd,
    Rain_macro = macro_rain
  )
  
  df_hourly %>%
    filter(as.Date(time) %in% date_seq) %>%
    mutate(date = as.Date(time)) %>%
    group_by(date) %>%
    summarise(
      Tmax_macro = max(Tair_macro, na.rm = TRUE),
      Wind_mean  = mean(Wind_macro, na.rm = TRUE),
      Rad_mean   = mean(Rad_macro, na.rm = TRUE),
      VPD_mean   = mean(VPD_macro, na.rm = TRUE),
      Rain_sum   = sum(Rain_macro * 3600, na.rm = TRUE),
      .groups    = "drop"
    )
}

#' Save a ggplot to disk and return the path invisibly.
save_plot <- function(p, path, width = 10, height = 7, dpi = 150) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, p, width = width, height = height, dpi = dpi)
  invisible(path)
}

#' Create the standard outputs/ sub-directory tree.
setup_outputs_dirs <- function(root = "outputs") {
  for (sd in c("clusters", "h2", "h1", "fpca", "gamm", "hobo", "s2_annex", "audit")) {
    dir.create(file.path(root, sd), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(root)
}

#' Read HOBO sub-canopy temperatures.
read_hobo_daily <- function(csv_file, date_seq, df_macro, ids_to_remove) {
  read.csv(csv_file) %>%
    mutate(
      datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      date     = as.Date(datetime)
    ) %>%
    filter(position_sensor == "a", date %in% date_seq,
           !(id_plot %in% ids_to_remove)) %>%
    group_by(id_plot, date) %>%
    summarise(Tmax_obs = max(t_hobo, na.rm = TRUE), .groups = "drop") %>%
    inner_join(df_macro, by = "date") %>%
    mutate(Delta_obs = Tmax_obs - Tmax_macro)
}

# ==============================================================================
# SECTION 3.  FOREST DATAFRAME CONSTRUCTION
# ==============================================================================

build_forest_dataframe <- function(r_stack) {
  df_all <- as.data.frame(r_stack, xy = TRUE, na.rm = FALSE) %>%
    filter(if_any(-c(x, y), ~ !is.na(.)))
  
  mat_lad  <- df_all %>% dplyr::select(starts_with("LAD_Layer_")) %>% as.matrix()
  row_sums <- rowSums(mat_lad, na.rm = TRUE)
  valid    <- row_sums > 0.5
  
  df_forest <- df_all[valid, ]
  mat_lad   <- mat_lad[valid, ]
  mat_lad[is.na(mat_lad)] <- 0
  
  list(df = df_forest, mat_lad = mat_lad, 
       z_breaks = as.numeric(gsub("LAD_Layer_", "", colnames(mat_lad))))
}

optimise_k_elbow <- function(df_scaled, k_range = 1:10) {
  wss <- vapply(k_range, function(k) {
    kmeans(df_scaled, centers = k, nstart = 25, iter.max = 100)$tot.withinss
  }, numeric(1))

  p1 <- c(k_range[1], wss[1]); p2 <- c(k_range[length(k_range)], wss[length(wss)])
  d <- vapply(seq_along(k_range), function(i) {
    p0 <- c(k_range[i], wss[i])
    abs((p2[2]-p1[2])*p0[1] - (p2[1]-p1[1])*p0[2] + p2[1]*p1[2] - p2[2]*p1[1]) /
      sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
  }, numeric(1))

  opt_k <- k_range[which.max(d)]

  p <- ggplot(data.frame(k = k_range, WSS = wss), aes(x = k, y = WSS)) +
    geom_line(colour = "grey50", linewidth = 1) +
    geom_point(size = 3, colour = "#31688e") +
    geom_vline(xintercept = opt_k, linetype = "dashed", colour = "#d8576b", linewidth = 1) +
    scale_x_continuous(breaks = k_range) +
    labs(title = "Elbow Method for K-means Stratification",
         subtitle = sprintf("Optimal k = %d", opt_k),
         x = "Clusters (k)", y = "Total Within Sum of Squares")

  list(opt_k = opt_k, wss = wss, k_range = k_range, plot = p)
}

label_clusters <- function(df_forest, k = NULL, vars = c("LAI", "Hmax", "fCover"), max_k = 10) {
  df_in <- df_forest %>% dplyr::select(all_of(vars))
  ok    <- complete.cases(df_in)
  df_scaled <- scale(df_in[ok, ])

  elbow_plot <- NULL
  if (is.null(k)) {
    cat("  [Auto-K] Determining optimal number of clusters via elbow method...\n")
    elbow_res  <- optimise_k_elbow(df_scaled, k_range = 1:max_k)
    k          <- elbow_res$opt_k
    elbow_plot <- elbow_res$plot
    cat(sprintf("  [Auto-K] Optimal k selected: %d\n", k))
  }

  km <- kmeans(df_scaled, centers = k, iter.max = 100, nstart = 25)
  df_forest$Cluster <- NA_integer_
  df_forest$Cluster[ok] <- km$cluster
  df_forest$Cluster <- factor(df_forest$Cluster)
  list(df = df_forest, elbow_plot = elbow_plot)
}

#' Profil LAD moyen par cluster — base d'interprétation typologique.
#'
#' Permet de nommer visuellement les clusters (bottom-heavy, top-heavy, etc.)
#' en affichant le profil moyen de chaque groupe dans l'espace de hauteur absolu.
#'
#' @param df_forest Dataframe avec colonnes x, y, Cluster (issu de label_clusters).
#' @param mat_lad   Matrice LAD alignée ligne-à-ligne avec df_forest.
#' @param z_breaks  Vecteur de hauteurs (m) correspondant aux colonnes de mat_lad.
plot_cluster_mean_profiles <- function(df_forest, mat_lad, z_breaks) {
  valid <- !is.na(df_forest$Cluster)
  clusters <- levels(df_forest$Cluster[valid])

  df_profiles <- lapply(clusters, function(cl) {
    idx <- which(df_forest$Cluster == cl & !is.na(df_forest$Cluster))
    mean_prof <- colMeans(mat_lad[idx, , drop = FALSE], na.rm = TRUE)
    data.frame(height = z_breaks, density = mean_prof, Cluster = cl, n = length(idx))
  })
  df_profiles <- bind_rows(df_profiles) %>%
    mutate(label = sprintf("Cluster %s\n(n = %d)", Cluster, n))

  ggplot(df_profiles, aes(x = density, y = height, colour = Cluster)) +
    geom_path(linewidth = 1.2) +
    facet_wrap(~ label, nrow = 1) +
    scale_colour_viridis_d(option = "turbo") +
    labs(
      title    = "Profil LAD moyen par cluster \u2014 base d\u2019interpr\u00e9tation typologique",
      subtitle = "K-means sur LAI, Hmax, fCover (forêt entière) \u2014 nommer chaque profil (bottom-heavy, top-heavy\u2026)",
      x        = "Densit\u00e9 foliaire LAD (m\u207b\u00b9)",
      y        = "Hauteur (m)"
    ) +
    theme(legend.position = "none")
}

# ==============================================================================
# SECTION 4.  cLHS SAMPLING
# ==============================================================================

sample_clhs_per_cluster <- function(df_forest, n_per_cluster = 100, vars = c("LAI", "Hmax", "fCover"), iter = 10000) {
  df_forest %>%
    filter(!is.na(Cluster)) %>%
    split(.$Cluster) %>%
    map_dfr(function(df_sub) {
      df_lhs   <- df_sub %>% dplyr::select(x, y, all_of(vars))
      n_target <- min(n_per_cluster, nrow(df_sub))
      lhs_res  <- clhs::clhs(df_lhs, size = n_target, iter = iter, simple = FALSE, progress = FALSE)
      df_sub[lhs_res$index_samples, ]
    })
}

# ==============================================================================
# SECTION 5.  LAD PROFILE CONSTRUCTORS
# ==============================================================================

make_lad_real <- function(plot_row) {
  lad_cols <- grep("LAD_Layer_", names(plot_row), value = TRUE)
  z        <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  d        <- as.numeric(plot_row[lad_cols]); d[is.na(d)] <- 0
  hmax     <- as.numeric(plot_row$Hmax)
  data.frame(height = z, density = d) %>% filter(height <= ceiling(hmax) + 1)
}

make_lad_uniform <- function(plot_row, canopy_base = 2) {
  hmax       <- as.numeric(plot_row$Hmax)
  lai        <- as.numeric(plot_row$LAI)
  max_layer  <- ceiling(hmax)
  if (max_layer <= canopy_base) return(data.frame(height = 1, density = 0))
  depth   <- max_layer - canopy_base
  density <- lai / depth
  heights <- 1:max_layer
  data.frame(height = heights, density = ifelse(heights >= canopy_base, density, 0))
}

make_lad_mean_factory <- function(df_sample) {
  lad_cols <- grep("LAD_Layer_", names(df_sample), value = TRUE)
  z        <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  mat      <- as.matrix(df_sample[, lad_cols]); mat[is.na(mat)] <- 0
  mean_prof<- colMeans(mat, na.rm = TRUE)
  mean_lai <- sum(mean_prof)
  
  function(plot_row) {
    hmax <- as.numeric(plot_row$Hmax)
    lai  <- as.numeric(plot_row$LAI)
    scale_fac <- if (mean_lai > 0) lai / mean_lai else 0
    data.frame(height = z, density = mean_prof * scale_fac) %>%
      filter(height <= ceiling(hmax) + 1)
  }
}

# ==============================================================================
# SECTION 6.  SCENARIO REGISTRY
# ==============================================================================

fn_real <- function(varname) { force(varname); function(plot_row) as.numeric(plot_row[[varname]]) }
fn_mean <- function(df_sample, varname) { m <- mean(df_sample[[varname]], na.rm = TRUE); function(plot_row) m }

scenario_reference <- function(df_sample) {
  list(name = "REF_all_real", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), 
       fcover_fn = fn_real("fCover"), lad_fn = make_lad_real)
}

scenarios_h2_uniform_vs_real <- function(df_sample) {
  list(
    REF_real_LAD = list(name = "H2_real_LAD", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_real),
    UNIFORM_LAD  = list(name = "H2_uniform_LAD", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_uniform)
  )
}

scenarios_h1_forward <- function(df_sample) {
  lai_mean <- fn_mean(df_sample, "LAI"); hmax_mean <- fn_mean(df_sample, "Hmax"); fcover_mean <- fn_mean(df_sample, "fCover"); lad_uniform <- make_lad_uniform
  list(
    Null_baseline   = list(name = "H1f_0_Null_baseline", lai_fn = lai_mean, hmax_fn = hmax_mean, fcover_fn = fcover_mean, lad_fn = lad_uniform),
    LAI_only        = list(name = "H1f_1_LAI_only", lai_fn = fn_real("LAI"), hmax_fn = hmax_mean, fcover_fn = fcover_mean, lad_fn = lad_uniform),
    LAI_Hmax        = list(name = "H1f_2_LAI_Hmax", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fcover_mean, lad_fn = lad_uniform),
    LAI_Hmax_fCover = list(name = "H1f_3_LAI_Hmax_fCover", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = lad_uniform),
    Full_real       = list(name = "H1f_4_Full_real", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_real)
  )
}

scenarios_h1_loo <- function(df_sample) {
  hmax_mean <- fn_mean(df_sample, "Hmax"); lai_mean <- fn_mean(df_sample, "LAI"); fcover_mean <- fn_mean(df_sample, "fCover"); lad_mean_fn <- make_lad_mean_factory(df_sample)
  list(
    drop_LAI_mean      = list(name = "H1l_dropLAI_mean", lai_fn = lai_mean, hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_real),
    drop_Hmax_mean     = list(name = "H1l_dropHmax_mean", lai_fn = fn_real("LAI"), hmax_fn = hmax_mean, fcover_fn = fn_real("fCover"), lad_fn = make_lad_real),
    drop_fCover_mean   = list(name = "H1l_dropfCover_mean", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fcover_mean, lad_fn = make_lad_real),
    drop_LAD_uniform   = list(name = "H1l_dropLAD_uniform", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_uniform),
    drop_LAD_meanShape = list(name = "H1l_dropLAD_meanShape", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = lad_mean_fn)
  )
}

scenarios_h1_factorial <- function(df_sample) {
  lai_mean <- fn_mean(df_sample, "LAI"); hmax_mean <- fn_mean(df_sample, "Hmax"); fcover_mean <- fn_mean(df_sample, "fCover"); lad_mean_fn <- make_lad_mean_factory(df_sample)
  grid <- expand.grid(LAI = c("real", "mean"), Hmax = c("real", "mean"), fCover = c("real", "mean"), LAD = c("real", "meanShape"), stringsAsFactors = FALSE)
  scenarios <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    nm <- sprintf("H1F_%s_%s_%s_%s", substr(g$LAI, 1, 1), substr(g$Hmax, 1, 1), substr(g$fCover, 1, 1), substr(g$LAD, 1, 1))
    scenarios[[i]] <- list(name = nm, lai_fn = if (g$LAI == "real") fn_real("LAI") else lai_mean,
                           hmax_fn = if (g$Hmax == "real") fn_real("Hmax") else hmax_mean,
                           fcover_fn = if (g$fCover == "real") fn_real("fCover") else fcover_mean,
                           lad_fn = if (g$LAD == "real") make_lad_real else lad_mean_fn, grid_row = g)
  }
  setNames(scenarios, vapply(scenarios, `[[`, "", "name"))
}

# ==============================================================================
# SECTION 7.  MuSICA BATCH RUNNER
# ==============================================================================

run_musica_one <- function(plot_row, scenario, out_nc_file, forcing_file, musica_cmd) {
  if (file.exists(out_nc_file)) return(invisible(out_nc_file))
  lai <- scenario$lai_fn(plot_row); hmax <- scenario$hmax_fn(plot_row)
  fcover <- scenario$fcover_fn(plot_row); allom <- scenario$lad_fn(plot_row)
  if (is.na(lai) || is.na(hmax) || is.na(fcover)) return(invisible(NULL))
  
  phenology <- calc_phenology(list.year = c(2021, 2022), nleafage = 1, budburst_date = 115, leaf_age_max_in = 0.56, relative_age_firstmax = 0.10, relative_age_lastmax = 0.75, LAI_max_per_cohort = lai)
  tryCatch({
    callmusica(
      musica.param = list(setupctl = list("clumping_factor" = fcover, "forcing_filename" = forcing_file, "history_filename" = out_nc_file, "forcing_height" = hmax + 2)),
      leaf.param = list(musica_veg1 = list("phenology" = phenology, "allometry" = allom, "leafphenologyctl" = list("lai_max_per_cohort" = lai),
                                           "leafallometryctl" = list("canopy_height_top" = hmax, "canopy_height_bottom" = 1), "leafmusicactl" = list("canopy_height_top" = hmax))),
      musica.cmd = musica_cmd, keep.tmp = FALSE, out.netcdf = TRUE
    )
  }, error = function(e) cat(sprintf("ERROR on %s: %s\n", basename(out_nc_file), e$message)))
  invisible(out_nc_file)
}

run_musica_scenario <- function(df_sample, scenario, out_dir, forcing_file = CFG$forcing_file, musica_cmd = CFG$musica_cmd, id_prefix = "Sim", force = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, pattern = "\\.nc$")
  if (!force && length(existing) >= nrow(df_sample)) return(invisible(out_dir))
  cat(sprintf("\n[scenario: %s] %d plots → %s (%d done)\n", scenario$name, nrow(df_sample), out_dir, length(existing)))
  for (i in seq_len(nrow(df_sample))) {
    plot_row <- df_sample[i, ]; sim_id <- sprintf("%s_%03d_X%d_Y%d", id_prefix, i, round(plot_row$x), round(plot_row$y))
    out_file <- file.path(out_dir, paste0("musica_out_", sim_id, ".nc"))
    if (file.exists(out_file)) next
    run_musica_one(plot_row, scenario, out_file, forcing_file, musica_cmd)
  }
  invisible(out_dir)
}

run_scenarios <- function(df_sample, scenarios, parent_dir, force = FALSE) {
  paths <- character(length(scenarios)); names(paths) <- names(scenarios)
  for (nm in names(scenarios)) {
    out_d <- file.path(parent_dir, scenarios[[nm]]$name)
    run_musica_scenario(df_sample, scenarios[[nm]], out_d, force = force)
    paths[nm] <- out_d
  }
  paths
}

# ==============================================================================
# SECTION 8.  POST-PROCESSING
# ==============================================================================

extract_deltatmax_one <- function(nc_file, df_macro, date_seq) {
  nc <- try(nc_open(nc_file), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  raw <- try(get_variable(nc, "Tair_z"), silent = TRUE); nc_close(nc); if (inherits(raw, "try-error") || is.null(raw)) return(NULL)
  raw %>% filter(nair == 1) %>% mutate(Tair_sim = Tair_z - 273.15, time = time - hours(2), date = as.Date(time)) %>%
    filter(date %in% date_seq) %>% group_by(date) %>% summarise(Tmax_micro = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
    inner_join(df_macro, by = "date") %>% mutate(Delta_Tmax = Tmax_micro - Tmax_macro)
}

extract_deltatmax_scenario <- function(scenario_dir, df_macro, date_seq) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE); if (length(nc_files) == 0) return(NULL)
  map_df(nc_files, function(f) {
    res <- extract_deltatmax_one(f, df_macro, date_seq); if (is.null(res)) return(NULL)
    res %>% mutate(x = as.numeric(str_extract(basename(f), "(?<=_X)\\d+")), y = as.numeric(str_extract(basename(f), "(?<=_Y)\\d+")), scenario = basename(scenario_dir))
  })
}

extract_all_scenarios <- function(scenario_paths, df_macro, date_seq) {
  map_df(scenario_paths, extract_deltatmax_scenario, df_macro = df_macro, date_seq = date_seq)
}

join_with_sample <- function(df_daily, df_sample) {
  df_meta <- df_sample %>% dplyr::select(x, y, Cluster, LAI, Hmax, fCover, VCI, starts_with("LAD_Layer_"))
  df_daily %>% inner_join(df_meta, by = c("x", "y"))
}

# ==============================================================================
# SECTION 9.  H2 — UNIFORM vs REAL LAD
# ==============================================================================

summarise_h2 <- function(df_h2) {
  df_h2 %>% pivot_wider(id_cols = c(x, y, date), names_from = scenario, values_from = Delta_Tmax) %>%
    rename(Real = matches("real_LAD"), Uniform = matches("uniform_LAD")) %>% mutate(diff = Real - Uniform)
}

plot_h2_distributions <- function(df_h2) {
  ggplot(df_h2, aes(x = Delta_Tmax, fill = scenario, colour = scenario)) + geom_density(alpha = 0.3, linewidth = 0.8) +
    scale_fill_manual(values  = c("#31688e", "#d8576b")) + scale_colour_manual(values = c("#31688e", "#d8576b")) +
    labs(title = "H2 — ΔTmax distributions: real vs uniform LAD", x = expression(Delta * T[max] ~ (degree*C)), y = "Density")
}

plot_h2_paired <- function(df_h2_wide) {
  r2_val   <- cor(df_h2_wide$Real, df_h2_wide$Uniform, use = "complete.obs")^2
  rmse_val <- sqrt(mean((df_h2_wide$Uniform - df_h2_wide$Real)^2, na.rm = TRUE))
  bias_val <- mean(df_h2_wide$Uniform - df_h2_wide$Real, na.rm = TRUE)
  stats_label <- sprintf("R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", r2_val, rmse_val, bias_val)
  
  axis_lims <- c(min(c(df_h2_wide$Real, df_h2_wide$Uniform), na.rm = TRUE), max(c(df_h2_wide$Real, df_h2_wide$Uniform), na.rm = TRUE))
  
  ggplot(df_h2_wide, aes(x = Real, y = Uniform)) + geom_hex(bins = 100) +
    scale_fill_gradient(low = "grey80", high = "midnightblue", trans = "log10", name = "Count\n(log)") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", colour = "#d73027", se = FALSE, linewidth = 1.2) + coord_fixed(xlim = axis_lims, ylim = axis_lims) +
    annotate("text", x = -Inf, y = Inf, label = stats_label, hjust = -0.1, vjust = 1.2, fontface = "bold", size = 5) +
    labs(title = "H2 — Per-plot, per-day ΔTmax: real vs uniform LAD", subtitle = "Paired comparison: impact of vertical profile SHAPE on cooling",
         x = expression("Real LAD" ~ Delta * T[max] ~ (degree*C)), y = expression("Uniform LAD" ~ Delta * T[max] ~ (degree*C))) +
    theme_bw(base_size = 14) + theme(legend.position = "right", plot.title = element_text(face = "bold"))
}

analyse_h2_distribution <- function(df_h2_wide, df_macro, threshold = 30) {
  df <- df_h2_wide %>% left_join(df_macro, by = "date") %>% mutate(period = ifelse(Tmax_macro >= threshold, "heatwave", "normal"))
  summ <- df %>% group_by(period) %>% summarise(n = n(), mean_diff = mean(diff, na.rm = TRUE), sd_diff = sd(diff, na.rm = TRUE), .groups = "drop")
  print(summ)
  p <- ggplot(df, aes(x = diff, fill = period)) + geom_histogram(bins = 80, alpha = 0.6, position = "identity") + geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c(normal = "#31688e", heatwave = "#d8576b")) + labs(title = "H2 — distribution of (Real - Uniform) ΔTmax", subtitle = "Split by macroclimate regime", x = "ΔTmax difference (°C)", y = "Count")
  print(p)
  invisible(p)
}

validate_hobo_paired <- function(hobo_res) {
  df <- hobo_res$daily
  key_counts <- df %>% count(id_plot, date, name = "n_sc")
  common_keys <- key_counts %>% filter(n_sc == length(unique(df$scenario)))
  df %>% inner_join(common_keys %>% select(id_plot, date), by = c("id_plot", "date")) %>% group_by(scenario) %>%
    summarise(n = n(), r2 = cor(Delta_obs, Delta_sim, use = "complete.obs")^2, rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)),
              mae  = mean(abs(Delta_obs - Delta_sim), na.rm = TRUE), bias = mean(Delta_sim - Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(rmse)
}

# ==============================================================================
# SECTION 10.  H1 — Forward / LOO
# ==============================================================================

score_scenarios_vs_reference <- function(df_all, ref_scenario_name) {
  df_ref <- df_all %>% filter(scenario == ref_scenario_name) %>% dplyr::select(x, y, date, Delta_ref = Delta_Tmax)
  df_all %>% inner_join(df_ref, by = c("x", "y", "date")) %>% mutate(diff = Delta_Tmax - Delta_ref) %>% group_by(scenario) %>%
    summarise(n = n(), mean_diff = mean(diff, na.rm = TRUE), mae = mean(abs(diff), na.rm = TRUE), rmse = sqrt(mean(diff^2, na.rm = TRUE)), .groups = "drop") %>% arrange(rmse)
}

plot_scenario_hierarchy <- function(df_scores, ref_scenario_name) {
  df_plot <- df_scores %>% filter(scenario != ref_scenario_name)
  ggplot(df_plot, aes(x = reorder(scenario, rmse), y = rmse, fill = mean_diff)) + geom_col() + coord_flip() +
    scale_fill_gradient2(low = "#31688e", mid = "white", high = "#d8576b", midpoint = 0, name = "Mean bias\n(°C)") +
    labs(title = "H1 — divergence from full-real reference", x = NULL, y = "RMSE vs reference (°C)")
}

plot_forward_curve <- function(df_scores, forward_order) {
  df_plot <- df_scores %>% filter(scenario %in% forward_order) %>% mutate(scenario = factor(scenario, levels = forward_order), step = as.integer(scenario))
  ggplot(df_plot, aes(x = step, y = rmse)) + geom_line(linewidth = 1) + geom_point(size = 3, colour = "#31688e") +
    scale_x_continuous(breaks = seq_along(forward_order), labels = forward_order) + labs(title = "H1 forward — information added incrementally", x = NULL, y = "RMSE vs full-real reference (°C)") + theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# ==============================================================================
# SECTION 11.  SHAPE METRICS (FPCA & Height of Median LAI)
# ==============================================================================

optimise_nbasis <- function(arg_vals, mat_lad_norm, nbasis_range = 4:15) {
  gcv <- vapply(nbasis_range, function(nb) {
    bb <- create.bspline.basis(c(min(arg_vals), max(arg_vals)), nbasis = nb)
    mean(smooth.basis(arg_vals, t(mat_lad_norm), bb)$gcv)
  }, numeric(1))

  p1 <- c(nbasis_range[1], gcv[1]); p2 <- c(nbasis_range[length(nbasis_range)], gcv[length(gcv)])
  d <- vapply(seq_along(nbasis_range), function(i) {
    p0 <- c(nbasis_range[i], gcv[i])
    abs((p2[2]-p1[2])*p0[1] - (p2[1]-p1[1])*p0[2] + p2[1]*p1[2] - p2[2]*p1[1]) /
      sqrt((p2[2]-p1[2])^2 + (p2[1]-p1[1])^2)
  }, numeric(1))
  opt_nb <- nbasis_range[which.max(d)]

  p <- ggplot(data.frame(nbasis = nbasis_range, GCV = gcv), aes(x = nbasis, y = GCV)) +
    geom_line(colour = "grey50", linewidth = 1) +
    geom_point(size = 3, colour = "#31688e") +
    geom_vline(xintercept = opt_nb, linetype = "dashed", colour = "#d8576b", linewidth = 1) +
    labs(title = "GCV Curve for FPCA Basis Functions",
         x = "Number of B-spline basis functions", y = "Mean GCV") +
    theme_bw()

  list(opt = opt_nb, gcv = gcv, range = nbasis_range, plot = p)
}

compute_fpca <- function(mat_lad, hmax_vec, z_breaks, n_harm = 3) {
  row_sums <- rowSums(mat_lad, na.rm = TRUE)
  mat_norm <- sweep(mat_lad, 1, pmax(row_sums, 1e-9), "/")
  n_bins   <- length(z_breaks); z_rel <- seq(0, 1, length.out = n_bins)
  mat_rel  <- matrix(0, nrow = nrow(mat_norm), ncol = n_bins)
  for (i in seq_len(nrow(mat_norm))) {
    h_i      <- max(hmax_vec[i], 1); z_i_rel <- z_breaks / h_i
    mat_rel[i, ] <- pmax(approx(z_i_rel, mat_norm[i, ], z_rel, rule = 2)$y, 0)
  }
  rs <- rowSums(mat_rel); rs[rs == 0] <- 1; mat_rel <- sweep(mat_rel, 1, rs, "/")
  opt    <- optimise_nbasis(z_rel, mat_rel)
  basis  <- create.bspline.basis(c(0, 1), nbasis = opt$opt)
  fd_obj <- Data2fd(z_rel, t(mat_rel), basis)
  fpca   <- pca.fd(fd_obj, nharm = n_harm)
  list(fpca = fpca, fd_obj = fd_obj, mat_rel = mat_rel, z_rel = z_rel,
       nbasis_opt = opt$opt, varprop = fpca$varprop, gcv_plot = opt$plot)
}

plot_fpc_harmonics <- function(fpca_res, mat_lad_real, z_breaks, hmax_vec) {
  scores <- fpca_res$fpca$scores; z_rel <- fpca_res$z_rel; mat_r <- fpca_res$mat_rel
  one <- function(k, label) {
    q <- quantile(scores[, k], c(0.05, 0.95), na.rm = TRUE)
    lo <- which(scores[, k] <= q[1]); hi <- which(scores[, k] >= q[2])
    df <- bind_rows(data.frame(rel_z = z_rel, lad = colMeans(mat_r[lo, , drop = FALSE]), grp = "Low 5%"), data.frame(rel_z = z_rel, lad = colMeans(mat_r), grp = "Mean"), data.frame(rel_z = z_rel, lad = colMeans(mat_r[hi, , drop = FALSE]), grp = "High 5%"))
    ggplot(df, aes(x = lad, y = rel_z, colour = grp)) + geom_path(linewidth = 1.1) + scale_colour_manual(values = c("#fca50a", "grey50", "#31688e")) + labs(title = sprintf("%s — %.1f%% var", label, 100 * fpca_res$varprop[k]), x = "Relative LAD", y = "Z / Hmax", colour = NULL) + theme(legend.position = "bottom")
  }
  one(1, "FPC1") | one(2, "FPC2") | one(3, "FPC3")
}

plot_fpc_loadings <- function(fpca_res) {
  harmonics <- fpca_res$fpca$harmonics; z_rel <- fpca_res$z_rel; vals <- eval.fd(z_rel, harmonics)
  df <- as.data.frame(vals); names(df) <- paste0("FPC", seq_len(ncol(df))); df$rel_z <- z_rel
  df_long <- pivot_longer(df, -rel_z, names_to = "FPC", values_to = "loading")
  ggplot(df_long, aes(x = loading, y = rel_z, colour = FPC)) + geom_path(linewidth = 1.1) + geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~ FPC, nrow = 1) + scale_colour_viridis_d(option = "turbo", end = 0.85) + labs(title = "Per-height contribution of each FPC (FPCA loadings)", x = "Loading", y = "Relative height (Z/Hmax)") + theme(legend.position = "none")
}

plot_fpc_scatters <- function(fpca_res, cluster_vec = NULL) {
  scores <- as.data.frame(fpca_res$fpca$scores[, 1:3])
  names(scores) <- c("FPC1", "FPC2", "FPC3")
  if (!is.null(cluster_vec)) scores$Cluster <- factor(cluster_vec)

  base <- function(xv, yv) {
    p <- ggplot(scores, aes(x = .data[[xv]], y = .data[[yv]]))
    if (!is.null(cluster_vec)) {
      p <- p + geom_point(aes(colour = Cluster), alpha = 0.5) + scale_colour_viridis_d(option = "turbo")
    } else {
      p <- p + geom_point(alpha = 0.4, colour = "#31688e")
    }
    p + labs(x = xv, y = yv)
  }
  (base("FPC1", "FPC2") | base("FPC1", "FPC3") | base("FPC2", "FPC3")) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
}

#' Compute the Height of Median LAI (Center of Gravity) for each plot.
compute_h_median <- function(mat_lad, z_breaks) {
  apply(mat_lad, 1, function(prof) {
    total_lai <- sum(prof, na.rm = TRUE)
    if (total_lai <= 0) return(0)
    cum_lai <- cumsum(prof)
    idx <- which(cum_lai >= (total_lai / 2))[1]
    return(z_breaks[idx])
  })
}

# ==============================================================================
# SECTION 12.  GAMM EMULATOR
# ==============================================================================

#' Prepare scaled inputs for the GAMM, explicitly integrating macroclimate and shape drivers.
prepare_gamm_data <- function(df_daily, df_sample, fpc_scores = NULL, h_median_vec = NULL) {
  df <- df_daily %>%
    mutate(plot_id = as.factor(paste0("X", x, "_Y", y)),
           date_factor = as.factor(date))
  
  if (!is.null(h_median_vec)) {
    df_sample$H_median <- h_median_vec
    df <- df %>% left_join(df_sample %>% dplyr::select(x, y, H_median), by = c("x", "y"))
  }
  
  if ("Cluster" %in% names(df)) df$Cluster <- as.factor(df$Cluster)

  scale_cols <- c("LAI", "Hmax", "fCover", "H_median", "Tmax_macro",
                  "Wind_mean", "Rad_mean", "VPD_mean", "Rain_sum")
  cols_to_scale <- intersect(scale_cols, names(df))
  for (nm in cols_to_scale) df[[paste0(nm, "_sc")]] <- as.numeric(scale(df[[nm]]))

  if (!is.null(fpc_scores)) {
    fpc_df <- as.data.frame(fpc_scores)
    names(fpc_df) <- paste0("FPC", seq_len(ncol(fpc_df)), "_sc")
    fpc_df[] <- lapply(fpc_df, function(v) as.numeric(scale(v)))
    fpc_df$x <- df_sample$x; fpc_df$y <- df_sample$y
    df <- df %>% left_join(fpc_df, by = c("x", "y"))
  }
  return(df)
}

#' Fit the reference GAMM, with toggles for macroclimate and shape predictor.
fit_reference_gamm <- function(df_gamm, shape_type = "NONE", use_macroclimate = TRUE) {
  base_terms <- c("s(LAI_sc)", "s(Hmax_sc)", "s(fCover_sc)")
  
  clim_terms <- if (use_macroclimate) {
    cat("  [GAMM] Including macroclimate drivers (Tmax, Wind, Rad, VPD, Rain)...\n")
    c("s(Tmax_macro_sc)", "s(Wind_mean_sc)", "s(Rad_mean_sc)", "s(VPD_mean_sc)", "s(Rain_sum_sc)")
  } else {
    cat("  [GAMM] Structural drivers only (Macroclimate excluded)...\n")
    character(0)
  }
  
  shape_terms <- switch(shape_type,
                        "CLUSTER" = {
                          cat("  [GAMM] Including profile Cluster type as categorical predictor...\n")
                          c("Cluster")   # terme paramétrique factoriel (directement interprétable)
                        },
                        "FPC" = {
                          cat("  [GAMM] Including FPC scores as shape predictors...\n")
                          c("s(FPC1_sc)", "s(FPC2_sc)", "s(FPC3_sc)")
                        },
                        "H_MEDIAN" = {
                          cat("  [GAMM] Including structural shape (H_median) as predictor...\n")
                          c("s(H_median_sc)")
                        },
                        {
                          cat("  [GAMM] Structural macro-metrics only (no specific shape variable)...\n")
                          character(0)
                        }
  )
  
  re_terms <- c("s(date_factor, bs='re')", "s(plot_id, bs='re')")
  rhs <- paste(c(base_terms, clim_terms, shape_terms, re_terms), collapse = " + ")
  form <- as.formula(paste("Delta_Tmax ~", rhs))
  
  bam(form, data = df_gamm, family = scat(), method = "fREML", discrete = TRUE)
}

#' Plot marginal effects of the GAMM predictors.
plot_gamm_marginal_effects <- function(gam_model, shape_type = "NONE", use_macroclimate = TRUE) {
  terms_to_plot <- c("LAI_sc", "Hmax_sc")

  if (shape_type == "CLUSTER") terms_to_plot <- c(terms_to_plot, "Cluster")
  if (shape_type == "FPC")     terms_to_plot <- c(terms_to_plot, "FPC1_sc")
  if (shape_type == "H_MEDIAN") terms_to_plot <- c(terms_to_plot, "H_median_sc")
  if (use_macroclimate) terms_to_plot <- c(terms_to_plot, "VPD_mean_sc", "Wind_mean_sc", "Rad_mean_sc")

  plots <- lapply(terms_to_plot, function(term) {
    pred <- ggpredict(gam_model, terms = term)
    plot(pred) +
      labs(title = paste("Effect of", gsub("_sc", "", term)),
           x = if (term == "Cluster") "Profile cluster type" else paste(gsub("_sc", "", term), "(Scaled)"),
           y = "Predicted \u0394Tmax (\u00b0C)") +
      theme_bw() + theme(plot.title = element_text(size = 10))
  })

  wrap_plots(plots, ncol = 2) + plot_annotation(title = "GAMM Marginal Effects", subtitle = "Holding all other variables at their mean")
}

diagnose_residuals <- function(gam_scat, df_gamm, plot = TRUE) {
  gam_gauss <- bam(formula(gam_scat), data = df_gamm, family = gaussian(), method = "fREML", discrete = TRUE)
  resid_g <- residuals(gam_gauss, type = "deviance")
  idx_sw  <- sample(length(resid_g), min(5000, length(resid_g))); sw_test <- shapiro.test(resid_g[idx_sw])
  excess_kurt <- mean((resid_g - mean(resid_g))^4) / sd(resid_g)^4 - 3
  aic_g <- AIC(gam_gauss); aic_s <- AIC(gam_scat)
  
  if (plot) {
    op <- par(mfrow = c(2, 2)); on.exit(par(op))
    qqnorm(resid_g, main = "QQ-plot (Gaussian residuals)", pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3)); qqline(resid_g, col = "red", lwd = 2)
    hist(resid_g, breaks = 80, freq = FALSE, main = "Residual distribution", xlab = "Residual", col = "lightblue", border = "white"); curve(dnorm(x, mean = mean(resid_g), sd = sd(resid_g)), add = TRUE, col = "red", lwd = 2)
    plot(fitted(gam_gauss), resid_g, pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3), xlab = "Fitted", ylab = "Residual", main = "Residuals vs fitted"); abline(h = 0, col = "red", lwd = 2); lines(lowess(fitted(gam_gauss), resid_g), col = "blue", lwd = 2)
    barplot(c(Gaussian = aic_g, `scat()` = aic_s), col = c("lightblue", "salmon"), main = "AIC comparison", ylab = "AIC")
  }
  
  data.frame(n_obs = length(resid_g), shapiro_W = unname(sw_test$statistic), shapiro_p = sw_test$p.value, excess_kurt = excess_kurt, aic_gauss = aic_g, aic_scat = aic_s, delta_aic = aic_g - aic_s, scat_preferred = (aic_g - aic_s) > 10)
}

# ==============================================================================
# SECTION 13.  DIAGNOSTICS — per-HOBO time-series overlay
# ==============================================================================

pick_representative_hobos <- function(df_hobo_daily) {
  df_summary <- df_hobo_daily %>% group_by(id_plot) %>% summarise(mean_delta = mean(Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(mean_delta) 
  n_total <- nrow(df_summary)
  res <- bind_rows(
    df_summary %>% slice(1) %>% mutate(rationale = "Strongest cooling (Max buffer)"),
    df_summary %>% slice(ceiling(n_total / 2)) %>% mutate(rationale = "Median cooling"),
    df_summary %>% slice(n_total) %>% mutate(rationale = "Weakest cooling (Min buffer)")
  ) %>% mutate(facet_label = sprintf("HOBO %s\n[%s : Mean ΔTmax = %.1f°C]", id_plot, rationale, mean_delta))
  return(res)
}

plot_timeseries_faceted <- function(df_selected_hobos, df_hobo, df_musica_real, df_musica_uniform, df_macro) {
  ids <- df_selected_hobos$id_plot
  df1 <- df_hobo %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_obs, source = "HOBO")
  df2 <- df_musica_real %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_micro, source = "MuSICA real LAD")
  df3 <- df_musica_uniform %>% filter(id_plot %in% ids) %>% transmute(id_plot, date, value = Tmax_micro, source = "MuSICA uniform LAD")
  df4 <- tidyr::expand_grid(id_plot = ids, df_macro) %>% transmute(id_plot, date, value = Tmax_macro, source = "ERA5 macro")
  
  df_all <- bind_rows(df1, df2, df3, df4) %>% inner_join(df_selected_hobos %>% dplyr::select(id_plot, facet_label), by = "id_plot") %>% mutate(facet_label = factor(facet_label, levels = df_selected_hobos$facet_label))
  
  ggplot(df_all, aes(x = date, y = value, colour = source)) + geom_line(linewidth = 0.6, alpha = 0.85) + facet_wrap(~ facet_label, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("HOBO" = "black", "MuSICA real LAD" = "#31688e", "MuSICA uniform LAD" = "#d8576b", "ERA5 macro" = "grey60")) +
    labs(title = "Time-series — Selected HOBOs vs MuSICA", subtitle = "Plots chosen by their average buffering capacity during summer", x = NULL, y = expression(T[max] ~ (degree*C)), colour = NULL) + theme(strip.text = element_text(size = 11, lineheight = 1.2))
}

# ==============================================================================
# SECTION 14.  HOBO VALIDATION (per-scenario)
# ==============================================================================

build_hobo_inputs <- function(hobo_geojson, r_stack, ids_to_remove) {
  hobo_pts <- st_read(hobo_geojson, quiet = TRUE)
  ext_vals <- terra::extract(r_stack, vect(hobo_pts), xy = TRUE)
  hobo_pts %>% st_drop_geometry() %>% dplyr::select(id_plot) %>% bind_cols(ext_vals) %>% filter(!(id_plot %in% ids_to_remove)) %>% drop_na(LAI, Hmax, fCover)
}

extract_deltatmax_hobo_scenario <- function(scenario_dir, df_macro, date_seq) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE); if (length(nc_files) == 0) return(NULL)
  map_df(nc_files, function(f) {
    id_str <- str_extract(basename(f), "(?<=HOBO_).*(?=\\.nc)"); res <- extract_deltatmax_one(f, df_macro, date_seq); if (is.null(res)) return(NULL)
    res %>% mutate(id_plot = id_str, scenario = basename(scenario_dir))
  })
}

validate_scenarios_at_hobos <- function(df_hobo_inputs, df_hobo_daily, scenarios, parent_dir, df_macro, date_seq, force = FALSE) {
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  hobo_paths <- character(length(scenarios)); names(hobo_paths) <- names(scenarios)
  for (nm in names(scenarios)) {
    sc <- scenarios[[nm]]; out_d <- file.path(parent_dir, sc$name); dir.create(out_d, recursive = TRUE, showWarnings = FALSE)
    existing <- list.files(out_d, pattern = "\\.nc$")
    if (!force && length(existing) >= nrow(df_hobo_inputs)) { hobo_paths[nm] <- out_d; next }
    cat(sprintf("\n[HOBO %s] %d HOBO plots → %s\n", sc$name, nrow(df_hobo_inputs), out_d))
    for (i in seq_len(nrow(df_hobo_inputs))) {
      plot_row <- df_hobo_inputs[i, ]; sim_id <- sprintf("HOBO_%s", plot_row$id_plot)
      out_file <- file.path(out_d, paste0("musica_out_", sim_id, ".nc")); if (file.exists(out_file)) next
      run_musica_one(plot_row, sc, out_file, CFG$forcing_file, CFG$musica_cmd)
    }
    hobo_paths[nm] <- out_d
  }
  df_sims <- map_df(hobo_paths, extract_deltatmax_hobo_scenario, df_macro = df_macro, date_seq = date_seq)
  df_obs <- df_hobo_daily %>% dplyr::select(id_plot, date, Delta_obs)
  df_daily <- df_sims %>% inner_join(df_obs, by = c("id_plot", "date")) %>% rename(Delta_sim = Delta_Tmax)
  df_metrics <- df_daily %>% group_by(scenario) %>% summarise(n = n(), r2 = cor(Delta_obs, Delta_sim, use = "complete.obs")^2, rmse = sqrt(mean((Delta_obs - Delta_sim)^2, na.rm = TRUE)), mae = mean(abs(Delta_obs - Delta_sim), na.rm = TRUE), bias = mean(Delta_sim - Delta_obs, na.rm = TRUE), .groups = "drop") %>% arrange(rmse)
  list(daily = df_daily, metrics = df_metrics)
}

plot_hobo_validation <- function(df_metrics, title = "HOBO validation — RMSE per scenario") {
  ggplot(df_metrics, aes(x = reorder(scenario, rmse), y = rmse, fill = bias)) + geom_col() + geom_text(aes(label = sprintf("R²=%.2f", r2)), hjust = -0.1, size = 3.2) + coord_flip() + scale_fill_gradient2(low = "#31688e", mid = "white", high = "#d8576b", midpoint = 0, name = "Mean bias\n(°C)") + labs(title = title, x = NULL, y = "RMSE vs observed HOBO ΔTmax (°C)") + expand_limits(y = max(df_metrics$rmse) * 1.15)
}

# ==============================================================================
# SECTION 15.  Sentinel-2 + FORMS-H ANNEX (H3)
# ==============================================================================

make_s2_formsh_scenario <- function(df_sample, in_dir, agg_factor = 2) {
  r_s2     <- rast(file.path(in_dir, "s2lai_summer_atbd_res_10_m.tif")); r_formsh <- rast(file.path(in_dir, "FORMS-H_Blois.tif"))
  r_s2_20  <- terra::aggregate(r_s2, fact = agg_factor, fun = "mean", na.rm = TRUE); r_fh_20  <- resample(r_formsh, r_s2_20, method = "bilinear")
  pts <- vect(as.matrix(df_sample[, c("x", "y")]), crs = crs(r_s2_20))
  df_sample$S2_LAI  <- as.numeric(terra::extract(r_s2_20, pts)[, 2]); df_sample$FORMS_H <- as.numeric(terra::extract(r_fh_20,  pts)[, 2])
  na_count <- sum(is.na(df_sample$S2_LAI) | is.na(df_sample$FORMS_H))
  if (na_count > 0) {
    df_sample$S2_LAI[is.na(df_sample$S2_LAI)] <- mean(df_sample$S2_LAI, na.rm = TRUE)
    df_sample$FORMS_H[is.na(df_sample$FORMS_H)] <- mean(df_sample$FORMS_H, na.rm = TRUE)
  }
  scenario <- list(name = "ANNEX_S2_FORMSH", lai_fn = function(plot_row) as.numeric(plot_row$S2_LAI), hmax_fn = function(plot_row) as.numeric(plot_row$FORMS_H), fcover_fn = function(plot_row) 1.0, lad_fn = make_lad_uniform)
  list(scenario = scenario, df_sample = df_sample)
}

# ==============================================================================
# SECTION 16.  main() — orchestrator
# ==============================================================================

main <- function() {

  cat("[0/7] Setup output directories...\n")
  setup_outputs_dirs()

  cat("[1/7] Loading LiDAR rasters...\n")
  rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)

  cat("[2/7] Building forest dataframe...\n")
  fd      <- build_forest_dataframe(rasters$stack)
  cl_res  <- label_clusters(fd$df, k = CFG$k_clusters)
  df_forest <- cl_res$df
  if (!is.null(cl_res$elbow_plot)) {
    save_plot(cl_res$elbow_plot, "outputs/clusters/kmeans_elbow.png", width = 8, height = 5)
    print(cl_res$elbow_plot)
  }

  df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)

  # ---- 2b. FPCA sur la forêt ENTIÈRE (avant échantillonnage) -----------------
  cat("[2b] FPCA sur la forêt entière (avant cLHS)...\n")
  fpca_res <- compute_fpca(fd$mat_lad, fd$df$Hmax, fd$z_breaks)
  fpc_var  <- fpca_res$varprop
  cat(sprintf("    FPC variance (forêt entière) : %s\n",
              paste(round(100 * fpc_var, 1), collapse = " / ")))

  writeLines(sprintf("FPC1: %.1f%%\nFPC2: %.1f%%\nFPC3: %.1f%%",
                     100 * fpc_var[1], 100 * fpc_var[2], 100 * fpc_var[3]),
             "outputs/fpca/fpca_variance.txt")

  save_plot(fpca_res$gcv_plot, "outputs/fpca/fpca_gcv_curve.png", width = 8, height = 5)
  print(fpca_res$gcv_plot)

  p_fpca_harm <- plot_fpc_harmonics(fpca_res, fd$mat_lad, fd$z_breaks, fd$df$Hmax)
  save_plot(p_fpca_harm, "outputs/fpca/fpca_harmonics.png", width = 12, height = 6)
  print(p_fpca_harm)

  p_fpca_load <- plot_fpc_loadings(fpca_res)
  save_plot(p_fpca_load, "outputs/fpca/fpca_loadings.png", width = 10, height = 5)
  print(p_fpca_load)

  p_fpca_scat <- plot_fpc_scatters(fpca_res, cluster_vec = df_forest$Cluster)
  save_plot(p_fpca_scat, "outputs/fpca/fpca_scatters.png", width = 12, height = 5)
  print(p_fpca_scat)

  p_clust_prof <- plot_cluster_mean_profiles(df_forest, fd$mat_lad, fd$z_breaks)
  save_plot(p_clust_prof, "outputs/clusters/cluster_mean_profiles.png", width = 12, height = 6)
  print(p_clust_prof)

  # ---- 3. cLHS sampling -------------------------------------------------------
  if (FLAGS$RUN_CLHS) {
    cat("[3/7] cLHS sampling...\n")
    df_sample <- sample_clhs_per_cluster(df_forest, CFG$n_per_cluster)
    saveRDS(df_sample, file.path(CFG$out_dir, "clhs_sample.rds"))
  } else {
    df_sample <- readRDS(file.path(CFG$out_dir, "clhs_sample.rds"))
    # Compatibilité : RDS anciens contiennent "Archetype" au lieu de "Cluster"
    if ("Archetype" %in% names(df_sample) && !"Cluster" %in% names(df_sample)) {
      df_sample <- df_sample %>% rename(Cluster = Archetype)
    }
  }
  cat(sprintf("    sample size: %d plots\n", nrow(df_sample)))

  # Projection des scores FPCA (forêt entière) sur le sous-ensemble cLHS.
  # Les scores sont alignés ligne-à-ligne avec fd$df / df_forest ; on retrouve
  # les indices par correspondance (x, y).
  idx_sample       <- match(paste(df_sample$x, df_sample$y),
                            paste(df_forest$x, df_forest$y))
  fpc_scores_sample <- fpca_res$fpca$scores[idx_sample, 1:3]
  
  ref_sc <- scenario_reference(df_sample)
  
  # ---- 4. H2 — uniform vs real LAD -----------------------------------------
  if (FLAGS$RUN_H2_UNIFORM_VS_REAL) {
    cat("[4/7] H2: uniform vs real LAD...\n")
    h2_scs   <- scenarios_h2_uniform_vs_real(df_sample)
    h2_paths <- run_scenarios(df_sample, h2_scs, CFG$out_h2)
    df_h2    <- extract_all_scenarios(h2_paths, df_macro, CFG$date_seq)
    df_h2_w  <- summarise_h2(df_h2)
    
    p_h2_dist <- plot_h2_distributions(df_h2)
    save_plot(p_h2_dist, "outputs/h2/h2_distributions.png")
    print(p_h2_dist)

    p_h2_pair <- plot_h2_paired(df_h2_w)
    save_plot(p_h2_pair, "outputs/h2/h2_paired_real_vs_uniform.png")
    print(p_h2_pair)

    cat(sprintf("    mean per-plot diff (real - uniform): %.3f °C\n", mean(df_h2_w$diff, na.rm = TRUE)))
    p_h2_hw <- analyse_h2_distribution(df_h2_w, df_macro, CFG$heatwave_thr)
    save_plot(p_h2_hw, "outputs/h2/h2_heatwave_histogram.png")
  }
  
  # ---- 5. H1 — forward / LOO / factorial -----------------------------------
  scenario_paths <- list()
  ref_path <- file.path(CFG$out_h1_forward, ref_sc$name)
  run_musica_scenario(df_sample, ref_sc, ref_path)
  scenario_paths[[ref_sc$name]] <- ref_path
  
  if (FLAGS$RUN_H1_FORWARD) {
    cat("[5a/7] H1 forward inclusion...\n")
    fw_scs <- scenarios_h1_forward(df_sample)
    scenario_paths <- c(scenario_paths, run_scenarios(df_sample, fw_scs, CFG$out_h1_forward))
  }
  if (FLAGS$RUN_H1_LOO) {
    cat("[5b/7] H1 leave-one-out...\n")
    loo_scs <- scenarios_h1_loo(df_sample)
    scenario_paths <- c(scenario_paths, run_scenarios(df_sample, loo_scs, CFG$out_h1_loo))
  }
  if (FLAGS$RUN_H1_FACTORIAL) {
    cat("[5c/7] H1 full factorial (HEAVY)...\n")
    fac_scs <- scenarios_h1_factorial(df_sample)
    scenario_paths <- c(scenario_paths, run_scenarios(df_sample, fac_scs, CFG$out_h1_factorial))
  }
  
  df_all_scenarios <- extract_all_scenarios(scenario_paths, df_macro, CFG$date_seq)
  df_scores        <- score_scenarios_vs_reference(df_all_scenarios, ref_sc$name)
  print(df_scores)
  write.csv(df_scores, "outputs/h1/h1_scores.csv", row.names = FALSE)

  p_h1_hier <- plot_scenario_hierarchy(df_scores, ref_sc$name)
  save_plot(p_h1_hier, "outputs/h1/h1_hierarchy.png")
  print(p_h1_hier)

  if (FLAGS$RUN_H1_FORWARD) {
    forward_order <- vapply(scenarios_h1_forward(df_sample), `[[`, "", "name")
    p_h1_fwd <- plot_forward_curve(df_scores, forward_order)
    save_plot(p_h1_fwd, "outputs/h1/h1_forward_curve.png")
    print(p_h1_fwd)
  }
  
  # ---- 6. GAMM emulator -----------------------------------------------------
  if (FLAGS$RUN_GAMM_EMULATOR) {
    cat("[6/7] GAMM emulator on reference scenario...\n")
    df_ref_daily <- df_all_scenarios %>% filter(scenario == ref_sc$name) %>% join_with_sample(df_sample)

    # Résolution des prédicteurs de forme (scores FPCA déjà calculés sur la
    # forêt entière et projetés sur l'échantillon à l'étape 2b/3).
    fpc_scores   <- if (FLAGS$GAMM_SHAPE_VAR == "FPC") fpc_scores_sample else NULL
    h_median_vec <- if (FLAGS$GAMM_SHAPE_VAR == "H_MEDIAN") {
      mat_s  <- df_sample %>% dplyr::select(starts_with("LAD_Layer_")) %>% as.matrix()
      mat_s[is.na(mat_s)] <- 0
      compute_h_median(mat_s, as.numeric(gsub("LAD_Layer_", "", colnames(mat_s))))
    } else NULL
    # CLUSTER : déjà présent dans df_gamm via join_with_sample → rien à ajouter.

    df_gamm <- prepare_gamm_data(df_ref_daily, df_sample,
                                  fpc_scores = fpc_scores, h_median_vec = h_median_vec)

    gam_ref <- fit_reference_gamm(df_gamm,
                                   shape_type      = FLAGS$GAMM_SHAPE_VAR,
                                   use_macroclimate= FLAGS$USE_MACROCLIMATE_GAMM)
    capture.output(summary(gam_ref), file = "outputs/gamm/gamm_summary.txt")
    print(summary(gam_ref))

    p_gamm_eff <- plot_gamm_marginal_effects(gam_ref,
                                              shape_type      = FLAGS$GAMM_SHAPE_VAR,
                                              use_macroclimate= FLAGS$USE_MACROCLIMATE_GAMM)
    save_plot(p_gamm_eff, "outputs/gamm/gamm_marginal_effects.png", width = 12, height = 8)
    print(p_gamm_eff)

    png("outputs/gamm/gamm_residuals_diagnostic.png", width = 1200, height = 900, res = 120)
    diag_res <- diagnose_residuals(gam_ref, df_gamm)
    dev.off()
    print(diag_res)
    write.csv(diag_res, "outputs/gamm/gamm_residuals_stats.csv", row.names = FALSE)
  }
  
  # ---- OPTIONAL: HOBO validation --------------------------------------------
  if (FLAGS$RUN_HOBO_VALIDATION) {
    cat("[HOBO] Running validation of forward scenarios at HOBO locations...\n")
    df_hobo_daily  <- read_hobo_daily(CFG$hobo_temp_csv, CFG$date_seq, df_macro, CFG$ids_to_remove)
    df_hobo_inputs <- build_hobo_inputs(CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove)
    cat(sprintf("    %d HOBO plots retained\n", nrow(df_hobo_inputs)))
    
    hobo_scs <- if (FLAGS$RUN_H1_FORWARD) c(list(REF = ref_sc), scenarios_h1_forward(df_sample)) else list(REF = ref_sc)
    
    hobo_res <- validate_scenarios_at_hobos(df_hobo_inputs, df_hobo_daily, hobo_scs, parent_dir = CFG$out_hobo, df_macro = df_macro, date_seq = CFG$date_seq)
    print(hobo_res$metrics)
    write.csv(hobo_res$metrics, "outputs/hobo/hobo_metrics.csv", row.names = FALSE)

    p_hobo_val <- plot_hobo_validation(hobo_res$metrics)
    save_plot(p_hobo_val, "outputs/hobo/hobo_validation_rmse.png")
    print(p_hobo_val)

    fw_names  <- vapply(scenarios_h1_forward(df_sample), `[[`, "", "name")
    rep_hobos <- pick_representative_hobos(df_hobo_daily)

    df_real <- hobo_res$daily %>% filter(scenario == fw_names["Full_real"]) %>% mutate(Tmax_micro = Delta_sim + Tmax_macro)
    df_unif <- hobo_res$daily %>% filter(scenario == fw_names["LAI_Hmax_fCover"]) %>% mutate(Tmax_micro = Delta_sim + Tmax_macro)

    p_ts <- plot_timeseries_faceted(rep_hobos, df_hobo_daily, df_real, df_unif, df_macro)
    save_plot(p_ts, "outputs/hobo/hobo_timeseries_representatives.png", width = 12, height = 10)
    print(p_ts)
  }
  
  # ---- OPTIONAL: Sentinel-2 + FORMS-H annex ---------------------------------
  if (FLAGS$RUN_S2_ANNEX) {
    cat("[annex] Running S2 + FORMS-H naive scenario...\n")
    annex    <- make_s2_formsh_scenario(df_sample, CFG$in_dir, CFG$agg_factor)
    p_annex  <- run_musica_scenario(annex$df_sample, annex$scenario, CFG$out_s2_annex)
    df_annex <- extract_deltatmax_scenario(p_annex, df_macro, CFG$date_seq)
    
    cat("    [annex] Computing metrics S2 vs Full LiDAR...\n")
    df_ref_compare <- df_all_scenarios %>% filter(scenario == ref_sc$name) %>% dplyr::select(x, y, date, Delta_ref = Delta_Tmax)
    df_s2_matched  <- df_annex %>% inner_join(df_ref_compare, by = c("x", "y", "date"))
    
    s2_r2   <- cor(df_s2_matched$Delta_Tmax, df_s2_matched$Delta_ref, use = "complete.obs")^2
    s2_rmse <- sqrt(mean((df_s2_matched$Delta_Tmax - df_s2_matched$Delta_ref)^2, na.rm = TRUE))
    s2_bias <- mean(df_s2_matched$Delta_Tmax - df_s2_matched$Delta_ref, na.rm = TRUE)
    
    stats_label_s2 <- sprintf("R² = %.2f\nRMSE = %.2f°C\nBias = %.2f°C", s2_r2, s2_rmse, s2_bias)
    
    df_s2_eval <- bind_rows(df_annex, df_all_scenarios %>% filter(scenario == ref_sc$name))
    df_plot_s2 <- df_s2_eval %>% mutate(scenario_label = ifelse(scenario == ref_sc$name, "LiDAR (Full 3D Real)", "Sentinel-2 + FORMS-H (2D Proxy)"))
    
    p_s2_density <- ggplot(df_plot_s2, aes(x = Delta_Tmax, fill = scenario_label, colour = scenario_label)) +
      geom_density(alpha = 0.3, linewidth = 0.8) +
      scale_fill_manual(values   = c("LiDAR (Full 3D Real)" = "#31688e", "Sentinel-2 + FORMS-H (2D Proxy)" = "#fca50a")) +
      scale_colour_manual(values = c("LiDAR (Full 3D Real)" = "#31688e", "Sentinel-2 + FORMS-H (2D Proxy)" = "#fca50a")) +
      labs(title = "Annex H3 — Optical Proxy vs Structural Reality (Density)",
           subtitle = sprintf("R² = %.2f | RMSE = %.2f°C | Bias = %.2f°C", s2_r2, s2_rmse, s2_bias),
           x = expression(Delta * T[max] ~ (degree*C)), y = "Density", fill = "Input Source", colour = "Input Source") +
      theme(legend.position = "bottom")
    
    save_plot(p_s2_density, "outputs/s2_annex/s2_density.png")
    print(p_s2_density)

    axis_lims_s2 <- c(min(c(df_s2_matched$Delta_ref, df_s2_matched$Delta_Tmax), na.rm = TRUE), max(c(df_s2_matched$Delta_ref, df_s2_matched$Delta_Tmax), na.rm = TRUE))
    
    p_s2_paired <- ggplot(df_s2_matched, aes(x = Delta_ref, y = Delta_Tmax)) +
      geom_hex(bins = 100) + scale_fill_gradient(low = "grey80", high = "midnightblue", trans = "log10", name = "Count\n(log)") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black", linewidth = 0.8) +
      geom_smooth(method = "lm", colour = "#fca50a", se = FALSE, linewidth = 1.2) +
      coord_fixed(xlim = axis_lims_s2, ylim = axis_lims_s2) +
      annotate("text", x = -Inf, y = Inf, label = stats_label_s2, hjust = -0.1, vjust = 1.2, fontface = "bold", size = 5, colour = "black") +
      labs(title = "Annex H3 — Per-plot, per-day ΔTmax: S2 Proxy vs LiDAR", subtitle = "Paired comparison: impact of 2D optical proxy on spatial prediction",
           x = expression("LiDAR (Full 3D Real)" ~ Delta * T[max] ~ (degree*C)), y = expression("Sentinel-2 + FORMS-H" ~ Delta * T[max] ~ (degree*C))) +
      theme_bw(base_size = 14) + theme(legend.position = "right", plot.title = element_text(face = "bold"))
    
    save_plot(p_s2_paired, "outputs/s2_annex/s2_paired.png")
    print(p_s2_paired)

    writeLines(sprintf("R2: %.4f\nRMSE: %.4f deg C\nBias: %.4f deg C", s2_r2, s2_rmse, s2_bias),
               "outputs/s2_annex/s2_metrics.txt")
  }

  cat("\n[done]\n")
  invisible(NULL)
}

if (interactive() || identical(Sys.getenv("RUN_MAIN"), "1")) {
  main()
}