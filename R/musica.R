# ==============================================================================
# Chapter 1 — MuSICA batch runner and post-processing (NetCDF → daily ΔTmax)
# ==============================================================================

#' Copy outputs from existing scenarios (forward/LOO) into the factorial folder.
#'
#' Avoids re-running simulations already available elsewhere.
#'
#' @param preload_map Named character vector: name = factorial target name,
#'   value = path to source directory (forward or LOO).
#' @param out_dir     Factorial output directory (CFG$out_h1_factorial).
#' @return Character vector of factorial names that were successfully preloaded.
preload_factorial_from_existing <- function(preload_map, out_dir) {
  # WARNING: this function copies NetCDF directories by name. If scenario
  # parameters (e.g. hmax_mean_val, lai_mean) change after the source was
  # produced, the copied NCs will be silently stale. The source directory name
  # encodes the scenario name only, not its parameter values. To detect
  # staleness: delete out_dir entirely before a full pipeline re-run whenever
  # scenario definitions change.
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  preloaded <- character(0)
  for (nm in names(preload_map)) {
    src <- preload_map[[nm]]
    dst <- file.path(out_dir, nm)
    if (dir.exists(dst)) {
      n_nc <- length(list.files(dst, pattern = "\\.nc$"))
      cat(sprintf("  [factorial preload] %s — déjà présent (%d NC), skip\n", nm, n_nc))
      preloaded <- c(preloaded, nm)
    } else if (dir.exists(src)) {
      n_src <- length(list.files(src, pattern = "\\.nc$"))
      ok <- file.copy(src, out_dir, recursive = TRUE)
      if (ok) {
        file.rename(file.path(out_dir, basename(src)), dst)
        cat(sprintf("  [factorial preload] %s ← copié depuis %s (%d NC)\n",
                    nm, basename(src), n_src))
        preloaded <- c(preloaded, nm)
      } else {
        cat(sprintf("  [factorial preload] ECHEC copie %s\n", nm))
      }
    } else {
      cat(sprintf("  [factorial preload] %s — source absente (%s), sera simulé\n",
                  nm, basename(src)))
    }
  }
  preloaded
}

#' Run a single MuSICA simulation for one plot × one scenario.
#'
#' @param plot_row     Single-row dataframe with structural covariates.
#' @param scenario     Scenario list (from scenarios.R).
#' @param out_nc_file  Path for the output NetCDF.
#' @param forcing_file Path to the ERA5 forcing NetCDF.
#' @param musica_cmd   Shell command to call MuSICA.
#' @return Invisible path to out_nc_file, or NULL on error.
run_musica_one <- function(plot_row, scenario, out_nc_file, forcing_file, musica_cmd,
                           musica_nml       = "musica.nml",
                           musica_variables = "./variables.csv",
                           extra_setup      = list()) {
  if (file.exists(out_nc_file)) return(invisible(out_nc_file))
  lai    <- scenario$lai_fn(plot_row)
  hmax   <- scenario$hmax_fn(plot_row)
  fcover <- scenario$fcover_fn(plot_row)
  if (is.na(lai) || is.na(hmax) || is.na(fcover)) return(invisible(NULL))
  # MuSICA convention: the model expects TOTAL (two-sided) leaf area = 2 × LAI.
  # Doubling `lai` here propagates to the LAD profile (lad_fn), the parametric
  # phenology (LAI_max_per_cohort) and lai_max_per_cohort. The S2-driven phenology
  # returns absolute Leaf_area and is doubled explicitly below.
  lai <- 2 * lai
  allom <- scenario$lad_fn(plot_row, hmax = hmax, lai = lai)

  # Phenology: S2-driven if phenology_fn provided, otherwise parametric fallback.
  # phenology_fn may return NULL for a given plot (missing S2 data) — also falls back.
  phenology <- NULL
  if (!is.null(scenario$phenology_fn)) {
    phenology <- scenario$phenology_fn(plot_row)
    if (!is.null(phenology)) {
      la_cols <- grep("^Leaf_area", names(phenology), value = TRUE)
      phenology[la_cols] <- lapply(phenology[la_cols], function(x) 2 * x)  # -> 2 × LAI
    }
  }
  if (is.null(phenology)) {
    phenology <- calc_phenology(
      list.year = c(2021, 2022), nleafage = 1,
      budburst_date = 115, leaf_age_max_in = 0.56,
      relative_age_firstmax = 0.10, relative_age_lastmax = 0.75,
      LAI_max_per_cohort = lai
    )
  }

  # forcing_height : si l'appelant le fournit via extra_setup on l'utilise,
  # sinon valeur dynamique hmax+2 (legacy). Permet de tester fh=30 constant.
  fh <- if (!is.null(extra_setup$forcing_height))
    extra_setup$forcing_height
  else
    hmax + 2

  setup_param <- c(
    list("clumping_factor"   = fcover,
         "forcing_filename"  = forcing_file,
         "history_filename"  = out_nc_file,
         "forcing_height"    = fh),
    extra_setup[setdiff(names(extra_setup), "forcing_height")]
  )

  tryCatch({
    callmusica(
      musica.file      = musica_nml,
      musica.variables = musica_variables,
      musica.param     = list(setupctl = setup_param),
      leaf.param       = list(musica_veg1 = list("phenology" = phenology, "allometry" = allom, "leafphenologyctl" = list("lai_max_per_cohort" = lai),
                                                  "leafallometryctl" = list("canopy_height_top" = hmax, "canopy_height_bottom" = 1), "leafmusicactl" = list("canopy_height_top" = hmax))),
      musica.cmd = musica_cmd, keep.tmp = FALSE, out.netcdf = TRUE, out.df = FALSE
    )
  }, error = function(e) cat(sprintf("ERROR on %s: %s\n", basename(out_nc_file), e$message)))
  invisible(out_nc_file)
}

#' Run all plots for a single scenario, writing NetCDFs to out_dir.
#'
#' @param df_sample      cLHS sample dataframe.
#' @param scenario       Scenario list (from scenarios.R).
#' @param out_dir        Output directory for this scenario's NetCDF files.
#' @param forcing_file   ERA5 forcing file path (default CFG$forcing_file).
#' @param musica_cmd     MuSICA shell command (default CFG$musica_cmd).
#' @param id_prefix      Prefix for simulation IDs (default "Sim").
#' @param force          Re-run even if output exists (default FALSE).
#' @param skip_existing  Skip the whole directory if any NC files exist (default FALSE).
#' @return Invisible out_dir.
run_musica_scenario <- function(df_sample, scenario, out_dir, forcing_file = CFG$forcing_file, musica_cmd = CFG$musica_cmd, id_prefix = "Sim", force = FALSE, skip_existing = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, pattern = "\\.nc$")
  if (!force && skip_existing && length(existing) > 0) {
    cat(sprintf("  [skip] %s : %d NC existants — répertoire conservé tel quel\n",
                scenario$name, length(existing)))
    return(invisible(out_dir))
  }
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

#' Run a list of scenarios sequentially, returning their output paths.
#'
#' @param df_sample      cLHS sample dataframe.
#' @param scenarios      Named list of scenario lists.
#' @param parent_dir     Parent directory; each scenario gets a sub-directory.
#' @param force          Re-run existing outputs (default FALSE).
#' @param skip_existing  Skip populated scenario directories (default FALSE).
#' @return Named character vector of output directory paths.
run_scenarios <- function(df_sample, scenarios, parent_dir, force = FALSE,
                          skip_existing = getOption("skip_existing_scenarios", FALSE)) {
  paths <- character(length(scenarios)); names(paths) <- names(scenarios)
  for (nm in names(scenarios)) {
    out_d <- file.path(parent_dir, scenarios[[nm]]$name)
    run_musica_scenario(df_sample, scenarios[[nm]], out_d,
                        force = force, skip_existing = skip_existing)
    paths[nm] <- out_d
  }
  paths
}

# ---- Post-processing: NetCDF → daily ΔTmax -----------------------------------

#' Interpolate Tair_z at a target absolute height (m above ground).
#'
#' MuSICA stocke Tair_z sur une grille verticale propre à chaque plot :
#' z(nair) = relative_height(nair) * Hmax. Comme nair=1 est à 0.028*Hmax,
#' sa hauteur absolue varie de 0.14 m (Hmax=5) à 1 m (Hmax=35) selon le plot.
#' Pour des comparaisons inter-plots cohérentes et alignées sur les HOBO,
#' on extrait Tair à une **hauteur absolue fixe** (1 m par défaut) en
#' interpolant linéairement entre les 2 couches encadrant cette hauteur.
#'
#' @param nc        Open ncdf4 handle (déjà ouvert par l'appelant).
#' @param z_target  Hauteur cible en mètres au-dessus du sol (défaut 1.0).
#' @return Tibble (time, Tair_sim en degrés C) ou NULL si lecture impossible.
get_tair_at_z <- function(nc, z_target = 1.0) {
  rh   <- try(ncdf4::ncvar_get(nc, "relative_height"),   silent = TRUE)
  hmax <- try(ncdf4::ncvar_get(nc, "veget_height_top"),  silent = TRUE)
  if (inherits(rh, "try-error") || inherits(hmax, "try-error")) return(NULL)

  raw <- try(get_variable(nc, "Tair_z"), silent = TRUE)
  if (inherits(raw, "try-error") || is.null(raw)) return(NULL)

  # Fallback nair=1 si Hmax non exploitable (Null_baseline canopée nulle, NaN, etc.)
  fallback_nair1 <- function() {
    raw %>% dplyr::filter(nair == 1) %>%
      dplyr::transmute(time, Tair_sim = Tair_z - 273.15)
  }

  hmax <- hmax[1]
  if (is.na(hmax) || !is.finite(hmax) || hmax <= 0) return(fallback_nair1())

  z_layers <- rh * hmax
  if (any(!is.finite(z_layers))) return(fallback_nair1())

  if (z_target <= z_layers[1]) {
    i_lo <- 1L; i_hi <- 1L; w <- 0.0
  } else if (z_target >= z_layers[length(z_layers)]) {
    i_lo <- length(z_layers); i_hi <- length(z_layers); w <- 0.0
  } else {
    i_lo <- max(which(z_layers <= z_target))
    i_hi <- i_lo + 1L
    w    <- (z_target - z_layers[i_lo]) / (z_layers[i_hi] - z_layers[i_lo])
  }

  if (i_hi == i_lo) {
    raw %>%
      dplyr::filter(nair == i_lo) %>%
      dplyr::transmute(time, Tair_sim = Tair_z - 273.15)
  } else {
    lo <- raw %>% dplyr::filter(nair == i_lo) %>% dplyr::select(time, Tair_lo = Tair_z)
    hi <- raw %>% dplyr::filter(nair == i_hi) %>% dplyr::select(time, Tair_hi = Tair_z)
    dplyr::inner_join(lo, hi, by = "time") %>%
      dplyr::mutate(Tair_sim = (1 - w) * Tair_lo + w * Tair_hi - 273.15) %>%
      dplyr::select(time, Tair_sim)
  }
}

#' Extract daily ΔTmax from a single MuSICA NetCDF output file.
#'
#' @param nc_file   Path to a MuSICA NetCDF output.
#' @param df_macro  Daily macroclimate (from extract_macro_daily).
#' @param date_seq  Dates of interest.
#' @param z_target  Hauteur absolue (m) à laquelle extraire Tair. Défaut
#'   `CFG$tair_target_height` si défini, sinon 1.0 m. Cohérent inter-plots
#'   et aligné sur la hauteur des sondes HOBO.
#' @return Tibble with date, Tmax_micro, Delta_Tmax, or NULL on failure.
# ============================================================================
# DEPRECATED for Chapter 1 — DO NOT USE IN ANY SHIPPED ANALYSIS (2026-07-31).
# This is the LEGACY ΔTmax extractor. It takes the sub-canopy daily maximum
# INDEPENDENTLY of the macroclimatic one and applies a -2 h clock shift. The
# chapter's canonical metric is the time-matched convention B of
# R/dtmax_convention.R (sub-canopy read AT the hour of the macro daily maximum,
# no clock shift): use macro_ref() + micro_hourly_at() + delta_tmax() instead.
# Kept only so the superseded h1_*/lovb_*/v9-v11 scripts still run. Fig. B1 was
# the last shipped figure on this extractor and was converted on 2026-07-31.
# ============================================================================
extract_deltatmax_one <- function(nc_file, df_macro, date_seq,
                                   z_target = if (exists("CFG") && !is.null(CFG$tair_target_height)) CFG$tair_target_height else 1.0) {
  nc <- try(nc_open(nc_file), silent = TRUE); if (inherits(nc, "try-error")) return(NULL)
  res <- get_tair_at_z(nc, z_target = z_target)
  nc_close(nc)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res %>%
    mutate(time = time - hours(2), date = as.Date(time)) %>%
    filter(date %in% date_seq) %>% group_by(date) %>%
    summarise(Tmax_micro = max(Tair_sim, na.rm = TRUE), .groups = "drop") %>%
    inner_join(df_macro, by = "date") %>% mutate(Delta_Tmax = Tmax_micro - Tmax_macro)
}

#' Extract hourly Tair_z at a fixed height + per-plot regression vs hourly ERA5.
#'
#' Reproduces the legacy `main_compare_era5_ign.R` buf/amp metric: fits
#' `Tair_sim ~ Tair_era5` on hourly series (filtered to date_seq), returns the
#' slope, log_slope = log(abs(slope)), and mean offset. slope < 1 ⇒ buffering,
#' slope > 1 ⇒ amplification.
#'
#' @param nc_file        Path to a MuSICA NetCDF output.
#' @param era5_hourly    Hourly ERA5 dataframe (cols: time POSIXct UTC, Tair_era5).
#' @param date_seq       Dates of interest.
#' @param z_target       Height above ground (m). Default `CFG$tair_target_height`.
#' @param use_nair1      If TRUE, use nair==1 (raw lowest air node) instead of
#'                       z_target interpolation, matching legacy
#'                       `main_compare_era5_ign.R` exactly.
#' @param time_shift_hr  Hours subtracted from MuSICA time. Legacy = 0 (UTC↔UTC).
#'                       Default 0 to recover the legacy 45/8 buf/amp split.
#' @return Tibble (slope, log_slope, dT_mean, n) or NULL if extraction fails.
extract_hourly_slope_one <- function(nc_file, era5_hourly, date_seq,
                                      z_target = if (exists("CFG") && !is.null(CFG$tair_target_height)) CFG$tair_target_height else 1.0,
                                      use_nair1 = TRUE,
                                      time_shift_hr = 0) {
  nc <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  on.exit(nc_close(nc), add = TRUE)
  if (isTRUE(use_nair1)) {
    raw <- tryCatch(musica.tools::get_variable(nc, "Tair_z"),
                    error = function(e) NULL)
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    res <- raw %>%
      dplyr::filter(nair == 1) %>%
      dplyr::transmute(time, Tair_sim = Tair_z - 273.15)
  } else {
    res <- get_tair_at_z(nc, z_target = z_target)
  }
  if (is.null(res) || nrow(res) == 0) return(NULL)
  df <- res %>%
    mutate(time = time - lubridate::hours(time_shift_hr),
           time = lubridate::floor_date(time, "hour"),
           date = as.Date(time)) %>%
    filter(date %in% date_seq) %>%
    select(time, Tair_sim)
  m <- merge(df, era5_hourly, by = "time")
  if (nrow(m) < 24) return(NULL)
  fit <- lm(Tair_sim ~ Tair_era5, data = m)
  sl <- as.numeric(coef(fit)[2])
  tibble::tibble(slope     = sl,
                 log_slope = log(abs(sl)),
                 dT_mean   = mean(m$Tair_sim - m$Tair_era5, na.rm = TRUE),
                 n         = nrow(m))
}

#' Build the hourly ERA5 Tair dataframe consumed by extract_hourly_slope_one.
#'
#' @param forcing_file  Path to ERA5 NetCDF forcing.
#' @param date_seq      Dates of interest.
#' @return Tibble with columns time (POSIXct UTC, hour-floored) and Tair_era5 (°C).
build_era5_hourly <- function(forcing_file, date_seq) {
  nc <- nc_open(forcing_file); on.exit(nc_close(nc), add = TRUE)
  tair <- ncvar_get(nc, "Tair") - 273.15
  t    <- force_utc_nc(forcing_file, "time")
  data.frame(time = lubridate::floor_date(t, "hour"), Tair_era5 = tair) %>%
    filter(as.Date(time) %in% date_seq) %>%
    distinct(time, .keep_all = TRUE)
}

#' Extract daily ΔTmax for all NetCDF files in a scenario directory.
#'
#' @param scenario_dir  Directory containing MuSICA NetCDF outputs.
#' @param df_macro      Daily macroclimate dataframe.
#' @param date_seq      Dates of interest.
#' @return Stacked tibble with x, y, scenario columns added.
extract_deltatmax_scenario <- function(scenario_dir, df_macro, date_seq) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE); if (length(nc_files) == 0) return(NULL)
  map_df(nc_files, function(f) {
    res <- extract_deltatmax_one(f, df_macro, date_seq); if (is.null(res)) return(NULL)
    res %>% mutate(x = as.numeric(str_extract(basename(f), "(?<=_X)\\d+")), y = as.numeric(str_extract(basename(f), "(?<=_Y)\\d+")), scenario = basename(scenario_dir))
  })
}

#' Extract daily ΔTmax across all scenarios.
#'
#' @param scenario_paths Named character vector of scenario directory paths.
#' @param df_macro       Daily macroclimate dataframe.
#' @param date_seq       Dates of interest.
#' @return Stacked tibble with all scenarios.
extract_all_scenarios <- function(scenario_paths, df_macro, date_seq) {
  map_df(scenario_paths, extract_deltatmax_scenario, df_macro = df_macro, date_seq = date_seq)
}

#' Join daily scenario results with the cLHS sample metadata.
#'
#' NC filenames store round(plot$x) as integers; terra returns cell centroids as
#' floats (often half-integers for 25m grids). round() on the df_sample side
#' aligns both coordinate sets so the inner_join matches all 400 plots.
#'
#' @param df_daily   Stacked scenario results (from extract_all_scenarios).
#' @param df_sample  cLHS sample dataframe.
#' @return Tibble with structural covariates attached.
join_with_sample <- function(df_daily, df_sample) {
  df_meta <- df_sample %>%
    mutate(x = round(x), y = round(y)) %>%
    dplyr::select(x, y, Cluster, LAI, Hmax, fCover, VCI, starts_with("LAD_Layer_"))
  df_daily %>% inner_join(df_meta, by = c("x", "y"))
}
