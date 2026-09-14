# ==============================================================================
# Chapter 1 — Vertical temperature profiles (H2 — by cluster median)
#
# NetCDF MuSICA structure:
#   Tair_z          : dims [nair=15 × time]
#   relative_height : vector [nair] of relative heights (0→1)
#   veget_height_top: absolute canopy height (time series, constant in practice)
#   Absolute Z axis = relative_height × veget_height_top[1]
# ==============================================================================

#' Extraire le profil vertical de Tmax pour un NetCDF MuSICA.
#'
#' @param nc_file   Chemin vers le fichier NetCDF.
#' @param date_seq  Dates d'intérêt (filtre sur l'été).
#' @param summary   "mean_summer" | "canicule_date" | "median_date"
#' @return Tibble (height_m, Tmax_z) ou NULL si échec.
extract_vertical_tmax_profile <- function(nc_file, date_seq,
                                           summary = "mean_summer") {
  if (is.na(nc_file) || !file.exists(nc_file)) return(NULL)

  nc <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)

  rel_h <- try(ncvar_get(nc, "relative_height"), silent = TRUE)
  h_top <- try(ncvar_get(nc, "veget_height_top"), silent = TRUE)
  raw   <- try(get_variable(nc, "Tair_z"), silent = TRUE)
  nc_close(nc)

  if (inherits(raw, "try-error") || is.null(raw)) return(NULL)

  hmax_nc <- if (!inherits(h_top, "try-error")) h_top[1] else NA_real_

  df <- raw %>%
    mutate(
      Tair_sim = Tair_z - 273.15,
      time     = time - hours(2),
      date     = as.Date(time)
    ) %>%
    filter(date %in% date_seq) %>%
    group_by(nair, date) %>%
    summarise(Tmax_daily = max(Tair_sim, na.rm = TRUE), .groups = "drop")

  df_out <- switch(summary,
    "mean_summer" = df %>%
      group_by(nair) %>%
      summarise(Tmax_z = mean(Tmax_daily, na.rm = TRUE), .groups = "drop"),

    "canicule_date" = {
      d_max <- df %>%
        group_by(date) %>%
        summarise(m = max(Tmax_daily, na.rm = TRUE), .groups = "drop") %>%
        slice_max(m, n = 1, with_ties = FALSE) %>%
        pull(date)
      df %>%
        filter(date == d_max) %>%
        dplyr::select(nair, Tmax_z = Tmax_daily)
    },

    "median_date" = {
      d_med <- df %>%
        group_by(date) %>%
        summarise(m = max(Tmax_daily, na.rm = TRUE), .groups = "drop") %>%
        arrange(m) %>%
        slice(ceiling(n() / 2)) %>%
        pull(date)
      df %>%
        filter(date == d_med) %>%
        dplyr::select(nair, Tmax_z = Tmax_daily)
    }
  )

  if (!inherits(rel_h, "try-error") && !is.na(hmax_nc)) {
    df_out$height_m <- rel_h[df_out$nair] * hmax_nc
  } else {
    n_layers <- max(df_out$nair, na.rm = TRUE)
    df_out$height_m <- seq(0, hmax_nc, length.out = n_layers)[df_out$nair]
  }

  df_out
}

#' Identifier le plot le plus proche de la médiane structurelle de chaque cluster.
#'
#' Distance euclidienne normalisée sur (LAI, Hmax, fCover) centrées
#' par la médiane du cluster.
#'
#' @param df_sample  cLHS avec colonnes (x, y, Cluster, LAI, Hmax, fCover).
#' @return Tibble 1 ligne × cluster : (Cluster, x, y, LAI, Hmax, fCover, plot_id).
select_cluster_median_plots <- function(df_sample) {
  if (!"Cluster" %in% names(df_sample)) {
    warning("[select_cluster_median_plots] colonne 'Cluster' absente de df_sample.")
    return(invisible(NULL))
  }

  df_sample %>%
    filter(!is.na(Cluster)) %>%
    group_by(Cluster) %>%
    mutate(
      sd_LAI    = sd(LAI,    na.rm = TRUE),
      sd_Hmax   = sd(Hmax,   na.rm = TRUE),
      sd_fCover = sd(fCover, na.rm = TRUE),
      d_med = sqrt(
        ((LAI    - median(LAI,    na.rm = TRUE)) / pmax(sd_LAI,    1e-9))^2 +
        ((Hmax   - median(Hmax,   na.rm = TRUE)) / pmax(sd_Hmax,   1e-9))^2 +
        ((fCover - median(fCover, na.rm = TRUE)) / pmax(sd_fCover, 1e-9))^2
      )
    ) %>%
    slice_min(d_med, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(plot_id = sprintf("X%d_Y%d", round(x), round(y))) %>%
    dplyr::select(Cluster, x, y, LAI, Hmax, fCover, plot_id)
}

#' Visualiser le profil LAD réel + ligne Hmax pour un plot donné.
#'
#' @param plot_row   Une ligne de tibble avec colonnes (x, y, Hmax, Cluster,
#'                   plot_id, fCover) et des colonnes LAD_Layer_*.
#' @param out_path   Chemin de sortie PNG.
#' @param annotate   Texte optionnel ajouté en sous-titre à la place du défaut.
#' @return Invisible ggplot object.
plot_lad_profile_for_plot <- function(plot_row, out_path, annotate = NULL) {
  lad_cols <- grep("^LAD_Layer_", names(plot_row), value = TRUE)
  if (length(lad_cols) == 0) {
    warning("[plot_lad_profile_for_plot] Aucune colonne LAD_Layer_ dans plot_row.")
    return(invisible(NULL))
  }

  hmax <- as.numeric(plot_row$Hmax[1])
  n    <- length(lad_cols)
  z    <- seq(0, hmax, length.out = n)
  lad  <- as.numeric(plot_row[1, lad_cols])

  df_lad <- tibble(z_m = z, LAD = lad)

  p <- ggplot(df_lad, aes(x = LAD, y = z_m)) +
    geom_path(linewidth = 1.2, colour = "#31688e") +
    geom_point(size = 2, alpha = 0.7, colour = "#31688e") +
    geom_hline(yintercept = hmax, linetype = "dashed",
               colour = "firebrick", linewidth = 0.9) +
    annotate("text",
             x     = max(lad, na.rm = TRUE) * 0.7,
             y     = hmax,
             label = sprintf("Hmax = %.1f m", hmax),
             vjust = -0.5, colour = "firebrick",
             fontface = "bold", size = 4) +
    labs(
      title    = sprintf("LAD profile — Plot %s (Cluster %s)",
                         plot_row$plot_id[1], plot_row$Cluster[1]),
      subtitle = if (!is.null(annotate)) annotate else
                 sprintf("LAI = %.2f | Hmax = %.1f m | fCover = %.2f",
                         plot_row$LAI[1], hmax, plot_row$fCover[1]),
      x = expression("LAD" ~ (m^2 ~ m^{-3})),
      y = "Height (m)"
    ) +
    theme_bw(base_size = 12)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 6, height = 7)
  cat(sprintf("  [plot_lad_profile_for_plot] Sauvegardé : %s\n", out_path))
  invisible(p)
}

#' Vérifier la robustesse du pic de Tmax proche du sol en Cluster 1
#' sur plusieurs plots (profils "canicule_date" Real vs Uniform LAD).
#'
#' @param df_sample       Tibble cLHS (Cluster, x, y, LAI, Hmax, fCover).
#' @param real_lad_dir    Dossier NetCDF H2_real_LAD.
#' @param uniform_lad_dir Dossier NetCDF H2_uniform_LAD.
#' @param date_seq        Dates d'intérêt.
#' @param out_dir         Dossier de sortie.
#' @param n_plots         Nombre de plots C1 à comparer (default 3).
#' @return Invisible list(plot, data).
verify_cluster1_peak <- function(df_sample, real_lad_dir, uniform_lad_dir,
                                  date_seq, out_dir = "outputs/h2",
                                  n_plots = 3) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!"Cluster" %in% names(df_sample)) {
    warning("[verify_cluster1_peak] Colonne 'Cluster' absente — abandon.")
    return(invisible(NULL))
  }

  c1_all <- df_sample %>%
    filter(as.character(Cluster) == "1") %>%
    mutate(
      sd_LAI    = sd(LAI,    na.rm = TRUE),
      sd_Hmax   = sd(Hmax,   na.rm = TRUE),
      sd_fCover = sd(fCover, na.rm = TRUE),
      d_med     = sqrt(
        ((LAI    - median(LAI,    na.rm = TRUE)) / pmax(sd_LAI,    1e-9))^2 +
        ((Hmax   - median(Hmax,   na.rm = TRUE)) / pmax(sd_Hmax,   1e-9))^2 +
        ((fCover - median(fCover, na.rm = TRUE)) / pmax(sd_fCover, 1e-9))^2
      )
    ) %>%
    arrange(d_med)

  if (nrow(c1_all) == 0) {
    warning("[verify_cluster1_peak] Aucun plot Cluster 1 trouvé — abandon.")
    return(invisible(NULL))
  }

  idx_sel <- c(2, 5, 10)
  idx_sel <- idx_sel[idx_sel <= nrow(c1_all)]
  if (length(idx_sel) == 0) idx_sel <- seq_len(min(n_plots, nrow(c1_all)))
  c1_sel  <- c1_all[idx_sel, ] %>%
    mutate(plot_id = sprintf("X%d_Y%d", round(x), round(y)))

  find_nc <- function(dir, x, y) {
    pat  <- sprintf("X%d_Y%d", round(x), round(y))
    hits <- list.files(dir, pattern = pat, full.names = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  c1_meta <- c1_sel %>%
    rowwise() %>%
    mutate(
      nc_real    = find_nc(real_lad_dir,    x, y),
      nc_uniform = find_nc(uniform_lad_dir, x, y)
    ) %>%
    ungroup() %>%
    filter(!is.na(nc_real) & !is.na(nc_uniform))

  if (nrow(c1_meta) == 0) {
    warning("[verify_cluster1_peak] Aucun NetCDF valide pour Cluster 1 — abandon.")
    return(invisible(NULL))
  }

  df_profiles <- purrr::pmap_dfr(
    c1_meta,
    function(Cluster, x, y, LAI, Hmax, fCover, plot_id, nc_real, nc_uniform, ...) {
      prof_real <- extract_vertical_tmax_profile(nc_real,    date_seq, "canicule_date")
      prof_unif <- extract_vertical_tmax_profile(nc_uniform, date_seq, "canicule_date")
      bind_rows(
        if (!is.null(prof_real)) prof_real %>% mutate(scenario = "Real LAD"),
        if (!is.null(prof_unif)) prof_unif %>% mutate(scenario = "Uniform LAD")
      ) %>%
        mutate(Cluster = as.character(Cluster),
               plot_id = plot_id,
               LAI     = LAI,
               Hmax    = Hmax)
    }
  )

  if (nrow(df_profiles) == 0) {
    warning("[verify_cluster1_peak] Aucun profil extrait pour Cluster 1.")
    return(invisible(NULL))
  }

  df_profiles <- df_profiles %>%
    mutate(panel_label = sprintf("Plot %s\nHmax=%.1fm | LAI=%.1f",
                                  plot_id, Hmax, LAI))

  df_hmax_lines <- df_profiles %>%
    group_by(panel_label) %>%
    summarise(Hmax = first(Hmax), .groups = "drop")

  p <- ggplot(df_profiles, aes(x = Tmax_z, y = height_m, colour = scenario)) +
    geom_path(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(data      = df_hmax_lines,
               aes(yintercept = Hmax),
               linetype  = "dashed", colour = "firebrick",
               linewidth = 0.8, inherit.aes = FALSE) +
    facet_wrap(~ panel_label, nrow = 1, scales = "free_x") +
    scale_colour_manual(
      values = c("Real LAD" = "#31688e", "Uniform LAD" = "#d8576b")
    ) +
    labs(
      title    = "Ground-peak check — Cluster 1 (heatwave date)",
      subtitle = sprintf("%d Cluster 1 plots ordered by distance to median",
                         nrow(c1_meta)),
      x      = expression(T[max] ~ "simulated" ~ (degree*C)),
      y      = "Height (m)",
      colour = "Scenario"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold"))

  png_path <- file.path(out_dir, "h2_verify_cluster1_peak.png")
  csv_path <- file.path(out_dir, "h2_verify_cluster1_peak.csv")

  save_plot(p, png_path, width = 14, height = 6)
  write.csv(df_profiles, csv_path, row.names = FALSE)
  cat(sprintf("  [verify_cluster1_peak] PNG : %s\n", png_path))
  cat(sprintf("  [verify_cluster1_peak] CSV : %s\n", csv_path))

  invisible(list(plot = p, data = df_profiles))
}

#' Produire la figure en 3 panneaux (un par cluster) : Tmax vs hauteur,
#' Real LAD vs Uniform LAD.
#'
#' @param df_sample       cLHS (besoin de Cluster, LAI, Hmax, fCover, x, y).
#' @param real_lad_dir    Dossier NetCDF du scénario H2_real_LAD.
#' @param uniform_lad_dir Dossier NetCDF du scénario H2_uniform_LAD.
#' @param date_seq        Dates d'intérêt (filtre été).
#' @param out_dir         Dossier de sortie PNG/CSV.
#' @param summary_mode    "mean_summer" | "canicule_date" | "median_date".
#' @return Invisiblement list(plot, data, medians).
plot_vertical_tmax_profiles <- function(df_sample, real_lad_dir,
                                         uniform_lad_dir, date_seq,
                                         out_dir = "outputs/h2",
                                         summary_mode = "mean_summer") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (!"Cluster" %in% names(df_sample)) {
    warning("[plot_vertical_tmax_profiles] colonne 'Cluster' absente — abandon.")
    return(invisible(NULL))
  }

  medians <- select_cluster_median_plots(df_sample)
  if (is.null(medians) || nrow(medians) == 0) {
    warning("[plot_vertical_tmax_profiles] Aucun plot médian sélectionné.")
    return(invisible(NULL))
  }
  cat(sprintf("  [vertical profile] %d plots médians :\n", nrow(medians)))
  print(as.data.frame(medians[, c("Cluster", "plot_id", "LAI", "Hmax", "fCover")]))

  find_nc <- function(dir, x, y) {
    pat  <- sprintf("X%d_Y%d", round(x), round(y))
    hits <- list.files(dir, pattern = pat, full.names = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  profiles_meta <- medians %>%
    rowwise() %>%
    mutate(
      nc_real    = find_nc(real_lad_dir,    x, y),
      nc_uniform = find_nc(uniform_lad_dir, x, y)
    ) %>%
    ungroup()

  missing <- profiles_meta %>%
    filter(is.na(nc_real) | is.na(nc_uniform))
  if (nrow(missing) > 0) {
    warning(sprintf(
      "[plot_vertical_tmax_profiles] NetCDF manquants pour %d plot(s) : %s",
      nrow(missing), paste(missing$plot_id, collapse = ", ")
    ))
    print(missing[, c("Cluster", "plot_id", "nc_real", "nc_uniform")])
    profiles_meta <- profiles_meta %>% filter(!is.na(nc_real) & !is.na(nc_uniform))
  }
  if (nrow(profiles_meta) == 0) {
    warning("[plot_vertical_tmax_profiles] Aucun NetCDF valide — abandon.")
    return(invisible(NULL))
  }

  df_profiles <- purrr::pmap_dfr(
    profiles_meta,
    function(Cluster, x, y, LAI, Hmax, fCover, plot_id, nc_real, nc_uniform, ...) {
      prof_real <- extract_vertical_tmax_profile(nc_real,    date_seq, summary_mode)
      prof_unif <- extract_vertical_tmax_profile(nc_uniform, date_seq, summary_mode)
      bind_rows(
        if (!is.null(prof_real))
          prof_real %>% mutate(scenario = "Real LAD"),
        if (!is.null(prof_unif))
          prof_unif %>% mutate(scenario = "Uniform LAD")
      ) %>%
        mutate(Cluster = as.character(Cluster),
               plot_id = plot_id,
               LAI     = LAI,
               Hmax    = Hmax)
    }
  )

  if (nrow(df_profiles) == 0) {
    warning("[plot_vertical_tmax_profiles] Aucun profil extrait.")
    return(invisible(NULL))
  }

  df_profiles <- df_profiles %>%
    mutate(panel_label = sprintf("Cluster %s\nHmax=%.1fm | LAI=%.1f (plot repr.)",
                                  Cluster, Hmax, LAI))

  df_hmax_lines <- df_profiles %>%
    group_by(panel_label) %>%
    summarise(Hmax = first(Hmax), .groups = "drop")

  p <- ggplot(df_profiles, aes(x = Tmax_z, y = height_m, colour = scenario)) +
    geom_path(linewidth = 1.2) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(data      = df_hmax_lines,
               aes(yintercept = Hmax),
               linetype  = "dashed", colour = "firebrick",
               linewidth = 0.8, inherit.aes = FALSE) +
    facet_wrap(~ panel_label, nrow = 1, scales = "free_x") +
    scale_colour_manual(
      values = c("Real LAD" = "#31688e", "Uniform LAD" = "#d8576b")
    ) +
    labs(
      title    = "Vertical Tmax profiles — Real LAD vs Uniform LAD",
      subtitle = sprintf(
        "Median plot per cluster  |  Mode: %s", summary_mode
      ),
      x      = expression(T[max] ~ "simulated" ~ (degree*C)),
      y      = "Height (m)",
      colour = "Scenario"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold"))

  png_path <- file.path(out_dir,
    sprintf("h2_vertical_tmax_profiles_%s.png", summary_mode))
  csv_path <- file.path(out_dir,
    sprintf("h2_vertical_tmax_profiles_%s.csv", summary_mode))

  save_plot(p, png_path, width = 14, height = 6)
  write.csv(df_profiles, csv_path, row.names = FALSE)
  cat(sprintf("  [vertical profile] PNG : %s\n", png_path))
  cat(sprintf("  [vertical profile] CSV : %s\n", csv_path))

  if (requireNamespace("patchwork", quietly = TRUE) &&
      length(unique(df_profiles$scenario)) == 2) {
    df_diff <- df_profiles %>%
      tidyr::pivot_wider(
        id_cols     = c(Cluster, nair, height_m, panel_label),
        names_from  = scenario,
        values_from = Tmax_z
      ) %>%
      dplyr::filter(!is.na(`Real LAD`) & !is.na(`Uniform LAD`)) %>%
      dplyr::mutate(Diff = `Real LAD` - `Uniform LAD`)

    if (nrow(df_diff) > 0) {
      p_diff <- ggplot(df_diff, aes(x = Diff, y = height_m)) +
        geom_path(linewidth = 1.0, colour = "#5ec962") +
        geom_vline(xintercept = 0, linetype = "dotted",
                   colour = "grey50", linewidth = 0.7) +
        geom_hline(data      = df_hmax_lines,
                   aes(yintercept = Hmax),
                   linetype  = "dashed", colour = "firebrick",
                   linewidth = 0.8, inherit.aes = FALSE) +
        facet_wrap(~ panel_label, nrow = 1, scales = "free_x") +
        labs(
          x = expression(Delta * T[max] ~ "(Real - Uniform)" ~ (degree*C)),
          y = "Height (m)"
        ) +
        theme_bw(base_size = 11) +
        theme(strip.text = element_blank())

      p_combined <- p / p_diff + patchwork::plot_layout(heights = c(2, 1))
      diff_path  <- file.path(out_dir,
        sprintf("h2_vertical_tmax_profiles_%s_with_diff.png", summary_mode))
      save_plot(p_combined, diff_path, width = 14, height = 8)
      cat(sprintf("  [vertical profile] PNG avec diff : %s\n", diff_path))
    }
  }

  invisible(list(plot = p, data = df_profiles, medians = medians))
}

#' Extraire le profil vertical de ΔTmax selon deux méthodes de définition du Tmax.
#'
#' Méthode 1 — indépendant : Tmax(z) = max_t(Tair_z(t)) par couche.
#' Méthode 2 — heure globale : profil à l'heure où la couche supérieure atteint son Tmax.
#' Retourne la moyenne sur toutes les dates de date_seq.
#'
#' @param nc_file  Chemin NetCDF.
#' @param df_macro Dataframe avec colonnes date et Tmax_macro (ERA5).
#' @param date_seq Dates de la période d'intérêt.
#' @return Tibble (height_m, nair, Delta_Tmax_mean, method) ou NULL.
extract_vertical_delta_both_methods <- function(nc_file, df_macro, date_seq) {
  if (is.na(nc_file) || !file.exists(nc_file)) return(NULL)

  nc <- try(nc_open(nc_file), silent = TRUE)
  if (inherits(nc, "try-error")) return(NULL)
  rel_h <- try(ncvar_get(nc, "relative_height"), silent = TRUE)
  h_top <- try(ncvar_get(nc, "veget_height_top"), silent = TRUE)
  raw   <- try(get_variable(nc, "Tair_z"), silent = TRUE)
  nc_close(nc)
  if (inherits(raw, "try-error") || is.null(raw)) return(NULL)

  hmax_nc <- if (!inherits(h_top, "try-error")) h_top[1] else NA_real_

  df_hourly <- raw %>%
    mutate(Tair_C = Tair_z - 273.15,
           time   = time - lubridate::hours(2),
           date   = as.Date(time)) %>%
    filter(date %in% date_seq) %>%
    inner_join(df_macro %>% dplyr::select(date, Tmax_macro), by = "date")

  top_nair <- max(df_hourly$nair, na.rm = TRUE)

  m1 <- df_hourly %>%
    group_by(nair, date, Tmax_macro) %>%
    summarise(Tmax_z = max(Tair_C, na.rm = TRUE), .groups = "drop") %>%
    mutate(Delta = Tmax_z - Tmax_macro,
           method = "Layer-independent")

  peak_times <- df_hourly %>%
    filter(nair == top_nair) %>%
    group_by(date) %>%
    slice_max(Tair_C, n = 1, with_ties = FALSE) %>%
    dplyr::select(date, peak_time = time)

  m2 <- df_hourly %>%
    inner_join(peak_times, by = c("date", "time" = "peak_time")) %>%
    group_by(nair, date, Tmax_macro) %>%
    summarise(Tmax_z = first(Tair_C), .groups = "drop") %>%
    mutate(Delta = Tmax_z - Tmax_macro,
           method = "Global Tmax time")

  df_out <- bind_rows(m1, m2) %>%
    group_by(nair, method) %>%
    summarise(Delta_Tmax_mean = mean(Delta, na.rm = TRUE), .groups = "drop")

  if (!inherits(rel_h, "try-error") && !is.na(hmax_nc)) {
    df_out$height_m <- rel_h[df_out$nair] * hmax_nc
  } else {
    n_lay <- max(df_out$nair)
    df_out$height_m <- seq(0, hmax_nc, length.out = n_lay)[df_out$nair]
  }

  df_out
}

#' Profils verticaux de ΔTmax par cluster — comparaison des deux méthodes Tmax.
#'
#' @param df_sample    cLHS avec Cluster.
#' @param real_lad_dir Dossier NetCDF du scénario full-real (H1f_4_Full_real ou REF).
#' @param df_macro     Dataframe ERA5 avec date et Tmax_macro.
#' @param date_seq     Dates de la période.
#' @param out_path     Chemin PNG de sortie.
#' @return Invisible list(plot, data).
plot_vertical_delta_by_cluster <- function(df_sample,
                                            real_lad_dir,
                                            df_macro,
                                            date_seq,
                                            out_path = "outputs/figures/06_h2_vertical_delta_by_cluster.png") {
  if (!"Cluster" %in% names(df_sample)) {
    warning("[plot_vertical_delta_by_cluster] colonne 'Cluster' absente — abandon.")
    return(invisible(NULL))
  }

  medians <- select_cluster_median_plots(df_sample)
  if (is.null(medians) || nrow(medians) == 0) return(invisible(NULL))

  find_nc <- function(dir, x, y) {
    pat  <- sprintf("X%d_Y%d", round(x), round(y))
    hits <- list.files(dir, pattern = pat, full.names = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  df_profiles <- medians %>%
    rowwise() %>%
    mutate(nc = find_nc(real_lad_dir, x, y)) %>%
    ungroup() %>%
    filter(!is.na(nc)) %>%
    purrr::pmap_dfr(function(Cluster, x, y, LAI, Hmax, plot_id, nc, ...) {
      prof <- extract_vertical_delta_both_methods(nc, df_macro, date_seq)
      if (is.null(prof)) return(NULL)
      prof %>% mutate(Cluster  = as.character(Cluster),
                      plot_id  = plot_id,
                      Hmax_plot = Hmax,
                      LAI_plot  = LAI)
    })

  if (is.null(df_profiles) || nrow(df_profiles) == 0) {
    warning("[plot_vertical_delta_by_cluster] aucun profil extrait — abandon.")
    return(invisible(NULL))
  }

  # Clip atmospheric layers above 1.5 × canopy height to avoid confounding
  # boundary-layer effects (very hot air above sparse, short canopies) with
  # within-canopy microclimate.
  df_profiles <- df_profiles %>%
    filter(height_m <= Hmax_plot * 1.5) %>%
    mutate(
      method      = factor(method, levels = c("Layer-independent",
                                              "Global Tmax time")),
      panel_label = sprintf("Cluster %s\nHmax=%.1fm | LAI=%.1f (plot repr.)",
                            Cluster, Hmax_plot, LAI_plot)
    )

  df_hmax <- df_profiles %>%
    group_by(panel_label) %>%
    summarise(Hmax = first(Hmax_plot), .groups = "drop")

  p <- ggplot(df_profiles,
              aes(x = Delta_Tmax_mean, y = height_m,
                  colour = method, linetype = method)) +
    geom_path(linewidth = 1.1) +
    geom_point(size = 1.8) +
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey50", linewidth = 0.6) +
    geom_hline(data = df_hmax, aes(yintercept = Hmax),
               linetype = "dashed", colour = "firebrick",
               linewidth = 0.8, inherit.aes = FALSE) +
    scale_colour_manual(
      values   = c("Layer-independent" = "#31688e",
                   "Global Tmax time"  = "#d8576b"),
      name     = "Tmax method"
    ) +
    scale_linetype_manual(
      values   = c("Layer-independent" = "solid",
                   "Global Tmax time"  = "dashed"),
      name     = "Tmax method"
    ) +
    facet_wrap(~ panel_label, nrow = 1, scales = "free_x") +
    labs(
      title    = "Vertical ΔTmax profiles by cluster — comparison of two methods",
      subtitle = "Blue = layer-independent Tmax | Red = profile at top-layer Tmax time\nRed horizontal dash = canopy height | ΔTmax = Tmax(layer) - Tmax ERA5",
      x        = expression(mean ~ Delta * T[max] ~ "summer mean" ~ (degree*C)),
      y        = "Height (m)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold"),
          plot.title      = element_text(face = "bold"))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 14, height = 7)
  csv_path <- sub("\\.png$", ".csv", out_path)
  write.csv(df_profiles, csv_path, row.names = FALSE)
  cat(sprintf("  [delta_cluster_profiles] PNG : %s\n", out_path))
  invisible(list(plot = p, data = df_profiles))
}

#' Profils verticaux Tmax multi-scénarios — figure triptyque comparatif.
#'
#' @param scenarios_dict Liste nommée : list("Etiquette" = "chemin/NetCDF", ...).
#' @param df_sample   cLHS (Cluster, x, y, LAI, Hmax, fCover).
#' @param date_seq    Dates d'intérêt.
#' @param summary_mode "mean_summer" | "canicule_date" | "median_date".
#' @param out_path    Chemin PNG de sortie.
#' @param palette     Vecteur de couleurs (longueur = nb scénarios).
#' @return Invisible ggplot.
plot_vertical_tmax_profiles_multi <- function(scenarios_dict, df_sample,
                                               date_seq,
                                               summary_mode = "canicule_date",
                                               out_path     = "outputs/figures/07_h2_vertical_3way.png",
                                               palette      = NULL,
                                               diff_from    = NULL) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  medians <- select_cluster_median_plots(df_sample)
  if (is.null(medians) || nrow(medians) == 0) {
    warning("[plot_vertical_tmax_profiles_multi] Aucun plot médian — abandon.")
    return(invisible(NULL))
  }

  find_nc <- function(dir, x, y) {
    pat  <- sprintf("X%d_Y%d", round(x), round(y))
    hits <- list.files(dir, pattern = pat, full.names = TRUE)
    if (length(hits) == 0) NA_character_ else hits[1]
  }

  df_all <- purrr::imap_dfr(scenarios_dict, function(nc_dir, sc_label) {
    purrr::pmap_dfr(medians, function(Cluster, x, y, LAI, Hmax, fCover, plot_id, ...) {
      nc_file <- find_nc(nc_dir, x, y)
      prof    <- extract_vertical_tmax_profile(nc_file, date_seq, summary_mode)
      if (is.null(prof)) return(NULL)
      prof %>% mutate(scenario  = sc_label,
                      Cluster   = as.character(Cluster),
                      plot_id   = plot_id,
                      LAI       = LAI,
                      Hmax      = Hmax)
    })
  })

  if (nrow(df_all) == 0) {
    warning("[plot_vertical_tmax_profiles_multi] Aucun profil extrait.")
    return(invisible(NULL))
  }

  df_all <- df_all %>%
    mutate(panel_label = sprintf("Cluster %s\nHmax=%.1fm | LAI=%.1f (plot repr.)",
                                  Cluster, Hmax, LAI))

  df_hmax <- df_all %>%
    group_by(panel_label) %>%
    summarise(Hmax = first(Hmax), .groups = "drop")

  sc_names <- names(scenarios_dict)
  if (is.null(palette)) {
    palette <- setNames(
      c("#31688e", "#d8576b", "#5ec962", "#fca50a", "#440154")[seq_along(sc_names)],
      sc_names
    )
  }

  annot_layer <- NULL

  # Differential mode: show T(scenario) - T(diff_from) instead of absolute T.
  # Makes near-zero effects visible when 3 lines would otherwise overlap.
  if (!is.null(diff_from) && diff_from %in% sc_names) {
    df_base <- df_all %>%
      filter(scenario == diff_from) %>%
      dplyr::select(panel_label, nair, Tmax_base = Tmax_z)
    df_all <- df_all %>%
      inner_join(df_base, by = c("panel_label", "nair")) %>%
      mutate(Tmax_z = Tmax_z - Tmax_base) %>%
      filter(scenario != diff_from)
    sc_names <- setdiff(sc_names, diff_from)
    palette  <- palette[sc_names]
    x_lab    <- bquote(Delta * T[max] ~ "vs" ~ .(diff_from) ~ (degree * C))
    vline    <- geom_vline(xintercept = 0, linetype = "dashed",
                           colour = "grey30", linewidth = 0.6)
    hline_obj <- geom_hline(data = df_hmax,
                             aes(yintercept = Hmax),
                             linetype = "dashed", colour = "firebrick",
                             linewidth = 0.8, inherit.aes = FALSE)
    title_str <- sprintf("Vertical ΔT profiles vs %s — %d scenarios", diff_from, length(sc_names))

    df_annot <- df_all %>%
      group_by(panel_label) %>%
      summarise(max_abs = max(abs(Tmax_z), na.rm = TRUE),
                y_pos   = max(height_m,    na.rm = TRUE) * 0.7,
                .groups = "drop") %>%
      mutate(label = sprintf("|ΔT|max = %.4f°C", max_abs))
    annot_layer <- geom_label(
      data        = df_annot,
      aes(y = y_pos, label = label),
      x           = -Inf,
      inherit.aes = FALSE,
      hjust       = -0.05, vjust = 0.5,
      size        = 3.2, fill = "lightyellow", colour = "grey20",
      label.size  = 0.4, fontface = "bold"
    )
  } else {
    x_lab     <- expression(T[max] ~ simulee ~ (degree * C))
    vline     <- NULL
    hline_obj <- geom_hline(data = df_hmax,
                             aes(yintercept = Hmax),
                             linetype = "dashed", colour = "firebrick",
                             linewidth = 0.8, inherit.aes = FALSE)
    title_str <- sprintf("Vertical Tmax profiles — comparison of %d scenarios",
                         length(sc_names))
  }

  p <- ggplot(df_all, aes(x = Tmax_z, y = height_m,
                            colour = scenario, linetype = scenario)) +
    geom_path(linewidth = 1.4) +
    geom_point(size = 2, alpha = 0.8) +
    vline +
    hline_obj +
    annot_layer +
    facet_wrap(~ panel_label, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = palette) +
    scale_linetype_manual(
      values = setNames(
        c("solid","dashed","dotted","dotdash","longdash")[seq_along(sc_names)],
        sc_names
      )
    ) +
    labs(
      title    = title_str,
      subtitle = sprintf("Median plot per cluster  |  Mode: %s", summary_mode),
      x        = x_lab,
      y        = "Height (m)",
      colour   = "Scenario", linetype = "Scenario"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          strip.text      = element_text(face = "bold", size = 12))

  save_plot(p, out_path, width = 16, height = 8)
  cat(sprintf("  [3-way profile] PNG : %s\n", out_path))
  invisible(p)
}
