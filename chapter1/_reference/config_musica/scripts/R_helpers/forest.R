# ==============================================================================
# Chapter 1 — Forest dataframe construction, clustering, and cLHS sampling
# ==============================================================================

#' Build the forest dataframe from a stacked LiDAR raster.
#'
#' Filters out pixels where the LAD profile sum is below a minimum threshold
#' and returns the pixel coordinates, scalar metrics, and LAD matrix.
#'
#' @param r_stack Terra SpatRaster stack (LAI, VCI, Hmax, fCover, LAD_Layer_*).
#' @return List with elements $df (dataframe), $mat_lad (matrix), $z_breaks (vector).
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

#' Find the optimal number of k-means clusters via the elbow method.
#'
#' Uses maximum perpendicular distance from the WSS curve to the line
#' connecting its endpoints.
#'
#' @param df_scaled Scaled numeric matrix or dataframe.
#' @param k_range   Integer range to test (default 1:10).
#' @return List with elements $opt_k, $wss, $k_range, $plot.
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

#' Assign k-means cluster labels to the forest dataframe.
#'
#' If k is NULL the optimal number is determined automatically via
#' optimise_k_elbow().
#'
#' @param df_forest Forest dataframe from build_forest_dataframe.
#' @param k         Fixed number of clusters, or NULL for automatic selection.
#' @param vars      Variables used for clustering (default LAI, Hmax, fCover).
#' @param max_k     Upper bound for automatic k selection (default 10).
#' @return List with elements $df (df_forest with Cluster column) and $elbow_plot.
label_clusters <- function(df_forest, k = NULL, vars = c("LAI", "Hmax", "fCover"), max_k = 10) {
  set.seed(42)
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

#' Plot mean LAD profile per cluster in relative height space.
#'
#' Profil LAD moyen par cluster — base d'interprétation typologique.
#' Permet de nommer visuellement les clusters (bottom-heavy, top-heavy, etc.)
#' en affichant le profil moyen de chaque groupe dans l'espace de hauteur absolu.
#'
#' @param df_forest Dataframe avec colonnes x, y, Cluster (issu de label_clusters).
#' @param mat_rel   Matrice LAD normalisée en hauteur relative (de compute_fpca).
#' @param z_rel     Vecteur de hauteurs relatives [0,1].
#' @return ggplot object.
plot_cluster_mean_profiles <- function(df_forest, mat_rel, z_rel) {
  valid    <- !is.na(df_forest$Cluster)
  clusters <- if (is.factor(df_forest$Cluster)) levels(droplevels(df_forest$Cluster[valid]))
              else sort(unique(as.character(df_forest$Cluster[valid])))

  df_profiles <- lapply(clusters, function(cl) {
    idx  <- which(df_forest$Cluster == cl & valid)
    mat  <- mat_rel[idx, , drop = FALSE]
    q25  <- apply(mat, 2, quantile, 0.25, na.rm = TRUE)
    q75  <- apply(mat, 2, quantile, 0.75, na.rm = TRUE)
    data.frame(
      rel_z   = z_rel,
      density = colMeans(mat, na.rm = TRUE),
      q25     = q25,
      q75     = q75,
      Cluster = cl,
      n       = length(idx),
      Hmax_mean = round(mean(df_forest$Hmax[idx], na.rm = TRUE), 1),
      LAI_mean  = round(mean(df_forest$LAI[idx],  na.rm = TRUE), 2)
    )
  })
  df_profiles <- bind_rows(df_profiles) %>%
    mutate(label = sprintf("Cluster %s  (n=%d)\nHmax moy.=%.1fm | LAI moy.=%.2f",
                           Cluster, n, Hmax_mean, LAI_mean))

  ggplot(df_profiles, aes(x = density, y = rel_z, colour = Cluster, fill = Cluster)) +
    geom_ribbon(aes(xmin = q25, xmax = q75), alpha = 0.15, colour = NA) +
    geom_path(linewidth = 1.3) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "firebrick",
               linewidth = 0.7) +
    annotate("text", x = Inf, y = 1, label = "mean Hmax",
             hjust = 1.05, vjust = -0.4, colour = "firebrick", size = 3) +
    facet_wrap(~ label, nrow = 1) +
    scale_colour_viridis_d(option = "turbo", end = 0.85) +
    scale_fill_viridis_d(option = "turbo",   end = 0.85) +
    labs(
      title    = "Mean LAD profile per cluster — typological interpretation basis",
      subtitle = "Relative height (z/Hmax) | Ribbon = IQR | Dash = z/Hmax = 1 (canopy top) | Hmax and LAI = cluster mean",
      x        = expression("Mean LAD" ~ (m^2 ~ m^{-3})),
      y        = "Relative height (z / Hmax)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          strip.text      = element_text(face = "bold"))
}

#' Plot ALL individual LAD profiles per cluster (with mean overlay).
#'
#' Vue d'inspection : superpose tous les profils du sample dans chaque cluster
#' (alpha = 0.15) pour voir la dispersion intra-cluster. Le profil moyen est
#' surimprime en epais pour reference.
#'
#' Complementaire de `plot_cluster_mean_profiles()` qui ne montre que la moyenne :
#' ici on voit aussi la variabilite (typique pour reperer les clusters "homogenes"
#' au sens VCI vs ceux qui melangent plusieurs sous-types).
#'
#' @param df  Dataframe avec colonnes x, y, Cluster, Hmax, LAI et LAD_Layer_* .
#'   Typiquement df_sample (cLHS, 100 plots/cluster) — plus lisible que df_forest.
#' @return ggplot object.
plot_cluster_all_profiles <- function(df) {
  lad_cols <- grep("^LAD_Layer_", names(df), value = TRUE)
  if (length(lad_cols) == 0) stop("Pas de colonnes LAD_Layer_* dans le dataframe.")

  df$plot_id <- seq_len(nrow(df))
  df_long <- df %>%
    dplyr::select(plot_id, Cluster, Hmax, LAI, all_of(lad_cols)) %>%
    tidyr::pivot_longer(cols = all_of(lad_cols),
                         names_to = "h_lab", values_to = "density") %>%
    mutate(height = as.numeric(gsub("LAD_Layer_", "", h_lab)),
           Cluster = factor(Cluster))

  df_mean <- df_long %>%
    group_by(Cluster, height) %>%
    summarise(density_mean = mean(density, na.rm = TRUE), .groups = "drop")

  df_stats <- df %>%
    group_by(Cluster) %>%
    summarise(n = dplyr::n(),
              Hmax_med = round(median(Hmax, na.rm = TRUE), 0),
              LAI_med  = round(median(LAI,  na.rm = TRUE), 2),
              .groups = "drop") %>%
    mutate(label = sprintf("Cluster %s  (n=%d)\nHmax_med=%dm | LAI_med=%.2f",
                            Cluster, n, Hmax_med, LAI_med))

  pal <- c("1" = "#440154", "2" = "#31688e", "3" = "#35b779", "4" = "#fde725")
  max_h <- max(df_long$height[df_long$density > 0], na.rm = TRUE) + 2
  max_d <- max(df_long$density, na.rm = TRUE) * 1.02

  ggplot(df_long, aes(x = density, y = height)) +
    geom_path(aes(group = plot_id, colour = Cluster),
              alpha = 0.15, linewidth = 0.35) +
    geom_path(data = df_mean,
              aes(x = density_mean, y = height),
              colour = "black", linewidth = 1.1) +
    geom_path(data = df_mean,
              aes(x = density_mean, y = height, colour = Cluster),
              linewidth = 0.85) +
    facet_wrap(~ Cluster, nrow = 1,
               labeller = labeller(Cluster = function(cl) {
                 sapply(cl, function(c) df_stats$label[df_stats$Cluster == c])
               })) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(0, max_d),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(breaks = seq(0, max_h, 5),
                       limits = c(0, max_h),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "All vertical LAD profiles per cluster — sample inspection",
      subtitle = "Thin lines: individual plots (alpha = 0.15) | Thick black + colour: cluster mean",
      caption  = "Heights in absolute meters. Cluster 1 = young/sparse; Cluster 2 = mature/dense bi-strata; Cluster 3 = low closed; Cluster 4 = intermediate.",
      x = expression("LAD density" ~ (m^2 ~ m^{-3})),
      y = "Height above ground (m)"
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text  = element_text(face = "bold", size = 9.5),
          plot.title  = element_text(face = "bold"),
          plot.caption= element_text(size = 8, colour = "grey35", hjust = 0,
                                      lineheight = 1.2))
}

#' Conditioned Latin Hypercube Sampling (cLHS) stratified per cluster.
#'
#' @param df_forest    Forest dataframe with Cluster column.
#' @param n_per_cluster Target sample size per cluster (default 100).
#' @param vars         Variables to include in the LHS design.
#' @param iter         Number of cLHS iterations (default 10000).
#' @return Dataframe of sampled plots.
sample_clhs_per_cluster <- function(df_forest, n_per_cluster = 100, vars = c("LAI", "Hmax", "fCover"), iter = 10000) {
  set.seed(42)
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
