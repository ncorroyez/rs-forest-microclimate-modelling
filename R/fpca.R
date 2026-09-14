# ==============================================================================
# Chapter 1 — Shape metrics: FPCA and Height of Median LAI
# ==============================================================================

#' Sélection automatique du nombre de FPCs à inclure dans le clustering / cLHS.
#'
#' Critère : gain marginal de variance expliquée par FPC_k. On inclut FPC_k
#' tant que varprop[k] > min_marginal.
#'
#' @param varprop      Vecteur de proportions de variance (issu de fpca$varprop).
#' @param min_marginal Seuil de gain marginal minimal (défaut 0.05 = 5 %).
#' @return Entier : nombre de FPCs à retenir (au moins 1).
select_n_fpc <- function(varprop, min_marginal = 0.05) {
  n <- sum(varprop >= min_marginal)
  max(1L, n)
}

#' Find the optimal number of B-spline basis functions via GCV elbow.
#'
#' @param arg_vals      Vector of argument values (e.g. relative heights).
#' @param mat_lad_norm  Normalised LAD matrix (rows = sites, cols = heights).
#' @param nbasis_range  Range of nbasis values to test (default 4:15).
#' @return List with elements $opt, $gcv, $range, $plot.
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

#' Compute FPCA on the LAD profiles after normalisation to relative height.
#'
#' Profiles are normalised by total LAI and interpolated to a common relative
#' height grid [0,1] before B-spline smoothing and PCA.
#'
#' @param mat_lad  LAD matrix (rows = sites, cols = height layers).
#' @param hmax_vec Vector of Hmax values aligned with mat_lad rows.
#' @param z_breaks Height breaks corresponding to the columns of mat_lad.
#' @param n_harm   Number of harmonics to retain (default 3).
#' @return List with elements $fpca, $fd_obj, $mat_rel, $z_rel,
#'   $nbasis_opt, $varprop, $gcv_plot.
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

#' Profile plots for each FPC: mean curve + Low 5% vs High 5% groups.
#'
#' @param fpca_res    Result of compute_fpca.
#' @param mat_lad_real Not used directly but kept for signature compatibility.
#' @param z_breaks    Height breaks (not used directly here).
#' @param hmax_vec    Hmax values (not used directly here).
#' @return patchwork ggplot of three harmonic panels.
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

#' Plot FPCA loading functions (contribution of each height to each FPC).
#'
#' @param fpca_res Result of compute_fpca.
#' @return ggplot object with facets per FPC.
plot_fpc_loadings <- function(fpca_res) {
  harmonics <- fpca_res$fpca$harmonics; z_rel <- fpca_res$z_rel; vals <- eval.fd(z_rel, harmonics)
  df <- as.data.frame(vals); names(df) <- paste0("FPC", seq_len(ncol(df))); df$rel_z <- z_rel
  df_long <- pivot_longer(df, -rel_z, names_to = "FPC", values_to = "loading")
  ggplot(df_long, aes(x = loading, y = rel_z, colour = FPC)) + geom_path(linewidth = 1.1) + geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~ FPC, nrow = 1) + scale_colour_viridis_d(option = "turbo", end = 0.85) + labs(title = "Per-height contribution of each FPC (FPCA loadings)", x = "Loading", y = "Relative height (Z/Hmax)") + theme(legend.position = "none")
}

#' Scatter plots of FPC score pairs, optionally coloured by cluster.
#'
#' @param fpca_res    Result of compute_fpca.
#' @param cluster_vec Optional factor vector of cluster labels.
#' @return patchwork ggplot with FPC1×FPC2, FPC1×FPC3, FPC2×FPC3.
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

#' Qualité de reconstruction FPCA : profils réels vs 3 premières composantes.
#'
#' Pour chaque site sélectionné, superpose le profil LAD normalisé réel
#' et sa reconstruction : mean(z) + Σ_{k=1}^{3} score_k × FPC_k(z).
#'
#' @param fpca_res    Résultat de compute_fpca() (contient mat_rel, z_rel, fpca).
#' @param idx_sites   Indices explicites. Si NULL : n_sites sites espacés sur FPC1.
#' @param n_sites     Nombre de sites si idx_sites est NULL (défaut 6).
#' @param cluster_vec Vecteur optionnel de labels Cluster pour annoter les panneaux.
#' @param out_path    Chemin PNG de sortie.
#' @return Invisible ggplot.
plot_fpca_reconstruction <- function(fpca_res, idx_sites = NULL, n_sites = 6,
                                      cluster_vec = NULL,
                                      out_path = "outputs/figures/10_fpca_reconstruction.png") {
  z_rel       <- fpca_res$z_rel
  mat_real    <- fpca_res$mat_rel
  scores      <- fpca_res$fpca$scores
  mean_curve  <- as.numeric(eval.fd(z_rel, fpca_res$fpca$meanfd))
  harm_matrix <- eval.fd(z_rel, fpca_res$fpca$harmonics)

  n_plots <- nrow(mat_real)
  if (is.null(idx_sites)) {
    fpc1_order <- order(scores[, 1])
    idx_sites  <- fpc1_order[round(seq(1, n_plots, length.out = n_sites))]
  }
  idx_sites <- idx_sites[idx_sites >= 1 & idx_sites <= n_plots]

  df_all <- purrr::imap_dfr(idx_sites, function(i, ii) {
    recon <- pmax(as.numeric(mean_curve + harm_matrix[, 1:3] %*% scores[i, 1:3]), 0)
    cl_label   <- if (!is.null(cluster_vec)) sprintf(" · Cl.%s", cluster_vec[i]) else ""
    site_label <- sprintf("Site #%d%s\nFPC1=%.2f", i, cl_label, scores[i, 1])
    bind_rows(
      data.frame(rel_z = z_rel, lad = mat_real[i, ],
                 source = "Real (LiDAR)",        site = site_label),
      data.frame(rel_z = z_rel, lad = recon,
                 source = "Reconstruit (3 CP)", site = site_label)
    )
  })

  var_str <- paste(sprintf("FPC%d=%.1f%%", 1:3, 100 * fpca_res$varprop[1:3]),
                   collapse = " | ")

  p <- ggplot(df_all, aes(x = lad, y = rel_z, colour = source, linetype = source)) +
    geom_path(linewidth = 1.1) +
    facet_wrap(~ site, nrow = 2, scales = "free_x") +
    scale_colour_manual(
      values = c("Real (LiDAR)" = "#31688e", "Reconstruit (3 CP)" = "#d8576b")
    ) +
    scale_linetype_manual(
      values = c("Real (LiDAR)" = "solid", "Reconstruit (3 CP)" = "dashed")
    ) +
    labs(
      title    = "FPCA reconstruction — real profiles vs 3 first principal components",
      subtitle = var_str,
      x        = "Normalised LAD (relative)",
      y        = "Relative height (z / Hmax)",
      colour = NULL, linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_plot(p, out_path, width = 14, height = 7)
  cat(sprintf("  [FPCA reconstruction] %d sites — %s\n", length(idx_sites), out_path))
  invisible(p)
}

#' Compute the Height of Median LAI (centre of gravity of the LAD profile).
#'
#' @param mat_lad  LAD matrix (rows = sites, cols = height layers).
#' @param z_breaks Height break values aligned with mat_lad columns.
#' @return Numeric vector of H_median values (one per site).
compute_h_median <- function(mat_lad, z_breaks) {
  apply(mat_lad, 1, function(prof) {
    total_lai <- sum(prof, na.rm = TRUE)
    if (total_lai <= 0) return(NA_real_)
    cum_lai <- cumsum(prof)
    idx <- which(cum_lai >= (total_lai / 2))[1]
    return(z_breaks[idx])
  })
}
