# ==============================================================================
# Chapter 1 — LAD profile constructors
# Chapter 3 — d_opt-truncated LAD constructor
# ==============================================================================

#' Rescale a LAD profile to a new height and leaf area, PRESERVING ITS SHAPE
#'
#' The shape-preserving rescaler behind the whole perturbation design. Height is changed
#' by stretching the profile in RELATIVE height (z / hmax) before re-interpolating onto
#' the 1 m grid; leaf area is changed by one scalar multiplying every layer. Neither
#' operation redistributes foliage vertically, which is what lets the real-versus-uniform
#' contrast be the ONLY perturbation in the chapter that changes shape (Section 2.6).
#'
#' @param z_real Heights of the measured profile, m.
#' @param d Leaf-area density at those heights.
#' @param hmax_ref Canopy height the measured profile belongs to, m.
#' @param hmax_tgt Target canopy height, m.
#' @param lai_tgt Target total leaf area, in the units `d` sums to (two-sided for MuSICA).
#' @return data.frame(height, density) on a 1 m grid from 1 m to ceiling(hmax_tgt).
.lad_rescale <- function(z_real, d, hmax_ref, hmax_tgt, lai_tgt) {
  z_rel  <- if (hmax_ref > 0) z_real / hmax_ref else z_real
  z_new  <- z_rel * hmax_tgt
  z_grid <- seq(1L, ceiling(hmax_tgt))
  d_new  <- approx(z_new, d, z_grid, rule = 2)$y
  d_new[is.na(d_new)] <- 0
  if (sum(d_new) > 1e-9) d_new <- d_new * (lai_tgt / sum(d_new))
  data.frame(height = z_grid, density = d_new)
}

#' Build a real LAD profile from a plot row (with optional Hmax/LAI override).
#'
#' @param plot_row Single-row dataframe with LAD_Layer_* columns and Hmax.
#' @param hmax     Optional Hmax override (m).
#' @param lai      Optional LAI override.
#' @return Dataframe with columns height (m) and density (m² m⁻³).
make_lad_real <- function(plot_row, hmax = NULL, lai = NULL) {
  lad_cols  <- grep("LAD_Layer_", names(plot_row), value = TRUE)
  z_real    <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  d         <- as.numeric(plot_row[lad_cols]); d[is.na(d)] <- 0
  hmax_real <- as.numeric(plot_row$Hmax)
  if (hmax_real <= 0)
    message(sprintf("[make_lad_real] hmax_real=0 for plot (x=%s,y=%s) — profile will be all-zero",
                    plot_row$x, plot_row$y))
  if (is.null(hmax)) hmax <- hmax_real
  if (is.null(lai))  lai  <- sum(d)
  .lad_rescale(z_real, d, hmax_real, hmax, lai)
}

#' Build a uniform LAD profile (constant density from canopy_base to Hmax).
#'
#' @param plot_row    Single-row dataframe with Hmax and LAI columns.
#' @param hmax        Optional Hmax override (m).
#' @param lai         Optional LAI override.
#' @param canopy_base Height of canopy base in metres (default 1).
#' @return Dataframe with columns height (m) and density (m² m⁻³).
make_lad_uniform <- function(plot_row, hmax = NULL, lai = NULL, canopy_base = 1) {
  if (is.null(hmax)) hmax <- as.numeric(plot_row$Hmax)
  if (is.null(lai))  lai  <- as.numeric(plot_row$LAI)
  max_layer <- ceiling(hmax)
  if (max_layer <= canopy_base) return(data.frame(height = 1L, density = 0))
  depth   <- max_layer - canopy_base + 1L
  density <- lai / depth
  heights <- seq_len(max_layer)
  data.frame(height = heights, density = ifelse(heights >= canopy_base, density, 0))
}

#' Factory: build a mean-shape LAD constructor from the full sample.
#'
#' Returns a function with the same signature as make_lad_real/uniform that
#' produces the sample-mean profile rescaled to each plot's (Hmax, LAI).
#'
#' @param df_sample cLHS sample dataframe with LAD_Layer_* columns.
#' @return Function(plot_row, hmax, lai) → data.frame(height, density).
make_lad_mean_factory <- function(df_sample) {
  lad_cols      <- grep("LAD_Layer_", names(df_sample), value = TRUE)
  z_real        <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  mat           <- as.matrix(df_sample[, lad_cols]); mat[is.na(mat)] <- 0
  mean_prof     <- colMeans(mat, na.rm = TRUE)
  hmax_mean_ref <- mean(df_sample$Hmax, na.rm = TRUE)

  function(plot_row, hmax = NULL, lai = NULL) {
    if (is.null(hmax)) hmax <- as.numeric(plot_row$Hmax)
    if (is.null(lai))  lai  <- as.numeric(plot_row$LAI)
    .lad_rescale(z_real, mean_prof, hmax_mean_ref, hmax, lai)
  }
}

#' Factory: build a cluster-type LAD constructor from the full sample.
#'
#' Returns a function that substitutes each plot's real shape with the mean
#' profile of its cluster, rescaled to its own (Hmax, LAI).
#'
#' @param df_sample cLHS sample with LAD_Layer_* and Cluster columns.
#' @return Function(plot_row, hmax, lai) → data.frame(height, density).
make_lad_cluster_type_factory <- function(df_sample) {
  lad_cols <- grep("LAD_Layer_", names(df_sample), value = TRUE)
  z_real   <- as.numeric(gsub("LAD_Layer_", "", lad_cols))

  cluster_profiles <- df_sample %>%
    filter(!is.na(Cluster)) %>%
    split(.$Cluster) %>%
    purrr::map(function(df_sub) {
      mat <- as.matrix(df_sub[, lad_cols]); mat[is.na(mat)] <- 0
      list(
        profile       = colMeans(mat, na.rm = TRUE),
        hmax_mean_ref = mean(df_sub$Hmax, na.rm = TRUE)
      )
    })

  function(plot_row, hmax = NULL, lai = NULL) {
    cl <- as.character(plot_row$Cluster)
    if (is.na(cl) || !(cl %in% names(cluster_profiles)))
      return(data.frame(height = z_real, density = rep(0, length(z_real))))
    cp <- cluster_profiles[[cl]]
    if (is.null(hmax)) hmax <- as.numeric(plot_row$Hmax)
    if (is.null(lai))  lai  <- as.numeric(plot_row$LAI)
    .lad_rescale(z_real, cp$profile, cp$hmax_mean_ref, hmax, lai)
  }
}

# ==============================================================================
# Chapter 3 — d_opt-truncated LAD constructor
# ==============================================================================

#' Build a real LAD profile truncated to d_opt height layers from the canopy top.
#'
#' @title d_opt-truncated LAD profile
#' @description Implements the physical interpretation of d_opt: Sentinel-2
#'   reflectance is sensitive only to the top `d_opt` metres of the canopy.
#'   This function zeroes out all LAD layers below `Hmax - d_opt` metres,
#'   retaining only the canopy layers that S2 actually sees, then rescales to
#'   `lai_tgt` so the total LAI matches the chosen LAI source.
#'
#'   Use this when pairing with `LAI_S2_DOPT` so that both the LAI amplitude
#'   (S2-corrected for depth) and the LAD shape (depth-truncated) are physically
#'   consistent with the effective canopy depth seen by S2.
#'
#' @param plot_row  Single-row dataframe with `LAD_Layer_*` columns and `Hmax`.
#' @param hmax      Optional Hmax override (m).
#' @param lai       Optional LAI override (applied to truncated layers only).
#' @param d_opt     Effective canopy depth seen by S2 in metres (default 6).
#' @return Dataframe with columns `height` (m) and `density` (m² m⁻³).
#'   Layers below `Hmax - d_opt` have density = 0.
#' @seealso [make_lad_real()]
make_lad_dopt <- function(plot_row, hmax = NULL, lai = NULL, d_opt = 6) {
  lad_cols  <- grep("LAD_Layer_", names(plot_row), value = TRUE)
  z_real    <- as.numeric(gsub("LAD_Layer_", "", lad_cols))
  d         <- as.numeric(plot_row[lad_cols]); d[is.na(d)] <- 0
  hmax_real <- as.numeric(plot_row$Hmax)
  if (is.null(hmax)) hmax <- hmax_real
  if (is.null(lai))  lai  <- sum(d)

  # First rescale to target hmax/lai (same as make_lad_real)
  lad_full  <- .lad_rescale(z_real, d, hmax_real, hmax, lai)

  # Then zero out layers below (Hmax - d_opt) — keep only top d_opt metres
  h_thresh  <- hmax - d_opt
  lad_full$density[lad_full$height <= h_thresh] <- 0

  # Re-normalise so the remaining layers integrate to lai_tgt
  tot <- sum(lad_full$density)
  if (tot > 1e-9) lad_full$density <- lad_full$density * (lai / tot)

  lad_full
}

#' Factory: build a make_lad_dopt closure with a fixed d_opt value.
#'
#' @title d_opt LAD factory
#' @description Returns a function with the same signature as `make_lad_real()`
#'   but using a fixed `d_opt` depth truncation. Useful for embedding d_opt into
#'   a scenario list without repeating the argument.
#'
#' @param d_opt  Effective canopy depth (m). Default 6 (Pareto Blois).
#' @return Function `f(plot_row, hmax, lai) → data.frame(height, density)`.
make_lad_dopt_factory <- function(d_opt = 6) {
  force(d_opt)
  function(plot_row, hmax = NULL, lai = NULL) {
    make_lad_dopt(plot_row, hmax = hmax, lai = lai, d_opt = d_opt)
  }
}
