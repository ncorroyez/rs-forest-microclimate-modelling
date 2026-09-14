# ==============================================================================
# Chapter 3 — Full scenario matrix (13 scenarios, 3 axes)
#
# Axis 1 — LAI source (6 variants)
#   ALS         : LiDAR full-depth LAI (reference)
#   ALS_DOPT    : LiDAR at d_opt layers from top (Pareto, Blois d_opt=6)
#   S2_ATBD     : S2 PROSAIL inversion, ATBD parameterisation
#   S2_DOPT     : S2 PROSAIL inversion, d_opt-optimised parameterisation
#   S2_RESCALED : S2_ATBD rescaled so peak matches ALS amplitude
#   RF          : Random-forest correction of S2_ATBD using structural covariates
#
# Axis 2 — Phenology mode (2 variants)
#   STATIC      : parametric calc_phenology() — same budburst/shape as Chapter 1
#   DYN         : S2-derived phenology (smoothed GAM, make_phenology_from_s2)
#
# Axis 3 — LAD shape (2 variants, only for DOPT scenarios)
#   LAD_FULL    : full-depth LiDAR LAD (make_lad_real)
#   LAD_DOPT    : top-d_opt-metres-only LiDAR LAD (make_lad_dopt)
#
# Naming convention: {MODE}_{LAI_SOURCE}[_LADOPT]
#   STATIC_ prefix = parametric phenology
#   DYN_    prefix = S2-driven phenology
#   _LADOPT suffix = depth-truncated LAD (only for DOPT scenarios)
#
# Scenario struct fields (backwards-compatible with Chapter 1):
#   $name         character  output directory name
#   $lai_fn       function(plot_row) → numeric
#   $hmax_fn      function(plot_row) → numeric
#   $fcover_fn    function(plot_row) → numeric
#   $lad_fn       function(plot_row, hmax, lai) → data.frame(height, density)
#   $phenology_fn function(plot_row) → data.frame | NULL  [NULL = parametric]
# ==============================================================================

# ---- Internal helpers --------------------------------------------------------

.fn_col <- function(col) function(plot_row) as.numeric(plot_row[[col]])

# Build a phenology_fn from a named list of per-plot time series
.make_pheno_fn <- function(ts_by_plot, list_year) {
  if (is.null(ts_by_plot)) return(NULL)
  make_phenology_fn_factory(ts_by_plot, list_year)
}

# Build one scenario list
.make_sc <- function(name, lai_col, lad_fn, pheno_fn = NULL) {
  list(
    name         = name,
    lai_fn       = .fn_col(lai_col),
    hmax_fn      = .fn_col("Hmax"),
    fcover_fn    = .fn_col("fCover"),
    lad_fn       = lad_fn,
    phenology_fn = pheno_fn
  )
}

# ---- Public: build full scenario matrix -------------------------------------

#' Build the Chapter 3 scenario matrix.
#'
#' @title Chapter 3 scenario matrix
#' @description Builds all scenarios covering the three experimental axes (LAI
#'   source, phenology mode, LAD shape) then filters to the requested mode.
#'
#'   `"minimal"` returns 5 scenarios: the Chapter 1 LiDAR reference
#'   (`STATIC_ALS`) and a clean 2×2 design (LAI correction: ATBD vs DOPT) ×
#'   (phenology: static vs dynamic).  This is the recommended mode for thesis
#'   defence with a 6-month timeline.
#'
#'   `"medium"` adds `DYN_ALS` (tests H2 shape-stability assumption) and
#'   `DYN_S2_DOPT_LADOPT` (isolates LAD vertical shape contribution).
#'
#'   `"full"` returns all 13–15 scenarios (original design).
#'
#' @param ts_list    Named list from `build_all_s2_ts()`:
#'   `$atbd, $dopt, $rescaled, $rescaled_dopt, $rf` — each a named list of
#'   per-plot smoothed time series `data.frame(doy, lai)`. NULL elements skip
#'   those DYN scenarios.
#' @param list_year  Integer vector of simulation years (default `c(2021, 2022)`).
#' @param d_opt      Effective canopy depth for LAD_DOPT scenarios (default 6 m,
#'   Pareto criterion for Blois).
#' @param mode       Scenario subset: `"minimal"` (default), `"medium"`, or
#'   `"full"`. Set via `CFG_C3$scenarios_mode`.
#' @return Named list of scenario lists.
#' @seealso [build_lai_correction_table()], [build_all_s2_ts()]
make_all_scenarios_c3 <- function(ts_list, list_year = c(2021, 2022), d_opt = 6,
                                   mode = "minimal") {
  lad_full <- make_lad_real
  lad_dopt <- make_lad_dopt_factory(d_opt = d_opt)

  pheno_atbd          <- .make_pheno_fn(ts_list$atbd,          list_year)
  pheno_dopt          <- .make_pheno_fn(ts_list$dopt,          list_year)
  pheno_rescaled      <- .make_pheno_fn(ts_list$rescaled,      list_year)
  pheno_rescaled_dopt <- .make_pheno_fn(ts_list$rescaled_dopt, list_year)
  pheno_rf            <- .make_pheno_fn(ts_list$rf,            list_year)
  pheno_rf_dopt       <- if (!is.null(ts_list$rf_dopt))
    .make_pheno_fn(ts_list$rf_dopt, list_year) else NULL
  pheno_rf_pheno      <- if (!is.null(ts_list$rf_pheno))
    .make_pheno_fn(ts_list$rf_pheno, list_year) else NULL
  pheno_rf_pheno_nodopt <- if (!is.null(ts_list$rf_pheno_nodopt))
    .make_pheno_fn(ts_list$rf_pheno_nodopt, list_year) else NULL
  pheno_dopt_pure       <- if (!is.null(ts_list$rescaled_dopt_pure))
    .make_pheno_fn(ts_list$rescaled_dopt_pure, list_year) else NULL
  # Real annual S2 phenology (19 dates, Whittaker + greenness-fraction × full
  # LiDAR magnitude) — captures observed leaf-out/senescence, not extrapolated.
  pheno_annual          <- if (!is.null(ts_list$annual))
    .make_pheno_fn(ts_list$annual, list_year) else NULL
  pheno_rf_clhs         <- if (!is.null(ts_list$rf_clhs))
    .make_pheno_fn(ts_list$rf_clhs, list_year) else NULL
  pheno_rf_clhs_dopt    <- if (!is.null(ts_list$rf_clhs_dopt))
    .make_pheno_fn(ts_list$rf_clhs_dopt, list_year) else NULL

  sc_list <- list(
    # ---- Static phenology (6 LAI sources × full LAD) -------------------------
    STATIC_ALS          = .make_sc("STATIC_ALS",         "LAI_ALS",         lad_full),
    STATIC_ALS_DOPT     = .make_sc("STATIC_ALS_DOPT",    "LAI_ALS_DOPT",    lad_full),
    STATIC_S2_ATBD      = .make_sc("STATIC_S2_ATBD",     "LAI_S2_ATBD",     lad_full),
    STATIC_S2_DOPT      = .make_sc("STATIC_S2_DOPT",     "LAI_S2_DOPT",     lad_full),
    STATIC_S2_RESCALED  = .make_sc("STATIC_S2_RESCALED", "LAI_S2_RESCALED", lad_full),
    # ---- S2-driven phenology (4 LAI sources × full LAD) ---------------------
    DYN_ALS             = .make_sc("DYN_ALS",            "LAI_ALS",         lad_full, pheno_atbd),
    DYN_S2_ATBD         = .make_sc("DYN_S2_ATBD",        "LAI_S2_ATBD",     lad_full, pheno_atbd),
    DYN_S2_DOPT         = .make_sc("DYN_S2_DOPT",        "LAI_S2_DOPT",     lad_full, pheno_dopt),
    DYN_S2_RESCALED     = .make_sc("DYN_S2_RESCALED",    "LAI_S2_RESCALED", lad_full, pheno_rescaled),
    # Sensor-consistency temporal correction: S2 phenological shape (d_opt
    # series) with amplitude anchored to LAI_ALS_dopt at the LiDAR date.
    DYN_S2_RESCALED_DOPT = .make_sc("DYN_S2_RESCALED_DOPT", "LAI_ALS_DOPT", lad_full, pheno_rescaled_dopt),
    # ---- d_opt-truncated LAD (paired scenarios: correction × LAD shape) ------
    # S2_DOPT + LAD_DOPT : fully coherent d_opt pairing (both amplitude and shape)
    STATIC_S2_DOPT_LADOPT = .make_sc("STATIC_S2_DOPT_LADOPT", "LAI_S2_DOPT",  lad_dopt),
    DYN_S2_DOPT_LADOPT    = .make_sc("DYN_S2_DOPT_LADOPT",    "LAI_S2_DOPT",  lad_dopt, pheno_dopt),
    # S2_ATBD + LAD_DOPT : tests if LAD shape alone corrects without amplitude fix
    STATIC_S2_ATBD_LADOPT = .make_sc("STATIC_S2_ATBD_LADOPT", "LAI_S2_ATBD",  lad_dopt),
    # ALS + LAD_DOPT + S2 pheno : decouples amplitude (LiDAR), shape-vertical (dopt),
    #                               and shape-temporal (S2) into separate contributions
    DYN_ALS_LADOPT        = .make_sc("DYN_ALS_LADOPT",         "LAI_ALS",      lad_dopt, pheno_atbd)
  )

  # Add RF scenarios only if RF time series is available (full mode only)
  if (!is.null(pheno_rf)) {
    sc_list$STATIC_RF <- .make_sc("STATIC_RF", "LAI_RF", lad_full)
    sc_list$DYN_RF    <- .make_sc("DYN_RF",    "LAI_RF", lad_full, pheno_rf)
  }
  # ML sensor-consistency correction (RF trained on d_opt LAI): ML analogue of
  # DYN_S2_RESCALED_DOPT. STATIC_RF_DOPT is its matched static counterpart.
  if (!is.null(pheno_rf_dopt)) {
    sc_list$STATIC_RF_DOPT <- .make_sc("STATIC_RF_DOPT", "LAI_RF_DOPT", lad_full)
    sc_list$DYN_RF_DOPT    <- .make_sc("DYN_RF_DOPT",    "LAI_RF_DOPT", lad_full, pheno_rf_dopt)
  }
  # Structure-conditioned S2 phenology emulator (target=LAI_S2). Magnitude is set
  # by the emulated series (~ATBD); compare against DYN_S2_ATBD (raw-smoothed S2)
  # and STATIC_S2_ATBD. lai_col only seeds lai_max_per_cohort (moot for DYN).
  if (!is.null(pheno_rf_pheno)) {
    sc_list$DYN_RF_PHENO <- .make_sc("DYN_RF_PHENO", "LAI_S2_ATBD", lad_full, pheno_rf_pheno)
  }
  if (!is.null(pheno_rf_pheno_nodopt)) {
    sc_list$DYN_RF_PHENO_NODOPT <- .make_sc("DYN_RF_PHENO_NODOPT", "LAI_S2_ATBD",
                                            lad_full, pheno_rf_pheno_nodopt)
  }
  # Strictly-consistent d_opt pair: S2_DOPT shape at LAI_ALS_dopt magnitude.
  if (!is.null(pheno_dopt_pure)) {
    sc_list$DYN_S2_DOPT_PURE <- .make_sc("DYN_S2_DOPT_PURE", "LAI_ALS_DOPT",
                                         lad_full, pheno_dopt_pure)
  }
  # Real annual S2 phenology × full LiDAR magnitude (the Part-1 pipeline output).
  if (!is.null(pheno_annual)) {
    sc_list$DYN_S2_ANNUAL <- .make_sc("DYN_S2_ANNUAL", "LAI_ALS", lad_full, pheno_annual)
  }
  # Naive baseline: a single, truly-constant LiDAR LAI all year (no phenology) —
  # physically wrong (LAI varies) but may win on flat-canopy windows.
  pheno_const <- if (!is.null(ts_list$const_als)) .make_pheno_fn(ts_list$const_als, list_year) else NULL
  if (!is.null(pheno_const)) {
    sc_list$CONST_ALS <- .make_sc("CONST_ALS", "LAI_ALS", lad_full, pheno_const)
  }
  # Pure-satellite naive baseline (Ch3): raw S2 LAI magnitude (saturated) ×
  # FORMS-H canopy height × UNIFORM LAD (no ALS vertical profile). Quantifies the
  # cost of having NO airborne LiDAR (magnitude + LAD shape + height all degraded).
  pheno_naive <- if (!is.null(ts_list$naive_s2)) .make_pheno_fn(ts_list$naive_s2, list_year) else NULL
  if (!is.null(pheno_naive)) {
    sc_list$NAIVE_S2_FORMSH <- list(
      name = "NAIVE_S2_FORMSH", lai_fn = .fn_col("LAI_S2_ATBD"),
      hmax_fn = .fn_col("FORMS_H"), fcover_fn = .fn_col("fCover"),
      lad_fn = make_lad_uniform, phenology_fn = pheno_naive)
  }
  # cLHS-trained transfer correction applied to HOBO (structure vs +d_opt).
  # lai_col seeds lai_max only (moot for DYN; magnitude set by the series).
  if (!is.null(pheno_rf_clhs))
    sc_list$DYN_RF_CLHS <- .make_sc("DYN_RF_CLHS", "LAI_S2_ATBD", lad_full, pheno_rf_clhs)
  if (!is.null(pheno_rf_clhs_dopt))
    sc_list$DYN_RF_CLHS_DOPT <- .make_sc("DYN_RF_CLHS_DOPT", "LAI_S2_ATBD", lad_full, pheno_rf_clhs_dopt)

  # Drop scenarios with NULL phenology_fn that came from missing ts_list elements
  sc_list <- Filter(Negate(is.null), sc_list)

  # ---- Filter by mode ----------------------------------------------------------
  keep <- switch(mode,
    minimal = c(
      "STATIC_ALS",
      "STATIC_S2_ATBD", "STATIC_S2_DOPT",
      "DYN_S2_ATBD",    "DYN_S2_DOPT"
    ),
    medium = c(
      "STATIC_ALS",
      "STATIC_S2_ATBD", "STATIC_S2_DOPT",
      "DYN_ALS",
      "DYN_S2_ATBD",    "DYN_S2_DOPT",
      "DYN_S2_DOPT_LADOPT"
    ),
    full = names(sc_list),
    stop(sprintf("Unknown scenarios mode '%s'. Use 'minimal', 'medium', or 'full'.", mode))
  )
  sc_list <- sc_list[intersect(keep, names(sc_list))]

  cat(sprintf("  %d Chapter 3 scenarios [mode='%s']: %s\n",
              length(sc_list), mode, paste(names(sc_list), collapse = ", ")))
  sc_list
}
