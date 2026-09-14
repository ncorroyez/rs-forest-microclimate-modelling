# ==============================================================================
# Chapter 1 — Scenario registry
#
# Each scenario is a named list:
#   $name       : character, used as directory name for outputs
#   $lai_fn     : function(plot_row) → numeric
#   $hmax_fn    : function(plot_row) → numeric
#   $fcover_fn  : function(plot_row) → numeric
#   $lad_fn     : function(plot_row, hmax, lai) → data.frame(height, density)
# ==============================================================================

# ---- Helper closures ---------------------------------------------------------

fn_real <- function(varname) { force(varname); function(plot_row) as.numeric(plot_row[[varname]]) }
fn_mean <- function(df_sample, varname) { m <- mean(df_sample[[varname]], na.rm = TRUE); function(plot_row) m }

# ---- Reference scenario (all inputs from LiDAR) ------------------------------

#' Build the full-real reference scenario (all inputs from LiDAR).
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named scenario list.
scenario_reference <- function(df_sample) {
  list(name = "REF_all_real", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"),
       fcover_fn = fn_real("fCover"), lad_fn = make_lad_real)
}

# ---- H2 scenarios (uniform vs real LAD) -------------------------------------

#' H2 scenarios: real vs uniform LAD, all other inputs held real.
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named list of two scenarios.
scenarios_h2_uniform_vs_real <- function(df_sample) {
  list(
    REF_real_LAD = list(name = "H2_real_LAD",    lai_fn = fn_real("LAI"),
                        hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"),
                        lad_fn = make_lad_real),
    UNIFORM_LAD  = list(name = "H2_uniform_LAD", lai_fn = fn_real("LAI"),
                        hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"),
                        lad_fn = make_lad_uniform)
  )
}

#' H2 scenarios: real vs cluster-type LAD.
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named list of two scenarios.
scenarios_h2_cluster_type <- function(df_sample) {
  lad_cluster_fn <- make_lad_cluster_type_factory(df_sample)
  list(
    REF_real_LAD     = list(
      name = "H2_real_LAD",
      lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"),
      fcover_fn = fn_real("fCover"), lad_fn = make_lad_real
    ),
    CLUSTER_TYPE_LAD = list(
      name = "H2_cluster_type_LAD",
      lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"),
      fcover_fn = fn_real("fCover"), lad_fn = lad_cluster_fn
    )
  )
}

# ---- H1 scenarios (forward inclusion) ----------------------------------------

#' H1 forward inclusion scenarios: SC0 → SC4 (null baseline to full-real).
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named list of 5 scenarios.
scenarios_h1_forward <- function(df_sample) {
  lai_mean      <- fn_mean(df_sample, "LAI")
  hmax_mean_val <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  hmax_mean     <- function(plot_row) hmax_mean_val
  fcover_absent <- function(plot_row) 1   # canopée fermée homogène, pas de clumping
  lad_uniform   <- make_lad_uniform
  list(
    Null_baseline   = list(name = "H1f_0_Null_baseline",       lai_fn = lai_mean,       hmax_fn = hmax_mean,       fcover_fn = fcover_absent, lad_fn = lad_uniform),
    LAI_only        = list(name = "H1f_1_LAI_only",             lai_fn = fn_real("LAI"), hmax_fn = hmax_mean,       fcover_fn = fcover_absent, lad_fn = lad_uniform),
    LAI_Hmax        = list(name = "H1f_2_LAI_Hmax",             lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fcover_absent, lad_fn = lad_uniform),
    LAI_Hmax_fCover = list(name = "H1f_3_LAI_Hmax_fCover",      lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = lad_uniform),
    Full_real       = list(name = "H1f_4_Full_real",             lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_real)
  )
}

# ---- H1 scenarios (leave-one-out) --------------------------------------------

#' H1 leave-one-out scenarios: each variable degraded in turn.
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named list of 5 scenarios.
scenarios_h1_loo <- function(df_sample) {
  lai_mean      <- fn_mean(df_sample, "LAI")
  hmax_mean_val <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  hmax_mean     <- function(plot_row) hmax_mean_val
  fcover_absent <- function(plot_row) 1   # cohérent avec baseline forward
  lad_mean_fn   <- make_lad_mean_factory(df_sample)
  list(
    drop_LAI_mean      = list(name = "H1l_dropLAI_mean",      lai_fn = lai_mean,       hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_real),
    drop_Hmax_mean     = list(name = "H1l_dropHmax_mean",     lai_fn = fn_real("LAI"), hmax_fn = hmax_mean,       fcover_fn = fn_real("fCover"), lad_fn = make_lad_real),
    drop_fCover_mean   = list(name = "H1l_dropfCover_mean",   lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fcover_absent,     lad_fn = make_lad_real),
    drop_LAD_uniform   = list(name = "H1l_dropLAD_uniform",   lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = make_lad_uniform),
    drop_LAD_meanShape = list(name = "H1l_dropLAD_meanShape", lai_fn = fn_real("LAI"), hmax_fn = fn_real("Hmax"), fcover_fn = fn_real("fCover"), lad_fn = lad_mean_fn)
  )
}

# ---- H1 full factorial (2^4 = 16 scenarios) ----------------------------------

#' H1 full factorial scenarios (2^4 lattice for Shapley computation).
#'
#' Generates all 16 combinations of {real, mean/absent/uniform} for
#' {LAI, Hmax, fCover, LAD}.
#'
#' @param df_sample cLHS sample dataframe.
#' @return Named list of 16 scenarios.
scenarios_h1_factorial <- function(df_sample) {
  lai_mean      <- fn_mean(df_sample, "LAI")
  hmax_mean_val <- round(mean(floor(df_sample$Hmax) + 1, na.rm = TRUE))
  hmax_mean     <- function(plot_row) hmax_mean_val
  fcover_absent <- function(plot_row) 1   # cohérent avec baseline forward/LOO
  grid <- expand.grid(LAI = c("real", "mean"), Hmax = c("real", "mean"),
                      fCover = c("real", "absent"), LAD = c("real", "uniform"),
                      stringsAsFactors = FALSE)
  scenarios <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g  <- grid[i, ]
    nm <- sprintf("H1F_%s_%s_%s_%s",
                  substr(g$LAI, 1, 1), substr(g$Hmax, 1, 1),
                  substr(g$fCover, 1, 1), substr(g$LAD, 1, 1))
    scenarios[[i]] <- list(
      name      = nm,
      lai_fn    = if (g$LAI    == "real")   fn_real("LAI")    else lai_mean,
      hmax_fn   = if (g$Hmax   == "real")   fn_real("Hmax")   else hmax_mean,
      fcover_fn = if (g$fCover == "real")   fn_real("fCover") else fcover_absent,
      lad_fn    = if (g$LAD    == "real")   make_lad_real     else make_lad_uniform,
      grid_row  = g
    )
  }
  setNames(scenarios, vapply(scenarios, `[[`, "", "name"))
}
