# ==============================================================================
# LOVB analysis — Module 2 : temporal aggregation (mean + P90)
#
# Aggregates daily Delta_Tmax to per-plot per-coalition scalars.
#   - Tmax_mean : mean over 122 summer days
#   - Tmax_P90  : 90th percentile using quantile(x, 0.9, type = 7)
#
# For cLHS, joins Cluster + structural traits (LAI, Hmax, fCover, FPC1) from
# clhs_sample.rds for context plots (Figure 2).
# For archetypes, derives Cluster from the archetype label.
# ==============================================================================

suppressMessages({
  library(data.table)
  library(here)
  library(cli)
})

# ---- 1. cLHS aggregation -----------------------------------------------------
lovb_aggregate_clhs <- function(DT_daily,
                                 clhs_sample_path = here("out_files/Sensitivity_Analysis/clhs_sample.rds"),
                                 cache_path       = here("outputs/lovb/data/DT_agg_cLHS.rds"),
                                 use_cache        = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached cLHS DT_agg from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Aggregate cLHS daily Delta_Tmax (mean + P90)")

  stopifnot(is.data.table(DT_daily),
            all(c("x", "y", "date", "bit_code", "Delta_Tmax") %in% names(DT_daily)))

  # Per-plot per-coalition mean + P90 (quantile type 7 = R default)
  DT_agg <- DT_daily[
    , .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
        Tmax_P90  = quantile(Delta_Tmax, 0.9, type = 7, na.rm = TRUE),
        n_days    = sum(!is.na(Delta_Tmax))),
    by = .(x, y, bit_code)
  ]

  # Join Cluster + structural traits from clhs_sample
  df_s <- as.data.table(readRDS(clhs_sample_path))
  if ("Archetype" %in% names(df_s) && !"Cluster" %in% names(df_s)) {
    setnames(df_s, "Archetype", "Cluster")
  }
  df_s[, x := as.integer(round(x))]
  df_s[, y := as.integer(round(y))]
  traits <- df_s[, .(x, y, Cluster, LAI, Hmax, fCover, FPC1)]
  setkey(traits, x, y)

  DT_agg[, x := as.integer(x)]
  DT_agg[, y := as.integer(y)]
  setkey(DT_agg, x, y)
  DT_agg <- traits[DT_agg, on = c("x", "y"), nomatch = NA]

  n_na_cluster <- sum(is.na(DT_agg$Cluster))
  if (n_na_cluster > 0L)
    cli_alert_warning("{n_na_cluster} rows have no Cluster join (coord mismatch?)")

  setkey(DT_agg, x, y, bit_code)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT_agg, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT_agg)}} rows)")
  DT_agg
}

# ---- 2. Archetypes aggregation -----------------------------------------------
lovb_aggregate_archetypes <- function(DT_daily,
                                       cache_path = here("outputs/lovb/data/DT_agg_archetypes.rds"),
                                       use_cache  = TRUE) {
  if (use_cache && file.exists(cache_path)) {
    cli_alert_info("Loading cached archetypes DT_agg from {.path {cache_path}}")
    return(as.data.table(readRDS(cache_path)))
  }
  cli_h1("Aggregate archetypes daily Delta_Tmax (mean + P90)")

  stopifnot(is.data.table(DT_daily),
            all(c("archetype", "date", "bit_code", "Delta_Tmax") %in% names(DT_daily)))

  DT_agg <- DT_daily[
    , .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
        Tmax_P90  = quantile(Delta_Tmax, 0.9, type = 7, na.rm = TRUE),
        n_days    = sum(!is.na(Delta_Tmax))),
    by = .(archetype, bit_code)
  ]
  # Derive Cluster from archetype label : "Arch_C1" -> "1"
  DT_agg[, Cluster := factor(sub("^Arch_C", "", archetype), levels = c("1", "2", "3", "4"))]
  setkey(DT_agg, archetype, bit_code)

  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(DT_agg, cache_path)
  cli_alert_success("Cached {.path {cache_path}} ({.val {nrow(DT_agg)}} rows)")
  DT_agg
}
