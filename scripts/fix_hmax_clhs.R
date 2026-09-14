# ==============================================================================
# Update Hmax in the cLHS sample with the buffer-MAX rule (= legacy Oct 2025).
#
# The current `clhs_sample.rds` has Hmax from the agg-MEAN 20 m raster, which
# systematically under-estimates canopy height. We keep the plot selection
# and cluster assignments unchanged (they were already validated), and only
# replace the Hmax values with the MAX of the raw 10 m raster within each
# 20 m cell (cell-MAX, equivalent to the buffer-MAX HOBO fix).
#
# Saves :
#   - out_files/Sensitivity_Analysis/clhs_sample_hmaxfix.rds      (updated)
#   - out_files/Sensitivity_Analysis/clhs_sample_floor05_hmaxfix.rds (updated)
#
# The original files are kept untouched as backup.
# ==============================================================================

suppressMessages({
  library(data.table); library(terra); library(sf); library(here); library(cli)
})

cli_h1("Updating Hmax in cLHS sample with cell-MAX aggregation")

# 1. Load raw Hmax raster (10 m)
r_hmax_raw <- terra::rast(here::here("in_files", "max_res_10_m.tif"))
cli_alert("Raw Hmax raster: {ncol(r_hmax_raw)}×{nrow(r_hmax_raw)} cells, res={res(r_hmax_raw)[1]} m")

# 2. Aggregate with MAX (factor 2 = 10 m → 20 m, same as cLHS native resolution)
r_hmax_max20 <- terra::aggregate(r_hmax_raw, fact = 2, fun = "max", na.rm = TRUE)
cli_alert("Aggregated Hmax 20 m (MAX): {ncol(r_hmax_max20)}×{nrow(r_hmax_max20)} cells, res={res(r_hmax_max20)[1]} m")

# 3. Compare to current MEAN aggregation (sanity check)
r_hmax_mean20 <- terra::aggregate(r_hmax_raw, fact = 2, fun = "mean", na.rm = TRUE)

# 4. For each cLHS sample, extract the new Hmax value
update_sample <- function(in_path, out_path) {
  if (!file.exists(in_path)) {
    cli_alert_warning("Skipping (not found): {.path {in_path}}")
    return(invisible(NULL))
  }
  df <- as.data.table(readRDS(in_path))
  cli_alert("Updating {.path {basename(in_path)}}: n={nrow(df)} plots")

  # Build point geometry from x/y (already in UTM31N)
  pts <- terra::vect(df, geom = c("x", "y"), crs = "EPSG:32631")

  # Extract new Hmax (cell-MAX) and old Hmax (cell-MEAN) for sanity check
  hmax_new  <- terra::extract(r_hmax_max20,  pts)[, 2]
  hmax_chk  <- terra::extract(r_hmax_mean20, pts)[, 2]

  # Compare
  m <- data.table(old = df$Hmax, recomputed_mean = hmax_chk, new = hmax_new)
  m[, diff_new := new - old]
  cli_alert("  Hmax stats (across {nrow(m)} cLHS plots):")
  cat(sprintf("    old (current MEAN) : mean=%.2f m, range [%.2f, %.2f]\n",
              mean(m$old, na.rm=TRUE),
              min(m$old, na.rm=TRUE), max(m$old, na.rm=TRUE)))
  cat(sprintf("    new (MAX agg 20 m) : mean=%.2f m, range [%.2f, %.2f]\n",
              mean(m$new, na.rm=TRUE),
              min(m$new, na.rm=TRUE), max(m$new, na.rm=TRUE)))
  cat(sprintf("    diff (new - old)   : mean=%+.2f m, n_increased=%d/%d\n",
              mean(m$diff_new, na.rm=TRUE),
              sum(m$diff_new > 0, na.rm=TRUE), nrow(m)))

  # Replace Hmax in-place
  df_new <- copy(df)
  df_new[, Hmax := hmax_new]

  saveRDS(df_new, out_path)
  cli_alert_success("Saved {.path {basename(out_path)}}")
  invisible(df_new)
}

sa_dir <- here::here("out_files/Sensitivity_Analysis")
update_sample(file.path(sa_dir, "clhs_sample.rds"),
              file.path(sa_dir, "clhs_sample_hmaxfix.rds"))
update_sample(file.path(sa_dir, "clhs_sample_floor05.rds"),
              file.path(sa_dir, "clhs_sample_floor05_hmaxfix.rds"))

cli_h2("Done")
cli_alert("Original files unchanged. New files saved with `_hmaxfix` suffix.")
cli_alert("To use them, point downstream scripts to the `_hmaxfix.rds` versions.")
