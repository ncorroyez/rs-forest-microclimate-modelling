# ==============================================================================
# V4 HOBO : use the LEGACY CSV inputs directly :
#   - metrics_results_25.csv  : pai (= LAI), max (= Hmax), clumping (= fCover)
#   - allometry_profiles.csv  : LAD profile per height (raw, 25 m point cloud)
#
# This bypasses raster extraction and reproduces exactly what
# `main_Blois_detailed.R` fed to MuSICA in February/October 2025.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(sf); library(terra); library(future); library(furrr)
  library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/lad.R"))
source(here::here("R/validation.R"))
source(here::here("R/h1_shapley_archetypes.R"))

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

cli_h1("V4 HOBO : LEGACY CSV inputs (metrics_results_25 + allometry_profiles)")

# 1. Load legacy CSVs
metr <- fread(here::here("out_files/radius_test/metrics_results_25.csv"))
allom <- fread(here::here("allometry_profiles.csv"))
allom <- allom[radius == 25]

# Filter to our 53 HOBO sensors
ids_keep <- setdiff(metr$id_plot, CFG$ids_to_remove)
cli_alert("HOBO sensors with legacy CSV data : {length(ids_keep)}")

# 2. Build plot rows with legacy values
df_hobo <- list()
for (id in ids_keep) {
  m <- metr[id_plot == id]
  a <- allom[id_plot == id]
  if (nrow(m) == 0 || nrow(a) == 0) next
  # LAD layers as columns
  lad_cols_raw <- grep("^h_", names(a), value = TRUE)
  z_real <- as.numeric(gsub("h_", "", lad_cols_raw))
  d <- as.numeric(a[1, ..lad_cols_raw])
  d[is.na(d)] <- 0

  # Build LAD_Layer_* columns (only non-zero entries)
  lad_named <- setNames(d, paste0("LAD_Layer_", z_real))

  # x, y from HOBO geojson
  hobo_pts <- sf::st_read(CFG$hobo_geojson, quiet = TRUE)
  c <- sf::st_coordinates(hobo_pts[hobo_pts$id_plot == id, ])

  row <- as.data.frame(c(
    list(id_plot = id, x = c[1, "X"], y = c[1, "Y"],
          LAI = m$pai, Hmax = m$max, fCover = m$clumping),
    as.list(lad_named)))
  df_hobo[[id]] <- row
}
df_hobo <- bind_rows(df_hobo)
# Floor fCover (per the floor05 convention used in the deck)
df_hobo$fCover[df_hobo$fCover < 0.5] <- 0.5
cli_alert("Built {nrow(df_hobo)} HOBO inputs from legacy CSVs")
cli_alert("Stats : LAI mean={round(mean(df_hobo$LAI),2)}, range [{round(min(df_hobo$LAI),2)}, {round(max(df_hobo$LAI),2)}]")
cli_alert("        Hmax mean={round(mean(df_hobo$Hmax),2)} m")

# 3. Forward order from V2 archetypes : LAD → fCover → Hmax → LAI
forward_order <- c("LAD", "fCover", "Hmax", "LAI")
var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}
cli_alert("Forward bits : {paste(bits_seq, collapse=' → ')}")

# 4. Scenarios with V2 archetype baselines
df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs <- build_factorial_scenarios_archetypes(df_floor)
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

out_root <- here::here("out_files/musica_hobo_v4_legacy")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (k in seq_along(bits_seq)) {
    bit <- bits_seq[k]
    sc_name <- bit2scn[bit]
    sc <- fac_scs[[sc_name]]
    out_dir <- file.path(out_root, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("Jobs : {length(hobo_jobs)}")

t_start <- Sys.time()
res <- future_walk(hobo_jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$pr, job$sc, job$out_nc,
                    CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("V4 sims done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")

# 5. Compute slope-based regime
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
DT_macro <- as.data.table(df_macro)[, .(date, Tmax_macro)]
regime_rows <- list()
for (f in list.files(file.path(out_root, "1111"), pattern="\\.nc$", full.names=TRUE)) {
  id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
  res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                    error=function(e) NULL)
  if (is.null(res) || nrow(res) == 0) next
  d <- as.data.table(res)
  if (!"Tmax_macro" %in% names(d)) d <- merge(d, DT_macro, by="date")
  fit <- lm(Tmax_micro ~ Tmax_macro, data = d)
  regime_rows[[id]] <- data.frame(id_plot = id, slope = coef(fit)[2],
                                     dT_mean = mean(d$Tmax_micro - d$Tmax_macro, na.rm = TRUE))
}
DT_reg <- as.data.table(do.call(rbind, regime_rows))
cli_h2("V4 (legacy CSV) slope-based regime")
cat(sprintf("  Buffer (slope < 1)    : %d / %d\n", sum(DT_reg$slope < 1), nrow(DT_reg)))
cat(sprintf("  Amplification (slope > 1) : %d / %d\n", sum(DT_reg$slope > 1), nrow(DT_reg)))
cat(sprintf("  slope mean=%.3f median=%.3f range=[%.3f, %.3f]\n",
            mean(DT_reg$slope), median(DT_reg$slope), min(DT_reg$slope), max(DT_reg$slope)))
cat(sprintf("  ΔTmax mean=%+.3f median=%+.3f range=[%+.3f, %+.3f]\n",
            mean(DT_reg$dT_mean), median(DT_reg$dT_mean),
            min(DT_reg$dT_mean), max(DT_reg$dT_mean)))
