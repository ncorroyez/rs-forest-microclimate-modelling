# ==============================================================================
# V3 HOBO : extract inputs (LAI, Hmax, fCover, LAD) from a 25 m radius
# buffer around each HOBO sensor (legacy `main_Blois_detailed.R` style).
# Hmax = buffer MAX, others = buffer MEAN.  LAI and LAD layers are × 2
# (PAD → LAI conversion).
#
# Runs 53 HOBO × 5 forward-selection coalitions = 265 sims.
# Then regenerates the HOBO forward-selection figure with the buffer25 data.
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

cli_h1("V3 HOBO : 25 m buffer extraction (legacy-style)")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
  hobo_buffer_mode = "buffer25",
  raw_raster_dir = CFG$in_dir))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
cli_alert("HOBO inputs (buffer25, n={nrow(df_hobo)}):")
cli_alert("  LAI    mean={round(mean(df_hobo$LAI),2)}, range [{round(min(df_hobo$LAI),2)}, {round(max(df_hobo$LAI),2)}]")
cli_alert("  Hmax   mean={round(mean(df_hobo$Hmax),2)} m, range [{round(min(df_hobo$Hmax),2)}, {round(max(df_hobo$Hmax),2)}]")
cli_alert("  fCover mean={round(mean(df_hobo$fCover),2)}, range [{round(min(df_hobo$fCover),2)}, {round(max(df_hobo$fCover),2)}]")

# Forward order : reuse the V2 archetype ranking (LAD → fCover → Hmax → LAI)
forward_order <- c("LAD", "fCover", "Hmax", "LAI")
cli_alert("Forward order : {paste(forward_order, collapse=' → ')}")

var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}
cli_alert("Forward bit sequence : {paste(bits_seq, collapse=' → ')}")

# Scenarios from V2 cLHS (with PAD-fixed)
df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
fac_scs_hobo <- build_factorial_scenarios_archetypes(df_floor)
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

out_root_hobo <- here::here("out_files/musica_hobo_v3")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (k in seq_along(bits_seq)) {
    bit <- bits_seq[k]
    sc_name <- bit2scn[bit]
    sc <- fac_scs_hobo[[sc_name]]
    out_dir <- file.path(out_root_hobo, bit)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_nc <- file.path(out_dir, sprintf("musica_out_HOBO_%s.nc", pr$id_plot))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    hobo_jobs[[length(hobo_jobs) + 1L]] <- list(
      pr = pr, sc = sc, out_nc = out_nc,
      lbl = sprintf("%s/%s", bit, pr$id_plot))
  }
}
cli_alert("HOBO jobs : {length(hobo_jobs)}")

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
cli_alert_success("HOBO V3 done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')), 1)} min")

# Compute slope-based regime
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
DT_macro <- as.data.table(df_macro)[, .(date, Tmax_macro)]

regime_rows <- list()
for (f in list.files(file.path(out_root_hobo, "1111"), pattern="\\.nc$", full.names=TRUE)) {
  id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
  res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                    error=function(e) NULL)
  if (is.null(res) || nrow(res) == 0) next
  d <- as.data.table(res)
  if (!"Tmax_macro" %in% names(d)) d <- merge(d, DT_macro, by="date")
  fit <- lm(Tmax_micro ~ Tmax_macro, data = d)
  regime_rows[[id]] <- data.frame(id_plot = id, slope = coef(fit)[2])
}
DT_reg <- as.data.table(do.call(rbind, regime_rows))
cli_h2("Slope-based regime (V3 buffer25 REF)")
cat(sprintf("  Buffer (slope < 1)    : %d / %d\n", sum(DT_reg$slope < 1), nrow(DT_reg)))
cat(sprintf("  Amplification (slope > 1) : %d / %d\n", sum(DT_reg$slope > 1), nrow(DT_reg)))
cat(sprintf("  slope mean=%.3f median=%.3f range=[%.3f, %.3f]\n",
            mean(DT_reg$slope), median(DT_reg$slope), min(DT_reg$slope), max(DT_reg$slope)))
