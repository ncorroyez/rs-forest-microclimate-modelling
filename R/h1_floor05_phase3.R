# ==============================================================================
# Phase 3 — Aggregate + merge daily Delta_Tmax for floor05 pipeline.
#
# Strategy : merge existing daily cache (factorial original) with new floor05
# NetCDFs for affected (plot × coalition) tuples only. No full re-read of 6400
# cLHS NCs needed.
#
# Outputs :
#   outputs/lovb_floor05/data/DT_daily_cLHS_floor05.rds  (16 coal × 400 plot)
#   outputs/lovb_floor05/data/DT_agg_cLHS_floor05.rds
#   outputs/lovb_floor05/data/DT_contrib_cLHS_floor05.rds  (wide 16 coalitions)
#   outputs/lovb_floor05/data/DT_daily_archetypes_floor05.rds
#   outputs/lovb_floor05/data/DT_contrib_archetypes_floor05.rds
#   outputs/lovb_floor05/data/DT_daily_HOBO_floor05.rds
#   outputs/lovb_floor05/data/DT_contrib_HOBO_floor05.rds
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli)
  library(ncdf4); library(tidyverse); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/musica.R"))
source(here::here("R/lovb_01_load.R"))
source(here::here("R/lovb_03_contrib.R"))

OUT_DIR <- here::here("outputs/lovb_floor05/data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cli_h1("Phase 3 — merge daily + build floor05 contribs")
t_start <- Sys.time()

# ---- 1. cLHS : merge original + floor05 affected -----------------------------
cli_h2("cLHS daily")
affected_clhs <- readRDS(here::here("out_files/Sensitivity_Analysis/floor05_affected_clhs.rds"))
affected_keys <- affected_clhs$plots[, paste(x, y, sep = "_")]
COAL_FCOVER_REAL <- c("0010","0011","0110","0111","1010","1011","1110","1111")

DT_orig_10 <- readRDS(here::here("outputs/lovb/data/DT_daily_cLHS.rds"))
DT_orig_6  <- readRDS(here::here("outputs/lovb/data/DT_daily_cLHS_missing6.rds"))
DT_orig <- rbind(DT_orig_10, DT_orig_6)
cli_alert("Original cLHS daily : {nrow(DT_orig)} rows ({length(unique(DT_orig$bit_code))} coalitions)")

# Load floor05 daily for affected coalitions only
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
load_one_dir <- function(scenario_dir, bit_code) {
  nc_files <- list.files(scenario_dir, pattern = "\\.nc$", full.names = TRUE)
  nc_files <- nc_files[file.size(nc_files) > 1e6]
  rows <- vector("list", length(nc_files))
  for (j in seq_along(nc_files)) {
    f <- nc_files[j]
    res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq),
                     error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0L) next
    dt <- as.data.table(res)
    bn <- basename(f)
    dt[, x := as.integer(sub(".*_X(\\d+)_Y\\d+\\.nc$", "\\1", bn))]
    dt[, y := as.integer(sub(".*_X\\d+_Y(\\d+)\\.nc$", "\\1", bn))]
    dt[, bit_code := bit_code]
    rows[[j]] <- dt
  }
  rbindlist(rows, fill = TRUE)
}

cli_alert("Loading floor05 affected daily ({length(COAL_FCOVER_REAL)} coalitions × {length(affected_keys)} plots) ...")
DT_floor_aff <- list()
for (b in COAL_FCOVER_REAL) {
  sc <- paste0("H1F_", lovb_bit_to_suffix(b))
  dir <- here::here("out_files/H1_factorial_floor05", sc)
  if (!dir.exists(dir)) { cli_alert_warning("Missing : {sc}"); next }
  DT_floor_aff[[b]] <- load_one_dir(dir, b)
}
DT_floor_aff <- rbindlist(DT_floor_aff, fill = TRUE)
DT_floor_aff <- DT_floor_aff[, .(x, y, date, bit_code, Delta_Tmax)]
cli_alert("floor05 affected daily : {nrow(DT_floor_aff)} rows")

# Merge : remove affected (plot × coalition) from orig, then bind floor05
DT_orig[, key := paste(x, y, sep = "_")]
DT_orig_keep <- DT_orig[!(bit_code %in% COAL_FCOVER_REAL & key %in% affected_keys)]
DT_orig_keep[, key := NULL]
DT_clhs_floor <- rbind(DT_orig_keep, DT_floor_aff)
setkey(DT_clhs_floor, x, y, bit_code, date)
cli_alert("Final cLHS daily floor05 : {nrow(DT_clhs_floor)} rows ({length(unique(DT_clhs_floor$bit_code))} coalitions)")
saveRDS(DT_clhs_floor, file.path(OUT_DIR, "DT_daily_cLHS_floor05.rds"))

# Aggregate to (plot, coalition)
DT_agg_clhs <- DT_clhs_floor[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
                                   Tmax_P90  = quantile(Delta_Tmax, 0.9,
                                                          na.rm = TRUE, type = 7),
                                   n_days = .N),
                              by = .(x, y, bit_code)]
saveRDS(DT_agg_clhs, file.path(OUT_DIR, "DT_agg_cLHS_floor05.rds"))
cli_alert("Aggregated : {nrow(DT_agg_clhs)} (plot, coalition) rows")

# Wide DT with all 16 coalitions
DT_w <- dcast(DT_agg_clhs, x + y ~ bit_code,
                value.var = c("Tmax_mean", "Tmax_P90"))
df_sample_floor <- as.data.table(readRDS(here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05.rds")))
DT_traits <- df_sample_floor[, .(x, y, Cluster, LAI, Hmax, fCover, FPC1)]
DT_full <- merge(DT_traits, DT_w, by = c("x","y"))
saveRDS(DT_full, file.path(OUT_DIR, "DT_contrib_cLHS_floor05_16coalitions.rds"))
cli_alert_success("Saved DT_contrib_cLHS_floor05_16coalitions.rds ({nrow(DT_full)} plots)")

# Also build the "10 coalitions" wide table (for LOVB/LVA via lovb_contrib_core)
# Map bit_code to label for compatibility with lovb_03_contrib code
DT_agg_for_contrib <- DT_agg_clhs[bit_code %in% c("0000","0001","0010","0100",
                                                     "1000","0111","1011","1101",
                                                     "1110","1111")]
DT_agg_for_contrib <- merge(df_sample_floor[, .(x, y, Cluster, LAI, Hmax, fCover, FPC1)],
                              DT_agg_for_contrib, by = c("x","y"))
DT_contrib_clhs <- .lovb_contrib_core(DT_agg_for_contrib,
                                         id_cols = c("x","y","Cluster","LAI",
                                                       "Hmax","fCover","FPC1"))
saveRDS(DT_contrib_clhs, file.path(OUT_DIR, "DT_contrib_cLHS_floor05.rds"))
cli_alert_success("Saved DT_contrib_cLHS_floor05.rds ({nrow(DT_contrib_clhs)} plots)")

# ---- 2. Archetypes : C1 floor05 + C2/C3/C4 original --------------------------
cli_h2("Archetypes daily")
DT_arch_orig <- readRDS(here::here("outputs/lovb/data/DT_daily_archetypes.rds"))
cli_alert("Original arch daily : {nrow(DT_arch_orig)} rows")
# Need to load C1 from floor05 for 16 coalitions
arch_c1_dir <- here::here("out_files/H1_archetypes_floor05/Arch_C1")
arch_files <- list.files(arch_c1_dir, pattern = "\\.nc$", full.names = TRUE)
arch_files <- arch_files[file.size(arch_files) > 1e6]
cli_alert("C1 floor05 : {length(arch_files)} NCs")
arch_rows <- vector("list", length(arch_files))
for (j in seq_along(arch_files)) {
  f <- arch_files[j]
  res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq),
                   error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) next
  dt <- as.data.table(res)
  # Get bit_code from filename : musica_out_Arch_C1_ARCH_<suff>.nc
  bn <- basename(f)
  suff <- sub(".*_ARCH_(.+)\\.nc$", "\\1", bn)
  # Suffix to bit_code reverse
  parts <- strsplit(suff, "_")[[1]]
  bit <- paste0(ifelse(parts[1] == "r", "1", "0"),
                  ifelse(parts[2] == "r", "1", "0"),
                  ifelse(parts[3] == "r", "1", "0"),
                  ifelse(parts[4] == "r", "1", "0"))
  dt[, archetype := "Arch_C1"]
  dt[, bit_code  := bit]
  arch_rows[[j]] <- dt
}
DT_arch_c1_floor <- rbindlist(arch_rows, fill = TRUE)
DT_arch_c1_floor <- DT_arch_c1_floor[, .(archetype, date, bit_code, Delta_Tmax)]
cli_alert("C1 floor05 daily : {nrow(DT_arch_c1_floor)} rows ({length(unique(DT_arch_c1_floor$bit_code))} coalitions)")
# Drop original C1 rows, replace with floor05
DT_arch_orig_keep <- DT_arch_orig[archetype != "Arch_C1"]
DT_arch_floor <- rbind(DT_arch_orig_keep, DT_arch_c1_floor, fill = TRUE)
setkey(DT_arch_floor, archetype, bit_code, date)
saveRDS(DT_arch_floor, file.path(OUT_DIR, "DT_daily_archetypes_floor05.rds"))
cli_alert_success("Saved DT_daily_archetypes_floor05.rds ({nrow(DT_arch_floor)} rows)")

DT_agg_arch <- DT_arch_floor[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
                                   Tmax_P90  = quantile(Delta_Tmax, 0.9,
                                                          na.rm = TRUE, type = 7),
                                   n_days = .N),
                              by = .(archetype, bit_code)]
# Cluster column for arch
DT_agg_arch[, Cluster := sub("Arch_C", "", archetype)]
DT_agg_arch[, Cluster := factor(Cluster, levels = c("1","2","3","4"))]
# Build wide
DT_contrib_arch <- .lovb_contrib_core(DT_agg_arch[bit_code %in% c("0000","0001","0010","0100",
                                                                       "1000","0111","1011","1101",
                                                                       "1110","1111")],
                                          id_cols = c("archetype","Cluster"))
saveRDS(DT_contrib_arch, file.path(OUT_DIR, "DT_contrib_archetypes_floor05.rds"))
cli_alert_success("Saved DT_contrib_archetypes_floor05.rds ({nrow(DT_contrib_arch)} archetypes)")

# ---- 3. HOBO : 4 sensors × 5 coalitions floor05 + rest original --------------
cli_h2("HOBO daily")
DT_hobo_orig <- readRDS(here::here("outputs/lovb/data/DT_daily_HOBO.rds"))
cli_alert("Original HOBO daily : {nrow(DT_hobo_orig)} rows ({length(unique(DT_hobo_orig$bit_code))} coalitions)")
affected_hobo <- readRDS(here::here("outputs/lovb/data/floor05_affected_hobo.rds"))
hobo_aff_ids <- affected_hobo$sensors$id_plot
hobo_aff_bits <- affected_hobo$coalitions

# Load floor05 HOBO daily
hobo_floor_root <- here::here("out_files/musica_hobo_validation_floor05")
hobo_rows <- list()
for (bit in hobo_aff_bits) {
  sc_dir <- file.path(hobo_floor_root, paste0("HOBO_", bit))
  nc_files <- list.files(sc_dir, pattern = "\\.nc$", full.names = TRUE)
  nc_files <- nc_files[file.size(nc_files) > 1e6]
  for (f in nc_files) {
    res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq),
                     error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0L) next
    dt <- as.data.table(res)
    bn <- basename(f)
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", bn)
    dt[, id_plot := id]; dt[, bit_code := bit]
    hobo_rows[[length(hobo_rows)+1L]] <- dt
  }
}
DT_hobo_floor_aff <- rbindlist(hobo_rows, fill = TRUE)
DT_hobo_floor_aff <- DT_hobo_floor_aff[, .(id_plot, date, bit_code, Delta_Tmax)]
cli_alert("HOBO floor05 affected daily : {nrow(DT_hobo_floor_aff)} rows")

# Merge : drop affected (id_plot × bit_code) from orig, bind floor05
DT_hobo_orig[, key := paste(id_plot, bit_code, sep = "_")]
floor_keys <- DT_hobo_floor_aff[, paste(id_plot, bit_code, sep = "_")]
DT_hobo_orig_keep <- DT_hobo_orig[!(key %in% floor_keys)]
DT_hobo_orig_keep[, key := NULL]
DT_hobo_floor <- rbind(DT_hobo_orig_keep, DT_hobo_floor_aff)
setkey(DT_hobo_floor, id_plot, bit_code, date)
saveRDS(DT_hobo_floor, file.path(OUT_DIR, "DT_daily_HOBO_floor05.rds"))
cli_alert_success("Saved DT_daily_HOBO_floor05.rds ({nrow(DT_hobo_floor)} rows)")

# HOBO contrib wide
DT_agg_hobo <- DT_hobo_floor[, .(Tmax_mean = mean(Delta_Tmax, na.rm = TRUE),
                                   Tmax_P90  = quantile(Delta_Tmax, 0.9,
                                                          na.rm = TRUE, type = 7),
                                   n_days = .N),
                              by = .(id_plot, bit_code)]
DT_contrib_hobo <- .lovb_contrib_core(DT_agg_hobo, id_cols = "id_plot")
saveRDS(DT_contrib_hobo, file.path(OUT_DIR, "DT_contrib_HOBO_floor05.rds"))
cli_alert_success("Saved DT_contrib_HOBO_floor05.rds ({nrow(DT_contrib_hobo)} sensors)")

cli_alert_success("Phase 3 done in {round(as.numeric(difftime(Sys.time(), t_start, units='mins')),1)} min")
