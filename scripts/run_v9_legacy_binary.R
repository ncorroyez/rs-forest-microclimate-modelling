# ==============================================================================
# V9 pipeline — re-run all sims with the TRUE LEGACY MuSICA binary
# (/home/corroyez/Documents/musica/musica, md5 53072278, Nov 2024) AND switch
# the buf/amp / forward-selection metric to HOURLY slope (Tair_sim ~ ERA5).
#
#   - Phase 1 : 4 archetypes × 16 coalitions = 64 sims (out_files/H1_archetypes_v9/)
#   - Phase 2 : determine forward-selection order from archetype LOO
#               (slope-based metric, log_slope < 0 = buffering)
#   - Phase 3 : 53 HOBO sensors × 5 forward coalitions = 265 sims
#               (out_files/musica_hobo_v9/)
#   - Phase 4 : regenerate LOO heatmap + HOBO forward-selection figure
#               using the hourly-slope metric
#
# CFG$musica_cmd has already been updated to point to the legacy binary.
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

stopifnot(file.exists(CFG$musica_cmd))
cli_h1("V9 pipeline — legacy binary + hourly-slope metric")
cli_alert("MuSICA binary : {CFG$musica_cmd}")
cli_alert("md5           : {tools::md5sum(CFG$musica_cmd)}")

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)

# ============================================================================
# Phase 1 : Archetype sims (4 × 16 = 64)
# ============================================================================
cli_h1("Phase 1 — Archetype sims")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds")))
df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
cli_alert("Archetypes :")
print(df_archetypes[, c("Cluster", "LAI", "Hmax", "fCover")])

fac_scs <- build_factorial_scenarios_archetypes(df_floor)
out_root_arch <- here::here("out_files/H1_archetypes_v9")
dir.create(out_root_arch, recursive = TRUE, showWarnings = FALSE)

jobs <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch <- df_archetypes[i, , drop = FALSE]
  arch_lbl <- sprintf("Arch_C%s", arch$Cluster)
  out_dir <- file.path(out_root_arch, arch_lbl)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (sc_nm in names(fac_scs)) {
    sc <- fac_scs[[sc_nm]]
    out_nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    if (file.exists(out_nc) && file.size(out_nc) > 1e6) next
    jobs[[length(jobs) + 1L]] <- list(arch = arch, sc = sc, out_nc = out_nc,
                                      lbl = paste(arch_lbl, sc_nm))
  }
}
cli_alert("Archetype jobs to run : {length(jobs)}")

t_arch <- Sys.time()
future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$arch, job$sc, job$out_nc,
                   CFG$forcing_file, CFG$musica_cmd),
    error = function(e) cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message)))
}, .options = furrr_options(seed = TRUE))
cli_alert_success("Archetype sims done in {round(as.numeric(difftime(Sys.time(), t_arch, units='mins')), 1)} min")

# ============================================================================
# Phase 2 : Forward-selection order from archetype LOO — HOURLY-SLOPE metric
# ============================================================================
cli_h1("Phase 2 — Forward order (hourly-slope LOO Δ)")

era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)
df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
LADMAP <- .ARCH_COAL_MAP

extract_slope <- function(nc_path) {
  if (!file.exists(nc_path) || file.size(nc_path) < 1e6) return(NA_real_)
  res <- tryCatch(extract_hourly_slope_one(nc_path, era5_hourly, CFG$date_seq,
                                           z_target = CFG$tair_target_height),
                  error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  as.numeric(res$slope)
}

arch_rows <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch_lbl <- sprintf("Arch_C%s", df_archetypes$Cluster[i])
  for (sc_nm in names(fac_scs)) {
    bit <- LADMAP[sc_nm]
    nc <- file.path(out_root_arch, arch_lbl,
                    sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    sl <- extract_slope(nc)
    arch_rows[[length(arch_rows) + 1L]] <- data.table(
      Cluster = df_archetypes$Cluster[i],
      bit_code = bit, slope = sl)
  }
}
DT_arch <- rbindlist(arch_rows)
DT_arch_w <- dcast(DT_arch, Cluster ~ bit_code, value.var = "slope")

# LOO Δ on hourly slope : slope(REF) - slope(LOO_v)
# negative Δ_v ⇒ removing v reduces buffering (i.e. v contributes to buffering)
VARS <- c("LAI","Hmax","fCover","LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
delta <- list()
for (v in VARS) delta[[v]] <- DT_arch_w[["1111"]] - DT_arch_w[[loo_bits[v]]]
DT_delta <- data.table(Cluster = DT_arch_w$Cluster, do.call(cbind, delta))
DT_delta_long <- melt(DT_delta, id.vars = "Cluster",
                      variable.name = "variable", value.name = "delta")
DT_avg <- DT_delta_long[, .(Avg = mean(delta, na.rm = TRUE)), by = variable]
cli_alert("Avg Δ_v (REF - LOO_v on hourly slope) per variable :")
print(DT_avg)

# Variables that strongly buffer have NEGATIVE Δ_v (removing them raises slope toward 1).
# Forward order should add the strongest buffer first → sort ASCENDING.
setorder(DT_avg, Avg)
forward_order <- as.character(DT_avg$variable)
cli_alert("Forward order (strongest buffer → weakest) : {paste(forward_order, collapse=' → ')}")

saveRDS(list(arch = DT_arch, delta = DT_delta_long, avg = DT_avg,
             forward_order = forward_order),
        here::here("outputs/v9_archetype_heatmap_data.rds"))

# ============================================================================
# Phase 3 : HOBO sims (53 × 5 coalitions)
# ============================================================================
cli_h1("Phase 3 — HOBO forward-selection sims")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
cli_alert("HOBO sensors : {nrow(df_hobo)}")
cli_alert("Hmax mean={round(mean(df_hobo$Hmax),2)} m, LAI mean={round(mean(df_hobo$LAI),2)}")

var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}
cli_alert("Forward bit sequence : {paste(bits_seq, collapse=' → ')}")

bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))
fac_scs_hobo <- build_factorial_scenarios_archetypes(df_floor)

out_root_hobo <- here::here("out_files/musica_hobo_v9")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (k in seq_along(bits_seq)) {
    bit <- bits_seq[k]
    sc_name <- bit2scn[bit]; if (is.na(sc_name)) next
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
cli_alert("HOBO jobs to run : {length(hobo_jobs)}")

t_hobo <- Sys.time()
future_walk(hobo_jobs, function(job) {
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
cli_alert_success("HOBO sims done in {round(as.numeric(difftime(Sys.time(), t_hobo, units='mins')), 1)} min")

cli_h1("V9 sims complete — see scripts/make_figs_v9_hourly.R for figures.")
saveRDS(list(forward_order = forward_order, bits_seq = bits_seq,
             arch_slopes = DT_arch_w),
        here::here("outputs/v9_pipeline_state.rds"))
