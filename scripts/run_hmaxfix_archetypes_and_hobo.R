# ==============================================================================
# Full re-run with corrected Hmax (buffer-MAX 25 m), current `musica.nml`.
#
# Phase 1 : 4 archetypes × 16 coalitions  =  64 sims
#           → outputs/H1_archetypes_hmaxfix/Arch_C{1..4}/musica_out_*.nc
# Phase 2 : 53 HOBO sensors × 4 new coalitions (NULL + 3 incremental)
#           → out_files/musica_hobo_hmaxfix/{0000,...}/musica_out_HOBO_*.nc
#           The REF coalition (1111) is already produced in REF/ from earlier
#           today.
# Phase 3 : Recompute Δ_v heatmap (LOO × 4 profiles) with corrected Hmax
# Phase 4 : Recompute HOBO forward selection figure with explicit labels
#
# Parallelism : 4 workers via {future} + {furrr}. Each MuSICA call gets its
# own tempdir, so concurrency is safe.
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
source(here::here("R/h1_shapley_archetypes.R"))  # for build_factorial_scenarios_archetypes

N_WORKERS <- 4
plan(multisession, workers = N_WORKERS)
cli_alert("Parallel plan : {N_WORKERS} workers")

# ----------------------------------------------------------------------------
# Phase 1 : Build archetype centroids from the Hmax-fixed cLHS sample
# ----------------------------------------------------------------------------
cli_h1("Phase 1 — Archetype simulations (4 × 16 = 64 sims, parallel)")

df_floor <- as.data.table(readRDS(here::here(
  "out_files/Sensitivity_Analysis/clhs_sample_floor05_hmaxfix.rds")))
df_archetypes <- make_synthetic_archetypes(as.data.frame(df_floor))
cli_alert("Archetype centroids (with corrected Hmax) :")
print(df_archetypes[, c("Cluster", "LAI", "Hmax", "fCover")])

# Build 16 coalitions
fac_scs <- build_factorial_scenarios_archetypes(df_floor)
cli_alert("Scenarios built : {length(fac_scs)}")

# Output root
out_root_arch <- here::here("out_files/H1_archetypes_hmaxfix")
dir.create(out_root_arch, recursive = TRUE, showWarnings = FALSE)

# Build job list (arch, scenario)
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
    jobs[[length(jobs) + 1L]] <- list(arch = arch, sc = sc, out_nc = out_nc, lbl = paste(arch_lbl, sc_nm))
  }
}
cli_alert("Archetype jobs to run : {length(jobs)}")

t_arch <- Sys.time()
res_arch <- future_walk(jobs, function(job) {
  suppressMessages({
    library(here); library(ncdf4); library(sf); library(terra)
    library(tidyverse); library(musica.tools); library(rmusica); library(lubridate)
  })
  source(here::here("R/config.R")); source(here::here("R/io.R"))
  source(here::here("R/musica.R")); source(here::here("R/lad.R"))
  tryCatch(
    run_musica_one(job$arch, job$sc, job$out_nc,
                    CFG$forcing_file, CFG$musica_cmd),
    error = function(e) {
      cat(sprintf("[%s] ERROR : %s\n", job$lbl, e$message))
      NULL
    })
}, .options = furrr_options(seed = TRUE))
cli_alert_success("Archetype sims done in {round(as.numeric(difftime(Sys.time(), t_arch, units='mins')), 1)} min")

# ----------------------------------------------------------------------------
# Phase 2 : Compute archetype LOO Δ_v matrix
# ----------------------------------------------------------------------------
cli_h1("Phase 2 — Compute archetype Δ_v and determine forward-selection order")

df_macro <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
# Coalition map: ARCH_<L>_<H>_<F>_<D> -> 4-bit (LAI|Hmax|fCover|LAD)
LADMAP <- .ARCH_COAL_MAP

extract_one <- function(nc_path) {
  if (!file.exists(nc_path) || file.size(nc_path) < 1e6) return(NA_real_)
  res <- tryCatch(extract_deltatmax_one(nc_path, df_macro, CFG$date_seq,
                                            z_target = CFG$tair_target_height),
                    error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

arch_rows <- list()
for (i in seq_len(nrow(df_archetypes))) {
  arch_lbl <- sprintf("Arch_C%s", df_archetypes$Cluster[i])
  out_dir <- file.path(out_root_arch, arch_lbl)
  for (sc_nm in names(fac_scs)) {
    bit <- LADMAP[sc_nm]
    nc <- file.path(out_dir, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    tair <- extract_one(nc)
    arch_rows[[length(arch_rows) + 1L]] <- data.table(
      Cluster = df_archetypes$Cluster[i],
      bit_code = bit, scenario = sc_nm,
      Tmax_micro = tair)
  }
}
DT_arch <- rbindlist(arch_rows)
DT_arch_w <- dcast(DT_arch, Cluster ~ bit_code, value.var = "Tmax_micro")
# Compute Δ_v = T_max(REF) - T_max(LOO_v)  (negative = v buffers, per user convention)
VARS <- c("LAI","Hmax","fCover","LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
delta <- list()
for (v in VARS) {
  delta[[v]] <- DT_arch_w[["1111"]] - DT_arch_w[[loo_bits[v]]]
}
DT_delta <- data.table(Cluster = DT_arch_w$Cluster, do.call(cbind, delta))
DT_delta_long <- melt(DT_delta, id.vars = "Cluster",
                          variable.name = "variable", value.name = "delta")
DT_avg <- DT_delta_long[, .(Avg = mean(delta)), by = variable]
cli_alert("Avg Δ_v per variable :")
print(DT_avg)

# Forward selection order : signed Avg from positive (anti-buffer) to negative (strong buffer)
setorder(DT_avg, -Avg)
forward_order <- as.character(DT_avg$variable)
cli_alert("Forward order (anti-buffer → strong buffer) : {.field {paste(forward_order, collapse=' → ')}}")

# Save heatmap data
saveRDS(list(arch = DT_arch, delta = DT_delta_long, avg = DT_avg,
              forward_order = forward_order),
         here::here("outputs/hmaxfix_archetype_heatmap_data.rds"))

# ----------------------------------------------------------------------------
# Phase 3 : HOBO sensors — run 4 new coalitions per sensor (NULL + 3 incremental)
# ----------------------------------------------------------------------------
cli_h1("Phase 3 — HOBO forward-selection sims (4 new coalitions × 53 sensors)")

rasters <- load_lidar_rasters(CFG$in_dir, CFG$agg_factor)
df_hobo <- as.data.frame(build_hobo_inputs(
  CFG$hobo_geojson, rasters$stack, CFG$ids_to_remove,
  hmax_raw_raster = file.path(CFG$in_dir, "max_res_10_m.tif"),
  hmax_buffer = 25))
df_hobo[df_hobo$fCover < 0.5, "fCover"] <- 0.5
cli_alert("n HOBO sensors : {nrow(df_hobo)}")

# Compute bit codes for forward steps
# Convention : variable order in bit string = (LAI, Hmax, fCover, LAD) = positions (1,2,3,4)
var_pos <- c(LAI = 1, Hmax = 2, fCover = 3, LAD = 4)
bits_seq <- character(5)
bits_seq[1] <- "0000"     # baseline
bs <- rep("0", 4)
for (i in seq_along(forward_order)) {
  bs[var_pos[forward_order[i]]] <- "1"
  bits_seq[i + 1L] <- paste(bs, collapse = "")
}
# bits_seq[5] should equal "1111" — REF, already done
cli_alert("Forward bit sequence : {paste(bits_seq, collapse=' → ')}")

# Build the scenarios for these bits
fac_scs_hobo <- build_factorial_scenarios_archetypes(df_floor)
# .ARCH_COAL_MAP gives ARCH_xxx → bit, invert
bit2scn <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))

out_root_hobo <- here::here("out_files/musica_hobo_hmaxfix")
hobo_jobs <- list()
for (i in seq_len(nrow(df_hobo))) {
  pr <- df_hobo[i, , drop = FALSE]
  for (k in seq_along(bits_seq)) {
    bit <- bits_seq[k]
    if (bit == "1111") {
      # REF already done in REF/ — skip
      next
    }
    sc_name <- bit2scn[bit]
    if (is.na(sc_name)) next
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
res_hobo <- future_walk(hobo_jobs, function(job) {
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

cli_h1("All sims complete. Run figure regeneration script next.")
cli_alert("Run : Rscript scripts/figs_hmaxfix.R")
