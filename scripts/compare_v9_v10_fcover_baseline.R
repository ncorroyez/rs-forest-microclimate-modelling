# ==============================================================================
# Compare V9 (fcov_b = 1) vs V10 (fcov_b = mean = 0.869) :
#   - Heatmap LOO Δ (archetypes)
#   - Inverse forward selection on HOBO
# Outputs side-by-side stats so we can judge whether Eva's suggested baseline
# changes the conclusions of the talk.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/h1_shapley_archetypes.R"))

OUT <- here::here("outputs/figs_MEB2026_final")
df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

# Returns mean Tmax_micro per (archetype, bit). For bits where fCover=0 (i.e.
# baseline-dependent), reads from the V10 directory; otherwise from V9.
archetype_tmax <- function(cl, bit, use_v10 = TRUE) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  sc_nm   <- setNames(names(.ARCH_COAL_MAP), unname(.ARCH_COAL_MAP))[bit]
  v10_path <- file.path("out_files/H1_archetypes_v10_fcovmean",
                         arch_lbl, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
  v9_path  <- file.path("out_files/H1_archetypes_v9",
                         arch_lbl, sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
  path <- if (use_v10 && file.exists(v10_path)) v10_path else v9_path
  res <- tryCatch(extract_deltatmax_one(path, df_macro, CFG$date_seq,
                                         z_target = CFG$tair_target_height),
                   error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0) return(NA_real_)
  mean(res$Tmax_micro, na.rm = TRUE)
}

# Build wide tables V9 and V10
build_wide <- function(use_v10) {
  rows <- list()
  for (cl in 1:4) {
    for (bit in unname(.ARCH_COAL_MAP)) {
      rows[[length(rows) + 1L]] <- data.table(
        Cluster = cl, bit = bit,
        Tmax = archetype_tmax(cl, bit, use_v10 = use_v10))
    }
  }
  dcast(rbindlist(rows), Cluster ~ bit, value.var = "Tmax")
}

cli_h1("Archetype Tmax means : V9 (fcov_b=1) vs V10 (fcov_b=mean=0.869)")
V9  <- build_wide(use_v10 = FALSE)
V10 <- build_wide(use_v10 = TRUE)
cat("\nV9 (fcov_b=1) wide :\n"); print(V9)
cat("\nV10 (fcov_b=mean) wide :\n"); print(V10)

# LOO Δ on Tmax
VARS <- c("LAI", "Hmax", "fCover", "LAD")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")

loo_delta <- function(W) {
  out <- list()
  for (v in VARS) out[[v]] <- W[["1111"]] - W[[loo_bits[v]]]
  data.table(Cluster = W$Cluster, do.call(cbind, out))
}
D_V9  <- loo_delta(V9)
D_V10 <- loo_delta(V10)
cat("\nLOO Δ_v V9  (°C) :\n");  print(D_V9)
cat("\nLOO Δ_v V10 (°C) :\n");  print(D_V10)

cat("\nAvg Δ_v across archetypes :\n")
cat("V9  :"); print(sapply(D_V9[, !"Cluster"], mean))
cat("V10 :"); print(sapply(D_V10[, !"Cluster"], mean))

# ---------- Inverse forward HOBO on V10 ------------------------
cli_h1("Inverse forward HOBO — V10 (fcov_b=mean)")

collect_step_v10 <- function(bit, metric = c("slope", "Tmax")) {
  metric <- match.arg(metric)
  subdir_v10 <- file.path("out_files/musica_hobo_v10_fcovmean", bit)
  subdir_v9  <- file.path("out_files/musica_hobo_v9", bit)
  # Use V10 if available (baseline-dependent), else V9
  subdir <- if (dir.exists(subdir_v10) && length(list.files(subdir_v10, "\\.nc$")) > 0) subdir_v10 else subdir_v9
  cat(sprintf("  bit %s -> %s\n", bit, subdir))
  nc_files <- list.files(subdir, pattern = "\\.nc$", full.names = TRUE)
  rows <- list()
  for (f in nc_files) {
    id <- sub("musica_out_HOBO_(.+)\\.nc$", "\\1", basename(f))
    if (metric == "slope") {
      res <- tryCatch(extract_hourly_slope_one(f, era5_hourly, CFG$date_seq,
                                                z_target = CFG$tair_target_height),
                       error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id, sim_metric = res$slope)
    } else {
      res <- tryCatch(extract_deltatmax_one(f, df_macro, CFG$date_seq,
                                             z_target = CFG$tair_target_height),
                       error = function(e) NULL)
      if (is.null(res) || nrow(res) == 0) next
      rows[[id]] <- data.table(id_plot = id, sim_metric = mean(res$Delta_Tmax, na.rm = TRUE))
    }
  }
  rbindlist(rows)
}

# HOBO observed (unchanged between V9 and V10)
hobo_raw <- as.data.table(read.csv(CFG$hobo_temp_csv)) %>%
  mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  filter(position_sensor == "a",
         as.Date(datetime) %in% CFG$date_seq,
         !id_plot %in% CFG$ids_to_remove) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(id_plot, time) %>%
  summarise(t_hobo = mean(t_hobo, na.rm = TRUE), .groups = "drop") %>%
  as.data.table()
hobo_with_era5 <- merge(hobo_raw, era5_hourly, by = "time")
hobo_slope_obs <- hobo_with_era5[, {
  fit <- lm(t_hobo ~ Tair_era5)
  .(slope_obs = as.numeric(coef(fit)[2]))
}, by = id_plot]

paths <- list(
  inv_slope = list(
    bits   = c("0000","0001","0101","1101","1111"),
    labels = c("Baseline","+LAD","+LAD+Hmax","+LAD+Hmax+LAI","REF (+fCover)"),
    metric = "slope", obs_col = "slope_obs")
)

make_path_fits <- function(p_cfg) {
  long <- list(); fits <- list()
  for (k in seq_along(p_cfg$bits)) {
    bit <- p_cfg$bits[k]
    d <- collect_step_v10(bit, p_cfg$metric)
    if (nrow(d) == 0) next
    d <- merge(d, hobo_slope_obs[, .(id_plot, obs_metric = slope_obs)], by = "id_plot")
    d[, step := k]; d[, label := p_cfg$labels[k]]; d[, bit := bit]
    long[[k]] <- d
    ok <- !is.na(d$sim_metric) & !is.na(d$obs_metric)
    x <- d$obs_metric[ok]; y <- d$sim_metric[ok]
    fits[[k]] <- data.table(step = k, label = p_cfg$labels[k],
                             r = suppressWarnings(cor(x, y)),
                             RMSE = sqrt(mean((y - x)^2)),
                             MAE  = mean(abs(y - x)), n = length(x))
  }
  list(long = rbindlist(long, fill = TRUE), fits = rbindlist(fits))
}

cli_h2("V10 inverse slope path :")
res_v10 <- make_path_fits(paths$inv_slope)
print(res_v10$fits)

# Compare with V9 fits (cached)
v9_fits <- fread(file.path(OUT, "tab_v9_HOBO_inverse_slope.csv"))
cli_h2("V9 inverse slope (for comparison) :")
print(v9_fits)

cat("\n=== V9 vs V10 inverse slope (r per step) ===\n")
cmp <- merge(v9_fits[, .(step, V9_r = r, V9_RMSE = RMSE)],
              res_v10$fits[, .(step, V10_r = r, V10_RMSE = RMSE)], by = "step")
cmp[, dr_v10_minus_v9 := V10_r - V9_r]
cmp[, dRMSE_v10_minus_v9 := V10_RMSE - V9_RMSE]
print(cmp)

fwrite(cmp, file.path(OUT, "tab_v9_vs_v10_fcover_baseline_comparison.csv"))
saveRDS(list(V9_arch = V9, V10_arch = V10, V9_LOO = D_V9, V10_LOO = D_V10,
              v9_HOBO_fits = v9_fits, v10_HOBO_fits = res_v10$fits, cmp = cmp),
         here::here("outputs/v10_fcover_baseline_comparison.rds"))
cli_h1("Comparison complete — see tab_v9_vs_v10_fcover_baseline_comparison.csv")
