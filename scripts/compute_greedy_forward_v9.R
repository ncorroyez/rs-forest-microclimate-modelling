# ==============================================================================
# True greedy forward feature selection from the V9 archetype factorial.
#
# At each step k in {1..4}, the variable v* added is the one whose inclusion
# minimises the mean RMSE (vs REF=1111) across the 4 archetypes — i.e. the
# variable that contributes the most to reaching REF.
#
# Compares with the LOO-ordered forward (current figure) and reports.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4)
  library(tidyverse); library(lubridate); library(musica.tools)
})

source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/musica.R"))
source(here::here("R/h1_shapley_archetypes.R"))

cli_h1("True greedy forward selection (V9 archetype factorial)")

df_macro    <- extract_macro_daily(CFG$forcing_file, CFG$date_seq)
era5_hourly <- build_era5_hourly(CFG$forcing_file, CFG$date_seq)

arch_root   <- here::here("out_files/H1_archetypes_v9")
clusters    <- 1:4
VARS        <- c("LAI", "Hmax", "fCover", "LAD")

# ----------------------------------------------------------------------------
# 1. Extract two metrics per (archetype, coalition) :
#    - mean Tmax_micro (daily)
#    - slope (hourly Tair_sim ~ Tair_era5)
# ----------------------------------------------------------------------------
rows <- list()
for (cl in clusters) {
  arch_lbl <- sprintf("Arch_C%s", cl)
  for (sc_nm in names(.ARCH_COAL_MAP)) {
    bit <- .ARCH_COAL_MAP[sc_nm]
    nc  <- file.path(arch_root, arch_lbl,
                      sprintf("musica_out_%s_%s.nc", arch_lbl, sc_nm))
    tmax_res  <- tryCatch(extract_deltatmax_one(nc, df_macro, CFG$date_seq,
                                                z_target = CFG$tair_target_height),
                          error = function(e) NULL)
    slope_res <- tryCatch(extract_hourly_slope_one(nc, era5_hourly, CFG$date_seq,
                                                   z_target = CFG$tair_target_height),
                          error = function(e) NULL)
    rows[[length(rows) + 1L]] <- data.table(
      Cluster   = cl, bit = bit, sc_nm = sc_nm,
      Tmax_mean = if (is.null(tmax_res)) NA_real_ else mean(tmax_res$Tmax_micro, na.rm = TRUE),
      slope     = if (is.null(slope_res)) NA_real_ else slope_res$slope)
  }
}
DT <- rbindlist(rows)
cli_alert("Extracted {nrow(DT)} (archetype × coalition) rows")

# ----------------------------------------------------------------------------
# 2. Decode bit codes (LAI|Hmax|fCover|LAD)
# ----------------------------------------------------------------------------
DT[, LAI    := as.integer(substr(bit, 1, 1))]
DT[, Hmax   := as.integer(substr(bit, 2, 2))]
DT[, fCover := as.integer(substr(bit, 3, 3))]
DT[, LAD    := as.integer(substr(bit, 4, 4))]

# ----------------------------------------------------------------------------
# 3. Greedy forward selection — criterion : minimise mean |Tmax_mean - REF|
# ----------------------------------------------------------------------------
REF_bit  <- "1111"
DT_ref   <- DT[bit == REF_bit, .(Cluster, Tmax_ref = Tmax_mean, slope_ref = slope)]
DT_full  <- merge(DT, DT_ref, by = "Cluster")
DT_full[, dTmax2  := (Tmax_mean - Tmax_ref)^2]
DT_full[, dSlope2 := (slope - slope_ref)^2]

# Helper : compute mean RMSE / slope-MSE across archetypes for a given bit pattern
bit_of_set <- function(set) {
  on_off <- ifelse(VARS %in% set, "1", "0")
  paste(on_off, collapse = "")
}
score_bit <- function(b, metric = "Tmax") {
  d <- DT_full[bit == b]
  if (nrow(d) == 0) return(NA_real_)
  if (metric == "Tmax") sqrt(mean(d$dTmax2, na.rm = TRUE))
  else                   sqrt(mean(d$dSlope2, na.rm = TRUE))
}

greedy_path <- function(metric = "Tmax") {
  active   <- character(0)
  remain   <- VARS
  path_log <- list()
  for (k in 0:length(VARS)) {
    if (length(remain) == 0) break
    candidates <- lapply(remain, function(v) c(active, v))
    scores <- sapply(candidates, function(s) score_bit(bit_of_set(s), metric))
    best_i <- which.min(scores)
    best_v <- remain[best_i]
    active <- c(active, best_v)
    remain <- setdiff(VARS, active)
    path_log[[length(path_log) + 1L]] <- data.table(
      step = k + 1L, added = best_v,
      cumulative = paste(active, collapse = "+"),
      bit = bit_of_set(active),
      rmse_vs_REF = scores[best_i])
  }
  rbindlist(path_log)
}

cli_h2("Greedy forward — criterion : minimise RMSE(Tmax_mean vs REF)")
path_tmax <- greedy_path("Tmax")
print(path_tmax)

cli_h2("Greedy forward — criterion : minimise RMSE(slope vs REF)")
path_slope <- greedy_path("slope")
print(path_slope)

# ----------------------------------------------------------------------------
# 4. Backward elimination (LOO) — current "forward order" approach
# ----------------------------------------------------------------------------
cli_h2("Reference : LOO Δ-ordered forward (current V9 figure)")
loo_bits <- c(LAI = "0111", Hmax = "1011", fCover = "1101", LAD = "1110")
DT_ref_w <- dcast(DT[, .(Cluster, bit, Tmax_mean)], Cluster ~ bit, value.var = "Tmax_mean")
loo_delta <- sapply(VARS, function(v)
  mean(DT_ref_w[["1111"]] - DT_ref_w[[loo_bits[v]]], na.rm = TRUE))
cat("\nLOO Δ_v on Tmax (REF - LOO_v, °C) :\n"); print(loo_delta)
loo_order_tmax <- names(sort(abs(loo_delta), decreasing = TRUE))
cat("LOO |Δ| order (largest first) :", paste(loo_order_tmax, collapse = " → "), "\n")

state <- readRDS(here::here("outputs/v9_pipeline_state.rds"))
cat("Slope-based forward order (current figure)  :",
    paste(state$forward_order, collapse = " → "), "\n")
cat("Greedy (Tmax)                              :",
    paste(path_tmax$added, collapse = " → "), "\n")
cat("Greedy (slope)                             :",
    paste(path_slope$added, collapse = " → "), "\n")

# ----------------------------------------------------------------------------
# 5. Save
# ----------------------------------------------------------------------------
saveRDS(list(DT_factorial = DT, greedy_tmax = path_tmax, greedy_slope = path_slope,
              loo_delta = loo_delta),
         here::here("outputs/v9_greedy_forward.rds"))
fwrite(path_tmax,  here::here("outputs/figs_MEB2026_final/tab_v9_greedy_forward_Tmax.csv"))
fwrite(path_slope, here::here("outputs/figs_MEB2026_final/tab_v9_greedy_forward_slope.csv"))

cli_h1("Greedy forward complete.")
