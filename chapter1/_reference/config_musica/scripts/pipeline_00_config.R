# ==============================================================================
# PIPELINE STAGE 00 — Config & conventions (single source of truth)
# Sourced by every stage. Loads the R/ library and freezes the FINAL decisions.
# ==============================================================================

suppressMessages({
  library(data.table); library(here); library(cli); library(ncdf4); library(sf)
  library(tidyverse); library(lubridate)
})
source(here::here("R/config.R"))          # CFG: paths, date_seq, musica_cmd, tair_target_height
source(here::here("R/io.R"))
source(here::here("R/forest.R"))
source(here::here("R/fpca.R"))
source(here::here("R/lad.R"))
source(here::here("R/musica.R"))
source(here::here("R/validation.R"))
source(here::here("R/cluster_relabel.R"))
source(here::here("R/h1_shapley_archetypes.R"))

# ---- branch: "" = legacy z1 (frozen) ; any suffix (z05, z1raw, …) = a LAD-source branch
BRANCH <- Sys.getenv("PIPE_BRANCH", "")
.musica_dir <- if (nzchar(BRANCH)) sprintf("out_files/musica_hobo_%s", BRANCH) else "out_files/musica_hobo_v10_fcovmean"
.fig_root   <- if (nzchar(BRANCH)) sprintf("outputs/figures_pipeline_%s", BRANCH) else "outputs/figures_pipeline"

# ---- FINAL conventions (do not drift — see pipeline/README.md) ---------------
PIPE <- list(
  BRANCH = BRANCH,
  strict = (BRANCH == ""),                          # enforce legacy regression asserts only on legacy branch
  # MuSICA attribution set
  FCOV_BASELINE   = "mean",                         # Eva's convention (v10), NOT 1 (v9)
  MUSICA_DIR_HOBO = here::here(.musica_dir),
  BITS            = c("0000","0001","0010","0011","0100","0101","0110","0111",
                      "1000","1001","1010","1011","1100","1101","1110","1111"),
  REF_BIT = "1111", BASE_BIT = "0000",
  # Per-METRIC extraction conventions (constraint #4 — never unify by accident)
  #   absolute ΔTmax / ΔVPDmax : fixed 1 m interpolation, -2 h shift (HOBO-comparable)
  #   relative micro/macro slope : ALSO fixed 1 m, no time shift (NOT nair==1; the
  #     legacy nair==1 / 45-8 / r=0.92 convention is superseded — all metrics at 1 m)
  Z_FIX = 1.0, TMAX_SHIFT_HR = 2L,
  SLOPE_USE_NAIR1 = FALSE, SLOPE_SHIFT_HR = 0L,   # slope at fixed 1 m (Z_FIX), like ΔTmax/ΔVPDmax
  P_HPA = 1013,                                      # for VPD reconstruction
  HOT_QUANTILE = 0.90,                               # 10% hottest days (macro Tmax)
  # Frozen clustering input (constraint #1 — NEVER silently re-cluster)
  CLUSTER_SAMPLE = here::here("out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds"),
  KMEANS_SEED = 42L,
  # Outputs (branch-specific)
  OUT_FIG = here::here(.fig_root),
  OUT_TAB = here::here(file.path(.fig_root, "tables")),
  OUT_DATA= here::here(file.path(.fig_root, "data")),
  # Regression acceptance targets (constraint #6)
  #   sum_phi_tmax : stage-05 additivity median (legacy HOBO Shapley).
  #   VAL : per-branch stage-07 validation targets. The PUBLISHED branch is z05
  #         (44/9, r=0.93); the legacy v10 branch is the historical guard (45/8).
  ASSERT = list(sum_phi_tmax = -0.311,
                VAL = list(legacy = list(r_pp = 0.92, buf = 45L, amp = 8L),
                           z05    = list(r_pp = 0.93, buf = 44L, amp = 9L)))
)
for (d in c(PIPE$OUT_FIG, PIPE$OUT_TAB, PIPE$OUT_DATA))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- small helpers -----------------------------------------------------------
esat_hpa <- function(Tc) 6.108 * exp(17.27 * Tc / (Tc + 237.3))
pipe_assert <- function(cond, msg) {
  if (isTRUE(cond)) cli::cli_alert_success(msg) else cli::cli_abort(paste("ASSERT FAILED:", msg))
}
pipe_near <- function(a, b, tol = 0.01) abs(a - b) <= tol

# binary sanity (legacy binary, md5 5307...; NOT model-3.2.3)
if (!file.exists(CFG$musica_cmd)) cli::cli_alert_warning("musica_cmd missing: {CFG$musica_cmd}")

cli::cli_alert_info("Pipeline config loaded — fCover baseline = {PIPE$FCOV_BASELINE}; ALL metrics at fixed 1 m (ΔTmax/ΔVPDmax -{PIPE$TMAX_SHIFT_HR}h shift, slope no shift).")
