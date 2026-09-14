# ==============================================================================
# LOVB extensions — single-command orchestrator
#
# Runs sequentially :
#   1) Module 6 : HOBO LOVB MuSICA sims (6 coalitions x 53 plots)         ~80 min
#   2) Module 8 : baseline Q10/Q90 MuSICA sims (10 coalitions x 52 plots) ~2 h
#   3) Module 7 : HOBO LOVB analysis + convergence with observed          ~30 s
#   4) Module 9 : baseline robustness analysis                            ~10 s
#
# Total : ~3.5 h on single thread (idle the machine).
# Safe to re-run : MuSICA sims skip existing NCs, analyses use cached .rds.
#
# Usage :
#   Rscript R/lovb_run_extensions.R 2>&1 | tee outputs/lovb/run_extensions.log
# ==============================================================================

suppressMessages({
  library(here)
  library(cli)
})

t0 <- Sys.time()

cli_h1("LOVB extensions pipeline — start")

# ---- 1. HOBO MuSICA sims (Module 6) ----------------------------------------
cli_h2("[1/4] HOBO MuSICA runs (~80 min)")
source(here::here("R/lovb_06_hobo_run.R"))
lovb_run_hobo(force = FALSE)
cli_alert_success("[1/4] HOBO runs done at {Sys.time()}")

# ---- 2. Baseline Q10/Q90 MuSICA sims (Module 8) ----------------------------
cli_h2("[2/4] Baseline Q10/Q90 MuSICA runs (~2 h)")
source(here::here("R/lovb_08_baseline_run.R"))
lovb_run_baselines(n_per_cluster = 13L)
cli_alert_success("[2/4] Baseline runs done at {Sys.time()}")

# ---- 3. HOBO analysis (Module 7) -------------------------------------------
cli_h2("[3/4] HOBO LOVB analysis + convergence")
source(here::here("R/lovb_07_hobo_analysis.R"))
lovb_run_hobo_analysis()

# ---- 4. Baseline robustness analysis (Module 9) ----------------------------
cli_h2("[4/4] Baseline robustness analysis")
source(here::here("R/lovb_09_baseline_analysis.R"))
lovb_run_baseline_analysis()

elapsed <- difftime(Sys.time(), t0, units = "mins")
cli_alert_success("LOVB extensions done in {round(as.numeric(elapsed), 2)} min")
cli_alert_info("HOBO data    : {.path outputs/lovb/data/DT_*HOBO*.rds}")
cli_alert_info("Baseline data: {.path outputs/lovb/data/DT_*baselines*.rds}")
cli_alert_info("Tables       : {.path outputs/lovb/tables/}")
cli_alert_info("Figures      : {.path outputs/lovb/figures/}")
