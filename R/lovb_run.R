# ==============================================================================
# LOVB analysis — orchestrator
#
# Single entry point that runs modules 1-5 in sequence :
#   1) extract daily Delta_Tmax (cLHS + archetypes, with cache)
#   2) aggregate to mean + P90 per plot per coalition
#   3) compute REF/LOVB/LVA contributions (wide format)
#   4) generate all 10 figures (5 types x 2 metrics)
#   5) generate all 5 tables (CSV + markdown to console)
#
# Usage (from project root) :
#   Rscript R/lovb_run.R
#
# First run extracts ~4000 cLHS NCs (10-15 min) + 40 archetype NCs (~15s),
# then all downstream is fast (<1 min).
# Subsequent runs use cached .rds in outputs/lovb/data/ -> <2 min total.
# ==============================================================================

suppressMessages({
  library(tidyverse)
  library(data.table)
  library(here)
  library(cli)
  library(ncdf4)
  library(ggplot2)
})

source(here::here("R/lovb_01_load.R"))
source(here::here("R/lovb_02_aggregate.R"))
source(here::here("R/lovb_03_contrib.R"))
source(here::here("R/lovb_04_figures.R"))
source(here::here("R/lovb_05_tables.R"))

t0 <- Sys.time()

# ------------------------------------------------------------------------------
cli_h1("LOVB pipeline — start")
# ------------------------------------------------------------------------------

# ---- 1. Load daily Delta_Tmax (cached if available) -------------------------
cli_h2("[1/5] Load daily Delta_Tmax")
DT_daily_clhs <- lovb_load_daily_clhs(use_cache = TRUE)   # 10-15 min first run
DT_daily_arch <- lovb_load_daily_archetypes(use_cache = TRUE)

# ---- 2. Aggregate -----------------------------------------------------------
cli_h2("[2/5] Temporal aggregation (mean + P90)")
DT_agg_clhs <- lovb_aggregate_clhs(DT_daily_clhs, use_cache = FALSE)
DT_agg_arch <- lovb_aggregate_archetypes(DT_daily_arch, use_cache = FALSE)

# ---- 3. Contributions -------------------------------------------------------
cli_h2("[3/5] Build contributions (REF/LOVB/LVA)")
DT_contrib_clhs <- lovb_contrib_clhs(DT_agg_clhs, use_cache = FALSE)
DT_contrib_arch <- lovb_contrib_archetypes(DT_agg_arch, use_cache = FALSE)

# ---- 4. Figures (10 total : 5 types x 2 metrics) ----------------------------
cli_h2("[4/5] Generate 10 figures")
FIG <- here::here("outputs/lovb/figures")

# Figure 1 — Type A absolu (cLHS, 4 panels x metric)
lovb_fig1_typeA(DT_contrib_clhs, "mean", file.path(FIG, "fig_LOVB_typeA_mean.png"))
lovb_fig1_typeA(DT_contrib_clhs, "P90",  file.path(FIG, "fig_LOVB_typeA_P90.png"))

# Figure 2 — Type B contextuel (cLHS, Delta_v vs trait, LOESS)
lovb_fig2_typeB(DT_contrib_clhs, "mean", file.path(FIG, "fig_LOVB_typeB_mean.png"))
lovb_fig2_typeB(DT_contrib_clhs, "P90",  file.path(FIG, "fig_LOVB_typeB_P90.png"))

# Figure 3 — Lecture D conditionnelle (cLHS, 4x4 grid)
lovb_fig3_lectureD(DT_contrib_clhs, "mean", file.path(FIG, "fig_LOVB_lectureD_mean.png"))
lovb_fig3_lectureD(DT_contrib_clhs, "P90",  file.path(FIG, "fig_LOVB_lectureD_P90.png"))

# Figure 4 — Barplot par archetype x variable
lovb_fig4_archetypes(DT_contrib_arch, "mean", file.path(FIG, "fig_LOVB_archetypes_mean.png"))
lovb_fig4_archetypes(DT_contrib_arch, "P90",  file.path(FIG, "fig_LOVB_archetypes_P90.png"))

# Figure 5 — REF vs LOVB scatter par archetype
lovb_fig5_archetypes_scat(DT_contrib_arch, "mean", file.path(FIG, "fig_LOVB_archetypes_REFvsLOVB_mean.png"))
lovb_fig5_archetypes_scat(DT_contrib_arch, "P90",  file.path(FIG, "fig_LOVB_archetypes_REFvsLOVB_P90.png"))

# ---- 5. Tables -------------------------------------------------------------
cli_h2("[5/5] Generate 5 tables")
lovb_tab1_global(DT_contrib_clhs)
lovb_tab2_per_archetype(DT_contrib_clhs)
lovb_tab3_convergence(DT_contrib_clhs)
lovb_tab4_archetypes(DT_contrib_arch)
lovb_tab5_convergence_archetypes(DT_contrib_arch)

elapsed <- difftime(Sys.time(), t0, units = "mins")
cli_alert_success("LOVB pipeline done in {round(as.numeric(elapsed), 2)} min")
cli_alert_info("Figures : {.path outputs/lovb/figures/}")
cli_alert_info("Tables  : {.path outputs/lovb/tables/}")
cli_alert_info("Cache   : {.path outputs/lovb/data/}")
