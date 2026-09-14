# =============================================================================
# run_chapter1.R — regenerate the Chapter 1 article analysis and figures from the
# MuSICA outputs that ALREADY EXIST on disk.
#
#   Rscript run_chapter1.R              # dry run: pre-flight + printed plan only
#   CH1_RUN=TRUE Rscript run_chapter1.R # execute stages A + B, then the ledger
#
# Optional switches (all default FALSE / off):
#   CH1_ALIGN_OBS=TRUE   also run c1_align_cluster_obs.R, the ONLY writer of the
#                        canonical observation table tab_hobo_perplot_cluster.csv
#   CH1_RUN_MUSICA=TRUE  also run stage B24, the ONE stage that simulates: the Appendix H
#                        boundary-layer-height insensitivity test (16 short runs, ~4 min)
#   CH1_SYNC_FIGS=TRUE   also copy/rename the figures into the manuscript dir
#   Rscript run_chapter1.R --sync-only   do ONLY the manuscript-dir sync
#
# NO STAGE RUNS MuSICA BY DEFAULT; only B24 does, and only under CH1_RUN_MUSICA=TRUE.
# See RUN_CHAPTER1.md for the full documentation, the
# list of scripts that must never be run, and the known hazards.
# =============================================================================

ROOT <- here::here()
if (is.na(ROOT) || !nzchar(ROOT)) ROOT <- getwd()
setwd(here::here())

ARGS      <- commandArgs(TRUE)
SYNC_ONLY <- "--sync-only" %in% ARGS
RUN       <- identical(toupper(Sys.getenv("CH1_RUN", "FALSE")), "TRUE") || SYNC_ONLY
ALIGN_OBS <- identical(toupper(Sys.getenv("CH1_ALIGN_OBS", "FALSE")), "TRUE")
# Opt-in stages name the environment flag that enables them. Before 2026-07-31 `optional`
# was hard-wired to CH1_ALIGN_OBS, so a second kind of opt-in stage could not be expressed.
RUN_MUSICA <- identical(toupper(Sys.getenv("CH1_RUN_MUSICA", "FALSE")), "TRUE")
gate_on <- function(g) switch(g, CH1_ALIGN_OBS = ALIGN_OBS, CH1_RUN_MUSICA = RUN_MUSICA, TRUE)
SYNC_FIGS <- identical(toupper(Sys.getenv("CH1_SYNC_FIGS", "FALSE")), "TRUE") || SYNC_ONLY

say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
rule <- function(ch = "-") say(strrep(ch, 78))

# ---------------------------------------------------------------------------
# 0. PRE-FLIGHT — the simulation assets that cannot be rebuilt without MuSICA
# ---------------------------------------------------------------------------
# 4th field = what the expected count means. "top" counts entries directly under the path
# (a set of run sub-directories); "nc" counts NetCDFs recursively. Made explicit 2026-07-31:
# the check used a single top-level count for everything, so H2_controlled_topheavy, whose
# 14 runs live under LAI6/ and LAI12/, reported "2/14" and the pre-flight declared itself
# INCOMPLETE while nothing was actually missing.
ASSETS <- list(
  list("forcing (canonical)", "in_files/FR-Blo_2021_v2.nc",                     NA,     NA),
  list("trait-perturbation design NetCDFs", "out_files/Chapter1/nc_sensitivity_perplot_units", 3200L, "top"),
  list("coalition sets (forward)", "out_files/musica_native20_forward",         16L,    "top"),
  list("full-model validation run", "out_files/musica_hobo_native20/1111",      NA,     NA),
  list("footprint-radius sweep", "out_files/radius_test_corr",                  7L,     "top"),
  list("cached Appendix C1 figure (MuSICA-derived, chapter lineage)",
       "out_files/H2_controlled_topheavy_chapter",                              42L,    "nc"),
  list("cLHS design sample (native20)",
       "out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds",       NA,     NA)
)

preflight <- function() {
  say("=== PRE-FLIGHT: simulation assets (NEVER DELETE THESE) ===")
  ok <- TRUE
  for (a in ASSETS) {
    lab <- a[[1]]; path <- a[[2]]; n <- a[[3]]; mode <- a[[4]]
    exists <- file.exists(path)
    detail <- ""
    if (exists && !is.na(n)) {
      got <- if (identical(mode, "nc"))
        length(list.files(path, pattern = "\\.nc$", recursive = TRUE))
      else length(list.files(path))
      detail <- sprintf(" (%d/%d entries)", got, n)
      if (got != n) exists <- FALSE
    }
    ok <- ok && exists
    say("  [%s] %-45s %s%s", if (exists) "x" else " ", lab, path, detail)
  }
  rule()
  ok
}

# ---------------------------------------------------------------------------
# 1. STAGES
# ---------------------------------------------------------------------------
S <- function(stage, label, script, env = character(0), optional = FALSE,
              gate = "CH1_ALIGN_OBS")
  list(stage = stage, label = label, script = script, env = env,
       optional = optional, gate = gate)

STAGES <- list(
  # ---- A. re-extraction from the existing NetCDFs -------------------------
  S("A1", "full-model validation vs 53 loggers (RETIRED -> Ch3)", "_archive_harmonisation_2026-09-14/c1_hobo_validation_to_ch3/c1_hobo_native20_validation_regen.R", optional = TRUE),
  S("A2", "align observed cols of the CANONICAL table", "chapter1/producers/c1_align_cluster_obs.R", optional = TRUE),
  S("A3", "trait-perturbation design (2 metrics x 2 periods)", "chapter1/producers/c1_sensitivity_units_extract2.R"),
  S("A4", "hot-day subset over the coalition sets",     "chapter1/producers/c1_hot_extract.R"),
  S("A5", "footprint-radius sweep, aligned clock",      "chapter1/producers/c1_radius_rescore.R"),
  S("A6", "correction levels (wind / scan angle / k)",  "chapter1/producers/c1_correction_levels_rescore.R"),
  S("A7", "bootstrap CIs on the trait sensitivities",   "chapter1/producers/c1_units_bootstrap_v2.R"),
  S("A8", "Table 1 (per-archetype traits, step/SD, profile swap)", "chapter1/producers/c1_make_table1.R"),
  S("A9", "cover-floor control (Appendix C)",            "chapter1/producers/c1_cover_floor_check.R"),
  S("A10","Table G1 bootstrap CIs (Appendix G)",          "chapter1/producers/c1_attribution_bootstrap_units.R"),
  S("A11","canopy-ratio control (Section 4.3)",           "chapter1/producers/c1_canopy_ratio_check.R"),
  S("A12","profile-centroid stats (Section 4.1, Appendix B)", "chapter1/producers/c1_profile_centroid_stats.R"),
  S("A13","raw scan-angle stats (Appendix F) [needs the external drive]",
          "chapter1/producers/c1_scanangle_raw_stats.R", optional = TRUE),

  # ---- B. figures ---------------------------------------------------------
  S("B1",  "Fig 1  typology",                     "chapter1/producers/c1_fig1_gril_native.R"),
  S("B2",  "Fig 3/4  attribution",                "chapter1/producers/c1_fig3_units.R"),
  S("B3",  "Fig 4/5  operating point",            "chapter1/producers/c1_operating_point_main.R"),
  S("B4",  "Fig 5/6 glass ceiling + Fig J1/I1 forward", "chapter1/producers/c1_fig5_j1_convB.R"),
  S("B5",  "Fig 6/8  observational corroboration","chapter1/producers/c1_fig6_obs_nested_native.R"),
  S("B6",  "Fig F2 / Fig 3  PCA recovering the archetypes", "chapter1/producers/c1_pca_recover_clusters.R"),
  S("B7",  "Fig A1 / Fig 7  footprint radius",    "chapter1/producers/c1_fig7_footprint_radius.R"),
  S("B8",  "Fig B1 / A1  ERA5 station bias",      "chapter1/producers/c1_era5_station_bias.R"),
  S("B9",  "Fig C2 / B2  within-canopy profiles (2 panels)", "scripts/make_vertical_profiles_native.R"),
  S("B10", "Fig C3 / B3  per-trait vertical gradient (per-plot, no baseline)", "scripts/make_trait_vertical_gradient_perplot.R"),
  S("B11", "Fig D1 / C1  residual vs topography", "scripts/make_residual_vs_topo.R"),
  S("B12", "Fig E1 / D1  HOBO coverage",          "chapter1/producers/c1_hobo_coverage_native.R"),
  S("B13", "Fig F3 / E2  observational profile shape", "chapter1/producers/c1_obs_profileshape_test.R"),
  S("B14", "Fig G1 / A2  wind log-profile correction", "chapter1/producers/c1_wind_profile_correction.R"),
  S("B15", "Fig G2 / F1  scan-angle control",     "scripts/fig_scanangle_control_csv.R"),
  S("B16", "Fig G3 / F2  LAI-correction comparison", "scripts/fig_lai_correction_compare.R"),
  S("B17", "Fig S1  hourly temperature scatter",  "chapter1/producers/c1_hourly_scatter.R"),
  S("B18", "Fig S2  attribution, hot days",       "chapter1/producers/c1_fig3_units.R",
           env = c(PERIOD = "hot", METRIC = "dt")),
  S("B19", "Figs S3-S5  hot-day companions",      "chapter1/producers/c1_hot_figures.R"),
  S("B20", "Fig F1 / E1  trait collinearity (native20)", "scripts/make_trait_collinearity_native20.R"),
  S("B21", "Fig C1 / B1  controlled H2 test (chapter lineage; runs MuSICA only if any of the 42 sims is absent)", "scripts/make_h2_controlled_topheavy.R"),
  S("B22", "Fig 2  pipeline schematic (graphviz)", "chapter1/producers/c1_fig2_methodo.R"),
  # Tables H1/H2 stay hand-written, but this stage prints every value they cite from the
  # namelists MuSICA actually stages, so a drift between the appendix and the model is one
  # command away rather than invisible. It runs nothing and writes only its own dump table.
  S("B23", "Tables H1/H2  MuSICA configuration check", "chapter1/producers/c1_dump_musica_config.R"),
  # Fig 1b's four representative plots have no generator; this audits the shipped file
  # against the canonical design sample so it cannot drift from the design it samples.
  S("B23b","Fig 1b  representative plots check",   "chapter1/producers/c1_typical_plots_check.R"),
  # Guards the one metric the whole chapter rests on: no shipped stage may fall back to
  # the legacy extractor, and the headline must stay on the aligned clock.
  S("B23c","ΔTmax convention guard",               "chapter1/producers/c1_check_dtmax_convention.R"),
  # Enforces the documentation contract: header, declared inputs/outputs, roxygen on
  # the sourced library functions. Conventions that are not checked decay silently.
  S("B23d","documentation contract",              "chapter1/producers/c1_check_docs.R"),
  # THE ONE STAGE THAT SIMULATES. Everything else in this runner reads archived outputs;
  # this backs the Appendix H boundary-layer-height claim and needs 16 short MuSICA runs
  # (~4 min). Off unless CH1_RUN_MUSICA=TRUE, so the runner's "no simulation" contract holds
  # by default. It writes only under out_files/Chapter1/nc_blh_insensitivity.
  S("B24", "Appendix H  boundary-layer-height insensitivity [RUNS MuSICA]",
           "chapter1/producers/c1_blh_insensitivity.R", optional = TRUE, gate = "CH1_RUN_MUSICA")
)

# ---------------------------------------------------------------------------
# 2. MANUSCRIPT FIGURE SYNC
#    Recovered by md5-matching every cited figure back to its generating script's
#    output. `short` and `full` are the two manuscript numbering schemes; NA = not
#    cited by that document.
# ---------------------------------------------------------------------------
# CH1_FIGDIR overrides the destination — use it to rehearse the sync into a scratch
# directory before letting it touch the manuscript figures.
FIGDIR <- Sys.getenv("CH1_FIGDIR", "chapter1/manuscript/figures/article_v323")
M <- function(src, short, full) list(src = src, short = short, full = full)

SYNC <- list(
  M("out_files/Chapter1/figures/Fig2_methodo.png",              "Fig2_methodo.png",              "Fig2_methodo.png"),
  M("outputs/figures_chap1/Fig1_typology_gril_native.png",       "Fig1_typology_gril_native.png", "Fig1_typology_gril_native.png"),
  # REMOVED 2026-07-30: this legacy hand-made diagram (2026-07-06) mapped to the SAME
  # destination as the graphviz render above and, being second, silently overwrote it on
  # every run — which is how the corrected Figure 2 kept reverting to "10 m" / "±1 SD".
  M("out_files/Chapter1/figures/Fig3_attribution_units.png",     "Fig3_attribution.png",          "Fig4_attribution.png"),
  M("outputs/figures_chap1/fig_operating_point_main.png",        "Fig4_operating_point.png",      "Fig5_operating_point.png"),
  M("out_files/Chapter1/figures/Fig5_obs_vs_sim_dumbbell.png","Fig5_obs_vs_sim_dumbbell.png","Fig6_obs_vs_sim_dumbbell.png"),
  M("out_files/Chapter1/figures/Fig6_obs_nested_native.png",     "Fig6_obs_corroboration.png",    "Fig7_obs_corroboration.png"),
  M("out_files/Chapter1/figures/FigSh_pca_clusters.png",         "FigF2_pca_recover_clusters_vci.png","Fig3_pca_recover_clusters.png"),
  M("out_files/Chapter1/figures/Fig7_footprint_radius.png",      "FigA1_footprint_radius_sensitivity.png","Fig8_footprint_radius.png"),
  # NB outputs/figures_chap1/FigX_era5_station_bias.png is a byte-identical STALE DUPLICATE
  # under the figure's old name; c1_era5_station_bias.R now writes B1_era5_station_bias.
  M("outputs/figures_chap1/B1_era5_station_bias.png",            "FigB1_era5_station_bias.png",      "FigA1_era5_station_bias.png"),
  M("outputs/figures_pipeline/annex/fig_h2_controlled_topheavy.png","FigC1_h2_controlled_topheavy.png","FigB1_h2_controlled_topheavy.png"),
  # C2/B2 is a MONTAGE, handled separately below.
  M("outputs/figures_chap1/fig_trait_vertical_gradient_perplot.png","FigC3_trait_vertical_gradient.png","FigB3_trait_vertical_gradient.png"),
  M("outputs/figures_pipeline/annex/fig_residual_vs_topo.png",   "FigD1_residual_vs_topo.png",       "FigC1_residual_vs_topo.png"),
  M("out_files/Chapter1/figures/FigSh_plot_trait_coverage_native.png", "FigE1_plot_trait_coverage.png",          "FigD1_plot_trait_coverage.png"),
  # F1/E1 PROVENANCE RESOLVED 2026-08-04. Stage B20 writes this exact path from the canonical
  # native20 sample, and re-running it reproduces the shipped png byte for byte (md5
  # 2faa73a8fd75afa21fe7b20d715c2504, unchanged). The former note, that the figure came from
  # pipeline/12_corrplot_traits.R on the superseded floor05_v2 sample with no reproducing
  # lineage, predates B20 and no longer holds.
  M("out_files/Chapter1/figures/fig_trait_collinearity_native20.png", "FigF1_trait_collinearity.png", "FigE1_trait_collinearity.png"),
  M("outputs/figures_pipeline/annex/fig_obs_profileshape.png",   "FigF3_obs_profileshape.png",       "FigE2_obs_profileshape.png"),
  M("outputs/figures_chap1/G1_wind_profile_correction.png",      "G1_wind_profile_correction.png","FigA2_wind_profile_correction.png"),
  M("outputs/figures_pipeline/annex/fig_scanangle_control.png",  "G2_scanangle_control.png",      "FigF1_scanangle_control.png"),
  M("out_files/Chapter1/figures/F2_lai_correction_compare.png",  "G3_lai_correction_compare.png", "FigF2_lai_correction_compare.png"),
  M("out_files/Chapter1/figures/FigJ1_forward_inclusion.png",    "FigJ1_forward_inclusion.png",   "FigI1_forward_inclusion.png"),
  M("out_files/Chapter1/figures/Fig3_attribution_units_hot.png", "FigS2_attribution_hot.png",      "FigS2_attribution_hot.png"),
  M("outputs/figures_chap1/FigS_obs_vs_sim_hot.png",          "FigS3_obs_vs_sim_hot.png",    "FigS3_obs_vs_sim_hot.png"),
  M("outputs/figures_chap1/FigS_forward_hot.png",                "FigS4_forward_hot.png",          "FigS4_forward_hot.png"),
  M("outputs/figures_chap1/FigS_obs_nested_hot.png",             "FigS5_obs_nested_hot.png",       "FigS5_obs_nested_hot.png"),
  M("outputs/figures_chap1/Fig_hourly_temp_scatter.png",         "FigS1_hourly_temp_scatter.png",  "FigS1_hourly_temp_scatter.png")
)

sync_figs <- function() {
  say("=== SYNC: working dirs -> %s ===", FIGDIR)
  dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
  n_ok <- 0L; n_miss <- 0L
  for (m in SYNC) {
    if (!file.exists(m$src)) { say("  [MISS] %s", m$src); n_miss <- n_miss + 1L; next }
    for (dst in unique(stats::na.omit(c(m$short, m$full)))) {
      file.copy(m$src, file.path(FIGDIR, dst), overwrite = TRUE)
      n_ok <- n_ok + 1L
    }
  }
  # --- C2 / B2 vertical gradient: montage of the two panels from stage B9 ---
  top <- "outputs/figures_chap1/fig_vertical_profiles_native_norm.png"
  bot <- "outputs/figures_chap1/fig_vertical_profiles_native_metres.png"
  cvt <- Sys.which("convert")
  if (file.exists(top) && file.exists(bot) && nzchar(cvt)) {
    for (dst in c("FigC2_vertical_gradient.png", "FigB2_vertical_gradient.png")) {
      system2(cvt, c(shQuote(top), shQuote(bot), "-append", shQuote(file.path(FIGDIR, dst))))
      n_ok <- n_ok + 1L
    }
    say("  [montage] C2/B2_vertical_gradient.png rebuilt from the 2 panels")
  } else {
    say("  [WARN] C2/B2_vertical_gradient.png NOT rebuilt (need both panels + ImageMagick `convert`);")
    say("         the existing file in %s is left untouched.", FIGDIR)
  }
  say("  %d files written, %d sources missing", n_ok, n_miss)
  rule()
}

# ---------------------------------------------------------------------------
# 3. DRIVER
# ---------------------------------------------------------------------------
plan <- function() {
  say("=== PLAN ===")
  st <- ""
  for (s in STAGES) {
    grp <- substr(s$stage, 1, 1)
    if (grp != st) { st <- grp; say("-- stage %s --", grp) }
    skip <- s$optional && !gate_on(s$gate)
    say("  %-4s %-52s %s%s%s", s$stage, s$label, s$script,
        if (length(s$env)) sprintf("  [%s]", paste(names(s$env), s$env, sep = "=", collapse = " ")) else "",
        if (skip) sprintf("   (SKIPPED: set %s=TRUE)", s$gate) else "")
  }
  say("-- stage C --")
  say("  C1   %-52s %s", "sync figures into the manuscript directory",
      if (SYNC_FIGS) "(enabled)" else "(SKIPPED: set CH1_SYNC_FIGS=TRUE)")
  say("  C2   %-52s %s", "figure ledger / presence check", "chapter1/ledger/make_article_figure_set.R")
  rule()
}

run_one <- function(s) {
  say("")
  rule("=")
  say(">>> [%s] %s", s$stage, s$label)
  say("    %s%s", s$script,
      if (length(s$env)) sprintf("   %s", paste(names(s$env), s$env, sep = "=", collapse = " ")) else "")
  rule("=")
  if (!file.exists(s$script)) stop(sprintf("[%s] script not found: %s", s$stage, s$script))
  t0 <- Sys.time()
  code <- system2("Rscript", s$script, env = if (length(s$env))
    paste(names(s$env), s$env, sep = "=") else character(0))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!identical(code, 0L))
    stop(sprintf("[%s] %s FAILED (exit %s) after %.0f s — halted.", s$stage, s$script, code, dt))
  say("<<< [%s] OK in %.0f s", s$stage, dt)
  dt
}

# ---- main -------------------------------------------------------------------
say("")
rule("=")
say("Chapter 1 — analysis + figures from the EXISTING MuSICA outputs (no simulation)")
say("root: %s", ROOT)
rule("=")

ok <- preflight()
if (!SYNC_ONLY) plan()

if (!RUN) {
  say("DRY RUN — nothing executed.")
  if (!ok) say("NOTE: pre-flight is INCOMPLETE; some steps would fail.")
  say("Execute with:  CH1_RUN=TRUE Rscript run_chapter1.R")
  quit(save = "no", status = 0)
}

if (SYNC_ONLY) { sync_figs(); quit(save = "no", status = 0) }
if (!ok) stop("Pre-flight failed: a required simulation asset is missing. See RUN_CHAPTER1.md section 0.")

T0 <- Sys.time()
for (s in STAGES) {
  if (s$optional && !gate_on(s$gate)) { say(""); say(">>> [%s] SKIPPED (%s) — set %s=TRUE to enable", s$stage, s$script, s$gate); next }
  run_one(s)
}
if (SYNC_FIGS) { say(""); sync_figs() }

say("")
rule("=")
say(">>> [C2] figure ledger")
rule("=")
if (!identical(system2("Rscript", "chapter1/ledger/make_article_figure_set.R"), 0L))
  stop("[C2] ledger check FAILED — halted.")

say("")
rule("=")
say("Chapter 1 regenerated in %.1f min.", as.numeric(difftime(Sys.time(), T0, units = "mins")))
say("Ledger must report 0 MISSING. Known hazards: see RUN_CHAPTER1.md section 6.")
rule("=")
