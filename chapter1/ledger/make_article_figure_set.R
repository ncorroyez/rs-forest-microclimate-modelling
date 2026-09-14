# ==============================================================================
# Chapter 1 figure ledger — provenance + presence checker
# ALIGNED ON THE LIVE MANUSCRIPT:
#   chapter1/manuscript/manuscript_chap1_FINAL_coherence_2026-09-09.md
#
# Rewritten 2026-09-14 (harmonisation, PRODUCER_FIGURE_MAP_2026-09-14.md). The previous
# ledger described a SUPERSEDED narrative (Fig3=attribution, Fig4=operating_point,
# Fig5=dumbbell, Fig6=obs_corroboration, plus glass-ceiling / forward-inclusion companions)
# tied to manuscript_chap1_EN_native20.md and a _SHORT.md with two numbering schemes. The
# live manuscript was restructured: it introduced the CHS41 hand-run chain (Fig 4-7, B1, S1)
# and DROPPED the logger-validation figures (validation moved to Chapter 3). Those dead
# entries are removed here.
#
# Live lineage (do NOT carry the old header values forward):
#   forcing        = out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc
#                    (CHS41 station + ERA5 + MERRA-2 h_sbl; NOT FR-Blo_2021_v2.nc, whose
#                     boundary layer is ~10x too low — see forcing_hsbl_trap)
#   baseline wind  = NO baseline wind correction (Appendix A reports the sensitivity only)
#   attribution    = on-manifold trait perturbation in fixed native steps (no Shapley, no per-SD)
#   ΔTmax          = CHS41 station 1.5 m, time-matched; MuSICA v3.2.3 iterative ABL ("yoyo"); k=0.5
#   design sample  = out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds (400 real 20 m pixels)
#
# This is a LEDGER, not a copier: for every figure the manuscript embeds it records the
# generating script, whether a sim re-run refreshes it, and checks the file is present. It
# writes tab_article_figures.csv + a MISSING list. Reorder / relabel by editing the lists.
#   Rscript chapter1/ledger/make_article_figure_set.R
#
# The manuscript embeds figures from TWO directories:
#   chapter1/manuscript/figures/article_v323/   (Fig1, Fig2, Fig3, A1, A2, C1, D1)
#   chapter1/manuscript/figures/                (Fig4-8, B1, S1 — the CHS41 chain)
# FIGDIR is the parent; the article_v323/ ones carry that prefix in the filename.
#
# KNOWN HAZARDS (see also RUN_CHAPTER1.md §6):
#   * TWO Fig6 files sit in the working dir out_files/Chapter1/figures/:
#       Fig6_h2_controlled_chs41.png  (c1_fig6_fig7_chs41.R) — NOT the embedded panel
#       Fig6_h2_gaussian_chs41.png    (c1_fig6v2_gaussian_heatmap.R) — the embedded one
#     A hand-copy must ship the gaussian panel, not the controlled one.
#   * Fig B1 and Fig S1 originally had NO plotting script: the shipped pngs were plotted
#     interactively on 2026-09-08 and the code was lost (searched all .R incl. archives + base
#     png() devices, all data readers, .Rhistory, RStudio history). RE-AUTHORED 2026-09-14 as
#     c1_b1_pervariable_plot_chs41.R (mirrors Fig7 on nc_b3_chs41/) and
#     c1_s1_attribution_hot_plot_chs41.R (mirrors Fig4 on hotdays_chs41.csv). Both reproduce the
#     shipped canvas exactly (3900x1200, 3900x1080) and the manuscript's numbers; the aesthetic
#     is a faithful mirror of the siblings, NOT byte-identical to the interactive originals, so
#     the shipped pngs are left in place until validated. FIG is env-overridable
#     (CH1_CHS41_FIGDIR) to rehearse into a scratch dir without touching the shipped files.
#   * The CHS41 chain (Fig 4-7, B1, S1) is NOT orchestrated by run_chapter1.R, and
#     run_chapter1.R --sync-only reaches only article_v323/, not figures/. Reproducing these
#     figures is manual. See PRODUCER_FIGURE_MAP_2026-09-14.md.
# ==============================================================================
suppressMessages({ library(here) })
FIGDIR <- here::here("chapter1/manuscript/figures")   # parent; article_v323/ figs carry the prefix

# list(id, filename (relative to FIGDIR), NA (kept for check_set arity), generating_script, refresh_on_sim_rerun)
# refresh=TRUE  -> depends on the native20 / CHS41 sims (values move on a re-run)
# refresh=FALSE -> structural / descriptive / diagnostic (independent of the sim values)
# Script labels are REPO-RELATIVE paths (from the project root) so the existence check
# below can resolve them with here::here(). Since the 2026-09-14 reorg every c1_* producer
# lives in chapter1/producers/; the two non-c1_ appendix producers stayed in scripts/.
# A label may chain two steps with "&&" or carry a "VAR=value" prefix; every ".R" token is
# checked. Free text in parentheses is ignored by the check.
P <- "chapter1/producers"
MAIN <- list(
  list("Fig1", "article_v323/Fig1_typology_gril_native.png", NA, file.path(P, "c1_fig1_gril_native.R"),                 FALSE),
  # graphviz diagram, source of truth scripts/c1_fig2_methodo.dot
  list("Fig2", "article_v323/Fig2_methodo.png",              NA, paste(file.path(P, "c1_fig2_methodo.R"), "(graphviz)"), FALSE),
  # seeded since 2026-07-29 (reads the native20 sample); was non-deterministic before
  list("Fig3", "article_v323/Fig3_pca_recover_clusters.png", NA, file.path(P, "c1_pca_recover_clusters.R"),             TRUE),
  # CHS41 chain: writes to out_files/Chapter1/figures/, hand-copied to figures/
  list("Fig4", "Fig4_attribution_chs41.png",                 NA, file.path(P, "c1_fig4_fig5_chs41.R"),                  TRUE),
  list("Fig5", "Fig5_operating_point_chs41.png",             NA, file.path(P, "c1_fig4_fig5_chs41.R"),                  TRUE),
  # HAZARD: c1_fig6_fig7_chs41.R writes a DIFFERENT Fig6 (h2_controlled) into the same dir
  list("Fig6", "Fig6_h2_gaussian_chs41.png",                 NA, file.path(P, "c1_fig6v2_gaussian_heatmap.R"),          TRUE),
  list("Fig7", "Fig7_vertical_Tprofile_chs41.png",           NA, file.path(P, "c1_fig6_fig7_chs41.R"),                  TRUE),
  list("Fig8", "Fig8_footprint_radius.png",                  NA,
       paste(file.path(P, "c1_radius_rescore.R"), "&&", file.path(P, "c1_fig7_footprint_radius.R")), TRUE)
)
APPENDIX <- list(
  list("A1", "article_v323/FigA1_era5_station_bias.png",       NA, file.path(P, "c1_era5_station_bias.R"),        FALSE),
  list("A2", "article_v323/FigA2_wind_profile_correction.png", NA, file.path(P, "c1_wind_profile_correction.R"),  FALSE),
  # plot re-authored 2026-09-14 (mirrors Fig7); sim step = c1_b3_pervariable_chs41.R (nc_b3_chs41/)
  list("B1", "FigB1_pervariable_gradient_chs41.png",          NA, file.path(P, "c1_b1_pervariable_plot_chs41.R"), TRUE),
  # five panels from the canonical native20 sample; reproduces the shipped png (RESOLVED 2026-08-04)
  list("C1", "article_v323/FigC1_trait_collinearity.png",     NA, "scripts/make_trait_collinearity_native20.R",   FALSE),
  list("D1", "article_v323/FigD1_scanangle_control.png",      NA, "scripts/fig_scanangle_control_csv.R",          FALSE)
)
SUPP <- list(
  # plot re-authored 2026-09-14 (mirrors Fig4); extract step = c1_hotdays_chs41.R (hotdays_chs41.csv)
  list("S1", "FigS1_attribution_hot_chs41.png",              NA, file.path(P, "c1_s1_attribution_hot_plot_chs41.R"), TRUE)
)
# NB tables (H1/H2 numbers, config) are inline markdown in the manuscript, not figures.

check_set <- function(items, kind) {
  do.call(rbind, lapply(items, function(it) {
    id <- it[[1]]; fn_s <- it[[2]]; fn_f <- it[[3]]; script <- it[[4]]; refresh <- it[[5]]
    do.call(rbind, lapply(list(c("live", fn_s), c("full", fn_f)), function(z) {
      if (is.na(z[2])) return(NULL)
      png <- file.path(FIGDIR, z[2]); ok <- file.exists(png)
      data.frame(kind = kind, doc = z[1], id = id, filename = z[2], script = script,
                 refresh_rerun = refresh, status = if (ok) "present" else "MISSING",
                 mtime = if (ok) format(file.info(png)$mtime, "%Y-%m-%d %H:%M") else NA_character_)
    }))
  }))
}
man <- rbind(check_set(MAIN, "main"), check_set(APPENDIX, "appendix"), check_set(SUPP, "supp"))
write.csv(man, file.path(FIGDIR, "tab_article_figures.csv"), row.names = FALSE)

cat(sprintf("\n=== Chapter 1 figure ledger (%d embedded figures, live manuscript) ===\n", nrow(man)))
print(man[, c("kind","id","filename","script","refresh_rerun","status","mtime")], row.names = FALSE)
miss <- man[man$status == "MISSING", ]
if (nrow(miss)) cat(sprintf("\n%d MISSING (regenerate):\n  %s\n", nrow(miss),
    paste(sprintf("%s = %s  [%s]", miss$id, miss$filename, miss$script), collapse = "\n  ")))
cat(sprintf("\n%d / %d entries refreshed by a sim re-run.\n", sum(man$refresh_rerun), nrow(man)))

# --- producer script existence check (distinct from the png presence check above) ------
# Split a script label into its repo-relative script paths: cut on "&&", drop leading
# "VAR=value" env prefixes, keep the first token of each command if it ends in ".R".
extract_script_paths <- function(label) {
  cmds <- trimws(strsplit(label, "&&", fixed = TRUE)[[1]])
  unlist(lapply(cmds, function(cmd) {
    toks <- strsplit(cmd, "[[:space:]]+")[[1]]
    toks <- toks[!grepl("^[A-Za-z_][A-Za-z0-9_]*=", toks)]     # strip PERIOD=... prefixes
    if (length(toks) && grepl("\\.R$", toks[1])) toks[1] else character(0)
  }))
}
scripts_tbl <- unique(do.call(rbind, lapply(seq_len(nrow(man)), function(i) {
  p <- extract_script_paths(man$script[i])
  if (!length(p)) return(NULL)
  data.frame(id = man$id[i], script_path = p, stringsAsFactors = FALSE)
})))
scripts_tbl$found <- file.exists(here::here(scripts_tbl$script_path))
not_found <- scripts_tbl[!scripts_tbl$found, ]
cat(sprintf("\n=== producer scripts: %d distinct paths referenced, %d found, %d NOT FOUND ===\n",
    length(unique(scripts_tbl$script_path)),
    length(unique(scripts_tbl$script_path[scripts_tbl$found])),
    length(unique(not_found$script_path))))
if (nrow(not_found)) cat(sprintf("SCRIPT NOT FOUND (fix the label or restore the file):\n  %s\n",
    paste(sprintf("%s -> %s", not_found$id, not_found$script_path), collapse = "\n  ")))
# entries whose label yields no checkable .R token (should be none; flag rather than hide)
unchecked <- man$id[!vapply(man$script, function(l) length(extract_script_paths(l)) > 0, logical(1))]
if (length(unchecked)) cat(sprintf("%d entries with no checkable script path: %s\n",
    length(unchecked), paste(unique(unchecked), collapse = ", ")))
