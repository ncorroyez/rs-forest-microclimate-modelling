# ==============================================================================
# run_chapter1_chs41.R — orchestrate the CHS41 figure chain of the live Chapter 1
# manuscript (Fig 4, 5, 6, 7, B1, S1). Companion to run_chapter1.R, which does NOT
# cover this chain. See RUN_CHAPTER1_CHS41.md for the full trace.
#
#   Rscript run_chapter1_chs41.R            # --check (default): verify assets + print plan, no side effects
#   Rscript run_chapter1_chs41.R --run      # extract -> figures -> sync (overwrites the six figures)
#   Rscript run_chapter1_chs41.R --run --no-sync
#   CH1_RUN_MUSICA=TRUE Rscript run_chapter1_chs41.R --run   # also rebuild missing NetCDF assets (heavy)
#
# RE-EXTRACTION + FIGURES only by default; NetCDFs/tables already exist on disk. A missing
# MuSICA asset with CH1_RUN_MUSICA unset is a hard error naming the script that rebuilds it.
# ==============================================================================
suppressMessages({ library(here) })
setwd(here::here())
ARGS      <- commandArgs(trailingOnly = TRUE)
DO_RUN    <- "--run"     %in% ARGS
NO_SYNC   <- "--no-sync" %in% ARGS
RUN_MUSICA<- identical(toupper(Sys.getenv("CH1_RUN_MUSICA")), "TRUE")
TBL <- "out_files/Chapter1/tables"
NCB <- "out_files/Chapter1"
FIGW<- "out_files/Chapter1/figures"                 # where the figure scripts write
FIGM<- "chapter1/manuscript/figures"               # manuscript dir (parent), where the .md reads

say <- function(...) cat(sprintf(...), "\n")
Rscript <- function(cmd) {
  say("  $ Rscript %s", cmd)
  st <- system2("Rscript", strsplit(cmd, " ")[[1]], stdout = TRUE, stderr = TRUE)
  code <- attr(st, "status"); if (!is.null(code) && code != 0) {
    cat(paste(st, collapse = "\n"), "\n"); stop("step failed: ", cmd, call. = FALSE) }
  invisible(st)
}

# MuSICA-produced assets: path -> the script that rebuilds it (only run with CH1_RUN_MUSICA=TRUE)
ASSETS <- list(
  list(path = file.path(TBL, "perturb_chs41_nowind.csv"),  musica = "chapter1/producers/c1_perturb_chs41_nowind.R 1 1"),
  list(path = file.path(NCB, "nc_perturb_chs41_nowind"),   musica = "chapter1/producers/c1_perturb_chs41_nowind.R 1 1"),
  list(path = file.path(NCB, "nc_archetype_chs41"),        musica = "chapter1/producers/c1_archetype_profiles_chs41.R"),
  list(path = file.path(NCB, "nc_b3_chs41"),               musica = "chapter1/producers/c1_b3_pervariable_chs41.R"),
  list(path = file.path(TBL, "h2_controlled_chs41.csv"),   musica = "chapter1/producers/c1_h2_controlled_chs41_clean.R"),
  list(path = file.path(TBL, "h2_gaussian_grid_chs41.csv"),musica = "chapter1/producers/c1_h2_gaussian_grid_chs41.R")
)
EXTRACT <- c("chapter1/producers/c1_perturb_chs41_slope_extract.R", "chapter1/producers/c1_hotdays_chs41.R")
FIGURES <- c("chapter1/producers/c1_fig4_fig5_chs41.R", "chapter1/producers/c1_fig6v2_gaussian_heatmap.R", "chapter1/producers/c1_fig6_fig7_chs41.R",
             "chapter1/producers/c1_b1_pervariable_plot_chs41.R", "chapter1/producers/c1_s1_attribution_hot_plot_chs41.R")
# manuscript figure  <-  file in FIGW to copy (Fig 6 = the gaussian panel, NOT h2_controlled)
SYNC <- c("Fig4_attribution_chs41.png", "Fig5_operating_point_chs41.png", "Fig6_h2_gaussian_chs41.png",
          "Fig7_vertical_Tprofile_chs41.png", "FigB1_pervariable_gradient_chs41.png",
          "FigS1_attribution_hot_chs41.png")

# ---- asset check ----
say("=== CHS41 chain: asset check ===")
missing <- Filter(function(a) !file.exists(a$path), ASSETS)
for (a in ASSETS) say("  [%s] %s", if (file.exists(a$path)) "X" else "MISSING", a$path)
if (length(missing)) {
  if (!RUN_MUSICA) {
    say("\n%d MuSICA asset(s) missing. Rebuild (heavy) with CH1_RUN_MUSICA=TRUE, or run:", length(missing))
    for (a in missing) say("  Rscript %s", a$musica)
    if (DO_RUN) stop("cannot --run with assets missing and CH1_RUN_MUSICA unset", call. = FALSE)
  } else if (DO_RUN) {
    say("\nCH1_RUN_MUSICA=TRUE: rebuilding %d missing asset(s) [HEAVY]", length(missing))
    for (a in missing) Rscript(a$musica)
  }
}

say("\n=== plan ===")
say("extract : %s", paste(EXTRACT, collapse = ", "))
say("figures : %s", paste(FIGURES, collapse = ", "))
say("sync    : %d figures  %s -> %s%s", length(SYNC), FIGW, FIGM,
    if (NO_SYNC) "  (SKIPPED: --no-sync)" else "")

if (!DO_RUN) { say("\n--check only (no side effects). Re-run with --run to execute."); quit(save = "no") }

say("\n=== extract ===");  for (s in EXTRACT) Rscript(s)
say("\n=== figures ===");  for (s in FIGURES) Rscript(s)

if (!NO_SYNC) {
  say("\n=== sync -> manuscript dir ===")
  for (fn in SYNC) {
    src <- file.path(FIGW, fn); dst <- file.path(FIGM, fn)
    if (!file.exists(src)) { say("  MISSING SOURCE %s", src); next }
    file.copy(src, dst, overwrite = TRUE, copy.date = TRUE); say("  %s -> %s", fn, FIGM)
  }
}
say("\n=== ledger ==="); Rscript("chapter1/ledger/make_article_figure_set.R")
say("\nDONE (CHS41 chain).")
