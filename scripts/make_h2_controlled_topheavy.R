# ==============================================================================
# R2bis — Test CONTRÔLÉ de H2 : à LAI et Hmax FIXÉS, la position verticale de la
# biomasse foliaire (bottom-heavy → top-heavy) change-t-elle le ΔTmax ?
#
# φ_LAD ≈ 0 dans les données réelles corrélées (Shapley) ; ici on isole l'effet
# PUR de la forme verticale : mêmes LAI & Hmax, seul varie le centre de masse du
# profil LAD (famille Beta). Si ΔTmax décroît quand la biomasse monte → H2 tient
# en conditions contrôlées. Si plat → confirme que la forme verticale n'agit pas.
#
# Reads :
#         out_files/H2_controlled_topheavy_chapter/LAI1s<lai>/musica_out_<profile>.nc
#           the 42 cached runs. MuSICA runs ONLY for a NetCDF that is missing or < 1 kB.
# Writes: outputs/figures_pipeline/annex/fig_h2_controlled_topheavy.{png,pdf}
#          + tab_h2_controlled_topheavy.csv
# ==============================================================================

suppressMessages({
  library(here); library(ncdf4); library(tidyverse); library(lubridate)
  library(rmusica); library(musica.tools)
})
setwd(here::here())                       # musica.nml / variables.csv sont à la racine
source(here::here("scripts/_article_style.R"))
source(here::here("R/config.R"))
source(here::here("R/io.R"))
source(here::here("R/lad.R"))
source(here::here("R/musica.R"))
source(here::here("R/dtmax_convention.R"))   # canonical time-matched ΔTmax (convention B)

OUT <- here::here("outputs/figures_pipeline/annex")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SIMDIR <- here::here("out_files/H2_controlled_topheavy")

# ---- LINEAGE. The 14 cached runs were produced with CFG's defaults: the LEGACY binary
# (~/Documents/musica/musica), the ERA5 forcing (musica_in_Blois.nc) and NO ABL coupling.
# That is NOT the chapter's configuration, and Fig. B1's caption used to claim otherwise.
# Set H2_LINEAGE=chapter to re-run the same 14 configurations on the chapter lineage
# (v3.2.3 binary, station-based FR-Blo_2021_v2.nc, abl_flag='iter') in a SEPARATE directory,
# leaving the cached legacy runs untouched, and see whether the conclusion survives.
# DEFAULT FLIPPED TO "chapter" 2026-07-31: every result in Chapter 1 must be on v3.2.3 +
# iterative ABL + CHS 41 station forcing. "legacy" is kept only to re-read the superseded
# cached runs for comparison; it must not produce a shipped figure.
LINEAGE <- tolower(Sys.getenv("H2_LINEAGE", "chapter"))
if (!LINEAGE %in% c("legacy", "chapter")) stop("H2_LINEAGE must be 'legacy' or 'chapter'")


# ---- design : LAI & Hmax fixés, profils Beta de centre de masse croissant -----
HMAX     <- 25                            # m, hauteur fixe
FCOVER   <- 0.87                          # clumping_factor fixe (≈ moyenne massif)
# LAI CONVENTION, FIXED 2026-07-31. These are ONE-SIDED leaf areas, the units the LiDAR
# retrieval and the whole chapter report. run_musica_one() doubles whatever lai_fn returns
# (R/musica.R: `lai <- 2 * lai`), so MuSICA receives the two-sided value it expects and no
# doubling is applied here. The previous LAI_LEVELS <- c(6, 12) were therefore fed as
# one-sided 6 and 12, not the "one-sided 3 and 6" the Fig. B1 caption claimed: a factor of
# two, and one-sided 12 is far outside the design, whose LAI runs 0.50 to 7.11 (median 3.48).
# The levels below span the real range instead.
LAI_LEVELS <- c(1, 2, 3, 4, 5, 6)         # ONE-SIDED LAI, capped at 6 (two-sided 12), the realistic ceiling for this oak forest (design max 7.11)

if (LINEAGE == "chapter") {
  SIMDIR <- here::here("out_files/H2_controlled_topheavy_chapter")
  source(here::here("R/wind_correction.R"))
  # The chapter baseline scales the forcing wind to just above the canopy. HMAX is fixed
  # here, so this is one constant factor, but it must be applied or the test would sit on a
  # different aerodynamic baseline from every other result.
  CFG$forcing_file <- windcorr_forcing(HMAX)
  CFG$musica_cmd   <- normalizePath("in_files/model-3.2.3/musica", mustWork = TRUE)
  EXTRA <- list("abl_flag" = '"iter"')
} else EXTRA <- list()
dir.create(SIMDIR, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("H2 lineage: %s | forcing %s | sims %s\n", LINEAGE,
            basename(CFG$forcing_file), basename(SIMDIR)))
# (a,b) Beta : a<b = bottom-heavy, a=b = centré, a>b = top-heavy
PROFILES <- list(
  list(lab = "very_bottom", a = 2, b = 6),
  list(lab = "bottom",      a = 2, b = 4),
  list(lab = "mid_low",     a = 2, b = 3),
  list(lab = "centered",    a = 2, b = 2),
  list(lab = "mid_high",    a = 3, b = 2),
  list(lab = "top",         a = 4, b = 2),
  list(lab = "very_top",    a = 6, b = 2)
)

# ---- constructeur de profil LAD synthétique (intègre à lai, grille 1 m) -------
#' Build a synthetic Beta LAD profile generator on a 1 m grid
#' @param a first Beta shape parameter (a < b bottom-heavy, a > b top-heavy)
#' @param b second Beta shape parameter
#' @param canopy_base lowest height (m) carrying foliage
#' @return a FUNCTION(plot_row, hmax, lai) returning data.frame(height, density) integrating to lai
make_lad_beta <- function(a, b, canopy_base = 1) {
  force(a); force(b); force(canopy_base)
  function(plot_row, hmax = NULL, lai = NULL) {
    if (is.null(hmax)) hmax <- as.numeric(plot_row$Hmax)
    if (is.null(lai))  lai  <- as.numeric(plot_row$LAI)
    max_layer <- ceiling(hmax)
    heights   <- seq_len(max_layer)
    u <- (heights - canopy_base + 0.5) / (max_layer - canopy_base + 1)   # ∈(0,1)
    w <- ifelse(heights >= canopy_base, dbeta(pmin(pmax(u, 1e-3), 1 - 1e-3), a, b), 0)
    w[!is.finite(w)] <- 0
    dens <- if (sum(w) > 0) lai * w / sum(w) else rep(0, max_layer)
    data.frame(height = heights, density = dens)
  }
}

# centre de masse (fraction de Hmax) d'un profil — indépendant du LAI
#' Centre of mass of a LAD profile, as a fraction of canopy height
#' @param lad_fn a profile generator as returned by make_lad_beta()
#' @param hmax canopy height (m), defaults to the fixed HMAX of the design
#' @return numeric in (0, 1); independent of LAI
com_frac <- function(lad_fn, hmax = HMAX) {
  p <- lad_fn(NULL, hmax = hmax, lai = 1)
  sum(p$density * p$height) / (sum(p$density) * hmax)
}

dummy <- data.frame(Hmax = HMAX, LAI = NA_real_, fCover = FCOVER)
# CONVENTION FIX 2026-07-31. This script used extract_deltatmax_one() from R/musica.R,
# the LEGACY extractor: it takes the sub-canopy daily maximum independently of the macro
# one and applies a -2 h clock shift. Every other stage of the chapter uses the canonical
# time-matched convention B of R/dtmax_convention.R (sub-canopy read at the hour of the
# macro daily maximum, no clock shift), so Fig. B1 was the one result on a different
# metric. Pure re-extraction from the cached NetCDFs; no simulation is re-run.
MREF <- macro_ref(CFG$forcing_file, CFG$date_seq)

# ---- run + extraction --------------------------------------------------------
rows <- list()
for (lai in LAI_LEVELS) {
  outdir <- file.path(SIMDIR, sprintf("LAI1s%g", lai))   # 1s = one-sided, vs the old ambiguous LAI6/LAI12
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (pf in PROFILES) {
    lad_fn <- make_lad_beta(pf$a, pf$b)
    sc <- list(name = pf$lab,
               lai_fn    = function(pr) lai,
               hmax_fn   = function(pr) HMAX,
               fcover_fn = function(pr) FCOVER,
               lad_fn    = lad_fn)
    out_nc <- normalizePath(file.path(outdir, sprintf("musica_out_%s.nc", pf$lab)),
                            mustWork = FALSE)
    # CACHE GUARD 2026-07-30: the 14 controlled-test NetCDFs are on disk under SIMDIR,
    # so this figure is re-makeable by pure re-extraction. Only a missing or truncated
    # output triggers MuSICA, which keeps the script honest if the sims are ever cleared.
    if (!file.exists(out_nc) || file.size(out_nc) < 1000)
      run_musica_one(dummy, sc, out_nc, CFG$forcing_file, CFG$musica_cmd, extra_setup = EXTRA)
    mh <- micro_hourly_at(out_nc, 1.0)
    if (is.null(mh) || nrow(mh) == 0) { cat(sprintf("  [skip] LAI%g %s\n", lai, pf$lab)); next }
    dd <- delta_tmax(mh, MREF)[date %in% CFG$date_seq]
    if (nrow(dd) == 0) { cat(sprintf("  [skip] LAI%g %s\n", lai, pf$lab)); next }
    rows[[length(rows) + 1L]] <- data.frame(
      LAI = lai, profile = pf$lab,
      com = com_frac(lad_fn),
      dTmax = mean(dd$Delta, na.rm = TRUE)
    )
    cat(sprintf("  LAI(1s)%-3g %-11s COM=%.2f  ΔTmax=%+.2f\n",
                lai, pf$lab, tail(rows, 1)[[1]]$com, tail(rows, 1)[[1]]$dTmax))
  }
}
D <- bind_rows(rows)
D$LAI <- factor(D$LAI, levels = LAI_LEVELS, labels = sprintf("%g", LAI_LEVELS))
write_csv(D, file.path(OUT, "tab_h2_controlled_topheavy.csv"))

# pente ΔTmax par unité de centre de masse (quantifie l'effet vertical pur)
slopes <- D %>% group_by(LAI) %>%
  summarise(slope = coef(lm(dTmax ~ com))[2],
            span  = max(dTmax) - min(dTmax), .groups = "drop")
print(as.data.frame(D)); print(as.data.frame(slopes))

# ---- figure ------------------------------------------------------------------
p <- ggplot(D, aes(com, dTmax, colour = LAI, group = LAI)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey45", linewidth = 0.5) +   # no-effect reference (committee)
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_colour_viridis_d(option = "D", end = 0.92, direction = -1,
                         name = "one-sided LAI") +
  labs(x = "LAD profile center of mass  (fraction of Hmax: 0 = bottom, 1 = top)",
       y = expression(bar(Delta*T[max])~"summer  ("*degree*"C)"),
       subtitle = sprintf("LAI, cover and Hmax (%g m) fixed; only the vertical shape varies  |  bottom-heavy to top-heavy", HMAX)) +
  theme_article() + legend_corner(0.01, 0.01)
ggsave_article(file.path(OUT, "fig_h2_controlled_topheavy"), p, 7.5, 5.5)

# Where does each LAI level put its optimum, and how big is the lever there?
opt <- D %>% group_by(LAI) %>%
  summarise(com_opt = com[which.min(dTmax)], dTmax_opt = min(dTmax),
            span = max(dTmax) - min(dTmax), .groups = "drop")
cat("\n=== optimum balance and lever size by one-sided LAI ===\n"); print(as.data.frame(opt))
write_csv(opt, file.path(OUT, "tab_h2_optimum_by_lai.csv"))
cli::cli_alert_success("Saved fig_h2_controlled_topheavy (png + pdf) + table")
