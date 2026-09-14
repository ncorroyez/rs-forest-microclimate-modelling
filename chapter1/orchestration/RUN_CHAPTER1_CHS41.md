# RUN_CHAPTER1_CHS41 — the CHS41 figure chain (Fig 4, 5, 6, 7, B1, S1)

Traced 2026-09-14 (harmonisation). This is the chain `RUN_CHAPTER1.md` does NOT cover: the six
figures of the live manuscript (`manuscript_chap1_FINAL_coherence_2026-09-09.md`) produced on the
CHS41-Rmerge / no-wind forcing. Driver: `run_chapter1_chs41.R`. See also the figure ledger
`chapter1/ledger/make_article_figure_set.R` and the map `chapter1/manuscript/PRODUCER_FIGURE_MAP_2026-09-14.md`.

**Like RUN_CHAPTER1.md, this covers RE-EXTRACTION + FIGURES only. No step runs MuSICA by default.**
The NetCDFs and design tables already exist on disk; every figure is rebuilt by re-reading them.

Canonical lineage:

| item | value |
|---|---|
| forcing | `out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc` (CHS41 station + ERA5 + MERRA-2 h_sbl; NOT FR-Blo_2021_v2.nc — see forcing_hsbl_trap) |
| baseline wind | none (Appendix A reports the sensitivity only) |
| design sample | `out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds` (400 real 20 m pixels) |
| model | MuSICA v3.2.3, iterative ABL ("yoyo"), k = 0.5 |

## 0. Simulation assets — REQUIRED ON DISK (rebuilt only with CH1_RUN_MUSICA=TRUE)

Each is produced by a MuSICA-running script (thousands of sims); never rebuild casually.

| asset | content | rebuilt by (MuSICA) |
|---|---|---|
| `out_files/Chapter1/nc_perturb_chs41_nowind/` + `tables/perturb_chs41_nowind.csv` | trait-perturbation design (per plot × arm) | `c1_perturb_chs41_nowind.R <chunk> <K>` |
| `out_files/Chapter1/nc_archetype_chs41/` (8) | per-archetype real/uniform vertical profiles | `c1_archetype_profiles_chs41.R` |
| `out_files/Chapter1/nc_b3_chs41/` (16) | per-variable within-canopy gradient sims | `c1_b3_pervariable_chs41.R` |
| `tables/h2_controlled_chs41.csv` | controlled H2 (ΔTmax vs centre of mass) | `c1_h2_controlled_chs41_clean.R` |
| `tables/h2_gaussian_grid_chs41.csv` | Gaussian LAD sweep, 5 μ × 5 σ × 5 LAI = 125 sims | `c1_h2_gaussian_grid_chs41.R` |

Auxiliary CHS41 tests (not among the six figures): `c1_h1_realvsreal_chs41.R`,
`c1_noise_floor_chs41.R`, `c1_bias_origin_chs41.R`, `c1_perturb_chs41_slope_extract.R` (extract only).

## 1. Extract (no MuSICA)

| step | command | reads | writes |
|---|---|---|---|
| E1 | `Rscript c1_perturb_chs41_slope_extract.R` | nc_perturb_chs41_nowind/ | `tables/perturb_chs41_nowind_SLOPE.csv` |
| E2 | `Rscript c1_hotdays_chs41.R` | nc_perturb_chs41_nowind/ | `tables/hotdays_chs41.csv` |

## 2. Figures (no MuSICA) — write to `out_files/Chapter1/figures/`

| fig | command | needs |
|---|---|---|
| Fig 4 + Fig 5 | `Rscript c1_fig4_fig5_chs41.R` | perturb_chs41_nowind.csv + _SLOPE.csv |
| Fig 6 | `Rscript c1_fig6v2_gaussian_heatmap.R` | h2_gaussian_grid_chs41.csv |
| Fig 7 | `Rscript c1_fig6_fig7_chs41.R` | h2_controlled_chs41.csv + nc_archetype_chs41/ |
| Fig B1 | `Rscript c1_b1_pervariable_plot_chs41.R` | nc_b3_chs41/ |
| Fig S1 | `Rscript c1_s1_attribution_hot_plot_chs41.R` | hotdays_chs41.csv |

**HAZARD.** `c1_fig6_fig7_chs41.R` also writes `Fig6_h2_controlled_chs41.png` into the same dir. The
manuscript's Fig 6 is the GAUSSIAN panel from `c1_fig6v2_gaussian_heatmap.R`. The sync (§3) copies the
gaussian file explicitly; never hand-copy `Fig6_h2_controlled_chs41.png`.

Only `c1_b1_*` and `c1_s1_*` honour `CH1_CHS41_FIGDIR` (rehearse into a scratch dir). The Fig 4/5/6/7
scripts hardcode `out_files/Chapter1/figures/`, so a full run overwrites those working-dir figures.

## 3. Sync into the manuscript directory

Copies the six CHS41 figures to `chapter1/manuscript/figures/` (parent, NOT article_v323/). Fig 8 is
NOT part of this chain (it comes from `run_chapter1.R`, radius sweep).

  Fig4_attribution_chs41.png  Fig5_operating_point_chs41.png  Fig6_h2_gaussian_chs41.png
  Fig7_vertical_Tprofile_chs41.png  FigB1_pervariable_gradient_chs41.png  FigS1_attribution_hot_chs41.png

## 4. Driver

```bash
Rscript run_chapter1_chs41.R            # --check (default): verify assets, print the plan, no side effects
Rscript run_chapter1_chs41.R --run      # extract -> figures -> sync (overwrites the six figures)
Rscript run_chapter1_chs41.R --run --no-sync   # regenerate but leave the manuscript dir untouched
CH1_RUN_MUSICA=TRUE Rscript run_chapter1_chs41.R --run   # additionally rebuild missing NetCDF assets (heavy)
```

A missing asset with `CH1_RUN_MUSICA` unset is a hard error naming the MuSICA script that rebuilds it.
