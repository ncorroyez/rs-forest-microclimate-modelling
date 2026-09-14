# Article Chap1 rerun — path/branch reconciliation (do this WHEN the station forcing arrives)

Companion to `run_article_v323.R`. The rerun is NOT a one-button job yet, because the figure
scripts read **canonical / v3.2.0** locations while the v3.2.3-iter sims write **`*_iter/` or the
`z05_v323` branch**. This file is the exact, ordered checklist to route everything through v3.2.3.
Everything here must be **run and verified with the real forcing** — do not pre-apply blind.

## The two branches (key fact)
- `outputs/figures_pipeline_z05/`  = **v3.2.0** (current canonical; what the manuscript numbers come from).
- `outputs/figures_pipeline_z05_v323/` = **v3.2.3** branch (ALREADY EXISTS, has its own `data/ref_validation.rds`).
- Per-plot sensitivity dirs: `sensitivity_perplot/` (v3.2.0) vs `sensitivity_perplot_iter/` and `sensitivity_perplot_v323/`.
- metrics6: `metrics6_v320/` vs `metrics6_iter/`  (the only chain already VER-parameterised, via `c1_deltadelta_perplot.R`).

## Recommended approach = PROMOTE iter → canonical (one file op, zero figure-script edits)
After the iter sims run, archive v3.2.0 and promote iter into the canonical names; then every downstream
generator + figure runs unchanged on v3.2.3. Concretely, in the DATA phase:

1. `build_forcing_pblh.R`  (new station forcing → `musica_in_Blois_pblh.nc`).
2. `c1_sensitivity_perplot_chunk_iter.R 1 1`  → `sensitivity_perplot_iter/`.
3. **PROMOTE:** `mv sensitivity_perplot sensitivity_perplot_v320archive && cp -r sensitivity_perplot_iter sensitivity_perplot`
   (archive-first so v3.2.0 is never lost). Do the same for the validation branch:
   point the figure set at `figures_pipeline_z05_v323/` (set `PIPE_BRANCH=z05_v323`) OR promote its
   `data/ref_validation.rds` into `figures_pipeline_z05/data/` (archive first).
4. Re-run the SIM-derived table generators so they rebuild from the now-canonical iter per-plot:
   - `c1_sensitivity_perplot_lad.R`   → `tab_sensitivity_perplot_lad.csv`  (runs MuSICA: real-vs-uniform LAD contrast — needs the iter binary+forcing).
   - `c1_sensitivity_percluster.R`    → `tab_sensitivity_percluster.csv`.
   - `c1_importance_sensitivity.R`    → `tab_importance_sensitivity.csv`  (only the inclusion ORDER for Fig 5; the bar figure itself is dropped).
   - `c1_metrics6_perplot.R 1 1 iter` → `metrics6_iter/`.
5. Station-ref + validation + VCI: `c1_deltatmax_station_ref`, `c1_slope_station_ref`,
   `c1_hobo_iter_validation`, `c1_vci_points_tilewise` + `c1_vci_finalize`.

## Figure → derived-table → iter-readiness map
| Fig | script | reads | iter-ready after step |
|-----|--------|-------|-----------------------|
| Fig2/3 | c1_deltadelta_perplot (VER=iter, NORM=sd) | `metrics6_iter/` | YES already (VER-param) |
| Fig4 validation | (z05_v323 branch) | `figures_pipeline_z05_v323/` | YES (branch exists) — set PIPE_BRANCH |
| Fig5 forward | c1_forward_percluster | `coal_metrics.rds`, `ref_validation.rds`, `tab_importance_sensitivity.csv` | after PROMOTE + step 4 |
| Fig6 vertical | make_vertical_profiles_wind_rh_vpd | sim profiles | needs an iter profile run |
| A1 potency/local | c1_sensitivity_figure | `tab_sensitivity_percluster.csv` | after step 4 |
| A2 lad/VCI | c1_lad_vci | `tab_sensitivity_perplot_lad.csv` + corrected VCI | after step 4 |
| B1 residual | make_residual_vs_topo | `ref_validation.rds` | after PROMOTE / branch |
| E1/E2 corrplot, E3 PCA, D1 dispersion, #14 | (read cLHS sample) | sample + corrected VCI | YES (version-independent) |
| F1/F2 obs-attribution, F3 coverage | c1_obs_attribution, c1_hobo_coverage | OBSERVED buffering + traits | YES (model-free, version-independent) |

## Caveats to verify on rerun
- `coal_metrics.rds` is a Shapley-era artefact (`pipeline/04_extract.R`); confirm `c1_forward_percluster` still
  needs it or refactor it out (forward order can come from the sensitivity importance directly).
- After PROMOTE, re-run `chapter1/ledger/make_article_figure_set.R` and check the manifest = 0 MISSING.
- Re-validate r vs the 53 HOBO under v3.2.3 + station forcing: prior v3.2.3 gave r≈0.63 under ERA5; the
  station forcing is the open question. If r stays low, the "r≈0.93 ranking" headline must be reworded.
- `c1_sensitivity_perplot_lad.R` / `c1_sensitivity_percluster.R` RUN MuSICA — they need the iter binary
  (`in_files/model-3.2.3/musica`) and the new forcing, not just a path swap.
