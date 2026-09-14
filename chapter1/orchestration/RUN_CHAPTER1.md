# RUN_CHAPTER1 — regenerating the Chapter 1 article analysis and figures

> **⚠️ SUPERSEDED IN PART (2026-09-14). This document describes the PRE-CHS41 lineage.**
> The stages below reproduce `manuscript_chap1_EN_native20.md` (+ `_SHORT.md`), NOT the live
> manuscript `manuscript_chap1_FINAL_coherence_2026-09-09.md`. The live manuscript was
> restructured: **Fig 4, 5, 6, 7, B1 and S1 come from a hand-run CHS41 chain this document
> does not cover** (`c1_fig4_fig5_chs41.R`, `c1_fig6v2_gaussian_heatmap.R`,
> `c1_fig6_fig7_chs41.R`, `c1_b3_pervariable_chs41.R`, `c1_hotdays_chs41.R`, and the upstream
> `c1_*_chs41.R` sims), and the logger-validation figures the stages below still build were
> DROPPED (validation moved to Chapter 3). Also: `run_chapter1.R --sync-only` reaches only
> `figures/article_v323/`, not `figures/` where the CHS41 figures live, so their copy is manual.
> The B1 and S1 figure-plotting producers are **not located** (files exist on disk, provenance
> unknown). The accurate figure→producer map for the live manuscript is
> `chapter1/manuscript/PRODUCER_FIGURE_MAP_2026-09-14.md`; the corrected ledger is
> `chapter1/ledger/make_article_figure_set.R`. Rewriting these stages onto the CHS41 narrative is a
> pending code chantier. Until then, treat the run order below as historical.

Derived 2026-07-29 by tracing `chapter1/ledger/make_article_figure_set.R` (the figure ledger),
the two manuscripts (`chapter1/manuscript/manuscript_chap1_EN_native20.md` and
`..._SHORT.md`), and the actual `ggsave` / `fwrite` targets of every script involved.

**This document covers RE-EXTRACTION AND FIGURES ONLY.** No step runs MuSICA by default. The single
exception is stage **B24** (`scripts/c1_blh_insensitivity.R`), which backs the Appendix H
boundary-layer-height claim with 16 short runs (about 4 minutes) and is skipped unless you set
`CH1_RUN_MUSICA=TRUE`. It writes only under `out_files/Chapter1/nc_blh_insensitivity` and its own
forcing cache, so it cannot disturb a canonical output.
Every number in the chapter is recovered by re-reading NetCDF outputs that already
exist on disk.

Canonical lineage of the chapter, for reference:

| item | canonical value |
|---|---|
| forcing | `in_files/FR-Blo_2021_v2.nc` (CHS 41 station based, solar time = UTC+1) |
| ΔTmax convention | `R/dtmax_convention.R`, time-matched, `OBS_CLOCK_OFFSET_H = -1` on **observations only** |
| design sample | `out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds` (400 real 20 m pixels) |
| archetype labels | `out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv`, column `P`, frozen |
| model | MuSICA v3.2.3, iterative ABL, canopy wind log-profile correction ON |

Superseded lineages that must **not** appear in a Chapter 1 number: `frblo`,
`z05`, `era5sd`, `v320`, `floor05_v2`, `convA`, `transfer/`,
`outputs/figures_pipeline_z05/`.

---

## 0. Simulation assets — REQUIRED ON DISK, NEVER DELETE

None of these can be rebuilt without thousands of CPU-hours of MuSICA. Every step
below reads them.

| path | content | size |
|---|---|---|
| `out_files/Chapter1/nc_sensitivity_perplot_units/` | trait-perturbation design, 8 tags × 400 plots = **3200 NetCDFs** | **8.6 GB** |
| `out_files/musica_native20_forward/` | 16 coalition sets (`0000`…`1111`) × 53 loggers | 2.6 GB |
| `out_files/musica_hobo_native20/1111` | full-model validation run, 53 loggers | — |
| `out_files/radius_test_corr/{5m,10m,12.5m,15m,20m,25m,50m}` | footprint-radius sweep (wind-corrected set) | — |
| `out_files/H2_controlled_topheavy/` | Appendix C1 controlled test, 14 NetCDFs | 58 MB |
| `in_files/FR-Blo_2021_v2.nc` | canonical station forcing | — |

Verify before starting:

```bash
cd /home/corroyez/Documents/z_Example_rmusica_31012025
ls out_files/Chapter1/nc_sensitivity_perplot_units | wc -l     # must be 3200
ls out_files/musica_native20_forward | wc -l                   # must be 16
ls out_files/radius_test_corr | wc -l                          # must be 7
ls -1 in_files/FR-Blo_2021_v2.nc out_files/musica_hobo_native20/1111 >/dev/null
```

## Scripts that MUST NOT be run

| script | why |
|---|---|

| `run_article_v323.R` | a stale one-command driver at repo root. Its DATA phase runs `c1_sensitivity_perplot_chunk_iter.R` and `c1_sensitivity_perplot_lad.R`, both of which invoke MuSICA; it also targets the pre-station forcing `in_files/musica_in_Blois.nc` and the superseded `floor05_v2` sample. Superseded by this document. |
| `pipeline/run_all.R` and stages `03*`–`11` | the Shapley / `floor05_v2` / `z05` pipeline. Superseded (see §5). Stage `03_musica.R` runs MuSICA. |
| `c1_hobo_direction.R` | historically overwrote the canonical observation table from a different branch. Now guarded (writes `tab_hobo_direction_z05branch.csv`), but it is not part of the article chain. |

## The one script that rewrites the canonical observation table

`c1_align_cluster_obs.R` is the **sole** writer of
`out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv`. It rewrites only
`dTmax_obs` and `slope_obs` from `tab_hobo_native20_validation.csv` and asserts the
`P` labels are unchanged. It is **opt-in** in the runner (`CH1_ALIGN_OBS=TRUE`) and
should only be run after step A1 has regenerated the validation table. Everything
else in the chapter reads this table and never writes it.

---

## 1. Stage A — re-extraction from the existing NetCDFs

Run in this order. All timings measured on this machine, 2026-07-29.

### A1. Full-model validation against the 53 loggers — ~30 s

```bash
Rscript c1_hobo_native20_validation_regen.R
```
* consumes `out_files/musica_hobo_native20/1111`, `in_files/FR-Blo_2021_v2.nc`, HOBO csv (`CFG$hobo_temp_csv`)
* produces `out_files/Chapter1/tables/tab_hobo_native20_validation.csv`
* expected: `dTmax r = 0.501 | bias = +0.98 °C | RMSE 1.64 | amplitude 15%`, `slope r = 0.782`

### A2. (opt-in) Align the observed columns of the canonical table — ~5 s

```bash
CH1_ALIGN_OBS=TRUE Rscript c1_align_cluster_obs.R    # only if A1 changed
```
* rewrites **only** `dTmax_obs` / `slope_obs` in `tab_hobo_perplot_cluster.csv`; asserts `P` unchanged
* idempotent: re-running with an unchanged A1 output is a no-op

### A3. Trait-perturbation design, 2 metrics × 2 periods — ~30 s

```bash
Rscript c1_sensitivity_units_extract2.R
```
* consumes the **3200 NetCDFs** in `out_files/Chapter1/nc_sensitivity_perplot_units/`
* produces `out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv`
* self-checks: recomputed ΔTmax must reproduce the original extraction (expect max deviation ~5e-15)

### A4. Hot-day (top 10%) subset over the coalition sets — ~1 min

```bash
Rscript c1_hot_extract.R
```
* consumes `out_files/musica_native20_forward/` (16 bits × 53 loggers)
* produces `out_files/Chapter1/hot_extract.rds`
* expected: `hot days: 13`, `obs 53 loggers, coal 848 rows (0 NA)`

### A5. Footprint-radius sweep, aligned clock — ~40 s

```bash
Rscript c1_radius_rescore.R
```
* consumes `out_files/radius_test_corr/*`
* produces `out_files/Chapter1/tables/tab_radius_sweep_aligned.csv`
* expected: `r` rising 0.414 (5 m) → 0.531 (50 m); P1 negative at every radius

### A6. Correction-level rescoring (wind / scan angle / k) — ~1 min

```bash
Rscript c1_correction_levels_rescore.R
```
* produces `tab_correction_levels.csv`, `tab_correction_levels_perplot.csv`, `tab_correction_levels_split.csv`
* expected biases `+0.92 / +1.03 / +0.98 °C`, `r = 0.53 / 0.49 / 0.50`, slope `r = 0.70 / 0.72 / 0.78`

### A7. Bootstrap CIs on the trait sensitivities (Appendix H) — ~2 min

```bash
Rscript c1_units_bootstrap_v2.R
```
* consumes `sensitivity_perplot_units_v2.csv`
* produces `out_files/Chapter1/tables/tab_attribution_units_bootstrap_v2.csv`

---

### A8. Table 1 — per-archetype traits, step/SD, profile swap — ~2 s

```bash
Rscript scripts/c1_make_table1.R
```

Emits `tab_table1_full.md` and `tab_table1_short.md` (the two manuscripts differ: the full one carries the
step/SD column). Paste-ready markdown, so Table 1 is no longer hand-maintained. **Warning:**
`out_files/Chapter1/tables/tab_cluster_meansd.csv` looks like the source but is the stale old-lineage version.

## 2. Stage B — figures

Each script is independent given Stage A. Order within the stage does not matter.

| # | command | produces (source path, no extension) | manuscript figure | time |
|---|---|---|---|---|
| B1 | `Rscript c1_fig1_gril_native.R` | `outputs/figures_chap1/Fig1_typology_gril_native` | Fig 1 | 8 s (self-contained: reads `in_files/fig1_crosssection_cache.rds`; the external drive is needed ONLY to rebuild that cache, i.e. if you delete it) |
| B2 | `Rscript scripts/c1_fig3_units.R` | `out_files/Chapter1/figures/Fig3_attribution_units` | Fig 3 (short) / Fig 4 (full) | 4 s |
| B3 | `Rscript scripts/c1_operating_point_main.R` | `outputs/figures_chap1/fig_operating_point_main` | Fig 4 (short) / Fig 5 (full) | 5 s |
| B4 | `Rscript scripts/c1_fig5_j1_convB.R` | `out_files/Chapter1/figures/Fig5_glass_ceiling_dumbbell` **and** `.../FigJ1_forward_inclusion` | Fig 5/6 + Fig J1/I1 | 13 s |
| B5 | `Rscript scripts/c1_fig6_obs_nested_native.R` | `out_files/Chapter1/figures/Fig6_obs_nested_native` | Fig 6 (short) / Fig 8 (full) | 6 s |
| B6 | `Rscript c1_pca_recover_clusters.R` | `out_files/Chapter1/figures/FigSh_pca_clusters` | Fig F2 (short) / Fig 3 (full) | 2 s |
| B7 | `Rscript scripts/c1_fig7_footprint_radius.R` | `out_files/Chapter1/figures/Fig7_footprint_radius` | Fig A1 (short) / Fig 7 (full) | 2 s |
| B8 | `Rscript c1_era5_station_bias.R` | `outputs/figures_chap1/B1_era5_station_bias` | Fig B1 / A1 | 3 s |
| B9 | `Rscript scripts/make_vertical_profiles_native.R` | `outputs/figures_chap1/fig_vertical_profiles_native_{norm,metres}` | Fig C2 / B2 (**after montage**, see §3) | 7 s |
| B10 | `Rscript scripts/make_trait_vertical_gradient_perplot.R` | `outputs/figures_chap1/fig_trait_vertical_gradient_perplot` | Fig C3 / B3 | 19 s |
| B11 | `Rscript scripts/make_residual_vs_topo.R` | `outputs/figures_pipeline/annex/fig_residual_vs_topo` | Fig D1 / C1 | 7 s |
| B12 | `Rscript scripts/c1_hobo_coverage_native.R` | `out_files/Chapter1/figures/FigSh_hobo_coverage_native` | Fig E1 / D1 | 6 s |
| B13 | `Rscript scripts/c1_obs_profileshape_test.R` | `outputs/figures_pipeline/annex/fig_obs_profileshape` | Fig F3 / E2 | 2 s |
| B14 | `Rscript c1_wind_profile_correction.R` | `outputs/figures_chap1/G1_wind_profile_correction` | Fig G1 / A2 | 2 s |
| B15 | `Rscript scripts/fig_scanangle_control_csv.R` | `outputs/figures_pipeline/annex/fig_scanangle_control` | Fig G2 / F1 | 2 s |
| B16 | `Rscript scripts/fig_lai_correction_compare.R` | `out_files/Chapter1/figures/F2_lai_correction_compare` | Fig G3 / F2 | 3 s |
| B17 | `Rscript c1_hourly_scatter.R` | `outputs/figures_chap1/Fig_hourly_temp_scatter` | Fig S1 | ~1 min |
| B18 | `PERIOD=hot METRIC=dt Rscript scripts/c1_fig3_units.R` | `out_files/Chapter1/figures/Fig3_attribution_units_hot` | Fig S2 | 4 s |
| B19 | `Rscript c1_hot_figures.R` | `outputs/figures_chap1/FigS_{glass_ceiling_hot,forward_hot,obs_nested_hot}` | Figs S3–S5 | 5 s |

Whole chain (stages A + B + ledger): **about 5 minutes** on this machine.

**Not regenerated here** (no script, or MuSICA-bound):

* **Fig 2 (methodo diagram)** — generated by stage B22 (`scripts/c1_fig2_methodo.R`) from the graphviz
  source `scripts/c1_fig2_methodo.dot`, which is the source of truth. It was hand-drawn until
  2026-08-03, and had drifted from the analysis (it showed 10 m traits and a +/-1 SD perturbation).
* **Fig C1 / B1 (controlled H2 test)** — cached at `outputs/figures_pipeline/annex/fig_h2_controlled_topheavy.png`; its generator runs MuSICA.
* **Fig F1 / E1 (trait collinearity)** — RESOLVED 2026-08-04. Stage B20
  (`scripts/make_trait_collinearity_native20.R`) writes it from the canonical native20 sample and
  reproduces the shipped png byte for byte. The former hazard note in §6 no longer applies.

---

## 3. Stage C — assemble the manuscript figure directory

The generating scripts write to three working directories
(`outputs/figures_chap1/`, `out_files/Chapter1/figures/`,
`outputs/figures_pipeline/annex/`). The manuscripts read
`chapter1/manuscript/figures/article_v323/` under **different names**, and under
**two different numbering schemes** (short doc / full doc). Historically this copy
was done by hand; the mapping below was recovered by md5-matching every cited
figure back to its source.

Rehearse it into a scratch directory first — this touches nothing in the manuscript:

```bash
CH1_FIGDIR=/tmp/figtest Rscript run_chapter1.R --sync-only
# expect: "43 files written, 0 sources missing" + "[montage] C2/B2 rebuilt"
```

Then, when satisfied:

```bash
Rscript run_chapter1.R --sync-only
```

`C2_vertical_gradient.png` / `B2_vertical_gradient.png` is **not** a plain copy: it
is a vertical montage of the two panels from B9 (normalized on top, metres below,
4800×3360 = 2 × 4800×1680). The runner rebuilds it with ImageMagick `convert`
(verified to reproduce the shipped dimensions exactly), and otherwise leaves the
existing file alone and warns.

## 4. Stage D — check the ledger

```bash
Rscript chapter1/ledger/make_article_figure_set.R
```

* writes `chapter1/manuscript/figures/article_v323/tab_article_figures.csv`
* must report **0 MISSING** across the 50 cited entries (25 per document)

---

## 5. Status of `pipeline/`

`pipeline/` is **superseded** for Chapter 1, with one exception.

Evidence:

* `pipeline/00_config.R:22` defaults the simulation directory to
  `out_files/musica_hobo_v10_fcovmean` (legacy v10 / MuSICA v3.2.0), while the
  chapter's figures come from `musica_native20_forward` and `musica_hobo_native20`.
* `pipeline/00_config.R:44` pins `CLUSTER_SAMPLE` to
  `clhs_sample_floor05_v2.rds`; the chapter's design sample is
  `clhs_sample_native20_floor05.rds`.
* `pipeline/00_config.R:56` asserts the published validation is `r_pp = 0.93`,
  split `44/9`; the chapter's headline validation is `r = 0.50`, bias `+0.98 °C`.
* The whole `05_shapley` / `05b` / `06` / `06b` / `10` / `10b` branch computes
  Shapley values. Shapley was dropped from Chapter 1 on 2026-06-20; neither
  manuscript cites a Shapley figure or table.

**Exceptions — do not delete `pipeline/`:**

1. `pipeline/00_config.R` is a **live dependency**. It is sourced by
   `c1_dumbbell_by_cluster.R`, `c1_radius_rescore.R`, `c1_hot_extract.R`,
   `c1_hourly_scatter.R` and others for `CFG$hobo_temp_csv`, `CFG$date_seq`,
   `CFG$ids_to_remove` and `PIPE$BITS`.
2. `pipeline/12_corrplot_traits.R` was long believed to be the true producer of the
   trait-collinearity figure shipped as Fig F1 / E1, on the superseded `floor05_v2` sample.
   RESOLVED 2026-08-04: stage B20 regenerates that exact path from the canonical native20
   sample and the result is byte-identical to the shipped figure, so the ledger is correct
   and the provenance hazard is closed.

Suggested (not applied): add a header comment to `pipeline/run_all.R` and to
stages `03`–`11` reading
`# SUPERSEDED for Chapter 1 (Shapley / floor05_v2 / z05 lineage). See RUN_CHAPTER1.md.`

---

## 6. Known hazards carried by this chain

Recorded here so the run order is not read as a clean bill of health. Details and
severity ranking are in the audit report.

1. **Fig 3 / Fig 4 (attribution) in the manuscript directory is a stale RENDERING —
   the numbers are intact.** `chapter1/manuscript/figures/article_v323/Fig3_attribution.png`
   (md5 `97f0c89f…`) does not match the current, deterministic output of its own
   generator `scripts/c1_fig3_units.R` (md5 `aa587326…`). However every value the
   generator prints reproduces the manuscript §3.2 text exactly (P4 leaf area
   −0.264 / cover −0.266, P1 0.049 / −0.021, |Hmax| ≤ 0.004, P4 asymmetry
   0.264 vs 0.224, profile −0.015 / −0.086 / −0.136 / +0.082; hot-day P4
   −0.356 as in the Fig. S2 caption). Stage C refreshes the rendering.
2. **`scripts/c1_operating_point_main.R:9`** reads the older per-part directory
   `out_files/Chapter1/tables/sensitivity_perplot_units/part_*.csv`, not the
   canonical `sensitivity_perplot_units_v2.csv` that Fig 3 uses.
3. **`c1_fig5_forward_2metric.R:8`** and **`scripts/make_residual_vs_topo.R:26`**
   read `tab_hobo_frblo_validation.csv` — the superseded **frblo** lineage.
4. **`scripts/make_trait_collinearity.R:19`** reads
   `transfer/musica_version_benchmark/run_simulations/inputs/fcover_per_plot.csv`.
5. **`c1_dumbbell_by_cluster.R:12`** reads `out_files/musica_hobo_z05_station_iter/1111`
   and a non-canonical forcing — the **z05** lineage. Not used by the current figures.
6. **Table 1 of the SHORT manuscript** is inline markdown with no generating script
   (its values are reproducible from `clhs_sample_native20_floor05.rds` +
   `R/cluster_relabel.R` and `sensitivity_perplot_units_v2.csv`, but nothing emits it).
7. **Appendix F collinearity: RESOLVED 2026-07-31.** The figure is now generated by
   `scripts/make_trait_collinearity_native20.R` (stage B20) from the canonical
   `clhs_sample_native20_floor05.rds`, in five panels as its caption promises, and the
   manuscript quotes exactly what that script prints (LAI-fCover 0.94 pooled; per archetype
   P1 0.79, P2 0.97, P3 0.93, P4 0.95). The former note claiming this blocked submission,
   and the three-way disagreement between shipped figure, ledger script and text, are retired.
8. **`c1_pca_recover_clusters.R`** is not deterministic (no seed; `ggrepel`),
   and reads `clhs_sample_floor05_v2.rds` rather than the native20 sample.
