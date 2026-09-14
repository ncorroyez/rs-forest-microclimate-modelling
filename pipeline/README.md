# Chapter 1 — clean reproducible pipeline

End-to-end, idempotent pipeline from LiDAR data to the attribution & validation
figures. Each stage reuses the vetted `R/` library and is skip-if-exists.

```
Rscript pipeline/run_all.R                       # full spine
PIPE_RECLUSTER=TRUE Rscript pipeline/run_all.R   # opt-in re-clustering check
Rscript pipeline/05_shapley.R                    # any stage alone (re-sources 00)
```

## DAG (default spine)

| Stage | Does | Key inputs | Outputs |
|---|---|---|---|
| `00_config` | conventions, paths, asserts (single source of truth) | `R/config.R` | `PIPE` list |
| `01_data` | LiDAR rasters → forest pixels → LAD-shape FPCA | `in_files/*.tif` | `data/forest_fpca.rds` |
| `02_clusters` | **frozen** P1–P4 typology + archetypes | `clhs_sample_floor05_v2.rds` | `data/archetypes.rds` |
| `03_musica` | HOBO × 16 coalitions (mean fCover baseline), legacy binary | cluster sample, rasters | `out_files/musica_hobo_v10_fcovmean/<bit>/` |
| `04_extract` | per-plot×coalition metrics (2 conventions, below) | stage 03 NCs, HOBO csv | `data/{coal_metrics,ref_validation,clusters}.rds` |
| `05_shapley` | exact Shapley φ per plot (ΔTmax, ΔVPDmax) | `coal_metrics` | `tables/tab_shapley_perplot.csv` |
| `05b_shapley_ci` | bootstrap 95% CI of median φ (robustness, NOT p-values) | `shapley_perplot` | `tables/tab_shapley_ci.csv` |
| `06_boxplots` | Shapley boxplots — 2 single figures, each with BOTH metrics (free-y per metric): `fig_shapley_pooled` (metric×period) + `fig_shapley_by_cluster` (metric×period rows × cluster cols, P1–P4 + All) | `shapley_perplot` | `fig_shapley_{pooled,by_cluster}` |
| `06b_ranking` | trait ranking by \|median φ\| — ordered Cleveland dot plot (median ± bootstrap CI, coloured by buffer/amplify/ns) + table | `tab_shapley_ci` | `fig_shapley_ranking`, `tab_shapley_ranking.csv` |
| `07_validation` | REF vs 53 HOBO (glass ceiling) | `ref_validation` | `fig_validation_*` |
| `07b_forward_validation` | forward validation in Shapley-rank order, per cluster + All | `coal_metrics`, `tab_shapley_ranking` | `fig_validation_forward_shapley`, `tab_validation_forward.csv` |
| `08_typology` | FPCA modes + cluster LAD profiles (cluster palette, no legend) + elbow + archetypes + map; per-cluster structure table | `forest_fpca`, frozen sample | `fig_fpca_*`, `fig_clusters_lad_profiles`, `tab_cluster_structure.csv`, … |
| `09_extra_metrics` | **ANNEX** — Shapley of ΔTmin / diurnal amplitude / T-stability (temperature-only, all-summer) | stage-03 NCs | `annex/fig_annex_shapley_extra_by_cluster`, `tab_annex_shapley_extra_ci.csv` |
| `10_conditional_shapley` | **ANNEX** — conditional Shapley for LAD (φ_LAD \| Hmax=real vs baseline): bounds the off-manifold artefact of correlated traits | `coal_metrics`, `shapley_perplot` | `annex/fig_annex_conditional_shapley_LAD`, `tab_annex_conditional_shapley_LAD.csv` |

Outputs land in `outputs/figures_pipeline/{,/tables,/data}`.

### cLHS attribution figures (article Fig 2 / 2b / 3)

| Stage | Does | Key inputs | Outputs |
|---|---|---|---|
| `11_clhs_attribution` | **Fig 2** (per-cluster mean\|φ\|) + **Fig 3** (global mean\|φ\|), canonical & LOCAL | floor05_v2 Shapley parts | `out_files/Chapter3/figures/FigSh_shapley_clhs_{cluster,global}.png` |

The 6400 cLHS sims (400 plots × 16 coalitions) are computed by `c3_shapley_chunk.R`
(per-cluster baseline) and `c1_shapley_global_chunk.R` (global baseline), both on the
**ARTICLE-canonical `clhs_sample_floor05_v2.rds`**, caching NetCDFs in
`out_files/Chapter3/nc_shapley2x[_global]/` and per-plot φ in
`out_files/Chapter3/tables/shapley_parts[_global]/`. Fig 2b (metric × period heatmap)
is `c1_metrics_chunk.R`/`c1_metrics_merge.R` from the same cached NetCDFs.
`c3_shapley_clhs_2x.R` (old `floor05` sample) is **DEPRECATED** (hard `stop()`), the
source of the former Fig 2/3 provenance ambiguity. All attribution outputs are LOCAL;
nothing writes to `~/Documents/NC_Full`.

## Frozen decisions (do NOT drift — see `00_config.R::PIPE`)

1. **Clustering is a frozen input.** P1–P4 come from `clhs_sample_floor05_v2.rds`.
   `02_clusters` never re-clusters by default (k-means random init would shift
   P1↔P4 and break every downstream figure). `PIPE_RECLUSTER=TRUE` re-runs with a
   pinned seed and asserts ARI > 0.95 vs the cache.
2. **fCover baseline = mean (0.869)** (Eva/v10), NOT 1 (v9). The 16-coalition set
   is self-contained in `musica_hobo_v10_fcovmean/` (stage 03 completes it).
3. **Never reuse `DT_daily_HOBO_floor05_*`** — it carries a third (validation)
   baseline (0000 → ΔTmax 1.28 °C), neither v9 nor mean.
4. **All metrics at fixed 1 m** (constraint #4 — `SLOPE_USE_NAIR1=FALSE`):
   - absolute ΔTmax / ΔVPDmax → fixed **1 m** interp + **−2 h** shift (HOBO-comparable);
   - relative micro/macro slope → **also fixed 1 m**, **no** shift.
   The legacy **nair==1 / 45-8 / r=0.92** slope convention is **superseded** (do not cite it).
5. **Idempotent / reuse-by-default.** No stage re-runs MuSICA unless an NC is missing.

## Regression acceptance (asserted in-pipeline)

- `05_shapley`: per-plot additivity Σφ = ΔTmax(REF) − ΔTmax(Base); **median Σφ(ΔTmax,all) ≈ −0.311 °C** (legacy branch).
- `07_validation` is **branch-aware** (`PIPE$ASSERT$VAL`), enforced on both branches:
  - **z05** (PUBLISHED, manuscript Fig 4/5): per-plot ΔTmax **r ≈ 0.93**; split **44/9**.
  - **legacy v10/z1** (historical guard): **r ≈ 0.92**; split **45/8**.

The published validation is the **z05** branch (`PIPE_BRANCH=z05`, understory-aware 0.5 m LAD).
Any drift in these = a convention or the clustering moved.

## Optional analyses (not in `run_all`)

```
Rscript pipeline/optional/10_canopy_ratio.R       # CR vs buffering (answers Bouwen)
Rscript pipeline/optional/11_vertical_profiles.R  # wind / RH / VPD profiles by cluster
Rscript pipeline/optional/12_height_robustness.R  # nair==1 vs 1 m robustness
```
