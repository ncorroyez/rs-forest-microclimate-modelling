# Chapter 3 — Dynamic mapping of microclimatic attenuation

## Scientific objective

Chapter 3 tests whether injecting a **Sentinel-2-derived LAI time series**
into MuSICA improves simulation of sub-canopy thermal buffering at HOBO
sensors in Blois (43 plots, summer 2021), and which LAI correction method
best explains the observed microclimate.

The central hypothesis is:

> **H2 — Shape stability**: The vertical LAD profile from LiDAR is temporally
> stable. The seasonal canopy can therefore be described as
> `LAD(z, t) = LAD_static(z) × (LAI_S2(t) / LAI_ALS)`

---

## Experimental design — 15 scenarios, 3 axes

| Axis | Levels |
|------|--------|
| **LAI source** | ALS, ALS_DOPT, S2_ATBD, S2_DOPT, S2_RESCALED, RF |
| **Phenology** | STATIC (`calc_phenology`) · DYN (S2 time series) |
| **LAD shape** | full depth · d_opt-truncated (cross-paired for DOPT) |

| Scenario name | LAI | LAD | Phenology | Scientific question |
|---|---|---|---|---|
| `STATIC_ALS` | LiDAR full | full | parametric | Chapter 1 reference |
| `STATIC_ALS_DOPT` | LiDAR d_opt | full | parametric | d_opt ALS amplitude |
| `STATIC_S2_ATBD` | S2 ATBD | full | parametric | raw S2 bias impact |
| `STATIC_S2_DOPT` | S2 d_opt | full | parametric | d_opt correction alone |
| `STATIC_S2_RESCALED` | S2 rescaled | full | parametric | amplitude fusion |
| `STATIC_RF` | RF-corrected | full | parametric | data-driven correction |
| `DYN_ALS` | LiDAR full | full | S2 ATBD ts | H2: shape stable, timing from S2 |
| `DYN_S2_ATBD` | S2 ATBD | full | S2 ATBD ts | raw S2 + temporal |
| `DYN_S2_DOPT` | S2 d_opt | full | S2 dopt ts | d_opt + temporal |
| `DYN_S2_RESCALED` | S2 rescaled | full | S2 rescaled ts | fusion + temporal |
| `DYN_RF` | RF | full | S2 RF ts | RF + temporal |
| `STATIC_S2_DOPT_LADOPT` | S2 d_opt | d_opt | parametric | coherent d_opt pairing |
| `DYN_S2_DOPT_LADOPT` | S2 d_opt | d_opt | S2 dopt ts | full d_opt coupling |
| `STATIC_S2_ATBD_LADOPT` | S2 ATBD | d_opt | parametric | LAD shape alone |
| `DYN_ALS_LADOPT` | LiDAR full | d_opt | S2 ATBD ts | decouple 3 contributions |

The last two scenarios isolate individual contributions: `STATIC_S2_ATBD_LADOPT`
tests if LAD shape truncation corrects without amplitude fix; `DYN_ALS_LADOPT`
decouples amplitude (LiDAR), vertical shape (d_opt), and temporal shape (S2).

---

## Two-script workflow

The chapter is split into two independent scripts to separate the
**LAI correction stage** (Chapter 2 knowledge) from the **microclimate
simulation stage** (MuSICA). Run them in order:

| Script | What it does | When to re-run |
|--------|-------------|----------------|
| `Chapter2bis_lai_corrections.R` | Loads all LAI sources, runs RF LOO-CV, smooths S2 time series, produces 5 diagnostic figures. Writes `lai_prep/` artefacts. | When LAI inputs change or `FLAGS_PREP` are modified. |
| `Chapter3_main.R` | Reads `lai_prep/`, builds 15 scenarios, runs MuSICA, validates, maps. | When running or extending simulations. |

Shared parameters live in `Chapter3_config.R`, sourced by both scripts.

```
z_Example_rmusica_31012025/
├── Chapter3_config.R               ← shared config (paths, d_opt, dates) — source first
├── Chapter3_01_lai_corrections.R   ← Step 1: LAI corrections + diagnostics
├── Chapter3_main.R                 ← Step 2: MuSICA + validation + mapping
├── R/
│   ├── lai_corrections.R     ← LAI loading, correction, LOO-CV, sensitivity
│   │     assemble_s2_ts_stack()      build_lai_correction_table()
│   │     compute_lai_s2_rescaled()   loo_cv_rf()
│   │     train_lai_rf()              compute_sensitivity_expectation()
│   │     build_all_s2_ts()           load_lai_prep()
│   ├── spatial_mapping.R     ← GAM emulator + domain mapping
│   │     fit_microclimate_emulator() predict_dtmax_domain()
│   │     map_correction_impact()     plot_correction_impact()
│   │     load_domain_rasters()       pick_heatwave_days()
│   ├── phenology_s2.R        ← make_phenology_from_s2(), make_phenology_fn_factory()
│   ├── scenarios_c3.R        ← make_all_scenarios_c3() (15 scenarios)
│   ├── lad.R                 ← make_lad_dopt(), make_lad_dopt_factory()
│   ├── musica.R              ← run_musica_one() + phenology_fn extension
│   ├── validation.R          ← validate_scenarios_at_hobos()
│   ├── io.R                  ← load_lidar_rasters(), read_hobo_daily()
│   └── sentinel.R            ← load_s2_ts(), smooth_s2_ts()
└── out_files/Chapter3/
    ├── lai_prep/             ← Step 1 outputs (df_plots_lai.rds, ts_by_plot.rds, …)
    ├── figs/                 ← PNG figures + GeoTIFF maps
    ├── nc/                   ← MuSICA NetCDF outputs per scenario
    └── tables/               ← metrics, sensitivity expectation CSVs
```

---

## Data requirements

### Already available in `in_files/`
| File | Used for |
|------|---------|
| `lai_z1_res_10_m.tif` | LAI_ALS |
| `s2lai_summer_atbd_res_10_m.tif` | LAI_S2_ATBD (summer snapshot) |
| `hmax_p95_res_10_m.tif` | Hmax |
| `fCover_res_10_m.tif` | fCover |
| `vci_res_10_m.tif`, `lcv_res_10_m.tif` | RF covariates |
| `musica_in_Blois.nc` | ERA5 forcing |
| `data_Blois_utm31n.geojson` | HOBO locations |
| `Blois_data_temperature.csv` | HOBO observations |

### Read-only from NC_Full (no copy needed)
| Path | Used for |
|------|---------|
| `03_RESULTS/Blois/Metrics/Deciduous_Only/s2lai_2021-*_atbd_res_10_m.tif` | S2_ATBD time series |
| `03_RESULTS/Blois/Metrics/Deciduous_Only/s2lai_2021-*_optim_Blois_res_10_m.tif` | S2_DOPT time series |
| `revision/output/intermediate/lai_als_dopt/Blois/LAI_ALS_dopt_per_site.tif` | LAI_ALS_DOPT |
| `revision/output/intermediate/sm6/Blois/s2lai_summer_opt_per_site_res_10_m.tif` | LAI_S2_DOPT |

---

## Running the pipeline

### Step 1 — LAI corrections and diagnostics (`Chapter3_01_lai_corrections.R`)

No MuSICA required. Run once before simulations.

```r
FLAGS_PREP$RELOAD_LAI  <- TRUE   # force recompute even if cache exists
FLAGS_PREP$RELOAD_TS   <- TRUE   # reprocess S2 time series
FLAGS_PREP$RUN_LOO_CV  <- TRUE   # RF LOO-CV (takes ~2 min)
FLAGS_PREP$SAVE_FIGURES <- TRUE
```

Produces (in `out_files/Chapter3/lai_prep/`):
- `df_plots_lai.rds` — per-plot table with all 6 LAI columns
- `ts_by_plot.rds` — smoothed S2 time series, 4 correction variants
- `rf_loo_cv.rds` — RF LOO-CV predictions + metrics

Figures (in `out_files/Chapter3/figs/`):
- `c3_lai_distribution.png` — LAI density per correction
- `c3_rf_loo_cv.png` — RF vs raw S2 vs LiDAR scatter
- `c3_phenology_sample.png` — seasonal curves × 4 sample plots
- `c3_lai_correction_map.png`, `c3_lai_scatter_all.png` — spatial diagnostics

### Step 2a — Dry run (`Chapter3_main.R`, 3 plots × 15 scenarios)

Verify MuSICA is callable and check per-plot runtime (~5 min for 45 runs).

```r
FLAGS_C3$DRY_RUN         <- TRUE
FLAGS_C3$RUN_SIMULATIONS <- TRUE
FLAGS_C3$RUN_SPATIAL_MAP <- FALSE
```

### Step 2b — Full run (`Chapter3_main.R`, 43 plots × 15 scenarios)

Expect ~3 h on a standard workstation.

```r
FLAGS_C3$DRY_RUN         <- FALSE
FLAGS_C3$RUN_SIMULATIONS <- TRUE
FLAGS_C3$RUN_SPATIAL_MAP <- TRUE
```

---

## Key design decisions

### RF validation — no data leakage

The RF correction is trained on HOBO plots and then used to drive MuSICA at
the same locations. To avoid data leakage:
- `loo_cv_rf()` runs before any MuSICA simulation
- LOO-CV R² and RMSE for RF vs raw S2 are reported upfront
- These metrics bound the maximum gain the RF scenario can legitimately claim

### d_opt temporal assumption

`LAD_DOPT` truncates the LAD to `d_opt = 6 m` from the canopy top, derived
from the Pareto criterion on summer data (DSM normalization). This depth is
held constant through the season. In spring/autumn with partial canopy fill,
S2 may see deeper — this is a documented simplification.

### Sensitivity expectation

Set `CH1_SENSITIVITY` (°C per m²/m² LAI) from the Chapter 1 GAMM marginal
effect of LAI on ΔTmax. The `compute_sensitivity_expectation()` function then
translates each LAI correction's mean delta-LAI into an expected RMSE gain
upper bound. If the observed gain in simulations matches this bound, the
correction's impact is physically explained; if not, other factors dominate.

### Assumption about MuSICA's internal scaling

H2 is cleanly implemented only if MuSICA scales the allometry at each timestep
by `phenology(t) / lai_max_per_cohort`. If MuSICA renormalises the phenology
curve internally, the `LAI_S2(t) / LAI_ALS` factor is altered. Verify against
MuSICA source before interpreting STATIC vs DYN scenario deltas.

---

## Expected outputs

| File | Description |
|------|-------------|
| `cache/lai_correction_table.rds` | Per-plot LAI table with 6 sources + RF model |
| `cache/s2_ts_all.rds` | Smoothed S2 time series (4 corrections × 43 plots) |
| `cache/rf_loo_cv.rds` | LOO-CV predictions and metrics |
| `tables/c3_metrics_all.csv` | RMSE, R², Bias, MAE for all 15 scenarios |
| `tables/c3_sensitivity_expectation.csv` | Theoretical RMSE gain per correction |
| `figs/c3_rf_loo_cv.png` | LOO-CV scatter: RF vs raw S2 vs LiDAR |
| `figs/c3_lai_distribution.png` | LAI distribution by correction |
| `figs/c3_phenology_sample.png` | Seasonal LAI curves by correction (4 plots) |
| `figs/c3_rmse_all.png` | RMSE bar chart (15 scenarios) |
| `figs/c3_scatter_matrix.png` | Obs vs sim ΔTmax scatter (faceted) |
| `figs/c3_seasonal_bias.png` | Seasonal bias per scenario (MAM/JJA/SON) |
| `figs/c3_best_h2h.png` | Best static vs best dynamic head-to-head |
| `figs/c3_map_{date}.png` | Predicted ΔTmax + correction impact map (3 heatwave days) |
| `figs/c3_dtmax_pred_{scenario}_{date}.tif` | GeoTIFF for GIS |
| `figs/c3_correction_impact_{date}.tif` | Correction impact GeoTIFF |

---

## Physical interpretation guide

**Static scenarios**: If `STATIC_S2_DOPT` beats `STATIC_S2_ATBD`, the d_opt
amplitude correction improves MuSICA even with a fixed phenology — the
summer-peak LAI bias matters.

**Static vs Dynamic**: If `DYN_*` systematically beats `STATIC_*` in
spring/autumn but not summer, it confirms that S2 phenological timing adds
value beyond amplitude correction alone.

**LAD_DOPT scenarios**: If `DYN_S2_DOPT_LADOPT` beats `DYN_S2_DOPT`, the
vertical structure truncation contributes beyond the amplitude fix. The
comparison `STATIC_S2_ATBD_LADOPT` vs `STATIC_S2_ATBD` isolates the LAD
shape contribution alone.

**`DYN_ALS_LADOPT`**: If this matches or beats `STATIC_ALS`, shape stability
(H2) holds — the temporal variation is captured by S2 without needing a
corrected LAI amplitude. If it loses to `STATIC_ALS`, the seasonal LAI
variation from S2 introduces noise rather than signal.

**Spatial maps**: Correction impact should be strongest in areas where
LAI_S2_ATBD diverges most from LAI_S2_DOPT — typically dense canopy where
S2 underestimates depth. These areas should also show the largest thermal
buffering in the reference map.
