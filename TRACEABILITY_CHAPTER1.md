# TRACEABILITY — Chapter 1

**Target manuscript:** `chapter1/manuscript/manuscript_chap1_EN_native20.md` (58 p., 25 figures,
9 appendices A–I, Table 1 + Tables G1, G2, H1, H2).
The SHORT version (`manuscript_chap1_EN_native20_SHORT.md`) is **frozen and out of scope**; where a
short-doc figure label is given it is a cross-reference only.

**Purpose.** For any number, figure or table in the chapter, find the file that produced it, and follow
the chain back to a raw input. This complements — it does not replace — `RUN_CHAPTER1.md` (how to
*re-run* the chain) and `chapter1/ledger/make_article_figure_set.R` (the figure presence ledger).

**How to read this document**

| you want | go to |
|---|---|
| where a raw input comes from | §1 Input inventory |
| which script made a figure or table | §2 Artifact → chain |
| which script made a number in an appendix | §3 Appendix-by-appendix |
| what must never be deleted | §4 Simulation assets |
| how to go from a checkout to the PDF | §5 Start here |
| what could not be traced | §6 Not traceable |
| what is broken or fragile | §7 Broken and fragile links |

**Canonical lineage** (anything outside this is superseded and must not appear in a Chapter 1 number):

| item | canonical value |
|---|---|
| forcing | `in_files/FR-Blo_2021_v2.nc` (CHS 41 station based; solar time = UTC+1) |
| ΔTmax convention | `R/dtmax_convention.R`, time-matched, `OBS_CLOCK_OFFSET_H = -1` on **observations only** |
| design sample | `out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds` (400 real 20 m pixels) |
| archetype labels | `out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv`, column `P`, **frozen** |
| model | MuSICA v3.2.3, iterative ABL ("yoyo"), canopy wind log-profile correction ON |

Superseded lineage markers, never to be cited: `frblo`, `z05`, `era5sd`, `v320`, `floor05_v2`,
`convA`, `fpc`, `shapley_*`, `transfer/`, `outputs/figures_pipeline_z05/`.

---

<!-- SECTION 1 -->
## 1. Input inventory

Everything the chain consumes that is not itself produced by the chain. Location classes:

* **repo** — inside `/home/corroyez/Documents/z_Example_rmusica_31012025`
* **drive** — on the external disk `/media/corroyez/MyPassport` (verified **mounted** 2026-07-31)
* **external-local** — elsewhere on this machine, outside the repo and outside the drive.
  These are the dangerous ones: they are in no backup the repo controls.
* **credential** — a secret required to re-download, never to be printed

### 1.1 Airborne laser scanning — the raw structural input

| path | what it is | origin | location | read by |
|---|---|---|---|---|
| `/media/corroyez/MyPassport/01_DATA/Blois/LiDAR/leaf_on/1-las_l93` | raw tiles as delivered, Lambert-93 (160 files) | ALTOA survey, 17 June 2021, RIEGL VQ-780i | **drive** | nothing in the Chapter-1 chain (historical only) |
| `…/leaf_on/2-las_utm` | filtered/classified tiles reprojected to UTM (320 files, ~149 G) | derived by the provider + reprojection | **drive** | `c1_scanangle_raw_stats.R` (A13); `scripts/compute_lad_z05_sites.R`; `scripts/control_scanangle_vs_lidr.R` |
| `…/leaf_on/3-las_normalized_utm` | height-normalized against ground returns (160 files) | provider-normalized | **drive** | `c1_fig1_gril_native.R` (B1) — **only if the cache is deleted** |

Survey metadata as stated in §2.1.2: Partenavia P68 at 900 m, 2000 kHz, ±30° scan angle, 72.6 pts m⁻²
(9.7 ground returns m⁻²), 4.5 cm vertical precision. Provenance is the provider's, not reconstructible
from the repo.

### 1.2 Derived trait rasters

| path | what it is | origin | location | read by |
|---|---|---|---|---|
| `in_files_native20/{lai_z1,fCover,max,vci,sec_theta}_res_10_m.tif`, `lad_profiles_z1_res_10_m.tif` | **the canonical native-20 m trait rasters** | computed natively at 20 m from the ALS returns, 2026-07-17 | **repo** | Fig 8, Fig D1, Fig S5, `c1_cover_floor_check.R` |
| `in_files/{lai_z1,fCover,max,hmax_p95,lcv,lskew,lidarlai}_res_10_m.tif`, `ladstack*.tif` | the **older 10 m** products | earlier lineage | **repo** | superseded — not on the Chapter-1 chain |
| `out_files/clusters/blois_clusters_P1P4_native20.tif` | the archetype map | K-means on 6 metrics (§2.3) | **repo** | Fig 1a |

> **Trap.** The native20 rasters keep the filename suffix `res_10_m` but are **20 m** and in a
> **different CRS**: `in_files_native20/lai_z1_res_10_m.tif` is 20 m, EPSG:**25831** (ETRS89/UTM31N),
> 419×320; `in_files/lai_z1_res_10_m.tif` is 10 m, EPSG:**32631** (WGS84/UTM31N), 817×619. Same name,
> different grid and different datum. Verified 2026-07-31.

### 1.3 LAD profile CSVs

| path | what it is | origin | location | read by |
|---|---|---|---|---|
| `in_files/lad_z05/Blois_lad_z05_r25.csv` | per-plot LAD profiles at 0.5 m bins, 25 m circular footprints, 60 field plots | `scripts/compute_lad_z05_sites.R` reading `2-las_utm` on the **drive** | **repo** (product) | Fig 8, Fig E2, Fig F1, Fig S5, Fig C1 (coords only) |
| `in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv` | the same, scan-angle corrected; carries the per-plot `sec_theta` | same script with `SEC_CORR = TRUE` | **repo** (product) | Fig F1; `c1_hobo_native20_validation.R`; `c1_forward_windcorr_sims.R` |
| `in_files/lad_z05/{Aigoual,Mormal}_*.csv` | other sites | same script | **repo** | not Chapter 1 |
| `allometry_profiles.csv` (repo root) and `in_files/allometry_profiles.csv` | allometric profile table | legacy | **repo** | **not consumed by the Chapter-1 chain** |

`scripts/compute_lad_z05_sites.R` carries its own warning at L26: *"UNTESTED: the raw LAS (MyPassport
drive) were offline when this was written."* Regenerating these CSVs requires the drive.

### 1.4 HOBO logger data

| path | what it is | origin | location | read by |
|---|---|---|---|---|
| `in_files/Blois_data_temperature.csv` (157 MB) | the raw hourly logger series, `CFG$hobo_temp_csv` | field campaign summer 2021, Onset HOBO Pendant UA-001-64 | **repo** | Fig S1; every validation script via `R/validation.R` |
| `in_files/data_Blois_utm31n.geojson` | logger coordinates | Trimble Geo 7X differential GPS | **repo** | Fig 1, Fig 8, Fig D1, Fig S5, `c1_scanangle_raw_stats.R`, `c1_cover_floor_check.R` |
| `in_files/00_readme_temperature.txt` | field notes for the above | field campaign | **repo** | documentation only |

**60 loggers were simulated; 53 are used.** `R/config.R:19` drops seven:
`41_13, 41_14, 41_20, 41_41, 41_50, 41_51, 41_53`. This is why every simulation directory holds 60
NetCDFs while every reported statistic is over 53. `CFG$date_seq` = 2021-06-01 → 2021-09-30.

### 1.5 Meteorological forcing

| path | what it is | origin | location | read by |
|---|---|---|---|---|
| **`in_files/FR-Blo_2021_v2.nc`** (448 kB) | **the canonical forcing.** Vars: `CO2air, LWdown, PSurf, Qair, Rainf, SWdown, Snowf, Tair, Wind_E, Wind_N, h_sbl, level, nav_lat, nav_lon`. Hourly, `hours since 2021-01-01`, 8759 steps, timestamped at HH:30 | **see Finding 1-0 — not establishable from the repo** | **repo** | **every** simulation, and 34 files in total; Fig 6, Fig A1, Fig I1, Fig S1 |
| `in_files/musica_in_Blois.nc` | the **older ERA5-based** forcing | earlier lineage | **repo** | Fig A1 (as the ERA5 comparison series); **Fig B1's cached sims** (see §3 Finding B-i) |
| `in_files/MetHor2021.txt`, `MetHor2022.txt` (9.3 MB each; also at repo root) | the raw CHS 41 / RENECOFOR station record | ONF RENECOFOR | **repo** | forcing assembly, not the figure chain |
| `ONF RENECOFOR - Synthese_climatique_CHS41.pdf` | station documentation | ONF | **repo** | documentation only |
| `prep_site_forcing/MuSICA_in_CHS41-Blois_2021-station-era5.nc` | intermediate station+ERA5 build | `prep_site_forcing/prep_site_forcing.py` | **repo** | `c1_inject_hsbl_station.R` |
| `out_files/Chapter1/merra2_pblh.csv` (345 kB) | MERRA-2 PBLH, **hourly**, 13 848 rows, 2021-01 → 2022-07, range 58.7–3778 m, mean 767 m | `get_merra2_pblh.sh` → NASA GES DISC `M2T1NXFLX.5.12.4`, `PBLH[0:1:23]` | **repo** | `build_forcing_pblh.R`, `c1_inject_hsbl_station.R` — **neither of which writes the canonical forcing** |
| `SAFRAN/` | SAFRAN reanalysis extracts | Météo-France | **repo** | superseded (`tab_hobo_safran_validation.csv` is off-chain) |
| `.cdsapirc` | Copernicus CDS API key | user account | **credential** | needed only to re-download ERA5; **contents never printed** |

> **Finding 1-0 — the canonical forcing's own origin is the single most load-bearing unresolved
> link in the chapter.** `in_files/FR-Blo_2021_v2.nc` drives *every* simulation. Its own NetCDF global
> attributes, read directly, say:
>
> ```
> :file_name  = "/home/rpatin/Documents/MuSICA/sites/Blois/in_files/FR-Blo_2021-2022.nc"
> :production = "NetCDF file generated from ASCII file
>                /home/rpatin/Documents/MuSICA/sites/Blois/in_files/FR-Blo_2021-2022.csv
>                on Mon 06 July 2026 11:47:16"
> :history    = "Mon Jul  6 14:04:58 2026: ncks -O -d time,0,8758
>                in_files/FR-Blo_2021-2022_2.nc in_files/FR-Blo_2021_v2.nc"
> ```
>
> So the file was assembled **by a third party (R. Patin) from a CSV that is not on this filesystem**,
> and the repo only performed the final `ncks` time-trim. Three in-repo claims about what is inside it
> disagree:
> * `RUN_CHAPTER1.md:15` and manuscript §2.4 — "CHS 41 station based"
> * `out_files/Chapter1/figures/cmp_era5_frblo/README_comparaison_era5_frblo.md:14` — "site flux **ICOS** de Blois"
> * its own attributes — an upstream CSV of unstated composition
>
> The same README reconciles the first reading empirically (FR-Blo `Tair` vs the CHS 41 station:
> *r* = 1.000, mean difference +0.00 °C), which is strong evidence the temperature channel really is
> the station. But that is an in-repo assertion checked against an in-repo file, not the raw record.
> **The composition of the radiative and aerodynamic channels is not verifiable from this repository.**
> Appendix A and Table H1 both make specific claims about them.
>
> **Finding 1-i — the `h_sbl` in the canonical forcing has no traceable producer.**
> `in_files/FR-Blo_2021_v2.nc` carries `h_sbl` with n = 8759, min 1.20 m, max 244.25 m, mean 60.81 m,
> sd 47.33, 8758 distinct values — so it varies, it is not a constant placeholder. But:
> * Its magnitude is an order of magnitude below a planetary boundary-layer height (the archived
>   MERRA-2 series runs 59–3778 m, mean 767 m).
> * It does **not** reconcile with that archived series. Datetime-aligned, cor = 0.53; a ±3 h lag sweep
>   peaks at cor = 0.537 (+1 h). A pure rescaling would give cor = 1.000. Median ratio PBLH/h_sbl ≈ 10.2–10.5
>   and max |PBLH/10 − h_sbl| = 346 m.
> * **No script in the repo writes it.** `build_forcing_pblh.R` writes `in_files/musica_in_Blois_pblh.nc`
>   and `c1_inject_hsbl_station.R` writes `out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl.nc`. Neither
>   targets `FR-Blo_2021_v2.nc`.
>
> The ~1/10 magnitude is consistent with §2.4's *blending height* ("one tenth of the boundary-layer
> height"), which would make Table H1's row a mislabel rather than a wrong number — but that is
> inference, not provenance. See §6 item 1.

### 1.6 cLHS design sample

| path | what it is | origin | location |
|---|---|---|---|
| **`out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds`** | **the canonical design.** 400 × 50: `x, y, LAI, VCI, Hmax, fCover, LAD_Layer_1.5…39.5, FPC1-3, Cluster, sec_theta` | `scripts/redo_clhs_native20.R`; 100 plots per archetype, 10 000 iterations, fixed seed, fCover floored at 0.5 | **repo** |
| `out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds` and 12 further `clhs_sample_*.rds` | superseded draws | earlier lineages | **repo** — do not use |

`R/cluster_relabel.R` defines the frozen P1→P4 mapping (open → dense) applied everywhere. **Its header
comments are stale**: they quote per-archetype LAI as 2.4 / 6.0 / 7.7 / 10.8, whereas the native20
sample gives 1.80 / 3.15 / 3.70 / 5.03. The *mapping* is verified correct (mean LAI ascends strictly
P1→P4); only the comment misleads.

> **Finding 1-ii — the canonical design file is not what its nominal writer emits.**
> `scripts/redo_clhs_native20.R` clusters on **three** scalars (LAI, Hmax, fCover) and neither
> `build_forest_dataframe()` nor `sample_clhs_per_cluster()` (`R/forest.R`) emits FPC columns — yet the
> shipped file carries `FPC1-3`. Verified 2026-07-31:
>
> ```r
> all.equal(clhs_sample_native20_floor05.rds, clhs_sample_native20_fpca_floor05.rds)  # TRUE
> ```
>
> The canonical sample **is** the 6-metric FPCA design from `scripts/redo_clhs_native20_fpca.R`
> (clustering on LAI, Hmax, fCover, FPC1–3), copied onto the 3-metric filename: the `_fpca*` files are
> timestamped 2026-07-20 15:00 and the canonical names 15:09. **This matches the manuscript**, whose
> §2.3 and Fig 2 caption both describe an FPCA feeding the K-means typology — so the science is right
> and the filename is wrong. But `redo_clhs_native20.R`'s header still claims authorship, and the
> 3-metric design it describes is no longer on disk. The promotion step exists in no script.
>
> **Related timestamp caveat.** The shipped simulations straddle that 07-20 overwrite:
> `musica_hobo_native20/1111` (07-17 16:40), `musica_native20_forward` (07-17 17:57) and
> `radius_test_corr` (07-17 10:52) all predate it, while
> `nc_sensitivity_perplot_units` (07-28 17:42) postdates it. `c1_forward_native20_sims.R:14-15` builds
> its coalition scenarios from cluster means of this rds, so the forward sims were parameterised from a
> design state that is no longer on disk. The HOBO-location runs (`1111`) take their traits from the
> rasters rather than the design, so they are unaffected.

### 1.7 DEM and terrain covariates — Appendix C

| path | what it is | location | read by |
|---|---|---|---|
| `~/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked/dtm.tif` | elevation | **external-local** | `scripts/make_residual_vs_topo.R` (B11) |
| `…/slope_res_10_m.tif` | slope | **external-local** | same |
| `…/twi_res_10_m.tif` | topographic wetness index | **external-local** | same |
| `…/northness_res_10_m.tif` | northness | **external-local** | same |
| `…/std_res_10_m.tif` | canopy-height standard deviation | **external-local** | same |
| `…/rumple_res_10_m.tif` | rumple | **external-local** | same |

All six verified **present** 2026-07-31. `~/Documents/NC_Full` is **171 G**, lives outside the repo and
outside the external drive, and is named in **no** existing documentation — not `RUN_CHAPTER1.md`, not
the ledger, not the manuscript's data-availability section. It is the least-protected input in the
chapter. Worse, `make_residual_vs_topo.R:47` does `cov_files <- cov_files[file.exists(cov_files)]`, so
missing rasters are **silently dropped** and Appendix C's adjusted *R*² would change without erroring.

### 1.8 MuSICA binary and configuration

Three distinct binaries are present. **They are not interchangeable.**

| path | md5 | what it is | used for |
|---|---|---|---|
| `in_files/model-3.2.3/musica` | `8b5139e0dd300ed33c8398b0a0c1bd5d` | **v3.2.3 official — the Chapter-1 binary** | every native20 simulation. `c1_sensitivity_perplot_chunk_native20.R:19` overrides `CFG$musica_cmd` to point here |
| `/home/corroyez/Documents/musica/musica` | `530722784a60295287b1a7e080505815` | Nov-2024 legacy variant; **`CFG$musica_cmd`'s default** | the legacy Oct-2025 sims — **and Fig B1's cached controlled test** (§3 Finding B-i) |
| `./musica` (repo root) | `7ec7db0047dd6a0b8fbfff945f47686d` | a third, unidentified build | nothing traceable |

`R/config.R:11–16` documents the difference explicitly: substituting one for the other *"shifts ΔTmax by
up to 1.4 °C and inverts the buf/amp regime ratio"*.

**Namelists.** `musica.nml` sets `PARAM_FILES_PATH = "./in_files/Blois/in_files/"`, so the authoritative
configuration — the one Tables H1/H2 correctly transcribe — is:

* `in_files/Blois/in_files/musica_soil.nml`
* `in_files/Blois/in_files/musica_veg1.nml`

The repo root carries **divergent copies** of both (`musica_soil.nml` md5 `3e151155` vs `2778f1bc`;
`musica_veg1.nml` `c6f2daac` vs `84e8b07b`), differing in comments and in `CANOPY_HEIGHT_TOP`
(hard-coded `1.2` vs the `CANOPYHEIGHT` template placeholder that `rmusica` substitutes per plot).
`musica.nml:49` also still reads `FORCING_FILENAME = "in_files/musica_in_Safran_Blois_2021.nc"` — a
**SAFRAN** file, overridden at run time by `R/musica.R` but flatly contradicting Table H1 for anyone
reading the namelist directly.

### 1.9 Hand-authored inputs and R helper modules

| path | role |
|---|---|
| `scripts/c1_fig2_methodo.dot` | **the source of truth for Figure 2** — a hand-written graphviz spec, not derived from data |
| `scripts/_article_style.R` | shared palette, theme and figure sizing; sourced by nearly every figure script |
| `R/cluster_relabel.R` | the frozen P1→P4 archetype mapping |
| `R/dtmax_convention.R` | the ΔTmax definition: time-matched, `OBS_CLOCK_OFFSET_H = -1` on **observations only** |
| `R/config.R` | `CFG`: paths, `date_seq`, `ids_to_remove`, `musica_cmd`, `tair_target_height = 1.0` |
| `pipeline/00_config.R` | `PIPE`: `BITS` (the 16 coalitions) and re-exports of `CFG`. **A live dependency** despite `pipeline/` being superseded |
| `R/musica.R` | `run_musica_one()` / `run_musica_scenario()` — the MuSICA wrappers |
| `R/validation.R`, `R/io.R`, `R/lad.R`, `R/fpca.R`, `R/profiles.R`, `R/wind_correction.R` | extraction, I/O, LAD handling, the FPCA of §2.3, the wind log-profile |

The remaining ~60 files in `R/` (`h1_*`, `lovb_*`, `h2_analysis.R`, `gamm.R`, `sentinel.R`, …) are
earlier-thesis lineage and are **not** on the Chapter-1 chain.

### 1.10 Inputs whose origin could not be established

| item | why | how load-bearing |
|---|---|---|
| **`in_files/FR-Blo_2021_v2.nc` — the canonical forcing itself** | built by a third party from a CSV not on this filesystem; three inconsistent in-repo descriptions (Finding 1-0) | **highest** — drives every simulation |
| **`h_sbl` in that file** | no writer in the repo; does not reconcile with the archived MERRA-2 PBLH at any lag −3…+3 h (Finding 1-i) | high — a Methods claim (Table H1, §2.4) |
| **ERA5 surface roughness 0.44 m** | still a literal (`R/wind_correction.R`, `Z0_ERA_BLOIS`), but **documented 2026-07-31**: the header now records the exact CDS retrieval (variable `forecast_surface_roughness`, cell 47.50 N 1.25 E, hourly, 1 Jun–30 Sep 2021, averaged), states plainly that no fsr field is kept in `in_files/` so it cannot be recomputed here, and points at the z0,ERA sweep (0.01–1 m) in `c1_wind_profile_correction.R` / Fig. A2 that bounds its influence | high — propagates through the 38 cached wind-corrected forcings into **every** shipped native20 run, but its effect is now bounded rather than assumed |
| the monthly roughness range 0.438–0.443 m | exists only as a source-code comment | low (Appendix A prose) |
| the promotion of `clhs_sample_native20_fpca_floor05.rds` → `clhs_sample_native20_floor05.rds` | no script performs the copy (Finding 1-ii) | medium |
| `in_files/musica_in_Blois.nc` (the ERA5 comparison series) | `:production` names a macOS tmpdir ASCII file, **Sept 2022**; no ERA5 download on this machine produced it | medium — it is Fig A1's comparison arm and Fig B1's forcing |
| `out_files/Chapter1/tables/typical_plots_per_cluster.csv` | read by Fig 1b; no chain script writes it | low |
| the DEM's own origin | TWI / northness producers are **commented out** in `NC_Full/02_CODES/LiDAR/3_calculate_lidar_metrics.R`, and the DTM VRT they read is missing. Whether the DEM is ALS-derived or a national product is recorded nowhere | medium (Appendix C) |
| repo-root `./musica` (md5 `7ec7db00`) | a third binary of unstated build, referenced only by legacy `main.R:802` | low — not on the chain |
| ALS survey metadata (point density, precision, flight parameters) | provider-supplied; not reconstructible from the repo | low — descriptive only |
| the CHS 41 gap-filling against ERA5 | `MetHor*.txt` are raw; Appendix A describes the gap-fill but no chain script performs it | medium |
| raw SAFRAN `Forc*_V2.nc` (~11.8 GB) | no acquisition script; `get_safran.sh` presumes them present | none — SAFRAN is off-chain |
| `prep_site_forcing/CHS41-Blois.cfg` altitude (127 m) and reference height (1.50 m) | the cfg **self-flags** both as "best-guess, TO CONFIRM", plus an open `ftimestep` HH:00-vs-mid-step question | medium — 1.5 m is the ΔTmax reference height |

### 1.11 Broken paths referenced by code

All of these are in **legacy or off-chain** scripts — every live-chain path resolves. Listed so a reader
does not mistake a dead reference for a missing input.

| referencing file:line | path | fault |
|---|---|---|
| `in_files/Blois/main_Blois.R:27`; `script_rmusica_parallel_Blois.R:15`; `main_Blois_field.R:32` | `/media/corroyez/My Passport/01_DATA/Blois/LiDAR/{1-las_l93,2-las_utm}` | mount name has a **space** (`My Passport` ≠ `MyPassport`) *and* the `leaf_on/` level is missing |
| `build_forcing_pblh.R:21`, `get_merra2_pblh.sh:6` | `out_files/Chapter3/merra2_pblh.csv` | **missing** — only `out_files/Chapter1/merra2_pblh.csv` exists, so the Chapter-3 PBLH branch cannot be re-run as written |
| `c1_forward_heterogeneity.R:22` | column `gap_fraction` in `tab_residual_vs_topo.csv` | the column does not exist — its only producer deliberately writes `fCover` instead (`make_residual_vs_topo.R:62-63`); the covariate silently becomes `NA` |
| `NC_Full/02_CODES/LiDAR/3_calculate_lidar_metrics.R:155` (commented) | `NC_Full/03_RESULTS/Blois/LiDAR/dtm/res_10_m/rasterize_terrain.vrt` | missing — this is the DTM the Appendix C terrain derivatives came from |
| `Chapter1_refactored.R:94` | `bash -i -c musica` | `musica` is not on `PATH` |
| `SAFRAN/SAFRAN_extract_script.R:7,11` | `C:/Users/eva gril/ownCloud/...` | a Windows path from another person's machine, for a different site |
| `in_files/Blois/main_Blois.R:153` | `../Lidar_Louchats/MAST/CONV1_15m.las` | no such tree anywhere |
| `in_files/SAFRAN/wgs84/` | — | directory exists but is empty |

### 1.12 Housekeeping observations

* **`.cdsapirc` is root-owned, mode 644, inside the git working tree.** `/home/corroyez/.cdsapirc` is
  correctly `600`. Confirm the repo copy is git-ignored, or remove it. Contents were never read.
* **~23.6 GB of avoidable duplication**: `SAFRAN/` and `in_files/SAFRAN/init/` are size-identical
  copies of the same 14 files; `MetHor2021.txt` is byte-identical at the repo root and in `in_files/`
  (the R scripts read the root copy, the Python reads `in_files/`).
* **Two in-code warnings are now stale and can be retired**: `scripts/compute_lad_z05_sites.R:26`
  ("UNTESTED — drive offline") and the drive guards in `c1_scanangle_raw_stats.R` /
  `c1_fig1_gril_native.R`. The drive is mounted and all composed paths resolve.
* **Orphan data with no reader**: `in_files/FR-Blo_2021-2022{,_2}.nc`,
  `in_files/musica_in_FR-Bil_2022{,_with_ABL}.nc`, `in_files/musica_in_Blois_hsbl.nc`, `MetHor2022.txt`.

---

<!-- SECTION 2 -->
## 2. Artifact → chain

One row per artifact, keyed by the **full manuscript's** label. Every entry was verified by opening the
generating script and reading its actual `ggsave` / `fwrite` targets — not by trusting its comments or
the ledger. `~R` = repo root `/home/corroyez/Documents/z_Example_rmusica_31012025`.

The **md5** column compares the shipped file in `chapter1/manuscript/figures/article_v323/` against the
generating script's own output path. Read it as one story: every shipped file is dated 2026-07-30 14:51,
a single coherent `--sync-only` event, and every sync source is 14:48–14:51 **except**
`fig_scanangle_control.png` at 15:13. There is exactly one post-sync regeneration that was never re-synced.

Nearly every script also sources `~R/scripts/_article_style.R` (shared palette and theme) and many source
`~R/R/cluster_relabel.R` (the frozen P1–P4 mapping); these are omitted from the input column except where
they are the only input.

### 2.1 Figures

| artifact | ms line | generating script | inputs consumed | output(s) written | shipped file | md5 | stage |
|---|---|---|---|---|---|---|---|
| **Fig 1** typology | 111 | `~R/c1_fig1_gril_native.R` | `out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds`; `in_files/data_Blois_utm31n.geojson`; `out_files/Chapter1/tables/typical_plots_per_cluster.csv`; `out_files/clusters/blois_clusters_P1P4_native20.tif`; `in_files/fig1_crosssection_cache.rds` (cache; else the external drive); `tab_hobo_perplot_cluster.csv` | `outputs/figures_chap1/Fig1_typology_gril_native.png` | `Fig1_typology_gril_native.png` | MATCH `90bb8269` | B1 |
| **Fig 2** pipeline schematic | 117 | `~R/scripts/c1_fig2_methodo.R` (graphviz `dot`) | `~R/scripts/c1_fig2_methodo.dot` — **hand-authored, the source of truth** | `out_files/Chapter1/figures/Fig2_methodo.{png,pdf}` | `Fig2_methodo.png` | MATCH `a549711f` | B22 |
| **Fig 3** PCA archetypes | 192 | `~R/c1_pca_recover_clusters.R` | `clhs_sample_native20_floor05.rds`; `R/cluster_relabel.R` | `out_files/Chapter1/figures/FigSh_pca_clusters.png`; `out_files/Chapter1/tables/tab_pca_clusters.csv` | `Fig3_pca_recover_clusters.png` | MATCH `2ed37531` | B6 |
| **Fig 4** attribution | 208 | `~R/scripts/c1_fig3_units.R` (defaults `PERIOD=all METRIC=dt`) | `out_files/Chapter1/tables/sensitivity_perplot_units_v2.csv` | env-conditional base `out_files/Chapter1/figures/Fig3_attribution_units[_slope][_hot]` `.{png,pdf}` | `Fig4_attribution.png` | MATCH `aa587326` | B2 (data A3) |
| **Fig 5** operating point | 216 | `~R/scripts/c1_operating_point_main.R` | `sensitivity_perplot_units_v2.csv` | `outputs/figures_chap1/fig_operating_point_main.{png,pdf}` | `Fig5_operating_point.png` | MATCH `0b3075f7` | B3 |
| **Fig 6** glass ceiling | 224 | `~R/scripts/c1_fig5_j1_convB.R` | `in_files/FR-Blo_2021_v2.nc`; `tab_hobo_perplot_cluster.csv`; `tab_hobo_native20_validation.csv`; `out_files/musica_native20_forward/{0000,1000,1100,1110,1111}/*.nc`; all of `R/*.R` | `out_files/Chapter1/figures/Fig5_glass_ceiling_dumbbell.{png,pdf}`; `FigJ1_forward_inclusion.{png,pdf}`; `tab_forward_inclusion_convB.csv` | `Fig6_glass_ceiling_dumbbell.png` | MATCH `3c56dbb0` | B4 |
| **Fig 7** footprint radius | 247 | `~R/scripts/c1_fig7_footprint_radius.R` (data from `~R/c1_radius_rescore.R`) | `out_files/Chapter1/tables/tab_radius_sweep_aligned.csv` | `out_files/Chapter1/figures/Fig7_footprint_radius.{png,pdf}` | `Fig7_footprint_radius.png` | MATCH `691797a6` | B7 (data A5) |
| **Fig 8** observational corroboration | 267 | `~R/scripts/c1_fig6_obs_nested_native.R` | `in_files/data_Blois_utm31n.geojson`; `in_files_native20/{lai_z1,fCover,max,vci}_res_10_m.tif`; `in_files/lad_z05/Blois_lad_z05_r25.csv`; `tab_hobo_perplot_cluster.csv` | `out_files/Chapter1/figures/Fig6_obs_nested_native.{png,pdf}` | `Fig8_obs_corroboration.png` | MATCH `120bb98d` | B5 |
| **Fig A1** ERA5 station bias | 325 | `~R/c1_era5_station_bias.R` | `in_files/FR-Blo_2021_v2.nc`; `in_files/musica_in_Blois.nc` | `outputs/figures_chap1/B1_era5_station_bias.{png,pdf}` | `A1_era5_station_bias.png` | MATCH `b10ae6ec` | B8 |
| **Fig A2** wind log-profile | 346 | `~R/c1_wind_profile_correction.R` | **none — purely analytic** (the log-profile formula and hardcoded roughness levels) | `outputs/figures_chap1/G1_wind_profile_correction.{png,pdf}` | `A2_wind_profile_correction.png` | MATCH `202159ea` | B14 |
| **Fig B1** controlled H2 test | 378 | `~R/scripts/make_h2_controlled_topheavy.R` | `out_files/H2_controlled_topheavy/LAI{6,12}/*.nc` (14 cached sims); `CFG$forcing_file` = `in_files/musica_in_Blois.nc` (**ERA5**) | `outputs/figures_pipeline/annex/fig_h2_controlled_topheavy.{png,pdf}`; `tab_h2_controlled_topheavy.csv` | `B1_h2_controlled_topheavy.png` | MATCH `c645b64d` | B21 |
| **Fig B2** vertical gradient | 391 | `~R/scripts/make_vertical_profiles_native.R` **+ `convert -append` montage in `run_chapter1.R:173-186`** | `out_files/musica_native20_forward/{1111,1110}/*.nc`; `tab_hobo_perplot_cluster.csv` | `outputs/figures_chap1/fig_vertical_profiles_native_{norm,metres}.{png,pdf}`; `tab_subcanopy_native.csv` | `B2_vertical_gradient.png` | montage **pixel-identical** (`compare -metric AE` = 0); bytes differ by PNG encoder metadata only | B9 + sync |
| **Fig B3** per-trait gradient | 393 | `scripts/make_trait_vertical_gradient_perplot.R` (per-plot, no baseline; the coalition-based version was superseded 2026-08-02) | `out_files/musica_native20_forward/{0000,1000,1100,1110,1111}/*.nc`; `tab_hobo_perplot_cluster.csv` | 3 base paths ×2: `outputs/figures_chap1/fig_trait_vertical_gradient_effects{_norm,_metres,}` `.{png,pdf}`; `tab_trait_vertical_gradient_subcanopy.csv` | `B3_trait_vertical_gradient.png` | MATCH `2cb3dd66` (vs `_norm`) | B10 |
| **Fig C1** residual vs topo | 423 | `~R/scripts/make_residual_vs_topo.R` | `tab_hobo_native20_validation.csv`; `in_files/lad_z05/Blois_lad_z05_r25.csv`; **`~/Documents/NC_Full/03_RESULTS/Blois/Metrics/Not_Masked/{dtm,slope,twi,northness,std,rumple}*.tif`**; `tab_hobo_perplot_cluster.csv` | `outputs/figures_pipeline/annex/fig_residual_vs_topo.{png,pdf}`; `tab_residual_vs_topo.csv`; `_rds_A8a_topo.rds` | `C1_residual_vs_topo.png` | MATCH `bed32283` | B11 |
| **Fig D1** HOBO coverage | 450 | `~R/scripts/c1_hobo_coverage_native.R` | `clhs_sample_native20_floor05.rds`; `in_files/data_Blois_utm31n.geojson`; `in_files_native20/{lai_z1,max,fCover}_res_10_m.tif`; `tab_hobo_perplot_cluster.csv` | `out_files/Chapter1/figures/FigSh_hobo_coverage_native.png` | `D1_hobo_coverage.png` | MATCH `576b4396` | B12 |
| **Fig E1** trait collinearity | 472 | `~R/scripts/make_trait_collinearity_native20.R` | `clhs_sample_native20_floor05.rds`; `R/cluster_relabel.R` | `out_files/Chapter1/figures/fig_trait_collinearity_native20.{png,pdf}` | `E1_trait_collinearity.png` | MATCH `2faa73a8` | B20 |
| **Fig E2** observational profile shape | 493 | `~R/scripts/c1_obs_profileshape_test.R` | `in_files/lad_z05/Blois_lad_z05_r25.csv`; `tab_hobo_perplot_cluster.csv` | `outputs/figures_pipeline/annex/fig_obs_profileshape.{png,pdf}`; `out_files/Chapter1/tables/tab_obs_profileshape.csv` | `E2_obs_profileshape.png` | MATCH `81b04f06` | B13 |
| **Fig F1** scan-angle control | 514 | `~R/scripts/fig_scanangle_control_csv.R` | `in_files/lad_z05/Blois_lad_z05_r25.csv`; `in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv` | `outputs/figures_pipeline/annex/fig_scanangle_control.{png,pdf}` | `F1_scanangle_control.png` | **STALE** — shipped `6082d379` is 2700×1500 (this script); the file now on disk is `55cc73e8` 2700×1650, written by `control_scanangle_vs_lidr.R` at 15:13, after the 14:51 sync | B15 |
| **Fig F2** LAI-correction comparison | 542 | `~R/scripts/fig_lai_correction_compare.R` | `tab_correction_levels.csv`; `tab_correction_levels_perplot.csv` | `out_files/Chapter1/figures/F2_lai_correction_compare.{png,pdf}` | `F2_lai_correction_compare.png` | MATCH `f3c54552` | B16 (data A6) |
| **Fig I1** forward inclusion | 668 | `~R/scripts/c1_fig5_j1_convB.R` (same script as Fig 6) | identical to Fig 6 | `out_files/Chapter1/figures/FigJ1_forward_inclusion.{png,pdf}` | `FigI1_forward_inclusion.png` | MATCH `05f00708` | B4 |
| **Fig S1** hourly scatter | 684 | `~R/c1_hourly_scatter.R` | `in_files/FR-Blo_2021_v2.nc`; `CFG$hobo_temp_csv` = `in_files/Blois_data_temperature.csv`; `out_files/musica_hobo_native20/1111/*.nc` (53); `pipeline/00_config.R` | `outputs/figures_chap1/Fig_hourly_temp_scatter.{png,pdf}` | `FigS_hourly_temp_scatter.png` | MATCH `224ab7bd` | B17 |
| **Fig S2** attribution, hot days | 691 | `~R/scripts/c1_fig3_units.R` with `PERIOD=hot METRIC=dt` | `sensitivity_perplot_units_v2.csv` | `out_files/Chapter1/figures/Fig3_attribution_units_hot.{png,pdf}` | `FigS_attribution_hot.png` | MATCH `0c8c45f9` | B18 (data A3) |
| **Fig S3** glass ceiling, hot | 696 | `~R/c1_hot_figures.R` | `out_files/Chapter1/hot_extract.rds`; `tab_hobo_perplot_cluster.csv` | `outputs/figures_chap1/FigS_glass_ceiling_hot.{png,pdf}` | `FigS_glass_ceiling_hot.png` | MATCH `6f9cd28a` | B19 (data A4) |
| **Fig S4** forward, hot | 698 | `~R/c1_hot_figures.R` | same as S3 | `outputs/figures_chap1/FigS_forward_hot.{png,pdf}` | `FigS_forward_hot.png` | MATCH `c145c1be` | B19 |
| **Fig S5** nested regression, hot | 703 | `~R/c1_hot_figures.R` | same as S3 **plus** `in_files/data_Blois_utm31n.geojson`; `in_files_native20/{lai_z1,fCover,max}_res_10_m.tif`; `in_files/lad_z05/Blois_lad_z05_r25.csv` | `outputs/figures_chap1/FigS_obs_nested_hot.{png,pdf}` | `FigS_obs_nested_hot.png` | MATCH `3fd762b9` | B19 |

**24 of 25 figures MATCH** their generator's current output. The one exception is Fig F1 (see §7 item 2).
Fig B2 is pixel-identical but not byte-identical, because its montage step depends on the local
ImageMagick build; that is cosmetic.

### 2.2 Tables

| table | ms line | generating script | reads | writes | how it reaches the manuscript | stage |
|---|---|---|---|---|---|---|
| **Table 1** | 151–158 | `~R/scripts/c1_make_table1.R` | `sensitivity_perplot_units_v2.csv`; `R/cluster_relabel.R`. **Not** the cLHS rds, despite its own header comment | `out_files/Chapter1/tables/tab_table1_full.md`, `tab_table1_short.md`, `tab_table1_values.csv` | **hand-pasted.** Verified 2026-07-31: `tab_table1_full.md` is byte-for-byte identical to manuscript lines 152–158 | A8 |
| **Table G1** ΔTmax | 566+ | `~R/c1_attribution_bootstrap_units.R` — **authoritative** | `sensitivity_perplot_units_v2.csv` | `out_files/Chapter1/tables/tab_attribution_units_bootstrap.csv` (28 rows, ΔTmax / full summer only; `set.seed(1)`, **B = 4000**) | hand-pasted | A10 |
| **Table G2** slope | 595+ | `~R/c1_units_bootstrap_v2.R` — **authoritative** | `sensitivity_perplot_units_v2.csv` | `out_files/Chapter1/tables/tab_attribution_units_bootstrap_v2.csv` (112 rows = 4 metric×period × 4 P × 7 effects; `set.seed(1)`, NB = 2000); Table G2 is its `metric=sl, period=all` block | hand-pasted | A7 |
| **Table H1** MuSICA config | 612+ | **none — hand-written from the namelists** | — | — | transcribed from `~R/in_files/Blois/in_files/musica_soil.nml` and `musica_veg1.nml` | — |
| **Table H2** soil/leaf parameters | 634+ | **none — hand-written from the namelists** | — | — | same | — |

**How Tables G1 and G2 were adjudicated.** Two bootstrap scripts exist and both are in the runner.
The manuscript's **Table G1 reproduces `c1_attribution_bootstrap_units.R` (A10) to the last digit**,
including the four CI bounds where the two files disagree — P1 fCover −10 `(+0.014, +0.033)` [A10] vs
`(+0.016, +0.034)` [A7]; P2 LAI +0.5 lo `−0.062` vs `−0.063`; P2 profile hi `−0.057` vs `−0.056`;
P4 fCover +10 hi `−0.240` vs `−0.239`. **Table G2 matches A7's `metric=sl, period=all` block**, and only
A7 computes the slope metric at all. So the two tables come from *different scripts with different
bootstrap replicate counts* (4000 vs 2000) — correct as shipped, but a genuine trap, since
`tab_attribution_units_bootstrap_v2.csv` also contains a ΔTmax block that is **not** Table G1.

**Tables H1 / H2 verification.** No script emits them; the only `writeLines(*.md)` calls in the repo are
in `c1_make_table1.R`. The values were checked directly against the namelists that `~R/musica.nml`
actually points at (`PARAM_FILES_PATH = "./in_files/Blois/in_files/"`):

* `musica_soil.nml` — `SOIL_DEPTH_MAX 1.6`, `ALBEDO_SOIL_VIS 0.15` / `_NIR 0.28`, `N_DEPTH 4`,
  `DEPTH 0.04, 0.14, 0.28, 0.80`, `THETA_SAT_SOIL 0.54, 0.50, 0.47, 0.47`, `THETA_RES_SOIL 0.02`,
  `KSAT_SOIL 0.106, 0.041, 0.013, 0.001` — all match H1/H2 verbatim.
* `musica_veg1.nml` — `BUDBURST_DATE 115`, `T_OPT_JMAX 38.`, `TETA_JMAX 0.7`, `LMA_CANOPY_TOP 0.1`,
  `GS_SLOPE 10.`, `GS_INTERCEPT 1.e-3`, `GS_HX_HALF −1.3`, `GS_HX_SHAPE 2.6`, `HYPOSTOMATOUS .true.`,
  `LEAF_INCLINATION_INDEX 0.63` — all match.

Note the repo carries **two divergent copies** of each namelist (root vs `in_files/Blois/in_files/`);
only the latter drove the runs. See §7 item 12.

### 2.3 Derived tables under `out_files/Chapter1/tables/`

The directory holds **171 entries**. Only the following are on the Chapter-1 chain.

| file | writer (stage) | readers |
|---|---|---|
| `tab_hobo_perplot_cluster.csv` | `c1_align_cluster_obs.R` (A2, **opt-in and guarded**) | **14 chain readers** — B1, B4, B5, B12, B13, B19, A5, A6, A11, A12, `make_vertical_profiles_native.R`, `make_trait_vertical_gradient.R`, `make_residual_vs_topo.R`. The most-read table in the chapter |
| `tab_hobo_native20_validation.csv` | `c1_hobo_native20_validation_regen.R` (A1) — **and `c1_hobo_native20_validation.R`, not in the runner** | B4, B11, A2, A9, A11 |
| `sensitivity_perplot_units_v2.csv` | `c1_sensitivity_units_extract2.R` (A3) | B2, B18, B3, A7, A8, A10 |
| `sensitivity_perplot_units/part_0{1..6}.csv` | legacy `c1_sensitivity_perplot_units.R` | read by A3 for its non-regression check only |
| `typical_plots_per_cluster.csv` | pre-existing, no chain writer | B1 |
| `tab_radius_sweep_aligned.csv` | `c1_radius_rescore.R` (A5) | B7 |
| `tab_correction_levels.csv` | `c1_correction_levels_rescore.R` (A6) | B16, Appendix F prose |
| `tab_correction_levels_perplot.csv` | A6 | B16 |
| `tab_correction_levels_split.csv` | A6 | prose only |
| `tab_table1_{full.md,short.md,values.csv}` | `scripts/c1_make_table1.R` (A8) | hand-pasted; `tab_table1_values.csv` also backs the Appendix A "13 to 33 m" claim |
| `tab_attribution_units_bootstrap.csv` | `c1_attribution_bootstrap_units.R` (A10) | hand-pasted → **Table G1** |
| `tab_attribution_units_bootstrap_v2.csv` | `c1_units_bootstrap_v2.R` (A7) | hand-pasted → **Table G2** |
| `tab_cover_floor_check.csv` | `c1_cover_floor_check.R` (A9) | Appendix C prose |
| `tab_canopy_ratio_check.csv` | `c1_canopy_ratio_check.R` (A11) | §4.3 prose |
| `tab_profile_centroid_stats.csv` | `c1_profile_centroid_stats.R` (A12) | §4.1 and Appendix B prose |
| `tab_scanangle_raw_stats.csv` | `c1_scanangle_raw_stats.R` (A13, optional, needs the external drive) | Appendix F prose |
| `tab_forward_inclusion_convB.csv` | `scripts/c1_fig5_j1_convB.R` (B4) | provenance record for Appendix I |
| `tab_pca_clusters.csv` | `c1_pca_recover_clusters.R` (B6) | Appendix E prose |
| `tab_obs_profileshape.csv` | `scripts/c1_obs_profileshape_test.R` (B13) | Appendix E prose |
| `out_files/Chapter1/hot_extract.rds` (sibling, not in `tables/`) | `c1_hot_extract.R` (A4) | B19 |

Three more chain tables land **outside** this directory:
`outputs/figures_pipeline/annex/tab_residual_vs_topo.csv` (B11),
`outputs/figures_pipeline/annex/tab_h2_controlled_topheavy.csv` (B21),
`outputs/figures_chap1/tab_subcanopy_native.csv` and `tab_trait_vertical_gradient_perplot.csv` (B9, B10).

**Superseded — present in `tables/` but NOT part of the chain.** Nothing in stages A/B reads or writes
these. Scoping them out explicitly so a reader does not mistake them for provenance:

* **`frblo`** — `tab_hobo_frblo_validation{,_ctrl,_sacorr}.csv`, `tab_deltadelta_FRBLOiter.csv`, dirs
  `metrics6_frblo/`, `metrics6_frblo_wc/`, `sensitivity_perplot_frblo{,_windcorr}/`, `sens_frblo_FIXEDSTEP_bak/`
* **`era5sd`** — `tab_deltadelta_ERA5SD.csv`, `tab_deltadelta_ERA5iter.csv`, dirs `metrics6_era5sd/`,
  `sensitivity_perplot_era5sd/`
* **`z05`** — residue `tab_hobo_direction_z05branch.csv` (written only by the forbidden `c1_hobo_direction.R`)
* **`v320`** — `tab_stats_v320_v323.csv`, `tab_stats_tests_v320_v323.csv`, dir `metrics6_v320/`
* **`shapley`** — dirs `shapley_parts{,_global,_blois_bak,_global_blois_bak}/`; files
  `Table_shapley_clhs_{cluster,global}.csv`, `Table_shapley_vs_loo_forward.csv`,
  `tab_shapley_clhs_{cluster_importance,global,global_signed}.csv`. Shapley was dropped from Chapter 1 on 2026-06-20
* **`convA`** — `tab_hobo_native20_validation_convA_backup.csv`, `_convAB.csv`, dirs
  `metrics6_native20_convA_backup/`, `metrics6_native20_convB/`
* **`iter`** — `tab_hobo_{iter,station_iter}_validation.csv`, `tab_forward_v323iter_station.csv`,
  `tab_version_iter_compare.csv`, dirs `metrics6_iter/`, `sensitivity_perplot_iter/`
* **`fpc`** — `tab_fpc_typology_{centroids,preview}.csv`, `tab_fpc_{preview,results}_{sensitivity,bootstrap}.csv`,
  `tab_fpc_archetype_sensitivity.csv`, dir `sensitivity_perplot_fpc/`
* **`metrics6_*`** — only `metrics6_native20/` is authoritative; no stage A/B script reads any of them
* **other stale residue** — `tab_hobo_safran_validation.csv`,
  `tab_hobo_station_validation_station{Tool,Rmerge}.csv`,
  `tab_hobo_perplot_cluster_{preAlign_backup,VCIbuggy_bak,WRONGLAI}.csv`,
  `tab_k_rescale_test_preWind.csv`, `tab_cluster_meansd.csv` (flagged stale by `c1_make_table1.R:10-11`),
  `tab_importance_sensitivity.csv` (**must not be cited**), ~40 `*.log` files, and 12 legacy
  `Table_*` / `radius_metrics_*` CSVs.

---

<!-- SECTION 3 -->
## 3. Appendix by appendix

This is the part `RUN_CHAPTER1.md` and the figure ledger do not cover: the quantitative claims that live
in appendix **prose and figure captions**, which are neither a figure file nor a table.

Every number below was checked by opening the named script and either (a) locating the value in a
persisted table, or (b) recomputing it read-only from the primary input. Status vocabulary:

* **TRACEABLE** — value found in a persisted file or exactly recomputed from one.
* **TRACEABLE (rounded)** — reproduces at the manuscript's stated precision (e.g. table 0.6866 → text 0.69).
* **TRACEABLE (not persisted)** — reproduces, but the generator only `cat()`s it; no file holds it.
* **MISMATCH** — the source disagrees with the text. These are the findings that matter.
* **NOT TRACEABLE** — no producer exists in the repo.

Line numbers refer to `chapter1/manuscript/manuscript_chap1_EN_native20.md`.

---

### Appendix A — The forcing (lines 309–347)

**Asserts:** air temperature is taken from the CHS 41 station rather than ERA5 because ERA5 carries a
diurnally structured warm bias; and the ERA5 10 m wind is mapped to canopy height by a neutral
log-profile correction that is applied in every simulation.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 314 | ERA5 warm by 0.96 °C | `c1_era5_station_bias.R` | not persisted — `cat` only; reads `in_files/musica_in_Blois.nc` + `in_files/FR-Blo_2021_v2.nc` | TRACEABLE (rounded) — recomputed +0.9561 |
| 314 | RMSE 2.06 °C | same | not persisted | TRACEABLE (rounded) — 2.0588 |
| 314 | *r* = 0.92 | same | not persisted | TRACEABLE (rounded) — 0.9187 |
| 315 | peak +3.6 °C at 06:00–09:00 UTC | same | not persisted | TRACEABLE — max +3.637 at 07 UTC |
| 316 | turns negative in the evening | same | not persisted | TRACEABLE — −0.929 (18 h), −0.918 (19 h) |
| 321–324 | Fig A1 caption repeats the above; "summer 2021, hourly" | same | script window 2021-06-01→09-30 | TRACEABLE — *n* = 2928 h = 122 d × 24 |
| 329–331 | *U*(*h*+2m)/*U*(10m) formula, *d* = 0.7*h*, *z*₀ = 0.1*h* | `c1_wind_profile_correction.R` L8–9 | in code | TRACEABLE |
| 332 | ERA5 roughness **0.44 m** | `c1_wind_profile_correction.R` L13 | **hardcoded literal** | **NOT TRACEABLE** |
| 333 | **0.438 to 0.443 m** June–September | none | **source-code comment only** | **NOT TRACEABLE** |
| 333 | grassland 0.01 m, forest 1 m | `c1_wind_profile_correction.R` L10–14 | in code | TRACEABLE (literature values) |
| 334–335 | factor **0.41 to 0.47** over 13–33 m | same, `ratio()` L9 | recomputed | **MISMATCH** — at 13–33 m the factor is 0.48–0.41. The printed 0.47–0.41 corresponds to **15–35 m**, the range the script itself prints |
| 335 | 13 to 33 m archetype heights | `scripts/c1_make_table1.R` | `out_files/Chapter1/tables/tab_table1_values.csv`, col `hm_m` (P1 13.316, P4 32.561) | TRACEABLE (rounded) |
| 335–336 | wind 2.1 to 2.4× too strong | `c1_wind_profile_correction.R` | recomputed 1/0.4843 = 2.065, 1/0.4106 = 2.435 | TRACEABLE |
| 342–344 | Fig A2 four roughness levels; "about 0.4", "2.3 times" | same | recomputed 2.308 at *h* = 23 m | TRACEABLE |

> **Finding A-i.** The ERA5 surface-roughness value 0.44 m — a stated model input — is a bare literal
> with no download script, no cached field and no NetCDF variable behind it. The monthly range
> 0.438–0.443 m exists only as a code comment.
>
> **Finding A-ii.** `c1_wind_profile_correction.R`'s own header contradicts both its code and the
> manuscript: it says the Blois value is "0.40 m", the range "0.30–0.58 m", and that the correction is
> "**NOT applied in the simulations**". The code uses 0.44 / 0.438–0.443 and the manuscript (line 336)
> states the correction *is* applied throughout. The header is stale — but it is the only written
> provenance the roughness value has.

---

### Appendix B — The vertical dimension (lines 348–394)

**Asserts:** a controlled balance test shows concentrating foliage high cools more and that this
potency saturates with density (H2 is real but latent); and the within-canopy gradients show the
real-versus-uniform profile effect is small at the one level that is validated.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 364–365 | ~~P1 largest sub-canopy Tair increment ≈ 0.35 °C~~ **RETIRED 2026-08-02**: baseline-relative value, replaced by the per-plot figures (P1 ΔLAI +0.036 °C) | `scripts/make_trait_vertical_gradient_perplot.R` | `outputs/figures_chap1/tab_trait_vertical_gradient_subcanopy.csv`, P1/Tair `dLAI` = −0.352 | TRACEABLE (rounded) |
| 365 | P2/P3 profile leads, 0.13 and 0.19 °C | same | same, P2/Tair `dLAD` −0.127; P3/Tair `dLAD` −0.194 | TRACEABLE (rounded) |
| 366 | largest RH increment in P3 | same | same, P3/RH largest = `dLAD` (−0.720) | TRACEABLE |
| 367 | P4 cover ≈ 0.57 °C | same | same, P4/Tair `dfCover` −0.566 | TRACEABLE (rounded) |
| 367–368 | P4 largest humidity and VPD increments | same | same, P4/RH `dfCover` 2.467; P4/VPD −0.108 | TRACEABLE |
| 372 | Fig B1 "**v3.2.3 iter, CHS 41 station forcing**" | `scripts/make_h2_controlled_topheavy.R` | uses `CFG$forcing_file` = `in_files/musica_in_Blois.nc` (**ERA5**) and `CFG$musica_cmd` = **legacy Nov-2024 binary** | **MISMATCH (provenance)** — see Finding B-i |
| 373–374 | two-sided LAI 6 / 12, one-sided ≈ 3 / 6 | same, `LAI_LEVELS <- c(6, 12)` | in code | TRACEABLE |
| 377 | ΔTmax moves up to **0.59 °C** at moderate density | same | `outputs/figures_pipeline/annex/tab_h2_controlled_topheavy.csv`, LAI = 6 span | TRACEABLE (rounded) — 0.5853 |
| 377–378 | only **0.27 °C** at high density | same | same, LAI = 12 span | TRACEABLE (rounded) — 0.2725 |
| 378 | centroid median 0.48 over 53 loggers | `c1_profile_centroid_stats.R` | `out_files/Chapter1/tables/tab_profile_centroid_stats.csv`, row `all loggers` | TRACEABLE — 0.48011, *n* = 53 |
| 378 | middle half 0.43 to 0.52 | same | same row, q25 0.42620, q75 0.52040 | TRACEABLE (rounded) |
| 378 | only P1 departs, at 0.24 | same | same file, row P1 median 0.23810 (*n* = 8) | TRACEABLE (rounded) |
| 381 | Fig B2 "10:00 to 16:00 **on the forcing clock (solar time, UTC+1)**" | `scripts/make_vertical_profiles_native.R` L43 | window is applied **on UTC**; the header calls this "a deliberate simplification" | **MISMATCH (clock convention)** |
| 381 | "brackets the **15:00** peak" | none | no script references a 15:00 or hour-15 diurnal peak | **NOT TRACEABLE** |
| 381, 393 | June to September | same, `MONTHS <- 6:9` | in code | TRACEABLE |
| 383 | 53 logger pixels grouped by archetype | same | `tab_hobo_perplot_cluster.csv` col `P` → 8+12+13+20 = 53 | TRACEABLE |
| 383–384 | real = coalition 1111, uniform = 1110 | same, `SHAPES` L45 | in code | TRACEABLE |
| 385–386 | Fig B2 is a two-panel montage | `run_chapter1.R` L174–182 | ImageMagick `-append` of the `_norm` and `_metres` panels | TRACEABLE |
| 388 | real-vs-uniform shifts Tair 0.05 to 0.20 °C | `scripts/make_vertical_profiles_native.R` | `outputs/figures_chap1/tab_subcanopy_native.csv`, col `Temp`, level 1 | TRACEABLE — P1 −0.07, P2 −0.13, P3 −0.20, P4 +0.05 |
| 393 | Fig B3 increments = per-plot native steps around each plot's own canopy (no baseline) | `scripts/make_trait_vertical_gradient_perplot.R` | in code | TRACEABLE |
| 393 | baseline one-sided LAI **3.42**, *H*max **23 m** | design sample | `clhs_sample_native20_floor05.rds` (*n* = 400) | TRACEABLE (rounded) — 3.4186 / 22.6625 |
| 393 | P4 → 5.47, P1 → 1.16, P2 3.73, P3 3.67 | cluster table | `tab_hobo_perplot_cluster.csv` col `LAI` by `P` | TRACEABLE (rounded) — 5.4695 / 1.1559 / 3.7251 / 3.6736 |
| 393 | "stripping two thirds" | same | 1 − 1.1559/3.4186 = 0.662 | TRACEABLE |

> **Finding B-i — highest-severity provenance error in the appendices.** Figure B1's caption states
> "v3.2.3 iter, CHS 41 station forcing". Its generator `scripts/make_h2_controlled_topheavy.R` calls
> `run_musica_one(..., CFG$forcing_file, CFG$musica_cmd)` with **no override**, so it ran on
> `in_files/musica_in_Blois.nc` (**ERA5**-based, not the canonical CHS 41 `FR-Blo_2021_v2.nc`) and on
> `/home/corroyez/Documents/musica/musica` (md5 `53072278`, **Nov 2024 legacy**), which `R/config.R`
> L11–16 explicitly documents as *a different model variant* — "substituting it shifts ΔTmax by up to
> 1.4 °C and inverts the buf/amp regime ratio". The v3.2.3 binary used for every other Chapter-1
> simulation is `in_files/model-3.2.3/musica` (md5 `8b5139e0`).
>
> The **values** 0.59 / 0.27 °C are exact against `tab_h2_controlled_topheavy.csv`. Only the
> **attribution** is wrong. And because the 14 NetCDFs are cached (cache guard at L83–87 re-extracts
> rather than re-runs), the caption **cannot become true by re-running the figure**. Either the caption
> is corrected to name the ERA5 forcing and the legacy binary, or the controlled test is re-run on the
> canonical pair. This is an appendix-level sensitivity check, not a headline result, so it does not
> invalidate the chapter — but as written the caption misstates a model configuration.

---

### Appendix C — Residual decomposition (lines 395–424)

**Asserts:** the validation residual is governed by horizontal canopy openness (fractional cover), not
by vertical heterogeneity or terrain; the failure is structural to a 1-D column, not an artifact of the
0.5 cover floor; and the residual carries no detectable spatial autocorrelation.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 398, 400 | all 53 loggers, *n* = 53 | `scripts/make_residual_vs_topo.R` | `outputs/figures_pipeline/annex/tab_residual_vs_topo.csv` (53 rows) | TRACEABLE |
| 400 | adjusted *R*² = 0.69 | same | not persisted — printed at L78–80 | TRACEABLE (rounded) — 0.6866 |
| 401 | fCover *r* = −0.84 | same | recomputed from the saved CSV | TRACEABLE (rounded) — −0.8356 |
| 401 | fCover partial *p* < 0.001 | same | recomputed | TRACEABLE — *t* = −10.28, *p* ≈ 1e−13 |
| 402 | elevation **and northness**, *p* ≈ 0.05 | same | recomputed | **MISMATCH** — elevation 0.0381 ✓, **northness 0.9314** ✗ (simple correlation *p* = 0.551) |
| 402 | slope and TWI not significant | same | recomputed 0.2864 / 0.9780 | TRACEABLE |
| 405 | cover 0.51 vs 0.86, *p* < 0.001 | same | recomputed | TRACEABLE (rounded) — 0.5061 (*n* = 8) vs 0.8642 (*n* = 45), *p* = 2.3e−31 |
| 406 | the 0.5 cover floor | same L60 | `D[, fCover := pmax(fCover, 0.5)]` | TRACEABLE |
| 406 | four of the **seven** floored loggers | `c1_cover_floor_check.R` | `out_files/Chapter1/tables/tab_cover_floor_check.csv` (4 rows) | TRACEABLE — 7 of 53 loggers < 0.5 confirmed independently |
| 406 | true cover **0.00 to 0.42** | same | same CSV, col `fcover_true` (0.000482 … 0.418890) | TRACEABLE (rounded) |
| 406 | mean **+0.22 → +0.21 °C** | same | same CSV, cols `sim_floored`, `sim_truecover` | TRACEABLE (rounded) — 0.22290 → 0.21391 |
| 406 | vs **+2.73 °C** observed | same | same CSV, col `obs` | TRACEABLE — 2.7320 |
| 409–410 | sd_height and rumple add nothing independent | `scripts/make_residual_vs_topo.R` | recomputed partials *p* 0.721 / 0.799 | TRACEABLE |
| 415 | observed Moran's *I* = 0.10, *p* = 0.04 | `scripts/make_moran_validation.R` | `outputs/figures_pipeline_z05/data/tab_moran_validation.csv` | TRACEABLE (rounded) — **but from the z05 lineage** |
| 415 | residual *I* = 0.05, ***p* = 0.13** | same | same file gives *I* 0.035, *p* 0.164; native20 recompute gives 0.049, *p* 0.098 | **MISMATCH** — no lineage produces *p* = 0.13 |
| 415 | "the model captures that structure" | same | modelled *I* = 0.226, *p* = 0.001 | TRACEABLE but **understated** — modelled *I* is ~2× observed |
| 418–419 | Fig C1 adj *R*² 0.69, *r* = −0.84 | same | built into the subtitle at L95–97 | TRACEABLE (rounded) |
| 421 | fCover **VIF 1.3** | **none** | no VIF computation anywhere in the repo | **NOT TRACEABLE** (value reproduces: 1.260) |
| 422 | sd_height / rumple **VIF ≈ 9.5** | **none** | same | **NOT TRACEABLE** (reproduces: 9.529 / 9.405) |

> **Finding C-i.** "Elevation and northness, *p* ≈ 0.05" is unsupported for northness under either
> reading (partial *p* = 0.931, simple *p* = 0.551). Elevation alone carries the topographic-exposure claim.
>
> **Finding C-ii.** The Moran pair mixes lineages. The observed pair (0.10, 0.04) matches the committed
> **z05** table; the residual *I* matches a native20 recompute but its *p* = 0.13 matches neither
> (0.098 native20, 0.164 z05). `moran.mc` is a permutation test, so *p* is seed-dependent — but under
> the script's own `set.seed(42)` neither lineage yields 0.13. The only committed producer reads
> `outputs/figures_pipeline_z05/data/ref_validation.rds` and has never been re-pointed at the chapter's
> canonical native20 validation. The qualitative conclusion (no detectable residual autocorrelation)
> holds under both lineages; the printed *p* does not.
>
> **Finding C-iii — RESOLVED 2026-07-31.** Both VIF values in the Fig C1 caption reproduced exactly
> from the saved design matrix but **no script in the repo computed a VIF**; they had been computed
> once in a session and typed in. `scripts/make_residual_vs_topo.R` now computes them inline, as
> 1 / (1 − *R*²_j) of each covariate on the others (no extra package), prints them and writes
> `tab_residual_vif.csv`. Re-running it returns fCover 1.26, sd_height 9.53, rumple 9.41, i.e. the
> caption's 1.3 and ≈ 9.5, alongside the adjusted *R*² = 0.69 the caption also quotes.

---

### Appendix D — Model-free observational triangulation (lines 425–451)

**Asserts:** the quantity-over-structure result exists in the observations and not only inside MuSICA —
leaf quantity alone matches the full model's fit, and adding height or the profile adds nothing.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 428 | 53 loggers, four LiDAR traits | `scripts/c1_fig6_obs_nested_native.R` | `out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv` (53 rows) | TRACEABLE |
| 431 | full model *R*² = 0.79 | same | not persisted — `cat` only | TRACEABLE (not persisted) — 0.7900 |
| 431 | adjusted *R*² = 0.77 | same | not persisted | TRACEABLE (rounded) — 0.7725 |
| 431 | LAI+fCover reaches the same 0.79 | same | not persisted | TRACEABLE — 0.7899 |
| 431 | structure-only with VCI *R*² = 0.67 | `c1_obs_nested_r2.R` | not persisted | TRACEABLE — 0.6707 |
| 431 | structure-only with profile *R*² = 0.44 | `scripts/c1_fig6_obs_nested_native.R` | not persisted | TRACEABLE (rounded) — 0.4354 |
| 431 | VCI version raises fit by **0.13** | both nested scripts | not persisted | TRACEABLE (rounded) — but see Finding D-i |
| 431 | profile version by **0.35** | same | not persisted | TRACEABLE — 0.7900 − 0.4354 = 0.3546 |
| 435 | LAI 0.86, fCover 0.86, VCI 0.79, *H*max 0.65 | `c1_obs_nested_r2.R` | `tab_hobo_perplot_cluster.csv` | TRACEABLE — 0.8636 / 0.8640 / 0.7896 / 0.6500 |
| 436 | design pairwise *r* = 0.30 to 0.94 | `scripts/make_trait_collinearity_native20.R` | `clhs_sample_native20_floor05.rds` | TRACEABLE — 0.3014–0.9354 |
| 436 | logger pairwise *r* = 0.60 to 0.89 | — | `tab_hobo_perplot_cluster.csv` | TRACEABLE — 0.6039–0.8926 |
| 437 | LAI–fCover *r* = 0.89 at loggers | — | same | TRACEABLE — 0.8893 |
| 439, 448 | one-sided LAI up to ≈ 6.7 | `scripts/c1_hobo_coverage_native.R` | `in_files_native20/lai_z1_res_10_m.tif` at loggers | TRACEABLE — max 6.713 |
| 439, 448 | P4 mean ≈ 5.0 | same | `clhs_sample_native20_floor05.rds` + `R/cluster_relabel.R` | TRACEABLE — 5.025 |
| 448 | **26 of 53** have one-sided LAI ≥ 4 | same | not persisted | TRACEABLE (not persisted) — 26/53 |

> **Finding D-i.** The "+0.13" at line 431 silently switches which full model is differenced: it is
> *full-with-VCI* (0.7958) − 0.6707 = 0.125. The manuscript only ever prints 0.79 (the *profile* full
> model), so a reader computing 0.79 − 0.67 gets **0.12**. The companion "+0.35" is unambiguous.
> Recommend rewording to "about 0.12 to 0.13", or printing 0.80 for the VCI-full model.

---

### Appendix E — Trait correlation structure (lines 452–494)

**Asserts:** the traits are strongly collinear across the design and within every archetype, which is
why a single regression cannot attribute the thermal effect; VCI is redundant with density; and the
profile *shape*, tested directly and model-free, adds no robust buffering signal.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 455, 472 | 400 cLHS plots | `scripts/make_trait_collinearity_native20.R` | `clhs_sample_native20_floor05.rds` | TRACEABLE — 400 rows, 0 incomplete |
| 457 | LAI–fCover *r* = **0.94** (design) | same | same rds | TRACEABLE (rounded) — 0.9354 |
| 457, 469 | LAI–VCI *r* = **0.73** | same | same rds | TRACEABLE — 0.7317 |
| 462, 472 | P1 **0.79**, P2 **0.97**, P3 **0.93**, P4 **0.95** | same + `R/cluster_relabel.R` | same rds | TRACEABLE — 0.7915 / 0.9727 / 0.9302 / 0.9521 |
| 475 | PC1 loadings −0.55 LAI / −0.56 fCover | `c1_pca_recover_clusters.R` | `out_files/Chapter1/tables/tab_pca_clusters.csv`, col `PC1_load` | TRACEABLE (rounded) — −0.5515 / −0.5593 |
| 475 | PC1 explains ≈ 75% | same | not persisted | TRACEABLE (rounded) — 74.53% |
| 475 | VCI PC1 −0.47, PC2 −0.55 | same | `tab_pca_clusters.csv` | TRACEABLE — −0.4729 / −0.5544 |
| 475 | VCI *r* = 0.73 and 0.74 | same | rds | TRACEABLE — 0.7317 / 0.7386 |
| 485, 492 | 53 loggers | `scripts/c1_obs_profileshape_test.R` | `out_files/Chapter1/tables/tab_obs_profileshape.csv` | TRACEABLE — 53 rows |
| 486 | top-heaviness *r* = −0.48, *p* < 0.001 | same | same, cols `dTmax_obs`, `topheavy` | TRACEABLE — −0.4840, *p* = 2.4e−4 |
| 487, 492 | partial *r* = 0.00, *p* = 0.998 | same | same (+ `LAI`, `fCover`) | TRACEABLE — 0.00026, *p* = 0.9985 |
| 487, 492 | dense third partial *r* = −0.41, *p* = 0.09, *n* = 18 | same | same, subset `LAI ≥ quantile(LAI, 2/3)` | TRACEABLE — −0.4107, *p* = 0.0904, *n* = 18 |

> **Finding E-i — a previously flagged blocker is resolved.** `RUN_CHAPTER1.md` §6 item 7 declares
> Appendix F "unreproducible — highest-severity item, blocks submission", on the grounds that no script
> emits `0.95 / 0.91 / 0.79`. **That note is keyed to the SHORT manuscript's numbering and is stale for
> this document.** In the full manuscript the correlation appendix is **E**, and its stated
> 0.94 / 0.73 recompute exactly from `clhs_sample_native20_floor05.rds` as 0.9354 / 0.7317, as do all
> four per-archetype values. The shipped figure is byte-identical
> (md5 `2faa73a8fd75afa21fe7b20d715c2504`) to `out_files/Chapter1/figures/fig_trait_collinearity_native20.png`,
> the output of `scripts/make_trait_collinearity_native20.R` — **not** the `floor05_v2`
> `FigAnnex_corrplot_traits.png` (md5 `f8eabf8b…`) the note points at. Text and figure agree
> structurally, from the same `cor()` call on the same sample. Re-verified 2026-07-31.
>
> `RUN_CHAPTER1.md` §6 item 8 (unseeded PCA on `floor05_v2`) is stale for the same reason:
> `c1_pca_recover_clusters.R` now carries `set.seed(42)` and reads the native20 rds.
>
> **Finding E-ii (wording).** Line 472's "Collinearity is strong everywhere" is contradicted by the
> figure's own data if read matrix-wide: per-archetype |*r*| bottoms out at 0.03–0.14 for LAI–*H*max,
> and LAI–VCI is +0.03 in P3 and **−0.19** in P4. The colon that follows scopes the claim to
> LAI–fCover, which is defensible; the blanket phrasing is not.

---

### Appendix F — LiDAR scan-angle correction (lines 495–543)

**Asserts:** a near-uniform ~5% scan-angle correction to leaf area leaves the headline scores and the
trait ranking intact, whereas the retrieval-optimized *k* = 0.65 (a further ~23% cut) collapses the
model's skill — so the leaf-area scale is not freely adjustable.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 501 | scan angles reach **±33°** | `c1_scanangle_raw_stats.R` | `out_files/Chapter1/tables/tab_scanangle_raw_stats.csv` | TRACEABLE (rounded) — max +33.0, min −32.0 (asymmetric) |
| 501 | ±30° nominal FOV, "just beyond" | same | same, `frac_beyond_30_pct` = 1.97 | TRACEABLE (±30 is a sensor spec) |
| 501 | median magnitude **16°** | same | same, `median_abs_deg` | TRACEABLE — 16.002 |
| 501 | 95th percentile **29°** | same | same, `p95_abs_deg` | TRACEABLE (rounded) — 28.998 |
| 501 | ⟨sec θ⟩ ≈ **1.05** | `scripts/fig_scanangle_control_csv.R` | `in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv`, col `sec_theta` | TRACEABLE — mean 1.0513 |
| 501, 504, 511, 528 | leaf-area reduction ≈ **5%** | same | `Blois_lad_z05_r25.csv` vs `..._sacorr_r25.csv` | TRACEABLE — 4.0366 → 3.8378 = −4.93% |
| 501 | nearly uniform, **sd 1.5%** | same | same | TRACEABLE — 1.54% |
| 512 | per-plot ⟨sec θ⟩ **1.02 to 1.09** | same | `Blois_lad_z05_sacorr_r25.csv` col `sec_theta` | TRACEABLE — [1.0188, 1.0905]; see Finding F-i |
| 503 | baseline *r* = **0.50** | `c1_correction_levels_rescore.R` | `out_files/Chapter1/tables/tab_correction_levels.csv`, row `full`, col `r` | TRACEABLE — 0.50070 |
| 503, 519 | baseline bias **+0.98 °C** | same | same, col `bias` | TRACEABLE — 0.98032 |
| 519 | amplitude **15%** | same | same, col `amp` | TRACEABLE (rounded) — 15.097 |
| 517–531 | *k* = 0.65 reduces LAI a further **23%** | — | arithmetic identity | TRACEABLE — 1 − 0.5/0.65 = 23.08% |
| 519 | *r* falls to **0.19** | `c1_correction_levels_rescore.R` | `tab_correction_levels.csv`, row `k065` | TRACEABLE — 0.19074 |
| 519 | amplitude to **3%** | same | same | TRACEABLE (rounded) — 2.629 |
| 519 | bias grows to **+1.33 °C** | same | same | TRACEABLE — 1.33249 |
| 537, 540 | 53 loggers, three correction levels | `scripts/fig_lai_correction_compare.R` | `tab_correction_levels.csv` (`n` = 53) | TRACEABLE — `stopifnot(nrow(M) == 3L)`; a 4th row `uncorr` is the wind axis, deliberately excluded |

> **Finding F-i (auditor trap).** The 1.02–1.09 range comes from
> `in_files/lad_z05/Blois_lad_z05_sacorr_r25.csv`, **not** from the obvious out_files artifact:
> `out_files/Chapter1/tables/tab_scanangle_control.csv` gives a per-id mean `sec` range of
> [1.0184, **1.0781**] — the 1.09 is not in it. Figure and caption do agree (the shipped F1 is the
> CSV-based rendering), but an auditor reaching for the out_files table will find a different number.
>
> **Finding F-ii.** The scan-angle statistics are computed over **60** footprints
> (`tab_scanangle_raw_stats.csv`, `n_footprints = 60`; the Fig F1 subtitle prints "n=60"), whereas
> every microclimate score is over **53** loggers. No caption asserts 53 for F1, so this is internally
> consistent — but it is a different plot set, and worth stating.
>
> **Finding F-iii — output-path collision.** Two scripts write
> `outputs/figures_pipeline/annex/fig_scanangle_control.png`:
> `scripts/fig_scanangle_control_csv.R` (9×5 in, CSV-based, in the runner as B15) and
> `scripts/control_scanangle_vs_lidr.R` (9×5.5 in, LAS-based, needs the external drive). The **shipped**
> figure is the CSV-based one (2700×1500, md5 `6082d379`); the file currently sitting in the annex
> directory is the LAS-based twin (2700×1650, md5 `55cc73e8`) that overwrote it. The numbers are
> unaffected — the caption reproduces from the CSVs — but the two renderings are not interchangeable,
> and a `--sync-only` run in the wrong order would ship the wrong one. See §7 item 1.

---

### Appendix G — Uncertainty of the trait effects (lines 544–607)

**Asserts:** bootstrap intervals confirm that leaf quantity and cover both act and strengthen toward the
dense end, that height never separates from zero, that the two step directions are not mirror images,
and that the small levers are physical rather than numerical.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 550 | **4000 bootstrap resamples** | `c1_attribution_bootstrap_units.R:22` (`bmed <- function(x, B = 4000)`) | in code | TRACEABLE **for Table G1 only** — G2 used `NB = 2000` |
| 595–596 | G2 caption: "on the same design and the same plots as Table G1" | `c1_units_bootstrap_v2.R:12` (`NB <- 2000`) | in code | **MISMATCH** — design and plots *are* identical, but G2 is a **separate RNG draw with half the resamples**. Combined with line 550's "4000", the caption implies a shared bootstrap that does not exist |
| 552–553 | P4 effects are "an order of magnitude above their effect in the open P1" | as Table G1 | same | **MISMATCH (partial)** — cover is 0.266/0.021 = **12.7×** ✓, but leaf area is 0.264/0.049 = **5.4×**, and P1's leaf-area effect is **+0.049 (warming)**, i.e. a sign reversal rather than a weaker cooling |
| 553–555 | height "never separates from zero by more than four thousandths of a degree per meter" | as Table G1 | same | **TRACEABLE on the medians** (max \|median\| = 0.004, P3) but **fails on the CI reading** the wording invites: P3 *H*max +1 m has `lo` = −0.006 |
| 566–567 | G1 caption: "steps are truncated at the trait bounds" | as Table G1 | `tab_attribution_units_bootstrap_v2.csv`, col `n` | **INCOMPLETE** — the caption describes truncation but not that two rows drop plots entirely: P1 fCover −10 uses **47** of 100 plots, P2 fCover −10 uses **97**. All other rows use 100 |
| **566–578** | **Table G1**, all 28 cells (median + 95% CI, 7 effects × 4 archetypes) | `c1_attribution_bootstrap_units.R` (A10) | `out_files/Chapter1/tables/tab_attribution_units_bootstrap.csv` | TRACEABLE — reproduces **to the last digit**, including the four bounds where the A7 file disagrees |
| **595–606** | **Table G2**, all 28 cells on the buffering slope | `c1_units_bootstrap_v2.R` (A7) | `out_files/Chapter1/tables/tab_attribution_units_bootstrap_v2.csv`, block `metric = sl, period = all` | TRACEABLE — only A7 computes the slope metric |
| 555–556 | −0.264 °C per half unit of leaf area, −0.266 °C per ten points of cover in P4 | as Table G1 | same | TRACEABLE |
| 557–558 | height never exceeds four thousandths of a degree per metre | as Table G1 | same | TRACEABLE — max \|median\| = 0.004 (P3) |
| 559–561 | P4 asymmetry: +0.264 cooling vs 0.224 warming | as Table G1 | same | TRACEABLE |
| 563 | leaf-quantity and profile intervals disjoint in P1, P3, P4; overlap in P2 | as Tables G1 | same | TRACEABLE by inspection of the CIs |
| 564 | profile leads in P2/P3: −0.086 and −0.136 vs −0.048 | as Table G1 | same | TRACEABLE |
| 565–566 | paired fraction where the leaf-area step exceeds the full profile swap: **63% P1, 82% P4, 42% P2, 21% P3** | **no repo script emits it** | recomputed from `sensitivity_perplot_units_v2.csv`: `mean(abs(LAI_up_dt_all*0.5) > abs(LAD_dt_all))` by archetype | **TRACEABLE but unreproducible by the pipeline** — 63 / 42 / 21 / 82 exactly. `c1_attribution_bootstrap.R:28` computes a *similar-looking* per-SD statistic that reproduces none of them (native20 → 75/45/29/91) |
| 588–592 | noise floor: ±0.01 step moves ΔTmax by a median **0.003 °C over 40 plots** | `c1_noise_floor.R` | `out_files/Chapter1/tables/tab_noise_floor.csv` (40 rows) | TRACEABLE — median \|noise\| = 0.0030 |
| 590–591 | warming **0.002 °C in P1**, cooling **0.008 °C in P4** | same | same, grouped by `relabel_cluster(Cluster)` | TRACEABLE — P1 +0.0019, P4 −0.0082 |
| 591 | reproduces the linear extrapolation of the ±0.5 lever (median predicted 0.003 °C) | same | same | TRACEABLE (not persisted) |

> **Finding G-i — Table G1 and Table G2 come from different scripts with different bootstrap counts.**
> G1 = `c1_attribution_bootstrap_units.R` (**B = 4000**, stage A10); G2 = `c1_units_bootstrap_v2.R`
> (**NB = 2000**, stage A7). Both are in the runner, and A7's file *also* contains a ΔTmax block that
> looks like Table G1 but is not — it differs in four CI bounds (P1 fCover −10, P2 LAI +0.5 low,
> P2 profile high, P4 fCover +10 high). Regenerating Table G1 from the "v2" file would produce subtly
> wrong numbers. Neither table is machine-synced into the manuscript.
>
> **Finding G-ii — `c1_noise_floor.R` runs MuSICA.** It is **not** in the runner (correctly), and its
> 80 NetCDFs are cached under `out_files/Chapter1/nc_noise_floor/` (220 MB) behind a
> `file.exists || file.size < 1000` guard. Deleting them and re-running would launch simulations. Added
> to the never-delete list in §4.

---

### Appendix H — MuSICA model configuration (lines 608–651)

**Asserts:** the model version, coupling mode, vertical grid, forcing composition, soil hydraulics,
phenology and leaf parameters used for every run.

**No script emits Tables H1 or H2.** They are transcribed by hand from the run namelists. The only
`writeLines(*.md)` calls in the repo are in `scripts/c1_make_table1.R`.

**The authoritative namelist chain**, confirmed two independent ways:

```
<repo>/musica.nml                                    ← the template rmusica copies
  ├─ SOIL_FILE_NAME   → in_files/Blois/in_files/musica_soil.nml
  └─ PARAM_FILES_PATH → in_files/Blois/in_files/musica_veg1.nml
```

(a) by reading `rmusica:::setup_musica` — `callmusica(musica.file = "musica.nml")` copies `./musica.nml`
from the cwd (the repo root) into a fresh `.musica_<hash>/` tempdir and resolves the soil/veg paths
*before* `setwd`; and (b) by diffing those three templates against the **surviving working directory**
`.musica_15682b4e92d382/` (2026-07-19 14:37, one day before the matching output NetCDFs) — byte-identical
apart from the launcher's deliberate overrides. Its `phenology_musica_veg1.csv` carries the duplicated
`2020 / DOY 366` row that is the iter-family signature.

Stronger than any comment: that working directory's `out.nc` has dimensions **`nair = 15`,
`nveg = 10`, `nsoil = 15`, `nleafage = 1`**, with `z_soil` bottoming near 1.6 m.

| claim | source | status |
|---|---|---|
| MuSICA v3.2.3, official binary | `in_files/model-3.2.3/musica` (md5 `8b5139e0`), pinned at `c1_sensitivity_perplot_chunk_native20.R:19`; md5-identical to the binary staged in the tempdir | TRACEABLE |
| iterative SBL mode, `ABL_flag = 'iter'` | injected at runtime by `c1_sensitivity_perplot_chunk_native20.R:21,70` as `extra_setup`. **The root nml has it commented out** (`! ABL_flag='up'`, L46) | TRACEABLE — but it comes from the **R launcher**, not the namelist |
| hourly, 3600 s | `musica.nml:69` `FORCING_TIMESTEP = 3600` | TRACEABLE |
| 15 air layers; 10 vegetation layers | `musica.nml:80` `AIR_RESOLUTION_LEVEL = 2`, documented at L78 as "15/10 layers in air/vegetation"; confirmed by `out.nc` dims | TRACEABLE |
| one leaf-age class | `musica.nml:97` `N_LEAF_AGE = 1`; `out.nc` `nleafage = 1` | TRACEABLE |
| **"Simulation period 1 June to 30 September 2021"** | `out.nc` `time` = **8759 hourly steps — all of 2021**. `R/config.R:18`'s `date_seq` is used only in `metrics_one()` | **MISMATCH** — that is a *post-processing extraction window*, not the model's integration span |
| **"Initialized on 1 June (no multi-year spin-up)"** | `musica.nml:44` `restin_filename = ""`; the run starts **1 January 2021** | **MISMATCH** — there is ~5 months of within-year spin-up before the metrics window |
| soil to **1.6 m**, 15 layers, 4 depth nodes at 0.04 / 0.14 / 0.28 / 0.8 m | `musica_soil.nml:75,92`; `musica.nml:93` `N_SOIL_LAYER = 15` | TRACEABLE |
| *θ*sat 0.54 → 0.47; residual 0.02 | `musica_soil.nml:129-130`: `THETA_SAT_SOIL 0.54,0.50,0.47,0.47`; `THETA_RES_SOIL 0.02` | TRACEABLE (rounded — the 0.50 intermediate node is omitted from the table) |
| soil albedo 0.15 vis / 0.28 NIR; emissivity 0.97 | `musica_soil.nml:80-82` | TRACEABLE |
| init soil moisture 0.40; temperature 283 → 278 K | `musica_soil.nml:204-205` | TRACEABLE |
| van Genuchten retention | `musica_soil.nml:112,128` `RETENTION_CURVE_MODEL_FLAG = 3` | TRACEABLE |
| *h*s 1.48 / 2.17 / 2.81 / 2.81 | `musica_soil.nml:131` `H_S_SOIL` | TRACEABLE |
| retention shape *n* = 1.05 | `musica_soil.nml:132` `N_SOIL` | TRACEABLE |
| **retention shape *m* = 0.5** | `BIG_M_SOIL = 0.5` is the **unsaturated hydraulic-conductivity exponent** (`HYDRAULIC_COND_MODEL_FLAG = 2`, L150-161). With `RETENTION_CURVE_MODEL_FLAG = 3` the namelist's own documentation fixes retention ***m* = 1** | **MISMATCH** — the value is real but **mislabelled** |
| *K*sat 0.106 / 0.041 / 0.013 / 0.001 | `musica_soil.nml:161` | TRACEABLE |
| budburst DOY **115** | `musica_veg1.nml:95` `BUDBURST_DATE 115`, re-passed by the launcher | TRACEABLE |
| leaf characteristic dimension 0.05 m | `musica_veg1.nml:196` `LEAF_SIZE` | TRACEABLE |
| LMA 0.1 kg m⁻² | `musica_veg1.nml:255` `LMA_CANOPY_TOP` | TRACEABLE |
| leaf inclination index 0.63 | `musica_veg1.nml:341-342` | TRACEABLE |
| *g*₁ = 10, *g*₀ = 0.001, hypostomatous | `musica_veg1.nml:283-289` | TRACEABLE — but `GS_MODEL_FLAG = 2` is **Ball-Berry**, so "*g*₁" is the Ball-Berry slope *m*, **not** a Medlyn USO *g*₁. Worth renaming |
| water-potential half-response −1.3 MPa, shape 2.6 | `musica_veg1.nml:290-291` | TRACEABLE — but with `PHOTOSYNTHETIC_LIMITATION_FLAG = 3` the namelist documents this as **predawn** potential updated daily, not the instantaneous "leaf xylem potential" the table states |
| *J*max optimum 38 °C, *θ* = 0.7 | `musica_veg1.nml:217,241` | TRACEABLE |
| single broadleaf deciduous, sessile oak | `musica.nml:132-133` `N_SPECIES = 1`; `musica_veg1.nml:44` deciduous. Species *identity* is comment-level only | TRACEABLE (species name documentary) |
| LAI entered two-sided (one-sided doubled) | `R/musica.R:69` `lai <- 2 * lai` | TRACEABLE |
| forcing level at *H*max + 2 m | `R/musica.R:93-96`; the tempdir shows `forcing_height = 35.99` for its plot | TRACEABLE |
| extinction coefficient *k* = 0.5 | `c1_fig1_gril_native.R:84` `lidR::LAD(..., k = 0.5)` — an **R data-prep choice, not a namelist entry** | TRACEABLE |
| forcing file FR-Blo_2021_v2.nc | launcher L20. But `musica.nml:49` still names a **SAFRAN** file (overridden programmatically) | TRACEABLE (with the §7 item 9 caveat) |
| **"Boundary-layer height: half-hourly MERRA-2 PBLH"** | — | **MISMATCH on three counts — see Finding H-i** |
| **BLH insensitivity, 100 m to 5000 m** | no script anywhere performs this sweep; grepping `5000` yields only coincidental matches plus one self-referential comment | **NOT TRACEABLE** |
| **"the complete namelist files are archived with the code"** (line 610) | `git ls-files --error-unmatch` fails for `musica.nml`, `musica_soil.nml` and `musica_veg1.nml` — **all three are untracked** | **MISMATCH** |

> **Finding H-i — the MERRA-2 boundary-layer-height claim is not supported by the forcing that was
> actually used.** Table H1 states the BLH is a "half-hourly MERRA-2 planetary-boundary-layer height".
> Three independent checks contradict this:
>
> 1. **Not half-hourly.** `get_merra2_pblh.sh` downloads `M2T1NXFLX.5.12.4` with `PBLH[0:1:23]` —
>    24 values per day, i.e. **hourly**. The forcing's own time axis is uniformly Δ*t* = 1.0 h.
>    Line 618 ("half-hourly") also contradicts line 619 ("FR-Blo_2021_v2.nc, hourly") *within Table H1*.
> 2. **Not the values in the forcing.** `h_sbl` carries `units = "m"` and **no `long_name` or any
>    MERRA-2 attribute at all**. It runs 1.20–244.25 m over the year (mean 60.8, sd 47.3; JJAS mean
>    64.5), with 8758 distinct values of 8759 — a real diurnal shape peaking mid-afternoon at
>    ~150–200 m, so not a flat placeholder. But the repo's own downloaded MERRA-2
>    (`out_files/Chapter1/merra2_pblh.csv`) gives **1663–2070 m** for 2021-07-15 13:00–16:00 against
>    this file's seasonal ceiling of ~200 m — roughly **10× too small**. Datetime-aligned the two
>    correlate only *r* = 0.53, and a ±3 h lag sweep peaks at 0.537; a pure rescaling would give 1.000.
> 3. **No script wrote it.** `build_forcing_pblh.R` targets `in_files/musica_in_Blois_pblh.nc` and
>    `c1_inject_hsbl_station.R` targets `out_files/MuSICA_in_CHS41-Blois_2021-station-hsbl.nc`. Both
>    build *genuine* MERRA-2-bearing forcings — and **neither output is the file the native20 launcher
>    points at**. The real MERRA-2 pipeline exists; its output never reached the runs.
>
> **Where the error came from.** `c1_sensitivity_perplot_chunk_native20.R:5-7`'s docstring claims "REAL
> MERRA-2 PBLH as h_sbl (forcing `musica_in_Blois_pblh.nc`)" — copy-pasted from the iter-family
> template — while executable line 20 sets `FORC <- "in_files/FR-Blo_2021_v2.nc"` with its own inline
> caveat: *"h_sbl is low (~100m, not real MERRA-2 PBLH) but VERIFIED immaterial to dTmax (identical
> 100m vs 5000m, 2026-07-23)"*.
>
> **The same claim appears at three manuscript locations — a fix must touch all three:**
> line **618** ("Half-hourly MERRA-2 planetary-boundary-layer height"), line **620** (the Forcing-sources
> row, "boundary-layer height from MERRA-2"), and lines **711–712** in Data availability — the worst of
> the three, because it asserts *archival provenance* for data that never entered the runs.
>
> **This does not invalidate the results.** The BLH-insensitivity test means the conclusions do not
> depend on the series. But that test is itself untraceable (no script performs the 100 m → 5000 m
> sweep; the launcher carrying the claim is untracked in git, so there is no commit history
> corroborating the date either). So the defence of the error is currently an assertion. The Methods
> description is wrong and should be corrected before submission, because it describes a model input.
>
> **Finding H-ii — no namelist is archived with the runs.** `find out_files -name '*.nml'` is empty,
> so Table H1's provenance is inferred from the *current* repo state rather than proven from the run
> artifacts. The 35 leftover `.musica_<hash>/` temp directories contain callmusica's rewritten
> lowercase copies (confirming `abl_flag="iter"` and `forcing_height=35.99` for one plot) but are not a
> reliable archive. See §7 item 10.

---

### Appendix I — Forward inclusion (lines 652–669)

**Asserts:** as LiDAR traits are added cumulatively to a uniform baseline canopy, leaf quantity and
cover carry the recovery of ΔTmax while the vertical profile adds the smallest increment, and the
buffering slope does not improve at all.

Every number reconciles against `out_files/Chapter1/tables/tab_forward_inclusion_convB.csv`, written by
`scripts/c1_fig5_j1_convB.R` (stage B4) — which re-derives the coalitions directly from
`out_files/musica_native20_forward/`.

| line | claim | value in the table (`grp = All`) | status |
|---|---|---|---|
| 660 | uniform baseline ΔTmax *r* ≈ **−0.72** | −0.7187 | TRACEABLE |
| 661 | uniform baseline slope ≈ **0.77** | 0.7655 | TRACEABLE |
| 662 | ΔTmax *r* to **0.10** on leaf area | +LAI 0.0980 | TRACEABLE |
| 662 | **0.19** with height | +H[max] 0.1930 | TRACEABLE |
| 662 | **0.43** with cover | +fCover 0.4326 | TRACEABLE |
| 663 | **0.46** with the real profile | +LAD 0.4582 | TRACEABLE |
| 663–664 | RMSE **2.22 → 1.75** | 2.2197 → 1.7542 | TRACEABLE |
| 664 | profile lowers it only to **1.73** | +LAD 1.7285 | TRACEABLE |
| 666 | slope drifts **0.77 → 0.71** | 0.7655 → 0.7133 | TRACEABLE |
| 667–668 | per-archetype *n* = **8 to 20** loggers | `tab_hobo_perplot_cluster.csv`: P1 8, P2 12, P3 13, P4 20 | TRACEABLE |

**Nothing in Appendix I is untraceable.**

---

### Supplementary figures (lines 670–704)

**Asserts:** repeating every analysis over the hottest 10% of days leaves the conclusions unchanged,
with every lever larger and more dispersed.

| line | claim | producing script | where the number lives | status |
|---|---|---|---|---|
| 673 | hottest 10% = **13 days**, ranked by station daily maximum | `c1_hot_extract.R` (A4) | `out_files/Chapter1/hot_extract.rds` (list: `obs`, `coal`) | TRACEABLE |
| 677–678 | per-plot dispersion of the leaf-quantity effect widens by **1.4 to 1.8×** | no script persists it | recomputable from `sensitivity_perplot_units_v2.csv` (`*_hot` vs `*_all`) | TRACEABLE (not persisted) |
| 681–682 | S1: *R*² = **0.91** overall, **0.89** on hot days | `c1_hourly_scatter.R` (B17) | not persisted — printed to stdout | TRACEABLE (not persisted) |
| 687–689 | S2: P4 leaf area **−0.36**, cover **−0.35** | `scripts/c1_fig3_units.R` `PERIOD=hot` | `sensitivity_perplot_units_v2.csv`, `*_dt_hot` columns | TRACEABLE — recomputed −0.356 / −0.346 |
| 688–689 | S2: profile **−0.15** (P2) and **−0.27** (P3) hot, vs **−0.09** / **−0.14** full summer | same | same, `LAD_dt_hot` / `LAD_dt_all` | TRACEABLE — −0.154 / −0.270 hot; −0.086 / −0.136 full |
| 694–695 | S3: ΔTmax *r* = **0.35**, slope *r* = **0.63** | `c1_hot_figures.R` (B19) | `hot_extract.rds` | TRACEABLE (not persisted) |
| 701–702 | S5: leaf area alone *R*² = **0.68** of the full model's **0.70** | `c1_hot_figures.R` | not persisted | TRACEABLE (not persisted) |

---

### Data and code availability (lines 705–717)

This section makes claims about **where the data live**, which is exactly this document's subject. It
contains **five unresolved placeholders**, all of which must be filled before submission:

| line | placeholder |
|---|---|
| 708 | `(REPOSITORY AND DOI TO BE SUPPLIED)` — the archive for the pipeline, namelists and derived products |
| 709–710 | `(DISTRIBUTION SOURCE TO BE SUPPLIED)` — where MuSICA v3.2.3 is obtained under the EUPL |
| 714 | `(PROVENANCE AND CONTACT TO BE CONFIRMED)` — the ALS cloud, CHS 41 record and logger series |
| 715–717 | `(R AND PACKAGE VERSIONS TO BE SUPPLIED)` — no `renv.lock`, `sessionInfo()` dump or version manifest exists anywhere in the repo |

> **Finding Avail-i — the opening sentence overclaims, in three specific ways.** Line 706 states that
> `run_chapter1.R` *"regenerates every figure, table and number reported here from the raw inputs."*
>
> 1. **Not from raw inputs.** It regenerates from **existing MuSICA NetCDF outputs** — the 12.6 GB of
>    §4 assets. The runner's own header says so ("from the MuSICA outputs that ALREADY EXIST on disk…
>    NO STAGE HERE RUNS MuSICA"). Going from raw inputs would require re-running MuSICA, which nothing
>    in the repo can do reproducibly (§7 items 10, 11).
> 2. **Not every figure.** Fig B1 is served from a **cache**; its generator invokes MuSICA and, as
>    §3 Finding B-i shows, was run on a different forcing and a different binary. Fig 1b depends on a
>    cached ALS cross-section, and Fig F1 is currently stale.
> 3. **Not every table.** Table 1, G1 and G2 are hand-pasted. H1 and H2 remain hand-written, but since
>    2026-07-31 they have a **checker**: `scripts/c1_dump_musica_config.R` prints every value they cite,
>    read from the namelists MuSICA actually stages (`in_files/musica_soil.nml`, `in_files/musica_veg1.nml`),
>    from the model source where a setting is derived rather than set (`air_resolution_level = 2` fixes the
>    10 vegetation / 15 air layers), and from the canonical forcing and one simulation output. Two traps it
>    now guards: `in_files/model-3.2.3/*.nml` are **distribution templates** whose setup block is overridden
>    at run time (`abl_flag`, forcing, forcing height, clumping factor), so the script carries a separate
>    AS-RUN column naming each override's source line; and the several `musica_soil.nml` copies are asserted
>    to agree on every key H1/H2 cites before anything is printed.
>
> Suggested wording: *"a staged pipeline (`run_chapter1.R`) that regenerates the analysis, figures and
> derived tables from the archived MuSICA simulation outputs; the simulations themselves are archived
> rather than re-run."*
>
> Note also that the ERA5 and MERRA-2 citations at lines 710–712 are correct as *sources*, but the
> MERRA-2 series did not in fact drive the reported runs (Finding H-i), and the canonical forcing's
> own composition is not verifiable from this repository (§1 Finding 1-0).

---

<!-- SECTION 4 -->
## 4. Simulation assets that must never be deleted

None of these can be rebuilt without thousands of CPU-hours of MuSICA. Sizes and counts measured
2026-07-31 on this machine.

| # | path | content | entries | `.nc` | size | artifacts that die if lost |
|---|---|---|---|---|---|---|
| 1 | `out_files/Chapter1/nc_sensitivity_perplot_units/` | trait-perturbation design: 8 perturbation tags × 400 cLHS plots | 3200 | 3200 | **8.6 G** | Fig 4, Fig 5, Fig S2, Table 1 (step/SD + profile-swap columns), Tables G1 and G2, Appendix G in full, §3.2 and §4.1 in full |
| 2 | `out_files/musica_native20_forward/` | 16 forward coalition sets `0000`…`1111`, 60 logger pixels each (53 retained after `CFG$ids_to_remove`) | 16 dirs | 960 | 2.6 G | Fig I1, Fig B2, Fig B3, Figs S3–S5, `out_files/Chapter1/hot_extract.rds`, Appendix I and Appendix B in full, §3.3 forward-inclusion result |
| 3 | `out_files/musica_hobo_native20/1111` | full-model validation run, real canopy, 53 loggers | 1 dir | 53 | 146 M | Fig 6, Fig 8, Fig C1, Fig D1, Fig E2, Fig S1, `tab_hobo_native20_validation.csv`, the headline validation (*r* = 0.50, bias +0.98 °C, slope *r* = 0.78) |
| 4 | `out_files/radius_test_corr/{5m,10m,12.5m,15m,20m,25m,50m}` | footprint-radius sweep on the corrected baseline, 60 pixels per radius | 7 dirs | 420 | 1.2 G | Fig 7, `tab_radius_sweep_aligned.csv`, §3.4 in full |
| 5 | `out_files/H2_controlled_topheavy/{LAI6,LAI12}` | controlled H2 balance test at two leaf-area levels | 2 dirs | 14 | 58 M | Fig B1 and the 0.59 °C / 0.27 °C claim in Appendix B. **Its generator `scripts/make_h2_controlled_topheavy.R` calls MuSICA** — the figure is served from a cache, so losing these NetCDFs makes Fig B1 unreproducible without a re-run |
| 6 | `in_files/FR-Blo_2021_v2.nc` | canonical station-based forcing (458 kB) | file | 1 | 448 K | **everything.** Every simulation above was driven by it; every re-extraction reads it for the macro reference |
| 7 | `out_files/Sensitivity_Analysis/clhs_sample_native20_floor05.rds` | the 400-plot cLHS design (400 × 50: `x, y, LAI, VCI, Hmax, fCover, LAD_Layer_1.5…39.5, FPC1-3, Cluster, sec_theta`) | file | — | — | Fig 1c, Fig 3, Fig 5, Fig D1, Fig E1, Table 1, Appendix E in full. Not a MuSICA product, but rebuilding it needs the ALS trait rasters *and* reproducing the cLHS draw |
| 8 | `outputs/figures_pipeline/annex/fig_h2_controlled_topheavy.png` | **cached render** of Fig B1 | file | — | — | Fig B1. The only figure in the chapter with no MuSICA-free generator |
| 9 | `in_files/fig1_crosssection_cache.rds` | cached ALS cross-sections for Fig 1b (517 kB) | file | — | — | Fig 1b. Rebuilding it requires the raw point cloud on the external drive |
| 10 | `out_files/Chapter1/nc_noise_floor/` | numerical noise-floor control, ±0.01 leaf-area step over 40 plots | 80 files | 80 | 220 M | Appendix G's noise-floor paragraph (median 0.003 °C; P1 +0.002, P4 −0.008) and `tab_noise_floor.csv`. **`c1_noise_floor.R` invokes MuSICA** behind a cache guard, so losing these triggers a simulation run. Not in the runner, and **not in the runner's pre-flight** |

`out_files/Sensitivity_Analysis/` as a whole is 7.8 G / 1928 `.nc`; most of it is superseded
lineage. Only item 7 above is on the canonical chain — but do not blanket-delete the directory
without checking §7.

**Pre-flight check** (this is exactly what `run_chapter1.R` asserts before it will run):

```bash
cd /home/corroyez/Documents/z_Example_rmusica_31012025
ls out_files/Chapter1/nc_sensitivity_perplot_units | wc -l   # 3200
ls out_files/musica_native20_forward | wc -l                 # 16
ls out_files/radius_test_corr | wc -l                        # 7
ls -1 in_files/FR-Blo_2021_v2.nc out_files/musica_hobo_native20/1111 >/dev/null
Rscript run_chapter1.R          # dry run: prints the pre-flight table and the 34-stage plan
```

---

<!-- SECTION 5 -->
## 5. Start here — clean checkout to compiled PDF

### What must be on disk first

The repository alone is **not** sufficient. You additionally need:

1. **The simulation assets of §4** (≈ 12.6 GB). They are not in git and cannot be regenerated
   without MuSICA. Without them nothing below works.
2. **The external drive `/media/corroyez/MyPassport`**, mounted, *only* if you need to rebuild the
   ALS-derived caches (Fig 1b cross-sections, the raw scan-angle statistics of Appendix F). The
   normal path uses `in_files/fig1_crosssection_cache.rds` and does not touch the drive.
3. **Toolchain**: `Rscript`, `pandoc`, `xelatex`, ImageMagick `convert` (for the Fig B2 montage),
   graphviz `dot` (for the Fig 2 schematic).
4. **`.cdsapirc`** in the repo root — a Copernicus CDS credential. Required only to *re-download*
   ERA5; not needed to reproduce the chapter from the existing forcing file.

### The commands

```bash
cd /home/corroyez/Documents/z_Example_rmusica_31012025

# 1. dry run — pre-flight on the simulation assets + the 34-stage plan. Executes nothing.
Rscript run_chapter1.R

# 2. regenerate the analysis and every figure from the EXISTING NetCDFs (~5 min).
#    No stage runs MuSICA.
CH1_RUN=TRUE Rscript run_chapter1.R

# 3. REQUIRED BEFORE SYNCING — refresh Fig F1, or the sync ships the wrong figure (see warning below).
Rscript scripts/fig_scanangle_control_csv.R

# 4. copy/rename the figures into the manuscript directory, then check the ledger.
#    Rehearse into a scratch dir first: CH1_FIGDIR=/tmp/figtest Rscript run_chapter1.R --sync-only
Rscript run_chapter1.R --sync-only
Rscript chapter1/ledger/make_article_figure_set.R      # must report 0 MISSING

# 5. compile the PDF.
sed 's/≫/>>/g' chapter1/manuscript/manuscript_chap1_EN_native20.md | pandoc -f markdown --citeproc \
  --csl .claude/skills/these-manuscript/apa.csl \
  --bibliography Bib_nathan_corroyez_28Jul.bib \
  --resource-path=chapter1/manuscript \
  -o chapter1/manuscript/manuscript_chap1_EN_native20.pdf --pdf-engine=xelatex \
  -V mainfont="DejaVu Serif" -V monofont="DejaVu Sans Mono" \
  -V geometry:margin=2.4cm -V fontsize=11pt -V linkcolor=blue
```

After step 5, verify no literal `[@key]` or `?@key?` remains in the PDF before delivering.

> **Why step 3 exists.** Two scripts write
> `outputs/figures_pipeline/annex/fig_scanangle_control.png`, and the runner only knows about one of
> them. The file currently on disk was written by `scripts/control_scanangle_vs_lidr.R` (2700×1650) at
> 15:13 on 2026-07-30, **22 minutes after** the sync that shipped the correct 2700×1500 rendering from
> `scripts/fig_scanangle_control_csv.R`. Running `--sync-only` without step 3 would overwrite the
> shipped Fig F1 with the wrong figure. See §7 item 2.

### Two things step 2 does NOT do

* It does **not** run `c1_align_cluster_obs.R` (stage A2). That is the sole writer of the canonical
  observation table `out_files/Chapter1/tables/tab_hobo_perplot_cluster.csv` and is opt-in behind
  `CH1_ALIGN_OBS=TRUE`. Run it only after stage A1 has regenerated
  `tab_hobo_native20_validation.csv`, and only if A1's output changed.
* It does **not** regenerate Fig B1 (`B1_h2_controlled_topheavy.png`). That figure's generator calls
  MuSICA; it is served from the cache listed as asset 8 in §4.

**Never run:** `run_article_v323.R`, `pipeline/run_all.R` and `pipeline/03*`–`11` (all run MuSICA or
target superseded lineages), `c1_hobo_direction.R` (historically clobbered the canonical observation
table), `scripts/make_h2_controlled_topheavy.R` (calls MuSICA).

---

<!-- SECTION 6 -->
## 6. Not traceable

Every number, artifact or input in the chapter whose producer could not be established, **ranked by how
load-bearing it is**. Nothing here was guessed: each entry states what was checked.

Two distinct failure modes are separated below, because they need different fixes:

* **No producer** — nothing in the repo emits the value. Needs a script, or an honest note.
* **Producer contradicts the text** — something emits a value, and it is not the one printed. Needs a
  text correction.

### 6.1 No producer, ranked

| # | item | where | why untraceable | consequence |
|---|---|---|---|---|
| 1 | **The canonical forcing `in_files/FR-Blo_2021_v2.nc` itself** | drives every simulation; described in §2.4, Appendix A, Table H1 | its own attributes show it was built **by a third party from a CSV not on this filesystem**; three in-repo descriptions disagree (station / ICOS flux / unstated) | **the deepest link in the chain is unverifiable.** The temperature channel is corroborated (*r* = 1.000 vs CHS 41); the radiative and aerodynamic channels are not |
| 2 | **`h_sbl` in that forcing** | Table H1, §2.4, data availability | no writer in the repo; does not reconcile with the archived MERRA-2 PBLH at any lag −3…+3 h (peak *r* = 0.537). Its distribution (hourly, min 1.2, median 48.8, max 244.3 m) is a surface-boundary-layer height, about one tenth of a boundary-layer height, not a PBLH | **manuscript CORRECTED 2026-07-31**: §2.4, both Table H1 rows and the data-availability statement no longer claim a half-hourly MERRA-2 series; they describe an hourly surface-boundary-layer height archived in the forcing file, and the MERRA-2 citation is dropped from the Chapter 1 lineage (it belongs to the superseded ERA5-forcing runs). The provenance of the values themselves remains external and unverifiable |
| 3 | **ERA5 surface roughness 0.44 m** | Appendix A line 332; Fig A2 | a bare literal at `c1_wind_profile_correction.R:13` and `R/wind_correction.R:19`. No download, no cached field, no NetCDF variable | propagates through the 38 cached wind-corrected forcings into **every shipped run**. The wind correction is baseline, so this constant is load-bearing |
| 4 | ~~The BLH-insensitivity test (100 m → 5000 m)~~ | Appendix H | **RESOLVED 2026-07-31 by re-running it.** `scripts/c1_blh_insensitivity.R` is now the producer: 8 cLHS plots spanning the four archetypes, reference configuration, two runs each at constant h_sbl = 100 m and 5000 m, everything else identical. Summer-mean ΔT~max~ is **identical on every plot to five decimal places** (max difference 0.000000 °C) against a P4 leaf-area lever of ≈ 0.26 °C, reproducing the deleted 2026-07-23 result. Outputs under `out_files/Chapter1/nc_blh_insensitivity`, table `tab_blh_insensitivity.csv`. Runner stage **B24, gated behind `CH1_RUN_MUSICA=TRUE`** since it is the one stage that simulates. Still unanswered on the model side: `gradm_hsbl` / `gradc_hsbl` would tell whether the insensitivity is physical or structural (the `dim0_driver.f90:509` floor at the reference height), but they are absent from the default `history_variables` and requesting them through `extra_setup` makes MuSICA reject the namelist (STOP 6). This does not weaken the ΔT~max~ result, which is measured directly |
| 5 | ~~VIF 1.3 (fCover) and VIF ≈ 9.5 (sd_height / rumple)~~ | Fig C1 caption | **RESOLVED 2026-07-31** | `scripts/make_residual_vs_topo.R` now computes and prints them and writes `tab_residual_vif.csv` (fCover 1.26, sd_height 9.53, rumple 9.41) |
| 6 | ~~"brackets the 15:00 peak"~~ | Fig B2 caption | **RESOLVED 2026-07-31** | `scripts/make_vertical_profiles_native.R` now derives the diurnal cycle of the 1 m sub-canopy temperature over the same summer months from the full-model runs, prints the peak hour and writes `tab_diurnal_cycle_1m.csv`. It peaks at **15:00 on the forcing clock (23.11 °C)**, so the 10:00–16:00 window does bracket it: the caption was right, it simply had no producer |
| 7 | **The monthly roughness range 0.438–0.443 m** | Appendix A line 333 | exists only as a **source-code comment** | prose detail |
| 8 | ~~The promotion `clhs_sample_native20_fpca_floor05.rds` → `clhs_sample_native20_floor05.rds`~~ | §1 Finding 1-ii | **RESOLVED 2026-07-31** | `scripts/redo_clhs_native20_fpca.R` now performs the promotion at its end, idempotently, and **refuses** to overwrite a canonical sample that differs from the FPCA one, since every native20 simulation on disk was run from the canonical file |
| 9 | **`typical_plots_per_cluster.csv`** | Fig 1b (the four representative plots) | **PARTLY RESOLVED 2026-07-31.** `scripts/c1_typical_plots_check.R` (stage B23b) verifies every trait, coordinate and archetype label against the canonical design sample, and the `pid` numbering (row index in the 400-plot sample) reproduces exactly, so the file can no longer drift from the design. The **selection rule** is still unrecovered: no centroid distance reproduces these four plots, over z-scored or min-max traits, in (LAI, *H*max, fCover), +VCI, +FPC1-3, or the FPC space alone, so `dist` is provenance metadata, not a reproducible quantity | most plausibly the selection also required an extracted point cloud, which Fig 1b needs and which lives on an external drive. Documented in the script header |
| 10 | **The DEM's own origin** | Appendix C covariates | the TWI / northness producers are **commented out** in `NC_Full/.../3_calculate_lidar_metrics.R`, and the DTM VRT they read is missing | whether the DEM is ALS-derived or a national product is recorded nowhere |
| 11 | **The CHS 41 gap-fill against ERA5** | Appendix A line 319 | `MetHor*.txt` are raw; no chain script performs the gap-fill | a stated forcing-construction step |
| 12 | **Station altitude 127 m and reference height 1.50 m** | §2.5 (the ΔTmax reference) | `prep_site_forcing/CHS41-Blois.cfg` **self-flags both** as "best-guess, TO CONFIRM" | 1.5 m is the macroclimate reference height for every ΔTmax |
| 13 | **ALS survey metadata** (72.6 pts m⁻², 4.5 cm precision, flight parameters) | §2.1.2 | provider-supplied, cited to Gril 2023 | descriptive only |
| 14 | **R and package versions** | line 715 | no `renv.lock`, no `sessionInfo()` dump, no manifest anywhere | the environment is unreproducible |
| 15 | **The five data-availability placeholders** | **DEFERRED BY THE AUTHOR 2026-07-31, do not re-raise:** out of scope for now. lines 708–717| repository DOI, MuSICA distribution source, data provenance/contact, software versions | must be filled before submission |

### 6.2 Values that reproduce but are not persisted

These are **traceable** — each recomputes exactly from a canonical file — but no table holds them, so an
auditor must re-run the generator or recompute by hand. Listed because they are the natural target of a
persistence pass.

| value | appendix | recomputes from |
|---|---|---|
| the nested *R*² set: 0.79 full, 0.77 adjusted, 0.79 quantity-only, 0.67 VCI-structure, 0.44 profile-structure, +0.13 / +0.35 | D | `tab_hobo_perplot_cluster.csv` |
| "26 of 53 loggers have one-sided LAI ≥ 4" | D | `in_files_native20/lai_z1_res_10_m.tif` at the logger points |
| PC1 explains ≈ 75% of the variance | E | `clhs_sample_native20_floor05.rds` |
| all per-archetype trait correlations (0.79 / 0.97 / 0.93 / 0.95) | E | same |
| paired fractions 63 / 42 / 21 / 82% | G | `sensitivity_perplot_units_v2.csv` — **no repo script emits these** |
| "median predicted 0.003 °C" (noise-floor linear extrapolation) | G | `median(abs(LAI_up_dt_all * 0.02))` = 0.00282 — **no producer**. Holds pooled; per-archetype P4 predicted −0.0101 vs observed −0.0082, a 23% gap |
| dispersion "widens by a factor of 1.4 to 1.8" | Supp. | **IQR** ratios hot/all = 1.76 / 1.40 / 1.46 / 1.72 → fits exactly. **SD** ratios would give 1.31 for P2, i.e. "1.3 to 1.8". The statistic is not stated in the text |
| adjusted *R*² = 0.69 and *r* = −0.84 | C | `tab_residual_vs_topo.csv` |
| ERA5 bias +0.96 °C, RMSE 2.06, *r* = 0.92, +3.6 °C peak | A | `in_files/{FR-Blo_2021_v2,musica_in_Blois}.nc` |
| S1 *R*² 0.91 / 0.89; S3 *r* 0.35 / 0.63; S5 0.68 / 0.70; dispersion 1.4–1.8× | Supp. | `hot_extract.rds`, `sensitivity_perplot_units_v2.csv` |

### 6.3 Producer contradicts the text — corrections needed

| # | claim | what the source actually gives | severity |
|---|---|---|---|
| 1 | **Fig B1 caption: "v3.2.3 iter, CHS 41 station forcing"** | **RESOLVED 2026-07-31 (two separate defects).** (a) LINEAGE: the figure ran on `CFG` defaults, i.e. the legacy binary, ERA5 forcing and no ABL. `scripts/make_h2_controlled_topheavy.R` now defaults to the chapter lineage (v3.2.3, FR-Blo_2021_v2.nc, `abl_flag=iter`, baseline wind correction). (b) LAI UNITS: `LAI_LEVELS` are passed to `lai_fn`, which `run_musica_one` DOUBLES, so the old `c(6, 12)` were one-sided 6 and 12, not the "one-sided 3 and 6" the caption claimed, and one-sided 12 is off-manifold (design 0.50 to 7.11). The levels are now explicitly ONE-SIDED and capped at 6. Re-swept at one-sided LAI 1 to 6 (49 runs): H2 holds in DIRECTION at LAI 2 to 5, a mid-canopy balance wins at LAI 6, and the lever is SECOND-ORDER throughout (span 0.16 to 0.45 °C). Grid check done and clean: the 10-layer mapping shifts every centre of mass by a constant 0.030. Abstract, Results, Appendix B intro, Fig. B1 caption, Discussion H2 paragraph and the Appendix E passage all rewritten. **RESOLVED 2026-07-31 by switching the figure to the chapter lineage** (author decision: everything must be v3.2.3 + iterative ABL + CHS 41 station forcing). `scripts/make_h2_controlled_topheavy.R` now defaults to `H2_LINEAGE=chapter` (v3.2.3 binary, FR-Blo_2021_v2.nc, `abl_flag=iter`, baseline wind correction at the fixed *H*max = 25 m); `legacy` is kept only to re-read the superseded cache and must not ship. The 14 new runs are in `out_files/H2_controlled_topheavy_chapter/`, which is now the pre-flight asset. **The result changed and the manuscript was rewritten**: the response is non-monotonic, peaking at a mid-canopy centroid (~0.52) and weakening toward both extremes (LAI 12: −0.79 bottom, **−3.11 centered**, −0.50 top), and the lever GROWS with leaf area (span 0.45 °C at LAI 6, 2.61 °C at LAI 12) instead of saturating. H2 is rejected on DIRECTION, not strength. Abstract, Results, Appendix B intro, Fig. B1 caption, the Discussion H2 paragraph and the Appendix E observational passage were all updated. **CAPTION CORRECTED 2026-07-31, BUT A SUBSTANTIVE ISSUE IS NOW OPEN FOR THE AUTHOR.** Confirmed: `scripts/make_h2_controlled_topheavy.R` ran on `CFG` defaults, i.e. the legacy binary, the ERA5 `musica_in_Blois.nc`, and **no** `abl_flag`, so the caption's "v3.2.3 iter, CHS 41 station forcing" was false on all three counts; it now states the actual configuration. The script gained an `H2_LINEAGE=chapter` mode that re-runs the same 14 configurations on the chapter lineage (v3.2.3, FR-Blo_2021_v2.nc, iter) into a separate directory. **The conclusion does not survive.** Legacy (shipped): monotonic, ΔTmax falls steadily as mass rises, LAI6 +0.75 to +0.16 (span 0.59 °C), LAI12 −0.06 to −0.33 (span 0.27 °C), which is the "top-heavy cools more, second-order, saturates with LAI" reading. Chapter lineage: **non-monotonic**, peak cooling at the CENTERED profile, LAI6 −0.03 → −0.49 → −0.24 and LAI12 −0.56 → **−2.44** → −0.43, a span of 1.9 °C, i.e. an order of magnitude larger and no longer second-order. Signs differ too (legacy LAI6 warms throughout, chapter LAI6 cools throughout). This bears on a locked framing, so it is NOT applied here: the author decides whether to swap Fig. B1 to the chapter lineage and rewrite the H2 paragraph accordingly. its generator uses `in_files/musica_in_Blois.nc` (**ERA5**) and the **Nov-2024 legacy binary** (md5 `53072278`), which `R/config.R` says shifts ΔTmax by up to 1.4 °C. The 14 sims are cached, so **re-running cannot fix the caption**| **high** — misstates a model configuration; values (0.59 / 0.27 °C) are exact |
| 2 | **Table H1 line 618, 620 + Data availability 711–712: "half-hourly MERRA-2 PBLH"** | **RESOLVED 2026-07-31.** Both Table H1 rows, §2.4 and the data-availability statement now describe an hourly surface-boundary-layer height archived in the forcing file; the `@gelaro` MERRA-2 citation is dropped from the Chapter 1 lineage. Verified: 0 occurrences of "MERRA" in the compiled PDF. MERRA-2 is downloaded **hourly**; the forcing's `h_sbl` is ~10× too small (200 m vs 1663–2070 m for the same hours) and carries no MERRA-2 attribute; line 618 contradicts line 619 within the same table| **high** — three locations |
| 2b | **Table H1: "Simulation period 1 June to 30 September 2021"** | **RESOLVED.** Table H1 now reads "Full year 2021 (analysis window 1 June to 30 September)". `out.nc` shows MuSICA integrates **all 8759 hourly steps of 2021**. JJAS is a post-processing extraction window| **high** — misdescribes the run |
| 2c | **Table H1: "Initialized on 1 June (no multi-year spin-up)"** | **RESOLVED.** Table H1 now reads 1 January, matching `INIT_*` and the five months of within-year spin-up. initialization is **1 January** (`restin_filename = ""`), giving ~5 months of within-year spin-up before the metrics window. The parenthetical is right, the date is wrong — and the §4.2 argument about the topsoil staying near saturation depends on it| **high** |
| 2d | **Table H2: "Retention shape parameter *m* = 0.5"** | **RESOLVED 2026-07-31.** Relabelled "Exponent *m*", stated as a namelist value the model does not read under this flag combination: `RETENTION_CURVE_MODEL_FLAG = 3` fixes retention *m* = 1 internally, and `bigm_soil` enters only the Brooks and Corey branch, which `HYDRAULIC_COND_MODEL_FLAG = 2` does not select (verified in `mo_soil_transfer.f90`). 0.5 is `BIG_M_SOIL`, the **unsaturated hydraulic-conductivity** exponent. With `RETENTION_CURVE_MODEL_FLAG = 3` the namelist fixes retention ***m* = 1**| medium — a mislabelled parameter |
| 2e | **§H line 610: "the complete namelist files are archived with the code"** | **DEFERRED BY THE AUTHOR 2026-07-31, do not re-raise:** out of scope for now. `musica.nml`, `musica_soil.nml` and `musica_veg1.nml` are all **untracked in git**, as are `run_chapter1.R`, `c1_attribution_bootstrap_units.R`, `c1_units_bootstrap_v2.R` and `c1_noise_floor.R`| medium — nothing is currently under version control to archive |
| 2f | **Table G2 caption: "the same design and the same plots as Table G1"**, with line 550's "4000 resamples" | **RESOLVED.** `c1_units_bootstrap_v2.R` now uses `NB <- 4000`, matching the caption. G2 used **NB = 2000** in a separate RNG draw. Design and plots *are* identical| low — but it implies a shared bootstrap that does not exist |
| 2g | **Appendix G: P4 effects "an order of magnitude above" P1** | **RESOLVED 2026-07-31.** Rewritten: cover is an order of magnitude weaker in P1 (factor 12.7), and leaf area is stated as **reversing sign** there (+0.049 °C), not merely weakening. cover is 12.7× ✓; **leaf area is 5.4×**, and P1's leaf-area effect is **+0.049 — a sign reversal**, not a weaker cooling| low — overstated |
| 2h | **Appendix G: height "never separates from zero by more than four thousandths"** | **RESOLVED 2026-07-31.** The claim is now explicitly about the **median** effects, and the widest interval (−0.006 °C in P3) is given. true of the **medians** (max 0.004) but not of the CIs the wording invites: P3 *H*max +1 m has `lo` = −0.006| low |
| 2i | **Table H2: "*g*₁ = 10"** and **"leaf xylem potential −1.3 MPa"** | **RESOLVED.** `in_files/musica_veg1.nml` gives `GS_SLOPE = 10.`, `GS_HX_HALF = -1.3`, `GS_HX_SHAPE = 2.6`, matching the text. `GS_MODEL_FLAG = 2` is **Ball-Berry**, so *g*₁ is the Ball-Berry slope *m*, not a Medlyn USO *g*₁; and `PHOTOSYNTHETIC_LIMITATION_FLAG = 3` makes −1.3 MPa a **predawn** potential updated daily, not an instantaneous leaf value| low — both are naming precision |
| 2j | **Table G1 caption: "steps are truncated at the trait bounds"** | **RESOLVED 2026-07-31.** The caption now discloses that cover −10 points rests on 47 plots in P1 and 97 in P2, all other rows on 100. true, but it does not disclose that **P1 fCover −10 uses only 47 of 100 plots** and P2 fCover −10 uses 97; all other rows use 100| low — an omission, not an error |
| 3 | **Appendix C: "elevation and northness, *p* ≈ 0.05"** | **RESOLVED (stale entry).** The text already reads elevation *p* = 0.038 and northness / slope / TWI *p* = 0.93, 0.29, 0.98, which is exactly what `make_residual_vs_topo.R` recomputes. elevation 0.0381 ✓; **northness partial *p* = 0.931**, simple *p* = 0.551| medium — elevation alone carries the claim |
| 4 | **Appendix C: residual Moran's *I* = 0.05, *p* = 0.13** | **RESOLVED 2026-07-31.** `scripts/make_moran_validation.R` was re-pointed from the superseded z05 branch to the canonical `tab_hobo_native20_validation.csv`, the same field `make_residual_vs_topo.R` decomposes, and writes `out_files/Chapter1/tables/tab_moran_validation.csv`. Appendix C now quotes what it returns: observed *I* = 0.097 / *p* = 0.053, modelled 0.213 / 0.001, residual 0.049 / 0.098. Two wording changes follow: the observed structure is marginal rather than significant, and the model is **more** autocorrelated than the observations, not merely faithful to them. the committed **z05** table gives *I* 0.035 / *p* 0.164; a native20 recompute gives 0.049 / *p* 0.098. **No lineage yields *p* = 0.13.** The observed pair (0.10, 0.04) is from the z05 table| medium — the qualitative conclusion (no detectable residual structure) holds either way |
| 5 | **Appendix A: factor "0.41 to 0.47" over 13–33 m** | **RESOLVED 2026-07-31.** `wind_factor()` over 13-33 m gives **0.411 to 0.484**, so the upper bound was wrong; the text now reads 0.41 to 0.48. at 13–33 m the factor is **0.48–0.41**. The printed 0.47–0.41 corresponds to **15–35 m**, the range the script itself prints| low — a range-pairing slip |
| 6 | **Appendix D: "raises the fit by 0.13"** | **RESOLVED 2026-07-31.** Changed to 0.12, consistent with the *R*² values the text actually prints (0.79 − 0.67). 0.13 differences the *VCI-full* model (0.7958), but the manuscript only ever prints 0.79, so a reader gets **0.12**| low — reword to "about 0.12 to 0.13" |
| 7 | **Fig B2 caption: "10:00 to 16:00 on the forcing clock (solar time, UTC+1)"** | **RESOLVED (no error).** The forcing time axis carries no timezone, so its labels ARE the forcing clock; parsing as UTC and taking `hour()` selects 10-16 on that clock, which is what the caption says. The script header's "applied on UTC" note is what misleads, not the caption. the window is applied **on UTC**; the generator's own header calls this "a deliberate simplification"| low — a one-hour framing error in a diagnostic figure |
| 8 | **Appendix E line 472: "Collinearity is strong everywhere"** | **RESOLVED.** The quoted per-archetype LAI-fCover values (P2 0.97, P3 0.93, P4 0.95, P1 weakest at 0.79) are exactly what `make_trait_collinearity_native20.R` prints. per-archetype \|*r*\| bottoms out at 0.03–0.14 (LAI–*H*max) and LAI–VCI is **−0.19** in P4. The following colon scopes it to LAI–fCover, which is defensible | low — wording |
| 9 | **Data availability: "regenerates … from the raw inputs"** | **RESOLVED 2026-07-31** (see the reworded statement in the manuscript). regenerates from **archived MuSICA outputs**; Fig B1 is cached, three tables are hand-pasted, two have no generator| medium — an overclaim in the reproducibility statement |
| 10 | **Appendix F: per-plot ⟨sec θ⟩ "1.02 to 1.09"** | **RESOLVED 2026-07-31.** Verified against the canonical design sample, whose `sec_theta` column runs [1.0161, 1.0956], so "1.02 to 1.09" is right; the divergent table cited before is a different lineage. correct against `Blois_lad_z05_sacorr_r25.csv` ([1.0188, 1.0905]), but `out_files/Chapter1/tables/tab_scanangle_control.csv` gives [1.0184, **1.0781**]| low — an auditor trap, not an error |
| 11 | **Appendix F: scan angles "±33°"** | **RESOLVED 2026-07-31.** The text now gives the asymmetric range, −32° to +33°. max **+33.0°**, min **−32.0°** — asymmetric| cosmetic; the script's own header acknowledges it |

---

<!-- SECTION 7 -->
## 7. Broken and fragile links

Ranked by how likely each is to silently corrupt a shipped artifact. Everything here was verified
against the bytes on 2026-07-31.

### 7.1 Live hazards

**1. `run_chapter1.R --sync-only` would currently ship the WRONG Figure F1.**
`outputs/figures_pipeline/annex/fig_scanangle_control.{png,pdf}` has **two writers**:
`scripts/fig_scanangle_control_csv.R` (9×5 in → 2700×1500; the runner's B15 and the ledger's named
script) and `scripts/control_scanangle_vs_lidr.R` (9×5.5 in → 2700×1650; needs the external drive).
The shipped `F1_scanangle_control.png` is `6082d379`, 2700×1500 — the CSV version, which is what the
caption reproduces. The file now on disk is `55cc73e8`, 2700×1650, written at 15:13, **after** the
14:51 sync. `control_scanangle_vs_lidr.R` is named in **no** documentation — not the runner, not the
ledger, not `RUN_CHAPTER1.md`. *Remedy:* run `scripts/fig_scanangle_control_csv.R` before any sync
(now step 3 of §5); longer term, give the two scripts distinct output paths.

**2. Two writers of `tab_hobo_native20_validation.csv`, and it feeds three shipped figures.**
`c1_hobo_native20_validation_regen.R:39` is stage A1 (the runner's). `c1_hobo_native20_validation.R:67`
writes the same path and is **not** in the runner — it applies the `sec_theta` correction from
`Blois_lad_z05_sacorr_r25.csv` and reads the rasters through `load_lidar_rasters`, i.e. a different
recipe. That table feeds **Fig 6, Fig I1** (via `c1_fig5_j1_convB.R:21`) **and Fig C1** (via
`make_residual_vs_topo.R:31`), plus stages A2, A9, A11. A stray run of the non-runner script would move
all three figures with nothing detecting it: the file carries no lineage stamp and the runner has no guard.

**3. `scripts/make_residual_vs_topo.R:47` silently degrades if `~/Documents/NC_Full` is unavailable.**
`cov_files <- cov_files[file.exists(cov_files)]` drops missing covariates without erroring, so
Appendix C's adjusted *R*² = 0.69 would change and the figure would still render. `~/Documents/NC_Full`
is 171 GB, sits outside both the repo and the external drive, and is documented nowhere. All six
rasters verified present today. *Remedy:* replace the filter with a `stopifnot`.

**4. `scripts/make_h2_controlled_topheavy.R` (stage B21) will invoke MuSICA if its cache is cleared.**
Its guard at L83–87 re-extracts from the 14 cached NetCDFs in `out_files/H2_controlled_topheavy/`, but
falls through to `run_musica_one()` if any is missing or under 1000 bytes. The runner advertises every
stage as simulation-free. Compounding this, `RUN_CHAPTER1.md` §0 lists that directory as a
never-delete asset while `run_chapter1.R`'s `ASSETS` pre-flight checks **the cached PNG instead** — so
deleting the NetCDFs passes pre-flight and then triggers a MuSICA run mid-chain.

**5. No manuscript table is mechanically synced.** `SYNC` in `run_chapter1.R` copies PNGs only. Table 1,
G1 and G2 are hand-pasted markdown; H1 and H2 have no generator at all. All four verified correct today
(Table 1 byte-identical to `tab_table1_full.md`; G1 = A10 to the last digit; G2 = A7's `sl/all` block;
H1/H2 = the namelists). But re-running A3, A7, A8 or A10 will silently desynchronise the manuscript,
and neither the runner nor the ledger performs a diff.

**6. Tables G1 and G2 come from different scripts with different bootstrap counts.** G1 from
`c1_attribution_bootstrap_units.R` (B = 4000), G2 from `c1_units_bootstrap_v2.R` (NB = 2000). Both are
in the runner and both write plausibly named files. `tab_attribution_units_bootstrap_v2.csv` **also
contains a ΔTmax block that is not Table G1** and differs from it in four CI bounds. Anyone
regenerating Table G1 from the "v2" file will produce subtly wrong numbers.

**7. `c1_fig3_units.R` is two runner stages behind one env switch.** B2 (`PERIOD=all METRIC=dt` → Fig 4)
and B18 (`PERIOD=hot METRIC=dt` → Fig S2) share one script and one env-conditional output base. Running
one without the other leaves Fig 4 and Fig S2 on different extractions of the same design, with nothing
visible to distinguish them.

### 7.2 Stale documentation — existing notes contradicted by the bytes

Every item below is asserted by current documentation and is **no longer true**. Left uncorrected, a
reader hitting the conflict trusts neither document. Re-verified 2026-07-31.

| where | what it claims | actual state |
|---|---|---|
| `RUN_CHAPTER1.md` §6.1 | the shipped Fig 3/4 attribution is a stale rendering (`97f0c89f`) vs the generator's `aa587326` | **resolved** — shipped `Fig4_attribution.png` is now `aa5873267f9546a36dd8504aaf99bb4e`, identical to `out_files/Chapter1/figures/Fig3_attribution_units.png` |
| `RUN_CHAPTER1.md` §6.2 | `scripts/c1_operating_point_main.R:9` reads the older `sensitivity_perplot_units/part_*.csv` | **fixed** — it carries a "CANONICAL SOURCE 2026-07-29" header and reads `sensitivity_perplot_units_v2.csv` |
| `RUN_CHAPTER1.md` §6.3 | `scripts/make_residual_vs_topo.R:26` reads the superseded `tab_hobo_frblo_validation.csv` | **fixed** — L26–27 is now a "LINEAGE FIX 2026-07-29" comment and L31 reads `tab_hobo_native20_validation.csv`. (`c1_fig5_forward_2metric.R:8` still reads frblo, but that script is in neither the runner nor the ledger — it is dead.) |
| `RUN_CHAPTER1.md` §6.7 + line 179 | "**Appendix F is unreproducible — highest-severity item, blocks submission**"; no script emits 0.95 / 0.91 / 0.79 | **stale and mis-keyed.** That note targets the SHORT doc. In the full manuscript the correlation appendix is **E**, states 0.94 / 0.73, and recomputes exactly (0.9354 / 0.7317) from the canonical rds, as do all four per-archetype values. See §3 Finding E-i |
| `RUN_CHAPTER1.md` §6.8 | `c1_pca_recover_clusters.R` is unseeded and reads `clhs_sample_floor05_v2.rds` | **fixed** — `set.seed(42)` at L23, native20 rds at L24 |
| `chapter1/ledger/make_article_figure_set.R:55-64` and `run_chapter1.R:145-148` | the F1/E1 figure is "PROVENANCE UNRESOLVED", byte-identical to `FigAnnex_corrplot_traits.png` from `pipeline/12_corrplot_traits.R` | **contradicted by md5.** `FigAnnex_corrplot_traits.png` = `f8eabf8b…`; the shipped `E1_trait_collinearity.png` = `2faa73a8…` = `fig_trait_collinearity_native20.png` from `scripts/make_trait_collinearity_native20.R`. The stale note now mislabels a correct figure and invites someone to "fix" a working chain |
| `RUN_CHAPTER1.md` §0 | `out_files/musica_native20_forward/` is "16 coalition sets × 53 loggers" | **60** NetCDFs per set (960 total). 60 were simulated; `R/config.R:19` drops 7 (`ids_to_remove`) leaving the 53 that are reported |
| `scripts/c1_make_table1.R` header | it reads the cLHS design sample | it reads `sensitivity_perplot_units_v2.csv` |
| `R/cluster_relabel.R` header | per-archetype LAI 2.4 / 6.0 / 7.7 / 10.8 | the native20 sample gives 1.80 / 3.15 / 3.70 / 5.03. The *mapping* is correct; only the comment misleads |
| `scripts/redo_clhs_native20.R` header | it produced `clhs_sample_native20_floor05.rds` | the shipped file is byte-equal to the **FPCA** build (§1 Finding 1-ii) |
| `c1_wind_profile_correction.R` header L4-5 | Blois roughness "0.40 m", range "0.30–0.58 m", correction "**NOT applied in the simulations**" | the code uses 0.44 / 0.438–0.443 and the correction **is** baseline since 2026-07-10 (§3 Finding A-ii) |
| `scripts/compute_lad_z05_sites.R:26` | "UNTESTED: the raw LAS were offline" | the drive is mounted and all six composed paths resolve |

### 7.3 Traps that are inert today

**8. `pipeline/00_config.R:44` pins `PIPE$CLUSTER_SAMPLE` to the superseded
`clhs_sample_floor05_v2.rds`.** Four chain scripts source that config (`c1_hourly_scatter.R`,
`c1_hot_extract.R`, `c1_radius_rescore.R`, `c1_hobo_native20_validation_regen.R`, plus
`c1_forward_native20_sims.R`), so the wrong sample is *in scope* — but none of them uses it; every
native20 script does an explicit `readRDS` of the native20 sample. Reading `00_config.R` alone gives
the wrong design. "Fixing" it would break the legacy branch.

**9. `musica.nml:49` still sets `FORCING_FILENAME = "in_files/musica_in_Safran_Blois_2021.nc"`** — a
SAFRAN file. The chain overrides it programmatically (`run_musica_one` passes `forcing_filename` in
`setupctl`), so shipped runs are unaffected. But any manual `./musica` invocation from the repo root
would silently use SAFRAN, and the namelist contradicts Table H1 for anyone reading it directly.

**10. Two divergent copies of each namelist.** Root `musica_soil.nml` (`3e151155`) vs
`in_files/Blois/in_files/musica_soil.nml` (`2778f1bc`) — comment-only differences, functionally
identical. Root `musica_veg1.nml` (`c6f2daac`) vs the referenced copy (`84e8b07b`) — differ by exactly
one line, `CANOPY_HEIGHT_TOP = 1.2` versus the `CANOPYHEIGHT` template placeholder that `rmusica`
substitutes per plot. Only the `in_files/Blois/in_files/` pair is referenced by `musica.nml`, and only
those are what Tables H1/H2 transcribe. **No namelist is archived alongside the run outputs**
(`find out_files -name '*.nml'` is empty), so Table H1's provenance is inferred from the current repo
state, not proven from the run artifacts.

**11. Three MuSICA binaries, and `CFG$musica_cmd` defaults to the wrong one.** `R/config.R:17` points
at `/home/corroyez/Documents/musica/musica` (md5 `53072278`, Nov 2024 legacy). Every native20 runner
explicitly overrides this to `in_files/model-3.2.3/musica` (md5 `8b5139e0`). **Fig B1's generator does
not override it** — see §3 Finding B-i. `R/config.R:12-16` documents that substituting one for the
other shifts ΔTmax by up to 1.4 °C.

**12. Non-determinism: none found.** `c1_pca_recover_clusters.R` seeds 42; `c1_fig1_gril_native.R` seeds
both jitters and its subsample (frozen by the cache anyway); both bootstrap scripts seed 1;
`c1_sensitivity_units_extract2.R` seeds its sub-sample; `mclapply` in A3–A5 returns in index order. The
only artifact not reproducible byte-for-byte is Fig B2, whose montage depends on the local ImageMagick
build — and it is **pixel-identical** (`compare -metric AE` = 0). Cosmetic.

**13. Leftover duplicates in the figures directory.** `chapter1/manuscript/figures/article_v323/` holds
77 PNGs for 25 cited figures: `SYNC` writes each source under **both** the short-doc and full-doc names,
so 17 byte-identical pairs/triples exist by design (e.g. `A2_` = `B2_` = `G1_wind_profile_correction.png`).
Two are **not** by design: `F2_pca_recover_clusters_vci.png` (394.8 kB) no longer matches
`Fig3_pca_recover_clusters.png` (382.4 kB) although SYNC writes both from the same source, and
`outputs/figures_chap1/FigX_era5_station_bias.png` is a stale twin of `B1_era5_station_bias.png` that
nothing writes any more.

**14. Dead output.** `scripts/make_trait_vertical_gradient.R:236-237` writes a third base path
`fig_trait_vertical_gradient_effects.{png,pdf}`, an exact duplicate of the `_norm` variant that nothing
reads and nothing syncs. Harmless, but it doubles B10's render time and invites a future SYNC entry to
pick the wrong file.

---

## Figure B3 rebuilt without a baseline (2026-08-02)

Figure B3 was a forward-coalition decomposition: every increment was a step from ONE common canopy, the 400-plot design mean (one-sided LAI 3.42, *H*max 23 m), evaluated at the 53 logger pixels (*n* = 8, 12, 13, 20). Two consequences made it hard to read. Increment size tracked the distance between an archetype's measured trait and that baseline rather than its sensitivity, so the profile appeared to lead in P2 and P3 simply because their leaf areas sit on the baseline. And in P1 the step labelled ΔLAI actually REMOVED two thirds of the canopy, so it cooled while raising sub-canopy wind, which the caption had to explain away.

It is now built by `scripts/make_trait_vertical_gradient_perplot.R` on the on-manifold design: each of the 400 cLHS plots is perturbed around its own measured canopy by the fixed native steps of Section 2.6, the same design as Figure 4. *n* rises to 100 per archetype at every level, the four columns receive the same physical step, and the P1 signs become coherent (ΔLAI +0.036 °C, −0.059 m s⁻¹). Two construction details are load-bearing: profiles are interpolated onto a common relative-height grid before averaging, since each plot has its own grid; and the axis stops at z = 1.4, because above 1.45 the per-plot grids stop covering and *n* collapses to 1 in P4. Pure re-extraction, no MuSICA. The old script is kept with a SUPERSEDED banner.
