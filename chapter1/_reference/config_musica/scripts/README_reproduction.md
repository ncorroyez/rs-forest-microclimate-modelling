# Scripts & fichiers utilisés — comparaison v3.2.0 vs v3.2.3 (reproductibilité)

Tout ce qui sert à produire la comparaison : scripts de run/analyse, helpers R, configs, données d'entrée.
Les **binaires + namelists + forçages** sont dans `../v3.2.0_legacy/` et `../v3.2.3/` ; les **figures/tables** dans `../../comparaison_versions/`.

## ⚠️ Chemins
Les scripts utilisent des **chemins relatifs à la RACINE du projet** (`out_files/`, `in_files/`, `R/`, `pipeline/00_config.R`).
Ce dossier est destiné au **dépôt / à l'inspection** (transparence). Pour **rejouer**, replacer les scripts à la racine du projet,
avec `R/` = `scripts/R_helpers/`, les `inputs/` dans `in_files/` (+ forçages), et les binaires accessibles.

## Dépendances
- **R ≥ 4** + packages : `rmusica`, `musica.tools` (wrappers MuSICA, fonction `callmusica` / `run_musica_one`), `ncdf4`, `data.table`, `dplyr`, `lubridate`, `ggplot2`, `terra`, `sf`, `here`.
- **Binaires** : `../v3.2.0_legacy/musica` (md5 5307…) et `../v3.2.3/musica` (8b5139e0).
- **pandoc** (pour le DOCX).

## Config
- `Chapter3_config.R` → `CFG_C3` (musica_cmd, forcing_file…).
- `pipeline_00_config.R` (= `pipeline/00_config.R`) → `CFG` (forcing_file, hobo_temp_csv, hobo_geojson, ids_to_remove, date_seq, in_dir, agg_factor).
- `R_helpers/` = toutes les fonctions sourcées (`musica.R` = `run_musica_one`/`get_tair_at_z`/`extract_macro_daily`/`build_era5_hourly`/`calc_phenology` ; `lad.R` = `make_lad_real`/`make_lad_uniform` ; `cluster_relabel.R` ; `io.R`, `forest.R`, `scenarios_c3.R`…).

## Données d'entrée (`../inputs/`)
| Fichier | Rôle |
|---|---|
| `clhs_sample_floor05_v2.rds` | les **400 plots cLHS** (traits + profils LAD), échantillon canonique |
| `Blois_lad_z05_r25.csv` | profils LAD **0,5 m** aux 53 HOBO (validation) |
| `Blois_data_temperature.csv` | températures **HOBO** horaires (obs) |
| `data_Blois_utm31n.geojson` | positions des HOBO |
| forçages | `../v3.2.0_legacy/forcing/musica_in_Blois.nc` ; `../v3.2.3/forcing/musica_in_Blois_pblh.nc` (avec `h_sbl`) |

## Scripts → rôle / figure
| Script | Rôle | Figure |
|---|---|---|
| `c1_sensitivity_perplot_chunk.R` | sensibilité per-unit 400 plots, **v3.2.0** (8 simus/plot) | 6, 9 |
| `c1_sensitivity_perplot_chunk_iter.R` | idem **v3.2.3 iter** (vraie BLH + fix phéno) | 6 |
| `c1_version_iter_compare.R` | importance per-unit v3.2.0 vs iter | 6 |
| `c1_metrics6_perplot.R` + `c1_metrics6_compare.R` | 3 métriques × 2 périodes (ré-extraction) | 7 |
| `c1_cmp_scatter_alltemps.R` | scatter horaire toutes-T | 8 |
| `c1_cmp_scatter_metrics.R` | ΔTmax & slope par plot (tous/chauds) | 9 |
| `c1_sensitivity_archetypes.R` | sensibilité « P1-P4 only » (profils moyens) | 10 |
| `c1_hobo_iter_validation.R` | validation HOBO en mode iter (53) | 5 |
| `c1_cmp_v320_v323iter.R` | validation + importance (figures) | 5, 6 |
| `c1_cmp_stats.R` | stats détaillées + tests (r±IC, Wilcoxon, bootstrap Δr) | 11, 12 |
| `get_merra2_pblh.sh` + `build_forcing_pblh.R` | télécharge MERRA-2 PBLH → injecte `h_sbl` dans le forçage v3.2.3 | — |
| `run_perplot_v323.sh`, `run_perplot_iter.sh`, `run_metrics6.sh` | lanceurs chunkés (12 process parallèles) | — |

## Pré-requis du mode yoyo (v3.2.3 `ABL_flag='iter'`)
1. Phéno du jour précédent (**2020/366**, ajoutée par les scripts iter) — sinon `STOP 16`.
2. Variable **`h_sbl`** dans le forçage (`musica_in_Blois_pblh.nc`) — sinon `STOP 6`.
