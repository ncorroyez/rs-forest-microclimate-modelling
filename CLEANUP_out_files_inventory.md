# out_files — inventaire pour ménage disque (2026-07-07)

Total `out_files` = **211 G**. Disque : **1.9 T, 98 % plein, ~53 G libres.**
Aucune suppression effectuée — ceci est un état des lieux pour que tu décides.
Colonnes : taille | date (mtime dossier) | rôle présumé | reco.

Reco : 🟢 DELETE-safe (backup/buggy explicite) · 🟡 OBSOLETE (0 réf récente, superseded) ·
🟠 REVIEW (ton jugement) · 🔒 KEEP (v3.2.3 courant ou référencé récemment).

---

## 🟢 DELETE-safe — backups / buggy explicites (≈ 62 G)

| Taille | Date | Dossier | Rôle | Reco |
|---|---|---|---|---|
| 47 G | 07-06 | `Chapter1/_WRONGLAI_bak` | runs 1×LAI (bug ×2) + v3.2.0 — invalides. Contient la seule copie `DYN_HYBRID_s12` mais inutilisable en v3.2.3 | 🟢 |
| 12 G | 07-06 | `Chapter1/_BASEFORC_bak` | backup ancien forçage (avant station forcing) | 🟢 |
| 2.1 G | 03-25 | `musica_option_b_ERA5_results_old` | « old » | 🟢 |
| 1019 M | 07-03 | `musica_hobo_radius_X2bak` | backup X2 (doublon de `musica_hobo_radius`) | 🟢 |
| 986 M | 07-03 | `Chapter1/nc_sens_frblo_FIXEDSTEP_bak` | backup | 🟢 |
| 516 K | 06-12 | `Chapter1/nc_shapley2x_v2_BAD` | « BAD » | 🟢 |
| 732 K | 06-08 | `Chapter1/lai_prep_COMMON_BACKUP` | backup prep | 🟢 |

---

## 🟡 OBSOLETE — 0 réf récente, remplacés par v3.2.3 (≈ 22 G)

| Taille | Date | Dossier | Rôle | Reco |
|---|---|---|---|---|
| 4.9 G | 04-20 | `H2_uniform_vs_real` | ancien test H2 (uniform vs real LAD), 0 réf récente | 🟡 |
| 3.5 G | 06-08 | `musica_hobo_z1raw` | z=1 brut, 0 réf (remplacé par z05) | 🟡 |
| 3.1 G | 06-07 | `Chapter3_d20` | Ch3 v3.2.0 alt-config d20, **0 réf** | 🟡 |
| 3.1 G | 06-05 | `Chapter3_per_site` | Ch3 v3.2.0 per-site, **0 réf** | 🟡 |
| 3.1 G | 06-05 | `Chapter3_common` | Ch3 v3.2.0 common-d_opt, **0 réf** | 🟡 |
| 1.7 G | 04-01 | `musica_option_b_ERA5_results_0104_20mf` | option_b ancien, 0 réf récente | 🟡 |
| 1.4 G | 03-26 | `musica_option_b_ERA5_results` | option_b ancien | 🟡 |
| 888 M | 03-26 | `musica_option_b_Safran_results` | option_b Safran ancien | 🟡 |

---

## 🟠 REVIEW — vieux runs / ton jugement requis (≈ 33 G)

| Taille | Date | Dossier | Rôle | Reco |
|---|---|---|---|---|
| 12 G | 01-27 | `TS` | time-series (janv. 2026, très ancien) ; 5 réf, 0 récente | 🟠 |
| 3.5 G | 06-05 | `musica_hobo_v10_fcovmean` | run versionné (fCover mean test) | 🟠 |
| 3.5 G | 06-02 | `musica_hobo_v11_model323` | run versionné v11 | 🟠 |
| 3.5 G | 06-15 | `musica_hobo_z05` | factoriel v3.2.0 (remplacé par z05_v323 ?) | 🟠 |
| 2.6 G | 05-19 | `musica_hobo_validation` | ancienne validation | 🟠 |
| 2.4 G | 05-31 | `musica_hobo_v9` | run versionné v9 | 🟠 |
| 2.3 G | 07-06 | `musica_hobo_z05_v323station_iter` | factoriel v3.2.3 station-iter (récent — peut servir) | 🟠 |
| 2.2 G | 06-15 | `radius_test_square` | test rayon carré (radius actif via `musica_hobo_radius`) | 🟠 |
| 2.2 G | 06-15 | `radius_test_circle` | test rayon cercle | 🟠 |
| 1.5 G | 02-03 | `radius_test` | test rayon (fév., ancien) | 🟠 |
| 1.1 G | 05-31 | `musica_hobo_v2` | run versionné v2 | 🟠 |
| 1.1 G | 05-31 | `musica_hobo_hmaxfix` | test Hmax fix | 🟠 |
| 874 M | 04-01 | `Annexe_S2_FORMS-H` | annexe FORMS-H (tiret ≠ `Annex_S2_FORMSH`) | 🟠 |
| 439 M | 05-17 | `version_comparison_hobo` | comparaison v3.2.0/v3.2.3 | 🟠 |
| 436 M | 05-12 | `version_comparison_hobo_fh30` | idem fh30 | 🟠 |
| 405 M | 2025-10 | `musica` | très ancien | 🟠 |
| ~1.3 G | 05/06 | `H1_archetypes_v2/v9/v10/v11/hmaxfix` (5×265 M) | archétypes versionnés anciens | 🟠 |
| ~650 M | 03-04 | `musica_hobo*_ERA5_results*` (3×219 M) | anciens ERA5 | 🟠 |
| ~1.2 G | 07-02/06 | `musica_hobo_{2yr_base,2yr_v2,REAL_v2,frblo_v2}` + `z05_station_*` (≈8×140-290 M) | tests forçage station récents (07-02/06) — voir si encore utiles | 🟠 |
| ~200 M | 05 | `musica_hobo_v5/v6/v7*`, `test_*`, `nmlpatch`, `before_Klara`, `modifs_Klara`, `tmp` | petits tests méthodo | 🟠 |

---

## 🔒 KEEP — v3.2.3 courant ou référencé récemment

| Taille | Date | Dossier | Rôle |
|---|---|---|---|
| 3.4 G | 07-07 | `Chapter3/nc` | **cache A canonique Ch3 v3.2.3 (2021)** — extraction en cours |
| 3.2 G | 06-22 | `Chapter1/nc_v323iter` | **cache B — §3.4/§3.5 v3.2.3 (_NM, timing, FUSION_H)** |
| 14 G | 07-04 | `Chapter1/nc_sensitivity_perplot_era5sd` | Shapley perplot (réf. récente) |
| 13 G | 06-20 | `Chapter1/nc_sensitivity_perplot_v323` | Shapley perplot v3.2.3 |
| 13 G | 06-30 | `Chapter1/nc_sensitivity_perplot_fpc` | Shapley perplot FPC (réf. récente) |
| 9.7 G | 06-19 | `Chapter1/nc_sensitivity_perplot` | Shapley perplot base (réf.) |
| 8.6 G | 07-06 | `Chapter1/nc_sensitivity_perplot_frblo` | Shapley perplot frblo (réf. récente) |
| 7.8 G | 07-03 | `Sensitivity_Analysis` | échantillons cLHS + analyses (89 réf. récentes) |
| 3.5 G | 06-09 | `musica_hobo_z05_v323` | factoriel v3.2.3 (9 réf. récentes) |
| 3.3 G | 06-22 | `Chapter1/nc_shapley2x` | Shapley 2× courant |
| 1.7 G | 04-22 | `Annex_S2_FORMSH` | cLHS Shapley 400 |
| 1018 M | 07-07 | `musica_hobo_radius` | sweep rayon actif (script 07-07) |
| 55–331 M | récents | `Chapter1/nc_sensitivity`, `nc_sensitivity_v323`, `nc_h2_v323`, `nc_archetypes`, `tables`, `figures` | sorties Ch1 courantes |

---

## Récap gains potentiels
- 🟢 seul : **≈ 62 G** (sûr) → passe de 53 G à ~115 G libres.
- 🟢 + 🟡 : **≈ 84 G**.
- 🟢 + 🟡 + 🟠 (tri au cas par cas) : jusqu'à **~117 G**.
