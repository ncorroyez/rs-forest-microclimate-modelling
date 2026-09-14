# Comparaison v3.2.0 vs v3.2.3 — extension microclimat fine

*Extension de la comparaison binaire existante (cf. `RECAP_travaux.md` §3,
`review/comparison_v320_v323iter.md`, et les `FigCmp_*` / `tab_*` de ce dossier).*

Le travail antérieur concluait, sur la validation HOBO ΔTmax à 1 m :
**v3.2.0 valide nettement mieux (r 0.93 → 0.63), v3.2.3 surpondère le LAD →
garder le legacy v3.2.0, yoyo écarté.** Cette extension reprend la même
opposition mais la **décompose finement** (métrique de Gril slope-and-equilibrium,
ΔTmax/ΔTmin/DTR, jour/nuit/global, JJAS/10 %-chauds, classement vs amplitude,
profils verticaux T/vent/RH/VPD) et **nuance le verdict** : ce n'est pas un
« v3.2.0 meilleur » en bloc, c'est un **arbitrage amplitude-diurne (v3.2.0) vs
classement-robuste + nocturne (v3.2.3)**.

Données communes : `STATIC_ALS` (LiDAR plein) des 53 plots HOBO de Blois, été
2021 ; macro = forçage free-air ; sims v3.2.0 = `out_files/Chapter3/nc/`,
v3.2.3-iter (ABL/yoyo, PBLH MERRA-2) = `out_files/Chapter3/nc_v323iter/`.

---

## Fiche 5-champs

**Méthodes.** (1) Métrique de Gril : `lmer(T_micro ~ T_macro + (1|Month))`,
slope = effet fixe (β<1 buffering), equilibrium = intercept/(1−slope). (2) ΔTmax,
ΔTmin, DTR, offset ΔT = mean(T_micro−T_macro). (3) Décomposition jour (10–15 h) /
nuit (22–03 h) / global ; période JJAS / 10 %-chauds (Tmax macro ≥ p90).
(4) Classement = Pearson r + Spearman ρ + IC bootstrap + test de Williams
(corrélations dépendantes) ; amplitude = SD(pred)/SD(obs) ; + biais, RMSE.
(5) Comparaison mécaniste (MuSICA forward) vs statistique (régression de Gril en
LOO sur canopée ± topo) ; attribution + partition de variance. (6) Profils
verticaux sur les 15 niveaux nc (RH/VPD via Tetens, P=101325 Pa).

**Pré-traitements.** HOBO position a, plots `ids_to_remove` exclus ; horaire
floor à l'heure ; topo extraite à 10 m (`03_RESULTS/.../Deciduous_Only/`) ;
seed phéno 2020/366 pour l'iter. 53/53 plots valides sur toutes les métriques.

**Données / figures.** voir tableau ci-dessous (script → figure → tables).

**Description.** Grille complète de comparaison des deux binaires sur le
microclimat de sous-bois et sa structure verticale, archétypes P1 (ouvert) → P4
(dense).

**Interprétation (verdict affiné).** voir « Bilan » plus bas.

---

## Scripts → figures → résultat

Scripts dans `scripts/` (à lancer **depuis la racine `z_Example_rmusica_31012025/`** :
`Rscript Chapitre1/comparaison_versions/scripts/02_analyses_microclimat/<x>.R`).
`01_generate_sims/` génère les nc v3.2.3 ; `02_analyses_microclimat/` produit
figures/tables.

| Script | Figure | Résultat clé |
|---|---|---|
| `c1_clhs_trait_dispersion` (racine) | `Fig 5` = `image_dispersion.png` (← `out_files/Chapter3/figures/FigAnnex_clhs_trait_dispersion.png`) | mise en place : dispersion des traits par archétype ; SD Hmax max en P1 (4–40 m), découplage hauteur–couvert (r=0.07) ; bornes cLHS = marginales empiriques |
| `cmp_v320_v323_direct` | `Fig 6` = `image_directcmp.png` (+ `tab_v320_v323_direct.csv`) | **comparaison directe modèle-à-modèle** (sans HOBO) : le yoyo amortit le cycle diurne — ΔTmax_chaud −0.48 °C, DTR −0.73 °C, ΔTmin +0.39 °C, biais inchangé (p=0.44) ; rang corrélé r 0.54–0.84 |
| `cmp_v320_v323_ampbuff` | `Fig 16` = `image_ampbuff.png` (+ `tab_v320_v323_ampbuff.csv`) | **split amplificateurs (ΔTmax>0, n=9) / tampons (n=44)** vs obs : chute Pearson v3.2.3 (0.92→0.60) = entièrement les 9 trouées (aplaties vers 0) ; sur les 44 tampons (refuges) v3.2.3 > v3.2.0 sur rang (ρ 0.79→0.84), Pearson (0.74→0.81) ET amplitude (41→51 %) ; Spearman global 0.86→0.90 |
| `cmp_v320_v323_complements` | `Fig 22-25` = `image_{cmp_byarchetype,gradient,taylor,divergence}.png` (+ `tab_v320_v323_robustness.csv`) | 5 compléments : #1 Δ par archétype = effet jour(ouvert P1)/nuit(dense P4) contrasté ; #2 v3.2.0 préserve la magnitude du gradient ΔTmax P1→P4, v3.2.3 l'aplatit ; #3 sans les 9 trouées, Pearson v3.2.3 (0.81) > v3.2.0 (0.74) → verdict s'inverse ; #4 Taylor (les 2 sous-dispersent ; v3.2.3 gagne l'amplitude nocturne) ; #5 divergence prédite (R²=0.81) par VCI<0/fCover<0/Hmax>0 = trouées ouvertes étagées |
| `plot_ts_lai_temp_archetypes` | `Fig_ts_lai_{temp,dTmax}_archetypes`, `Fig_monthly_lai_dTmax_archetypes` | séries LAI + temp/ΔTmax par scénario × 4 HOBO ; sims tassent le gradient |
| `plot_slope_evolution_archetypes` | `Fig_slope_evolution_archetypes` | pente β mensuelle ; sims épinglées ~1, obs déploie 0.5–1.5 |
| `rf_topo_lai_phase1` | — (`rf_topo_lai_*csv`) | topo n'améliore PAS la prédiction du LAI (ΔR²_LOO = −0.03) |
| `rf_topo_microclimate_direct` | `Fig_microclim_Q1_gradient`, `Fig_microclim_Q2_topo_residual` | empirique recale le gradient ; topo = northness sur ΔTmax seulement |
| `gril_mechanistic_vs_statistical` | `Fig_gril_mechanistic_vs_statistical_slope` | slope : v3.2.0 r=0.96, v3.2.3 r=0.79, Gril LOO R²=0.87 ; équilibre : v3.2.0 inverse (ρ −0.81) |
| `gril_recap_advanced_JJAS` | `Fig_recap_JJAS_advanced` | IC bootstrap + Williams (v3.2.0>v3.2.3 p<1e-4 ; mécaniste≈statistique en rang) ; buffering = canopée (84 % variance unique) |
| `fig_v323_amplitude_vs_ranking` | `Fig_v323_amplitude_vs_ranking` | slope : amplitude identique (31 %), v3.2.3 range mieux ; ΔTmax : v3.2.0 amplitude (44 vs 20 %) |
| `compare_binaries_extra_axes` | `Fig_compare_binaries_extra_axes` | nuit : v3.2.3 gagne (slope SDrec 53→113 %) ; ΔTmin froid : r 0.47→0.82 |
| `compare_binaries_full_grid` | `Fig_compare_full_grid` | grille rang/amplitude × jour-nuit-global × JJAS/chaud : v3.2.3 ρ ≥ partout ; v3.2.0 inverse (équil, jour-chaud) |
| `compare_vertical_profiles[_full]` | `Fig_vertical_profiles_{midday,full}` | profil ouvert P1 jour : sous-bois +3 °C/VPD +0.35 (v3.2.0) aplati par v3.2.3 ; inversion nocturne + forte en v3.2.3 ; cisaillement vent + raide en v3.2.3 |
| `cmp_v320_v323iter_ch3[_rmse]` | `Fig_v320_v323iter_{monthly_slopeR2,rmse_grid}` | pente R² mensuelle + RMSE (antérieur, 13 scénarios) |

**Figures numérotées (jeu livrable).** Les 25 figures du document, nommées par leur numéro narratif `Fig1_…` à `Fig25_…`, sont dans `figures/numbered/` (le manuscrit `Comparaison_v320_v323_RESUME.md` les référence depuis ce dossier). Les sources brutes restent dans `figures/` (`Fig_*.png`) et `figures/doc_media/` (`image*.png`).

Tables associées : `*.csv` dans `tables/` (préfixes `recap_JJAS_`, `compare_*`,
`gril_`, `vertical_profiles_`, `microclim_direct_`, `rf_topo_`, `ts_archetypes_`,
`slope_evolution_`, `Table_v320_v323iter_ch3*`).

---

## Bilan — arbitrage, pas « meilleur/pire »

- **Classement (rang) inter-plots** : v3.2.3 **≥ v3.2.0 presque partout** en
  Spearman (slope 0.92 vs 0.90 ; ΔTmax été 0.90 vs 0.86). La « chute » de Pearson
  de v3.2.3 est **pilotée par 2–3 plots ouverts extrêmes**, pas par l'ordre.
- **Amplitude** : se sépare par moment. **Jour** → v3.2.0 récupère plus
  (ΔTmax 44 % vs 20 %). **Nuit** → v3.2.3 récupère plus (slope nuit 113 % vs
  53 % ; ΔTmin froid r 0.82 vs 0.47). Sur la **slope** l'amplitude est identique
  (31 %), donc v3.2.3 y domine (même amplitude, meilleur rang).
- **Équilibre** : v3.2.0 **inverse** le classement (ρ −0.81, artefact du pôle
  slope→1 sur les plots amplifiants) ; v3.2.3 n'inverse pas (ρ ≈ 0).
- **Profils verticaux** : v3.2.3 **amortit le chauffage convectif du sous-bois
  ouvert le jour** (sa faiblesse Tmax) et **renforce l'inversion/découplage
  nocturne + le cisaillement du vent** (sa force) — signature physique cohérente
  du couplage ABL.
- **Mécaniste vs statistique** : MuSICA (v3.2.0) **classe aussi bien** que la
  régression de Gril (Williams ns) ; le statistique ne gagne que la *calibration*
  d'amplitude (R² 0.87 vs 0.37). Le buffering est un phénomène de **canopée**
  (LAI + fCover, 84 % de variance unique) ; la **topo** (northness/insolation)
  n'apporte qu'une correction *secondaire et métrique-spécifique* (ΔTmax, pas la
  slope), et **rien** à la prédiction du LAI.

**Conséquence pour le choix de binaire.** La décision §3 (garder v3.2.0) reste
défendable **si la cible est l'amplitude des extrêmes chauds diurnes**. Si la
cible est la **carte ordonnée des refuges** (rang, transférabilité) ou le
**régime nocturne**, **v3.2.3 est au moins aussi bon, voire meilleur**. À
trancher selon la question thermique mise en avant — ce n'est plus un rejet en
bloc de v3.2.3.
