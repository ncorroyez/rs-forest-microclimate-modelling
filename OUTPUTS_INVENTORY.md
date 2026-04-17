# Inventaire des outputs — Chapter1_refactored.R

*Généré le 2026-04-17 — commit court : `ed9ab37`*
*Branche : `docs/outputs-inventory`*

---

## Statistiques globales

| Métrique | Valeur |
|---|---|
| Fichiers totaux | 45 |
| PNG | 30 |
| CSV | 12 |
| TXT | 3 |
| Taille totale | ~4.2 Mo |
| Dernier mtime | 2026-04-17 10:26:52 (s2_annex/s2_metrics.txt) |

---

## Arborescence

```
outputs/
├── audit/       ( 9 fichiers, ~532 Ko)
│   ├── concurvity_full.csv
│   ├── concurvity_pairwise_worst.csv
│   ├── h2_amplitude_histogram.png
│   ├── h2_boxplot_archetype.png          ← vestige branche précédente
│   ├── h2_boxplot_cluster.png
│   ├── h2_diff_correlations.csv
│   ├── h2_scatter_metrics.png
│   ├── h2_structure_by_archetype.csv     ← vestige branche précédente
│   └── h2_structure_by_cluster.csv
├── clusters/    ( 2 fichiers, ~156 Ko)
│   ├── cluster_mean_profiles.png
│   └── kmeans_elbow.png
├── fpca/        ( 6 fichiers, ~1.2 Mo)
│   ├── fpca_gcv_curve.png
│   ├── fpca_harmonics.png
│   ├── fpca_loadings.png
│   ├── fpca_reconstruction.png
│   ├── fpca_scatters.png
│   └── fpca_variance.txt
├── gamm/        ( 4 fichiers, ~200 Ko)
│   ├── gamm_marginal_effects.png
│   ├── gamm_residuals_diagnostic.png
│   ├── gamm_residuals_stats.csv
│   └── gamm_summary.txt
├── h1/          ( 3 fichiers, ~140 Ko)
│   ├── h1_forward_curve.png
│   ├── h1_hierarchy.png
│   └── h1_scores.csv
├── h2/          (14 fichiers, ~1.3 Mo)
│   ├── c1_median_lad_profile.png
│   ├── c2_median_lad_profile.png
│   ├── c3_median_lad_profile.png
│   ├── h2_distributions.png
│   ├── h2_heatwave_histogram.png
│   ├── h2_paired_real_vs_uniform.png
│   ├── h2_verify_cluster1_peak.csv
│   ├── h2_verify_cluster1_peak.png
│   ├── h2_vertical_tmax_profiles_3way_canicule.png
│   ├── h2_vertical_tmax_profiles_canicule_date.csv
│   ├── h2_vertical_tmax_profiles_canicule_date.png
│   ├── h2_vertical_tmax_profiles_mean_summer.csv
│   ├── h2_vertical_tmax_profiles_mean_summer.png
│   ├── h2_vertical_tmax_profiles_median_date.csv
│   └── h2_vertical_tmax_profiles_median_date.png
├── hobo/        ( 3 fichiers, ~536 Ko)
│   ├── hobo_metrics.csv
│   ├── hobo_timeseries_representatives.png
│   └── hobo_validation_rmse.png
└── s2_annex/    ( 3 fichiers, ~312 Ko)
    ├── s2_density.png
    ├── s2_metrics.txt
    └── s2_paired.png
```

---

## Par sous-dossier

### `outputs/audit/`

#### `audit/concurvity_full.csv`

- **Colonnes** : (index), para, s(LAI_sc), s(Hmax_sc), s(fCover_sc), s(date_factor), s(plot_id) — 3 lignes (worst / observed / estimate)
- **Head** :

| | para | s(LAI_sc) | s(Hmax_sc) | s(fCover_sc) | s(date_factor) | s(plot_id) |
|---|---|---|---|---|---|---|
| worst | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| observed | 1.000 | 1.000 | 1.000 | 1.000 | 0.005 | 0.031 |
| estimate | 1.000 | 1.000 | 1.000 | 1.000 | 0.012 | 0.100 |

- **Producteur R** : `run_concurvity_audit()` ligne 1198
- **Description** : Audit de concurvité global (full=TRUE) — worst-case par terme smooth. Les smooth structurels (LAI, Hmax, fCover) montrent une concurvité maximale entre eux (=1.0 worst), mais s(date_factor) et s(plot_id) ont des valeurs observed très faibles (0.005 et 0.031), confirmant leur indépendance par rapport aux termes structurels.
- **Usage recommandé** : Diagnostic interne

---

#### `audit/concurvity_pairwise_worst.csv`

- **Colonnes** : (index), para, s(LAI_sc), s(Hmax_sc), s(fCover_sc), s(date_factor), s(plot_id) — 6 lignes (matrice pairwise)
- **Head** :

| | para | s(LAI_sc) | s(Hmax_sc) | s(fCover_sc) | s(date_factor) | s(plot_id) |
|---|---|---|---|---|---|---|
| para | 1.000 | 0.010 | 0.012 | 0.318 | 1.000 | 1.000 |
| s(LAI_sc) | 0.010 | 1.000 | 0.440 | 0.792 | 0.010 | 1.000 |
| s(Hmax_sc) | 0.012 | 0.440 | 1.000 | 0.200 | 0.012 | 1.000 |
| s(fCover_sc) | 0.318 | 0.792 | 0.200 | 1.000 | 0.318 | 1.000 |

- **Producteur R** : `run_concurvity_audit()` ligne 1200
- **Description** : Matrice pairwise worst — révèle la concurvité élevée LAI ↔ fCover (0.792). Tous les termes montrent une concurvité maximale avec s(plot_id), attendue (plot_id capte la variation intra-site).
- **Usage recommandé** : Diagnostic interne

---

#### `audit/h2_amplitude_histogram.png`

- **Objet** : Histogramme de l'amplitude temporelle de l'effet LAD par plot (q95 - q05 de la différence ΔTmax Real - Uniform), coloré par Cluster. Médiane globale et % de plots avec amplitude > 1°C annotés.
- **Dimensions** : 1500 × 900 px
- **Taille** : 51 Ko
- **Producteur R** : `analyse_h2_structure()` ligne 847 — `ggsave(file.path(out_dir, "h2_amplitude_histogram.png"), p4, width=10, height=6, dpi=150)`
- **Usage recommandé** : Annexe / diagnostic

---

#### `audit/h2_boxplot_archetype.png`

- **Objet** : Identique à `h2_boxplot_cluster.png` — boxplot + jitter de la différence moyenne par plot (Real - Uniform ΔTmax) groupé par archetype/cluster. Vestige d'une ancienne nomenclature (`archetype` → `cluster` lors du refactor `refactor/rename-archetype-to-cluster`).
- **Dimensions** : 1500 × 900 px
- **Taille** : 101 Ko
- **Producteur R** : `analyse_h2_structure()` (version antérieure) — produit par `ggsave` avec `out_dir` pointant vers `outputs/audit/`
- **Usage recommandé** : **Obsolète — peut être supprimé** (doublon de `h2_boxplot_cluster.png`). # TODO: à clarifier avec Nathan

---

#### `audit/h2_boxplot_cluster.png`

- **Objet** : Boxplot + jitter de la différence par plot (mean ΔTmax Real - Uniform LAD) groupé par Cluster (3 clusters). Sous-titre : médiane par cluster et taille d'échantillon. Titre : "H2 — Structure de l'effet LAD par cluster (per-plot mean)".
- **Dimensions** : 1500 × 900 px
- **Taille** : 99 Ko
- **Producteur R** : `analyse_h2_structure()` ligne 727 — `ggsave(file.path(out_dir, "h2_boxplot_cluster.png"), p1, width=10, height=6, dpi=150)`
- **Usage recommandé** : Slide principale H2

---

#### `audit/h2_diff_correlations.csv`

- **Colonnes** : metric, pearson, spearman, n — 3 lignes (LAI, Hmax, fCover)
- **Head** :

| metric | pearson | spearman | n |
|---|---|---|---|
| LAI | -0.067 | -0.066 | 300 |
| Hmax | -0.088 | -0.108 | 300 |
| fCover | -0.056 | -0.075 | 300 |

- **Producteur R** : `analyse_h2_structure()` ligne 797
- **Description** : Corrélations Pearson et Spearman entre la différence moyenne per-plot (Real - Uniform ΔTmax) et les métriques structurelles (LAI, Hmax, fCover). Toutes très faibles (|r| < 0.11) — confirme que l'effet LAD n'est pas simplement lié aux métriques 2D classiques.
- **Résumé** : pearson ∈ [-0.088, -0.056] ; spearman ∈ [-0.108, -0.066] ; n = 300 pour toutes.
- **Usage recommandé** : Données brutes pour ré-analyse / narratif H2

---

#### `audit/h2_scatter_metrics.png`

- **Objet** : Nuage de points facetté (un panneau par métrique : LAI, Hmax, fCover) — diff mean ΔTmax (Real - Uniform) en Y vs valeur de la métrique en X, avec régression linéaire et annotation r de Pearson. Titre : "H2 — Diff LAD (real-uniform) vs métriques structurelles par plot".
- **Dimensions** : 1800 × 1200 px
- **Taille** : 248 Ko
- **Producteur R** : `analyse_h2_structure()` ligne 822 — `ggsave(file.path(out_dir, "h2_scatter_metrics.png"), p3, width=12, height=8, dpi=150)`
- **Usage recommandé** : Annexe H2 / diagnostic

---

#### `audit/h2_structure_by_archetype.csv`

- **Colonnes** : Archetype, period, n, median_diff, mean_diff, sd_diff, q05, q95 — 6 lignes (3 archetypes × 2 périodes)
- **Head** (3 premières lignes) :

| Archetype | period | n | median_diff | mean_diff | sd_diff |
|---|---|---|---|---|---|
| 1 | heatwave | 1100 | -0.008 | -0.008 | 0.717 |
| 1 | normal | 11100 | 0.011 | 0.037 | 0.860 |
| 2 | heatwave | 1100 | 0.013 | 0.006 | 0.198 |

- **Producteur R** : `analyse_h2_structure()` (version antérieure avec nomenclature `archetype`) — vestige de la branche `refactor/rename-archetype-to-cluster`.
- **Description** : Identique à `h2_structure_by_cluster.csv` (contenu numérique strictement identique). Vestige terminologique.
- **Usage recommandé** : **Obsolète** — utiliser `h2_structure_by_cluster.csv` # TODO: à clarifier avec Nathan

---

#### `audit/h2_structure_by_cluster.csv`

- **Colonnes** : Cluster, period, n, median_diff, mean_diff, sd_diff, q05, q95 — 6 lignes (3 clusters × 2 périodes)
- **Head** :

| Cluster | period | n | median_diff | mean_diff | sd_diff | q05 | q95 |
|---|---|---|---|---|---|---|---|
| 1 | heatwave | 1100 | -0.008 | -0.008 | 0.717 | -0.400 | 0.366 |
| 1 | normal | 11100 | 0.011 | 0.037 | 0.860 | -0.377 | 0.514 |
| 2 | heatwave | 1100 | 0.013 | 0.006 | 0.198 | -0.227 | 0.220 |

- **Producteur R** : `analyse_h2_structure()` ligne 749
- **Description** : Statistiques descriptives de l'effet LAD (Real - Uniform ΔTmax) stratifiées par Cluster et période macroclimatique (heatwave : Tmax_macro ≥ 30°C / normal). Cluster 2 a la plus petite sd_diff (0.198 vs 0.717 pour C1) — profil en tige uniforme, effet LAD très resserré.
- **Résumé** : n heatwave = 1100 / cluster, n normal = 11100 / cluster ; sd_diff ∈ [0.20, 0.86] ; médiane globale ≈ 0.
- **Usage recommandé** : Données brutes pour ré-analyse / table manuscrit

---

### `outputs/clusters/`

#### `clusters/cluster_mean_profiles.png`

- **Objet** : Profils LAD moyens par cluster (3 courbes) en hauteur relative (z/Hmax), avec annotation pour nommer chaque profil-type (tige uniforme, exponentiel dégressif, bimodal). Titre : "Profil LAD moyen par cluster — base d'interprétation". Dimensions 1800 × 900 px.
- **Dimensions** : 1800 × 900 px
- **Taille** : 99 Ko
- **Producteur R** : `plot_cluster_mean_profiles()` → appelé `main()` ligne 1929 — `save_plot(p_clust_prof, "outputs/clusters/cluster_mean_profiles.png", width=12, height=6)`
- **Usage recommandé** : Slide principale — figure d'introduction aux clusters

---

#### `clusters/kmeans_elbow.png`

- **Objet** : Courbe du coude (Elbow Method) pour la sélection du nombre optimal de clusters K-means sur les métriques LiDAR (LAI, Hmax, fCover). Titre : "Elbow Method for K-means Stratification". Optimal k annoté en sous-titre.
- **Dimensions** : 1200 × 750 px
- **Taille** : 52 Ko
- **Producteur R** : `optimise_k_elbow()` → `main()` ligne 1888 — `save_plot(cl_res$elbow_plot, "outputs/clusters/kmeans_elbow.png", width=8, height=5)`
- **Usage recommandé** : Annexe méthodologique / diagnostic

---

### `outputs/fpca/`

#### `fpca/fpca_gcv_curve.png`

- **Objet** : Courbe GCV (Generalized Cross-Validation) pour la sélection du nombre de fonctions de base B-spline de la FPCA. Titre : "GCV Curve for FPCA Basis Functions". Permet de justifier le nombre de bases retenu.
- **Dimensions** : 1200 × 750 px
- **Taille** : 43 Ko
- **Producteur R** : `optimise_nbasis()` → `main()` ligne 1905 — `save_plot(fpca_res$gcv_curve, "outputs/fpca/fpca_gcv_curve.png", width=8, height=5)`
- **Usage recommandé** : Annexe méthodologique

---

#### `fpca/fpca_harmonics.png`

- **Objet** : Figure multi-panneaux (1 panneau par FPC, nrow=1) montrant l'harmonique ±1 écart-type de chaque composante principale fonctionnelle (FPC1, FPC2, FPC3) avec la courbe moyenne et les ±perturbations. Titre par panneau : "FPC k — XX.X% var". Titre général inféré de `plot_fpc_harmonics()`.
- **Dimensions** : 1800 × 900 px
- **Taille** : 156 Ko
- **Producteur R** : `plot_fpc_harmonics()` → `main()` ligne 1909 — `save_plot(p_fpca_harm, "outputs/fpca/fpca_harmonics.png", width=12, height=6)`
- **Usage recommandé** : Slide principale FPCA — figure clé pour l'interprétation des composantes

---

#### `fpca/fpca_loadings.png`

- **Objet** : Loadings (contributions par hauteur) des 3 premières FPC en hauteur relative (Z/Hmax), un panneau par FPC. Titre : "Per-height contribution of each FPC (FPCA loadings)". Couleurs viridis.
- **Dimensions** : 1500 × 750 px
- **Taille** : 89 Ko
- **Producteur R** : `plot_fpc_loadings()` → `main()` ligne 1913 — `save_plot(p_fpca_load, "outputs/fpca/fpca_loadings.png", width=10, height=5)`
- **Usage recommandé** : Annexe / slide FPCA secondaire

---

#### `fpca/fpca_reconstruction.png`

- **Objet** : Reconstruction FPCA pour 6 sites représentatifs (échantillonnés le long du gradient FPC1) — comparaison profil réel LiDAR vs reconstruit par les 3 premières CP. Sous-titre : "FPC1=41.1% | FPC2=30.9% | FPC3=15.2%". Titre : "Reconstruction FPCA — profils réels vs 3 premières composantes".
- **Dimensions** : 2100 × 1050 px
- **Taille** : 171 Ko
- **Producteur R** : `plot_fpca_reconstruction()` → `main()` ligne 1024 — `save_plot(p, out_path, width=14, height=7)`
- **Usage recommandé** : Slide principale FPCA — valide la qualité de la représentation (87.2% de variance)

---

#### `fpca/fpca_scatters.png`

- **Objet** : Nuages de points FPC1 vs FPC2, FPC1 vs FPC3, FPC2 vs FPC3 (3 panneaux) colorés par Cluster. Permet de voir la séparation des clusters dans l'espace fonctionnel.
- **Dimensions** : 1800 × 750 px
- **Taille** : 658 Ko
- **Producteur R** : `plot_fpc_scatters()` → `main()` ligne 1917 — `save_plot(p_fpca_scat, "outputs/fpca/fpca_scatters.png", width=12, height=5)`
- **Usage recommandé** : Slide principale / annexe H2 — lien FPCA ↔ clustering

---

#### `fpca/fpca_variance.txt`

- **Contenu** : Variance expliquée par les 3 premières FPC. FPC1: 41.1% / FPC2: 30.9% / FPC3: 15.2% (total: 87.2%)
- **Taille** : 36 octets
- **Producteur R** : `main()` ligne 1901-1903 — `writeLines(sprintf(...), "outputs/fpca/fpca_variance.txt")`
- **Usage recommandé** : Diagnostic interne / reporting automatique

---

### `outputs/gamm/`

#### `gamm/gamm_marginal_effects.png`

- **Objet** : Grille de panneaux des effets marginaux du GAMM de référence (via `ggpredict`) — un panneau par terme smooth (LAI_sc, Hmax_sc, fCover_sc, Cluster) + forçages macroclimatiques (VPD, Vent, Rad) si `use_macroclimate=TRUE`. Titre global : "GAMM Marginal Effects". Sous-titre : "Holding all other variables at their mean".
- **Dimensions** : 1800 × 1200 px
- **Taille** : 97 Ko
- **Producteur R** : `plot_gamm_marginal_effects()` → `main()` ligne 2176 — `save_plot(p_gamm_eff, "outputs/gamm/gamm_marginal_effects.png", width=12, height=8)`
- **Usage recommandé** : Slide principale GAMM — figure de résultats centraux

---

#### `gamm/gamm_residuals_diagnostic.png`

- **Objet** : 4 panneaux de diagnostic des résidus gaussiens (via `diagnose_residuals()`) : QQ-plot, histogramme des résidus + courbe normale, résidus vs fitted, comparaison AIC Gaussien vs scat(). Produit via `png()` base R (pas `ggsave`).
- **Dimensions** : 1200 × 900 px
- **Taille** : 86 Ko
- **Producteur R** : `diagnose_residuals()` → `main()` ligne 2179 — `png("outputs/gamm/gamm_residuals_diagnostic.png", width=1200, height=900)`
- **Usage recommandé** : Diagnostic interne — justifie le choix de la famille `scat()` (ΔAIC = 32 163)

---

#### `gamm/gamm_residuals_stats.csv`

- **Colonnes** : n_obs, shapiro_W, shapiro_p, excess_kurt, aic_gauss, aic_scat, delta_aic, scat_preferred — 1 ligne
- **Head** :

| n_obs | shapiro_W | shapiro_p | excess_kurt | aic_gauss | aic_scat | delta_aic | scat_preferred |
|---|---|---|---|---|---|---|---|
| 36600 | 0.689 | 3.0e-70 | 31.95 | 71177 | 39015 | 32163 | TRUE |

- **Producteur R** : `diagnose_residuals()` → `main()` ligne 2183
- **Description** : Métriques de diagnostic résiduel : Shapiro-Wilk W=0.689 (p≈0, forte non-normalité), excès de kurtosis=31.9, ΔAIC=32 163 en faveur de la famille scaled-t (`scat()`).
- **Usage recommandé** : Données brutes pour ré-analyse / table supplémentaire

---

#### `gamm/gamm_summary.txt`

- **Contenu** : Sortie textuelle de `summary(gam_ref)` — GAMM `bam` famille scat(3, 0.249), n=36600, R²(adj)=0.827, Deviance=73.5%. Termes : s(LAI_sc) edf=6.65 ***, s(Hmax_sc) edf=4.97 ***, s(fCover_sc) edf=3.67 ***, s(date_factor) edf=120.84 ***, s(plot_id) edf=254.68 n.s.
- **Taille** : 1 041 octets
- **Producteur R** : `main()` ligne 2170 — `capture.output(summary(gam_ref), file="outputs/gamm/gamm_summary.txt")`
- **Usage recommandé** : Référence manuscrit / annexe statistique

---

### `outputs/h1/`

#### `h1/h1_forward_curve.png`

- **Objet** : Courbe d'inclusion progressive (Forward inclusion) — RMSE vs référence Full_real en Y, scénarios en X (ordre d'inclusion : LAI, Hmax, fCover, Full_real). Titre : "H1 forward — information added incrementally". Montre la contribution marginale de chaque paramètre.
- **Dimensions** : 1500 × 1050 px
- **Taille** : 73 Ko
- **Producteur R** : `plot_forward_curve()` → `main()` ligne 2145 — `save_plot(p_h1_fwd, "outputs/h1/h1_forward_curve.png")`
- **Usage recommandé** : Slide principale H1 — figure centrale pour la hiérarchie des paramètres

---

#### `h1/h1_hierarchy.png`

- **Objet** : Barplot horizontal RMSE par scénario (tous les scénarios H1f et H1l) classés par RMSE croissant. Titre : "H1 — divergence from full-real reference". Montre le score de chaque scénario d'omission/inclusion par rapport au scénario de référence full-real.
- **Dimensions** : 1500 × 1050 px
- **Taille** : 52 Ko
- **Producteur R** : `plot_scenario_hierarchy()` → `main()` ligne 2139 — `save_plot(p_h1_hier, "outputs/h1/h1_hierarchy.png")`
- **Usage recommandé** : Annexe H1 / slide secondaire

---

#### `h1/h1_scores.csv`

- **Colonnes** : scenario, n, mean_diff, mae, rmse — 11 lignes (scénarios H1f et H1l)
- **Head** (5 premières lignes) :

| scenario | n | mean_diff | mae | rmse |
|---|---|---|---|---|
| H1f_4_Full_real | 36600 | 0.000 | 0.000 | 0.000 |
| REF_all_real | 36600 | 0.000 | 0.000 | 0.000 |
| H1l_dropLAD_meanShape | 36600 | 0.031 | 0.243 | 0.664 |
| H1f_3_LAI_Hmax_fCover | 36600 | -0.019 | 0.232 | 0.699 |
| H1l_dropLAD_uniform | 36600 | -0.019 | 0.232 | 0.699 |

- **Producteur R** : `score_scenarios_vs_reference()` → `main()` ligne 2136
- **Description** : Scores de divergence (RMSE vs référence Full_real) pour tous les scénarios H1 forward et LOO. Le scénario Null_baseline (RMSE=1.605°C) vs Full_real (RMSE=0°C) borne la plage. Drop fCover (RMSE=0.868) > Drop LAD (RMSE=0.664) > Drop Hmax (RMSE=0.791).
- **Résumé** : RMSE ∈ [0.000, 1.605] ; MAE ∈ [0.000, 1.188] ; n ∈ [30378, 36600] selon les scénarios.
- **Usage recommandé** : Données brutes pour ré-analyse — table H1 manuscrit

---

### `outputs/h2/`

#### `h2/c1_median_lad_profile.png`

- **Objet** : Profil LAD individuel du plot médian de Cluster 1 (hauteur absolue en m en Y, LAD m²·m⁻³ en X). Hmax annoté par ligne dashed rouge. Sous-titre : LAI, Hmax, fCover du plot. Titre : "Profil LAD — Plot {id} (Cluster 1)".
- **Dimensions** : 900 × 1050 px
- **Taille** : 43 Ko
- **Producteur R** : `plot_lad_profile_for_plot()` → `main()` ligne 1463 — `save_plot(p, out_path, width=6, height=7)`
- **Usage recommandé** : Slide principale H2 — illustration du profil-type de chaque cluster

---

#### `h2/c2_median_lad_profile.png`

- **Objet** : Profil LAD individuel du plot médian de Cluster 2.
- **Dimensions** : 900 × 1050 px
- **Taille** : 56 Ko
- **Producteur R** : `plot_lad_profile_for_plot()` → `main()` ligne 1463
- **Usage recommandé** : Slide principale H2

---

#### `h2/c3_median_lad_profile.png`

- **Objet** : Profil LAD individuel du plot médian de Cluster 3.
- **Dimensions** : 900 × 1050 px
- **Taille** : 58 Ko
- **Producteur R** : `plot_lad_profile_for_plot()` → `main()` ligne 1463
- **Usage recommandé** : Slide principale H2

---

#### `h2/h2_distributions.png`

- **Objet** : Distribution de densité des ΔTmax pour les deux scénarios (Real LAD vs Uniform LAD) superposées. Titre : "H2 — ΔTmax distributions: real vs uniform LAD".
- **Dimensions** : 1500 × 1050 px
- **Taille** : 68 Ko
- **Producteur R** : `plot_h2_distributions()` → `main()` ligne 1965 — `save_plot(p_h2_dist, "outputs/h2/h2_distributions.png")`
- **Usage recommandé** : Slide principale H2 — vue d'ensemble de la distribution de l'effet LAD

---

#### `h2/h2_heatwave_histogram.png`

- **Objet** : Histogramme de la différence ΔTmax (Real - Uniform) stratifié par période macroclimatique (normal vs heatwave, Tmax_macro ≥ 30°C). Titre : "H2 — distribution of (Real - Uniform) ΔTmax". Couleurs bleu/rouge par régime.
- **Dimensions** : 1500 × 1050 px
- **Taille** : 41 Ko
- **Producteur R** : `analyse_h2_distribution()` → `main()` ligne 1974 — `save_plot(p_h2_hw, "outputs/h2/h2_heatwave_histogram.png")`
- **Usage recommandé** : Slide principale H2 — montre la sensibilité de l'effet LAD aux canicules

---

#### `h2/h2_paired_real_vs_uniform.png`

- **Objet** : Comparaison par plot/jour (scatter ou paired plot) des ΔTmax Real LAD vs Uniform LAD. Titre : "H2 — Per-plot, per-day ΔTmax: real vs uniform LAD". Sous-titre : "Paired comparison: impact of vertical profile SHAPE on cooling".
- **Dimensions** : 1500 × 1050 px
- **Taille** : 231 Ko
- **Producteur R** : `plot_h2_paired()` → `main()` ligne 1969 — `save_plot(p_h2_pair, "outputs/h2/h2_paired_real_vs_uniform.png")`
- **Usage recommandé** : Slide principale H2 — figure centrale de l'hypothèse H2

---

#### `h2/h2_verify_cluster1_peak.csv`

- **Colonnes** : nair, Tmax_z, height_m, scenario, Cluster, plot_id, LAI, Hmax, panel_label — 180 lignes (3 plots × ~30 couches air × 2 scénarios)
- **Head** (3 premières lignes) :

| nair | Tmax_z | height_m | scenario | Cluster | plot_id |
|---|---|---|---|---|---|
| 1 | 33.163 | 0.499 | Real LAD | 1 | X370920_Y5270300 |
| 2 | 33.175 | 1.556 | Real LAD | 1 | X370920_Y5270300 |
| 3 | 33.199 | 2.740 | Real LAD | 1 | X370920_Y5270300 |

- **Producteur R** : `verify_cluster1_peak()` ligne 1595
- **Description** : Profils verticaux de Tmax (Real vs Uniform LAD, mode canicule_date) pour 3 plots médians de Cluster 1, destinés à vérifier la robustesse du pic de température proche du sol observé dans ce cluster. LAI et Hmax inclus.
- **Résumé** : Tmax_z ∈ [~30°C, ~35°C] selon les hauteurs ; height_m ∈ [0.5 m, ~Hmax].
- **Usage recommandé** : Données brutes pour ré-analyse / vérification robustesse

---

#### `h2/h2_verify_cluster1_peak.png`

- **Objet** : Figure facettée (1 panneau par plot Cluster 1) : profils verticaux Tmax (hauteur en Y, Tmax_z en X), Real LAD vs Uniform LAD comparés, ligne Hmax annotée. Titre : "Vérification pic sol — Cluster 1 (canicule_date)". Permet de confirmer que le pic de température proche du sol en Cluster 1 est robuste sur plusieurs plots.
- **Dimensions** : 2100 × 900 px
- **Taille** : 136 Ko
- **Producteur R** : `verify_cluster1_peak()` ligne 1594 — `save_plot(p, png_path, width=14, height=6)`
- **Usage recommandé** : Annexe H2 / diagnostic robustesse Cluster 1

---

#### `h2/h2_vertical_tmax_profiles_3way_canicule.png`

- **Objet** : Profils verticaux Tmax multi-scénarios (3 scénarios : Real LAD, Uniform LAD, Cluster-type LAD) pour le plot médian de chaque cluster, mode canicule_date. Titre dynamique : "Profils verticaux Tmax — comparaison 3 scénarios". Sous-titre : mode canicule_date.
- **Dimensions** : 2100 × 1050 px
- **Taille** : 128 Ko
- **Producteur R** : `plot_vertical_tmax_profiles_multi()` → `main()` ligne 1866 — `save_plot(p, out_path, width=14, height=7)` avec `out_path="outputs/h2/h2_vertical_tmax_profiles_3way_canicule.png"`
- **Usage recommandé** : Slide principale H2 — figure comparative des 3 scénarios, clé pour H2+H3

---

#### `h2/h2_vertical_tmax_profiles_canicule_date.csv`

- **Colonnes** : nair, Tmax_z, height_m, scenario, Cluster, plot_id, LAI, Hmax, panel_label — 180 lignes
- **Head** (3 premières lignes) :

| nair | Tmax_z | height_m | scenario | Cluster | plot_id |
|---|---|---|---|---|---|
| 1 | 33.163 | 0.499 | Real LAD | 1 | X370920_Y5270300 |
| 2 | 33.175 | 1.556 | Real LAD | 1 | X370920_Y5270300 |
| 3 | 33.199 | 2.740 | Real LAD | 1 | X370920_Y5270300 |

- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="canicule_date") → `main()` ligne 1727
- **Description** : Profils verticaux de Tmax pour les plots médians de chaque cluster (mode : canicule_date — moyennée sur les jours de forte chaleur). 3 clusters × 2 scénarios (Real/Uniform) × ~30 couches.
- **Usage recommandé** : Données brutes pour ré-analyse

---

#### `h2/h2_vertical_tmax_profiles_canicule_date.png`

- **Objet** : Profils verticaux Tmax pour chaque cluster (3 panneaux), Real LAD vs Uniform LAD, mode canicule_date. Titre : "Profils verticaux de Tmax — Real LAD vs Uniform LAD". Sous-titre : mode canicule_date.
- **Dimensions** : 2100 × 900 px
- **Taille** : 126 Ko
- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="canicule_date") → `main()` ligne 1726 — `save_plot(p, png_path, width=14, height=6)`
- **Usage recommandé** : Slide principale H2 — profils pendant canicule

---

#### `h2/h2_vertical_tmax_profiles_mean_summer.csv`

- **Colonnes** : nair, Tmax_z, height_m, scenario, Cluster, plot_id, LAI, Hmax, panel_label — 180 lignes
- **Head** :

| nair | Tmax_z | height_m | scenario | Cluster | plot_id |
|---|---|---|---|---|---|
| 1 | 25.382 | 0.556 | Real LAD | 1 | X367840_Y5266380 |
| 2 | 25.031 | 1.733 | Real LAD | 1 | X367840_Y5266380 |
| 3 | 24.930 | 3.053 | Real LAD | 1 | X367840_Y5266380 |

- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="mean_summer") → `main()` ligne 1727
- **Description** : Profils verticaux Tmax moyennés sur tout l'été. Tmax_z ∈ [~22°C, ~27°C], moins chaud que canicule_date.
- **Usage recommandé** : Données brutes pour ré-analyse

---

#### `h2/h2_vertical_tmax_profiles_mean_summer.png`

- **Objet** : Profils verticaux Tmax Real vs Uniform, mode mean_summer (moyenne sur tous les jours d'été JJA). Titre : "Profils verticaux de Tmax — Real LAD vs Uniform LAD". Sous-titre : mode mean_summer.
- **Dimensions** : 2100 × 900 px
- **Taille** : 125 Ko
- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="mean_summer") → `main()` ligne 1726
- **Usage recommandé** : Annexe H2 — référence "conditions moyennes" vs canicule

---

#### `h2/h2_vertical_tmax_profiles_median_date.csv`

- **Colonnes** : nair, Tmax_z, height_m, scenario, Cluster, plot_id, LAI, Hmax, panel_label — 180 lignes
- **Head** : Identique en structure à canicule_date.csv (même colonnes, valeurs de Tmax_z intermédiaires)
- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="median_date") → `main()` ligne 1727
- **Description** : Profils verticaux pour la date médiane de l'été.
- **Usage recommandé** : Données brutes pour ré-analyse

---

#### `h2/h2_vertical_tmax_profiles_median_date.png`

- **Objet** : Profils verticaux Tmax Real vs Uniform, mode median_date (date médiane de l'été). Titre : "Profils verticaux de Tmax — Real LAD vs Uniform LAD". Sous-titre : mode median_date.
- **Dimensions** : 2100 × 900 px
- **Taille** : 133 Ko
- **Producteur R** : `plot_vertical_tmax_profiles()` (mode="median_date") → `main()` ligne 1726
- **Usage recommandé** : Annexe H2 — condition médiane pour comparaison

---

### `outputs/hobo/`

#### `hobo/hobo_metrics.csv`

- **Colonnes** : scenario, n, r2, rmse, mae, bias — 6 lignes (scénarios H1f validés aux HOBO)
- **Head** (4 premières lignes) :

| scenario | n | r2 | rmse | mae | bias |
|---|---|---|---|---|---|
| H1f_3_LAI_Hmax_fCover | 6337 | 0.435 | 2.628 | 2.247 | 2.089 |
| H1f_4_Full_real | 6459 | 0.482 | 2.637 | 2.252 | 2.081 |
| REF_all_real | 6459 | 0.482 | 2.637 | 2.252 | 2.081 |
| H1f_1_LAI_only | 4995 | 0.388 | 2.742 | 2.393 | 2.166 |

- **Producteur R** : `validate_scenarios_at_hobos()` → `main()` ligne 2199
- **Description** : Métriques de validation des scénarios H1 forward vs capteurs HOBO terrain (ΔTmax observé). Biais systématique de +2.08°C constant sur tous les scénarios — confirme le biais résiduel MuSICA non-lié à la structure forestière. R² ∈ [0.004, 0.482]. Le scénario Null a r²=0.004.
- **Résumé** : RMSE ∈ [2.628, 2.942] ; MAE ∈ [2.247, 2.561] ; Bias ∈ [1.934, 2.166] ; r² ∈ [0.004, 0.482].
- **Usage recommandé** : Données brutes pour ré-analyse — table de validation manuscrit

---

#### `hobo/hobo_timeseries_representatives.png`

- **Objet** : Séries temporelles estivales pour les HOBO représentatifs (choisis par leur capacité de tamponnage médiane) — comparaison HOBO observé / MuSICA Real / MuSICA Uniform / Macroclimate ERA5. Titre : "Time-series — Selected HOBOs vs MuSICA". Figure facettée (un panneau par HOBO).
- **Dimensions** : 1800 × 1500 px
- **Taille** : 475 Ko
- **Producteur R** : `plot_timeseries_faceted()` → `main()` ligne 2212 — `save_plot(p_ts, "outputs/hobo/hobo_timeseries_representatives.png", width=12, height=10)`
- **Usage recommandé** : Slide principale validation HOBO — figure d'illustration qualitative du biais

---

#### `hobo/hobo_validation_rmse.png`

- **Objet** : Barplot horizontal RMSE par scénario H1 forward (axes : scénario reordonné par RMSE, RMSE vs HOBO en X) avec couleur de remplissage = biais moyen (gradient bleu→blanc→rouge). Titre : "HOBO validation — RMSE per scenario". Annotation r² par barre.
- **Dimensions** : 1500 × 1050 px
- **Taille** : 49 Ko
- **Producteur R** : `plot_hobo_validation()` → `main()` ligne 2202 — `save_plot(p_hobo_val, "outputs/hobo/hobo_validation_rmse.png")`
- **Usage recommandé** : Slide principale validation — diagnostic rapide des scénarios vs terrain

---

### `outputs/s2_annex/`

#### `s2_annex/s2_density.png`

- **Objet** : Densités de ΔTmax pour LiDAR Full 3D Real vs Sentinel-2 + FORMS-H (2D Proxy) superposées. Titre : "Annex H3 — Optical Proxy vs Structural Reality (Density)". Sous-titre avec R², RMSE, Bias calculés.
- **Dimensions** : 1500 × 1050 px
- **Taille** : 85 Ko
- **Producteur R** : `main()` ligne 2245 — `save_plot(p_s2_density, "outputs/s2_annex/s2_density.png")`
- **Usage recommandé** : Annexe H3 — figure d'introduction à "l'illusion optique"

---

#### `s2_annex/s2_metrics.txt`

- **Contenu** : R2: 0.3343 / RMSE: 1.5544 deg C / Bias: -0.8396 deg C
- **Taille** : 50 octets
- **Producteur R** : `main()` ligne 2263-2264 — `writeLines(sprintf(...), "outputs/s2_annex/s2_metrics.txt")`
- **Description** : Métriques de l'annexe Sentinel-2 (S2 + FORMS-H proxy vs LiDAR full 3D). Biais négatif de -0.84°C (S2 sous-estime ΔTmax), R²=0.33 — confirme "l'illusion optique" : la proxy 2D capte 33% de la variance spatiale du tamponnage.
- **Usage recommandé** : Reporting / narratif H3

---

#### `s2_annex/s2_paired.png`

- **Objet** : Scatter hexbin per-plot/per-day : LiDAR ΔTmax (X) vs Sentinel-2+FORMS-H ΔTmax (Y), droite y=x dashed, régression lm orange. Annotation R²/RMSE/Bias. Titre : "Annex H3 — Per-plot, per-day ΔTmax: S2 Proxy vs LiDAR". Sous-titre : "Paired comparison: impact of 2D optical proxy on spatial prediction".
- **Dimensions** : 1500 × 1050 px
- **Taille** : 216 Ko
- **Producteur R** : `main()` ligne 2260 — `save_plot(p_s2_paired, "outputs/s2_annex/s2_paired.png")`
- **Usage recommandé** : Annexe H3 — figure principale de l'annexe Sentinel-2

---

## Par catégorie d'usage

### Figures pour slides principales

| Fichier | Section | Description courte |
|---|---|---|
| `clusters/cluster_mean_profiles.png` | Stratification | Profils LAD moyens par cluster — clé pour interpréter les 3 types |
| `fpca/fpca_harmonics.png` | FPCA | Harmoniques FPC1-3 — valide la décomposition fonctionnelle |
| `fpca/fpca_reconstruction.png` | FPCA | Qualité de reconstruction (87.2% variance) |
| `fpca/fpca_scatters.png` | FPCA | Séparation des clusters dans l'espace fonctionnel |
| `gamm/gamm_marginal_effects.png` | GAMM | Effets marginaux des prédicteurs structurels et climatiques |
| `h1/h1_forward_curve.png` | H1 | Courbe d'inclusion progressive — hiérarchie LAI > Hmax > fCover > profil |
| `h2/c1_median_lad_profile.png` | H2 | Profil-type Cluster 1 (plot médian) |
| `h2/c2_median_lad_profile.png` | H2 | Profil-type Cluster 2 |
| `h2/c3_median_lad_profile.png` | H2 | Profil-type Cluster 3 |
| `h2/h2_distributions.png` | H2 | Distribution ΔTmax Real vs Uniform |
| `h2/h2_heatwave_histogram.png` | H2 | Effet LAD amplifié en canicule |
| `h2/h2_paired_real_vs_uniform.png` | H2 | Comparaison per-plot/jour Real vs Uniform — figure H2 centrale |
| `h2/h2_vertical_tmax_profiles_3way_canicule.png` | H2 | Profils 3 scénarios — figure comparative clé |
| `h2/h2_vertical_tmax_profiles_canicule_date.png` | H2 | Profils verticaux en canicule |
| `hobo/hobo_timeseries_representatives.png` | HOBO | Séries temporelles qualitatives HOBO vs MuSICA |
| `hobo/hobo_validation_rmse.png` | HOBO | RMSE par scénario vs terrain |
| `audit/h2_boxplot_cluster.png` | H2 | Structure de l'effet LAD par cluster |

### Figures pour annexe

| Fichier | Section | Description courte |
|---|---|---|
| `clusters/kmeans_elbow.png` | Méthodes | Choix du nombre de clusters (elbow) |
| `fpca/fpca_gcv_curve.png` | FPCA | Choix du nombre de bases B-spline (GCV) |
| `fpca/fpca_loadings.png` | FPCA | Loadings par hauteur |
| `h1/h1_hierarchy.png` | H1 | Hiérarchie complète de tous les scénarios |
| `h2/h2_verify_cluster1_peak.png` | H2 | Vérification robustesse pic sol Cluster 1 |
| `h2/h2_vertical_tmax_profiles_mean_summer.png` | H2 | Profils en conditions moyennes estivales |
| `h2/h2_vertical_tmax_profiles_median_date.png` | H2 | Profils à la date médiane |
| `audit/h2_amplitude_histogram.png` | H2 | Amplitude temporelle de l'effet LAD |
| `audit/h2_scatter_metrics.png` | H2 | Corrélations diff ↔ métriques structurelles |
| `s2_annex/s2_density.png` | H3 | Densités ΔTmax LiDAR vs Sentinel-2 |
| `s2_annex/s2_paired.png` | H3 | Scatter hexbin S2 vs LiDAR (figure principale annexe H3) |

### Diagnostics internes

| Fichier | Description courte |
|---|---|
| `audit/concurvity_full.csv` | Concurvité GAMM — worst-case par terme smooth |
| `audit/concurvity_pairwise_worst.csv` | Matrice pairwise concurvité |
| `gamm/gamm_residuals_diagnostic.png` | 4-panel diagnostic résidus Gaussien |
| `gamm/gamm_residuals_stats.csv` | Métriques diagnostic (Shapiro, kurtosis, ΔAIC) |
| `fpca/fpca_variance.txt` | Variance expliquée par les 3 FPC |
| `s2_annex/s2_metrics.txt` | Métriques R²/RMSE/Bias S2 |

### Données brutes pour ré-analyse

| Fichier | Colonnes clés | Lignes |
|---|---|---|
| `audit/h2_diff_correlations.csv` | metric, pearson, spearman, n | 3 |
| `audit/h2_structure_by_cluster.csv` | Cluster, period, n, median_diff, sd_diff, q05, q95 | 6 |
| `gamm/gamm_summary.txt` | Résumé complet du GAMM (R²=0.827, n=36600) | — |
| `h1/h1_scores.csv` | scenario, n, mae, rmse | 11 |
| `hobo/hobo_metrics.csv` | scenario, n, r2, rmse, mae, bias | 6 |
| `h2/h2_verify_cluster1_peak.csv` | nair, Tmax_z, height_m, scenario, Cluster, plot_id | 180 |
| `h2/h2_vertical_tmax_profiles_*.csv` (×3) | nair, Tmax_z, height_m, scenario, Cluster | 180 chacun |

### Fichiers obsolètes / à clarifier

| Fichier | Statut |
|---|---|
| `audit/h2_boxplot_archetype.png` | Doublon de `h2_boxplot_cluster.png` — vestige avant `refactor/rename-archetype-to-cluster` |
| `audit/h2_structure_by_archetype.csv` | Doublon de `h2_structure_by_cluster.csv` — même contenu, ancienne nomenclature |
