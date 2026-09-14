# Ch1 — carte producteur → figure (à rebours du manuscrit vivant, 2026-09-14)

Manuscrit : `manuscript_chap1_FINAL_coherence_2026-09-09.md`. Établi en partant des **14 figures
réellement embarquées**, puis en remontant au script qui écrit chacune. Aucun `mv` fait à partir de
cette carte (sauf les 12 backups datés, déjà archivés).

## 1. Les 14 figures embarquées et leur producteur RÉEL

| Figure manuscrit | Fichier | Producteur vivant |
|---|---|---|
| Fig 1 typologie | `figures/article_v323/Fig1_typology_gril_native.png` | `c1_fig1_gril_native.R` |
| Fig 2 méthodo | `figures/article_v323/Fig2_methodo.png` | `scripts/c1_fig2_methodo.R` (graphviz) |
| Fig 3 PCA clusters | `figures/article_v323/Fig3_pca_recover_clusters.png` | `c1_pca_recover_clusters.R` |
| Fig 4 attribution | `figures/Fig4_attribution_chs41.png` | `c1_fig4_fig5_chs41.R` |
| Fig 5 operating point | `figures/Fig5_operating_point_chs41.png` | `c1_fig4_fig5_chs41.R` |
| Fig 6 H2 gaussien | `figures/Fig6_h2_gaussian_chs41.png` | `c1_fig6v2_gaussian_heatmap.R` |
| Fig 7 profil vertical T | `figures/Fig7_vertical_Tprofile_chs41.png` | `c1_fig6_fig7_chs41.R` |
| Fig 8 rayon empreinte | `figures/Fig8_footprint_radius.png` | `c1_radius_rescore.R` → `scripts/c1_fig7_footprint_radius.R` |
| Fig A1 biais station | `figures/article_v323/FigA1_era5_station_bias.png` | `c1_era5_station_bias.R` |
| Fig A2 corr. vent | `figures/article_v323/FigA2_wind_profile_correction.png` | `c1_wind_profile_correction.R` |
| Fig B1 gradient/var | `figures/FigB1_pervariable_gradient_chs41.png` | `c1_b3_pervariable_chs41.R` |
| Fig C1 collinéarité | `figures/article_v323/FigC1_trait_collinearity.png` | `scripts/make_trait_collinearity_native20.R` |
| Fig D1 scan-angle | `figures/article_v323/FigD1_scanangle_control.png` | `scripts/fig_scanangle_control_csv.R` |
| Fig S1 attribution hot | `figures/FigS1_attribution_hot_chs41.png` | chaîne `_chs41` hot (à confirmer : `c1_hotdays_chs41.R`) |

**Producteurs vivants (à GARDER) :** les 13 ci-dessus + la chaîne `_chs41` amont qui alimente Fig 4-7/B1/S1
(`c1_archetype_profiles_chs41`, `c1_perturb_chs41_*`, `c1_h2_gaussian_*_chs41`, `c1_bias_origin_chs41`,
`c1_h1_realvsreal_chs41`, `c1_noise_floor_chs41`, `c1_h2_controlled_chs41*`, `c1_hotdays_chs41`).

## 2. Défaut concret : le ledger `make_article_figure_set.R` est PÉRIMÉ

Son mapping décrit une narration **pré-chs41** qui ne correspond plus au manuscrit :

| Ledger dit | Manuscrit embarque |
|---|---|
| Fig3 = `Fig4_attribution.png` (`c1_fig3_units.R`) | Fig3 = `Fig3_pca_recover_clusters.png` |
| Fig4 = operating_point (`c1_operating_point_main.R`) | Fig4 = `Fig4_attribution_chs41.png` |
| Fig5 = dumbbell (`c1_fig5_j1_convB.R`) | Fig5 = `Fig5_operating_point_chs41.png` |
| Fig6 = obs_corroboration | Fig6 = `Fig6_h2_gaussian_chs41.png` |
| (rien) | Fig7 chs41, Fig8, B1/S1 chs41 |

Le ledger **ne connaît aucune figure `_chs41`**. C'est la même maladie que les orchestrateurs
(`run_chapter1.R` ≠ `RUN_CHAPTER1.md`, aucun ne nomme la chaîne chs41) : **l'infra décrit un autre
chapitre.** Le dossier `figures/article_v323/` contient d'ailleurs ~17 PDF non embarqués (glass_ceiling,
forward_inclusion, dumbbell…) = sorties de cette narration morte.

→ **Correctif prioritaire (chantier code, pas rangement) :** réécrire le ledger `make_article_figure_set.R`
et un orchestrateur unique sur la narration réelle (Fig 1-8 + A1/A2/B1/C1/D1/S1, chaîne chs41 incluse), de
sorte que « la procédure documentée reproduise le Ch1 qui existe ».

## 3. Les 70 orphelins : résultat et pourquoi je ne les déplace PAS en masse

- **Aucun des 70 ne produit une figure embarquée** (intersection vide). C'est un signal fort *mais pas
  suffisant* : un orphelin peut écrire une **donnée** lue *par chemin* par un producteur vivant
  (ex. `c1_hot_extract.R` → `hot_extract.rds`), et il paraîtrait alors mort à tort.
- **git est inutilisable** (`c1_*` non suivis, mtimes mass-touchés) → pas de départage par récence.
- Le tri auto entrée/sortie de données par grep est **non fiable** (les fichiers canoniques partagés —
  `clhs_sample_native20_floor05.rds`, forçage `.nc`, `tab_hobo_perplot_cluster.csv` — sont lus par presque
  tout ; l'extraction confond inputs et outputs). Vérifié : le classeur auto donne des faux « KEEP ».

**Catégories des 70 (par nom, indicatif) :** slope ×6, diag ×6, era5 ×5, hobo ×5, cmp ×5, pointcloud ×5,
frblo ×4, _test ×4, v323 ×3, windcorr ×2, fpc ×2, glass_ceiling ×1, response_curve ×1, autres.

**Recommandation :** ne pas balayer par motif. Faire une passe *revue* (lire chaque orphelin : que
produit-il vraiment, quelqu'un le lit-il), par lots de ~10, Nathan valide chaque lot. Candidats de plus
haute confiance à traiter en premier (produisent une figure NON embarquée d'une narration abandonnée) :
`c1_fig4_glass_ceiling.R`, `c1_fig1_gril_style.R`, `c1_fig3_attribution_native.R`, `c1_cmp_era5_frblo*.R`,
`c1_safran_vs_era5_scatter.R` — à confirmer par lecture, pas par le nom.
