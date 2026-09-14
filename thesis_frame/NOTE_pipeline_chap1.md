# Chapitre 1 — note pipeline (très brève)

## But
Déterminer **quelle dimension de la structure de canopée gouverne le tampon microclimatique estival** de
sous-bois : la **quantité de feuilles** (LAI, couverture) ou leur **arrangement vertical** (profil LAD).
Approche : traits LiDAR (ALS) → modèle biophysique **MuSICA** → **attribution Shapley exacte** de ΔTmax,
validée sur **53 loggers HOBO** (chênaie de Blois, été 2021).

## Hypothèses
- **H1** : le profil LAD module les transferts radiatifs/turbulents → paramètre clé pour une simulation fidèle.
- **H2** : à LAI et Hmax fixés, un profil **top-heavy** atténue davantage le cycle diurne (interception haute).

## Méthodes (résumé)
- **2.1** Site Blois ; 53 HOBO à 1 m (boucliers PVC, été 2021) ; ΔTmax = max micro − max macro (journalier).
- **2.2** Profils LAD 0.5 m depuis ALS leaf-on ; variables LAI, Hmax, fCover, VCI (VCI → SupMat) ; grille 20 m.
- **2.3** FPCA du profil (double-normalisé → forme pure) ; K-means sur **LAI, Hmax, fCover, FPC1–3** (~87 % var)
  → **4 archétypes P1–P4** ; échantillon **cLHS 100/cluster = 400 plots**.
- **2.4** MuSICA multicouche ; **16 scénarios/plot** (chaque variable réelle vs baseline) = 6400 sims ; baseline
  scalaires = moyenne, LAD = uniforme (préserve LAI et Hmax, n'aplatit que la forme) ; 2 baselines : **par
  cluster** (principal) et **global 400 plots** (complémentaire).
- **2.5** Métriques **toutes à 1 m** : (a) ΔTmax/ΔVPDmax (max journalier, −2 h) ; (b) pente micro-macro
  (régression horaire, no shift ; log(slope)<0 = tampon, >0 = amplification).
- **2.6** **Shapley exact** (φ, additivité à la précision machine) ; importance = mean|φ| (baseline cluster),
  direction = φ signé (baseline global) ; contraste **partenaire-fixé** pour borner l'artefact off-manifold
  (LAD|Hmax, LAI|fCover).
- **2.7** Validation vs 53 HOBO (ΔTmax, pente) ; forward-inclusion ; Moran's I ; résidu ~ topo/hétérogénéité.

## Résultats + comment chaque figure est faite

| Fig | Message | Script / source |
|---|---|---|
| **1** Archétypes | 4 profils LAD moyens P1→P4 (ouvert→dense) | `08_typology` / `outputs/figs_MEB2026_final/fig_archetypes_profiles.png` |
| **2** Importance par cluster | **density-dependent** : P1 LAI 0.41 ≫ ; P4 **LAD 0.13 mène** (LAI/fCover saturent) | `pipeline/11_clhs_attribution.R` ← parts `c3_shapley_chunk.R` (floor05_v2, 6400 nc cachés) |
| **2b** Heatmap 3 métriques × 2 périodes | bascule LAI→LAD tient partout | `c1_metrics_chunk.R` → `c1_metrics_merge.R` |
| **3** Importance paysage (baseline global) | LAI 0.54 > fCover 0.31 ≫ LAD 0.09 (mean\|φ\|) | `pipeline/11_clhs_attribution.R` ← `c1_shapley_global_chunk.R` |
| **4** Validation ΔTmax | r=0.93, biais +0.9 °C, split **44/9**, pente 0.95 vs 0.86 | `pipeline/07_validation.R` (branche **z05**) |
| **5** Forward-inclusion (grille P1–P4 + All) | trait dominant de chaque type construit l'ajustement (r entre loggers) | `c1_forward_percluster.R` (branche z05) |
| **6** Profils verticaux (vent/RH/VPD) | illustration mécaniste | `outputs/figures_pipeline/annex/` |
| **A1** Shapley vs forward vs LOO | Shapley borné par les 2 limites ; justifie l'estimand | `c1_loo_forward_chunk.R` → `c1_loofwd_merge.R` |
| **B1** Résidu vs gap_fraction | plafond de verre : adj R²=0.62, gap r=0.76 (53 plots) | `scripts/make_residual_vs_topo.R` (covariables Not_Masked) ; Moran : `scripts/make_moran_validation.R` |
| **C1** Sensibilité au rayon | fit insensible au footprint 5–50 m (cercle/carré), 12.5 m justifié | `c1_radius_{circle,square}.R` → `c1_radius_metric.R` → `c1_radius_figure.R` |
| **D** (table) Dispersion par archétype | SD LAI plate (1.25–1.52) → effondrement = saturation, pas échantillonnage | `outputs/figures_pipeline_z05/tables/tab_trait_dispersion_by_archetype.csv` |
| **E1/E2** Corrplots traits (+ FPC1–3) | colinéarité forte au paysage, **s'effondre dans P4** ; FPCs quasi indépendants des scalaires | `pipeline/12_corrplot_traits.R` ; `pipeline/12b_corrplot_traits_fpc.R` |

**Conclusion** : le contrôle est **density-dependent** — quantité de feuilles dans l'ouvert, arrangement
vertical dans le dense (par saturation, effet petit en absolu) ; MuSICA reproduit le **classement** (r=0.93)
mais bute sur le **niveau absolu** (plafond de verre = trouées/topo hors d'un modèle 1-D). Pont vers le Chap 2
(Sentinel-2 : illusion optique).
