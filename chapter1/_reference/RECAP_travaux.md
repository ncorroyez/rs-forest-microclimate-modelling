# Chapitre 1 — récapitulatif des travaux

*Synthèse de tout ce qui a été fait sur cette période.*

## 1. Grand pivot : abandon du Shapley → analyse de sensibilité
- **Décision** : retrait **total** du Shapley du chapitre (manuscrit, figures, mémoires).
- **Pourquoi** : « le LAD mène dans le dense » était un résultat **spécifique au Shapley** (importance pondérée par la dispersion, off-manifold). Aucune lecture de sensibilité ne le reproduit → claim **adouci** : *le LAI domine partout, sa domination s'érode P1→P4, le profil vertical devient COMPARABLE (co-mène) en P4 dense sans jamais dépasser le LAI*.
- **Nouvelle attribution = sensibilité par perturbation** (sur la variété des données, sans baseline ni coalitions) :
  - `c1_importance_sensitivity.R` → **Fig 2** (importance = |sensibilité médiane par unité| × dispersion intra-cluster ; LAD = |réel-vs-uniforme|). P1 LAI 0.57 → P4 LAI 0.12 ≈ LAD 0.11.
  - **Fig 3** = sensibilités par plot (unités physiques) ; **Fig 5** forward-inclusion régénérée avec l'ordre sensibilité (LAI 1er partout).
  - Appendice A : potency (universelle) + sensibilité locale (saturation ×4) — remplace l'ancienne annexe Shapley-vs-LOO (supprimée).
- **Manuscrit entièrement réécrit** (titre, abstract, §1, §2.4, §2.6, §3.2, §3.3, §4.1–4.4, §5, annexes) ; mémoires mises à jour (Shapley marqué historique).

## 2. Intégration du VCI pour le terme LAD
- Confirmé : **VCI = 1 ⟺ profil uniforme** (Shannon-evenness).
- Le contraste LAD réel-vs-uniforme = « VCI_réel → 1 » ; intra-archétype, plus le profil est concentré (VCI bas) plus le contraste est grand : **ρ partiel = −0.33** (|LAI). Le **centre de masse échouait** (≈0.5 partout) — le VCI réussit.
- `c1_lad_vci.R` → **Fig A2** ; intégré en §2.6/§3.2/§4.3 (VCI = descripteur/normalisateur, **pas** un input MuSICA ; nécessaire-mais-pas-suffisant → profil complet via FPCA).

## 3. Comparaison binaire v3.2.0 (validé) vs v3.2.3 (yoyo / `ABL_flag='iter'`)
- **Demande** : tester le « nouveau musica » + le mode yoyo.
- **Mise en route du yoyo** (diagnostic via build de debug puis binaire officiel restauré) : 2 pré-requis trouvés —
  - phénologie du **jour précédent** (2020/366, `calc_phenology` ne fait que 365 j/an) → sinon `STOP 16` ;
  - variable **`h_sbl`** (hauteur de couche limite) dans le forçage → sinon `STOP 6`.
- **Vraie BLH** récupérée : **MERRA-2 PBLH** (NASA Earthdata, après autorisation app GES DISC), cellule Blois, horaire 2021-01→2022-07, injectée → `in_files/musica_in_Blois_pblh.nc` (`get_merra2_pblh.sh` + `build_forcing_pblh.R`). *(ERA5/CDS inaccessible : pas de clé.)*
- **Résultats** (3200 simus iter + validation HOBO) :
  - **iter ≈ none** : importance identique (Δ≤0.11 °C), validation r=0.633 ≈ none 0.629 → **le yoyo est immatériel** pour le ΔT_max à 1 m (confirmé avec vraie BLH).
  - **v3.2.3 valide nettement moins bien** : r **0.93 → 0.63** ; et surpondère le LAD (mène en P3/P4).
  - **Robuste aux 3 métriques × 2 périodes** (ΔTmax, pente, ΔVPDmax × tous/10%-chauds) : contraste v3.2.0(LAI)/v3.2.3(LAD-dense) tient partout.
- **Décision (sur preuve)** : **garder le legacy v3.2.0 validé** ; yoyo documenté comme contrôle de robustesse, écarté.
- Livrables : `comparaison_versions/` (doc + `FigCmp_hobo_validation`, `FigCmp_importance`, `FigCmp_metrics6` + tables).

## 4. Liste de la réunion — statut
- **Tout le calculable est fait** : profils 0.5 m (validation), sensibilité ±LAI/ΔHmax/ΔfCover, VCI + boxplot, ACP récupère P1-P4, Lmom/Hmed, LAD↔VCI, normalisation range/pas, cohérence MuSICA↔HOBO, microclimat↔clustering, valeur du LiDAR, HOBO dans le récit, clustering FPC, cLHS-vs-moyenne, Fig 1 corrigée (homothétie H_max).
- **Caduc** (Shapley abandonné) : tests médiane/moyenne + côté baseline cLHS-vs-moyenne.
- **Reste** : (1) lien écologique stade/gestion *quantitatif* = données ONF externes ; (2) 4 trous factuels manuscrit (modèle HOBO, specs LiDAR, réfs Ramsay&Silverman + Hersbach) ; (3) attribution 0.5 m (Opt 2) **parquée**.

## 5. Lien écologique des 4 profils — inséré + validé à 1 m (§4.4)
- Ouverture **affirmative** : l'ACP {LAI, fCover, H_max, VCI} **récupère** les 4 clusters (PC1≈75%) ; le **microclimat mesuré ordonne** les types (**KW p≈3×10⁻⁶**) → typologie LiDAR 1 m prédit le régime thermique réel ⇒ sens écologique établi **à 1 m**, sans 0.5 m.
- Lecture stade/gestion = hypothèse testable (nécessite inventaire). Figures supplémentaires `FigS1_pca_clusters`, `FigS2_microclimate_by_cluster` ajoutées + documentées.
- **Décision** : tout reste **à 1 m** (attribution) ; 0.5 m en validation ; robustesse documentée.

## 6. Nettoyage disque (~95 G libérés : 47 G → 115 G)
- **Shapley devenu inutile (~49 G)** : `nc_shapley2x_global` (26 G) + les 5600 coalitions ≠1111/1110 de `nc_shapley2x` (~23 G). **Conservé** : `*_1111.nc` + `*_1110.nc` (réutilisés par la sensibilité).
- **H1_* pré-refonte (~46 G)** : `H1_factorial/forward/loo/lovb_baselines` (sur ta validation).
- Reste intact (sur ta demande) : option_b ERA5, vieux HOBO v2-4/v9-11, variantes Chapter3, TS, Annex_S2, nc de comparaison (~30 G).

## 7. Organisation projet + mémoires
- **`Chapitre1/`** créé : `redaction/` (manuscrit prose + **version bullet points** + fiche `README_figures_methodes.md` détaillée par figure + figures + bib) et `comparaison_versions/` (fiche + figures + tables) + `README.md` index.
- **Convention** : tout livrable s'accompagne d'une fiche 5-champs (méthodes · pré-traitements · données/fig · description · interprétation).
- Mémoires créées/mises à jour : pivot Shapley→sensibilité, recette iter/yoyo + BLH, binaire validé, structure de dossier + convention.

## 8. État final
- **Article** : manuscrit prose réécrit (sensibilité, co-lead, VCI, plafond de verre, lien éco validé à 1 m) + version bullet points, dans `chapter1/manuscript/`.
- **Robustesse** : binaire (v3.2.0 vs v3.2.3/yoyo), résolution (1 m vs 0.5 m), 3 métriques × 2 périodes — tout testé et documenté.
- **Prêt pour la réunion / Nextcloud**, modulo les 4 blancs factuels à remplir et la décision (déjà prise) de rester à 1 m.
