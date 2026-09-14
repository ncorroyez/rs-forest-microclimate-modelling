# Comparaison MuSICA v3.2.0 vs v3.2.3 (yoyo / `ABL_flag='iter'`)

*Chapitre 1 — Blois, été 2021. Objet : tester si le couplage couche-limite atmosphérique
(le « yoyo » de MuSICA, disponible uniquement dans la v3.2.3) change l'attribution des traits
et/ou améliore la validation terrain, par rapport au binaire validé v3.2.0.*

*Numérotation : Figures 1–4 = contexte (typologie / création des 400 plots) ; Figures 5–12 = comparaison.*

---

## 1. Méthode (en bullet points)

**Les deux binaires comparés**
- **v3.2.0 « legacy »** : `/home/corroyez/Documents/musica/musica` (md5 `5307…`, nov. 2024). Binaire **validé** des chapitres ; **antérieur à `ABL_flag`** → **pas de yoyo**.
- **v3.2.3** : `in_files/model-3.2.3/musica` (md5 `8b5139e0`, mai 2026), lancé avec **`ABL_flag='iter'`** = yoyo itératif (`itermax=2` par pas de temps).

**Mise en route du yoyo (deux pré-requis, sinon le modèle plante)**
- **Phénologie du jour précédent** : au 1ᵉʳ pas, le yoyo lit la veille → **2020/366** (bissextile). `calc_phenology` ne fait que 365 j/an → ligne dormante ajoutée, sinon `STOP 16`.
- **`h_sbl`** (hauteur de couche limite) dans le forçage → sinon `STOP 6`. **Vraie BLH** = **MERRA-2 PBLH** (NASA/Earthdata), Blois, horaire 2021-01→2022-07, injectée dans `in_files/musica_in_Blois_pblh.nc`. *(ERA5/CDS inaccessible.)*

**Métriques (3) × périodes (2)** — toutes à **1 m fixe**, été 2021 (JJAS)
- **ΔT_max** = max diurne T micro − T macro, décalage **−2 h**. Négatif = tampon.
- **slope micro/macro** = pente régression horaire (T micro ~ T macro), sans décalage. <1 = tampon.
- **ΔVPDmax** = max diurne du VPD reconstruit à 1 m. *(Non validable HOBO — capteurs T seuls.)*
- **Périodes** : tous les jours / **10 % les plus chauds** (T_max macro ≥ p90). Les 6 combos sont ré-extraits des nc (aucune simu nouvelle).

**Création des 400 plots cLHS (chaîne LiDAR → typologie → échantillon)** — *cf. Figures 1–4*
- **LiDAR ALS leaf-on** → profils **LAD** (inversion du profil de retours) + traits LAI, H_max, fCover, VCI. Grille 20 m, LAD binné 1 m.
- **FPCA** sur les profils LAD double-normalisés (forme pure) → **FPC1–3** (~87 % de la variance de forme).
- **K-means** sur **6 métriques** (LAI, H_max, fCover, FPC1–3) → **4 archétypes P1→P4**, gelés avant simu.
- **cLHS** : 100 plots/archétype → **400 plots**. fCover planché à 0,5.

**Analyse de sensibilité — comment**
- **Sur la variété des données** : chaque trait bougé d'un petit pas natif autour du point réel de chaque plot (partenaires au réel) → réponse partielle locale, **jamais off-manifold** (pas de baseline globale, pas de coalitions).
- **LAD** = seule exception (input fonctionnel) : effet = **contraste réel-vs-uniforme** (uniforme = VCI maximal).
- 400 plots : base + 6 perturbations (**LAI ±1 one-sided, fCover ±0,1, H_max ±5 m**) + 1 contraste LAD = **8 simus/plot** (3200/binaire). Mêmes 400 plots pour les deux.
- **Importance** = |pente locale médiane par unité| × **dispersion intra-archétype** (LAD = |médiane réel−uniforme|). *(Readout cLHS ; voir aussi « P1–P4 only », Fig 10.)*

---

# Figures 1–4 — Contexte : typologie & création des 400 plots

## Figure 1 — Archétypes LAD P1–P4
`FigCtx1_archetypes.png`
- **Données** : profils LAD (1 m) des 400 plots cLHS + labels de cluster.
- **Process** : profil LAD **moyen par cluster**, ramené au H_max moyen du cluster par **homothétie** (`.lad_rescale`) avant moyennage (correction V3).
- **Description** : P1 ouvert (LAI 1,2 ; fCover 0,55 ; H_max 20,7 ± 10,7 m) → P4 dense (5,4 ; 1,0 ; 33,7 m). H_max **non monotone** (P2 le plus court, 14,8 m).
- **Interprétation** : gradient structural densité/fermeture ; base de l'analyse densité-dépendante.

## Figure 2 — Modes fonctionnels du profil (FPCs)
`FigCtx2_fpc_harmonics.png`
- **Données** : profils LAD double-normalisés (hauteur + densité) des 400 plots.
- **Process** : FPCA → harmoniques **FPC1–3** (modes de variation de la *forme* du profil).
- **Description** : 3 modes capturent ~87 % de la variance de forme (FPC1 ≈ mode dominant de répartition verticale ; FPC2–3 = modes secondaires).
- **Interprétation** : la forme du profil (fonction) est résumée par 3 axes **orthogonaux**, largement **indépendants des scalaires** (LAI/H_max/fCover) → entrées propres du clustering (pourquoi un scalaire seul ne suffit pas).

## Figure 3 — ACP : le clustering est récupéré (validation structurale)
`FigCtx3_pca_clusters.png`
- **Données** : traits scalaires (LAI, fCover, H_max) + **VCI** sur les 400 plots, labels de cluster.
- **Process** : ACP non supervisée ; projection des plots coloriés par archétype.
- **Description** : les **4 clusters se séparent** le long de PC1 (~75 % de variance = axe densité/fermeture portant LAI+fCover+VCI) ; PC2 = H_max. **VCI colinéaire** à l'axe densité ; P3/P4 se chevauchent.
- **Interprétation** : une ACP **indépendante** retrouve les clusters → typologie **non arbitraire**. Mais un **scalaire (VCI) ne crée pas d'axe séparateur** indépendant → justifie le profil complet (FPCA).

## Figure 4 — Boxplots du VCI par cluster
`FigCtx4_vci_boxplot.png`
- **Données** : VCI (indice de complexité verticale) par plot, groupé par archétype.
- **Process** : boxplot du VCI par cluster (P1–P4).
- **Description** : **gradient monotone** médiane P1 ≈ 0,33 → P4 ≈ 0,60 ; dispersion qui se resserre de P1 (large) à P4 (étroite).
- **Interprétation** : les clusters capturent bien un **axe d'inhomogénéité verticale** (VCI) — métrique de forme **nécessaire mais pas suffisante** (d'où le clustering sur le profil complet).

---

# Figures 5–12 — Comparaison v3.2.0 vs v3.2.3 (yoyo/iter)

## Figure 5 — Validation HOBO : simulé vs observé (ΔTmax)
`FigCmp_hobo_validation_v320_v323iter.png`
- **Données d'entrée** : ΔT_max **observé** (53 HOBO) ; **simulé** à la config réelle (z05) par chaque binaire.
- **Process** : un point = un logger (obs en x, sim en y) ; 1:1, ajustement, r/biais/RMSE sur 53 plots.
- **Description** : v3.2.0 nuage serré, **r = 0,93**, biais +0,92, RMSE 1,47. v3.2.3 iter **quasi plat**, **r = 0,63**, biais +0,61, RMSE 1,67.
- **Interprétation** : v3.2.0 **reproduit le classement entre-plots** (base du glass ceiling). v3.2.3 **sur-lisse** (perte de structure spatiale).

## Figure 6 — Importance des traits par unité, par archétype
`FigCmp_importance_v320_v323iter.png`
- **Données d'entrée** : importance par trait (LAI, fCover, LAD, Hmax) et archétype, par binaire (perturbations per-unit, 400 plots).
- **Process** : barres groupées ; importance = |pente médiane par unité| × dispersion intra-archétype (LAD = |réel-vs-uniforme|).
- **Description** : **v3.2.0** LAI domine, s'érode P1→P4 (0,57 → 0,12), LAD co-mène en P4 (0,12 ≈ 0,11). **v3.2.3 iter** P1 s'effondre, **LAD mène P3/P4** (P4 : LAD 0,80 ≫ LAI 0,39).
- **Interprétation** : le **yoyo ne change quasi rien** (iter ≈ none) ; le **vrai écart est le binaire** : v3.2.3 surpondère le vertical — sur le binaire qui valide le moins bien (Fig 5).

## Figure 7 — Trait dominant sur 3 métriques × 2 périodes (heatmap)
`FigCmp_metrics6_v320_v323iter.png`
- **Données d'entrée** : importance per-unit des 4 traits pour les **6 combos** métrique × période, par archétype et version.
- **Process** : trait dominant (couleur) + part du LAD (étiquette) ; lignes = combos, colonnes = archétypes, facettes = version. (même formule que Fig 6.)
- **Description** : **v3.2.0** LAI dominant ~partout (LAD ne passe en tête qu'en P4 chaud) ; **v3.2.3 iter** LAD mène en P3/P4 sur les 6 combos (42–66 %). Jours chauds renforcent.
- **Interprétation** : le contraste v3.2.0(LAI)/v3.2.3(LAD) **n'est pas un artefact d'une métrique unique** — il tient sur ΔTmax, slope et ΔVPDmax.

## Figure 8 — Scatter horaire toutes les températures
`FigCmp_scatter_hourly_v320_v323iter.png`
- **Données d'entrée** : **toutes** les T air horaires à 1 m, sim vs obs, été 2021 (**155 016 points**), par version.
- **Process** : hexbins (densité log), 1:1 + ajustement ; R²/RMSE/biais.
- **Description** : **v3.2.0 R² 0,78** / RMSE 2,39 / biais +1,13 ; **v3.2.3 iter R² 0,77** / 2,41 / +1,09 — quasi identiques.
- **Interprétation** : sur la T brute (cycle diurne+saisonnier) les deux sont équivalents ; la différence de version est dans la **structure fine entre-plots** (Fig 5), pas le suivi grossier.

## Figure 9 — Validation par plot : ΔTmax & slope (tous / 10 % chauds)
`FigCmp_scatter_metrics_v320_v323iter.png`
- **Données d'entrée** : par plot (53), ΔTmax et slope, obs vs sim, tous jours et 10 % chauds.
- **Process** : R²/RMSE par panneau ; **panneaux carrés**, x = y range (1:1 vraie diagonale). 8 panneaux (4 combos × 2 versions).
- **Description** (R² | RMSE) : ΔTmax tous **0,87**/1,47 vs 0,40/1,67 ; ΔTmax chauds 0,75/2,06 vs 0,32/2,12 ; slope tous **0,92**/0,14 vs 0,61/0,14 ; slope chauds 0,87/0,22 vs 0,50/0,22.
- **Interprétation** : **v3.2.0 valide mieux sur les 4 combos** ; le slope (amplitude) est le mieux reproduit (R² 0,92).

## Figure 10 — Sensibilité « P1-P4 only » (4 profils moyens)
`FigCmp_sensitivity_archetypes_v320_v323iter.png`
- **Données d'entrée** : les **4 profils MOYENS d'archétype** (au lieu des 400 plots). 8 simus/archétype × 2 versions = 64.
- **Process** : effet par unité au **point de fonctionnement moyen** de chaque archétype (sensibilité brute, pas d'importance×dispersion). = le « 2/ moy » du « 1/ cLHS 2/ moy » de la réunion.
- **Description** : v3.2.0 effets modérés **sauf P1 où LAD = 1,12** ; v3.2.3 iter effets bien plus grands (LAD 1,2–2,0).
- **Interprétation** : **P1 LAD = 1,12 = artefact de Jensen** (le profil moyen d'un archétype hétérogène n'est pas représentatif ; par plot, LAD P1 ≈ 0,02 et s'annule) → confirme le **readout cLHS (400 plots) comme principal** ; le « P1-P4 only » illustre justement le biais de la moyenne.

## Figures 11–12 — Statistiques & tests détaillés
Skill (r ± IC bootstrap 95 %) :
`FigCmp_stats_skill_v320_v323iter.png`
Erreur absolue appariée (Wilcoxon) :
`FigCmp_stats_error_v320_v323iter.png`
Tables : `tab_stats_v320_v323.csv`, `tab_stats_tests_v320_v323.csv`.
- **Données / process** : par plot (53). r/R²/RMSE/MAE/biais + **IC bootstrap 95 %** (3000 rééch.). Tests appariés : **Wilcoxon** sur |erreur|/plot ; **bootstrap Δr**. Script `c1_cmp_stats.R`.
- **Tableau (r [IC95] · RMSE · biais)** :
  | Métrique | v3.2.0 r [IC95] | v3.2.3 r [IC95] | Δr [IC95] (p) | RMSE 320/323 | biais 320/323 |
  |---|---|---|---|---|---|
  | ΔTmax · tous | 0.93 [0.88–0.97] | 0.63 [0.57–0.72] | **+0.30 [0.19–0.37]** (p<0.001) | 1.47 / 1.67 | +0.92 / +0.61 |
  | ΔTmax · chauds | 0.86 [0.76–0.93] | 0.57 [0.48–0.69] | **+0.30 [0.12–0.41]** (p=0.005) | 2.06 / 2.12 | +1.56 / +1.10 |
  | slope · tous | 0.96 [0.93–0.98] | 0.78 [0.72–0.85] | **+0.18 [0.11–0.23]** (p<0.001) | 0.144 / 0.137 | +0.082 / +0.030 |
  | slope · chauds | 0.93 [0.88–0.96] | 0.71 [0.65–0.78] | **+0.22 [0.15–0.27]** (p<0.001) | 0.222 / 0.224 | +0.133 / +0.085 |
  - **T horaire 1 m** : R² ≈ 0,89, RMSE ≈ 1,91 °C, r/plot ≈ 0,945 pour **les deux** ; Wilcoxon RMSE/plot **p = 0,61 (aucune différence)**.
- **Tests — deux axes distincts** :
  - **Ranking spatial (r)** : v3.2.0 **significativement meilleur** sur les 4 — IC bootstrap **disjoints** (Fig 11), Δr > 0 à p ≤ 0,005.
  - **Niveau absolu (|erreur|, biais)** : c'est **v3.2.3 qui a l'erreur la plus faible** (biais chaud moindre) — significatif pour slope (p≈0,003/0,004) & ΔTmax-chauds (p=0,009) ; ΔTmax-tous ns.
  - **T horaire brute** : aucune différence (p=0,61).
- **Interprétation** : les binaires se départagent sur des **axes différents** — v3.2.0 = meilleure **structure spatiale** (r) ; v3.2.3 = meilleur **niveau** (biais, dominé par le forçage). L'attribution du chapitre reposant sur le **classement entre-plots**, **v3.2.0 est retenu** ; l'avantage de niveau de v3.2.3 ne compense pas la perte de structure (et ne vient pas du yoyo, immatériel).

---

## Conclusion

- **Le yoyo est immatériel** pour le ΔT_max à 1 m : iter ≈ none (sensibilité **et** validation), confirmé avec une **vraie BLH** (MERRA-2).
- « Utiliser le yoyo » = **adopter v3.2.3**, qui **valide moins bien en ranking** (r 0,93 → 0,63) et **change la conclusion** (LAD mène dans le dense).
- **Robuste aux 3 métriques × 2 périodes** (Fig 7) : pas un artefact d'une métrique unique.
- **« Plus complet » n'est pas « plus juste » ici** : le couplage ABL valide moins bien la structure spatiale.
- **Décision (sur preuve)** : **v3.2.0 (validé, sans yoyo)** = binaire de l'article ; le yoyo (vraie BLH) **testé et assumé** en robustesse. Réponse solide au reviewer « pourquoi ignorer le couplage ABL ? ».
