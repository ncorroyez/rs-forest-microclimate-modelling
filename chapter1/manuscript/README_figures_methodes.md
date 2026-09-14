# Chapitre 1 — Méthodes, pré-traitements et figures (fiche détaillée)

*Compagnon de `manuscript_chap1_EN.md`. Pour CHAQUE figure : données d'entrée · pré-traitements ·
process (script) · description · interprétation. Objet du chapitre : attribuer le tampon microclimatique
estival de sous-bois aux traits de structure de canopée (LiDAR → MuSICA), validé sur 53 HOBO. Blois, été 2021.*

---

## 0. Conventions transversales

- **Binaire MuSICA** : **legacy v3.2.0** (`/Documents/musica/musica`, md5 5307…), le seul **validé** (r=0.93 vs HOBO). Le couplage ABL « yoyo » (v3.2.3) a été testé et écarté → voir `../comparaison_versions/`.
- **LAI** : le LiDAR donne le LAI *one-sided* ; MuSICA attend le total *two-sided* → ×2 appliqué **une seule fois** dans `run_musica_one` (R/musica.R).
- **Métrique ΔT_max** : (max diurne T micro à **1 m** interpolé) − (max diurne T macro), été 2021 (JJAS), décalage **−2 h** (aligne les pics simulés/observés). Négatif = tampon.
- **Métriques complémentaires** : pente micro/macro (régression horaire, sans décalage ; <1 = tampon) ; ΔVPDmax (1 m, reconstruit de T + `wair_z`). **Périodes** : tous les jours / 10% les plus chauds (T_max macro ≥ p90).
- **Échantillon** : `out_files/Sensitivity_Analysis/clhs_sample_floor05_v2.rds` — 400 plots cLHS (100/archétype), fCover planché à 0.5.

## 1. Pré-traitements communs (chaîne LiDAR → MuSICA)

- **LiDAR ALS leaf-on** → profils de **densité de surface foliaire (LAD)** par inversion du profil de retours [@macarthur1969; @bouvier2015] ; traits scalaires LAI, H_max, fCover, VCI. Agrégation grille **20 m**, LAD **binné à 1 m**. *(Validation HOBO : reconstruction fine **0.5 m** « z05 » depuis le nuage brut, capte le sous-bois <2 m.)*
- **FPCA** sur les profils LAD (double-normalisés en hauteur ET densité → forme pure, sans échelle) → **FPC1–3** (~87% de la variance de forme).
- **K-means** sur 6 métriques (LAI, H_max, fCover, FPC1–3) → **4 archétypes P1→P4** (densité croissante). Clustering **gelé avant toute simulation**.
- **cLHS** : 100 plots/archétype → **400 plots** couvrant l'espace des traits.
- **Design sensibilité (par plot, 8 simus)** : config réelle (base) + perturbations **±1 LAI (one-sided)**, **±0.1 fCover**, **±5 m H_max** (un trait à la fois, les autres au réel) + contraste **LAD réel-vs-uniforme**. Sur la **variété des données** (jamais de baseline globale off-manifold).
- **Importance d'un trait dans un archétype** = |pente locale médiane par unité| × **dispersion intra-archétype** du trait ; LAD = |médiane du contraste réel-vs-uniforme|.
- **Validation** : 53 HOBO, run MuSICA dédié par logger (LAD 0.5 m), ΔT_max sim vs obs.

---

# Figures principales

## Fig. 1 — Quatre archétypes structuraux (profils LAD)
- **Données** : profils LAD (1 m) des 400 plots cLHS + labels de cluster.
- **Pré-traitements** : typologie FPCA + K-means (§1).
- **Process** : profil LAD **moyen par cluster**, chaque profil **ramené au H_max moyen du cluster par homothétie** (`.lad_rescale` : étirement axe hauteur + renormalisation densité à LAI constant) **avant** moyennage (correction « V3 » du bug d'étirement). LAI affiché one-sided. Script `scripts/make_fig_archetypes_profiles_v3.R` → `Fig01_archetypes.png`.
- **Description** : P1 ouvert (LAI 1.2, fCover 0.55, H_max 20.7 ± 10.7 m) → P4 dense (LAI 5.4, fCover 1.0, H_max 33.7 m) ; H_max **non monotone** (P2 le plus court, 14.8 m).
- **Interprétation** : gradient structural densité/fermeture ; la hauteur n'ordonne pas la densité (P1 = classe ouverte hétérogène, pas « jeune peuplement »). Justifie de lire l'attribution le long du gradient de densité.

## Fig. 2 — Importance des traits par archétype (sensibilité par unité)
- **Données** : simulations per-unit des 400 plots (`out_files/Chapter3/nc_sensitivity_perplot` + base/uniforme `nc_shapley2x/*_{1111,1110}.nc`).
- **Pré-traitements** : 8 simus/plot (§1).
- **Process** : importance = |sensibilité médiane par unité| × dispersion intra-cluster (LAD = |médiane réel-vs-uniforme|). Scripts `c1_sensitivity_perplot_chunk.R` (+ `_lad`) → `c1_importance_sensitivity.R` → `Fig02_importance_sensitivity.png`.
- **Description** : P1 LAI 0.57 ≫ fCover 0.24 > H_max 0.18 > LAD 0.02 → P4 **LAI 0.12 ≈ LAD 0.11** (LAD a dépassé cover/hauteur). Érosion monotone du LAI P1→P4.
- **Interprétation** : **le LAI domine partout, mais sa domination s'érode** ; dans le dense saturé, le profil vertical devient **comparable** (co-mène) — sans jamais dépasser le LAI.

## Fig. 3 — Sensibilités par plot, en unités physiques (distributions)
- **Données** : mêmes simulations per-unit.
- **Pré-traitements** : effet per-unit (perturbation − base)/pas, par plot, par direction.
- **Process** : boxplots par archétype des effets ± (LAI/fCover/H_max) + contraste LAD réel-vs-uniforme. Bornes fCover gérées (saturation ≥0.95 → 0 ; plancher 0.5 → exclu). Script `c1_sensitivity_perplot_figure.R` → `Fig03_sensitivity_perplot.png`.
- **Description** : ajouter 1 LAI refroidit ~0.47 °C (ouvert) → 0.17 °C (dense) ; **asymétrie** ouvert : retirer 1 LAI = +1.18 °C vs ajouter −0.47 °C (×2.5) ; LAD réel-vs-uniforme +0.02 (P1) → +0.11 (P4).
- **Interprétation** : **saturation** (le pouvoir refroidissant de la quantité s'effondre vers le dense) + **asymétrie** (perdre du feuillage coûte plus que densifier ne rapporte, dans l'ouvert).

## Fig. 4 — Validation sur 53 HOBO (simulé vs observé)
- **Données** : ΔT_max observé (53 HOBO, été 2021) + simulé à la config réelle, LAD **0.5 m**, un run MuSICA par logger (`out_files/musica_hobo_z05`).
- **Pré-traitements** : LAD reconstruit à 0.5 m (capte sous-bois <2 m) ; ΔT_max 1 m / −2 h.
- **Process** : sim vs obs par plot, r/biais ; partition tampon/amplification. Pipeline `03z`/`07`.
- **Description** : **r = 0.93**, biais chaud +0.9 °C, partition 44/9 (obs 45/8), pente médiane 0.95 vs obs 0.86.
- **Interprétation** : le modèle **capture le classement spatial** mais **plafonne sur le niveau absolu** et comprime l'amplitude → le **« plafond de verre »**.

## Fig. 5 — Validation forward (inclusion cumulative des traits)
- **Données** : coalitions HOBO (`musica_hobo_z05`) + **ordre d'importance** (`tab_importance_sensitivity.csv`).
- **Pré-traitements** : à partir de la canopée baseline de chaque archétype, ajout cumulatif des traits dans **l'ordre d'importance sensibilité**.
- **Process** : à chaque étape, r entre ΔT_max simulé moyen et observé, par archétype. Script `c1_forward_percluster.R` → `Fig05_forward_inclusion.png`.
- **Description** : **LAI ajouté en 1er dans tous les archétypes** (r 0.84–0.98 seul) ; en P3/P4 le **LAD ajouté en 2e** relève le fit (P4 : r 0.86 → 0.90).
- **Interprétation** : le trait qui domine l'attribution est aussi celui qui **récupère le terrain en premier** — pont entre attribution interne et validation.

## Fig. 6 — Profils verticaux simulés (illustration mécaniste)
- **Données** : nc MuSICA, **15 couches** (seule figure hors lecture 1 m).
- **Process** : profils verticaux T air / vent / HR / VPD par archétype. `fig_vertical_profiles_legacy.png`.
- **Description / Interprétation** : illustre l'atténuation verticale ; pas de référence verticale terrain (illustratif).

---

# Annexes

## Fig. A1 — Potency intrinsèque (a) vs sensibilité locale normalisée (b)
- **Données** : simulations per-unit (a : pleine gamme, partenaires à la moyenne du cluster ; b : fenêtre intra-cluster p10–p90, normalisée par le range plein-échantillon).
- **Process** : `c1_sensitivity_percluster.R` → `FigA1a_potency.png`, `FigA1b_local_sensitivity.png`.
- **Description** : (a) potency **cluster-invariante** : LAI −2.6 ≫ fCover −1.2 > H_max −0.25 > LAD +0.1. (b) sensibilité locale du LAI **s'effondre** ×4 : −4.2 (P1) → −1.0 (P4) ; idem fCover.
- **Interprétation** : la **puissance physique est universelle** (LAI toujours le plus potent) mais la **marge locale** varie → le trait dominant par archétype = potency × range disponible. Justifie l'analyse **par archétype**.

## Fig. A2 — Effet LAD vs écart à l'uniforme (VCI)
- **Données** : contraste per-plot dT_LAD (réel − uniforme) + VCI par plot (`tab_sensitivity_perplot_lad.csv`, `tab_lad_vci.csv`).
- **Process** : dT_LAD vs (1 − VCI), régressions par archétype, corrélation partielle | LAI. Script `c1_lad_vci.R` → `FigA2_lad_vci.png`.
- **Description** : intra-archétype, profil plus concentré (VCI bas) → contraste plus grand (**ρ partiel = −0.33** | LAI ; net en P1/P4). Le centre de masse, lui, est ~uniforme (0.51–0.54 H_max) et n'explique rien.
- **Interprétation** : l'effet du profil est une prise scalaire **via le VCI** (écart à l'uniforme), pas via le top/bottom ; le VCI reste **nécessaire mais pas suffisant** (colinéaire densité, un-à-plusieurs avec la forme) → d'où la typologie sur profil complet.

## Fig. B1 — Décomposition du résidu de validation
- **Données** : résidu (obs − modélisé ΔT_max, 53 HOBO) + covariables terrain/hétérogénéité (élévation, pente, TWI, northness, σ hauteur, rumple, **gap fraction**).
- **Process** : régression OLS ; Moran's I. `fig_residual_vs_topo.png`.
- **Description** : **adj R² = 0.62** ; dominé par **gap fraction** (r = 0.76 ; plots amplifiants gap 0.52 vs 0.01) ; résidu peu autocorrélé spatialement.
- **Interprétation** : ce que le modèle 1-D ne voit pas = **trouées sous-pixel** (hétérogénéité horizontale), pas la structure verticale → délimite le plafond de verre.

## Fig. C1 — Sensibilité au rayon d'empreinte
- **Données** : traits ré-extraits à 9 rayons (5–50 m) × 2 géométries (cercle/carré), runs HOBO (`radius_test_circle/square`).
- **Process** : R²/RMSE/MAE horaire vs HOBO. `FigAnnex_radius_sensitivity.png`.
- **Description** : R² **plat** 0.642–0.648, géométries indiscernables ; erreur minimale ~10–12.5 m.
- **Interprétation** : le résultat **ne dépend pas** du rayon/forme d'empreinte → choix 12.5 m justifié.

## Fig. D1 — Dispersion des traits par archétype + découplage hauteur–couvert
- **Données** : 400 plots cLHS (`clhs_sample_floor05_v2.rds`), LAI one-sided, fCover (planché à 0.5, = valeur passée à MuSICA), H_max ; labels P1–P4 (`relabel_cluster`).
- **Process** : violons + boxplots par trait × archétype (haut) ; nuage H_max vs fCover coloré par archétype avec corrélation intra-archétype annotée (bas). Script `c1_clhs_trait_dispersion.R` → `FigAnnex_clhs_trait_dispersion.png` (copie livrable `FigD1_trait_dispersion.png`).
- **Description** : SD de LAI quasi constante le long du gradient (0.62–0.76) ; fCover **sature** vers 1 dans le dense ; **SD de H_max de loin maximale en P1** (10.7 m vs 3.4–4.8 m), où H_max est **décorrélé du couvert** (r = 0.07) et couvre toute la gamme de hauteur (4–40 m) ; en P3/P4 hauteur et couvert sont **verrouillés ensemble**.
- **Interprétation** : le cLHS **n'impose aucune borne externe** — les étendues par trait sont les marginales empiriques intra-cluster, fidèlement reproduites. Le **découplage hauteur–densité au bout ouvert** (un couvert clair peut être jeune-bas OU mature-haut-dégradé) est une propriété de la forêt, pas un artefact ; il permet de sonder le levier hauteur **indépendamment de la densité** là où la densité le laisse libre. (Compagnon visuel de la Table D1.)

## Fig. E1–E2 — Structure de corrélation des traits (+ FPCs)
- **Données** : traits (LAI, fCover, H_max, VCI) et FPC1–3, sur 400 plots et par archétype.
- **Process** : matrices de corrélation. `FigE1_corrplot_traits.png`, `FigE2_corrplot_fpc.png`.
- **Description** : forte colinéarité globale (LAI-fCover 0.81, LAI-VCI 0.87) qui **s'effondre en P4** (LAI-fCover 0.17) ; FPCs faiblement/moyennement corrélés aux scalaires.
- **Interprétation** : motive l'analyse de sensibilité par perturbation (vs régression) ; la colinéarité densité-dépendante = le miroir de la saturation ; le profil (FPC) est un axe largement indépendant.

---

# Supplémentaires — validation de la typologie à 1 m (lien écologique, §4.4)

## Fig. S1 — ACP récupère les 4 archétypes
- **Données** : traits scalaires (LAI, fCover, H_max) + VCI sur les 400 plots, labels de cluster.
- **Process** : ACP non supervisée ; projection des plots, coloriés par archétype. Script `c1_pca_recover_clusters.R` → `FigS1_pca_clusters.png`.
- **Description** : les 4 clusters **se séparent** le long de PC1 (~75% de variance = axe densité/fermeture portant LAI+fCover+VCI) ; PC2 = H_max. VCI **colinéaire** à l'axe densité (un scalaire ne crée pas d'axe séparateur indépendant ; P3/P4 se chevauchent).
- **Interprétation** : la typologie n'est pas arbitraire — une ACP indépendante des traits scalaires la **retrouve**, à la résolution **1 m**. Et le chevauchement P3/P4 = miroir descriptif de la saturation (justifie le profil complet vs un scalaire).

## Fig. S2 — Le microclimat mesuré se retrouve dans le clustering
- **Données** : ΔT_max et pente **observés** (53 HOBO), groupés par archétype LiDAR.
- **Process** : distribution observée par cluster + Kruskal–Wallis. Script `c1_hobo_direction.R` → `FigS2_microclimate_by_cluster.png`.
- **Description** : le ΔT_max observé **ordonne les types** — amplifiant dans l'ouvert → tampon dans le dense ; **KW p ≈ 7 × 10⁻⁵** (ΔT_max), ≈ 3 × 10⁻⁷ (pente) *(clusters canoniques `clusters.rds`, après correction du bug d'unités logger one-sided vs cLHS two-sided)*.
- **Interprétation** : **validation forte et indépendante** — la typologie structurale (profils LiDAR 1 m) **prédit le régime thermique réellement mesuré** au sol. Son sens écologique est établi **à 1 m**, sans recourir au 0.5 m. *(Support direct du §4.4 : les 4 profils sont réels, pas une abstraction.)*

---

## Reproductibilité (figure → script → données)
| Fig | Script | Données |
|---|---|---|
| 1 | make_fig_archetypes_profiles_v3.R | LAD 1 m + clusters |
| 2 | c1_importance_sensitivity.R | nc_sensitivity_perplot (+1111/1110) |
| 3 | c1_sensitivity_perplot_figure.R | idem + tab_sensitivity_perplot_lad |
| 4 | pipeline 03z/07 | musica_hobo_z05 + HOBO obs |
| 5 | c1_forward_percluster.R | musica_hobo_z05 + tab_importance_sensitivity |
| 6 | (pipeline annex) | nc 15 couches |
| A1 | c1_sensitivity_percluster.R | nc per-unit |
| A2 | c1_lad_vci.R | tab_sensitivity_perplot_lad + VCI |
| B1 | make_residual_vs_topo.R | résidu + covariables terrain |
| C1 | c1_radius_*.R | radius_test_circle/square |
| D1 | c1_clhs_trait_dispersion.R | clhs_sample_floor05_v2 |
| E1/E2 | make_corrplot_* | traits + FPC |
