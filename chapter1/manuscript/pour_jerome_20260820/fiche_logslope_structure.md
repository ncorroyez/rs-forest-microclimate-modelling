# Fiche — log(pente micro-macro) vs métriques structurales, avant / après l'été

Demande de Jérôme, 20 août 2026 : retrouver la figure log(slope) vs Hmax, VCI et PAI
(observé et modélisé), celle du poster SilviLaser, pour voir comment les résultats ont
évolué entre les simulations d'avant et d'après l'été.

Script : `scripts/c1_logslope_vs_structure.R` (un seul `Rscript`, aucune simulation MuSICA).

## Méthodes

Pente = coefficient d'une régression OLS de la température horaire sous canopée à 1 m
sur la température horaire du macroclimat, une régression par logger, JJAS 2021.
C'est exactement la pente de validation du chapitre (`c1_hobo_native20_validation_regen.R`),
donc la métrique du panneau (b) de la Figure 6 du manuscrit. `log(pente) < 0` = tamponnement,
`> 0` = amplification. Les 53 pentes sont positives, le log est défini partout.

Chaque lignée de simulation est régressée sur **son propre** forçage macro : `musica_in_Blois.nc`
pour les runs d'avant l'été (`out_files/musica_hobo_ERA5_results/`, 31 mars 2026),
`in_files/FR-Blo_2021_v2.nc` pour les runs actuels (`out_files/musica_hobo_native20/1111/`,
17 juillet 2026). C'est ce qui rend chaque série cohérente en interne. L'axe des x est en
revanche le **tableau de traits actuel** dans les deux cas, de sorte que les deux courbes
modèle reposent sur la même abscisse et que seule la lignée de simulation bouge.

Ajustement quadratique par panneau et par série (`y ~ poly(x, 2)`), comme sur le poster.

Les figures ne portent que la légende des courbes, avec la version du modèle dans le nom de chaque
série. Le détail des corrections reste ici :

- **Avant l'été** (run du 31 mars 2026) : MuSICA v3.2.0, forçage `musica_in_Blois.nc`, pas de
  correction de vent, traits antérieurs à la refonte native 20 m.
- **Après l'été** (run du 17 juillet 2026) : MuSICA v3.2.3, couplage ABL `abl_flag = 'iter'`,
  forçage `FR-Blo_2021_v2.nc` avec correction de vent en profil logarithmique par placette,
  LAI corrigé de l'angle de visée, traits calculés nativement sur carrés de 20 m,
  couvert fractionnel plancher à 0.5.

Les numéros de version viennent de l'attribut `version` des NetCDF de sortie, pas d'une note ;
le `abl_flag` et le forçage de la lignée actuelle viennent de `c1_hobo_native20_validation.R:13`
et du dump `tab_musica_config_dump.csv` (colonne AS-RUN).

## Pré-traitements

- Traits recalculés nativement à 20 m (lignée native20), pas d'agrégation 10 m → 20 m.
- VCI recalculé sur hauteurs normalisées au sol. C'est pourquoi il court ici de 0.05 à 0.93
  alors que le poster affichait 0 à 0.73 : l'ancien raster n'était pas normalisé.
- Le troisième panneau porte le **LAI une face** utilisé dans le chapitre (0.02 à 6.71),
  pas le PAI brut du poster qui montait au-delà de 12. Seule la forme est comparable,
  pas l'abscisse.
- Correction de vent en profil logarithmique intégrée au forçage des runs actuels ;
  MuSICA v3.2.3 avec `abl_flag = iter`.

## Données et figures

| Fichier | Contenu |
|---|---|
| `Fig_logslope_vs_structure_native20.png` (+ pdf) | Deux séries, capteurs de terrain et MuSICA actuel. L'analogue direct du panneau (c) du poster. |
| `Fig_logslope_vs_structure_evolution.png` (+ pdf) | Trois séries, avec les simulations d'avant l'été en violet. |
| `Fig_logslope_AVANT_ETE_2026-04.png` | La figure d'avant l'été telle qu'elle avait été produite (avril 2026, `Plots/Silvilaser/compfull.png`), avec ses propres traits et son propre VCI. |
| `Poster_SilviLaser_2025.png` | Le poster d'origine, pour mémoire. |
| `out_files/Chapter1/tables/tab_logslope_vs_structure.csv` | Table par logger : pentes obs, avant, après, plus LAI, Hmax, VCI. |
| `out_files/Chapter1/tables/tab_logslope_presummer.csv` | Cache des pentes d'avant l'été (évite de relire les 53 NetCDF). |

## Description

Amplitude de log(pente) sur les 53 loggers :

| Série | Étendue | Part de l'étendue observée |
|---|---|---|
| Capteurs de terrain | -0.363 à +0.312 (0.675) | — |
| MuSICA avant l'été | -0.048 à +0.112 (0.160) | 24 % |
| MuSICA après l'été | -0.330 à +0.015 (0.345) | 51 % |

Corrélation avec la pente observée : 0.934 avant, 0.782 après.
Loggers à pente > 1 : 8 observés, 31 avant l'été, 8 après.

Ajustements quadratiques, R² (Spearman entre parenthèses) :

| Panneau | Terrain | Avant | Après |
|---|---|---|---|
| Hmax | 0.623 (-0.81) | 0.666 (-0.78) | 0.727 (-0.90) |
| VCI | 0.681 (-0.74) | 0.675 (-0.69) | 0.391 (-0.72) |
| LAI | 0.770 (-0.89) | 0.843 (-0.88) | 0.750 (-0.84) |

## Interprétation

Le gradient structural du modèle a plus que doublé. Avant l'été la courbe MuSICA était
quasi plate et collée à zéro : le modèle tamponnait à peine, et il amplifiait sur 31 des
53 loggers alors que le terrain n'en amplifie que 8. Après l'été la courbe descend
franchement avec la hauteur et avec le LAI, elle couvre la moitié de l'étendue observée
au lieu du quart, et le partage tamponnement / amplification tombe juste, 8 contre 8.

Le plafond de verre reste entier, et se lit directement sur la figure. Le modèle plafonne
à une pente de 1.02 quand le terrain monte à 1.37 : il n'amplifie toujours pas vraiment.
Sur les placettes ouvertes, à Hmax faible ou LAI faible, les points verts montent bien
au-dessus des points orange.

Un point à ne pas cacher : la corrélation avec la pente observée **baisse**, de 0.93 à 0.78.
L'ancienne lignée classait mieux les placettes tout en écrasant massivement l'amplitude,
un classement presque parfait sur un intervalle quatre fois trop étroit. La nouvelle
restitue l'amplitude au prix d'un peu de rang. C'est le compromis habituel entre justesse
et précision, et c'est la restitution d'amplitude qui compte pour l'argument du chapitre,
puisque c'est elle qui porte le ΔTmax.

Le panneau VCI est celui qui bouge le plus dans sa forme, et c'est un artefact de mesure
plus qu'un résultat : la normalisation des hauteurs a tassé les placettes denses entre
0.80 et 0.93, si bien que l'ajustement quadratique y travaille sur un nuage très concentré
à droite et une seule placette à 0.05. Le R² du modèle y chute à 0.39. Ne pas lire ce
panneau comme une perte de performance.

## Provenance des axes, et écart avec le poster (question de Jérôme, 20 août 2026)

**Footprint des métriques.** Poster : disque de rayon 25 m, `out_files/radius_test/metrics_results.csv`
filtré sur `radius == 25` (`main_compare_era5_ign.R:98`). Aujourd'hui : valeur de la maille
native 20 m contenant le logger, rasters `in_files_native20/` extraits au point (vérifié :
écart maximal de 5e-15 sur LAI, Hmax, fCover et VCI contre `tab_hobo_perplot_cluster.csv`).
400 m² contre environ 1960 m².

| Métrique | Poster, disque r = 25 m | Maille native 20 m |
|---|---|---|
| Hmax | 6.5 - 43.5, médiane 29.5 | 2.88 - 37.56, médiane 26.72 |
| VCI | 0.19 - 0.66, médiane 0.58 | 0.052 - 0.929, médiane 0.846 |
| PAI / LAI | 0.02 - 12.65, médiane 8.04 | 0.022 - 6.71, médiane 3.93 |

**Hauteur d'extraction du modèle.** Le poster échantillonnait `nair == 1`
(`main_compare_era5_ign.R:125`), le premier niveau d'air. `relative_height[1]` vaut 0.02849
dans tous les fichiers (vérifié sur les 53), soit 0.07 m sous une canopée de 2.5 m et 1.07 m
sous 37.5 m. A/B sur les mêmes NetCDF de mars, avec exactement la convention de la série
violette (macro `musica_in_Blois.nc`, JJAS 2021), en ne changeant que la hauteur :
`nair == 1` donne log(pente) de -0.067 à +0.291, 1 m fixe interpolé donne -0.048 à +0.112
(cette seconde ligne est la série violette de la figure). Quatre placettes dépassent +0.2
sous `nair == 1`, toutes à Hmax < 7.1 m. Les log(pente) très positives du poster en milieu
ouvert sont donc un effet de hauteur, pas un effet de version.

**Incohérence du poster.** Les simulations de mars ont été paramétrées avec un footprint
proche de r = 5 à 10 m : leur `veget_height_top` s'écarte au maximum de 2 à 3 m du balayage
r = 5 / r = 10, contre 27 m pour r = 25, alors que le poster les trace contre les métriques
r = 25. `veget_height_top` est quantifié au demi-mètre, donc r = 5 et r = 10 ne se
départagent pas. La maille native 20 m est plus proche mais pas identique : écart maximal de
9.46 m, écart moyen de 1.00 m. L'abscisse de la figure n'est donc pas celle avec laquelle la
série violette a tourné, seulement une bien meilleure approximation que le r = 25.

**Macro observé.** `musica_in_Blois.nc` (poster) contre `FR-Blo_2021_v2.nc` (maintenant),
r = 0.92 entre les deux séries horaires, écart moyen de 0.96 °C. L'étendue de log(pente)
observée passe de 0.740 à 0.675, le maximum de +0.365 à +0.312.
