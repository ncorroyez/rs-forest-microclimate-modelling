# Ce qui a changé entre le poster SilviLaser et la figure d'aujourd'hui

Tout ce qui suit est vérifié sur les fichiers (attributs NetCDF, scripts, rasters), pas de mémoire.

## Tableau

| | Avant (poster, sims du 31 mars 2026) | Après (chapitre, sims du 17 juillet 2026) |
|---|---|---|
| Modèle | MuSICA v3.2.0 | MuSICA v3.2.3, `abl_flag = 'iter'` |
| Sorties | `out_files/musica_hobo_ERA5_results/` | `out_files/musica_hobo_native20/1111/` |
| Forçage | `in_files/musica_in_Blois.nc` | `in_files/FR-Blo_2021_v2.nc` |
| Vent | aucune correction | profil log de canopée par placette, intégré au forçage |
| LAI | brut | corrigé de l'angle de visée (sec θ) |
| Traits qui ont **paramétré** les sims | valeur du **pixel 10 m** sous le logger | **MAX sur un disque de rayon 10 m** des rasters 10 m |
| Traits tracés en **x** | disque de **rayon 25 m** (`radius_test`, `radius == 25`) | **maille native 20 m** (`in_files_native20/`) |
| VCI | raster 10 m, élévation absolue, 0.19-0.66 | raster natif 20 m, hauteurs normalisées au sol, 0.05-0.93 |
| 3e axe | PAI deux faces, 0.02-12.65 | LAI une face, 0.02-6.71 |
| Hauteur d'extraction | `nair == 1` (hauteur relative 0.0285, soit 0.07 à 1.07 m) | 1 m fixe interpolé, la hauteur des HOBO |
| Période | 1 juin - 23 sept 2021 | 1 juin - 30 sept 2021 |
| Macro des pentes | `musica_in_Blois.nc` | `FR-Blo_2021_v2.nc` |

## Comment le footprint des sims a été identifié

`veget_height_top` des NetCDF, comparé aux candidats d'extraction sur les 53 loggers
(le résidu de 0.5 à 1 m est l'arrondi au demi-mètre que MuSICA applique) :

| Candidat | écart max / moyen contre les sims d'avant | contre les sims d'après |
|---|---|---|
| pixel 10 m | **1.11 / 0.29 m** | 9.82 / 1.31 m |
| maille 20 m | 9.46 / 1.15 m | 2.10 / 0.73 m |
| MAX disque r = 10 m (rasters 10 m) | 9.49 / 0.97 m | **0.99 / 0.49 m** |
| MAX disque r = 25 m | 26.74 / 3.17 m | - |

## Deux décalages entre ce qui a tourné et ce qui est tracé

**Le poster : 27 m.** Les sims ont tourné sur le pixel 10 m, la figure les trace contre le
disque de 25 m. Vingt fois la surface. C'est ce qui explique l'essentiel du déplacement
horizontal des points entre les deux figures.

**Ma figure : 2.1 m.** Les sims d'après ont tourné sur le MAX d'un disque de 10 m, l'axe des
x porte la maille native 20 m, qui est la table de traits canonique du chapitre. L'écart est
petit mais pas nul, et il concerne aussi les autres figures du chapitre qui croisent une
sortie MuSICA avec cette table.

## Effets sur le résultat

- Étendue de log(pente) du modèle : 0.160 avant, 0.345 après, contre 0.675 observés.
- Placettes à pente > 1 : 31 avant, 8 après, 8 observées.
- Corrélation avec la pente observée : 0.934 avant, 0.782 après.
- Les log(pente) > 0.2 du poster en canopée ouverte : effet de `nair == 1`. Sur les mêmes
  sims de mars, en ne changeant que la hauteur, +0.291 à `nair == 1` contre +0.112 à 1 m fixe,
  et les 4 placettes au-dessus de +0.2 sont toutes à Hmax < 7.1 m.

## Variante : chaque série sur les traits qui l'ont paramétrée

`scripts/c1_logslope_vs_structure_ownx.R` -> `Fig_logslope_vs_structure_ownx.png`.
Avant l'été sur le pixel 10 m, après l'été sur le MAX en disque de 10 m avec le LAI corrigé
de l'angle de visée (facteur global 1.0513). Les capteurs de terrain, qui ne sont
paramétrés par rien, sont tracés sur le footprint d'après l'été. Le VCI n'est un input
d'aucune des deux lignées, seul son footprint est aligné, sur le même raster normalisé.

Seul l'axe des x bouge, les pentes sont inchangées, donc toutes les statistiques de la note
ci-dessus (étendues, corrélation avec l'observé, nombre de placettes amplificatrices)
restent identiques. Ce qui bouge, ce sont les ajustements quadratiques :

| Panneau | Série | R² abscisse canonique | R² abscisse propre |
|---|---|---|---|
| Hmax | terrain | 0.623 | 0.646 |
| Hmax | avant | 0.666 | 0.742 |
| Hmax | après | 0.727 | 0.726 |
| VCI | terrain | 0.681 | 0.638 |
| VCI | avant | 0.675 | 0.496 |
| VCI | après | 0.391 | 0.329 |
| LAI | terrain | 0.770 | 0.863 |
| LAI | avant | 0.843 | 0.908 |
| LAI | après | 0.750 | 0.822 |

Hauteur et LAI se resserrent (le footprint colle enfin à ce qui a tourné), le VCI se
dégrade, ce qui est attendu puisqu'il n'est pas un input et que son alignement est
conventionnel. Conclusion : le récit avant/après ne dépend pas du choix d'abscisse.
