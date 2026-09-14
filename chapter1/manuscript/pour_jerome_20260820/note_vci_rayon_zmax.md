# VCI : pourquoi pas de points près de 0, contrairement à la Fig. 4 de Gril et al. 2023

Question de Jérôme, 20 août 2026. Script : `scripts/c1_vci_radius_zmax_test.R`
(clip LiDAR aux 53 loggers, 3 rayons, 2 normalisations, sans le garde-fou `zmax < 2` de
`c1_vci_tile.R`. Sans effet en pratique : `lidR::entropy` applique de toute façon
`if (zmax < 2*by) return(NA)`, et aucune des 53 placettes n'a un max sous 2 m, donc les
53 valeurs sont finies dans les six combinaisons.)

## Même indice, même fonction, une seule différence

Elle et moi utilisons le même van Ewijk et al. 2011 via `lidR::VCI`, entropie de Shannon des
retours par tranches de 1 m, normalisée. Le code de `lidR` est explicite : `VCI(z, zmax, by)`
écarte d'abord `z >= zmax`, puis `entropy()` découpe `seq(0, ceiling(zmax/by)*by, by)`, compte
l'occupation de chaque tranche et divise l'entropie par celle du cas uniforme, soit
log(ceiling(zmax)). La seule différence entre nos deux calculs est donc l'argument `zmax`, qui
fixe le nombre de tranches et donc ce dénominateur.

Deux détails du même code, à garder en tête : les retours sol sont inclus (on filtre `z >= 0`),
donc la tranche [0, 1[ pèse lourd et tire le VCI vers le bas ; et `zmax = max(z)` fait que le
filtre `z < zmax` retire exactement le point le plus haut.

Je passe **le max de hauteur de la placette elle-même**. Un peuplement ouvert de 3 m a donc
3 tranches : ses retours peuvent s'y répartir, et il atterrit en milieu de gamme. Avec un
**plafond fixe à l'échelle du massif**, la même placette a 40 tranches dont 3 occupées,
l'entropie est minuscule devant log(40), et le VCI s'effondre vers 0. L'article ne dit pas
quel `zmax` a été passé ; les 40.5 m qu'il mentionne sont le plafond des tranches de PAD.

## Le rayon n'est pas la cause

| Rayon | VCI own-max min / médiane / max | < 0.2 |
|---|---|---|
| 5 m (celui d'Eva) | 0.072 / 0.829 / 0.935 | 1 |
| 10 m | 0.047 / 0.847 / 0.924 | 1 |
| 25 m | 0.127 / 0.859 / 0.915 | 1 |
| maille 20 m (chapitre) | 0.052 / 0.846 / 0.929 | 1 |

Passer au rayon de 5 m ne crée aucun point près de 0. La distribution est la même partout.

## Le plafond de normalisation, si

| Rayon | VCI plafond fixe 40.5 m, min / médiane / max | < 0.2 | < 0.1 |
|---|---|---|---|
| 5 m | 0.027 / 0.708 / 0.866 | 4 | 2 |
| 10 m | 0.018 / 0.730 / 0.872 | 3 | 3 |
| 25 m | 0.061 / 0.755 / 0.903 | 3 | 1 |

Le plancher descend et quelques points près de 0 apparaissent, la signature que Jérôme
cherche. La médiane descend aussi, de 0.83 à 0.71. On ne reproduit pas pour autant toute
l'étendue de sa Fig. 4 : son maximum est vers 0.65, le nôtre reste à 0.87. Il reste donc
d'autres choix en jeu (largeur de tranche, inclusion ou non des retours sol, seuil de
hauteur) que l'article ne précise pas.

## Et un argument pour garder notre choix

| Rayon | r(VCI, Hmax) own-max | r(VCI, Hmax) plafond fixe |
|---|---|---|
| 5 m | 0.504 | 0.891 |
| 10 m | 0.606 | 0.906 |
| 25 m | 0.637 | 0.887 |

Avec un plafond fixe, le VCI est à 0.89 - 0.91 de corrélation avec la hauteur maximale : c'est
une métrique de hauteur déguisée, puisqu'un peuplement court ne peut structurellement pas être
« réparti » sur 40 tranches. Il prédit d'ailleurs mieux la pente observée pour cette raison
(r = -0.89 contre -0.69 à 5 m), mais en portant Hmax avec lui. Notre normalisation par le max
de la placette mesure ce qu'on veut vraiment, l'homogénéité du remplissage à l'intérieur de la
gamme de hauteur du peuplement, et reste bien moins confondue (r = 0.50 - 0.64).

## Le rayon, pour mémoire

Les simulations actuelles ont tourné sur un disque de rayon 10 m, donc dans la borne
recommandée par Gril et al. (≤ 10 m). L'axe des x des figures porte la maille native 20 m,
écart maximal de 2.1 m sur Hmax, testé sans effet. Et le balayage de rayons 5 à 50 m
(420 simulations) est déjà dans le manuscrit, Fig. 8.
