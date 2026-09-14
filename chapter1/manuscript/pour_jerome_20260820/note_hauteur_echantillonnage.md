# Lire le modèle près du sol : ce que ça change

Demande de Jérôme, 20 août 2026 : refaire les figures à 0.1 m plutôt qu'à 1 m, les HOBO étant
contre le tronc et donc à l'abri des masses d'air, peut-être plus proches d'une température de
sol. Script : `scripts/c1_sampling_height_test.R` et `scripts/c1_logslope_sampling_height.R`.
Aucune simulation relancée.

## 0.1 m n'existe pas dans le modèle

La grille d'air de MuSICA n'a aucun nœud sous 0.0285 × Hmax. Sur les 53 loggers ça donne
0.085 m dans le peuplement le plus court et 1.083 m dans le plus haut, médiane 0.769 m.
`micro_hourly_at()` **plafonne** au nœud 1 au lieu d'extrapoler, donc demander 0.1 m renvoie
le nœud 1 sur toutes les placettes de plus de 3.5 m de haut. La seule expérience bien définie
est donc `nair == 1`, le nœud le plus bas.

À noter au passage : sur 7 placettes (Hmax 36 à 38 m) le nœud 1 est déjà à 1.026 - 1.083 m,
donc la lecture « à 1 m » y renvoie en fait 1.03 à 1.08 m. Écart maximal de 8 cm, sans
conséquence, mais « toutes les placettes à 1 m fixe » est une approximation à 8 cm près et
non une exactitude.

## Le piège : ce n'est pas un décalage uniforme

L'hypothèse de Jérôme est un décalage vers le bas identique partout, les loggers étant tous à
1 m. `nair == 1` fait autre chose : il descend à 0.26 m en moyenne dans le tiers le plus
ouvert et reste à 1.00 m dans le tiers le plus haut. La manipulation est donc confondue avec
la hauteur de canopée, qui est l'un des axes de la figure. Tout ce qui suit est découpé par
tercile de Hmax pour que ce confondant reste visible.

## Résultats

Réseau entier, n = 53 :

| Lecture | ΔTmax r | biais | amplitude | pente r | pente moyenne |
|---|---|---|---|---|---|
| 1 m | 0.501 | +0.98 °C | 15 % | 0.782 | 0.913 |
| `nair == 1` | 0.577 | +0.98 °C | 18 % | 0.847 | 0.898 |

Par tercile de hauteur :

| Tercile | n | Hmax moyen | nœud bas | biais 1 m | biais bas | ΔTmax r 1 m → bas | pente r 1 m → bas |
|---|---|---|---|---|---|---|---|
| court | 18 | 8.8 m | 0.26 m | +0.02 °C | +0.04 °C | 0.664 → 0.757 | 0.835 → 0.837 |
| moyen | 17 | 24.7 m | 0.71 m | +1.34 °C | +1.32 °C | 0.727 → 0.756 | 0.614 → 0.710 |
| haut | 18 | 34.7 m | 1.00 m | +1.60 °C | +1.59 °C | 0.700 → 0.701 | 0.716 → 0.718 |

Déplacement par placette : ΔTmax bouge de 0.037 °C en moyenne, 0.19 °C au maximum ; la pente
de 0.017 en moyenne, 0.083 au maximum. L'étendue de log(pente) passe de -0.330..+0.015 à
-0.330..+0.025, et les mêmes 8 placettes amplifient.

## Interprétation

**Le biais chaud ne bouge pas.** +0.98 °C dans les deux cas, et le déplacement par tercile est
de ±0.02 °C. Dans la lignée actuelle MuSICA n'a quasiment pas de gradient thermique sous 1 m.
L'hypothèse du tronc, prise comme un simple effet de hauteur, n'explique donc pas le biais.

**Et surtout, elle n'est pas testable là où le biais vit.** Le biais est nul dans le tiers
ouvert (+0.02 °C) et vaut +1.34 et +1.60 °C dans les tiers moyen et haut. Or c'est précisément
dans ces peuplements-là que le nœud le plus bas est déjà à 1.00 m : on ne peut pas y descendre
le capteur virtuel. Le biais chaud du modèle est un problème de canopée fermée, pas de canopée
ouverte.

**Les corrélations s'améliorent, mais exactement là où le nœud descend.** ΔTmax r passe de
0.664 à 0.757 dans le tiers court, où le nœud descend à 0.26 m, et de 0.700 à 0.701 dans le
tiers haut, où il ne descend pas. Le gain global de 0.501 à 0.577 est donc porté par les
placettes ouvertes, ce qui est cohérent avec un effet de hauteur, mais aussi indissociable du
confondant Hmax. Je ne le présenterais pas comme un argument en faveur de l'hypothèse.
