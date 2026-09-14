Réponse à Jérôme, 20 août 2026 (3e échange) : test de la hauteur d'échantillonnage (0.1 m) et
calcul du VCI (rayon et normalisation).
Pièce jointe : Fig_logslope_sampling_height.png
Détails : note_hauteur_echantillonnage.md, note_vci_rayon_zmax.md

---

Bonjour Jérôme,

J'ai fait les deux tests.

**La hauteur d'échantillonnage.** Un 0.1 m fixe n'existe pas dans le modèle : la grille d'air
de MuSICA n'a aucun nœud sous 0.0285 × Hmax, soit 0.085 m dans le peuplement le plus court
mais 1.083 m dans le plus haut, médiane 0.77 m. Le mieux que je puisse faire est de lire le
nœud le plus bas, ce que fait la figure jointe. Attention à la lecture : ce n'est pas un
décalage uniforme, on descend à 0.26 m en moyenne dans le tiers ouvert et on reste à 1.00 m
dans le tiers haut. La manipulation est donc confondue avec la hauteur de canopée.

Résultat : **le biais chaud ne bouge pas du tout**, +0.98 °C dans les deux cas, et le
déplacement du ΔTmax par placette est de 0.04 °C en moyenne. Dans la lignée actuelle MuSICA
n'a quasiment pas de gradient thermique sous 1 m.

Et surtout, ton hypothèse n'est pas testable là où le problème se trouve. Le biais est nul
dans le tiers ouvert (+0.02 °C) et vaut +1.34 et +1.60 °C dans les tiers moyen et haut. Or
c'est exactement dans ces peuplements que le nœud le plus bas est déjà à 1.00 m : impossible
d'y descendre le capteur virtuel. Le biais chaud est un problème de canopée fermée.

Les corrélations, elles, s'améliorent : ΔTmax r de 0.501 à 0.577, pente r de 0.782 à 0.847.
Mais le gain est concentré là où le nœud descend (tiers court, r de 0.664 à 0.757) et nul là
où il ne descend pas (tiers haut, 0.700 à 0.701), donc je ne le lirais pas comme un appui à
l'hypothèse. Si tu veux un vrai 0.1 m partout, il faudrait doubler la résolution verticale
(`AIR_RESOLUTION_LEVEL = 1`, 20 couches de végétation au lieu de 10) et relancer. Dis-moi si
ça t'intéresse, un run test suffit à voir où tombe le nouveau nœud bas.

**Le VCI.** Même indice, van Ewijk 2011, et la même fonction `lidR::VCI`. La différence est
l'argument `zmax`, qui fixe le nombre de tranches et donc le dénominateur. Je passe le maximum
de hauteur de la placette elle-même : un peuplement ouvert de 3 m a 3 tranches, ses retours
s'y répartissent, il atterrit en milieu de gamme. Avec un plafond fixe à l'échelle du massif
il aurait 40 tranches dont 3 occupées, et son VCI s'effondrerait vers 0. C'est ça qui explique
tes points près de 0, pas le rayon.

Je l'ai vérifié en recalculant aux 53 loggers. À 5 m, ton rayon, avec ma normalisation :
0.072 à 0.935, une seule placette sous 0.2, c'est-à-dire la même distribution qu'à 10 m, 25 m
ou sur la maille de 20 m. Avec un plafond fixe à 40.5 m et le même rayon de 5 m : 0.027 à
0.866, quatre placettes sous 0.2 et deux sous 0.1. Le plafond fait le plancher, le rayon non.

Un argument pour garder ma version, cela dit. Avec un plafond fixe, le VCI corrèle à 0.89 avec
Hmax : un peuplement court ne peut structurellement pas être réparti sur 40 tranches, donc
l'indice devient une métrique de hauteur déguisée. Il prédit d'ailleurs mieux la pente
observée (r = -0.89 contre -0.69), mais en portant Hmax avec lui. Avec ma normalisation la
corrélation avec Hmax tombe à 0.50 et l'indice mesure ce qu'on lui demande, l'homogénéité du
remplissage à l'intérieur de la gamme de hauteur du peuplement. Je préfère garder ça, quitte
à ne pas être directement superposable à ta Fig. 4.

Sur le rayon enfin, je suis dans ta borne : les simulations tournent sur un disque de 10 m,
et le balayage 5 à 50 m est déjà dans le chapitre (Fig. 8, 420 simulations).

Bien à toi,

Nathan
