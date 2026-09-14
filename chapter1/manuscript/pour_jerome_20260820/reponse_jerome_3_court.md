Version courte du 3e échange (détails : note_hauteur_echantillonnage.md, note_vci_rayon_zmax.md).
Pièce jointe : Fig_logslope_sampling_height.png

---

Bonjour Jérôme,

J'ai fait les deux tests.

**La hauteur.** Un 0.1 m fixe n'existe pas dans le modèle : la grille d'air de MuSICA n'a
aucun nœud sous 0.0285 × Hmax, soit 0.09 m dans le peuplement le plus court mais déjà 1.08 m
dans le plus haut. J'ai donc lu le nœud le plus bas, c'est la figure jointe. Le biais chaud ne
bouge pas d'un iota, +0.98 °C dans les deux cas : sous 1 m le modèle n'a quasiment pas de
gradient thermique. Et ton hypothèse n'est pas testable là où le problème se trouve, puisque
le biais est nul en canopée ouverte (+0.02 °C) et vaut +1.3 à +1.6 °C en canopée fermée, là
où justement le nœud le plus bas est déjà à 1 m. Si tu veux un vrai 0.1 m partout, il faut
doubler la résolution verticale du modèle et relancer ; dis-moi si ça t'intéresse.

**Le VCI.** Même indice et même fonction que toi, van Ewijk 2011 via `lidR::VCI`. La seule
différence est le `zmax` passé en argument, qui fixe le nombre de tranches. Je prends le
maximum de hauteur de la placette elle-même, si bien qu'un peuplement de 3 m a 3 tranches et
peut très bien être « complexe » ; avec un plafond fixe à l'échelle du massif il aurait
40 tranches dont 3 occupées et tomberait vers 0. C'est ça qui fait tes points près de 0, pas
le rayon : à 5 m, ton rayon, avec ma normalisation, j'obtiens exactement la même distribution
qu'à 10 ou 25 m.

Je garde ma version, pour une raison précise : avec un plafond fixe le VCI corrèle à 0.89 avec
Hmax, il devient une métrique de hauteur déguisée. Avec la mienne il tombe à 0.50 et mesure
vraiment l'homogénéité du remplissage. Quitte à ne pas être superposable à ta Fig. 4.

Sur le rayon, je suis dans ta borne : les simulations tournent sur un disque de 10 m, et le
balayage 5 à 50 m est déjà dans le chapitre.

Bien à toi,

Nathan
