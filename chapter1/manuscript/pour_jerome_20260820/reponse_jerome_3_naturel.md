Version naturelle du 3e échange. Pièce jointe : Fig_logslope_sampling_height.png

---

Bonjour Jérôme,

J'ai regardé les deux, et le premier point m'a appris quelque chose.

Pour la hauteur, je ne peux pas descendre à 0.1 m : la grille d'air de MuSICA n'a pas de nœud
sous 0.0285 × Hmax, ce qui fait 9 cm dans le peuplement le plus court mais déjà 1.08 m dans le
plus haut. J'ai donc pris le nœud le plus bas, c'est la figure jointe. Le biais chaud ne bouge
pas, +0.98 °C avant comme après, tout simplement parce que le modèle n'a presque pas de
gradient thermique sous 1 m. Le plus embêtant, c'est qu'on ne peut pas vraiment tester ton
idée là où elle compterait : le biais est nul en canopée ouverte et monte à +1.3 ou +1.6 °C en
canopée fermée, or c'est justement là que le nœud le plus bas est déjà à 1 m. Il faudrait
doubler la résolution verticale du modèle et tout relancer pour aller plus bas, dis-moi si tu
penses que ça vaut le coup.

Pour le VCI, on utilise bien le même, van Ewijk et la fonction de lidR. La seule différence
est le zmax qu'on passe en argument, qui fixe le nombre de tranches. Je prends le maximum de
hauteur de la placette elle-même, donc un peuplement de 3 m a 3 tranches et peut très bien
être réparti dessus ; avec un plafond fixe à l'échelle du massif, il aurait 40 tranches dont
3 occupées et s'effondrerait vers 0. C'est ça qui fait tes points près de zéro. J'ai vérifié
en recalculant à ton rayon de 5 m et je retrouve exactement la même distribution qu'à 10 ou
25 m, donc le rayon n'y est pour rien.

Je préfère garder ma normalisation, pour une raison assez nette : avec le plafond fixe, le VCI
corrèle à 0.89 avec la hauteur maximale, il devient une métrique de hauteur déguisée, alors
qu'avec la mienne il tombe à 0.50 et mesure ce qu'on lui demande. Quitte à ne pas être
superposable à ta figure 4.

Et sur le rayon lui-même je reste dans ta borne, les simulations tournent sur un disque de
10 m, avec le balayage 5 à 50 m déjà dans le chapitre.

Bien à toi,

Nathan
