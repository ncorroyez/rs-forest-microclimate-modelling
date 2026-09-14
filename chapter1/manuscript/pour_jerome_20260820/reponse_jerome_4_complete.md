Réponse complète et courte : hauteur d'échantillonnage (grille doublée) + VCI (rayon et zmax).
Pièce jointe : Fig_logslope_airres1.png

---

Bonjour Jérôme,

Ton idée sur la hauteur m'a poussé plus loin que prévu.

Je ne pouvais pas descendre sous 1.08 m dans les grands peuplements parce que la grille d'air
de MuSICA n'a que 10 couches. Je l'ai passée à 20, c'est juste un paramètre de namelist, et le
nœud le plus bas tombe à 0.26 m. Pas 0.1 m, il faudrait 28 couches et le modèle n'en propose
que 10 ou 20, mais assez pour tester en canopée fermée, ce qui était impossible jusqu'ici.

Sur le biais tu as partiellement raison : lire à 0.2 m au lieu de 1 m en enlève 0.16 à
0.20 °C en canopée fermée. Réel, mais ça ne fait que 12 % d'un biais de 1.6 °C.

Le vrai gain est ailleurs, et je ne l'attendais pas. La corrélation sur le ΔTmax passe de 0.50
à 0.72, celle sur la pente de 0.78 à 0.88, et l'amplitude restituée de 15 à 22 %. Sur la
figure jointe, le modèle lu près du sol suit bien mieux le nuage des capteurs. Le plafond de
verre tient, 22 % reste loin de 100 %, mais il se fissure. Réserve quand même : lire le modèle
à 0.2 m face à des capteurs posés à 1 m ne se justifie que si ton mécanisme est le bon, le
gain pourrait aussi venir de ce qu'on échantillonne un niveau à plus fort contraste. Ça rend
ton hypothèse crédible, ça ne la démontre pas.

Pour le VCI, on utilise bien le même, van Ewijk et la fonction de lidR. La seule différence est
le zmax passé en argument, qui fixe le nombre de tranches. Je prends le maximum de hauteur de
la placette elle-même, donc un peuplement de 3 m a 3 tranches et peut très bien être réparti
dessus ; avec un plafond fixe à l'échelle du massif il aurait 40 tranches dont 3 occupées et
s'effondrerait vers 0. C'est ça qui fait tes points près de zéro, pas le rayon : j'ai recalculé
à ton rayon de 5 m et je retrouve exactement la même distribution qu'à 10 ou 25 m.

Je préfère garder ma normalisation, pour une raison nette : avec le plafond fixe le VCI corrèle
à 0.89 avec la hauteur maximale, il devient une métrique de hauteur déguisée, alors qu'avec la
mienne il tombe à 0.50 et mesure ce qu'on lui demande. Quitte à ne pas être superposable à ta
figure 4. Et sur le rayon je reste dans ta borne, les simulations tournent sur un disque de
10 m, avec le balayage 5 à 50 m déjà dans le chapitre.

Je garde la grille doublée en lignée parallèle pour l'instant. Dis-moi si tu penses qu'elle
doit entrer dans le chapitre, ça changerait la validation et une partie de la discussion.

Bien à toi,

Nathan
