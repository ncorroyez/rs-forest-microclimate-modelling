Version courte du 4e échange. Pièce jointe : Fig_logslope_airres1.png

---

Bonjour Jérôme,

Ton idée m'a poussé plus loin que prévu.

Je ne pouvais pas descendre sous 1.08 m dans les grands peuplements parce que la grille d'air
de MuSICA n'a que 10 couches. Je l'ai passée à 20, c'est juste un paramètre de namelist, et le
nœud le plus bas tombe à 0.26 m. Pas 0.1 m, il faudrait 28 couches et le modèle n'en propose
que 10 ou 20, mais assez pour tester en canopée fermée, ce qui était impossible jusqu'ici.

Sur le biais tu as partiellement raison : lire à 0.2 m au lieu de 1 m en enlève 0.16 à
0.20 °C en canopée fermée. Réel, mais ça ne fait que 12 % d'un biais de 1.6 °C.

Le vrai gain est ailleurs, et je ne l'attendais pas. La corrélation sur le ΔTmax passe de 0.50
à 0.72, celle sur la pente de 0.78 à 0.88, et l'amplitude restituée de 15 à 22 %. Sur la
figure jointe, le modèle lu près du sol suit bien mieux le nuage des capteurs. Le plafond de
verre tient, 22 % reste loin de 100 %, mais il se fissure.

Une réserve quand même : lire le modèle à 0.2 m face à des capteurs posés à 1 m ne se justifie
que si ton mécanisme est le bon. Le gain pourrait aussi venir de ce qu'on échantillonne un
niveau à plus fort contraste. Ça rend ton hypothèse crédible, ça ne la démontre pas.

Je garde ça en lignée parallèle pour l'instant. Dis-moi si tu penses que ça doit entrer dans
le chapitre, ça changerait la validation et une partie de la discussion.

Bien à toi,

Nathan
