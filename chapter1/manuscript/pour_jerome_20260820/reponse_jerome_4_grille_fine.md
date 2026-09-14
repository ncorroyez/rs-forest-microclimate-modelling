Réponse à Jérôme, 20 août 2026 (4e échange) : grille verticale doublée, lecture près du sol.
Pièce jointe : Fig_logslope_airres1.png. Détails : tab_airres1_compare.csv.

---

Bonjour Jérôme,

Ton idée m'a poussé plus loin que prévu, et ça valait le coup.

Je ne pouvais pas descendre sous 1.08 m dans les grands peuplements parce que la grille d'air
de MuSICA n'a que 10 couches de végétation. Je l'ai passée à 20, ce qui n'est qu'un paramètre
de namelist, et le nœud le plus bas tombe à 0.26 m au lieu de 1.08 m sur la placette la plus
haute. Un vrai 0.1 m resterait hors d'atteinte, il faudrait 28 couches et le modèle n'en
propose que 10 ou 20, mais 0.26 m suffit pour tester ton hypothèse là où elle compte, en
canopée fermée, ce qui était impossible jusqu'ici.

J'ai séparé les deux effets, parce que la grille fine redécoupe aussi la canopée et n'est donc
pas une simple relecture plus basse de la même simulation. À 1 m dans les deux cas, la grille
seule n'apporte rien : le biais passe de +0.98 à +1.03 °C et l'amplitude se comprime un peu
plus. C'est en descendant que ça bouge.

Sur le biais, tu as partiellement raison. Dans les tiers moyen et haut, lire à 0.2 m au lieu
de 1 m enlève 0.16 à 0.20 °C. C'est réel, mais ça ne représente que 12 % d'un biais de
1.6 °C : la hauteur d'échantillonnage n'explique pas le réchauffement du modèle en canopée
fermée.

En revanche le gain sur la structure spatiale est net, et je ne l'attendais pas. Contre la
baseline du chapitre, la corrélation sur le ΔTmax passe de 0.50 à 0.72, celle sur la pente de
0.78 à 0.88, et l'amplitude restituée de 15 à 22 %. La figure jointe le montre bien : lu près
du sol, le modèle suit le nuage des capteurs beaucoup plus fidèlement, surtout sur la hauteur
et le LAI. Le plafond de verre tient toujours, 22 % reste très loin de 100 %, mais il se
fissure.

Une réserve que je me fais à moi-même. Lire le modèle à 0.2 m pour le comparer à des capteurs
physiquement posés à 1 m n'est légitime que si ton mécanisme est le bon, l'abri du tronc.
L'amélioration pourrait tout aussi bien venir de ce qu'on échantillonne un niveau à plus fort
contraste spatial. Ça rend ton hypothèse nettement plus crédible, ça ne la démontre pas. Il
faudrait un logger en air libre à côté d'un logger sur tronc pour trancher.

Pour l'instant je laisse ça comme une lignée parallèle, la baseline du chapitre reste à 1 m
sur 10 couches. Mais si tu penses que ça mérite d'entrer dans le chapitre, dis-le moi : ça
changerait la validation et une partie de la discussion.

Bien à toi,

Nathan
