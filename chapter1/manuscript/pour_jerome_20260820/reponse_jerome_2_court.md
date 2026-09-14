Réponse à Jérôme, 20 août 2026 (2e échange). Détail complet dans note_avant_apres.md.
Pièces jointes : Fig_logslope_vs_structure_native20.png, Fig_logslope_vs_structure_evolution.png.

---

Bonjour Jérôme,

Bon courage pour le rapport DR1.

Tu as raison sur les deux points.

**Le footprint.** Le poster trace les métriques d'un disque de rayon 25 m, la nouvelle figure
la maille native de 20 m qui contient le logger. D'où un Hmax qui descend beaucoup plus bas
(2.9 m au minimum au lieu de 6.5 m) : un disque de 25 m attrape presque toujours un gros
arbre, une maille de 20 m peut tomber dans une trouée. Deux autres choses ont bougé sur ces
axes : le VCI est maintenant calculé sur hauteurs normalisées au sol, ce qui l'étale de
0.19-0.66 à 0.05-0.93, et le troisième axe est le LAI une face et non le PAI, d'où le facteur
2 sur l'abscisse. Enfin l'axe des y observé a bougé un peu aussi, le macro de référence étant
passé de `musica_in_Blois.nc` à `FR-Blo_2021_v2.nc` : l'étendue de log(pente) observée passe
de 0.740 à 0.675.

**Les log(pente) > 0.2 en milieu ouvert.** C'est la hauteur d'extraction. Le poster
échantillonnait le modèle à `nair == 1`, le premier niveau d'air, qui n'est pas 1 m : sa
hauteur relative vaut 0.0285, donc 7 cm sous une canopée de 2.5 m et 1.07 m sous 37.5 m. Sur
les mêmes simulations de mars, en ne changeant que ça, log(pente) monte à +0.291 à
`nair == 1` et plafonne à +0.112 à 1 m fixe, et les quatre placettes au-dessus de +0.2 sont
toutes à Hmax inférieur à 7.1 m. Le chapitre interpole maintenant tout le monde à 1 m, la
hauteur des HOBO, donc l'amplification en canopée ouverte n'est plus gonflée par un point de
mesure posé au sol. Ce n'est pas la version de MuSICA qui les a fait disparaître.

Un dernier point, trouvé en vérifiant tout ça. Les simulations du poster avaient tourné sur
la valeur du pixel 10 m sous le logger, pas sur le disque de 25 m contre lequel la figure les
trace : c'est la première cause du déplacement horizontal des points. Ma figure a le même
défaut mais très réduit, 2.1 m sur Hmax, les simulations actuelles ayant tourné sur un MAX en
disque de 10 m alors que l'axe porte la maille de 20 m. J'ai refait la figure en mettant
chaque série sur les traits qui l'ont réellement paramétrée : les nuages ne bougent
pratiquement pas et aucun des chiffres ci-dessus ne change. Je te l'envoie si tu veux la
voir, mais ça ne change rien à la lecture.

Bien à toi,

Nathan
