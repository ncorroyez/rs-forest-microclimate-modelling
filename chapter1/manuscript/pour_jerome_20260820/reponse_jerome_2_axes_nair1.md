Réponse à Jérôme, 2e échange du 20 août 2026 : pourquoi les données ont changé entre le
panneau (c) du poster SilviLaser et la figure log(slope) régénérée, et où sont passées les
log(pente) > 0.2 en milieu ouvert.

---

Bonjour Jérôme,

Bon courage pour le rapport DR1, et merci d'avoir regardé d'aussi près : les deux points que
tu soulèves sont justes, et le second m'a fait trouver autre chose.

D'abord une vérification, pour que tu saches que la reconstruction est fidèle. En remettant
l'ancien macro de référence, je retrouve l'observé du poster presque exactement : log(pente)
de -0.376 à +0.365, contre environ -0.39 à +0.35 lus sur la figure. Ce qui suit n'est donc
pas une reconstruction approximative, ce sont bien les mêmes données.

**1. Oui, le footprint a changé, mais ce n'est pas le seul changement.**

Le poster porte les métriques d'un disque de rayon 25 m. La figure d'aujourd'hui porte la
valeur de la seule maille native de 20 m qui contient le logger, soit 400 m² contre environ
1960 m², cinq fois plus petit. Ça se voit surtout sur la hauteur : Hmax allait de 6.5 à
43.5 m (médiane 29.5), il va maintenant de 2.9 à 37.6 m (médiane 26.7). Un disque de 25 m
attrape presque toujours un gros arbre, une maille de 20 m peut tomber entièrement dans une
trouée.

Le VCI, lui, a changé de définition et pas seulement de rayon. Il est maintenant calculé sur
des hauteurs normalisées au sol ; l'ancien raster travaillait sur l'élévation absolue et
comprimait d'environ 0.68×. De 0.19 à 0.66 (médiane 0.58) auparavant, contre 0.05 à 0.93
(médiane 0.85) aujourd'hui.

Le troisième axe est un simple changement d'unité. Le PAI du poster allait de 0.02 à 12.65
(médiane 8.04), le LAI une face d'aujourd'hui de 0.02 à 6.71 (médiane 3.93) : un rapport de
médianes de 2.04, cohérent avec le passage de deux faces à une face. Le rayon n'y joue
presque rien.

Et l'axe des y observé a bougé un peu lui aussi, parce que le macro de référence est passé
de `musica_in_Blois.nc` à `FR-Blo_2021_v2.nc` (r = 0.92 entre les deux séries horaires,
écart moyen de 0.96 °C). L'étendue de log(pente) observée passe de 0.740 à 0.675 et le
maximum de +0.365 à +0.312.

**2. Les log(pente) > 0.2 en milieu ouvert : c'est la hauteur d'extraction, pas la lignée.**

Le poster échantillonnait le modèle à `nair == 1`, le premier niveau d'air. Ce n'est pas 1 m :
sa hauteur relative vaut 0.02849 dans tous les fichiers, soit 0.07 m dans une placette
ouverte de 2.5 m et 1.07 m dans un peuplement de 37.5 m. En canopée ouverte le point de
mesure se retrouvait quasiment au sol, là où le modèle amplifie fort.

Je l'ai vérifié sur les mêmes NetCDF de mars, avec exactement la convention de la série
violette de ma figure (macro `musica_in_Blois.nc`, JJAS 2021), en ne changeant que la hauteur
d'échantillonnage :

- à `nair == 1` : log(pente) de -0.067 à **+0.291**
- à 1 m fixe interpolé : log(pente) de -0.048 à +0.112

Quatre placettes dépassent +0.2 sous `nair == 1`, et toutes ont un Hmax inférieur à 7.1 m.
Ce sont exactement tes milieux ouverts. Elles ne disparaissent donc pas parce que la version
de MuSICA a changé, mais parce que le chapitre interpole désormais toutes les placettes à 1 m
fixe, la hauteur des HOBO, pour que la comparaison se fasse à la même hauteur physique
partout. L'amplification en canopée ouverte était gonflée par un point de mesure posé à 7 cm
du sol.

**3. Un point que j'ai trouvé en vérifiant.**

Le panneau du poster n'était pas cohérent entre ce qui a été simulé et ce qui a été tracé.
Les simulations de mars ont été paramétrées avec un footprint proche de 5 à 10 m de rayon :
leur `veget_height_top` s'écarte au maximum de 2 à 3 m du balayage r = 5 et r = 10, contre
27 m pour r = 25. La hauteur est quantifiée au demi-mètre dans les sorties, donc r = 5 et
r = 10 ne se départagent pas, et 5 à 10 m est le plafond honnête. Or le poster les trace
contre les métriques du disque de 25 m.

La maille native 20 m d'aujourd'hui en est plus proche, mais pas identique : écart maximal de
9.5 m, écart moyen de 1.0 m. Je ne prétends donc pas que l'abscisse de ma figure soit celle
avec laquelle la série violette a tourné, seulement qu'elle en est nettement plus proche que
le r = 25. Si tu préfères la stricte comparabilité visuelle avec le poster, je peux retracer
les x sur le disque de 25 m, c'est une ligne à changer.

Bien à toi,

Nathan
