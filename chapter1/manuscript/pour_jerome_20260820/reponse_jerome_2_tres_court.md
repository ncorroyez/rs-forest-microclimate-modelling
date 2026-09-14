Version très courte (les chiffres sont dans note_avant_apres.md).

---

Bonjour Jérôme,

Bon courage pour le rapport DR1.

Tu as vu juste sur les deux points.

Les x ont changé de footprint : le poster utilisait un disque de rayon 25 m, la nouvelle
figure la maille de 20 m qui contient le logger. Le VCI a aussi changé de définition, il est
maintenant calculé sur hauteurs normalisées au sol, et le troisième axe est le LAI une face
et non le PAI, d'où le facteur 2.

Les log(pente) > 0.2 en milieu ouvert venaient de la hauteur d'extraction. Le poster prenait
le modèle à `nair == 1`, le premier niveau d'air, qui vaut 0.0285 × Hmax, soit 7 cm dans une
placette ouverte. J'ai refait tourner les simulations de mars en ne changeant que ça :
+0.291 à `nair == 1`, +0.112 à 1 m fixe, et les quatre placettes concernées sont toutes en
dessous de 7 m de haut. Le chapitre échantillonne maintenant tout le monde à 1 m, la hauteur
des HOBO, donc l'amplification en canopée ouverte n'est plus gonflée par un capteur virtuel
posé au sol. Ce n'est pas la version de MuSICA.

Bien à toi,

Nathan
