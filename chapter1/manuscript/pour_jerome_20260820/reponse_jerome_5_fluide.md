Version fluide. Pièce jointe : Fig_logslope_iter_noiter.png

---

Bonjour Jérôme,

Pas de souci pour les aller-retours, ils m'ont fait trouver des choses que je n'aurais pas
cherchées.

La figure est en pièce jointe, 1 m avec itération et 1 m sans, tout le reste identique.
Prends-la pour ton rapport si elle te sert.

Tu as raison sur l'amplitude, et c'est moi qui me suis mal exprimé. Mon 22 % était la pente
d'une régression du ΔTmax simulé sur le ΔTmax observé, donc une tout autre variable que la
pente micro-macro. Si je fais le même calcul sur log(pente), je trouve 44 % pour la baseline et
53 % avec la lecture basse, et 39 % sur la pente non transformée. Tes 50 et 40 % tombent donc
exactement juste, je parlais simplement d'autre chose.

Ton intuition sur l'itération marche, mais d'un seul côté. Sans elle le modèle décomprime
nettement, l'amplitude sur log(pente) passe de 44 à 54 %. Seulement tout le gain est du côté
tamponnant : en canopée dense le log(pente) moyen va chercher -0.209 au lieu de -0.171, contre
-0.259 observé, alors qu'en canopée ouverte il ne bouge pas d'un cheveu. Sur la figure les deux
courbes se séparent à droite et restent collées à gauche.

Sur ce bout ouvert, justement, j'ai testé tout ce que j'ai pu. Le maximum de log(pente) y vaut
+0.015 en baseline, +0.020 sans itération, +0.023 sans correction de vent, +0.016 en rendant
aux placettes leur vrai couvert au lieu du plancher, et +0.034 en lisant à 4 cm. Le terrain
monte à +0.312. Aucun levier ne rattrape le dixième de l'écart, et je crois maintenant que
c'est structurel : sur ces placettes le modèle reçoit un LAI de 0.00 à 0.36, donc une colonne
quasi vide forcée cinq mètres au-dessus d'un sol nu. L'air à 1 m est asservi au forçage, la
pente se colle à 1, et c'est bien ce qu'il rend, +0.013 quelle que soit la placette. Le
terrain, lui, mesure +2 à +4 °C d'amplification, produite par un sol nu qui chauffe une couche
de quelques décimètres. Ce n'est pas une sous-amplification à corriger, c'est un mécanisme que
le 1-D n'a pas.

Et la différence avec SilviLaser ne vient pas de la hauteur de lecture. Sur les simulations de
mars, même macro et même période, `nair == 1` donnait +0.291 et 1 m fixe +0.112 : la hauteur y
comptait beaucoup. Aujourd'hui elle ne compte presque plus. Ce qui a disparu entre les deux,
c'est la sensibilité du modèle près du sol en canopée ouverte, pas ma façon de le lire.

Un mot enfin sur la Fig B2, parce que la lecture s'y retourne. Son point le plus bas est à la
hauteur relative 0.0285, soit environ 1 m en absolu pour une P4 de 35 m : elle ne dit donc rien
sous 1 m, et le réchauffement que tu y vois est entre 3 m et 1 m. Sous 1 m, avec la grille
fine, le profil s'inverse et refroidit de 0.20 °C en canopée dense, ce qui est cohérent avec un
sol ombragé à midi. C'est d'ailleurs pour ça que descendre le capteur virtuel gagne un peu en
canopée fermée et rien du tout en canopée ouverte.

Bien à toi,

Nathan
