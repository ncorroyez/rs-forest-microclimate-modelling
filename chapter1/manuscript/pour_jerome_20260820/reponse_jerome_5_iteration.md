Réponse à Jérôme, 21 août 2026 : figure itération / sans itération, définition de l'amplitude,
et pourquoi le modèle n'amplifie pas en canopée ouverte.
Pièce jointe : Fig_logslope_iter_noiter.png

---

Bonjour Jérôme,

Pas de souci pour les aller-retours, ils m'ont fait trouver des choses que je n'aurais pas
cherchées.

La figure que tu demandes est en pièce jointe, 1 m avec itération et 1 m sans, tout le reste
identique. Prends-la pour ton rapport si elle te sert.

**Tu as raison sur l'amplitude, et je me suis mal exprimé.** Mon 22 % est la pente d'une
régression du ΔTmax simulé sur le ΔTmax observé, donc une autre variable que la pente
micro-macro. Sur log(pente) je trouve bien ce que tu lis :

| lecture | amplitude sur ΔTmax | sur log(pente) | sur la pente brute |
|---|---|---|---|
| baseline, 1 m | 15 % | 44 % | 39 % |
| grille fine, nœud bas | 22 % | 53 % | 48 % |

Tes 50 % et 40 % tombent exactement sur la baseline. Je préciserai la métrique la prochaine
fois plutôt que de dire « l'amplitude ».

**L'itération : ton intuition marche, mais d'un seul côté.** Sans elle le modèle décomprime
nettement, l'amplitude sur log(pente) passe de 44 à 54 % et l'étendue de 51 à 65 % de
l'observé. Mais tout le gain est du côté tamponnant : en P4 le log(pente) moyen va chercher
-0.209 au lieu de -0.171, contre -0.259 observé, alors qu'en P1 il ne bouge pas, +0.016 contre
+0.013. Sur la figure les deux courbes se séparent à droite et restent collées à gauche.

**Sur le bout ouvert, j'ai testé tout ce que je pouvais.** Maximum de log(pente) en P1 :
+0.015 en baseline, +0.020 sans itération, +0.023 sans correction de vent, +0.016 en rendant
aux placettes leur vrai couvert au lieu du plancher à 0.5, +0.034 en lisant à 4 cm. Le terrain
monte à +0.312. Aucun levier ne dépasse le dixième de l'écart.

La raison est structurelle et je la crois maintenant solide. Sur ces placettes le modèle reçoit
un LAI une face de 0.00 à 0.36 : la colonne est quasi vide, forcée à Hmax + 2 m, soit environ
5 m au-dessus d'un sol nu. L'air à 1 m est alors asservi au forçage, donc pente ≈ 1 et
log(pente) ≈ 0, et c'est exactement ce qu'il rend, +0.013 quelle que soit la placette. Le
terrain, lui, mesure +2 à +4 °C d'amplification du ΔTmax, produite par un sol nu qui chauffe
une couche de quelques décimètres. Ce n'est pas une sous-amplification à corriger, c'est un
mécanisme que le 1-D n'a pas.

**Sur SilviLaser, la différence n'est pas la hauteur de lecture.** Sur les simulations de mars,
même macro et même période, `nair == 1` donnait +0.291 et 1 m fixe +0.112 : la hauteur y
comptait beaucoup. Aujourd'hui elle ne compte presque plus, +0.013 à 1 m contre +0.018 au nœud
bas. Ce qui a disparu entre les deux, c'est la sensibilité du modèle près du sol en canopée
ouverte, pas la façon dont je le lis.

**Sur la Fig B2, la lecture se retourne.** Son point le plus bas est à la hauteur relative
0.0285, ce qui fait environ 1 m en absolu pour une P4 de 35 m. B2 ne dit donc rien sous 1 m :
le réchauffement que tu y vois est entre 3 m et 1 m. Sous 1 m, avec la grille fine, le profil
s'inverse en canopée dense et refroidit de 0.20 °C, ce qui est cohérent avec un sol ombragé à
midi. C'est d'ailleurs pour ça que descendre le capteur virtuel réduit un peu le biais en
canopée fermée et rien du tout en canopée ouverte.

Bien à toi,

Nathan
