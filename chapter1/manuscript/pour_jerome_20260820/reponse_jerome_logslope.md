Réponse à Jérôme, demande du 20 août 2026 : la figure log(slope) vs Hmax / VCI / PAI,
observée et modélisée, pour comparer les simulations d'avant et d'après l'été.
Pièces jointes : Fig_logslope_vs_structure_native20.png, Fig_logslope_vs_structure_evolution.png,
Fig_logslope_AVANT_ETE_2026-04.png.

---

Bonjour Jérôme,

Merci, ça avance bien de mon côté.

Oui, je l'avais en stock, mais seulement dans l'ancienne version. Je l'ai donc régénérée
sur la lignée actuelle, et je t'envoie les deux.

La première pièce jointe est l'analogue direct du panneau du poster : capteurs de terrain
en vert, MuSICA en orange, log(pente) contre Hmax, VCI et LAI. La seconde superpose en
violet les simulations d'avant l'été, sur la même abscisse, pour que la comparaison se
fasse d'un coup d'oeil. Je joins aussi la version d'avril pour mémoire.

Ce qui a bougé, en trois chiffres. L'étendue de log(pente) simulée passe de 0.16 à 0.35,
contre 0.68 observés : le modèle couvrait le quart du gradient observé, il en couvre
maintenant la moitié. Le partage tamponnement / amplification tombe juste, 8 loggers à
pente supérieure à 1 dans le modèle comme sur le terrain, alors qu'avant l'été le modèle
en amplifiait 31 sur 53. Et la courbe orange descend franchement avec la hauteur et avec
le LAI, là où la violette restait collée à zéro.

Le plafond de verre, lui, est intact, et il se lit mieux qu'avant sur cette figure. Le
modèle plafonne à une pente de 1.02 quand le terrain monte à 1.37. Sur les placettes
ouvertes, à faible Hmax ou faible LAI, les points verts restent nettement au-dessus des
oranges : le modèle ne sait toujours pas amplifier.

Un point que je préfère te donner tout de suite plutôt que tu le trouves : la corrélation
avec la pente observée baisse, de 0.93 à 0.78. L'ancienne lignée classait mieux les
placettes, mais sur un intervalle quatre fois trop étroit, donc un classement presque
parfait d'une amplitude qui n'existait pas. La nouvelle restitue l'amplitude et perd un
peu de rang. Pour l'argument du chapitre c'est l'amplitude qui compte, puisque c'est elle
qui porte le ΔTmax, mais la baisse est réelle et je la dis.

Deux réserves de lecture, sur les axes. Le VCI est recalculé sur hauteurs normalisées au
sol, donc il court ici de 0.05 à 0.93 alors que le poster affichait 0 à 0.73 ; l'effet est
que les placettes denses se tassent entre 0.80 et 0.93 et que l'ajustement quadratique y
travaille sur un nuage très concentré. C'est pour ça que le R² du modèle chute à 0.39 dans
ce panneau seulement : c'est un effet de la métrique, pas une perte de performance. Et le
troisième panneau porte le LAI une face du chapitre, de 0.02 à 6.71, pas le PAI brut du
poster qui montait au-delà de 12 ; seule la forme est comparable, pas l'abscisse.

Tout est régénérable d'un coup avec scripts/c1_logslope_vs_structure.R, et la table par
logger est à côté si tu veux refaire tes propres ajustements.

Bien à toi,

Nathan
