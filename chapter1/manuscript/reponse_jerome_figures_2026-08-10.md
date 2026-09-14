Réponse à Jérôme, retours figures du 9 août 2026 (Chapitre 1, manuscrit AoFS).
Points traités : Fig 4, 5, 6, 8, B1, B3. Deux corrections apportées à ses hypothèses (Fig 6, normalisation B3).

---

Bonjour Jérôme,

Merci pour ces retours, très utiles. J'ai traité ce qui était actionnable tout de suite, et j'ai vérifié les
questions de méthode dans le code plutôt que de te répondre de mémoire.

**Fig 4 et 5.** Parfait, je ne touche à rien. Ta lecture de la Fig 4 est celle que je fais aussi : le contraste
entre P1 et P2-P4 tient à ce que le LAI de P1 est concentré là où on estime Tmax.

**Fig 6.** C'est fait, les deux panneaux partagent maintenant l'ordre du panneau (a), le ΔTmax observé, donc un
logger occupe la même ligne à gauche et à droite. Mais ton intuition ne se confirme qu'à moitié : les deux ordres
observés sont en fait très concordants, Spearman rho = 0.93. Le panneau (b) reste donc quasi monotone dans ce
nouvel ordre, il n'y a pas de réarrangement spectaculaire à lire.

Sur les sites que le modèle prédirait bien sur les deux métriques, une mise en garde de méthode d'abord. Les
erreurs des deux panneaux sont elles-mêmes corrélées (rho = 0.76), donc ce n'est pas un sous-ensemble
particulier. Et surtout |erreur ΔTmax| corrèle à +0.85 avec |ΔTmax observé| : comme le modèle comprime
l'amplitude à 15 % avec un biais de +0.98 °C, il « réussit » mécaniquement là où le tamponnement observé est
faible. Les dix meilleurs sur les deux critères incluent deux P1 à fCover 0.23 et 0.01 dont le ΔTmax observé est
quasi nul. J'ai donc contrôlé par |ΔTmax observé|, plus LAI et fCover, avant de tester tes deux pistes.

Sur l'erreur en ΔTmax (panneau a), aucune des deux ne ressort : variabilité spatiale du LAI, partiel r = +0.16
(p = 0.27) ; VCI, partiel r = +0.12 (p = 0.41).

Sur l'erreur en pente (panneau b), c'est plus intéressant. La variabilité spatiale du LAI semble forte en brut
(r = +0.52) mais s'effondre au contrôle (partiel r = -0.13, p = 0.36) : c'était l'ouverture déguisée, le CV
spatial du LAI corrèle à -0.77 avec le fCover. Ta seconde piste, en revanche, tient dans la direction que tu
devinais : à quantité de feuilles et couvert fixés, un profil verticalement moins dispersé est mieux prédit
(partiel VCI / erreur = +0.39, p = 0.004). Attention, le signe s'inverse entre le brut (-0.48) et le contrôlé,
donc il ne se lit qu'avec le contrôle explicite.

Je ne le mettrais pas dans le chapitre en l'état. L'intervalle bootstrap sur ce partiel est [-0.12, +0.69], il
traverse zéro, et le VCI est très concentré chez nous (médiane 0.85, trois quarts des placettes entre 0.81 et
0.93) avec une seule à 0.05 qui pèse lourd. Disons une piste cohérente avec ton intuition, à retester ailleurs,
plutôt qu'un résultat.

**Fig 8.** Oui, exactement, une simulation par rayon : 7 rayons (5, 10, 12.5, 15, 20, 25, 50 m) x 60 runs
loggers = 420 simulations MuSICA, avec les traits ré-extraits à chaque rayon. 53 loggers entrent dans chaque
corrélation, le réseau en compte 60 mais 53 ont une observation exploitable.

**Fig B1.** J'ai fait la figure, et elle apporte quelque chose, mais pas ce qu'on espérait. En brut,
r = -0.48 (p < 0.001), les peuplements top-heavy tamponnent davantage. Une fois LAI et fCover contrôlés, le
partiel tombe à r = 0.00 (p = 1.00) : la pente apparente est de la quantité de feuilles déguisée,
r(LAI, hauteur pondérée) = +0.45, et ça se voit à l'oeil sur le dégradé de couleur. J'ai tracé l'axe des x sur
0-1, la plage que B1 balaie en simulation, et l'IQR observé (0.43-0.52) n'en occupe que 9 %. La platitude est
donc le résultat, pas un échec : le levier existe, les peuplements réels n'en parcourent pas assez pour
l'exercer. Je l'ai mise en annexe (Fig. B4), en complément de l'added-variable plot déjà présent, qui lui
masquait la plage réellement occupée.

**Fig B3.** Trois réponses.

Oui, les pas scalaires sont tous positifs : le script ne lit que les runs LAI+, fCover+ et Hmax+, et la
différence est (perturbé - base). Une nuance, ΔLAD n'est pas un pas mais le contraste réel-moins-uniforme, donc
il n'est pas sur le même pied que les trois autres.

Pour Hmax, la correction que tu proposes est déjà en place, et c'est vérifiable directement : sur une placette à
Hmax 32 m et sa perturbée à 33 m, les dix couches de végétation ont des hauteurs relatives strictement
identiques. Les couches MuSICA sont proportionnelles à la hauteur de canopée, donc chaque profil est bien
normalisé par son propre Hmax, Hmax+ΔHmax compris.

D'où la réponse à ta dernière question : l'appariement est relatif, pas absolu. En absolu les niveaux diffèrent,
la couche 5 est à 10.15 m contre 10.47 m. Concrètement la courbe orange répond « à la même fraction de hauteur
de canopée », pas « à 1 m fixe pendant que la canopée grandit ». Ça se referme sur la Fig 4, qui mesure ΔTmax à
1 m absolu fixe et y trouve Hmax tout aussi négligeable (±0.04 °C), donc le classement « Hmax dernier partout »
ne dépend pas de ce choix. J'ai ajouté tout ça dans la légende de B3, ta question montrait qu'un relecteur la
poserait.

Un défaut résiduel, mineur : au-dessus du sommet de canopée les niveaux sont à un décalage métrique fixe (+2.02,
+6.07 m), donc leur hauteur relative diffère légèrement entre runs, alors que le tracé étiquette tout avec la
grille du run de base. Sur la plage affichée l'écart est inférieur à 0.01 unité relative, soit environ 30 cm.
Sous la canopée, il est nul.

Bonne lecture du chapitre, et merci encore.

Nathan
