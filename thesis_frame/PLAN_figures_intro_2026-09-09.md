# Figures de l'introduction générale — inventaire et propositions (2026-09-09)

Réponse à : « mettre des figures dans l'intro, explorer les dossiers pour voir les
figures existantes, proposer d'autres à faire (en codage, pas gen par IA) ».

Rien n'a été modifié dans `intro_generale_these_EN.md`. Ce fichier est la liste à
valider / raturer avant que je code quoi que ce soit.

---

## 0. État actuel

Cinq figures, **aucun tableau**. Toutes conceptuelles, toutes produites par un script
`fig_intro_*.R` sourçant `_intro_fig_style.R` :

| N° | Section | Ce qu'elle dit | Script |
|----|---------|----------------|--------|
| 1.1 | §1.2.2 | le LAI est une intégrale, et une intégrale jette de l'information | `fig_intro_lai_integral.R` |
| 1.2 | §1.3.2 | multicouche vs big-leaf : seul le multicouche accepte un profil | `fig_intro_musica.R` |
| 1.3 | §1.5.2 | quelle part de la colonne chaque capteur atteint | `fig_intro_sensor_column.R` |
| 1.4 | §1.5.3 | le compromis structure 3D / revisite, et la case vide | `fig_intro_tradeoff.R` |
| 1.5 | §1.6 | l'enchaînement des trois gaps → trois chapitres | `fig_intro_roadmap.R` |

Deux contraintes héritées :

- **Règle de style.** `_intro_fig_style.R` déclare en toutes lettres que les figures de
  l'intro ne portent **aucune valeur mesurée et aucun label d'axe numérique**. C'est une
  décision, pas un oubli.
- **Numérotation en dur.** Les légendes contiennent `**Figure 1.n.**` écrit à la main.
  **Vérifié :** aucun renvoi `Figure 1.x` dans le corps de l'intro, ni dans
  `discussion_generale_these_EN.md`, ni dans `manuscript_chap1_EN_native20.md`, ni dans
  `Chapter3_article_standalone_EN.md`. Insérer une figure au milieu ne casse donc que les
  numéros des légendes suivantes — renumérotation mécanique, en une passe, sans risque.

**Sections sans figure** : 1.1.1, 1.1.2, 1.1.3, 1.1.4, 1.2.1, 1.3.1, 1.3.3, 1.4.1,
1.4.2, 1.5.1, 1.7.

---

## 1. Ce qui reste ouvert de la passe précédente

Trois propositions de la session du 26 août n'ont jamais été construites. Je les remets
en tête parce qu'elles priment sur tout ce que j'ajoute aujourd'hui.

### 1.1 — §1.5.1 · Le verrou d'attribution *(c'était ma recommandation principale)*
Deux panneaux : à gauche, les peuplements réels alignés en bande étroite dans l'espace
quantité × arrangement, où une régression n'ajuste qu'une direction et non un effet ; à
droite, la croix de perturbations qui déplace un trait à la fois autour d'un point de
cette bande. C'est le pivot méthodologique qui justifie l'existence du Ch1, et il est
toujours invisible. **Reste ma priorité n°1.**

### 1.2 — §1.2.2 · Tableau du vocabulaire de la surface foliaire
LAI vrai / LAI effectif / PAI / fCover / indice d'agrégation × ce que l'instrument
renvoie × convention une face ou deux faces × ce que MuSICA attend. Ton texte affirme
que le désaccord entre produits LAI « is definitional before it is physical » ; sans
tableau, la phrase reste une pétition de principe. C'est aussi le seul endroit où figer
proprement le piège du ×2.

### 1.3 — §1.7 · Tableau des trois chapitres
Sites, capteurs, échelle, modèle, validation, sortie. Argument spécifique : tu as écarté
une section Matériel & Méthodes commune, donc les faits concrets (trois sites, 53
capteurs, grille 10 m, été 2021) ne sont regroupés nulle part.

---

## 2. Inventaire : ce qui existe ailleurs dans tes dossiers

2 668 fichiers figure dans les deux dépôts. Filtre appliqué : **je ne retiens que ce dont
je peux nommer le script producteur** — sinon ce n'est pas « du codage », c'est un PNG
orphelin dont on ne peut ni vérifier ni régénérer le contenu.

| Figure existante | Producteur | Verdict |
|---|---|---|
| `outputs/figures_chap1_pointclouds/side_P{1..4}.png` — nuages vus de côté, échelle z commune 0–38 m | `c1_pointcloud_sideview.R` | **Recyclable comme brique** dans la figure 3.1. ✔ disque `MyPassport` monté, vérifié : régénérable depuis le dépôt. |
| `pc3d_typical_P{1..4}.png` — vues 3D obliques | `c1_pointcloud_typical_fig1.R` | Alternative. Plus joli, moins lisible pour lire un profil vertical. |
| `out_files/Chapter1/figures/FigStation_archetype_map.png` — 53 loggers, Blois | `c1_archetype_map.R` | **Re-coder, pas réutiliser** : en français, et colorée par archétype P1–P4 — ce qui divulgue la typologie du Ch1 dès l'intro. |
| `paper1_Figures/Fig_1.png` — France + 3 sites + rasters LAI_ALS/LAI_S2 | *aucun script trouvé* | Producteur introuvable, et les rasters côte à côte **sont déjà le résultat du Ch2**. À re-coder dépouillé si 3.3 est retenue. |

**Écarté volontairement** — figures de résultats de chapitre : `Fig01_archetypes`,
`Fig1_saturation_styled`, `Fig4_dopt_mechanism`, `Fig0_lai_sensitivity`, tout
`manuscripts/ch4/figures/`, tout `figs_MEB2026_final/`. Les mettre dans l'intro générale
grille la chute des chapitres et crée un doublon dans le manuscrit assemblé. Je préfère
que ce soit écrit plutôt que tu croies que je ne les ai pas vues.

---

## 3. Nouvelles propositions

### Priorité haute

#### 3.1 — §1.4.1 · La chaîne de restitution LiDAR *(conceptuelle)*
`fig_intro_als_chain.R`

Quatre panneaux : nuage normalisé → fraction de trouées couche par couche → inversion
Beer-Lambert → profil de densité foliaire, LAI intégré comme aire sous la courbe. Les
hypothèses annotées **là où elles entrent** : le coefficient d'extinction *k* sur la
flèche d'inversion, la distribution aléatoire du feuillage sur la couche.

C'est le point que §1.4.1 martèle — « a LiDAR LAI is a model-based estimate resting on
stated assumptions, not a direct measurement » — sans rien pour le soutenir.

*Vérifié en ouvrant la Fig 1.3* : elle trace bien un profil de densité foliaire, mais
pour dire **quelle part de la colonne chaque capteur atteint**. Elle ne montre nulle part
la fraction de trouées ni l'inversion. Recouvrement visuel (même idiome de profil),
pas de recouvrement d'argument. Panneau (a) : un vrai `side_P{3,4}.png` plutôt qu'un
dessin — seul écart à la règle de style, et un écart honnête (une image, pas une valeur).

#### 3.2 — §1.4.2 · L'inversion PROSAIL et son caractère mal posé *(conceptuelle)*
`fig_intro_prosail_chain.R`

Pendant exact de 3.1, pour l'autre capteur : réflectance observée → LUT de spectres
simulés → régresseur entraîné → LAI. Un encart montre **deux jeux de paramètres distincts
produisant des spectres quasi superposés** — l'ill-posedness rendue visible, et la raison
pour laquelle les *priors* sont le troisième facteur du Ch2.

*Encart saturation — à arbitrer, mais la ligne de partage est nette.* La prose de §1.4.2
s'engage déjà sur le **phénomène** (« exhausted somewhere above five », avec citations) :
un schéma sans graduations ne fait que redire le texte. Ce qui appartient au Ch2, c'est le
**mécanisme** — d_opt, la surface foliaire enfouie sous la profondeur optique, les 6 à 10 m
site par site. À noter quand même : la Fig 1.3 engage **déjà** « effective optical depth »
dans une figure de l'intro. Le seuil est donc plus haut que je ne l'avais supposé, et
l'encart passe sans difficulté.

#### 3.3 — §1.1.4 · Le système d'étude et la limite spatiale des loggers *(données)*
`fig_intro_study_system.R`

**Seule figure porteuse de données que je propose**, et je propose de l'assumer : une
deuxième classe, avec sa propre convention de légende (source + année), plutôt que de
casser en douce la règle du fichier de style.

Trois panneaux : (a) carte de France, trois sites ; (b) Blois, MNH en fond, les 53
loggers **d'une seule couleur** (pas d'archétype : la typologie appartient au Ch1) ;
(c) la surface couverte par les 53 placettes rapportée à celle du massif.

**Deux réserves que je dois poser franchement.**

*Je me contredis.* En août j'ai écarté un visuel « capteur ponctuel vs champ continu »
en §1.1.4 comme « une évidence déjà dite en deux phrases ». Ce qui a changé : je propose
maintenant une carte du **système d'étude**, pas une illustration de l'argument. La
question de l'évidence reste valable pour le panneau (c) seul.

*Problème d'emplacement, pas de contenu.* §1.1.4 est une section de **littérature** —
SoilTemp, blindage des HOBO, plan d'échantillonnage, tout cité. Y déposer une carte la
transforme en section de méthodes au milieu d'un paragraphe, et l'intro n'a aucune section
« sites d'étude » pour l'accueillir. À décider **avant** de coder : nouvelle §1.1.5 ?
§1.7 ? Ça change la légende, le jeu de panneaux, et si (c) est le bon troisième panneau.
La plus chère du lot (`sf`, MNH, geojson).

#### 3.4 — §1.1.2 · Tampon et amplification : le changement de signe *(conceptuelle)*
`fig_intro_buffering_sign.R`

ΔTmax en ordonnée, ouverture de la canopée en abscisse, la courbe traversant le zéro :
fermé → tampon, ouvert → amplification. Encart jour/nuit : le même feuillage refroidit le
jour et réchauffe la nuit, et seule la branche diurne atténue les extrêmes.

ΔTmax est **la métrique de toute la thèse**, elle est définie dans ce paragraphe, et elle
n'a aucune figure. Le changement de signe prépare aussi le résultat densité-dépendant du
Ch1 — le plan FR le note explicitement.

### Priorité moyenne

#### 3.5 — §1.3.1 · Corrélatif vs mécaniste *(conceptuelle)* — **re-proposition**
`fig_intro_stat_vs_mech.R`

En août je l'avais proposée en **tableau** et classée dernière, en disant qu'elle
risquait de paraphraser la prose. Je la reprends en **figure** parce qu'un tableau
paraphrase, alors qu'un schéma montre ce que la prose ne peut que décrire : à gauche,
des traits colinéaires entrent et les contributions individuelles ressortent
indissociables ; à droite, le bilan d'énergie est calculé, donc **un trait peut être bougé
seul**.

Réserve, après avoir ouvert la Fig 1.2 : celle-ci porte déjà en pied de panneau
« *the arrangement can be changed at fixed leaf area* ». L'idée de perturbation est donc
déjà effleurée. Ce qui reste propre à 3.5 est le côté gauche — l'impasse corrélative.
**Et cette moitié-là est exactement le panneau gauche de la figure 1.1 (§1.5.1).**
Concrètement : si tu prends 1.1, ne prends pas 3.5. Elles se recouvrent à moitié.

#### 3.6 — §1.3.3 · De la placette au paysage *(conceptuelle)*
`fig_intro_plot_to_landscape.R`

Une colonne MuSICA → la même colonne répliquée sur une grille de structure télédétectée.
Dit pourquoi §1.4 arrive juste après. Vérifié : ni la Fig 1.2 ni la Fig 1.4 ne montrent
cette réplication. Reste néanmoins la proposition la plus dispensable de la liste.

#### 3.7 — Tableau des capteurs — **je déconseille**
La Fig 1.4 dit déjà exactement ça, et mieux : elle place en plus la case vide « ce dont un
modèle mécaniste sur paysage aurait besoin ». Un tableau ne ferait que la paraphraser.

---

## 4. Courbe de phénologie

Tu avais dit « courbe de phéno à voir pk pas, j'ai les données ». Elle n'a jamais été
construite. Une réserve, puis trois options.

**Ne pas tracer `blois_s2_lai_ts_2021.rds`.** Deux raisons indépendantes, établies dans ce
dépôt : (1) `NC_Full/DIAGNOSTIC_s2_summer_decline.md` attribue la chute de −39 % entre juin
et août à la **géométrie d'illumination post-solstice**, pas à la phénologie — cette série
contient un artefact identifié ; (2) **aucun script du dépôt ne produit ces `.rds`**.
Le document le plus lu de la thèse est le pire endroit où mettre une série non traçable
contenant un artefact connu.

Trois options acceptables, par ordre de préférence :
- **(a) Courbe schématique**, sans graduations, intégrée à 3.2 : trajectoire saisonnière
  continue de S2 contre le trait unique de la campagne ALS. Zéro risque.
- **(b) Phénologie prescrite MuSICA** (`phenology_musica_veg1.csv`, `Leaf_area_1yr` par
  jour julien) en §1.3.2 : illustre ce que le modèle attend en entrée. Traçable, mais
  c'est une entrée de modèle, pas une observation — à dire dans la légende.
- **(c) Fenêtre débourrement → maximum seulement** de la série S2, série et limite nommées
  dans la légende. Défendable techniquement ; je ne le ferais pas dans l'intro générale.

---

## 5. Ordre proposé

1. Tu rayes ce qui ne t'intéresse pas. **Le seul arbitrage qui bloque : 1.1 ou 3.5** (elles
   se recouvrent à moitié), et **l'emplacement de 3.3**.
2. Je code **1.1** (verrou d'attribution) — priorité n°1, purement conceptuelle.
3. Puis **3.1 et 3.2** (les deux chaînes de restitution) : les deux plus gros trous du
   texte, conceptuelles, `_intro_fig_style.R` réutilisé tel quel.
4. Puis **3.4** (changement de signe).
5. Les deux tableaux (1.2, 1.3), bon marché.
6. **3.3** en dernier : seule à demander des données et une décision de style de ta part.
7. Renumérotation des légendes en une passe, une fois la liste gelée.

Si tout est retenu, on arrive à **onze figures et deux tableaux** pour ~10 000 mots.
C'est au-dessus de ce que je te disais tenable en août (cinq figures, deux tableaux pour
35 pages). À toi de trancher : pour serrer, l'ordre de sacrifice est 3.6, puis 3.7, puis
le panneau (c) de 3.3.

---

## 6. Ajout du 2026-09-09 (2ᵉ passe) — concepts généraux, MaCCMic, MuSICA, technologies

Passe demandée sur les figures « générales / globales ». J'ai fouillé `Presentations/`,
`MaCCMic_Imprint/`, `Admin/ED_GAIA/` et `musica/documentation/`, que la première passe
n'avait pas ouverts. Cinq trouvailles, dont deux qui changent des décisions plus haut.

### 6.1 — La figure MuSICA existe, et elle est de qualité impression
`Admin/ED_GAIA/CSI2/MuSICA_general_diagram.jpg` — **2321 × 1615 px**, deux panneaux :
(A) la discrétisation (couches de végétation, espèces / classes d'âge, types de feuilles
éclairées-ombragées-mouillées-sèches, couches sol/litière, profils racinaires) ;
(B) les processus (transfert radiatif, interception pluie/neige, photosynthèse et
respiration, transpiration/évaporation, transfert turbulent et microclimat, transferts
d'eau et de chaleur).

C'est très probablement la figure que tu avais en tête en août (« j'en ai un normalement
dans mes dossiers qqpart »). Je ne l'avais pas trouvée et j'avais codé la Fig 1.2 à la
place.

**Recommandation : garder les deux, elles ne disent pas la même chose.** La Fig 1.2
argumente *pourquoi* le schéma de canopée doit être multicouche (multicouche vs big-leaf,
« le profil ne peut pas être accepté » à droite). Celle-ci dit *ce que MuSICA contient*.
La première est un argument, la seconde une description. Si tu n'en veux qu'une en §1.3.2,
garde la Fig 1.2 : c'est elle qui porte l'argument du Ch1.

Variante économique : ne reprendre que le **panneau B**, qui est celui où « microclimat »
apparaît explicitement comme sortie du transfert turbulent.

⚠️ **Attribution.** Figure de Jérôme Ogée (auteur de MuSICA — la doc `MuSICA_Short_Description.pdf`
le confirme comme auteur). Il faut son accord et un crédit en légende. C'est une ligne de
mail, mais il faut la faire : c'est la seule figure du lot qui n'est pas de toi.

### 6.2 — La figure MaCCMic : je recommande de la re-coder, pas de l'importer
`Presentations/figure_MaCCMic_web_EN-768x342.jpg` (EN) et `maccmic_fr.png` /
`Admin/ED_GAIA/CSI3/figures/fig00_maccmic_schema.png` (FR, 1019 × 426).

Elle dit exactement ce que dit §1.1.2 : deux mâts, T_micro sous couvert et T_macro en
zone ouverte, la couverture forestière atténuant les extrêmes, et la biodiversité du
sous-bois répondant au microclimat. C'est littéralement la définition de ΔTmax de ta thèse,
en image.

**Quatre raisons de ne pas l'importer telle quelle :**
1. **Résolution.** La version EN est en 768 px de large — environ 122 dpi sur 16 cm. La FR
   monte à 162 dpi. Aucune n'atteint 300 dpi, et je n'ai pas trouvé de source vectorielle
   (fouillé les PDF MaCCMic et le pptx de nov. 25).
2. **Langue.** La version haute résolution est la française ; le manuscrit est en anglais.
3. **Contenu.** Conifères et joggeurs : c'est une figure de communication de projet, pas
   une figure de thèse sur les chênaies-hêtraies décidues.
4. **Attribution.** Figure du projet ANR, crédit obligatoire.

**Recommandation : fusionner avec la proposition 3.4.** Une seule figure en §1.1.2, deux
panneaux : (a) les deux mâts, T_micro / T_macro, ΔTmax défini — le cadrage MaCCMic ;
(b) le changement de signe le long de l'ouverture de canopée — tampon puis amplification.
Ça règle résolution, langue, contenu et attribution d'un coup, et ça donne à §1.1.2 la
figure d'ouverture qui lui manque. **3.4 devient de fait ma priorité n°2, après 1.1.**

### 6.3 — `Presentations/france.geojson` : l'objection de coût sur 3.3 tombe
Le fond de carte France existe déjà en geojson dans ton dossier de présentations. Le panneau
(a) de la figure « système d'étude » (3.3) devient trivial. Il reste le MNH de Blois pour le
panneau (b), mais ce n'était pas la partie chère. **La question qui bloque sur 3.3 reste
l'emplacement, pas le coût.**

### 6.4 — La comparaison LiDAR / S2 « sur les technologies » : un tableau, pas une figure
Tu la demandes explicitement. Analyse honnête de ce qui est déjà couvert :

- La **Fig 1.3** trace déjà les impulsions laser aller-retour jusqu'au sol et le trajet
  optique du satellite. La géométrie d'acquisition y est donc **à moitié faite**.
- La **Fig 1.4** donne déjà revisite × information structurelle × emprise.

Ce qui n'est nulle part, c'est la **physique de la mesure** : actif contre passif, ce que
l'instrument mesure réellement (temps de vol contre luminance réfléchie), ce que renvoie une
acquisition (retours discrets géoréférencés contre une valeur par bande et par pixel), la
source d'énergie, et ce qui limite chacun (plan de vol et densité d'impulsions contre
couverture nuageuse et angle solaire).

C'est un **tableau**, pas une figure : ce sont des attributs appariés, pas une géométrie.

**Je reviens donc sur mon « je déconseille » du 3.7**, et je dis ce qui a changé : je l'avais
écarté en le lisant sur l'axe de la Fig 1.4 (compromis spatio-temporel), où il faisait
effectivement doublon. Sur l'axe instrumental que tu proposes, il ne recoupe rien. Le
tableau reste à écrire, la Fig 1.4 reste à sa place.

Une figure de géométrie d'acquisition (avion + fauchée laser à gauche, satellite + pixel
10 m à droite) reste possible en ouverture de §1.4, mais elle recouperait la Fig 1.3 à
peu près de moitié. Je ne la recommande pas tant que le tableau n'existe pas.

Matière disponible pour ce tableau, et pour un second sur les bandes S2 :
`Presentations/s2bands.png` est une **capture d'écran d'une page web** — inutilisable telle
quelle, mais c'en est le contenu : B02–B12, longueur d'onde, résolution native 10/20/60 m,
et la colonne qui compte pour §1.4.2 — quelles bandes l'inversion utilise réellement. Ton
texte insiste sur ce compromis 10 m contre red-edge et SWIR ; un petit tableau le fige.

### 6.5 — Deux images secondaires, pour information
- `Presentations/lidardata.png` — grand nuage de points ALS coloré par hauteur, layons et
  hétérogénéité de canopée bien visibles. **Palette arc-en-ciel**, non sûre pour les
  daltoniens : à régénérer avec une rampe correcte via la machinerie `c1_pointcloud_*.R`
  si tu la veux. Bonne image d'ouverture pour §1.4.1 si tu préfères du réel au schéma.
- `Presentations/Hemispherical_photo1.jpg` — photographie hémisphérique. Illustrerait bien
  la ligne « LAI effectif / mesure indirecte » du tableau de vocabulaire (1.2). Décoratif
  seul, utile en vignette dans le tableau.

### 6.6 — Ordre révisé
1. **1.1** — verrou d'attribution (§1.5.1). Inchangé, priorité n°1.
2. **3.4 fusionnée avec 6.2** — ΔTmax et le changement de signe (§1.1.2). Monte en priorité :
   c'est la figure d'ouverture du manuscrit et elle absorbe le cadrage MaCCMic.
3. **3.1 et 3.2** — les deux chaînes de restitution (§1.4.1, §1.4.2).
4. **6.1** — décider si le diagramme MuSICA d'Ogée rejoint la Fig 1.2 en §1.3.2 (et si oui,
   lui écrire).
5. **Les tableaux** : vocabulaire du LAI (1.2), les trois chapitres (1.3), technologies
   LiDAR/S2 (6.4), bandes S2 (6.4). Tous bon marché.
6. **3.3** — système d'étude, une fois l'emplacement tranché.

Les arbitrages qui bloquent sont maintenant au nombre de trois : **1.1 ou 3.5** ;
**l'emplacement de 3.3** ; **une ou deux figures MuSICA en §1.3.2**.

---

## 7. Ajout du 2026-09-09 (3ᵉ passe) — figures publiées trouvées dans ta bibliothèque

Décisions enregistrées :
- **MaCCMic : gardée telle quelle**, plus de re-codage. Choix concret à faire : la version EN
  fait 768 px de large (≈ 122 dpi sur 16 cm), la FR 1019 px (≈ 162 dpi). **Je recommande la
  EN, placée en largeur réduite** (12 cm plutôt que pleine page) pour remonter le dpi effectif.
  Conséquence : ma proposition **3.4 (le changement de signe de ΔTmax) redevient une figure
  codée à part entière**, elle n'absorbe plus MaCCMic.
- **Comparaison LiDAR / S2 : on cherche une figure publiée**, pas un tableau. Recherche en
  cours, état ci-dessous.

### 7.1 Quatre figures publiées, confirmées en les ouvrant

| Source | Figure | Section visée | Ce qu'elle apporte |
|---|---|---|---|
| **Zellweger et al. 2019**, *TREE* 34(4) | **Fig. 1** — réseau de capteurs dans un peuplement scanné au LiDAR (coupe du nuage de points, S1 en ouvert / S2 sous couvert), courbes Tmin–Tmax, cartes de hauteur de canopée et de topographie, puis **« calibrer et valider des modèles statistiques et prédire le microclimat sur le paysage »** → cartes de Tmax et de VPD | §1.3.1 ou §1.5.1 | **La trouvaille la plus utile de la session.** C'est l'image canonique de l'approche **statistique** que ta thèse conteste. §1.3.1 et §1.5.1 argumentent aujourd'hui contre une méthode que le lecteur n'a jamais vue. |
| **Zellweger et al. 2019**, encadré, **Fig. I(C)** | macroclimat vs microclimat en température maximale, l'**offset** et sa **tendance** en aires colorées | §1.1.2 | Complément exact de MaCCMic : MaCCMic donne le schéma de paysage, celle-ci donne la série temporelle. Les deux se répondent. |
| **Bramer et al. 2018**, *Adv. Ecol. Res.* 58 | **Fig. 1** — gradient vertical macroclimat → microclimat, et **une flèche par facteur forçant, dont la longueur est l'étendue verticale sur laquelle il agit** (propriétés du sol, évapotranspiration, ombrage de canopée, topographie, drainage d'air froid, rayonnement, précipitation, type de végétation, vent, proximité de l'eau, mélange atmosphérique) | §1.2.1 | §1.2.1 est de la prose de bilan d'énergie pure, sans aucune figure. Celle-ci la rend lisible d'un coup. Bémol : niveaux de gris, photo-composite, style un peu daté. |
| **Bramer et al. 2018** | **Fig. 5** — résolution spatiale × résolution temporelle des capteurs satellitaires, avec une bande **« Microclimate / Topo-Local / Meso »** sur l'axe spatial | §1.4.2 ou §1.3.3 | Montre que le microclimat exige < 10 m, et que peu de capteurs y sont. **Attention : aucun LiDAR dedans — ce n'est pas la comparaison de technologies demandée.** Recoupe ta Fig 1.4 sur l'axe revisite, en diffère sur l'autre (résolution spatiale vs information structurelle). |

**Conséquence sur les propositions existantes.** Si la Fig. 1 de Zellweger entre, ma
proposition **3.5 (corrélatif vs mécaniste) change de forme ou disparaît** : la figure
publiée *est* le panneau gauche, en mieux et en vrai. Ce qui resterait à coder est le
panneau droit seul — le modèle mécaniste où un trait peut être bougé seul.

### 7.2 La comparaison LiDAR / S2 sur les technologies — état de la recherche

Ouvert et écarté : Laslier et al. 2023 (microclimat Sentinel-1 + LiDAR + optique — mais
figures = site, workflow, biplots) ; Gavilan-Acuna et al. 2025 (ALS vs satellite optique en
sylviculture de précision — la Fig. 1 oppose bien les deux capteurs, mais sur des questions
de fertilisation, d'éclaircie et de récolte, hors sujet ici) ; Li et al. 2023 (actif vs
passif, mais étude de cas sur mangrove) ; Zhang & Lin 2017 (fusion en photogrammétrie).

**Diagnostic.** La littérature de fusion et d'application suppose que le lecteur sait déjà
ce qu'est un laser et ce qu'est un spectromètre : elle dessine des chaînes de traitement,
pas de la physique de la mesure. Une figure qui oppose le temps de vol actif à la
réflectance passive vit dans les **manuels et les supports d'enseignement**, pas dans les
articles de RSE. La recherche se poursuit de ce côté.

### 7.3 Permissions — une seule liste, à faire une fois
Elsevier (Bramer 2018 ; Zellweger 2019 dans *TREE*) et Springer (*Annals of Forest Science*)
autorisent tous la reproduction en thèse via RightsLink, en général gratuitement, mais il
faut la demander et créditer. Plus l'accord de Jérôme Ogée pour le diagramme MuSICA (6.1),
et le crédit ANR MaCCMic. À grouper en une session de dix minutes le jour où la liste est gelée.

### 7.4 Deux sections qui ne veulent pas de figure
**§1.1.1** (changement climatique et extrêmes) demanderait des données climatiques externes
aux deux dépôts — hors périmètre, et la section est courte. **§1.1.3** (pourquoi la
variabilité fine compte) est un argument de conséquences écologiques : une figure ne ferait
que le décorer. Je préfère te le dire que te fabriquer deux candidats faibles.

---

## 8. Ajout du 2026-09-09 (4ᵉ passe) — thèse de Bouwen, et verdict sur la figure « technologies »

### 8.1 La thèse de Klara Bouwen : quatre figures directement réutilisables
`~/Zotero/storage/WAPSAPHQ/` — son introduction est le voisin le plus proche du tien, et
elle a déjà fait le travail de sélection.

| Figure | Source d'origine | Section visée | Verdict |
|---|---|---|---|
| **Fig. I.1** — profils verticaux de température dans le peuplement : (a) ciel clair 3 h vs 13 h, (b) ciel couvert, (c) **minuit, canopée claire vs dense**, (d) **midi, canopée claire vs dense** | **De Frenne et al. 2021**, *Global Change Biology* | **§1.1.2** | **La meilleure trouvaille de cette passe.** Les panneaux (c) et (d) sont exactement ma proposition 3.4 — asymétrie jour/nuit **et** contraste clair/dense — mais publiés, avec de vrais profils, et De Frenne 2021 est déjà cité dans ton intro. Ça **remplace 3.4** ou la réduit à un complément. Bonus : ce sont des profils *verticaux*, ce qui prépare l'argument multicouche de §1.3.2. |
| **Fig. I.2** — interactions climat / peuplement | Aussenac 2000 | §1.2.1 | Classique de la littérature forestière française. Alternative à Bramer Fig. 1. |
| **Fig. I.3** — (A) échecs de plantation et mortalité des semis **en France**, (B) semis | — | **§1.1.3** | **Corrige ce que j'ai écrit en 7.4.** Je disais que §1.1.3 ne voulait pas de figure ; celle-ci donne à l'argument régénération une base chiffrée et française, au lieu de le décorer. À vérifier : quelle est sa source primaire. |
| **Fig. I.5** — (A) flux d'énergie de surface, (B) hydrologie, (C) cycle du carbone | **Bonan 2019** | §1.2.1 | Le **panneau A seul** est la figure de bilan d'énergie de §1.2.1. Les panneaux B et C débordent du sujet. Concurrente directe de Bramer Fig. 1 : Bramer est plus microclimat (chaque facteur avec son étendue verticale), Bonan est plus bilan d'énergie (les flux nommés). **Bramer si tu veux les facteurs, Bonan A si tu veux les flux.** |
| Fig. I.6 — schéma de résolution du modèle multicouche | — | §1.3.2 | Alternative au diagramme d'Ogée (6.1). Trois figures MuSICA candidates pour une seule section : il faut en éliminer. |

### 8.2 La comparaison LiDAR / S2 « sur les technologies » — verdict
Cherché sérieusement : la bibliothèque Zotero (~1 000 PDF, filtrée sur les 80 entrées
combinant LiDAR et optique), la thèse de Bouwen, et le web.

**Ouvert et écarté :**
- *Fassnacht et al. 2024*, **Forestry** 97(1), en libre accès CC BY — la grande revue du
  domaine. Vérifié : **aucune figure de comparaison de capteurs**, tout est en prose.
- *NASA ARSET, « The Fundamentals of LiDAR »* (Podest), diapositive 4 « Active and Passive
  Sensors » — domaine public, donc pas de permission à demander, mais **trop élémentaire** :
  un satellite, le soleil, la Terre. Rien sur la canopée. Indigne d'une intro de thèse.
- *Gavilan-Acuna et al. 2025*, **Annals of Forest Science**, Fig. 1 — oppose bien ALS et
  satellite optique par des flèches de couleur, mais sur des questions de fertilisation,
  d'éclaircie et de récolte. Structure correcte, domaine faux.
- *Bramer et al. 2018* Fig. 5 — résolution spatiale × temporelle, mais **satellites
  uniquement, aucun LiDAR**.
- *Laslier et al. 2023*, *Li et al. 2023*, *Zhang & Lin 2017*, *Landry et al. 2020* — chaînes
  de traitement et études de cas, pas de physique de la mesure.

**Pourquoi ça ne se trouve pas.** La littérature de fusion et d'application suppose que le
lecteur sait ce qu'est un laser et ce qu'est un spectromètre : elle dessine des workflows.
La figure qui oppose le temps de vol actif à la réflectance passive vit dans les manuels
(Lillesand & Kiefer, Jensen, Bonan) et les supports d'enseignement — donc soit trop
élémentaire, soit sous copyright de manuel, ce qui est le pire cas pour une permission.

**Ce que je propose à la place, et c'est mon avis honnête : tu l'as déjà.** Ta **Fig. 1.3**
trace les impulsions laser aller-retour jusqu'au sol contre le trajet optique du satellite,
sur un vrai profil de densité foliaire, avec la profondeur optique effective. Aucune figure
publiée que j'ai ouverte ne fait ça aussi bien pour *ton* argument. Ce qui lui manque, ce
n'est pas une figure concurrente, c'est **la ligne de texte qui nomme la différence
physique** : actif contre passif, temps de vol contre luminance réfléchie, retours discrets
géoréférencés contre une valeur par bande et par pixel.

Trois issues possibles, à toi de trancher :
- **(a)** Ajouter ces attributs en encart dans la Fig. 1.3 — coût quasi nul, et la figure
  devient la comparaison de technologies.
- **(b)** Coder une figure de géométrie d'acquisition en ouverture de §1.4 (avion et fauchée
  laser à gauche, satellite et pixel 10 m à droite). Recoupe la Fig. 1.3 à moitié, mais
  ouvre proprement la section.
- **(c)** Reprendre la Fig. 1 de Gavilan-Acuna 2025 en la reciblant sur *tes* questions
  (structure verticale, dynamique saisonnière, emprise) plutôt que sur la sylviculture de
  précision — c'est une adaptation, donc du codage, avec citation de l'original.

Je recommande **(a)**, puis **(b)** si tu veux vraiment une figure séparée.

### 8.3 Où on en est, section par section
| Section | Figure | Statut |
|---|---|---|
| §1.1.1 | — | Pas de figure (données climatiques hors dépôts) |
| §1.1.2 | **MaCCMic** (gardée telle quelle) + **De Frenne 2021 / Bouwen I.1** + Zellweger encadré I(C) | Trois candidates, en garder deux |
| §1.1.3 | **Bouwen I.3** (mortalité des semis en France) | Nouveau — corrige mon 7.4 |
| §1.1.4 | 3.3 système d'étude *(emplacement à trancher)* | À coder |
| §1.2.1 | **Bramer Fig. 1** *ou* **Bonan / Bouwen I.5A** | Choisir |
| §1.2.2 | Fig 1.1 intégrale (existante) + tableau vocabulaire LAI | Tableau à faire |
| §1.3.1 | **Zellweger Fig. 1** (l'approche statistique, en image) | Le meilleur ajout du lot |
| §1.3.2 | Fig 1.2 (existante) + **diagramme Ogée** *ou* **Bouwen I.6** | Trois candidates, en garder une de plus au maximum |
| §1.3.3 | 3.6 placette → paysage | Dispensable |
| §1.4.1 | 3.1 chaîne LiDAR | À coder |
| §1.4.2 | 3.2 chaîne PROSAIL + tableau bandes S2 | À coder |
| §1.5.1 | **1.1 verrou d'attribution** | Priorité n°1, à coder |
| §1.5.2 | Fig 1.3 (existante) — **+ encart technologies, option (a)** | Petit ajout |
| §1.5.3 | Fig 1.4 (existante) | Fait |
| §1.6 | Fig 1.5 feuille de route (existante) | Fait |
| §1.7 | Tableau des trois chapitres | À faire |

### 8.4 Permissions — liste consolidée
Elsevier : Bramer 2018, Zellweger 2019. Wiley : De Frenne 2021 (via RightsLink, gratuit en
thèse). Springer : Gavilan-Acuna 2025 si option (c). Bonan 2019 : éditeur du manuel.
Aussenac 2000 : *Annals of Forest Science*. Plus l'accord de Jérôme Ogée (diagramme MuSICA)
et le crédit ANR MaCCMic. **Reprendre une figure via la thèse de Bouwen ne dispense pas de
la permission de la source d'origine** — c'est l'éditeur d'origine qu'il faut solliciter,
pas Klara.

---

## 9. Fait — construit et inséré le 2026-09-09

Quatre figures codées, une figure importée, insérées et renumérotées. L'intro compile
en **39 pages, 10 figures, 0 citation non résolue**. Sauvegarde du markdown avant
insertion : `review/.bak_intro_2026-09-09.md` (fichier caché — `ls -a` pour le voir).

| N° | Section | Figure | Script / source |
|----|---------|--------|-----------------|
| **1.1** | §1.1.2 | **MaCCMic**, gardée telle quelle (version EN, 88 % de largeur) | `figures/Fig_intro_maccmic.jpg` |
| 1.2 | §1.2.2 | LAI est une intégrale | `fig_intro_lai_integral.R` |
| 1.3 | §1.3.2 | multicouche vs big-leaf | `fig_intro_musica.R` |
| **1.4** | §1.4 (ouverture) | **les deux capteurs comme instruments** | `fig_intro_sensor_instruments.R` *(nouveau)* |
| **1.5** | §1.4.1 | **chaîne de restitution LiDAR** | `fig_intro_als_chain.R` *(nouveau)* |
| **1.6** | §1.4.2 | **chaîne PROSAIL et ses deux limites** | `fig_intro_prosail_chain.R` *(nouveau)* |
| **1.7** | §1.5.1 | **le verrou d'attribution** | `fig_intro_attribution_lock.R` *(nouveau)* |
| 1.8 | §1.5.2 | ce que chaque capteur atteint | `fig_intro_sensor_column.R` |
| 1.9 | §1.5.3 | compromis spatio-temporel | `fig_intro_tradeoff.R` |
| 1.10 | §1.6 | feuille de route | `fig_intro_roadmap.R` |

Chaque nouveau script se termine par un `stopifnot`, comme les précédents, de sorte que
la figure vérifie sa propre affirmation :
- **1.4** — l'inversion de Beer-Lambert boucle exactement (2 × 10⁻¹⁶) et le profil est
  reconstruit à 0,45 % du pic ; la fraction de trouées reste dans [0, 1] et croît avec
  la hauteur.
- **1.5** — les deux canopées « différentes » ne diffèrent que de 4 % de l'étendue
  spectrale (c'est le caractère mal posé, rendu mesurable) ; la pente tombe à 8 % de sa
  valeur en canopée ouverte (c'est la saturation).
- **1.7** — le nuage a une dispersion transversale de 16 % de sa dispersion longitudinale
  (colinéarité), et **chacune** des quatre croix de perturbation atteint au moins un point
  hors de la bande occupée par les placettes.
- **1.4** — comparaison appariée : les quatre attributs sont renseignés pour les deux
  capteurs, aucune case vide.

### Décisions prises, plus à trancher
- **3.5 (corrélatif vs mécaniste) est absorbée dans 1.6.** Son panneau gauche *est* le
  panneau (a) du verrou d'attribution. Je ne la porte plus.
- **§1.1.2 : MaCCMic + De Frenne I.1**, et j'abandonne Zellweger encadré I(C) — c'est le
  même argument d'offset que les deux mâts de MaCCMic, avec moins d'apport.
- **La comparaison « technologies » est faite**, mais **en figure autonome (1.4) plutôt qu'en
  encart de la Fig 1.8** : celle-ci est déjà dense (deux panneaux, quatre entrées de légende,
  et elle porte « effective optical depth », terme qui appartient au Ch2). En ouverture de
  §1.4, la comparaison instrumentale donne aussi à la section la figure d'entrée qui lui
  manquait, sans toucher à une figure qui fonctionne.

### Convention Beer-Lambert — vérifiée
La Fig 1.5 écrit `P_gap(z) = exp(−k·PAI(z))`, `PAI(z) = −ln(P_gap(z))/k`, « puis dériver
avec la hauteur ». Vérifié : c'est **mot pour mot** la formule du Ch3 §2
(`LAI = −ln(P_gap)/k`) et de `R/lidar_lai.R` (`PAI = -log(P) / k`). La dérivation avec la
hauteur est l'étape MacArthur-Horn que le Ch1 §2 décrit. Aucune divergence de convention
entre la figure et les chapitres. À noter tout de même : **le Ch1 utilise k = 0,5 et le
Ch3 k = 0,65** — la figure ne donne aucune valeur, donc pas de conflit, mais l'écart
mérite une phrase quelque part.

### Reste à faire, par ordre
1. Les figures sous permission : **De Frenne 2021 / Bouwen I.1** (§1.1.2),
   **Zellweger 2019 Fig. 1** (§1.3.1), **Bramer Fig. 1 ou Bonan I.5A** (§1.2.1),
   **Bouwen I.3** (§1.1.3), **diagramme MuSICA d'Ogée** (§1.3.2).
3. Les tableaux : vocabulaire du LAI (§1.2.2), bandes S2 (§1.4.2), les trois chapitres (§1.7).
4. La carte du système d'étude (§1.1.4 ou nouvelle section), une fois l'emplacement tranché.

### Convention pour les figures empruntées
La légende de la Fig 1.1 porte **`[permission pending: ANR MaCCMic]`** en gras. Toute
figure importée entre avec ce marqueur, de sorte qu'une compilation ne puisse jamais
produire silencieusement un manuscrit contenant une figure non autorisée. Le marqueur se
retire quand l'accord est obtenu.

---

## 10. Fait — deuxième vague (2026-09-09, suite)

L'intro compile en **42 pages, 11 figures, 3 tableaux, 0 citation non résolue**.

### 10.1 Nouveautés
| N° | Section | Contenu | Source |
|----|---------|---------|--------|
| **Tableau 1.1** | §1.2.2 | vocabulaire de la surface foliaire : LAI vrai, PAI effectif ALS, profil, fCover, LAI optique, surface deux faces | Ch1 §2 (l. 142-144, 164), `R/musica.R:65`, `R/validation.R:96` |
| **Tableau 1.2** | §1.4.2 | les dix bandes Sentinel-2, longueur d'onde S2A, résolution native | @druschSentinel2ESAsOptical2012 + Ch2 (discussion) |
| **Fig 1.11** | §1.7 | **le système d'étude** : France + trois forêts ; Blois + les 53 loggers | `fig_intro_study_system.R` *(nouveau)* |
| **Tableau 1.3** | §1.7 | ce que chaque chapitre utilise : question, forêt, capteurs, modèle, résolution, unité d'analyse, juge | Ch1 §2.1/2.3, Ch2 §2.1, Ch3 §2 |

### 10.2 L'emplacement de la carte : tranché
Elle n'est **pas** allée en §1.1.4. Cette section est de la littérature (SoilTemp, blindage
des HOBO, plan d'échantillonnage) et y déposer une carte la transformait en méthodes au
milieu d'un paragraphe. Elle ouvre **§1.7**, avec le Tableau 1.3 : le lecteur y rencontre
pour la première fois « Blois », « les 53 loggers » et « trois sites », et la section
devient « voici la thèse, voici où et avec quoi ». L'argument « un logger mesure un point »
reste textuel en §1.1.4, où il est à sa place.

### 10.3 Deux entorses assumées
- **La Fig 1.11 porte des positions mesurées**, contrairement aux dix autres. Sa légende le
  dit en toutes lettres (« Unlike the other figures of this chapter, this one carries
  measured positions »), et son script est le seul à dépendre de `sf`. C'est la deuxième
  classe annoncée en 3.3, pas une dérive du style.
- **Il y a désormais un renvoi croisé dans le texte** : « Figure 1.11 locates the forests »
  en tête de §1.7. C'était zéro jusqu'ici. Toute insertion **après** la feuille de route
  (§1.6) doit donc mettre ce renvoi à jour à la main.

### 10.4 Provenance des 53 loggers, vérifiée
`in_files/data_Blois_utm31n.geojson` contient **60** placettes. Les 53 retenues viennent de
`outputs/figures_pipeline/data/clusters.rds`. Vérifié : ce jeu de 53 identifiants est
**identique** dans les trois lignées de pipeline présentes (`figures_pipeline`,
`figures_pipeline_z05`, `figures_pipeline_z05_v323`), donc le filtre n'engage la figure sur
aucune d'elles. Les sept retirées : 41_13, 41_14, 41_20, 41_41, 41_50, 41_51, 41_53.
Le script vérifie `nrow(plots) == 53` et que les trois sites tombent dans le contour de la
France métropolitaine.

### 10.5 ⚠️ Une tension entre chapitres, à trancher côté science
En construisant le Tableau 1.2 j'ai buté sur une contradiction que je n'invente pas :

- **Ch2 (discussion, RSE r2)** : « We restricted the inversion to the 10-m spectral bands,
  as finer resolution is critical in heterogeneous forests […]. Additional investigations
  are necessary to identify the relevance of lower spatial resolution red-edge and SWIR
  spectral bands (20 m). » → l'inversion tourne sur **quatre** bandes.
- **Ch3 §2** : « The red-edge and short-wave-infrared bands, acquired natively at 20 m, were
  therefore resampled to 10 m. These bands carry much of the leaf-area information least
  affected by saturation. » → se lit comme si elles **entraient** dans la restitution.

Les deux ne peuvent pas être vraies de la même chaîne. Ch3 dit par ailleurs prendre les
résultats du Ch2 « as given », donc soit Ch3 a relancé sa propre inversion avec un jeu de
bandes différent, soit sa phrase décrit seulement la grille de sortie et est trompeuse.

**Ce que j'ai fait en attendant** : la légende du Tableau 1.2 attribue la restriction au
**Ch2 seul** (« Chapter 2 restricts the inversion to the four 10 m bands »), ce qui est
exact et n'engage pas le Ch3. À toi de trancher, puis de rectifier la phrase du Ch3 ou la
légende.

### 10.6 Le coefficient d'extinction
Le Ch1 utilise *k* = 0,5 (MacArthur-Horn, distribution sphérique) et le Ch3 *k* = 0,65.
La légende du Tableau 1.1 pose le principe sans donner les valeurs : le coefficient est un
choix déclaré, pas une propriété mesurée, et comme il rééchelonne la restitution par une
constante il laisse le classement entre placettes inchangé. **Les valeurs elles-mêmes
relèvent des chapitres** ; l'écart 0,5 / 0,65 mérite une phrase quelque part, mais pas dans
l'introduction générale.

### 10.7 Permissions
Liste consolidée, modèles de courriel et cases à cocher dans
**`review/PERMISSIONS_figures_intro.md`**. Quatre figures retenues attendent un accord et
**ne sont pas insérées** ; la seule importée (Fig 1.1) porte son marqueur.

---

## 11. Correctifs de dernière passe, et ce qu'ils ont révélé

### 11.1 Les longueurs d'onde du Tableau 1.2 — corrigées
La première version venait d'une capture d'écran (`Presentations/s2bands.png`) : c'étaient les
valeurs **S2A en vol** (492, 704, 741, 833, 1614, 2202). La légende cite pourtant Drusch et
al. 2012. J'ai ouvert le Tableau 3 de Drusch et repris **ses** valeurs de conception
(490, 705, 740, 842, 1610, 2190). Six lignes sur dix changeaient. Le tableau et sa source
disent maintenant la même chose, et plus aucun chiffre du manuscrit ne remonte à une capture
d'écran. Ajouté aussi la clause qui manquait : les trois bandes à 60 m servent la correction
atmosphérique et l'écran nuageux, d'où dix lignes pour treize bandes.

### 11.2 ⚠️ L'inversion PROSAIL n'utilise pas quatre bandes, mais **trois**
En vérifiant la tension signalée en 10.5, j'ai remonté le code plutôt que la prose :

- `02_CODES/Sentinel_2/Main_02_produceLAI.R:26` → `S2BandSelect <- list('lai' = c('B3','B4','B8'))`
- `02_CODES/Sentinel_2/3_train_predict_prosail.R:194-197` → `bands_10m <- c('B3','B4','B8')`,
  puis `bands_select <- list(bands_10m)` qui **écrase** la ligne précédente offrant aussi
  `bands_20m <- c('B3','B4','B5','B6','B7','B8A','B11','B12')`.

Donc l'inversion tourne sur **B3 (vert), B4 (rouge), B8 (PIR)**. La bande bleue B2 n'entre
pas, alors qu'elle est à 10 m. Le jeu à huit bandes existe dans le code, désactivé, et il a
servi à une étude comparative (`03_RESULTS/*/Plots/@atbd_bands_noise_study/`,
`bands_3_4_8` contre `bands_3_4_5_6_7_8A_11_12`).

**Conséquence pour le Ch2** : sa phrase « We restricted the inversion to the 10-m spectral
bands » est imprécise, c'est trois des quatre. À resserrer si l'occasion se présente.

### 11.3 ⚠️ Et la tension du 10.5 est tranchée : la phrase du Ch3 est fausse
Les scripts du Ch3 lisent
`output/intermediate/sm6/<site>/s2lai_summer_atbd_T_res_10_m.tif`
(`c3_hybrid_lai.R:16`, `c3_s2_improvement.R:14`) — c'est-à-dire **le produit de la chaîne du
Ch2**, sans suffixe de jeu de bandes, donc le défaut `bands_3_4_8`. Le Ch3 ne relance pas
d'inversion.

Donc, en §2 du Ch3, « The red-edge and short-wave-infrared bands, acquired natively at 20 m,
were therefore resampled to 10 m. These bands carry much of the leaf-area information least
affected by saturation » **décrit une chose qui n'a pas lieu** : ces bandes n'entrent pas du
tout dans la restitution que le Ch3 utilise. La phrase décrit sans doute la grille de travail
à 10 m et a glissé. Elle dit aujourd'hui au lecteur quelque chose de faux sur la méthode du
chapitre, et la discussion du Ch3 s'appuie dessus (« We did not test a retrieval at the
native 20 m resolution followed by super-resolution to 10 m »).

**C'est une correction à faire dans le Ch3, pas dans l'introduction.** Le seul point qui
touchait l'intro, la légende du Tableau 1.2, est déjà neutre : elle dit « Chapter 2 restricts
the inversion to the 10 m bands », ce qui reprend la formulation du chapitre sans avancer un
nombre. J'avais écrit « the four 10 m bands » : retiré.

Les dix bandes du Ch3 apparaissent bien quelque part, mais ailleurs : c'est la forêt
aléatoire sur réflectances brutes (`c3_archetype_scores.R:96`), un diagnostic qui contourne
la restitution. Aucun rapport avec l'inversion.

---

## 12. Passe finale — figures empruntées insérées, et le Ch3 corrigé

L'intro compile en **47 pages, 15 figures, 3 tableaux, 0 citation non résolue**, plus aucun
marqueur de permission. Les autorisations sont ramenées à une formalité d'avant-dépôt,
listée dans `PERMISSIONS_figures_intro.md` (réduit à une page).

### 12.1 Quatre figures publiées, extraites et insérées
| N° | Section | Figure | Extraction |
|----|---------|--------|-----------|
| **1.2** | §1.1.2 | profils verticaux de température, clair/couvert et clair/dense, minuit/midi — **De Frenne et al. 2021, Fig. 4** | `pdfimages` p24, image intégrée 912 × 912 |
| **1.3** | §1.1.3 | échecs de plantation et mortalité des semis en France, 2006-2023 — **panneau A de Bouwen 2025 Fig. I.3, d'après le DSF** | rendu 300 dpi p17, recadré sur le panneau A seul |
| **1.5** | §1.3.1 | la route statistique vers une carte de microclimat — **Zellweger et al. 2019, Fig. 1** | rendu 300 dpi p8 puis recadrage : `pdfimages` renvoyait l'image **sans le texte**, qui est vectoriel |
| **1.6** | §1.3.2 | ce que MuSICA contient, discrétisation et processus — **diagramme de Jérôme Ogée** | fichier d'origine 2321 × 1615 |

Piège rencontré : sur le PDF de Zellweger, `pdfimages` extrait un raster **muet**, tout le
texte étant vectoriel par-dessus. Il faut passer par `pdftocairo -r 300` puis recadrer.

### 12.2 Bramer écartée, après l'avoir regardée
`Bramer et al. 2018` Fig. 1 (les facteurs forçants et leur étendue verticale) est extraite et
lisible, mais c'est un photo-montage en niveaux de gris de style daté : elle jurerait avec les
onze figures dessinées, et la prose de §1.2.1 tient sans elle. Le fichier reste dans
`figures/_src/bramer-000.png`.

### 12.3 Ce que la nouvelle numérotation donne
1.1 MaCCMic · **1.2 profils De Frenne** · **1.3 régénération** · 1.4 intégrale ·
**1.5 route statistique** · **1.6 MuSICA d'Ogée** · 1.7 multicouche vs big-leaf ·
1.8 instruments · 1.9 chaîne LiDAR · 1.10 chaîne PROSAIL · 1.11 verrou d'attribution ·
1.12 colonne · 1.13 compromis · 1.14 feuille de route · 1.15 système d'étude.
Le renvoi croisé de §1.7 a suivi : « Figure 1.15 locates the forests ».

### 12.4 ⚠️ J'ai modifié un manuscrit de chapitre
`manuscripts/ch3/Chapter3_article_standalone_EN.md`, deux phrases. Sauvegarde :
`.bak_ch3_2026-09-09.md`.

**§2, avant** : « All analyses were run on the 10 m grid to match the microclimate
resolution. The red-edge and short-wave-infrared bands, acquired natively at 20 m, were
therefore resampled to 10 m. »
**après** : « All analyses were run on the 10 m grid to match the microclimate resolution,
and the inversion therefore used only the bands acquired natively at that resolution: green,
red and near infrared. The red-edge and short-wave-infrared bands, acquired at 20 m, do not
enter the retrieval. »

**§4, avant** : « …which resamples the red-edge and short-wave-infrared bands from their
native 20 m » → **après** : « …which leaves the red-edge and short-wave-infrared bands out of
the inversion ». Et la phrase suivante devient « a retrieval at the native 20 m resolution,
**using them**, followed by super-resolution », pour rester cohérente.

Base : `Main_02_produceLAI.R:26` → `c('B3','B4','B8')`, et
`3_train_predict_prosail.R:194-197` où `bands_select <- list(bands_10m)` écrase l'option à
huit bandes. Le Ch3 lit le produit du Ch2 (`c3_hybrid_lai.R:16`) et ne relance pas
d'inversion.

### 12.5 Trois retouches dans l'intro, d'après tes notes
- **§1.5.3** : une phrase était dupliquée mot pour mot (« They have to be addressed in
  order. »). Supprimée.
- **§1.4.2** : la saisonnalité S2 était vendue comme son avantage principal sans réserve.
  Ajouté que la réflectance est directionnelle, que la géométrie soleil-visée change au fil
  de la saison indépendamment de la canopée [@maignanBidirectionalReflectanceEarth2004 ;
  @breonAnalysisHotSpot2002], et qu'une trajectoire restituée peut donc dériver après le
  solstice pour des raisons étrangères à la chute des feuilles. C'est ton point « fiabilité
  de la saisonnalité, qu'est-ce que ça traduit réellement », et c'est exactement ce que
  `NC_Full/DIAGNOSTIC_s2_summer_decline.md` a diagnostiqué.
- **RQ3** : l'objectif listait les quatre scénarios sans jamais dire quel contraste organise
  le chapitre. Ajouté : surface foliaire satellite qui décrit la partie haute contre surface
  foliaire LiDAR qui intègre toute la colonne, donc le déficit attendu en dense est
  exactement la surface foliaire située sous la profondeur que le signal optique atteint, et
  les archétypes du Ch1 séparent les canopées où ce déficit compte de celles où il ne compte
  pas.

### 12.6 Ce que je n'ai pas touché, et pourquoi
Le reste de tes notes est **déjà dans l'intro** : quatre scénarios seulement, la pente
comme métrique primaire et l'offset en secondaire, la stratification par archétype, le
mécanisme d_opt en §1.7, le lien au Ch1. Je ne réécris pas ce qui concorde.
**FORMS-H et la modélisation de capteur restent hors de l'intro**, tu les places toi-même en
Discussion générale.
