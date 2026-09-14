# Chapitre 1 — Crib de défense DÉTAILLÉ (le « pourquoi » à fond)
2026-07-28. Pour chaque choix : **Quoi · Pourquoi (mécaniste) · Alternative rejetée · Si on me pousse.**

## A. PIPELINE

### 1. Traits ALS natifs à 20 m (pas 10 m agrégé)
- **Pourquoi 20 m :** MuSICA est une colonne 1-D qui représente un canopée homogène ; le HOBO intègre une source de quelques dizaines de mètres. 20 m = l'échelle où un peuplement est cohérent ET comparable à l'empreinte du capteur. Plus fin (10 m) = bruit sous-placette qu'un modèle 1-D ne sait pas exploiter.
- **Pourquoi « natif » (calculé direct à 20 m) et pas agrégé depuis 10 m :** agréger un trait **non-linéaire** (profil LAD, VCI) depuis des tuiles 10 m ≠ le calculer sur **tous les retours** de la cellule 20 m — moyenner des profils lisse la structure. Le natif = la vraie distribution des retours.
- **Si poussé :** recette natif-vs-agrégé validée (LAI r 0.97, Hmax 0.99, LAD 0.95) → le choix ne distord pas, il est plus propre ; et le fit est **insensible au rayon d'empreinte de 5 à 50 m** (Annexe A).

### 2. Les 4 traits : LAI, fCover, Hmax, profil LAD
- **Pourquoi ces 4 :** ils couvrent les axes qui *pourraient* piloter le buffering — **quantité** (LAI, interception), **fermeture** (fCover, ombrage/ciel vu), **taille** (Hmax, profondeur d'interception + rugosité), **arrangement** (profil LAD, où le feuillage intercepte). Les 3 premiers, tout capteur les voit ; le profil est la contribution **unique** du LiDAR — toute la question est de savoir s'il « mérite sa place ».
- **Alternative rejetée :** ajouter rumple, gap fraction, etc. → colinéaires avec ces 4, pas de levier mécaniste nouveau. Parcimonie.

### 3. VCI = entropie de Shannon van Ewijk, sur hauteurs sol-relatives, descripteur seulement
- **Pourquoi cette formule :** standard, 0–1, **comparable à Gril (même site)** ; mesure l'**évenness** de la distribution verticale des retours.
- **Pourquoi hauteurs sol-relatives et PAS z/Hmax :** VCI doit capter l'étalement vertical **en mètres absolus** (un canopée de 30 m étalé sur 25 m ≠ un de 10 m sur 8 m). Rescaler par Hmax effacerait l'info de hauteur que VCI encode en partie.
- **Pourquoi descripteur et JAMAIS input :** MuSICA ingère le **profil LAD complet** (toutes les couches), qui contient strictement plus d'info qu'un scalaire. VCI est un résumé 1-nombre pour la typologie/report, et surtout **colinéaire à la densité** (r 0.79 avec LAI) → le donner « comme la forme » serait faux (c'est exactement le piège que le harsh review a trouvé sur la Fig 6).

### 4. LAI ×2 avant MuSICA
- **Pourquoi :** le transfert radiatif de MuSICA compte la surface foliaire **totale (bilatérale)** pour l'interception ; l'inversion MacArthur-Horn du LiDAR rend du **one-sided (projeté)**. Ne pas doubler → interception divisée par 2 → bilan d'énergie faux.
- **Pourquoi ce N'EST PAS du clumping :** le clumping (Ω < 1, feuillage aggloméré) est une correction **séparée**. Le ×2 est purement la conversion d'unités one-sided → two-sided. Les confondre = double comptage.

### 5. k = 0.5 (coefficient d'extinction, inversion Beer-Lambert)
- **Pourquoi 0.5 :** distribution sphérique des angles foliaires → G = 0.5 ; défaut conventionnel feuillu quand l'angle est inconnu.
- **Pourquoi c'est un choix load-bearing :** k **fixe le LAI retrouvé** (LAI ∝ 1/k dans Beer-Lambert). k = 0.65 (optimisé pour la *précision de rétrieval* dans une étude compagnon) → LAI **−23 %** → **efface le skill** de buffering.
- **Parade :** on a pris la valeur **physiquement conventionnelle**, pas celle optimisée sur le résultat ; et on **reporte le résultat comme conditionnel à k = 0.5** (Annexe G). On ne peut pas être accusé d'avoir tuné l'input dont on teste l'importance.

### 6. FPCA sur profils LAD double-normalisés → FPC1-3
- **Pourquoi FPCA :** le profil LAD est une **fonction** (courbe), pas un scalaire ; la FPCA est la façon naturelle de décomposer des données fonctionnelles en modes orthogonaux (Ramsay). Réduit ~15 couches à quelques scores de forme interprétables.
- **Pourquoi double-normaliser (z/Hmax + intégrale unité) :** pour obtenir la **forme PURE**, scale-free en hauteur ET en quantité — sinon FPC1 ne ferait que récupérer Hmax et LAI (qu'on a déjà en scalaires). La double-normalisation retire ça → les composantes ne captent **que l'arrangement**.
- **Pourquoi FPC1-3 :** 92 % de la variance de forme ; au-delà = bruit.

### 7. k-means, 4 archétypes, sur 6 métriques
- **Pourquoi clusteriser :** avoir des archétypes **interprétables** pour la présentation et le report par archétype. **PAS pour piloter l'analyse** — chaque perturbation est par placette.
- **Pourquoi 6 métriques (3 scalaires + 3 FPC) :** pour que les archétypes diffèrent sur **tous** les axes (quantité, taille, fermeture, forme), pas seulement la densité.
- **Pourquoi 4 :** coude du WCSS près de 4, interprétable (ouvert→dense). **Assumer la silhouette 0.27 :** clustering faible → **partition de commodité d'un continuum**, pas des types naturels.
- **Parade « vos clusters sont arbitraires » :** oui, dispositif de présentation ; la science est **par placette et continue** — Fig 4 le prouve (chaque levier = fonction continue du point de fonctionnement, pas un saut aux frontières).

### 8. cLHS = 400 pixels RÉELS (100/cluster, 10 000 itérations, seed 42)
- **Pourquoi échantillonner :** simuler toute la population × 8 perturbations = trop coûteux ; il faut un sous-ensemble représentatif.
- **Pourquoi cLHS et pas aléatoire :** l'aléatoire sous-échantillonne les **queues** (rares peuplements ouverts) ; cLHS reproduit les **marginales ET les corrélations** des traits, couvre toute la gamme uniformément.
- **Pourquoi SÉLECTIONNER des vrais pixels (pas un LHS synthétique) :** pour rester **on-manifold** — les vrais canopées ont des combinaisons de traits physiquement cohérentes ; un LHS synthétique générerait des combos impossibles (LAI fort + cover faible) que MuSICA simulerait mais qui n'existent pas.
- **Pourquoi 100/cluster :** assez pour les distributions intra-archétype + le bootstrap.
- **Assumer le plancher fCover 0.5 :** sous 0.5 les hypothèses sol/sous-étage de MuSICA et le rétrieval deviennent peu fiables → **ça exclut le régime ouvert** (limite honnête, visible sur Fig 1c où des loggers tombent sous le plancher).

### 9. MuSICA v3.2.3, couplage ABL itératif (« yoyo »)
- **Pourquoi itératif :** l'analyse de sensibilité **change la structure de canopée**, ce qui change les flux vers l'atmosphère **au-dessus**. Un forçage fixe appliquerait la **même** météo au-dessus à un canopée dense et à un ouvert = **physiquement incohérent** (l'argument exact de Bouwen). Le yoyo laisse l'air au-dessus répondre.
- **Pourquoi c'est immatériel ici :** testé BLH 100→5000 m → **ΔTmax inchangé**, parce que le ΔTmax à 1 m est fixé par le bilan radiatif/turbulent **dans** la canopée, pas par le sommet de la couche limite. → représentation plus complète, mais **robustesse, pas driver** (donc on ne s'en sert pas pour se différencier de Bouwen).

### 10. Forçage : station CHS41 + ERA5 + MERRA-2
- **Pourquoi Tair = station et pas ERA5 :** ERA5-Land a un **biais chaud à structure diurne** vs la station (jusqu'à **+3.6 °C** tôt le matin) ; l'utiliser hériterait ce biais dans l'offset. La station = le vrai macroclimat local.
- **Pourquoi ERA5 pour rad/vent/pression :** la station ne les mesure pas ; ERA5 = meilleure source grillée.
- **Pourquoi référencer l'offset à la MÊME station :** ainsi forçage et référence **partagent la météo** → le biais chaud +0.85 °C est **généré dans le schéma de canopée, pas hérité du forçage**. (Défense clé du plafond de verre.)

### 11. Correction de vent (log neutre, dans le forçage)
- **Pourquoi corriger :** ERA5 donne le vent 10 m **en terrain ouvert** ; MuSICA a besoin du vent **juste au-dessus de chaque canopée**. Profil log neutre (déplacement d = 0.7h, rugosité z0 = 0.1h) → le vent brut 10 m est **~2.3× trop fort** au-dessus des grands canopées.
- **Pourquoi offline dans le forçage et figé à Hmax de base :** transparent, par placette, et **figé** pendant les perturbations → perturber Hmax **ne confond pas** le levier hauteur avec le vent.
- **Pourquoi la garder alors qu'elle dégrade le fit (0.52→0.42) :** c'est l'input **physiquement correct** ; la dégradation **expose** l'erreur structurelle du modèle avec de meilleurs inputs — elle ne dit pas que la correction est fausse.

### 12. Sensibilité = ±1 SD on-manifold, par placette
- **Pourquoi ±1 SD intra-archétype :** mesurer la leverage **réalisée** de chaque trait au point de fonctionnement réel, pas un changement unitaire abstrait. ±1 SD = la variation que le trait montre **vraiment** dans l'archétype.
- **Pourquoi par placette et pas sur les profils moyens :** **Jensen** — MuSICA est non-linéaire, f(canopée moyen) ≠ moyenne de f(canopées) ; perturber chaque vraie placette et agréger est **non biaisé**.
- **Pourquoi partenaires corrélés tenus réels :** isole l'effet d'un trait **sur le manifold** (interventionnel mais on-manifold), évite les chimères hors-manifold.
- **Pourquoi le profil en swap réel-vs-uniforme brut (pas ±SD, pas (1−VCI)) :** le profil est **fonctionnel**, on ne le bouge pas d'un « SD » scalaire ; l'intervention naturelle = « remplacer la forme réelle par l'uniforme à LAI/Hmax fixés » = **changement maximal = borne haute**. Pas rescalé (1−VCI) car **un VCI = une infinité de profils**.
- **Pourquoi différence centrale symétrique :** réduit le biais de courbure, capte la pente locale au point de fonctionnement.

### 13. Extraction des métriques (ΔTmax à 1 m, slope)
- **Pourquoi 1 m fixe :** les loggers sont à 1 m → comparer du comparable. Réf station à 1.5 m = hauteur météo standard (et Gril) ; le décalage 0.5 m est petit et connu.
- **Pourquoi ΔTmax (offset du max journalier) :** l'extrême est ce qui compte écologiquement (stress thermique) et c'est là que le buffering est le plus grand.
- **Pourquoi aussi la slope :** sépare « à quel point ça tamponne » (amplitude) de « position sur le gradient » (slope) ; le modèle récupère mieux la slope (r 0.70) que l'amplitude (r 0.42).

### 14. Validation (53 HOBO)
- **Pourquoi des runs MuSICA dédiées par logger (pas interpoler les cLHS) :** obs-moins-modèle défini au **vrai canopée** de chaque logger → résidu propre, pas une interpolation.
- **Pourquoi forward inclusion :** montre la **valeur incrémentale** de chaque trait ajouté dans un ordre fixe → le profil ajoute le plus petit incrément.
- **Pourquoi la régression nested model-free :** corroborer **sans MuSICA** → la domination de la **quantité** ne dépend pas du modèle. (**ATTENTION :** elle teste quantité vs métriques **scalaires** Hmax/VCI, PAS l'arrangement — VCI ≠ forme.)
- **Pourquoi la décomposition du résidu :** vérifier que le résidu est structuré par la **fermeture** (ouverture), pas par de la structure verticale non modélisée ou du terrain → l'échec est la limite 1-D connue.

## B. FIGURES — mécanique détaillée
Pour chacune : **Données · Comment c'est construit · Ce que ça montre · Pourquoi cette forme · À assumer.**

### Fig 1 — Site, typologie, distribution des traits  (`c1_fig1_gril_native.R`)
- **Données :** raster de clusters natif 20 m ; 400 pixels cLHS ; 53 loggers ; 4 placettes représentatives ; nuages LAS bruts.
- **Comment :** 3 blocs assemblés en patchwork. **(a)** carte : raster d'archétypes (trous intérieurs comblés par vote majoritaire des voisins, 10 itérations — les gaps sous-canopée/routes ; l'extérieur non-forêt est laissé intact), points cLHS, croix HOBO, cercles a–d, barre d'échelle 1 km + flèche nord. **(b)** pour chaque archétype : nuage LAS vu de côté (coloré par hauteur, rampe brun→vert) **aligné sur le même axe 0–40 m** que son profil LAD, les deux dérivés **du MÊME clip LAS** (r=12 m, MacArthur-Horn dz=1, k=0.5) → garantit que nuage et profil décrivent le même canopée ; Hmax/LAI/VCI annotés dans le coin. **(c)** demi-violons (ggdist) des 4 traits sur les 400 cLHS + points bruts + médiane, avec les **53 loggers en croix noires**.
- **Ce que ça montre :** la typologie couvre tout le gradient structurel ; les loggers échantillonnent le même espace de traits **sauf à l'ouvert** (plusieurs croix sous le plancher fCover 0.5 et à faible VCI).
- **Pourquoi cette forme :** style « Gril » (carte + coupes) pour ancrer visuellement ; nuage+profil côte à côte = le lecteur voit que le profil LAD *est* le nuage.
- **À assumer :** les 4 placettes représentatives sont **illustratives** (elles précèdent l'échantillon natif), d'où le recalcul depuis le LAS plutôt qu'une ligne de l'échantillon.

### Fig 2 — Schéma du pipeline (fait main)
- **Ce que ça montre :** ALS → 4 traits → FPCA/typologie → cLHS 400 → MuSICA (v3.2.3 iter, forçage station) → attribution ±1 SD → validation 53 loggers.
- **Pourquoi :** une image qui rend lisible l'enchaînement des choix ; les flèches y sont légitimes (notation de schéma).

### Fig 3 — Attribution : sensibilité des traits par archétype  (`scripts/c1_fig3_attribution_native.R`)
- **Données :** `metrics6_native20` (les mêmes valeurs par placette que le bootstrap / Table H1). Aucune simulation nouvelle.
- **Comment :** sensibilité **par placette** = différence centrale `(add − rem)/2` pour LAI/fCover/Hmax ; profil = `dT_LAD` **brut** (swap réel−uniforme). Pour la métrique slope, on passe en Δlog(β) ≈ Δβ/β_base. Violons `scale="width"` = **distribution des 100 placettes** de chaque archétype, dodgés par cluster, + boîtes blanches (médiane + IQR). 2 lignes = ΔTmax et Δlog(slope).
- **Ce que ça montre :** LAI domine et **se renforce vers le dense** ; fCover second ; Hmax petit (max en P1) ; le profil (borne haute) est le plus gros levier en **P2/P3** seulement.
- **Pourquoi cette forme :** violons et non barres → on voit la **dispersion inter-placettes**, pas juste une moyenne (demande du comité) ; garder le **signe** rend visibles les inversions.
- **À assumer :** les **assises diffèrent par construction** (scalaires = ±1 SD réalisé ; profil = swap complet = borne haute) → non interchangeables ; P1 lu avec prudence (le modèle n'y reproduit pas le terrain).

### Fig 4 — Le levier est fonction du point de fonctionnement, pas du cluster  (`scripts/c1_operating_point_main.R`)
- **Données :** les 400 placettes, labels de cluster **jetés** dans le modèle.
- **Comment :** 4 panneaux (a) levier LAI, (b) levier Hmax, (c) levier fCover (plancher 0.5 exclu), (d) swap profil — chacun tracé **contre SON propre point de fonctionnement** (LAI de base, Hmax de base, fCover de base ; le profil n'ayant pas d'axe scalaire, il est tracé vs densité). Points colorés par cluster **a posteriori**, lissage loess robuste. Annotation en coin = **partition de variance** : R² du levier vs point de fonctionnement **continu** contre R² vs **étiquette de cluster** (0.88 vs 0.72 ; 0.75 vs 0.28 ; 0.66 vs 0.21 ; 0.63 vs 0.34).
- **Ce que ça montre :** chaque levier est une **courbe lisse et continue** ; les archétypes tombent comme des bacs le long de cette courbe et **traversent les frontières sans marche**. Le panneau (d) montre le swap profil **négatif au milieu du gradient puis repassant positif** au fort LAI (l'inversion de signe en P4).
- **Pourquoi cette forme :** c'est **la** réponse à « vos clusters fabriquent le résultat ». La **continuité elle-même est la preuve** ; la partition de variance la chiffre.
- **À assumer :** la partition est **illustrative** (l'étiquette de cluster est un binning grossier du même point de fonctionnement) — on ne la vend pas comme un test.
- **NB :** une variante séparée à 8 panneaux (`c1_operating_point_response.R`) montre **+1 SD ET −1 SD ET la valeur réelle non perturbée** (gris pointillé) — faite exprès pour la question du comité « il n'y a que +1 SD ? la valeur initiale est-elle incluse ? ». À dégainer si on te repose la question.

### Fig 5 — Plafond de verre, logger par logger  (`c1_dumbbell_by_cluster.R`)
- **Données :** 53 HOBO horaires vs 53 **runs MuSICA dédiées** ; macro = station 1.5 m.
- **Comment :** pour chaque logger, on calcule la même paire de métriques des deux côtés (ΔTmax = max journalier micro − max macro ; slope = régression horaire micro~macro). Puis **dumbbell** : cercle = observé, triangle = simulé, **segment = l'erreur du modèle**, couleur = archétype, loggers **triés par la valeur observée**. 4 variantes = 2 métriques × (tous les jours / 10 % les plus chauds).
- **Ce que ça montre :** **(a)** en ΔTmax le modèle **comprime toute la gamme vers zéro** (sous-tamponne le dense, ne reproduit pas le réchauffement des ouverts) ; **(b)** en slope les écarts sont bien plus petits → il **place** les placettes sur le gradient sans en restituer l'amplitude.
- **Pourquoi cette forme :** le dumbbell rend l'erreur **par placette** visible d'un coup d'œil (un scatter la moyennerait) ; trier par l'observé fait apparaître la compression comme un éventail qui se referme.

### Fig 6 — Corroboration model-free  (`scripts/c1_fig6_obs_nested_native.R`)
- **Données :** traits natifs 20 m **extraits au pixel de chaque logger** + buffering **observé** (dTmax_obs). **Aucun MuSICA.**
- **Comment :** 4 régressions linéaires emboîtées sur n=53 : LAI seul ; **Quantité** (LAI+fCover) ; **Structure** (Hmax+VCI) ; **Complet** (les 4). R² calculés **dynamiquement** dans le script, barres + ΔR² annotés.
- **Ce que ça montre :** Quantité seule = **0.80** = le modèle complet ; ajouter la structure **n'ajoute rien** (ΔR² < 0.01) ; ajouter la quantité à la structure ajoute **+0.13**.
- **Pourquoi cette forme :** preuve **sans modèle de process et sans partitionner** la variance partagée de traits colinéaires → le résultat central ne repose pas sur le skill de MuSICA.
- **À assumer (le point trouvé par le harsh review) :** le bras « structure » est **Hmax + VCI**, or **VCI est un descripteur de densité, pas de forme** → cette figure **ne teste PAS l'arrangement vertical**. Le seul test model-free de la forme est la **top-heaviness** (Fig F3), **non significative** (partial r = −0.23, p = 0.10). Le texte le dit maintenant explicitement.

### Annexe C2 / C3 — Gradients verticaux  (`make_vertical_profiles_native.R`, `make_trait_vertical_gradient_perplot.R`)
- **Données :** coalitions forward natives déjà simulées — **1111 = profil réel**, **1110 = profil uniforme** (LAI/Hmax/fCover réels). 53 pixels-loggers **groupés par archétype**. **Aucune simulation nouvelle.**
- **Comment (C2) :** moyenne JJAS, fenêtre diurne 10–16 h, par **niveau** (15 couches d'air) ; RH et VPD **dérivés** de Tair_z/wair_z (Tetens) car MuSICA ne les sort pas ; **deux vues** : normalisée z/h_canopée (haut) et **mètres absolus** (bas, avec le sommet propre à chaque archétype).
- **Ce que ça montre :** l'écart réel−uniforme est **petit partout** (0.05–0.20 °C au niveau sous-canopée) face à un écart **inter-archétypes** bien plus grand (~0.55 °C) piloté par la densité → H1 rendue visible.
- **C3 :** effet marginal **de chaque trait** sur le gradient. Le trait dominant **change avec la densité** : P1 → LAI (−0.35 °C), P2/P3 → **profil** (−0.13 / −0.19), P4 → fCover (−0.57).
- **À assumer :** diagnostic **interne au modèle** — un seul niveau (≈1 m) est ancré au terrain ; et la ligne de base de C3 est le **canopée moyen global**, pas la base on-manifold par placette de la Fig 3 (donc les signes se lisent comme des écarts au canopée moyen).

- **Provenance complète :** `chapter1/ledger/make_article_figure_set.R` = registre des **25 figures citées**, chacune → son script + un drapeau « rafraîchie par un re-run » (15/25).

## C. HEADLINES à asserter
1. **La densité horizontale (quantité + fermeture) gouverne le buffering réalisé** ; levier ↑ avec la densité (LAI −0.33 °C/SD en P4).
2. **Profil vertical = puissant mais non réalisé** (borne haute, domine P2/P3, mais quasi invariant dans les vrais peuplements) → **H1 non réalisée**.
3. **Modèle : rang oui (r=0.42), amplitude non** (~11 %, +0.85 °C), **échoue à l'ouvert** (r=−0.37) = plafond de verre.
4. **La domination de la quantité est model-indépendante** (régression model-free + tout Beer-Lambert) ; **le négatif profil, lui, repose sur le modèle** (cohérent terrain, non confirmé indépendamment).

## D. LIMITES à OWN (les dire en premier)
- **r=0.42 faible** → attribution lue **en canopée fermée** où le modèle marche = **limite de portée assumée**.
- **LAI vs cover non séparable** (r 0.91 loggers / 0.95 cLHS) → dire **« axe quantité »**, pas « LAI devant cover ».
- **Négatif H1 sous-puissant** (forward r 0.39→0.42 sur n=53, dans le bruit) → « pas d'incrément **détectable** », pas « n'existe pas ».
- **Mécaniste (R²≈0.18) < statistique de Gril (0.91) même site** → **complémentaires** : statistique = *où*, mécaniste = *quel trait*.
- **Négatif profil :** le terrain ne peut ni confirmer ni réfuter (Annexe E) ; seul test model-free de la forme (Fig F3) = **non-significatif** (partial r=−0.23, p=0.10).
- **Un site, une espèce, un été feuillé.**
- **k=0.5** conditionne l'amplitude du levier quantité ; **BLH** corrompue dans le forçage mais **prouvée immatérielle** (test 100–5000 m).

## E. PIÈGES DE RÉUNION (réponses courtes)
- *« Pourquoi pas Shapley / GAMM ? »* → colinéarité (r 0.91) + concurvity ; le GAMM V1 était une **boîte noire** (effet « Date »). La sensibilité on-manifold évite les **chimères hors-manifold** d'une attribution interventionnelle.
- *« Vos archétypes sont-ils réels ? »* → non, commodité (silhouette 0.27) ; la densité-dépendance est **continue** (Fig 4).
- *« Pourquoi croire le mécaniste s'il perd contre le statistique ? »* → il n'a pas le même job : Gril prédit **où**, MuSICA explique **quel trait** ; et le résultat quantité tient **aussi sans MuSICA** (Fig 6).
- *« Et Sentinel-2 ? »* → puisque la 3D verticale n'ajoute rien au **réalisé**, un proxy 2-D « quantité » pourrait suffire → **pont Chap 2** (illusion optique / saturation).
- *« Le biais +0.85 °C invalide tout ? »* → non : un **offset fixe** s'annule dans la différence ±SD → **le rang des leviers survit** ; seul le *gain* est affecté, et on ne lit l'attribution que là où le modèle est fidèle.
