# Discussion Générale — carte des matériaux (passe Fable, 2026-08-13)

Synthèse de trois relectures complètes (intro générale + Ch1 canonique native20,
Ch2 r3 août 2026, Ch3 standalone) croisées avec les perspectives déjà construites
(GEDI/Annexe G, test 2022, fenêtre vernale, nuit, design loggers). Document de
travail FR ; la Discussion elle-même sera EN.

---

## 1. Les promesses de l'intro que la Discussion DOIT honorer

- **RQ1** (quantité vs arrangement vertical, attribution mécaniste malgré la
  collinéarité) — fermée par Ch1, avec le résidu avoué : « What the design cannot
  provide is a direct empirical test of the attribution itself ».
- **RQ2** (réconcilier S2 et LiDAR : consistance inter-capteurs, PAS accuracy) —
  fermée par Ch2 ; les deux caveats promis (domaine restreint ; « better agreement
  in leaf area does not automatically propagate into better simulated
  microclimate ») sont exactement ce que Ch3 arbitre.
- **RQ3** (quel forçage, dans quel régime, déployable hors-LiDAR) — fermée par Ch3
  (fusion gatée hauteur), avec l'hypothèse d'entrée confirmée (« no single sensor
  wins everywhere »).
- **Le séquençage** : l'intro affirme que les 3 gaps « have to be addressed in
  order » — la Discussion doit montrer que chaque chapitre a bien rétréci la
  question du suivant (c'est le cas, cf. §2.1 ci-dessous).
- **La doctrine triple-axe** (ranking / dispersion / biais rapportés séparément)
  — à réutiliser comme grille de lecture dans toute la Discussion.
- **La frontière nommée** : « This thesis does not produce a wall-to-wall
  microclimate map. It establishes what a map would have to be built on…
  Naming that boundary is not a concession made late, it is the reason the
  General Discussion turns where it does. » → la Discussion doit s'ouvrir là.
- **Les deux avenues explicitement réservées à la Discussion** (§1.7 intro) :
  (a) LiDAR spatial (GEDI) — ⚠️ l'intro le DÉGRADE (« sampling geometry and
  geolocation uncertainty nonetheless keep it short, for now, of what plot-scale
  mechanistic microclimate modeling requires ») ; (b) modèles 3D/bords
  (ARPS-MuSICA [REF: Grulois 2026 — placeholder NON résolu dans l'intro],
  ForEdgeClim).

## 2. Les fils de synthèse inter-chapitres (le squelette)

### 2.1 LE fil unificateur : « l'exactitude de récupération n'est pas la fidélité fonctionnelle »
La plus belle phrase transversale est déjà dans Ch1 : « The leaf area that is most
accurate for retrieval is therefore not the leaf area that best couples to the
microclimate. » La chaîne complète :
- Ch1 Appendix F : une réduction UNIFORME de 23 % du LAI efface le skill
  (r 0.50→0.19 ; amplitude 15→3 % ; biais +0.98→+1.33) — et le k = 0.65 issu du
  travail inter-capteurs (Ch2) fait précisément cela.
- Ch2 : le LUT « opt » réduit RMSE de 1.2 et biais de 1.5 m²/m² (le meilleur
  produit de RÉCUPÉRATION) ;
- Ch3 : ce même produit est le PIRE forçage microclimat (pooled R² 0.00,
  cor espace-LAI +0.19, signe inversé) et la troncature d_opt = le pire biais
  (+1.14 °C).
→ Section de Discussion à part entière : *ce qu'on optimise détermine ce qu'on
obtient* ; la consistance inter-capteurs (Ch2) et le couplage fonctionnel (Ch3)
sont deux fonctions-objectif différentes. Le « anti-corollaire » de l'intro
(1.5.2) l'avait promis.

### 2.2 Un seul gradient, trois manifestations
Le gradient de densité organise TOUTE la thèse : Ch1 — le modèle n'a de skill
qu'en canopée fermée (dense r = 0.69, ouverte r = +0.20 n.s.) et le levier
quantité est max en dense (−0.26 °C/0.5 LAI en P4) ; Ch2 — la saturation S2 est
un phénomène de canopée dense (63 % du LAI sous d_opt en dense vs 19 % en
ouvert) ; Ch3 — le renversement stratifié (dense LiDAR 0.60/S2 0.02 ; ouvert S2
0.71/LiDAR ~0). → La Discussion peut superposer les trois seuils : transition
microclimat LAI ≈ 2–3 (Ch3 segmented), crossover capteurs 3.86, saturation
optique 5–6 : trois lignes sur le même axe, pas une coïncidence mais une
géométrie (là où le buffering se joue, l'optique ne voit plus).

### 2.3 « LiDAR requis » ≠ « la 3D requise » — désambiguïser
Ch1 borne le détail vertical (profil = potent-but-unrealized, ΔR² = 0.000 ;
Hmax ≤ 0.004 °C m⁻¹ comme LEVIER) ; Ch3 rend le LiDAR indispensable en dense —
mais pour sa MAGNITUDE non-saturée, et la hauteur y sert de CLASSIFICATEUR de
régime (gate 17.8 m), pas de levier physique. Sans ce distinguo explicite, un
lecteur colle Ch1 (« vertical detail unnecessary ») contre Ch3 (« LiDAR
required ») ; l'intro a déjà la formule-pont (« the forcing itself still draws
canopy height and profile shape from LiDAR in both regimes ») — la reprendre.

### 2.4 L'épistémologie assumée : pas de vérité terrain, nulle part
Aucun LAI terrain (Ch1, Ch2, Ch3), pas de « ground truth » : la thèse tient par
un triangle — attribution mécaniste (Ch1), consistance inter-capteurs (Ch2),
arbitrage par 53 loggers (Ch3). L'argument Woodgate (même les instruments
terrain divergent) est déjà dans Ch2 r3 : en faire la position épistémique de
la thèse entière, pas une excuse locale. Y rattacher la circularité bornée
(≤17 % prior, scoring non-borné) et le k-dependence (« only the variable
ranking is claimed to carry »).

### 2.5 Le plafond du modèle — ce que AUCUN forçage ne répare
À rassembler en une sous-section : biais chaud +0.98 °C (scheme-level, survit
aux corrections vent/scan), ~15 % d'amplitude capturée (Ch1), SDrec ≤ 0.28
(Ch3), biais nuit +1.5–1.7 °C (Ch3), skill dégradé les jours les plus chauds
(r = 0.35, Ch1) et biais accru en 2022 (+1.16, graine A). Message : la marge de
progrès est passée du FORÇAGE au MODÈLE (1-D, sans bords) → transition naturelle
vers l'avenue 3D/ARPS-MuSICA/ForEdgeClim.

## 3. Liste rouge — les tensions à désamorcer (phrase par phrase)

1. **Ch3 App. A/B : « the value of the fusion is spatial…, not temporal » et
   « the deficit is structural, not temporal »** — catégoriels ; à re-gloser en
   énoncés *fenêtre-été* sinon le gain d'automne GEDI (Annexe G, §G.4) les
   contredit. Le corps de Ch3 est déjà prudent (« not testable beyond summer
   with these data ») : citer le corps, pas les annexes.
2. **Intro §1.7 dégrade GEDI** : notre perspective doit expliciter que le rôle
   démontré est une exigence PLUS FAIBLE que celle que l'intro écarte (timing
   saisonnier au niveau site + magnitude moyenne, PAS structure 3D par plot) —
   l'Annexe G le dit déjà (pattern footprint-bound) ; ajouter la phrase-pivot
   dans la sous-section perspectives.
3. **Ch2 : « systematic underestimation of LAI_S2 » (H1) et « underestimated
   LAI beyond these values » (intro Ch2)** — ne JAMAIS reprendre tel quel : à
   k = 0.65 les biais moyens sont −0.32/+1.11/+0.03 (S2 comparable ou
   AU-DESSUS) ; réutiliser les formulations sûres de Ch2 même (« mechanically
   compresses the LAI_S2 range and depresses the slope »).
4. **Ch2 §4.5 : « LAI_S2_opt may serve as an operational input… forest
   microclimate modelling is one such application »** — clash frontal avec Ch3
   (opt = pire forçage). La Discussion DOIT traiter cette phrase de front :
   c'est l'anti-corollaire prédit par l'intro, vérifié par Ch3. Ne pas la fuir,
   la retourner (c'est un des plus beaux résultats de la thèse).
5. **Hmax** : négligeable comme levier (Ch1) vs gate de fusion (Ch3) — toujours
   dire « height as classifier, not as lever ».
6. **ΔTmin** : Ch3 dit « the leaf-area field ranks [night cooling] only
   weakly » alors que ΔTmin poolé Ch3 = 0.68–0.73 et notre graine donne
   cor(ΔTmin, LAI) = +0.83. Réconcilier : le signe (dense = nuits plus CHAUDES),
   la mécanique intro (piégeage longwave) et le garde-fou Ch1 (« a distinct
   mechanism we do not examine ») → présenter la nuit comme territoire NOUVEAU,
   compression bilatérale d'amplitude, avec le double tranchant écologique
   (nuits chaudes ≠ protectrices, déficit hydrique, intro §1.1.1).
7. **Le piège de collinéarité** : les corrélations observationnelles des graines
   (2022 r = −0.88 ; nuit +0.83 ; vernal +0.48) sont des CHECKS DE COHÉRENCE de
   l'attribution mécaniste Ch1, pas des attributions — le dire une fois,
   clairement, sinon la Discussion ré-ouvre le piège que la thèse a construit
   pour éviter.
8. **2022** : le « deepening » (−0.73→−1.22) va contre la littérature
   sécheresse citée par l'intro (buffering s'affaiblit quand les sols sèchent)
   → le formuler comme observation site/année (chênaie de plaine, nappe ?) qui
   NUANCE la littérature, pas qui la réfute ; et rappeler que le skill absolu
   dégrade en chaud extrême (Ch1 r = 0.35 ; biais 2022 +1.16).
9. **Cohérences numériques** : d_opt Aigoual = 8 m (Table 3 r3) mais « 7 m »
   deux fois dans Ch2 §4.2 → standardiser (8 m ; « 6–10 m selon les sites, 7 m
   combiné ») ; [REF: Grulois 2026] non résolu dans l'intro ; r = 0.50 (jamais
   0.46 du forward-inclusion) ; « quantité » = LAI+cover ensemble, jamais
   classés l'un contre l'autre.

## 4. Les perspectives, hiérarchisées (avec état des preuves)

**Étage 1 — seedées et chiffrées (prêtes, documents compilés)** :
- GEDI = dimension temporelle LiDAR (Annexe G complète ; produit composé
  FUSION_GEDIMAX = meilleur forçage évalué ; EOS +25–34 j ; réplication
  Cotrina). Positionner contre l'intro (cf. §3.2).
- 2022 / microrefugia : couplage tient (−0.88), buffering se creuse, fusion
  transfère hors-année ; caveat biais modèle croissant.
- Fenêtre vernale en degrés (+0.48) + nul hiver→SOS (0.05) : la flèche mesurée
  va de la phéno vers le microclimat — paire propre avec Wu 2024.
- Nuit : compression bilatérale (jour −0.885 / nuit +0.83) — candidat à monter
  dans la SYNTHÈSE plutôt qu'en perspective (résultat de fond).
- Design : 15–20 loggers suffisent pour établir le couplage pente~LAI
  (bootstrap) — une phrase d'implication pratique pour la réplication.

**Étage 2 — hooks des chapitres, non testés (paragraphes courts)** :
- Modèles 3D/bords : ARPS-MuSICA [REF Grulois], ForEdgeClim — l'avenue réservée
  par l'intro ET le débouché du §2.5 (plafond modèle).
- Assimilation séquentielle S2 (EnKF/particule) — « the one untested fusion
  family » (Ch3 §4.3) ; se marie naturellement avec l'arbitre GEDI.
- Sous-bois < 2 m : re-traiter les nuages de points (hook Ch3, mécanisme ouvert
  du skill S2 en ouvert).
- Champagnes LAI terrain standardisées (hook Ch2 #3) — fermer l'axe accuracy.
- Conifères / mixtes / leaf-off (hook Ch2 #2 ; Aigoual = cas d'école où gate et
  GEDI échouent tous deux).
- Hyperspectral (CHIME, S2-NG) et S1 SAR tout-temps (Soudani 2021) pour le
  timing des épaules.
- Sylviculture : buffering et récolte (hook Ch1, « directional corollary »).

**Étage 3 — à ne PAS promettre** : carte wall-to-wall microclimat (frontière
nommée par l'intro) ; toute claim d'accuracy absolue ; toute attribution depuis
des corrélations observationnelles.

## 5. Squelette proposé pour la Discussion Générale

1. Retour sur les trois questions (RQ1→RQ3, une synthèse par chapitre, le
   séquençage démontré).
2. Un gradient, trois lignes (2.2) + le fil « accuracy ≠ fidélité » (2.1) —
   les deux sections de fond.
3. « LiDAR requis, mais pas pour sa 3D » (2.3) + épistémologie sans vérité
   terrain (2.4).
4. Le plafond du modèle (2.5) — pivot vers les perspectives.
5. Perspectives étage 1 (GEDI temporel en tête, puis 2022/refugia, vernal,
   nuit) puis étage 2, fermées par la frontière assumée (pas de carte — encore).
