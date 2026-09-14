# Synthèse de la passe croisée Fable — Discussion Générale vs reste de la thèse

*2026-08-18. Cinq revues indépendantes (intro, Ch1, Ch2, Ch3, Annexe G/perspectives),
rapports détaillés : `disc_crosscheck_{intro,ch1,ch2,ch3,annexG}.md`. Sauvegardes
pré-édition dans `_bak_2026-08-18/`. Ce fichier = ce qui a été trouvé, ce qui a été
appliqué, ce qui reste à arbitrer.*

## 1. Verdict global

L'arc tient : chaque RQ de l'intro trouve sa réponse, la logique de gradient (D.2)
et la posture épistémique (D.5) sont celles de l'intro, tous les chiffres Ch1/Ch2/Ch3
cités par la Discussion sont exacts **contre les versions canoniques**. Les problèmes
sont de trois natures : (a) deux erreurs de données côté Annexe G (SOS phénologique
calculé sur la série d'été ; comparaison 2022 mal appariée) ; (b) une dizaine de
glissements de cadrage entre la Discussion et les chapitres (caveat attribué au Ch2 qui
ne l'énonce pas, « best retrieval », rôle de S2 « summer timing », automne « non
testable » dans le Ch3, hétérogénéité intra-pixel absente, Bouwen jamais refermé) ;
(c) deux contradictions dont la faute est **côté intro**, pas côté Discussion.

## 2. Ce qui a été corrigé aujourd'hui (données)

- **SOS par placette (D.7.3 / §G.7)** — l'ancien `sos_s2` (Table23c/24c) était un
  premier franchissement brut sur `tb$atbd`, la série S2 *d'été* du Ch3 extrapolée à
  plat sur l'année : 15–20 placettes « débourraient » en janvier. Recalcul par
  `c4_disc_seeds3_sos_fix.R` : TRS50 lissé (whit2, semainier, dernier franchissement
  avant le pic) sur la série pleine-année `annual` (47 placettes) → SOS DOY 106–134.
  Nouveaux chiffres : printemps r = +0.55, partiel|LAI +0.40, partiel|LAI+Hmax +0.33
  (n = 33) ; hiver r = +0.40, partiel +0.19 (n = 35). Sur `naive_s2` (53) : +0.51/+0.30
  et +0.34/+0.03. Le résultat vernal survit affaibli ; l'asymétrie hiver/printemps
  survit sur les partiels. Tables : `Table23c_c4_pheno_loop_v2_annual.csv`,
  `Table24c_c4_winter_sos_v2_annual.csv`, `Table25_c4_sos_fix_summary_{annual,naive_s2}.csv`.
  Note : SOS corrèle à −0.67 avec LAI_ALS (denses = plus tôt) → ce sont les partiels
  qu'il faut lire ; c'est écrit dans §G.7.
- **2022 (D.7.2 / §G.6 / Table G4)** — le texte opposait « tous jours 2021 » (−0.73)
  à « jours chauds 2022 » (−1.22). À fenêtre égale (p90) : −1.13 vs −1.22, r(slope,LAI)
  −0.882 vs −0.881, pente 0.86 vs 0.93. Table G4 passe à 4 colonnes avec la pente ;
  §G.6 perd « buffers most exactly when it matters most » (déjà rétracté côté
  Discussion, l'annexe traînait derrière).

## 3. Ce qui a été corrigé aujourd'hui (texte de la Discussion, 31 remplacements)

D.1 : « ranks far better than reproduces » → « recovers the ordering… far better than
the amplitude » ; 23 % défini comme ratio k 0.5/0.65 ; paragraphe Ch2 réécrit
(plafond 6–7 + overshoot sparse, k glosé, d_opt = proxy empirique, **hétérogénéité
intra-pixel ajoutée**, le Ch2 *propose* le produit, il ne « flagge » pas le caveat ; « the
caveat became the result » → « the offer became the test ») ; « among the products
tested » ; « near zero » scopé à l'offset ; switch « statistically tied with S2 alone on
the pooled ranking » ; forme opérationnelle : seul le gate est LiDAR-free ; ajout du
poids inégal des deux moitiés (dense = model-free r −0.74 ; open = via le modèle).
D.2 : « crossover ≈ 3.9 » → « reversal, split médian 3.86, stable 2.0–4.5 » ; plafond
6–7 du Ch2 + « près de 4 » attribué au pool Blois k = 0.5 ; conventions k nommées ;
63 % scopé aux placettes denses de Blois ; « per unit leaf area » → +0.41 °C de biais ;
la tension Ch1 (échec structurel en ouvert) vs Ch3 (S2 classe l'ouvert à 0.71) est
nommée et lue comme association proxy ; « genuine resolution » → « division of labor ».
D.3 : titre → « Inter-sensor consistency is not functional fidelity » ; « best retrieval »
→ « most inter-sensor-consistent », gain mesuré contre le LiDAR tronqué, prior Blois =
LAI_ALS_dopt ; « worst forcing » → « no better than the untuned retrieval (0.07 vs
0.01) » ; « the route not the tuning » → contrainte en amont du tuning (raw reflectance
R² 0.35).
D.4 : increment profil chiffré (+0.026, CI) ; hauteur = borne inférieure (voie
aérodynamique fixée) ; direction H2 énoncée ; « excellent classifier » → « serviceable
(85 %) » ; « stop resolving beyond 5–6 » → ranking perdu avant le plafond ; **paragraphe
Bouwen ajouté** (profil bêta ≠ cause de l'under-buffering ; la colonne l'est).
D.5 : « no field LAI measured » → « used… none concurrent available » ; loggers
« independent of both the model and the LAI retrievals » ; résidu Ch1 réassigné
(ranking corroboré model-free, magnitudes model-internal) ; circularité d_opt ajoutée ;
4e axe model-free ; phrase VPD.
D.6 : « at or below 0.28 » ; nuit attribuée au Ch3 seul ; direction non testée
(k < 0.5 / clumping) avant de conclure « margin shifted to the model » ; sub-pixel gaps
+ queue amplifiante ; « unshielded » → non-aspirated shields + 1 m vs 1.5 m ; nuit :
« warms the night » → fermé tient la référence, ouvert tombe 2–4 °C sous ; Ch3 nuit
recadré (R² poolé 0.68–0.73, le manque = amplitude/biais, pas ranking) ; macro = fichier
de forçage.
D.7.1 : « at one site » ; automne = archive au-delà du Ch3, 47 series plots ; « same
53 plots » → « same logger network » ; matrice capteurs : S2 = ranking ouvert + timing
leaf-on (nul sur le plateau, Ch3) ; footprint n = 11–18 ; Aigoual chiffré et renvoyé à
§G.9.
D.7.2 : fenêtres appariées (cf. §2).  D.7.3 : nouveaux chiffres (cf. §2).
D.7.4 : Grulois sans tiret, phrase dé-dupliquée de l'intro ; sub-2 m relié au 0.5 m du
Ch1 et au h_min 2–5 m du Ch2 ; loggers « within this network… untested elsewhere ».
D.8 : « airborne LiDAR wherever the canopy is dense » ; « a domain of validity » ; « first
two exploratory trials » ; « a network of 15–20 loggers ».

Annexe G (14 remplacements + légendes) : titre/légendes sans tiret ; « nine
acquisitions, eight composite points » ; EOS « about 25 days (34 with the spline
metric) » ; Table G3 rendue (ligne vide) ; ≤ 0.18 restreint aux temporels + statiques
0.32/0.49 nommés ; footprint 100–150 m n = 11–18 ; §G.6/Table G4 (cf. §2) ; §G.7
réécrit avec méthode + retrait explicite des anciens chiffres ; §G.8 nuit + forçage ;
loggers ±25–30 % vs ±11 % à n = 40, « within this network » ; §G.9 : Aigoual ajouté
(7 caveats), « 47 series plots », test hiver déclaré fait.

PDF recompilés : `discussion_generale_these_EN.pdf`, `annexe_GEDI_EN.pdf` ; 0 citation
non résolue ; `[REF: Grulois 2026]` toujours en attente de l'entrée Zotero.

## 4. Ce qui reste — À TOI d'arbitrer (hors périmètre Discussion)

### 4.1 Intro générale (2 sévères, côté intro)
- **I808-811 et I882-883** : « leaf quantity saturates… profile co-leads in dense
  canopy » — inversé par rapport au Ch1 canonique (levier quantité *se raidit* vers le
  dense ; profil = plus grand levier en P3, borne sup.) et « co-leads » est sur la
  liste rouge du skill. Réécriture proposée : `disc_crosscheck_intro.md` C-intro 1.
- **I899-909** : « outside summer the loggers do not exist » — faux (archive 2020-07 →
  2023-08) et contredit par D.7.1/D.7.2/D.7.3. Scoper le *chapitre*, pas l'archive :
  C-intro 2.
- Moindres : I864-866 « accuracy » → « consistency » (C-intro 3) ; annoncer le
  « field test » model-free du Ch1 (C-intro 4) ; un seul nom pour la pente
  (« buffering slope ») repris par la Discussion (C-intro 5) ; I906-908 peut rester tel
  quel maintenant que D.7.1 est reformulé.

### 4.2 Chapitre 1
- **`review/manuscript_chap1_EN.md` est la copie périmée du 25 juin** (r 0.93, +0.9,
  « Shapley-style », pas de k = 0.65). La canonique est
  `chapter1/manuscript/manuscript_chap1_EN_native20.md` (14 août). À mettre en
  quarantaine avant assemblage (je ne supprime rien sans ton feu vert).
- App. F du Ch1 dit k = 0.65 « more accurate » — contredit la posture D.5 / Ch2 r3.
  Suggestion : « optimized for inter-sensor consistency ».
- Écart +0.98 (Ch1) vs +0.73 °C (Ch3) de biais chaud pour un forçage LiDAR — deux
  fenêtres/étalons différents, une phrase de réconciliation quelque part serait bienvenue.

### 4.3 Chapitre 3
- **App. A** : « loggers cover the summer only… not testable beyond summer with these
  data » → « the summer window retained for this chapter » (l'automne est utilisé en
  Annexe G).
- §3.3 : si tu veux garder « per unit leaf area » quelque part, ajouter le test
  β_below ≈ β_top (il existe dans les tables projet, pas dans le manuscrit) ; sinon la
  Discussion s'en passe désormais.
- Vérifier que le jeu de figures de thèse tire de `figK_hybrid_ch2.csv` (0.573/0.536)
  et non de `figK_hybrid_switch.csv` (0.550/0.540).

### 4.4 Chapitre 2
- Rien de bloquant. Points que la Discussion reprend maintenant et que le Ch2 doit
  bien contenir en r3 : h_min 2–5 m à Blois, résidu +0.7 à Blois, prior Blois =
  LAI_ALS_dopt (Table A.16).

### 4.5 Fichier `disc_generale_perspectives_GEDI_EN.md`
- Recommandation des cinq revues : **archiver** (redondant avec D.7). Quatre phrases à
  porter dans D.7 si tu y tiens : accord SOS à 3 jours (déjà en D.7.1 via §G.3),
  séparation feuille/bois, caveat marcescence, complémentarité FORMS. Je ne l'ai pas
  supprimé.

### 4.6 Non vérifiable par les agents (inchangé)
- Représentativité de la Tmin macro du fichier pblh (cold pooling) — dit dans le texte.
- Deux jeux d'amplitudes GEDI dans le dossier (Table13 vs Table21b) ; l'annexe utilise
  Table13 ; annoter la provenance avant soutenance.

## 5. Chiffres re-vérifiés dans cette passe (tous OK contre les sources canoniques)
Ch1 native20 : r 0.50, +0.98, ~15 %, dense 0.69/open +0.20, −0.26/0.5 LAI en P4,
≤ 0.004 °C m⁻¹, +0.026 [−0.007, +0.070], 23 % → r 0.19. Ch2 : d_opt 6–10/7 m, tops
17–32 m, −1.2 m² m⁻², ≤ 17 %, 0.25–0.80, k 0.5→0.65 = −23 %. Ch3 : 0.60/0.02,
0.71, 0.573→0.536, +0.67/+0.71, 1.60/1.62, 17.8 m, 85 %, 63 %, +1.14/+0.73, SDrec
≤ 0.28, r −0.92, +1.5–1.7 °C nuit. Annexe G : Table G3 (40 cellules) vs Table20,
G2 vs 15b/18d/22b, 1.609/0.464/1.568, Table23b, nuit +0.83/−0.89, Table24d, EOS 25 j,
SOS 122/125.

## 6. Suite « change en fct » (2026-08-18, après-midi) — §4 appliqué

- **Intro** (8 édits, sauvegarde `_bak_2026-08-18/`) : I808-811 réécrit (levier quantité se
  raidit, profil = plus grand levier en canopée intermédiaire, borne sup., plus de « co-leads ») ;
  I882-883 idem ; I899-909 scope le chapitre et nomme l'archive 2022 ; I864-866 « accuracy » →
  « consistency » ; check model-free annoncé (I815) ; « buffering slope » introduit à I346 et
  substitué à « thermal-coupling slope » (I855). PDF recompilé, 0 citation cassée.
- **Ch1 native20** App. F : « more accurate » / « retrieval accuracy » → « maximizes inter-sensor
  consistency » (2 édits). PDF recompilé (9.9 MB).
- **Ch3 standalone** App. A : « loggers cover the summer only » → « summer window retained for
  the logger validation of this chapter » ; « not testable beyond summer with these data » →
  « not tested beyond summer in this chapter ». PDF recompilé.
- **Discussion D.7.1** : quatre phrases portées depuis le fichier autonome (plancher hivernal =
  bois, séparation feuille/bois par saison ; marcescence ; complémentarité FORMS). PDF recompilé.
- **Quarantaine** (déplacés, pas supprimés) dans `review/_archive_stale/` avec README :
  `manuscript_chap1_EN.md/.pdf` (copie périmée) et `disc_generale_perspectives_GEDI_EN.md/.pdf`.
- Non fait (jugement à toi) : réconciliation +0.98 (Ch1) vs +0.73 °C (Ch3) ; test β_below ≈ β_top
  dans Ch3 §3.3 ; provenance Table13 vs Table21b ; entrée Zotero Grulois.
