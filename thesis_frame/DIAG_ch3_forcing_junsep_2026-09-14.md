# Ch3 — quel forçage porte les chiffres de juin-septembre ? (diagnostic, 14/09/2026)

Produit par `c3_chs41_junsep_diag.R` (pure ré-extraction de NetCDF existants, aucune simulation).
Table par placette : `out_files/Chapter3_CHS41/tables/c3_chs41_junsep_diag_perplot.csv`.

## Le constat

`c3_junsep_full_recompute.R` (07/09) lit `out_files/Chapter3/nc_genuine53_windcorr`, l'arbre du
20 juillet, **avec correction de vent**. L'arbre `out_files/Chapter3_CHS41/nc` (24-28 août), sur
`MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc` et **sans** correction de vent, existait déjà et
n'a pas été utilisé. Le Ch1 est sur CHS41-Rmerge sans vent. Les deux chapitres ne sont donc pas sur
la même convention.

`gen_opt_CHS41_nowind.R` le dit explicitement : « NO wind correction: Chapter 3 runs on the CHS41
forcing as-is. The per-plot log-profile wind scaling is kept out of the chapter and illustrated only
in Chapter 1. » C'est l'intention ; l'extraction de septembre ne l'a pas suivie.

## Ce que change le forçage (macro = `musica_in_Blois_pblh.nc` dans les deux cas, juin-septembre)

| Scénario | | pente *r* | ΔTmax *r* | biais ΔTmax | dispersion |
|---|---|---|---|---|---|
| LiDAR fixe | windcorr (publié) | 0,79 | 0,49 | +0,73 | 0,24 |
| LiDAR fixe | **CHS41 sans vent** | **0,79** | **0,63** | **+0,30** | 0,23 |
| S2 seul | windcorr (publié) | 0,88 | 0,74 | +0,82 | 0,16 |
| S2 seul | CHS41 sans vent (n = 45) | 0,64 | 0,28 | +1,01 | 0,08 |
| S2 « opt » | windcorr (publié) | 0,85 | 0,05 | +1,17 | 0,08 |
| S2 « opt » | CHS41 sans vent | 0,78 | 0,24 | +0,72 | 0,05 |
| Combinaison | windcorr (publié) | 0,78 | 0,49 | +0,75 | 0,22 |
| Combinaison | CHS41 sans vent (n = 45) | 0,83 | 0,62 | +1,00 | 0,27 |

La pente du LiDAR fixe ne bouge pas (0,79). Le ΔTmax du LiDAR s'améliore nettement (0,49 → 0,63) et
son biais tombe de moitié. Les deux scénarios Sentinel-2 bougent beaucoup et **dans l'autre sens**.

**L'arbre CHS41 est incomplet** : `DYN_S2_RESCALED` (= Combinaison) n'a que 47 fichiers sur 53, et
seules 45 placettes sont exploitables. Une comparaison à quatre scénarios sur CHS41 n'est donc pas
encore possible sans relancer ces placettes.

## Les trois analyses orphelines, régénérées en juin-septembre

**1. Compétence de rang par strate de densité** (forçage LiDAR, coupure à la médiane du LAI_ALS, 3,80) :

| | dense (n = 26) | clair (n = 27) |
|---|---|---|
| windcorr, ΔTmax *r* | +0,77 | **+0,01** |
| CHS41 sans vent, ΔTmax *r* | +0,78 | **+0,82** |

C'est le point qui s'inverse. Sous la convention publiée du Ch3, la compétence se concentre bien dans
le dense, comme le disait l'ancienne phrase. Sous la convention vivante du Ch1, **le contraste
disparaît**. Par archétype : windcorr P1 −0,10 / P2 +0,61 / P3 +0,61 / P4 +0,71 ; CHS41 P1 +0,39 /
P2 +0,69 / P3 +0,51 / P4 +0,72.

**2. Triplet de validation** (forçage LiDAR, juin-septembre) : windcorr pente *r* = 0,79, ΔTmax
*r* = 0,49, biais +0,73 °C, dispersion 0,24. CHS41 sans vent : 0,79, 0,63, +0,30 °C, 0,23.
L'ancienne « amplitude d'environ 15 % » correspond à la dispersion, qui vaut ici 0,23-0,24.

**3. Incrément du profil dans l'ajustement model-free** (observations seules, donc indépendant du
forçage ; ΔTmax observé de juin-septembre, n = 53, LAI de forçage `LAI_ALS`) :

R² = 0,783 (LAI) → 0,865 (+fCover) → 0,872 (+Hmax) → 0,885 (+profil VCI).
**Incrément du profil ΔR² = +0,013, test de permutation p = 0,034.**
Sur la pente : 0,846 → 0,899 → 0,901 → 0,903, **ΔR² = +0,002, p = 0,257**.

Avec le LAI du clustering au lieu du LAI de forçage : +0,002 sur ΔTmax, +0,000 sur la pente.

Nuance par rapport à l'ancienne formule « indistinguishable from zero (+0,026, IC −0,007 à +0,070) » :
sur la **pente** la conclusion tient, sur **ΔTmax** l'incrément est petit (1,3 % de variance) mais
détectable. À ne pas écrire dans un manuscrit sans décision : c'est un résultat nouveau.

**4. Le rescaling de 23 %** n'est pas régénérable par ré-extraction. Le passage k = 0,5 → 0,65 est
précisément la réduction de 23 %, et le Ch3 tourne déjà en k = 0,65 ; il faudrait 53 simulations
MuSICA à k = 0,5 sous CHS41 pour retrouver un chiffre. La Discussion n'en a plus besoin : elle porte
désormais la formulation vivante du Ch1 §4.3, qui est une revendication de résolution, sans chiffre.

## Biais nocturne : le chiffre de D.6 est bon

`+1,5 à 1,7 °C` n'était écrit nulle part, mais il se vérifie. ΔTmin observé sur les 53 placettes en
juin-septembre = **−0,66 °C**, exactement la valeur que l'ancien Ch3 IMRADC publiait. Biais nocturne
simulé − observé, arbre publié (windcorr) : LiDAR fixe **+1,72**, S2 seul +1,64, Combinaison +1,65,
S2 « opt » +1,24 °C. Donc « +1,5 à 1,7 °C » décrit bien trois scénarios sur quatre, et c'est le
quatrième (« opt », +1,24) qui sort de la fourchette.

Sous CHS41 sans vent, le biais nocturne tombe à **+0,14 à +0,60 °C**. La phrase de D.6 dépend donc,
elle aussi, de la convention retenue.

---

# Ré-extraction complète sur CHS41 (14/09, après décision de Nathan)

## Ce qui a été produit

- **16 simulations manquantes complétées** (`gen_fix_CHS41_dyn.R`) : les 8 placettes « prob8 »
  (41_17, 41_18, 41_19, 41_27, 41_30, 41_39, 41_47, 41_49) pour `DYN_S2_ATBD` et pour
  `DYN_S2_RESCALED`. Elles échouaient parce que l'arbre portait des fichiers de **0 octet** que
  MuSICA ne réécrit pas. Les nouvelles sorties sont écrites dans
  `out_files/Chapter3_CHS41/**nc_fix**/`, rien n'a été supprimé ; l'extraction lit `nc` puis
  `nc_fix` en repli. Les 8 fichiers vides sont toujours là, à toi de dire s'il faut les écraser.
- **Ré-extraction complète** (`c3_junsep_full_recompute_CHS41.R`) : 53 placettes × 4 scénarios,
  aucun NA, deux références macro. Tables dans `out_files/Chapter3_CHS41/tables/` :
  `junsep_chs41_{perplot,by_archetype,paired_contrasts}_macro-{chs41,pblh}.csv`.
- **Figures 1 à 4 régénérées** (`c3_reframe_figs_junsep_CHS41.R`) sous
  `Fig{1,2,3,4}_reframe_junsep_chs41.png`. La Figure 5 (densité de points) ne dépend pas du forçage
  et reste `Fig5_reframe_junsep.png`.

## Choix de la référence macro

ΔTmax = maximum sous-couvert moins maximum macro. Le macro retenu est **la station CHS 41**, celle
qui force les simulations et celle du Ch1, pour les observations comme pour les simulations. Le
tableau `macro-pblh` est conservé comme test de sensibilité : les corrélations ne bougent
pratiquement pas (pente 0,79 et ΔTmax 0,63 dans les deux cas), seuls les offsets se décalent.

## Publié (windcorr) → CHS41 sans vent

| | publié | CHS41 |
|---|---|---|
| Prémisse, pente vs LAI : obs / sim | −0,92 / −0,91 | **−0,92 / −0,91** |
| Placette à placette, pente / ΔTmax | 0,79 / 0,49 | **0,79 / 0,63** |
| P2, pente : S2 vs LiDAR | 0,71 / 0,63 | **0,60 / 0,67** |
| P2, ΔTmax : S2 vs LiDAR | 0,81 / 0,61 | **0,75 / 0,69** |
| P4, pente : LiDAR / Comb / S2 / opt | 0,71 / 0,70 / 0,34 / 0,12 | **0,72 / 0,68 / 0,32 / 0,09** |
| Poolé, pente : S2 vs LiDAR | 0,88 / 0,79 | **0,86 / 0,79** |
| Contraste S2−LiDAR, pente, P2 | −0,024 | **0,000** [−0,009 ; 0,011] |
| Contraste S2−LiDAR, pente, P4 | +0,048 | **+0,063** [0,046 ; 0,082] |
| Combinaison au plus près du LiDAR | 0,013 | **0,039** (P4) |
| Biais ΔTmax, LiDAR fixe | +0,73 °C | **+0,30 °C** |
| Biais nocturne, 4 scénarios | +1,24 à +1,72 °C | **+0,14 à +0,60 °C** |

**Ce qui tient :** la prémisse, la fidélité de rang sur la pente, la domination du LiDAR en P4,
la faiblesse de P3, le paradoxe de Simpson en poolé, l'échec de « opt » en dense.

**Ce qui change de sens :** en P2 le Sentinel-2 brut **ne bat plus** le LiDAR sur la pente (0,60
contre 0,67) ; il le bat toujours sur ΔTmax (0,75 contre 0,69). Et la substitution déplace le
couplage **le plus en canopée dense**, pas en canopée intermédiaire : le contraste P2 est
exactement nul. Le titre de la §3.2 a donc été rendu neutre (« How far the forcing choice moves the
simulated coupling ») et le Key message reformulé métrique par métrique, en attendant ton arbitrage.

## Les trois analyses orphelines, sur la convention vivante

| | publié (windcorr) | CHS41 sans vent |
|---|---|---|
| Rang, moitié dense (n = 26) | +0,77 | +0,78 |
| Rang, moitié claire (n = 27) | **+0,01** | **+0,82** |

Le contraste dense/clair disparaît, ce qui confirme l'avertissement du Ch1 : c'était un effet de la
correction de vent. La clause « its validation skill concentrates there » reste donc retirée de D.3,
remplacée par le fait mesuré.

Ajustement model-free (observations seules, ΔTmax juin-septembre, n = 53, LAI de forçage) :
R² = 0,783 (LAI) → 0,865 (+fCover) → 0,872 (+Hmax) → 0,885 (+profil). **ΔR² du profil = +0,013,
permutation p = 0,034** ; sur la pente **+0,003, p = 0,263**. Réinséré en D.4.

Le triplet de validation devient pente *r* = 0,79, ΔTmax *r* = 0,63, biais +0,30 °C, dispersion
0,23. Le rescaling de 23 % reste non régénérable sans 53 simulations à k = 0,5 ; la Discussion
n'en a plus besoin.

---

## Addendum (14/09, soir) — le coefficient d'extinction, et le test des 23 %

**Découverte.** Le Ch3 ne tourne pas à *k* = 0,65 comme ses Méthodes l'affirmaient. Le raster branché,
`in_files/lai_z1_res_10_m.tif`, est byte-identique à `03_RESULTS/Blois/Metrics/Raw/lai_z1_res_10_m.tif`,
produit par `myPAI(Z, zmin = 1)` dont le défaut est **`k = 0.5`** (`02_CODES/libraries/functions_lidar2.R:8`),
jamais surchargé. Contrôle indépendant : le ratio médian `df_plots_lai.rds$LAI_ALS` / raster aux
53 placettes vaut 1,0009. **Ch1 et Ch3 sont donc sur le même *k*.** Le 0,65 est réel mais appartient au
Ch2, où il sert au *d*~opt~ ; il avait été transplanté dans le §2.2 du Ch3.

**Le test.** `gen_k065_CHS41.R` relance STATIC_ALS sur les 53 placettes avec la surface foliaire
× 0,769 (= 0,5/0,65), forme du profil, *H*max et fCover tenus fixes. 53/53 en 1,8 min, écriture dans
`Chapter3_CHS41/nc_k065/`, rien de supprimé. Extraction par `c3_k065_extract_CHS41.R`, macro CHS 41,
juin-septembre, z = 1 m. Sorties : `tables/k065_{perplot,paired}_macro-chs41.csv`.

**Résultat — l'affirmation du Ch1 est confirmée.**

| Quantité | *k* = 0,5 | *k* = 0,65 | Décalage [IC 95 %] |
|---|---|---|---|
| Pente simulée | — | — | **+0,029** [+0,023 ; +0,036] |
| ΔTmax simulé | — | — | **+0,238 °C** [+0,178 ; +0,305] |
| Décalage pente par archétype | P1 +0,000 | P2 +0,031 | P3 +0,011 · **P4 +0,052** |
| Décalage ΔTmax par archétype | P1 +0,001 °C | P2 +0,194 °C | P3 +0,058 °C · **P4 +0,477 °C** |
| *r*(sim, obs) pente | 0,792 | 0,802 | classement préservé |
| *r*(sim, obs) ΔTmax | 0,625 | 0,559 | recul de 0,066 |
| Biais chaud ΔTmax | +0,303 °C | +0,541 °C | l'échelle suit la feuille |

Moins de feuillage, tampon plus faible : la pente remonte vers 1 et le sous-bois se réchauffe. L'effet
est nul en trouée, où il n'y a rien à rééchelonner, et maximal en dense.

**Nuance à trancher.** Le **classement** survit au décalage sur la pente (0,79 → 0,80). La prescription
du Ch1, « an optical proxy should keep its leaf-area error below 23 % », vaut donc pour le couplage
absolu et non pour le rang. Le Ch1 ne fait pas la distinction.

**Écrit.** Ch3 §2.2 (le *k* corrigé), §2.6 (troisième test de robustesse), Annexe B *Extinction
coefficient* (l'affirmation remplacée par la mesure). Ch1 l. 599 : le renvoi faux « (Appendix D) »
retiré. L'intro et la discussion générales ne portaient plus les 23 %, rien à y faire.

---

## Les trois arbitrages, tranchés (14/09, soir)

Nathan a pris les trois. Le chiffre qui les a tranchés est un **bootstrap apparié sur la différence
des corrélations de classement** entre deux scénarios, à l'intérieur d'un archétype (B = 4000,
`c3_rank_diff_bootstrap_CHS41.R`, sortie `tables/k_rankdiff_S2_vs_LiDAR_macro-chs41.csv`). Il
n'existait pas avant : le chapitre comparait des *r* sans jamais porter la taille d'échantillon
dans la comparaison.

| Archétype | *n* | pente : S2 − LiDAR [IC 95 %] | ΔTmax : S2 − LiDAR [IC 95 %] | résolu ? |
|---|---|---|---|---|
| P1 | 8 | −0,002 [−0,44 ; +0,42] | +0,216 [−0,12 ; +0,91] | non |
| **P2** | 12 | **−0,070 [−0,47 ; +0,34]** | **+0,059 [−0,40 ; +0,36]** | **non** |
| P3 | 13 | −0,161 [−0,36 ; −0,02] | −0,957 [−1,27 ; −0,57] | **oui** |
| P4 | 20 | −0,402 [−0,66 ; −0,18] | −0,467 [−0,96 ; −0,10] | **oui** |

**Ce que ça change.** La scission par métrique en P2 était du bruit dans les deux sens : ni l'avance
sur ΔTmax ni le retard sur la pente ne survivent à 12 placettes. Le chapitre affirmait une avance
qu'il ne pouvait pas mesurer. À l'inverse, **P3 est séparable** et dans le même sens que P4, ce que le
chapitre passait sous silence en le réduisant à « la strate faible ». Le résultat défendable est donc :
le LiDAR bat le S2 partout où l'échantillon peut le dire (P3 et P4), et les deux sont indiscernables
là où il ne le peut pas (P1 et P2).

**Écrit** — Ch3 : Key message, résumé (Results et Conclusion), §2.5 (le bootstrap en Méthodes),
titre et corps du §3.2 (« The forcing choice moves the coupling most in dense canopy », P4 en tête),
§3.3 (P2 non séparable, P3 ajouté), §4.1 (dense en tête, P2 rétrogradé), §4.2 (le « 0,01 » périmé
devient 0,04, et la règle de rejet mentionne aussi ce que l'échantillon ne résout pas), §5.
Ch1 : la prescription des 23 % gagne sa clause (« wherever the absolute coupling matters, although
the between-plot ranking survives a shift of that size »).

**Contrôles :** Ch3 24 p., 0 tiret cadratin, 0 espace sécable avant unité (49 insécables), phrase la
plus longue 35 mots, 5 figures résolues. Ch1 38 p., 0 citation non résolue.

**Deux pièges rencontrés, pour mémoire.** (1) `P == get("P")` dans un `i` de data.table se référence
la colonne : le bootstrap renvoyait n = 53 pour les quatre archétypes. (2) Le re-rewrap des
paragraphes détruit les espaces insécables ; il faut les réappliquer mécaniquement après chaque
passe, et `chr(160)` plutôt qu'un littéral, qui se fait normaliser en chemin.

---

## Passe de cohérence finale (14/09, nuit) — ce que les passes précédentes avaient laissé

Nathan a demandé si l'Intro, le Ch3 et la Discussion étaient bons, plans compris. Non, ils ne
l'étaient pas. Cinq résidus, tous trouvés en recalculant depuis les tables plutôt qu'en relisant.

**1. Trois légendes de figure sur cinq étaient périmées.** La passe CHS41 avait régénéré les PNG
sans relire les légendes du markdown. Fig. 1 donnait ΔTmax simulé *r* = −0,71 là où les deux macros
donnent **−0,84** ; Fig. 2 disait que le S2 brut « lowers the coupling slope in P2 », alors que le
contraste y est **0,000** ; Fig. 3 disait « Raw Sentinel-2 leads LiDAR in P2 », le verdict qu'on
venait de retirer du corps. Corrigées, les cinq légendes tiennent sous 60 mots.

**2. La Discussion surestimait l'effet des 23 %, deux fois.** Elle écrivait « moves the simulated
coupling by more than Chapter 3 resolves against the loggers ». La mesure dit l'inverse : le décalage
vaut +0,029 sur la pente et +0,238 °C sur ΔTmax, soit **moins** que la dispersion résiduelle
intra-archétype (0,033 à 0,044 sur la pente, 0,35 à 0,49 °C sur ΔTmax), et le classement survit.
Réécrit sur ce que la mesure dit : la grandeur bouge, l'ordre non.

**3. Le ratio de dispersion « at or below 0.23 » était faux.** Il vaut 0,23 sur ΔTmax poolé mais
**0,38 sur la pente**, et surtout la compression n'est pas uniforme : elle est concentrée dans les
trouées (0,015 en P1) et **disparaît en dense** (1,03 sur la pente, 0,79 sur ΔTmax en P4). Réécrit
avec les deux valeurs poolées et la localisation.

**4. Le « 23 % cut in Chapter 1 »** appartient désormais au Ch3, qui porte le test. Réattribué.

**5. Le bloc « attend ton arbitrage »** traînait, identique, dans les trois plans FR, et le plan de
l'intro n'avait pas reçu la section du verdict P2. Les trois sont alignés, docx régénérés.

**État final vérifié :** Ch3 24 p., Intro 45 p., Discussion 10 p. Partout 0 tiret cadratin,
0 espace sécable avant unité, 0 citation non résolue. Phrase la plus longue nette de citations :
35 mots au Ch3, 34 en Discussion, 36 en Intro (la question RQ3, exception acceptée de longue date).
Chaque chiffre du §3 du Ch3 a été recalculé depuis `junsep_chs41_*_macro-chs41.csv` et concorde.
