# Plan détaillé — Article Chapitre 1
### (calé sur la « Histoire de la thèse » officielle + résultats réels du pipeline)

> ⚠️ **OBSOLÈTE sur la méthode d'attribution (mis à jour 2026-06-20).** Le Shapley a été **entièrement
> retiré** du Chapitre 1 ; l'attribution est désormais une **analyse de sensibilité par perturbation de
> traits** (sur la variété des données, sans baseline ni coalitions). Le claim est adouci : **le LAI domine
> partout, sa domination s'érode P1→P4, le profil vertical devient COMPARABLE (co-mène, ne dépasse pas) le LAI
> en P4 dense**. Voir mémoire `density_dependent_shapley_main.md` et `review/manuscript_chap1_EN.md`. Le reste
> du plan (angle éco, glass ceiling r≈0.93, annexe H3/S2, périmètre) reste valable.
>
> ⚠️ **MISE À JOUR 2026-06-27 — moteur & référence.** L'article passe sur **MuSICA v3.2.3 en mode
> yoyo/iter** (couplage ABL) avec un **nouveau forçage basé sur la station locale CHS 41** (+ gap-fill
> ERA5 ou SAFRAN, choix à faire), à venir sous quelques jours. La **référence ΔTmax devient la station
> (1,5 m)** (bloc A comité réu 26/06), le biais forçage↔station passe en annexe. ⚠️ **Tous les chiffres
> ci-dessous (r=0.93, +0.9 °C, +0.2 °C ERA5, R²=0.62, splits, sensibilités) sont des runs v3.2.0/ERA5 et
> doivent être REFAITS au re-run** ; le r=0.93 en particulier est à re-tester (v3.2.3 donnait 0.63 sous
> forçage ERA5). Voir mémoire `article_v323_station_forcing.md` et la bannière en tête de
> `manuscript_chap1_EN.md`.

**Angle** : écologique (contribution relative LAD vertical vs LAI total au buffering),
rigueur méthodo valorisée (attribution sous colinéarité, off-manifold). **Périmètre** :
mono-site (Blois, chênaie), été 2021, 53 HOBO — énoncé tôt. **H3/S2** : annexe + perspective.

ΔTmax = écart thermique diurne micro − macro (macro = météo régionale milieu ouvert).

---

## ⚠️ DEUX TENSIONS À TRANCHER (narration officielle ↔ pipeline réel)

**T1 — Méthode : GAMM (narration) vs Shapley exact (pipeline).**
La narration officielle dit *« GAMMs + évaluation de la concurvité pour isoler l'effet
partiel du vertical »*. Mais le pipeline réel **remplace** le GAMM par une **attribution
Shapley exacte** sur 16 coalitions MuSICA — ce qui EST la résolution du verrou
concurvity/colinéarité (Shapley moyenne sur toutes les coalitions, robuste là où le GAMM
souffrait de concurvity et où les partial plots évaluent hors-variété).
→ **Recommandation** : l'article utilise **Shapley** ; le GAMM = approche V1 dépassée,
mentionnée en intro comme le verrou levé. À acter dans le texte (sinon incohérence).

**T2 — Résultat ↔ hypothèses : H1/H2 misent sur le vertical (LAD) ; le résultat dit LAI.**
H1/H2 posent que **le profil vertical LAD est le paramètre clé**. Or l'attribution Shapley
donne **LAI > fCover > LAD ≈ 0 > Hmax** (φ_LAD(ΔTmax) ≈ 0, change même de signe selon
branche). **Le LAD apporte peu au-delà du LAI total** dans la chênaie de Blois.
→ **Recommandation** : reframe honnête. L'histoire devient *« on s'attendait à ce que la
distribution verticale soit déterminante ; on montre qu'à l'échelle du buffering ΔTmax,
c'est le LAI total qui domine et que la forme verticale n'ajoute presque rien »*. C'est un
résultat **publiable et fort** (nuance le dogme 3-D), pas un échec. H2 (top-heavy →
refroidit plus) peut être testé spécifiquement en sous-analyse même si φ_LAD global ≈ 0.

---

## Titre (travail)
> *Total leaf area, not its vertical arrangement, governs simulated summer understory
> buffering: an exact-Shapley attribution of LiDAR canopy traits, and its glass ceiling
> against field sensors.*

## Questions de recherche (officielles)
1. Sensibilité des T° sous couvert simulées (MuSICA) aux variations de distribution
   **verticale ET horizontale** de la végétation.
2. **Contribution relative du LAD (profil vertical) vs LAI (surface totale)** sur ΔTmax.

## Hypothèses (officielles — à confronter aux résultats)
- **H1** — La structure verticale module les transferts radiatifs/turbulents → le profil
  LAD est un paramètre clé pour simuler fidèlement le microclimat. *(partiellement infirmée
  : voir T2)*
- **H2** — À LAI et hauteur constants, les profils **top-heavy** (biomasse concentrée en
  haut) induisent une plus forte atténuation diurne (rayonnement intercepté plus haut,
  sources de chaleur déplacées vers le sommet). *(à tester spécifiquement)*
- **H3 (perspective, annexe)** — Le signal optique 2-D (S2) donne une illusion de buffering
  (proxy de structure) → motive le Ch.2.

---

## 1. Introduction
1. Refuges thermiques / découplage microclimatique (De Frenne, Lembrechts, Zellweger).
2. Structure de canopée comme moteur — mais **quelle dimension** : LAI total vs **profil
   vertical 3-D** ? Positionnement vs **Bouwen** (LiDAR+MuSICA voisin, seuil de densité).
3. Verrou V1 : GAMM monolithique, **colinéarité / concurvity**, effet « Date » boîte noire.
   → besoin d'une attribution qui explore **toutes les combinaisons**, robuste à la
   colinéarité, sans partial plots off-manifold.
4. Approche : **scénarios MuSICA** + **Shapley exact** (2⁴), validé contre 53 HOBO.
5. Questions / H1–H2 (+ H3 perspective).

## 2. Matériel & Méthodes
> Garde-fous reviewer : périmètre mono-site énoncé tôt ; **×2 LAI = aire foliaire deux
> faces** (pas « clumping ») ; gradient vertical = illustration mécaniste.

- **2.1 Site & capteurs** — Blois chênaie, été 2021, 53 HOBO ~1 m. *Paragraphe de périmètre.*
- **2.2 LiDAR & traits** — profils LAD 0.5 m (capte le sous-étage) ; traits LAI, Hmax,
  fCover, **VCI**, forme du profil ; LAI une-face → **×2 = aire deux faces** (MuSICA).
- **2.3 Typologie** — **K-Means** sur métriques ALS (LAI, Hmax, VCI, fCover) →
  **archétypes P1–P4** ; (+ FPCA de la forme LAD en appui) ; clustering gelé, cLHS.
- **2.4 MuSICA & plan de simulation** — choisir/présenter :
  - *Option A (in-silico contraint)* : plan factoriel LAI×hauteur × gabarits LAD
    d'archétypes, **filtré par l'enveloppe structurelle** du massif (rasters complets).
  - *Option B (in-situ)* : profils LAD réels par pixel tirés en strates (covariance
    naturelle préservée).
  - *(Pipeline réel ≈ B + coalitions ; T1)* : 16 coalitions (2⁴), baseline fCover moyen
    0.869, LAI_b cLHS ≈ 6.73, binaire v3.2.0.
- **2.5 Extraction — 2 conventions** : (a) ΔTmax/ΔVPDmax absolus → 1 m interp + −2 h ;
  (b) pente micro–macro (split buf/amp) → nair==1, sans shift.
- **2.6 Attribution** — **Shapley exact** φ par trait/placette (ΔTmax, ΔVPDmax), additivité
  vérifiée, CI bootstrap ; **isole φ_LAD vs φ_LAI** (= la question 2, en remplacement du
  GAMM). **Shapley conditionnel** (φ_LAD | Hmax) borne l'off-manifold.
- **2.7 Validation** — REF vs 53 HOBO ; ΔTmax journalier ; pente micro–macro (45/8).

## 3. Résultats (`affirmation → figure → chiffre`)

### R1 — Typologie en 4 archétypes
→ `fig_clusters_lad_profiles`, `fig_archetype_profiles`, `fig_fpca_*`, `tab_cluster_structure`
→ **4 archétypes P1 (clair) → P4 (dense)**.

### R2 — LAI domine ; la distribution verticale ajoute peu *(cœur, reframe T2)*
→ `fig_shapley_ranking`, `fig_shapley_pooled`, `tab_shapley_ranking.csv`
→ **LAI > fCover > LAD ≈ 0 > Hmax** ; φ_LAI(ΔTmax) ≈ −0.17 °C, fCover −0.08,
  **LAD ≈ 0**, Hmax −0.02 ; Σφ ≈ −0.311 °C ; même ordre sur VPD.
→ **Message : à l'échelle du ΔTmax, le LAD vertical n'apporte quasi rien au-delà du LAI.**

### R2bis — Test CONTRÔLÉ in-silico de H2 → **(va en DISCUSSION, pas en Résultats)**
> Article = **focus terrain** : les Résultats restent empiriques (R2/R2ter). R2bis est le
> support *mécaniste théorique* qui explique, en Discussion, le résultat terrain (D1bis).
→ `annex/fig_h2_controlled_topheavy.png`, `annex/tab_h2_controlled_topheavy.csv`
→ À **LAI & Hmax fixés**, profils Beta bottom→top-heavy : ΔTmax **décroît** quand la
  biomasse monte (top-heavy refroidit plus) → **H2 vérifiée mécaniquement**.
  - LAI=6 : ΔTmax +0.75 → +0.16 °C (amplitude **0.59 °C**, pente −1.18 °C/COM).
  - LAI=12 : ΔTmax −0.06 → −0.33 °C (amplitude **0.27 °C** → l'effet **diminue de moitié
    quand le couvert est dense** : saturation de l'interception).
→ **Réconcilie T2** : la forme verticale agit dans le bon sens (H2 vraie) mais son ampleur
  (~0.3–0.6 °C) est petite et **sature avec le LAI** → elle ne ressort PAS comme facteur
  attribuable (φ_LAD ≈ 0) une fois le LAI total pris en compte dans les couverts réels
  corrélés. Message : *« le profil compte mécaniquement, mais c'est un effet de second
  ordre, dominé par le LAI total et qui s'efface quand la canopée est dense ».*

### R2ter — **Validation empirique de H2 sur le terrain : effet vertical NUL une fois les traits contrôlés**
→ `annex/fig_h2_field_validation.png`, `annex/tab_h2_field_validation.csv`
  (script `scripts/make_h2_field_validation.R` — **données réelles uniquement** : profils LAD
   ALS z0.5 → centre de masse COM, buffering ΔTmax obs HOBO, LAI/Hmax/fCover, clusters P1–P4)
→ **Brut** : COM corrèle avec le buffering (r=−0.54) — mais c'est le gradient de densité
  (P1 sparse/bas-COM amplifie → P4 dense/haut-COM tamponne).
→ **À LAI seul contrôlé** : r partiel = **−0.43 (p=0.001)** — *semble* valider H2…
→ **À LAI+Hmax+fCover contrôlés** (= conditions du Shapley) : r partiel = **−0.04 (p=0.79)**,
  **NUL**. Le COM ne faisait que **proxy fCover (r=0.62) et Hmax (r=0.47)**.
→ **Conclusion** : en conditions réelles, la position verticale de la biomasse **n'a aucun
  effet indépendant** sur le buffering une fois les autres traits pris en compte →
  **confirme empiriquement φ_LAD ≈ 0**. Unifie Shapley + R2bis + R2ter : le mécanisme H2
  existe (R2bis, 2ᵈ ordre, sature) mais est **indiscernable de zéro sur le terrain** dès
  qu'on contrôle les traits dominants corrélés.
→ **Bonus méthodo (→ D3)** : le contraste r=−0.43 (naïf) vs r=−0.04 (contrôlé) est une
  démonstration directe de **pourquoi il faut le Shapley** (contrôle TOUS les traits) plutôt
  qu'une corrélation partielle / un effet partiel GAMM naïf qui confond les traits corrélés.

### R3 — La 3-D distribue la variance, mais plafond de verre
→ `fig_validation_perplot`, `fig_validation_by_cluster`, `fig_validation_forward_shapley`
→ **r spatial ≈ 0.93** (classement bien capté) ; **biais +0.9 °C** ; split terrain **45/8**,
  modèle **41/12** ; **~8 placettes amplificatrices non reproduites** (pente HOBO →1.44).
> Phrase d'attaque : *« on attribue le buffering simulé ; modèle validé sur le classement
> spatial (r≈0.93), plafond connu sur le niveau absolu et les amplificateurs. »*

### R4 — Mécanisme : gradient vertical (illustration)
→ `annex/fig_vertical_profiles_wind_rh_vpd_V9vsV11` (⚠️ illustration, pas de réf verticale).

## 4. Discussion
- **D1** — Pourquoi le LAI domine et pas le profil : interception/extinction quasi saturée
  en chênaie dense ; cohérent avec le seuil de Bouwen ; **Canopy Ratio (Starck) n'ajoute
  rien** au-delà du LAI (`optional/10`). Nuance le dogme « 3-D indispensable ».
- **D1bis — La différence théorie ↔ terrain (clôt le débat « le vertical compte-t-il ? »).**
  *Résultat terrain* (R2ter) : à LAI+Hmax+fCover contrôlés, la position verticale (COM) n'a
  **aucun effet** sur le buffering observé (r=−0.04, p=0.79) — cohérent avec φ_LAD≈0.
  *Explication mécaniste* (R2bis, in-silico) : le mécanisme de H2 **existe** (top-heavy
  refroidit plus) mais il est **petit et sature** avec le LAI. Quantification signal/bruit :
  - sensibilité modèle ≈ 0.5–1.2 °C par unité de COM ;
  - mais la variation de COM **réellement disponible** à structure fixée est étroite
    (SD résiduel ≈ 0.09, vs span 0.50 imposé artificiellement dans R2bis) ;
  - → effet ΔTmax attendu sur le terrain ≈ **0.05–0.10 °C**, soit **~10× sous le bruit
    inter-placettes (SD ≈ 0.65 °C)** → **non identifiable**.
  **Donc pas de désaccord modèle↔réalité** : le profil module bien le ΔTmax *dans le bon
  sens*, mais son amplitude sur la variabilité architecturale réelle est d'un ordre de
  grandeur sous le signal inter-placettes → **ce n'est pas un levier opérant en peuplement
  réel**. (Caveat : COM ALS bruité en basses strates → dilution additionnelle.)
- **D2** — Plafond de verre : classement OK, niveau absolu non ; +2 °C → **hypothèse** ET
  manquante (+ alternatives : hauteur forçage, humidité sol, radiatif capteur).
- **D3** — Apport méthodo (valorisé) : Shapley exact vs GAMM (concurvity, boîte noire
  Date→forçages ERA5) ; toutes les coalitions ; robuste colinéarité ; conditionnel borne
  l'off-manifold. Limite : attribution **interne au modèle, validée spatialement**.
- **D4** — Robustesse (classement identique v3.2.0/v3.2.3, annexe ; v3.2.0 r 0.93 > 0.63) +
  périmètre mono-site.
- **D5** — Perspective Ch.2 : illusion optique S2 (annexe). **Baseline naïve S2+FORMS-H →
  MuSICA = scénario du Ch.3** (renvoi).

## 5. Conclusion
Le **LAI total** est le premier levier du buffering estival simulé ; **la forme verticale du
profil ajoute peu** (nuance forte vs H1/H2). La 3-D distribue bien la variance entre
placettes (r≈0.93) mais le modèle plafonne (niveau absolu + amplificateurs) → cadre honnête
« plafond de verre », ouverture Ch.2.

---

## Figures principales
1. `fig_clusters_lad_profiles` (+ `fig_archetype_profiles`)
2. `fig_shapley_ranking` (phare)
3. `fig_shapley_by_cluster`
4. `fig_validation_perplot`
5. `fig_validation_forward_shapley`
6. `fig_vertical_profiles_wind_rh_vpd_V9vsV11` (illustration)

## Annexes
A1 FPCA · A2 Shapley conditionnel LAD (off-manifold) · A3 métriques extra (Tmin/amp/stabilité)
· A4 robustesse v3.2.0/v3.2.3 · A5 hauteur nair==1 vs 1 m · A6 LAD z0.5 vs z1 · A7 S2/illusion optique (→Ch.2)

## Checklist pré-emption Reviewer 2
- [ ] Périmètre mono-site énoncé tôt.
- [ ] Validation = **r≈0.93 (spatial)** en chiffre phare ; +0.9 °C + 8 amplificateurs = le plafond.
- [ ] Attribution = **interne au modèle, validée spatialement** (phrase d'attaque Résultats).
- [ ] **T1 tranché** : Shapley assumé comme méthode (GAMM = V1 levée), cohérent partout.
- [ ] **T2 assumé** : LAI domine, LAD ≈ 0 — hypothèses reframées, pas présentées comme confirmées.
- [ ] ×2 = aire deux faces (jamais « clumping »).
- [ ] Gradient vertical = illustration ; ET = hypothèse parmi d'autres ; confusion cluster↔LAI reconnue.
