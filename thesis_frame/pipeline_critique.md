# Revue critique du pipeline — Article Chapitre 1

Audit de bout en bout du code produisant les figures et chiffres de
`review/manuscript_chap1_EN.md` (pipeline `pipeline/` + scripts `c1_*`/`c3_*` +
librairie `R/`). Lecture du code réel (pas des commentaires), vérifications sur
les artefacts primaires (tables, NetCDF, échantillon cLHS).

Classement en trois tiers :
- **Tier 1 — science à risque** : menace une conclusion du papier.
- **Tier 2 — reproductibilité / provenance** : le résultat existe mais n'est pas
  traçable/régénérable proprement.
- **Tier 3 — hygiène / introduit cette session**.

Pour chaque point : verdict, preuve (fichier:ligne ou artefact), correctif.

---

## TIER 1 — Science à risque

### 1.1 [BLOCKER] La décomposition du résidu (Fig 5, R²=0.51) tourne sur n=44 et exclut 4 des 5 plus gros amplificateurs

**Preuve.** `scripts/make_residual_vs_topo.R` régresse le résidu (obs − modèle ΔTmax)
sur 7 covariables (elevation, slope, twi, northness, sd_height, rumple,
gap_fraction). Dans `outputs/figures_pipeline/annex/tab_residual_vs_topo.csv`,
**9 plots ont NA sur toutes les covariables topo/hétérogénéité** :
`41_17, 41_18, 41_19, 41_27, 41_30, 41_39, 41_47, 41_49, 41_55`. `lm()` les
supprime → la régression tourne sur **n=44, pas 53** (warnings "Removed 9 rows"
dans `run_full_chap1.log`).

Les plots supprimés ne sont pas aléatoires : `41_18` (résidu **+3.00**),
`41_39` (**+2.50**), `41_49` (**+2.15**), `41_19` (**+1.03**) sont **les plus
gros résidus positifs**, c.-à-d. exactement les plots amplificateurs que la
section 3.3/4.2 prétend expliquer par la topographie et les trouées. C'est un
défaut de **couverture de covariables** (ces 9 plots tombent hors des tuiles
LiDAR externes), pas une intention — mais la conséquence est la même : la
décomposition est **muette sur les amplificateurs qu'elle invoque**.

**Aggravant.** `gap_fraction`, présenté comme dominant (r=0.58), est **non nul
pour seulement 3 plots retenus** : `41_07` (0.05), `41_56` (0.06), `41_59` (0.02).
Le r marginal et le β partiel (le seul significatif, p=0.001) reposent sur ~3
points contre une distribution massée à zéro, dont `41_07` qui est à la fois le
plus gros amplificateur retenu (+1.99) et un point gap non nul → point à très
fort levier. "Résidu dominé par gap_fraction" = un quasi-indicateur binaire sur
3 plots, ajusté in-sample, sans rééchantillonnage.

**Aussi.** R²=0.51 est **in-sample, non ajusté** (R² ajusté = 0.41 est calculé
mais non rapporté), 7 prédicteurs pour n=44 (≈6 obs/prédicteur), aucune
validation croisée.

**Correctif (obligatoire, le résultat est porteur pour §3.3 et §4.2).**
1. Récupérer la couverture topo/CHM pour les 9 plots manquants et **re-tourner
   sur 53**. Si impossible, **dire explicitement n=44 et lister les 4
   amplificateurs exclus** dans le texte et la légende.
2. Rapporter le **R² ajusté (0.41)** comme chiffre de tête, avec n.
3. Pour `gap_fraction` : reconnaître que l'effet repose sur 3 plots, ou le
   traiter en indicateur binaire, ou retirer la prétention de "dominance".

### 1.2 [MAJEUR — réfuté par les données, à blinder dans le texte] Le glissement de mean|φ| LAI→LAD est-il un artefact d'échantillonnage ?

**L'objection Reviewer 2.** mean|φ| mesure l'effet de déplacer un trait de son
baseline-cluster à sa valeur réalisée → il croît avec (dispersion intra-cluster
du trait) × (sensibilité physique). Si l'ouvert (P1) a une large gamme de LAI et
le dense (P4) une gamme étroite, alors "LAI domine l'ouvert / s'effondre dans le
dense" serait en partie un effet de **gamme d'échantillonnage**, pas de physique.

**Vérification (dans `clhs_sample_floor05_v2.rds`, mapping densité P1→P4).**

| Archetype | LAI moy | LAI SD | fCover SD |
|---|---|---|---|
| P1 (ouvert, LAI≈2.35) | 2.35 | **1.39** | 0.07 |
| P2 | 6.05 | 1.52 | 0.07 |
| P3 | 7.72 | 1.25 | 0.06 |
| P4 (dense, LAI≈10.8) | 10.80 | **1.27** | 0.06 |

**La dispersion intra-cluster du LAI est quasi constante (SD 1.25–1.52).** Donc
l'effondrement de l'importance LAI en P4 (0.41 → 0.10 °C) se produit à
**dispersion d'échantillonnage comparable** → c'est de la **saturation
physique**, pas un artefact de gamme. Même conclusion pour fCover (SD ~0.06–0.07
partout, importance 0.21 → 0.01).

**Action.** Le résultat tient, mais l'objection est trop évidente pour être
laissée implicite. Ajouter une phrase + le tableau de dispersion en annexe, et
**lier explicitement** à l'expérience contrôlée H2 (`make_h2_controlled_topheavy.R`,
mémoire `h2_controlled_topheavy`) : à LAI et Hmax **fixés**, le profil top-heavy
refroidit plus → preuve sans dispersion que la saturation est physique. C'est la
réponse blindée à donner d'avance.

### 1.3 [À CONFIRMER, impact à scoper] Double-doublement du LAI dans la branche z05 (qui alimente la validation Fig 4)

**Preuve (tracée).** `pipeline/03z_musica_z05.R:41` `df[, LAI := LAI * LAI_X2]`
double le LAI scalaire (log : "LAI med now 8.49"). Puis `run_musica_one`
(`R/musica.R:69`) `lai <- 2*lai` double **encore** → ~17 atteint MuSICA
(non physique pour du chêne ; LAI deux-faces réaliste ≈ 8–14). De plus les
**baselines viennent de `df_clhs`** (`03z:48`, floor05_v2) qui ne passe **pas**
par la ligne 41 → le toggle réel-vs-baseline du LAI est sur des **échelles
incohérentes**. Le NetCDF de sortie n'expose pas le LAI d'entrée, donc je ne peux
pas lire 17 vs 8.5 directement.

**Scope (important).** La branche z05 alimente **uniquement la Fig 4 (validation)**,
pas l'attribution Shapley (qui est le run externe floor05_v2). Direction de
l'effet : un LAI gonflé → canopée plus dense → **plus** de buffering. Or le
modèle **sous-buffer quand même** (biais chaud +0.9 °C) → corriger le bug rend le
biais chaud **plus grand** et le plafond de verre **plus marqué**. **La conclusion
qualitative survit** ; seuls bougent les chiffres 44/9 et la pente 0.95.

**Correctif.** Logger le LAI réellement passé au binaire pour un plot z05 1111.
Si ~17 : supprimer `03z:41` (laisser le seul ×2 de `musica.R:69`) et
re-générer la validation z05. Re-rapporter 44/9, pente, biais.

---

## TIER 2 — Reproductibilité & provenance

### 2.1 [MAJEUR] Le résultat principal (Fig 2/3, attribution density-dependent) est HORS du pipeline propre, et sa provenance d'échantillon est ambiguë

- Les stages `pipeline/03→05` calculent un Shapley à **baseline GLOBAL** avec
  **reporting par cluster** (`05_shapley.R:39` merge le cluster *après* φ), PAS un
  baseline par cluster. L'analyse par-cluster du papier vit dans les scripts ad-hoc
  `c3_shapley_chunk.R` (+ `c3_shapley_merge.R`), hors `run_all.R`.
- **Ambiguïté d'échantillon.** `c3_shapley_chunk.R:24` utilise `floor05_v2`
  (canonique, LAI_b≈6.73). Mais `c3_shapley_clhs_2x.R:20` utilise `floor05`
  (l'échantillon que tes notes marquent FAUX, LAI~3.1). **Les deux scripts
  écrivent le MÊME fichier** `NC_Full/.../FigSh_shapley_clhs_clusters.png` → celui
  qui tourne en dernier gagne.
- La figure citée par le manuscrit, `out_files/Chapter3/figures/FigSh_shapley_clhs_cluster.png`
  (singulier, 52589 o, 13 juin), est une **copie orpheline** : **aucun script
  suivi n'écrit ce chemin exact**. Impossible de certifier qu'elle vient de
  floor05_v2 et non de floor05.

**Correctif.** Porter le calcul par-cluster (v2) dans une étape numérotée du
pipeline ; faire écrire à `c3_shapley_merge.R` un nom **versionné/non
collisionnant** ; supprimer ou neutraliser `c3_shapley_clhs_2x.R` (floor05) ;
régénérer Fig 2/3 depuis v2 et figer le chemin cité.

### 2.2 [MAJEUR] Moran's I (I=0.30, p<0.001) n'existe pas dans le dépôt

`grep` global (`moran|nb2listw|knearneigh|spdep|ape::|listw`) → **rien**. Aucune
définition de poids spatiaux, aucun n, aucun test de permutation. Le chiffre
n'est **pas reproductible depuis le dépôt** (vit peut-être dans un notebook hors
dépôt). Les coordonnées et le résidu sont dans `tab_residual_vs_topo.csv`, donc le
test est faisable.

**Correctif.** Committer le code Moran (poids, n, méthode), ou retirer la
prétention. Ne pas la présenter comme acquise tant que non reproductible.

### 2.3 [MAJEUR] Les figures principales mélangent TROIS provenances

- **Fig 4 (validation)** : `outputs/figures_pipeline_z05/` → branche **z05**
  (LAD 0.5 m). Source des chiffres r=0.93 / +0.9 °C / 44-9 / pente 0.95.
- **Fig 5 (résidu) & Fig 6 (profils verticaux)** : `outputs/figures_pipeline/annex/`
  → branche **z1 legacy** (`v10_fcovmean`). `fig_residual_vs_topo.png` n'existe
  PAS dans la branche z05.
- **Fig 2/2b/3/A1/B1 (attribution)** : `out_files/Chapter3/figures/FigSh_*` → run
  externe `c3_*` (échantillon ambigu, cf. 2.1).

Donc la validation (z05), la décomposition du résidu (z1) et l'attribution
(externe) ne sont **pas le même jeu de simulations**. Le "plafond de verre" (z05)
et l'attribution (floor05) ne proviennent pas de la même source LAD.

**Correctif.** Choisir une branche LAD de référence (z05 recommandée si c'est la
validation publiée) et **régénérer toutes les figures principales depuis cette
branche**, ou justifier explicitement le mélange en Méthodes.

### 2.4 [MOYEN] Dérive de documentation : pente "nair==1 / split 45-8"

`00_config.R:39` `SLOPE_USE_NAIR1 = FALSE` (pente calculée à 1 m fixe ;
"legacy 45/8 superseded"). Mais : `00_config.R:37` (commentaire) et `:65`
(bannière au chargement), `04_extract.R:5`, `07_validation.R:2-4` et
`pipeline/README.md:45,51` disent encore "nair==1, 45/8, r=0.92". L'ASSERT
`slope_buf=45L/slope_amp=8L` (`00_config.R:50`) **ne plante pas**, mais seulement
parce qu'il n'est exécuté que sur la branche par défaut v10 (`strict`), qui
**donne 45/8** — alors que **le manuscrit rapporte la branche z05 (44/9)**.
L'assert teste donc une **autre branche que celle publiée**.

Le manuscrit lui-même est interne-cohérent (Méthodes 2.5 = 1 m, no shift ;
Résultats = 44/9, pente 0.95 vs 0.86). Ce sont la **doc et l'assert** qui ont
dérivé.

**Correctif.** Mettre à jour README + commentaires + bannière (1 m, pas nair==1) ;
pointer l'assert sur la branche réellement publiée (z05, 44/9) ou documenter
pourquoi deux branches coexistent.

### 2.5 [MOYEN] Deux r distincts (0.92 et 0.93) présentés comme "le r du classement"

- `07_validation.R:18` : r=0.92 = Pearson sur **ΔTmax moyen par plot** (n=53).
- `scripts/make_slope_sim_vs_obs_groups.R:54` : r=0.93 = Pearson sur la **pente
  micro/macro par plot** (n=53). Quantité différente.

Le papier/README les emploie de façon interchangeable. **Correctif** : préciser
quelle figure rapporte quel r, et sur quelle métrique.

### 2.6 [MOYEN] Validation forward (Fig B1) : r sur n=8 sans n affiché ni garde-fou

`07b_forward_validation.R:49-66` calcule et stocke `n` par (cluster, step) mais
**n n'est jamais dessiné**. P1 (n=8, LAI r≈0.84) a le même poids visuel que All
(n=53). Aucun IC, aucun caveat small-n. Un r sur n=8 est quasi non informatif.

**Correctif.** Annoter n sur chaque cellule, mettre un caveat small-n (déjà
partiellement dans le texte de l'annexe B, à renforcer dans la figure).

### 2.7 [MOYEN] Non reproductible depuis un checkout propre

Stage 03 est **skip-if-exists** sur des NetCDF cachés (848 sims HOBO déjà
présents). Les **6400 sims cLHS** de l'attribution ne sont **pas dans `run_all.R`**
(disque externe `NC_Full/`, script `c3_*`). La régénération exige le **binaire
legacy hard-codé** `/home/corroyez/Documents/musica/musica` (md5 5307…) ; le
binaire alternatif `model-3.2.3` décale ΔTmax jusqu'à 1.4 °C et inverse le ratio
buf/amp (`R/config.R:13-16`). Donc science non reproductible sans ce binaire exact
+ le disque externe.

**Correctif.** Documenter ces dépendances dans le README ("non régénérable sans
binaire X + disque Y") ; idéalement versionner le binaire (ou son hash + source)
et rapatrier les NetCDF cLHS sous un chemin du dépôt.

---

## TIER 3 — Points solides & corrigé cette session

### Solide (crédit)
- **Formule Shapley exacte correcte** : poids `s!(n−s−1)!/n!`, additivité
  Σφ = v(1111)−v(0000) assertée à l'identité machine (`05_shapley.R:15-25,51`).
- **LAD uniforme conserve l'aire foliaire totale** (change la forme, pas l'intégrale) :
  `R/lad.R:13,43-52` renormalise à la même LAI → le toggle LAD réel↔uniforme est propre.
- **Doublement LAI legacy correct** : PAD→LAI (`io.R:44`) puis one-sided→two-sided
  (`musica.R:69`), deux ×2 distincts appliqués une fois chacun sur la voie legacy z1.
- **Typologie gelée garantie** : `02_clusters.R` ne re-cluster pas par défaut ;
  `PIPE_RECLUSTER=TRUE` re-tourne au seed 42 et asserte ARI>0.95 vs cache.
- **Matching validation propre** : `04_extract.R` lit des MuSICA dédiés AUX 53
  positions logger (`musica_hobo_*`), jamais interpolés depuis les plots cLHS.
- **Extraction 1 m sans off-by-one**, conventions abs (1 m + −2 h) et pente
  (1 m, no shift) bien séparées (`04_extract.R:33-59`).

### Corrigé cette session
- **Phrase "footprint 12.5 m" en Méthodes 2.2** : les traits principaux viennent
  de **pixels 20 m** (`R/config.R:20` `agg_factor=2`, 10 m→20 m), PAS d'un footprint
  12.5 m. L'étude de rayon (Annexe C) ré-extrait depuis le **nuage de points brut**
  (voie différente). Méthodes 2.2 corrigé : traits agrégés à 20 m + Annexe C
  reframée en simple test d'insensibilité au footprint.

### Placeholders restants
- `:125` `[REF: Ramsay & Silverman]` (FPCA) ; `:291` `[REF: Hersbach et al. 2020, ERA5]`.

---

## RÉSOLUTIONS (session 2026-06-15)

- **1.1 RÉSOLU.** Covariables résidu basculées Oak_Only → **Not_Masked** (les 9 NA
  étaient les placettes hors masque chêne = les amplificatrices). Régression sur
  **n=53**, **adj R²=0.62** (vs 0.41/n=44), `gap_fraction` **r=0.76** (vs 0.58),
  partial p<0.001 ; amplificatrices gap 0.52 vs 0.01 (p=0.006). Inclure les
  amplificatrices RENFORCE le résultat. `scripts/make_residual_vs_topo.R` + garde NA.
- **1.2 RÉSOLU.** Tableau dispersion (Appendice D) + défense H2 ajoutés au manuscrit.
  LAI SD plat (1.25–1.52) → effondrement importance = saturation physique.
- **1.3 RÉSOLU.** Double-doublement supprimé dans `03z` (×2 uniquement dans
  `musica.R`). z05 re-tourné (848 sims). Validation **inchangée** (r=0.931, biais
  +0.92, 44/9, pente 0.945 vs 0.855) : le LAI réel du chêne est déjà saturé, donc
  8.5 vs 17 ne change pas ΔTmax — cohérent avec la thèse de saturation. Inputs
  désormais physiquement corrects.
- **2.1 RÉSOLU.** `pipeline/11_clhs_attribution.R` (canonique, local, floor05_v2)
  régénère Fig 2 (reproduit exactement : P1 LAI 0.41 → P4 LAD 0.13) et Fig 3.
  `c3_shapley_clhs_2x.R` (floor05) DÉPRÉCIÉ (stop()). Merge redirigé hors NC_Full.
  **Découverte :** l'ancienne Fig 3 signée "≈-0.2" venait du mauvais floor05 ;
  sous floor05_v2 le rang d'importance tient (LAI 0.54 > fCover 0.31 ≫ LAD 0.09)
  mais le φ signé enjambe zéro (off-manifold) → Fig 3 passée en mean|φ|, texte corrigé.
- **2.2 RÉSOLU.** `scripts/make_moran_validation.R` (spdep, kNN=8, 999 perms).
  **Correction :** observé I≈0.10 (pas 0.30) ; modèle PLUS autocorrélé (I=0.23),
  pas "amorti" → le modèle sur-lisse (cohérent glass ceiling) ; résidu faiblement/
  pas autocorrélé. Texte 3.3 réécrit.
- **2.3 RÉSOLU.** Fig 2/3 → local floor05_v2 ; Fig 4/5 → z05 (Fig 5 lisait déjà z05).
- **2.4-7 RÉSOLU.** "nair==1 / 45-8" purgé de README, config (bannière+commentaires),
  en-têtes de stages. Numéros validation z05 (44/9, r=0.93) confirmés = manuscrit.
  Assert legacy 45/8 conservé (garde la branche v10 ; branche publiée = z05, documenté).

## Top 5 à traiter avant le jury / soumission (ÉTAT)

1. **Résidu (1.1)** : re-tourner sur 53 ou avouer n=44 + amplificateurs exclus ;
   R² ajusté ; gap_fraction sur 3 points.
2. **Moran's I (2.2)** : committer le code ou retirer la prétention.
3. **Provenance attribution (2.1, 2.3)** : régénérer Fig 2/3 depuis floor05_v2,
   figer le chemin, supprimer le script floor05 collisionnant.
4. **z05 double-LAI (1.3)** : confirmer (~17) et corriger ; re-rapporter 44/9.
5. **Dispersion / mean|φ| (1.2)** : déjà réfuté par les données — ajouter le
   tableau + le lien H2 contrôlé comme défense écrite.
