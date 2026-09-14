# Plan de retouches — Chapitre 1 (2026-09-07)

Capture des décisions et analyses issues de la revue du 07/09. Le Ch1 n'est **pas
finalisé** : prévoir explicitement du temps de retouche (absent du planning actuel),
faire relire par les encadrants, et soumettre les arbitrages de structure au CSI.

Version courante = `manuscript_chap1_EN_native20_csi3pass_2026-09-02.md`
(5 figures main : Fig 1–5 ; 6 annexes : A, B, E, F, G, H).
Version complète archivée = `manuscript_chap1_EN_native20.md` (8 main, 9 annexes)
— source pour rapatrier ce qui a été coupé (radius, glass ceiling, model-free).

## 1. Recadrage éditorial

- **Emphase déplacée** : le Ch1 n'insiste pas sur l'accord mesure-modèle mais sur
  l'**analyse de sensibilité du modèle aux variables de structure** (Fig 4/5 = cœur).
- **Validation HOBO (ex-Fig 6 glass ceiling, ex-Fig 7 model-free) → Chapitre 3.**
  Ce sont des figures de validation contre les 53 loggers ; leur place est au Ch3.
  Déplacement à confirmer avec le CSI (« éventuellement déplaçable »).
- Conséquence : la métrique radius/footprint doit être refaite **sans** cible
  observationnelle (voir §4).

## 2. Renvois vers le Chapitre 3

- Aller plus loin, mono-site, **à partir du microclimat** (réf : un article de
  Maclean à retrouver) : générer des températures type thermocouple HOBO,
  approche **bilan d'énergie**.
- Les figures de validation HOBO déplacées du Ch1 y sont réintégrées.

## 3. Analyses scientifiques à mener (Ch1, profils verticaux)

1. **Milieu ouvert (P1) — modèle vs mesure.** Trancher si l'absence de skill en
   canopée ouverte vient du modèle (1-D, near-neutral) ou des mesures terrain
   (HOBO contaminés sol/rayonnement en trouée). Prudence d'interprétation requise.
2. **Profils verticaux de température par LAI.** Ils doivent différer selon le LAI.
   Question directe et testable : **se rejoignent-ils à 1 m** (hauteur logger) ?
3. **Ne pas s'arrêter à la validation HOBO.** Sortir les profils simulés :
   T haut de canopée, T air, T feuillage (moyenne), T sol. À défaut, rester
   **prudent** et le discuter explicitement.
4. **Position relative du pic de LAI vs capteur.** Ajouter du LAI **au-dessus**
   vs **en-dessous** du capteur (1 m) et mesurer l'effet sur ΔTmax.
5. **Normalisation des boxplots.** Tester une normalisation par **ΔTmax / ΔVar**
   (révèle-t-elle plus de structure ?). Et **dénormaliser la Fig B3**.

## 4. Radius / footprint (ex-Fig 8) — refonte sur cLHS

Ancienne version : corrélation ΔTmax vs 53 loggers à 7 rayons. Remplacée par :

- **Base = 400 plots cLHS** (les 53 HOBO partent au Ch3).
- **Référence = footprint 20 m** (pixel natif actuel).
- Pour chaque rayon r ∈ {5, 10, 12.5, 15, 25, 50 m} : re-clip du nuage →
  recompute des 4 variables (LAI, fCover, Hmax, profil LAD) → re-run MuSICA →
  **écart ΔTmax(r) − ΔTmax(20 m)** par plot. Des **écarts de ΔTmax**, pas des
  corrélations. (20 vs 20 = 0, référence.)
- Livrables complémentaires :
  - **gamme des 4 variables structurelles par rayon** (stabilité de l'extraction) ;
  - re-run de la **perturbation Fig 4/5 par rayon** (robustesse de l'attribution
    au choix de footprint).

⚠️ Charge : 400 × 6 rayons ≈ 2400 runs MuSICA (+ perturbations × rayons).
**Dry-run obligatoire** : ~2 rayons × ~20 plots d'abord, estimer le temps, valider,
puis full run.

### Plots écartés aux petits rayons (à discuter, pas à imputer)
Au rayon 5 m, 4 plots sur 400 ne donnent pas de canopée exploitable, pour deux
raisons distinctes (vérifié par re-clip) :
- **1 trouée écologique réelle** (clhs_243, P1) : points présents mais aucun retour
  de végétation au-dessus de 2 m, donc LAD = 0 ; rien à forcer dans MuSICA.
- **3 trous de couverture LiDAR** (clhs_264, 287, 380) : le centre du plot tombe
  dans un gap de la donnée (ligne de vol / no-data, ~7 à 15 m de rayon) ; 0 point
  jusqu'à 10 m, données pleines à 20 m. Leur valeur « 20 m » vient de la fraction
  couverte de la cellule. Ce n'est pas de l'écologie, c'est la donnée.

Traitement : la comparaison inter-rayons (robustesse de l'attribution) se fait sur
le **jeu commun de plots présents à tous les rayons** (échantillon constant), et on
**reporte le n par rayon**. Pas d'imputation. Effet confiné aux petits footprints ;
aux rayons ≥ 10 m (sauf clhs_264) les 400 plots passent. À mentionner en une phrase
dans le chapitre.

## 5. Logistique

- [ ] Prévoir un créneau de retouche Ch1 (actuellement non planifié).
- [ ] Lecture encadrants.
- [ ] Soumettre au CSI l'arbitrage « validation HOBO déplacée au Ch3 ».
- [ ] Retrouver la référence Maclean (bilan d'énergie / génération T HOBO).

## Ordre d'attaque proposé

a. Radius sur cLHS (le plus cadré) — §4.
b. Profils verticaux par LAI + « se rejoignent-ils à 1 m » — §3.2/3.3.
c. LAI au-dessus/en-dessous du capteur — §3.4.
d. Normalisation boxplots + dénorm B3 ; ouvert modèle-vs-mesure — §3.5/3.1.
