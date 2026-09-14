# Triage des retours encadrants / CSI (collés 2026-09-10)

Les commentaires collés mélangent Ch1, Ch3, analyses futures et planning. Triés ci-dessous.
Seuls les items **[Ch1-TEXTE]** sont traités dans la passe Fable en cours ; les analyses attendent
un dry-run et ton go (règle CLAUDE.md : dry-run obligatoire avant calcul lourd). Les items **[Ch3]**
sont à reporter dans le magasin mémoire Ch3 (dossier NC_Full).

> Fait de structure qui tranche plusieurs items : **le Ch1 n'a PAS de Fig B3** (figures = 1–7, A1/A2,
> B1, C1, D1) et **aucun contenu HOBO/logger** (déplacé au Ch3). Donc « fig B3 », « validation HOBO »,
> « figure validation HOBO » = Ch3, pas Ch1.

---

## [Ch1-TEXTE] — traités dans la passe Fable (édition maintenant)

1. **« pas d'emphase accord mesure/modèle, + sur l'analyse de sensibilité du modèle aux variables de
   structure ».** = la consigne la plus nette. Le Ch1 n'a aucune validation terrain ; toute phrase qui
   se lit comme « le modèle est d'accord avec l'observation » est à corriger. Points chauds : abstract,
   conclusion, §4.3 (réfs étude compagnon 2026a), Annexe A (biais ERA5-vs-station RMSE 2,06 / r 0,92 =
   validation du FORÇAGE, pas du modèle → risque de lecture erronée). → **vérifié par les agents Fable**.
3. **(volet Ch1) milieu ouvert : prudence P1.** Le §4.3 dit déjà « colonne 1D à fidélité limitée en
   canopée ouverte » et « pas de validation terrain ». Vérifier que l'interprétation ailleurs ne va pas
   plus loin que cet aveu. Le fond « modèle vs mesures » est surtout un sujet Ch3 (validation).
7. **(volet prudence) « dommage de s'arrêter à la validation HOBO → sinon être prudent et en discuter ».**
   Le Ch1 n'ayant pas de HOBO, la version Ch1 = rester prudent sur le fait que la métrique est un ΔTmax
   à 1 m sans profil vertical ni température terrain derrière. Vérifier le ton du §4.3 et de la discussion.

---

## [Ch1-ANALYSE] — analyses futures, dry-run + go requis (NON lancées)

4. **Profils de température verticaux selon le LAI ; se rejoignent-ils à 1 m ?** Relié à l'Annexe B
   (« within-canopy gradient, per variable »). Analyse : sortir les profils T°(z) simulés par archétype/
   LAI et regarder la convergence à 1 m (la métrique). Valide/qualifie le choix du capteur à 1 m.
5. **Normaliser les boxplots par ΔTmax/ΔVar** (Fig 4/5) pour « voir + de choses ». Idée figure : passer
   les effets en fraction de la variation de la variable → comparabilité inter-archétypes. À prototyper.
7. **(volet analyse) Analyser le profil complet** : T° haut de canopée, air, feuillage (moyenne), sol,
   pas seulement le ΔTmax à 1 m. Renforcerait l'interprétation mécaniste (ou impose la prudence si non fait).
8. **Position relative du pic de LAI vs capteur ; ajouter du LAI au-dessus / en-dessous.** Extension
   directe du test H2 gaussien (sweep μ,σ) : ici cibler explicitement le rapport pic–capteur (1 m) et une
   redistribution asymétrique. Recoupe l'item 4.
11. **Rayon d'emprise (radius) c1 : corréler à une simulation de référence** (20 m par défaut) ; comparer
    5 m vs 20 m, 10 m vs 20 m ; **exprimer en ΔTmax, pas en corrélations**. Teste la sensibilité du
    résultat à l'échelle de la fenêtre LiDAR.
12. **Rayon c1 : passer des 53 HOBOs aux 400 placettes cLHS** pour ce test (si les HOBOs migrent au Ch3).
    Cohérent avec le périmètre Ch1 = 400 cLHS uniquement.

---

## [Ch3] — à reporter dans le magasin Ch3 (PAS Ch1)

2. **Aller plus loin sur un seul site, à partir du microclimat** (cf un article Maclean) : générer des
   températures type thermocouple HOBO et un bilan d'énergie. = programme Ch3.
3. **(volet Ch3) milieu ouvert : est-ce le modèle ou les mesures terrain le problème ?** Question de
   validation Ch3 (placettes ouvertes P1).
6. **Fig B3 : dénormaliser ?** = figure Ch3 (pas de B3 au Ch1).
9. **Figure de validation HOBO : dans le chapitre 3 ?** Placement, Ch3.

---

## [PLANNING] — à noter, pas d'action de rédaction

10. Prévoir du temps pour retoucher le Ch1 (actuellement non prévu). Voir avec le CSI leur avis.
    Validation HOBO éventuellement déplaçable. **Ch1 pas finalisé, lecture encadrants à faire.**

---

## Résultat de la passe extensive écriture/cohérence/citation/consignes (2026-09-10)

Passe demandée « Fable extensive ». Les 4 agents Fable ont échoué (limite de session sur le modèle
`claude-fable-5-1`, reset 19:20 Paris) → passe faite par Opus avec le skill `nathan-scientific-voice`
chargé en entier, section par section + global, mêmes consignes que les agents (POLISH, before/after
ciblés, jamais réécriture de paragraphe ; citations = style only, vérif plein-texte déjà faite).

**Verdict : le chapitre est propre.** Rien à éditer qui vaille le risque sur un texte déjà passé 4 fois.
- **Voix (balayages mécaniques sur toute la prose)** : 0 cadratin ; 0 verbe portentous (underscore/
  highlight/reveal) ; 0 métaphore positionnelle (« sits/stands » = tous des « stands » = peuplements) ;
  0 participe traînant ; 0 clivage de révélation ; 0 phrase de prose > 35 mots (les lignes longues =
  légendes de tableaux auto-suffisantes, règle Féret, ou retours-ligne). Anglais britannique intact.
- **Cohérence sections + global** : Key message/Abstract/Méthodes/Résultats/Discussion/Conclusion/
  Annexes A–F relus ; cadrages verrouillés tenus (quantité domine, profil second-ordre, jamais LAI-vs-
  cover, « P2 to P4 », σ axe contrôlé dominant sauf LAI 4, LAI 6 = bord, Tables 2/E1 médianes, pas de
  HOBO). L'édit A1 (« primacy often attributed ») recouple bien avec §4.1 L478.
- **Consigne « pas d'emphase accord mesure/modèle »** : RESPECTÉE. Le texte désavoue explicitement la
  validation terrain (§2.6 L308, §4.1 L519, §4.3 L575/L627, Annexe B L703, L709-710). Annexe A r=0,92/
  RMSE 2,06 = biais ERA5-vs-station (forçage), légende claire. Abstract dit « simulated buffering ».
- **Consigne « prudence milieu ouvert / ΔTmax 1 m seul »** : RESPECTÉE. §3.2 L354-355 dit P1
  « weakly constrained here » + renvoi étude compagnon 2026a ; Annexe B L703 « model constructs, with no
  multi-height measurement to anchor them ».

**Non-édité volontairement** : légendes Tables 1 (103 mots) et 2 (116 mots) dépassent le budget légende
figure (≤60 mots) mais sont des légendes de TABLEAU auto-suffisantes (convention de signe, troncature,
comptes de placettes) → règle Féret « la légende porte sa propre définition » > budget. Ne pas raccourcir.

**Option** : relancer les 4 agents Fable après 19:20 pour un second regard indépendant du modèle Fable.
