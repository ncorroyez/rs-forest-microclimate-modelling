---
title: "Point encadrants pré-CSI3 (enrichi) — 26 août 2026"
date: 2026-08-26
attendees: Nathan Corroyez, Sylvie Durrieu, Jérôme Ogée, Jean-Baptiste Féret
duration: 9 h 00 – 11 h 00 (visio, lien envoyé par Nathan)
---

# Point encadrants pré-CSI3 (enrichi) — 26 août 2026

## 0. Statut en une ligne

Rapport CSI3 v3 (10 p), slides (20 diapos, notes orateur sur les 20), manuscrit WIP (134 p, Ch1 + Ch2 pleins, intro et Ch3 en juste milieu, discussion en outline) sont prêts. Ch1 AoFS figé chez encadrants, Ch2 r2 chez éditeur RSE, Ch3 draft structurellement complet, en attente lecture encadrants (gel 30/09). La réunion sert à trancher 4 décisions bloquantes avant le CSI du 1ᵉʳ septembre.

---

## 1. Ordre du jour timeboxé (2 h)

| Temps | Sujet | Livrable |
|---|---|---|
| 09:00 – 09:10 | Ce que je veux sortir avec (§9) | Cadrage |
| 09:10 – 09:30 | D1 RSE fallback + D2 Ch3 saisonnier vs régime-dépendant | Position commune CSI |
| 09:30 – 09:50 | D3 réinscription 4A + D4 financement 2 mois ANR | Feu vert Direction |
| 09:50 – 10:15 | Ch1 AoFS : retours consolidés + fenêtre de soumission | Date cible |
| 10:15 – 10:40 | Ch3 : plan de finalisation Ch3 avant 30/09 + réponse à la remarque Jérôme sur S2/PAI | Priorités |
| 10:40 – 11:00 | Rapport + slides : dernières remarques, ouvertures manuscrit | Version à figer |

---

## 2. Les 4 décisions à trancher

### D1 — Stratégie article RSE si retour défavorable
- **Position Nathan** : ne pas re-formater avant la réponse RSE, mais pré-identifier ce matin un journal de repli (ISPRS J. Photogramm. RS, IEEE JSTARS, Int. J. Applied Earth Obs. Geoinf.) pour économiser 2-3 semaines si rejet.
- **Pushback Jérôme (mail 25/08)** : discussion prématurée, rien à décider avant retour éditeur.
- **Fallback si Jérôme confirme sa position** : accepter, mais me faire préciser publiquement ce matin quel journal Jérôme viserait le jour où on saurait ; ça vaut identification implicite.

### D2 — Ch3, régime-dépendant vs saisonnier S2 (le point Jérôme)
- **Position Nathan** : le régime-dépendant (LiDAR dense, S2 open, split LAI 3.86) reste la réponse déployable ; le saisonnier S2 est **testé en Annexe A** du Ch3 et n'ajoute rien sur ΔTmax à Blois en été (null result documenté).
- **Remarque Jérôme (mail 25/08)** : les écarts LiDAR-S2 dans MuSICA viennent probablement du PAI S2 plus faible → mécaniquement meilleur accord en ouvert, pire en dense. Suggère : utiliser S2 pour le saisonnier.
- **Ce qu'il faut clarifier avec lui ce matin** : sa suggestion vise (a) une extension inter-annuelle Blois multi-années, (b) un gap-filling entre deux vols LiDAR, ou (c) autre chose. Le Ch3 tel qu'écrit n'a pas testé (a) ni (b), qui sont deux perspectives valides. Cadrer sa remarque comme perspective de discussion générale, pas comme refonte Ch3.
- **Fallback** : si Jérôme veut une refonte Ch3, chiffrer le coût en semaines. Impact direct sur D3 (gel 30/09 devient intenable).

### D3 — Réinscription en 4A
- **Position Nathan** : réinscription nécessaire, trois causes documentées (rapport §7.1) + deux décès familiaux au printemps 2026.
- **Attendu des encadrants** : validation explicite pour la présenter comme telle au CSI (« la Direction soutient la demande »).
- **Fallback** : si Sylvie ou Jérôme veut nuancer publiquement, obtenir au moins un accord de principe et une phrase-clé qu'ils accepteront de dire au comité.

### D4 — Financement ANR MaCCMic post-décembre 2026
- **Position Nathan** : besoin des 2 mois (janv+févr 2027) pour tenir le scénario 1 (dépôt mi-déc, soutenance février). Sans cela, scénario 2 ou 3, 5 semaines de manuscrit en moins.
- **Attendu de la Direction** : engagement à ouvrir la démarche cette semaine (contact avec le PI ANR MaCCMic).
- **Fallback** : si pas possible, cadrer le CSI directement sur le scénario 2 (dépôt fin novembre) et anticiper la charge d'écriture.

---

## 3. Anticipation des questions et objections

### 3.1 Sylvie (issues des annotations `_sd` et de son style)

| Question probable | Réponse à donner | Où c'est |
|---|---|---|
| « Tu as bien remplacé "world" par "fields" en §1.1 ? » | Oui, rapport §1.1 revu au stop-slop. | Rapport §1.1 |
| « L'espèce du chêne à Blois est bien précisée ? » | Blois = *Quercus petraea* (sessile), Mormal = *Q. robur* + *Fagus sylvatica*. | Rapport §2 |
| « Les séries temporelles S2 du Ch3 sont dans la description des données ? » | Oui, §2 mentionne « full 2021–2022 dynamic Sentinel-2 LAI series at Blois ». | Rapport §2 + slide 5 |
| « Les archétypes P1-P4 sont introduits dans le Design de Ch1 avant d'être utilisés ? » | Oui, §3 Design décrit K-means sur 6 métriques ALS + cLHS 100 plots/archétype = 400 plots. | Rapport §3 + slide 7 |
| « Tu as ajouté les conférences (EGU, MT180) ? » | Oui, §8 les liste. | Rapport §8 + slide 18 |
| « Le calendrier reflète les cycles d'AR encadrants + co-auteurs ? » | Oui, Table 4 caption dit explicitement que chaque freeze wrappe un cycle supervisor + co-auteur. | Rapport §7.1 |
| « GEDI c'est bien pour FORMS-H seulement ? » | Non, dans le rapport GEDI-based comparison est distincte de FORMS-H (qui n'utilise GEDI que pour calibration). Reformulé en §7.3. | Rapport §7.3 |
| « Tu as bien nuancé la tolérance ED GAIA ? » | Oui, §7.2 dit « usually covers... not a fully guaranteed line ». | Rapport §7.2 |
| « Ton contrat n'est validé que sur Ch1 par tes encadrants ; tu as intégré le temps d'AR ? » | Oui, note explicite sous Table 4. | Rapport §7.1 |
| « Le style est encore un peu robotique. » | Reconnaître, expliquer que la passe stop-slop a tourné cette semaine sur intro + Ch1 + Ch3 + discussion (Ch2 laissé intact) ; scores 38-40/50. | Voir /tmp/.../fable_*_report.md |
| Sylvie sur R²=0 LiDAR Ch3 : « Tu as réintégré les pixels non couvrants ? » | Oui, les 53 loggers couvrent fCover 0.4-1.0. Le 0.04 (dense S2) reflète la saturation PROSAIL, pas une exclusion. Note ajoutée sous Table 3. | Rapport Table 3 caption |

### 3.2 Jérôme (mail 25/08)

| Question probable | Réponse à donner | Où c'est |
|---|---|---|
| « Le style est difficile à comprendre. » | Reconnaître, dire que la passe stop-slop a tourné. Ne pas défendre le style. | — |
| « Les archétypes P1-P4 sont décrits avant d'être utilisés ? » | Oui, ajouté §3 Design. | Rapport §3 + slide 7 |
| « Boundary-layer height n'est pas utile même en iter, pourquoi c'est là ? » | Retiré. Le forçage §3 Ch1 dit maintenant simplement « CHS 41 station + ERA5 secondaires ». | Rapport §3 |
| « C'est quoi "wind-corrected native-resolution product" ? » | Reformulé : « native 20 m ALS product forced with per-plot wind-corrected macroclimate (log-profile 10m → sommet canopée) ». | Rapport §3 |
| « Table 2, la corrélation ΔTmax est model-observation pas between-plot. » | Relabelée : `Model-observation ΔTmax correlation (n = 53)`. L'amplitude reste `Between-plot`. | Rapport Table 2 |
| « Je n'arrive pas à réconcilier tes r Ch1 et R² Ch3. » | Ajouté note sous Table 2 : R² ≈ r², r=0.69 dense Ch1 → R²=0.48 ≈ 0.56 Ch3 dense LiDAR. Différence = scénario LAI testé (STATIC_ALS vs Ch1 native), pas la métrique. | Rapport Table 2 caption |
| « Ne discutons pas fallback RSE avant retour éditeur. » | Accepter. Question 6 déjà reformulée : « no decision needed before RSE editor's answer ». Ne pas insister. Rester ferme sur le principe d'identification implicite (voir D1). | Rapport §4 + Q6 |
| « Les différences LiDAR-S2 dans MuSICA sont dues au PAI S2 plus faible. Utilise plutôt S2 pour le saisonnier. » | (a) Le PAI S2 plus faible est un fait diagnostiqué au Ch2 (d_opt = 7 m pour Blois → S2 voit moins). C'est bien pour ça que le régime-dépendant marche : LiDAR pour la vraie densité en dense, S2 en ouvert où sa saturation ne joue pas. (b) Le saisonnier S2 a été testé en Annexe A Ch3 → null result à Blois en été (dynamique ≈ statique). (c) Sa suggestion est-elle sur (i) le multi-années, (ii) le gap-filling entre vols LiDAR, ou (iii) une autre piste ? À clarifier. Cadrer comme perspective de la discussion générale, pas comme refonte Ch3. | Ch3 §4 + Annexe A |

### 3.3 Jean-Baptiste Féret (silencieux depuis 5/07)

| Question probable | Réponse à donner | Où c'est |
|---|---|---|
| « Où en est le r2 RSE côté PROSAIL priors ? » | 225 configurations testées (grille 5 × 5 × 3 × 3 sur ALA, LAI, LMA, BROWN). Partial-attribution table borne la circularité prior LiDAR à ≤ 17 % du RMSE reduction. Chez éditeur depuis 10/08. | Rapport §4 + Ch2 r2 |
| « La retrieval "opt" forestière du Ch2 est bien intégrée au Ch3 ? » | Oui, `STATIC_S2_DOPT` et `DYN_S2_DOPT_LADOPT` testent explicitement l'opt. L'opt améliore l'accord numérique LAI vs LAI-tronqué, mais **inverse** le ranking (r LAI-space +0.18 vs −0.18 pour ATBD brut). C'est le contre-résultat du Ch3. | Ch3 §3.4 + rapport §5 |
| « Le preprocS2 est bien cité ? » | Oui, §2 rapport + Ch2 méthodes. | Rapport §2 |
| « Le forçage MuSICA en Ch1 est basé sur quelle station ? » | CHS 41 (RENECOFOR) pour Tair/RH/precip + ERA5 gap-fill pour radiation, vent, pression. | Rapport §3 + Ch1 §2.4 |
| « Silvilaser 2025, MEB 2026, EGU 2025 sont dans le rapport ? » | Oui, §8 les liste + MT180. | Rapport §8 |
| « Ta stratégie de rédaction sur les 5 semaines Ch3 ? » | Gel 30/09. Reste à finaliser : Declarations/Author contributions/Funding à compléter, passe éditoriale, cycle d'AR encadrants. : finaliser §3.4 (fusion scenarios) + §4 discussion + §5 conclusion + appendices A-D. Cycle d'itération : draft complet mi-septembre, retour encadrants fin septembre. | Rapport §7.1 + Table 4 |

---

## 4. Chiffres à avoir sous la main

| Chapitre | Grandeur | Valeur |
|---|---|---|
| Ch1 | *r* ranking (n=53) | 0.50 (CI 0.26-0.68) |
| Ch1 | Biais chaud | +0.98 °C |
| Ch1 | Dense stratum (n=27) | *r* = 0.69 (CI 0.42-0.85) |
| Ch1 | Sparse stratum (n=26) | *r* = 0.20 (n.s.) |
| Ch1 | Slope correlation | *r* = 0.78 |
| Ch1 | Amplitude recovered | ≈ 15 % |
| Ch1 | Design | K-means (6 métriques ALS) → 4 archétypes P1-P4 ; cLHS 100/archétype → 400 plots |
| Ch1 | Extinction | k = 0.5 |
| Ch2 | d_opt (Aigoual / Blois / Mormal) | 8 / 6 / 10 m |
| Ch2 | Circularité prior LiDAR | ≤ 17 % du RMSE reduction |
| Ch2 | LUT PROSAIL | 225 configurations (5 × 5 × 3 × 3) |
| Ch2 | Hétérogénéité horizontale | DSM_SD (plus fort covariant du résidu) |
| Ch2 | Extinction | k = 0.65 |
| Ch2 | 5 sensibilités reviewer | k, angle zénithal, h_min veg, fCover, d_thr ∈ {10,15,20 m} |
| Ch3 | Split LAI 3.86 | dense n=25, open n=28 |
| Ch3 | Dense LiDAR vs S2 | R² 0.56 vs 0.01 |
| Ch3 | Open LiDAR vs S2 | R² 0.44 vs 0.63 (CIs overlap) |
| Ch3 | Pooled LiDAR vs S2 | R² 0.38 vs 0.49 (indistinguables) |
| Ch3 | Régime-dépendant | pooled R² 0.55, biais +0.84 °C, RMSE 1.41 °C |
| Ch3 | FORMS-H switch 17.8 m | R² 0.51 in-sample |
| Ch3 | d_opt Blois utilisé | 7 m (cross-site combined, dans la fenêtre 6-10 m de Ch2) |
| Ch3 | Extinction | k = 0.65 |
| Général | Sites | Aigoual (hêtre), Blois (chêne sessile), Mormal (chêne pédonculé + hêtre) |
| Général | Loggers HOBO | 180 posés (60/site), 53 utilisables Blois été 2021 |
| Général | Forçage MuSICA | CHS 41 (Tair/RH/precip) + ERA5 (radiation/vent/pression) |
| Général | Convention ΔTmax | time-matched (Bouwen 2025) |

---

## 5. Follow-up des retours antérieurs

| Item | Origine | État | Note |
|---|---|---|---|
| D1 5/7 — contradiction modèle/obs en canopée ouverte | 5/7 | Fait | Assumé (Option A) : Ch1 §3.3-4.3 documente, sparse r=0.20 n.s. dit dans rapport CSI3. |
| D2 5/7 — H1 rejetée (profil vertical négligeable) | 5/7 | Fait | Profil « puissant mais non réalisé » ; Ch1 abstract + discussion refondus. |
| D3 5/7 — design de perturbation pas-fixe | 5/7 | Fait | Pas natif LAI ±0.5, fCover ±10 pts, Hmax ±1 m, profil swap réel-vs-uniforme. |
| Ch2 recadrage accuracy → consistency | 5/7 | Fait | Recadrage r2, resoumis 10/08. |
| Passe Fable coherence CSI3 | interne 25/08 | Fait | 3 fixes high-prio + 3 med-prio appliqués. |
| Sylvie _sd : 21 corrections | Sylvie 25/08 | Fait | Toutes appliquées dans rapport v3. |
| Jérôme mail 25/08 : Table 2 label | Jérôme 25/08 | Fait | Renommée `Model-observation ΔTmax correlation`. |
| Jérôme mail 25/08 : wind-corrected clarifié | Jérôme 25/08 | Fait | §3 reformulée. |
| Jérôme mail 25/08 : RSE fallback assoupli | Jérôme 25/08 | Fait | §4 + Q6 reformulés « no decision before RSE answer ». |
| Jérôme mail 25/08 : boundary-layer height retiré | Jérôme 25/08 | Fait | §3 nettoyée. |
| Jérôme mail 25/08 : réconciliation r Ch1 / R² Ch3 | Jérôme 25/08 | Fait | Note sous Table 2 caption. |
| Jérôme mail 25/08 : style robotique | Jérôme 25/08 | Fait | Passe stop-slop sur rapport + intro + Ch1 + Ch3 + disc. |
| Jérôme mail 25/08 : Ch3 saisonnier vs régime | Jérôme 25/08 | Ouvert | D2 aujourd'hui. |
| Poster manuscrit WIP sur Nextcloud | Sylvie 5/7 | En attente | À faire 27/08. |

---

## 6. Fichiers à ouvrir par point d'agenda

| Point d'agenda | Ouvrir |
|---|---|
| Cadrage + statut | `/home/corroyez/Documents/Admin/ED_GAIA/CSI3/V3_Corroyez_Nathan_CSI_3A_4A_Report.pdf` §1-2 |
| D1 RSE fallback | Rapport §4 + §9 Q6 |
| D2 Ch3 saisonnier | Rapport §5, table 3 + `/home/corroyez/Documents/NC_Full/manuscripts/ch3/Chapter3_article_standalone_EN.md` §3.4 et Appendix A |
| D3 réinscription | Rapport §7.1 + §9 Q6 (renseigner la formulation retenue) |
| D4 financement | Rapport §7.2 + slide 17 |
| Ch1 AoFS | `/home/corroyez/Documents/z_Example_rmusica_31012025/chapter1/manuscript/manuscript_chap1_EN_native20_AoFS.docx` |
| Ch3 plan 5 semaines | Ch3 article md + Table 4 rapport |
| Rapport / slides | Rapport v3 + `/home/corroyez/Documents/Admin/ED_GAIA/CSI3/Corroyez_Presentation_CSI_3A_4A_01_09_26.pptx` |
| Manuscrit WIP | `/home/corroyez/Documents/NC_Full/thesis_current_state_wip.pdf` (134 p) |

---

## 7. Actions attendues à l'issue

| Action | Owner | Date cible |
|---|---|---|
| Poster manuscrit WIP sur Nextcloud, envoyer lien au comité | Nathan | 27/08 |
| Version finale rapport + slides après retours réunion | Nathan | 29/08 |
| Envoi rapport + slides + lien Nextcloud au comité | Nathan | 30/08 |
| CSI 3A → 4A | Comité | 01/09 |
| Ch1 AoFS retours consolidés par les 3 encadrants | Sylvie, JBF, Jérôme | Date à fixer ce matin |
| Ch3 gap 40 % rédigé | Nathan | 30/09 |
| Ch3 relectures encadrants | Sylvie, JBF, Jérôme | 30/09 → 15/10 |
| Démarche sécurisation 2 mois ANR MaCCMic | Direction (Sylvie + Jérôme) | Cette semaine si accord ce matin |
| Intro générale + discussion générale gelées | Nathan | 22/11 |
| Dépôt manuscrit ADUM (scénario 1) | Nathan | 15/12 |

---

## 8. Blockers / risques ouverts

- **RSE** : verdict pas garanti avant le dépôt manuscrit. Fallback identifié implicitement mais pas engagé.
- **Ch3** : 40 % restant à écrire en 5 semaines, avec 3 relectures encadrants dans la même fenêtre. Si Jérôme demande une refonte sur la piste saisonnière S2, calendrier compromis.
- **Discussion générale** : encore à consolider, gel 22/11 tenable si Ch3 gèle le 30/09.
- **Financement 2 mois** : sans démarche ANR cette semaine, on tombe en scénario 2 (dépôt fin novembre) : 5 semaines de manuscrit en moins.
- **Ch1 AoFS soumission** : encore à fixer. Si soumission avant dépôt, ajoute un cycle reviewer en pleine fenêtre d'assemblage.

---

## 9. Ce que Nathan doit demander explicitement au début de la réunion

1. **Accord de principe sur la réinscription en 4A** (D3), formulé de sorte que je puisse le rapporter au CSI comme validé par la Direction.
2. **Feu vert pour la démarche 2 mois ANR MaCCMic** (D4), avec un porteur nommé (Sylvie ou Jérôme) et une échéance cette semaine.
3. **Une position commune sur la remarque Jérôme "S2 saisonnier"** (D2), soit acceptée comme perspective de discussion générale, soit cadrée comme refonte Ch3 avec calendrier révisé.

Sortir avec ces trois choses est le seuil de succès de la réunion.
