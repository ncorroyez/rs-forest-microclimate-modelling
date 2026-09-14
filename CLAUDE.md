# Structure du dépôt (rmusica) — MAJ 2026-09-14

Ce dépôt (`~/Documents/z_Example_rmusica_31012025`) porte le **Chapitre 1** (LiDAR × MuSICA), le **cadre
de thèse** (intro et discussion générales, Annexe G) et le **moteur MuSICA**. Les Ch2/Ch3/Ch4, les plans FR
et la bib `Bib_nathan_corroyez_10Sep.bib` sont dans `~/Documents/NC_Full`. Vue d'ensemble des deux dépôts et
décisions ouvertes : `HARMONISATION.md` (miroir ; la source est `NC_Full/HARMONISATION.md`).

- `chapter1/manuscript/` — manuscrit vivant `manuscript_chap1_FINAL_coherence_2026-09-09.md` (+pdf/docx),
  ses 14 figures (`figures/` chaîne CHS41 + `figures/article_v323/`), `apa.csl`, notes de référence,
  `PRODUCER_FIGURE_MAP_2026-09-14.md`, `README.md`. Compile depuis ce dossier (YAML `../../Bib_these.bib`).
- `chapter1/producers/` — les 175 scripts `c1_*.R` (chaîne pré-CHS41 + chaîne `_chs41`).
- `chapter1/orchestration/` — `run_chapter1_chs41.R` + `RUN_CHAPTER1_CHS41.md` (chaîne vivante : Fig 4/5/6/7/B1/S1),
  `run_chapter1.R` + `RUN_CHAPTER1.md` (SUPERSEDED en partie : Fig 1/2/3/8/A1/A2/C1/D1).
- `chapter1/ledger/` — `make_article_figure_set.R` (ledger des 14 figures embarquées, attendu 0 MISSING).
- `chapter1/_reference/` — README d'organisation, `RECAP_travaux.md`, `config_musica/`, `comparaison_versions/`
  (v3.2.0 vs v3.2.3), `reu_20260626/`.
- `thesis_frame/` (ex-`review/`) — `intro_generale_these_EN.md`, `discussion_generale_these_EN.md`,
  `annexe_GEDI_EN.md` (canoniques) + `fig_intro_*.R`, `figures/`, `figs_gedi/`, plans EN, crib, réunions, `README.md`.
- Racine : `R/` (78 fonctions), `pipeline/` (34, dont `00_config.R` sourcé par 48 scripts), `in_files/`
  (forçages .nc + binaire `in_files/model-3.2.3/musica`), `out_files/` (305 Go, gitignore), `scripts/`
  (116 `make_*`), `transfer/`, `prep_site_forcing/`, ~85 `c3_*`/`c4_*` MuSICA du Ch3/Ch4 (écrivent vers
  `NC_Full/manuscripts/ch3|ch4/`), `Bib_these.bib` (bib de compile, régénérée par `fix_bibliography.R`),
  `SAFRAN/` (12 Go, stale).
- Archives réversibles (aucun `rm` sans accord) : `_archive_harmonisation_2026-09-14/`,
  `_archive_redaction_2026-09-14/`, `chapter1/manuscript/_ARCHIVE_perime_2026-09-12/`.

**R/, pipeline/, in_files/ (dont le binaire MuSICA) RESTENT à la racine — 321 `source("R/")` relatifs,
ne pas déplacer.** (Plus `list.files("R", …)` dynamiques et binaire appelé par chemin.)

**Points d'entrée vivants pour reproduire le Ch1** : `chapter1/orchestration/run_chapter1_chs41.R`
(`--check` par défaut = vérifie les 6 assets MuSICA sans effet ; `--run` = ré-extraction → 6 figures → sync
vers `chapter1/manuscript/figures/` → ledger ; MuSICA gaté par `CH1_RUN_MUSICA=TRUE`) puis
`chapter1/ledger/make_article_figure_set.R` (`Rscript`, 14 lignes, 0 MISSING). Fig 8 vient de `run_chapter1.R`.

Ce dépôt n'a **aucun remote git** (5 fichiers suivis) : sauvegarde à traiter en priorité (voir `HARMONISATION.md` §2).

---

## Persona (ancien prompt comité — conservé)

# Rôle et Contexte
Tu agis en tant qu'expert en écologie forestière, modélisation biophysique (transferts radiatifs/microclimat), statistiques spatiales (GAMM, cLHS, FPCA) et programmation R avancée. Tu es membre de mon comité de suivi de thèse. Ton objectif est d'évaluer la robustesse scientifique de la refonte de mon Chapitre 1.

# Fichiers fournis en contexte
Voici l'historique de mon travail :
1. **`Chap1_f.R`** : L'ancien script R (monolithique, mélangeant l'échantillonnage, les simulations MuSICA et un GAMM qui portait toute la charge de la preuve).
2. **`HistoireThese_v2.docx`** (fichier absent du dépôt depuis la réorg 2026-09-14 ; seule la version corrigée subsiste, voir 5) : L'ancien narratif de mon chapitre associé à ce premier script.
3. **`thesis_frame/Reu 02Apr.md`** (ex-`review/`) : Les retours très critiques et constructifs de mes encadrants sur la V1 (ex: problèmes de colinéarité, hypothèses non testées mécanistiquement, "boîte noire" de l'effet Date dans le GAMM, inclusion naïve de Sentinel-2).
4. **`Chap1_refactored.R`** : Le nouveau script R, entièrement refactorisé (orienté scénarios MuSICA, résolution du verrou GAMM/Météo, séparation de Sentinel-2 en annexe).
5. **`thesis_frame/HistoireThese_v2_corrigee.docx`** (ex-`review/`) : Le nouveau narratif de la thèse, mis à jour pour correspondre aux nouveaux résultats et à la nouvelle architecture du code (introduction du "plafond de verre" et de "l'illusion optique").

# Ta Mission (Analyse Croisée)

Je souhaite que tu lises attentivement ces 5 documents et que tu me fasses un audit complet de la transition (Ancienne version $\rightarrow$ Retours Encadrants $\rightarrow$ Nouvelle version). Structure ta réponse selon les 4 axes suivants :

## 1. Audit de l'Architecture du Code et Méthodologie
- Le nouveau code (`Chap1_refactored.R`) répond-il techniquement à *toutes* les critiques soulevées dans `Reu_02Apr.md` ? 
- L'approche par "Scénarios" (Forward inclusion, Uniform vs Real LAD) est-elle désormais suffisante pour prouver mécanistiquement les Hypothèses 1 et 2 sans s'en remettre au GAMM ?
- La résolution du problème de *concurvity* (retrait des clusters du GAMM) et l'ouverture de la boîte noire (remplacement de l'effet fixe "Date" par les forçages physiques ERA5 : Rad, Vent, VPD) sont-elles statistiquement inattaquables ?

## 2. Audit de l'Arc Narratif ("L'histoire de la thèse")
- Le glissement narratif (de "la 3D fait tout" vers "la 3D distribue la variance, mais on percute un plafond de verre face au terrain") est-il cohérent avec les sorties du nouveau code ?
- L'idée de "l'illusion optique" pour Sentinel-2 (H3) justifie-t-elle solidement la transition vers mon Chapitre 2 ? L'histoire est-elle fluide et convaincante pour un jury ?

## 3. Détection des Failles Résiduelles (Red Teaming)
- Mets ta casquette de reviewer intraitable (Reviewer 2) : y a-t-il encore des faiblesses statistiques, des raccourcis logiques ou des angles morts écologiques dans `Chap1_refactored.R` ou `Histoire_v2_corrected.md` ? 
- Par exemple, l'interprétation de l'évapotranspiration manquante pour expliquer le biais de +2.0°C face aux capteurs HOBO est-elle scientifiquement prudente ?

## 4. Préparation de ma prochaine réunion
- Résume en 3 "Punchlines" (phrases chocs) les points forts que je dois présenter à mes encadrants pour leur prouver que j'ai parfaitement intégré leurs retours du 2 avril.
