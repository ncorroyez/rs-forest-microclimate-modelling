# Ch1 — manuscrit canonique

- **LE fichier vivant : `manuscript_chap1_FINAL_coherence_2026-09-09.md`** (+ `_2026-09-11.pdf/.docx`).
  `manuscript_chap1_EN.md` (dans `thesis_frame/`, en .docx) est l'ancienne lignée (Shapley) — ne pas y toucher.
- Compile **depuis ce dossier** (le YAML déclare `bibliography: ../../Bib_these.bib` = racine rmusica ;
  `apa.csl` est ici ; pas de `--bibliography` à ajouter, sinon voir mémoire `ch1_compile_bib_trap`) :
  `sed 's/≫/>>/g' manuscript_chap1_FINAL_coherence_2026-09-09.md | pandoc -f markdown --citeproc --csl apa.csl --resource-path=. -o manuscript_chap1_FINAL_coherence_$(date +%F).pdf --pdf-engine=xelatex -V mainfont="DejaVu Serif" -V monofont="DejaVu Sans Mono" -V geometry:margin=2.4cm -V fontsize=11pt -V linkcolor=blue`
- Les 14 figures embarquées : `figures/` (7, chaîne CHS41 : Fig 4/5/6/7/B1/S1 + Fig 8) et `figures/article_v323/` (7).
  Elles sont ÉCRITES par `../orchestration/run_chapter1_chs41.R` (sync) et `../ledger/make_article_figure_set.R`
  (14 lignes, 0 MISSING attendu) ; producteur de chacune : `PRODUCER_FIGURE_MAP_2026-09-14.md`.
- Tout le reste ici = support / historique / sauvegardes : notes de référence (`CHS41_numbers_ref_*`, `FINAL_markers_map_*`,
  `CITATION_AUDIT_*`, `SUPERVISOR_COMMENTS_*`, `README_figures_methodes.md`, `note_wind_correction_FR.*`), plans et
  réponses aux encadrants, `Bib_nathan_corroyezJune13.bib` (copie PÉRIMÉE, ne pas compiler dessus), `tools/`, `tmp/`,
  `csi3pass/`, `analyses_2026-09-10/`, `pour_jerome_20260820/`, et `_ARCHIVE_perime_2026-09-12/` (201 artefacts, `rm` en attente).
- Vue d'ensemble des deux dépôts : `../../HARMONISATION.md` (miroir ; source dans NC_Full).
