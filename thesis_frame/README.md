# thesis_frame — cadre général de la thèse (intro, discussion, Annexe G)

- **Canoniques** : `intro_generale_these_EN.md` (intro générale, 16 figures dans `figures/`),
  `discussion_generale_these_EN.md` (discussion générale), `annexe_GEDI_EN.md` (Annexe G, 4 figures dans `figs_gedi/`).
  Chacun a son `.pdf` à côté. Les plans FR correspondants sont dans `NC_Full/plans/PLAN_{intro,discussion}_generale_FR.md`.
- Compile **depuis ce dossier**. L'intro porte `bibliography: ../Bib_these.bib` dans son YAML ; la discussion et
  l'Annexe G n'ont pas de YAML → ajouter `--bibliography ../Bib_these.bib`. CSL = celui du skill `these-manuscript` :
  `sed 's/≫/>>/g' <fichier>.md | pandoc -f markdown --citeproc --csl ~/.claude/skills/these-manuscript/apa.csl --bibliography ../Bib_these.bib --resource-path=. -o <fichier>.pdf --pdf-engine=xelatex -V mainfont="DejaVu Serif" -V monofont="DejaVu Sans Mono" -V geometry:margin=2.4cm -V fontsize=11pt -V linkcolor=blue`
  (`--resource-path` obligatoire : le `.md` passe par stdin, une figure manquante disparaît sans erreur).
- Figures de l'intro : produites par `fig_intro_*.R` + `_intro_fig_style.R` (écrivent `thesis_frame/figures/`) ;
  autorisations/emprunts dans `PERMISSIONS_figures_intro.md`, numérotation dans `PLAN_figures_intro_2026-09-09.md`.
- Matériel de cadre général conservé ici : plans EN (`plan_intro_generale_these.md`, `plan_article_chap1*.md`),
  `defense_crib_chap1.*`, comptes rendus de réunions (`meeting_*`, `Reu*`), `disc_crosscheck_*`, audits pipeline,
  `HistoireThese_v2_corrigee.docx`, `planning_these_juin_sept_2026.md`.
- Tout le reste = support / historique / sauvegardes (`*_rewrite_*`, `*_voicepass_*`, `.bak_*` du jour, `*.bak.md`,
  `manuscript_chap1_EN.docx` = ancienne lignée Ch1). Ex-dossier `review/` (renommé le 2026-09-14) ; certains `.md`
  anciens citent encore `review/` ou `manuscripts/ch3/`. Vue d'ensemble : `../HARMONISATION.md` (miroir de NC_Full).
