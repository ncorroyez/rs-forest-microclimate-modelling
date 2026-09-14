# Chapitre 1 — organisation du projet

Deux sous-dossiers, chacun **auto-documenté** (un fichier détaille méthodes, pré-traitements,
données par figure, description et interprétation) :

## `redaction/`  — rédaction de l'article Chapitre 1
- `manuscript_chap1_EN.md` — le manuscrit prose (IMRaD, article pour soumission) ; **`Manuscrit_chap1.docx`** = version DOCX, figures intégrées (citations en clair `[@clé]` — `--citeproc` bloqué par une clé `.bib` malformée).
- `manuscript_chap1_bullets.md` — **version « masse de bullet points »** du même contenu (relecture rapide / Nextcloud) : tous les chiffres et claims en puces, par section IMRaD.
- **`README_figures_methodes.md`** — fiche détaillée : pour CHAQUE figure (Fig 1–6 + annexes A–E),
  les données d'entrée · pré-traitements · process (script) · description · interprétation, + les
  conventions transversales et la chaîne LiDAR→MuSICA.
- `figures/` — figures du manuscrit (noms narratifs Fig01…FigE2).
- `Bib_nathan_corroyezJune13.bib` — bibliographie.

## `comparaison_versions/`  — comparaison binaire v3.2.0 vs v3.2.3 (yoyo/iter)
- **`README_comparaison.md`** + **`Comparaison_v320_v323_yoyo.docx`** (version DOCX, figures intégrées) — fiche détaillée :
  méthode (bullets) + pour chaque figure (validation HOBO, importance per-unit, heatmap 3 métriques × 2 périodes,
  scatter horaire toutes-T, scatter par-plot ΔTmax & slope tous/chauds) données · process · description · interprétation.
- `figures/` — `FigCmp_*`. `tables/` — `tab_version_*`, `tab_metrics6_*`, `tab_hobo_iter_validation`.

## `config_musica/`  — bundles pour faire tourner MuSICA (Nextcloud)
- **`v3.2.0_legacy/`** (binaire validé `53072278`, sans yoyo) et **`v3.2.3/`** (`8b5139e0`, avec mode yoyo `ABL_flag='iter'`).
- Chacun : `musica` + namelists (`musica.nml`, `musica_soil.nml`, `musica_veg1.nml`, `variables.csv`) + `forcing/` + `README.md`.
- v3.2.3 inclut `forcing/musica_in_Blois_pblh.nc` (h_sbl = PBLH MERRA-2) requis pour le yoyo.

**Conclusion de la comparaison** : le yoyo (couplage ABL) est immatériel (iter ≈ none) et n'existe
que dans v3.2.3, qui valide moins bien (r 0.93→0.63) → **binaire validé v3.2.0 retenu pour l'article** ;
le yoyo documenté comme contrôle de robustesse.

---
*Convention : tout nouveau livrable (analyse ou figure) s'accompagne d'une fiche détaillée
(méthodes · pré-traitements · données · description · interprétation) dans son sous-dossier.*
*Les scripts générateurs restent à la racine du projet (ils écrivent dans `out_files/`) ; les figures/tables
ici sont des copies « livrables » — voir le tableau « figure → script → données » de chaque fiche pour régénérer.*
