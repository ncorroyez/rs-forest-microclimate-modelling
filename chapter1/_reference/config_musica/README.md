# Configurations MuSICA — bundles autonomes (Nextcloud)

Deux dossiers self-contained pour faire tourner chaque binaire (binaire + namelists + `variables.csv` + forçage + README).

- **`v3.2.0_legacy/`** — binaire **validé** (md5 `53072278`, r=0.93), **sans yoyo**. Mode du chapitre.
- **`v3.2.3/`** — binaire `8b5139e0`, valide moins bien (r≈0.63) ; contient le **mode yoyo `ABL_flag='iter'`** (2 pré-requis : phénologie 2020/366 + `h_sbl`=PBLH ; forçage `musica_in_Blois_pblh.nc` fourni).

Chaque dossier a son `README.md` détaillé (contenu, lancement, conventions, pièges). Les .nml sont des **gabarits** : dans l'article ils sont réécrits par plot via le wrapper R `rmusica`/`musica.tools`.

- **`scripts/`** — tous les scripts R/bash de run & comparaison v3.2.0 vs v3.2.3 (+ `R_helpers/`, `Chapter3_config.R`, `pipeline_00_config.R`) + **`README_reproduction.md`** (dépendances, mapping script→figure, caveat chemins).
- **`inputs/`** — données d'entrée : `clhs_sample_floor05_v2.rds` (400 plots), `Blois_lad_z05_r25.csv` (LAD 0,5 m), `Blois_data_temperature.csv` (HOBO), `data_Blois_utm31n.geojson` (positions).
