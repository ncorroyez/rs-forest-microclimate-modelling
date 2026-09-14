# Plan de restructuration Ch1 (2026-09-08)

Décisions Nathan : **B1 (test H2) → figure main ; G1 (incertitude) → Table 2 ;
B2 (gradients verticaux) → figure main (fort)**. Fichier cible = la version de prose
courante `manuscript_chap1_EN_native20_csi3pass_stylepass_2026-09-07.md`.
À appliquer APRÈS que tous les chiffres/figures soient sur CHS41-Rmerge/sans-vent.

## 1. Déplacements

| Bloc actuel | Devient | Contenu |
|---|---|---|
| Fig **B1** (controlled test H2, beta shapes) | **Figure 6** (main) | test H2 contrôlé, LAI/Hmax fixés |
| Fig **B2** (within-canopy T/wind/RH/VPD) | **Figure 7** (main) | gradients verticaux réel vs uniforme |
| Table **G1** (incertitude ΔTmax par archétype) | **Table 2** (main) | médiane + IC des effets |
| Fig **B3** (per-variable within-canopy gradient) | reste annexe → **Figure B1** | inchangé de rôle |
| Table **G2** (incertitude pente) | reste annexe → **Table G1** | inchangé de rôle |

## 2. Nouvelle numérotation

- **Figures main** : 1 typologie · 2 pipeline · 3 PCA · 4 attribution · 5 point de
  fonctionnement · **6 test H2 contrôlé (ex-B1)** · **7 gradients verticaux (ex-B2)**.
- **Tables main** : 1 moyennes par archétype · **2 incertitude des effets ΔTmax (ex-G1)**.
- **Annexe B** : ne garde que l'ex-B3 → renommée **Figure B1** ; titre à réduire
  (« The within-canopy gradient, per variable » ; retirer « controlled balance test »
  et « within-canopy gradients » qui partent en main).
- **Annexe G** : ne garde que l'ex-G2 → **Table G1** ; titre inchangé.
- **Option (revue Fable)** : descendre **Fig 3 (PCA) en annexe** pour tenir le main à 6
  figures. NON décidé — à confirmer avec Nathan. Par défaut on garde 7 figures main.

## 3. Placement narratif

- **Figure 6 (H2)** : dans la section qui discute H2 (profil), là où le texte renvoie
  aujourd'hui à « Appendix B » pour le test contrôlé (§3.2 / §4.1).
- **Figure 7 (gradients)** : comme figure mécanistique — pourquoi l'effet profil est petit
  au readout 1 m (les profils divergent en hauteur, **convergent à 1 m**). §4 (mécanisme).

## 4. Cross-références à mettre à jour (checklist, lignes du fichier stylepass)

- `Fig. B1` / `Figure B1` (test H2) → **Figure 6** : L245, L367, L371.
- `Fig. B2` / `Figure B2` (gradients) → **Figure 7** : L168, L259, L369, L373, L375.
- `Figure B3` (per-variable) → **Figure B1** (annexe renommée) : L369, L375.
- `Table G1` (incertitude ΔTmax) → **Table 2** : L214, L392, L458, L480, L512.
- `Table G2` (incertitude pente) → **Table G1** : L216, L493, L511.
- `Appendix B` : L168, L174, L212, L245, L259, L296, L365 → repointer vers Figure 6/7
  ou l'annexe B slim selon le contexte.
- `Appendix G` : L214, L216, L288, L455, L458 → repointer (Table 2 vs annexe G slim).

## 5. Statut régénération CHS41-Rmerge/sans-vent (préalable aux déplacements)

| Figure/Table | Source CHS41 | Statut |
|---|---|---|
| Fig 4 attribution + Table 2 (ex-G1) + Table G1 (ex-G2) | `perturb_chs41_nowind.csv` | **rerun en cours** |
| Fig 6 (ex-B1, test H2) | `scripts/c1_h2_controlled_chs41.R` | **prêt, à lancer** (42 runs) |
| Fig 7 (ex-B2, gradients) | sims archétype réel/uniforme sur CHS41 | **à créer** (B2 lit du legacy ; `make_vertical_profiles_wind_rh_vpd_metres.R` à re-pointer) |
| Fig B1 (ex-B3, per-variable gradient) | idem sims archétype CHS41 | **à créer** avec Fig 7 |

## 6. Ordre d'exécution

1. Fin du rerun Fig 4/H1 (en cours) → régénérer Fig 4, Table 2, Table G1.
2. Lancer B1 CHS41 (`c1_h2_controlled_chs41.R`) → Figure 6.
3. Créer les sims archétype CHS41 (réel + uniforme, 4 archétypes) → Figures 7 et B1.
4. Appliquer les déplacements + renumérotation + checklist §4 sur le fichier stylepass.
5. MAJ Table H1 (CHS41-Rmerge, h_sbl 550 m, sans vent), §2.4 (vent → annexe), Bouwen.
