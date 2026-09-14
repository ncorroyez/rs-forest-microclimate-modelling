# Chiffres CHS41-Rmerge / sans-vent pour la MAJ du manuscrit Ch1 (2026-09-08)

Source : `out_files/Chapter1/tables/{perturb_chs41_nowind,h1_realvsreal_chs41,h2_controlled_chs41}.csv`.
Convention : ΔTmax time-matched, sans shift ; fenêtre juin–sep ; T à 1 m.

## Table 2 (ex-G1) — ΔTmax par pas natif, par archétype, moyenne [IC 95 %]

| Pas | P1 | P2 | P3 | P4 |
|---|---|---|---|---|
| LAI +0.5 | +0.024 (+0.021, +0.026) | −0.057 (−0.067, −0.047) | −0.056 (−0.066, −0.046) | −0.203 (−0.215, −0.191) |
| LAI −0.5 | +0.020 (+0.017, +0.023) | −0.033 (−0.041, −0.026) | −0.040 (−0.049, −0.032) | −0.174 (−0.187, −0.161) |
| fCover +10 | −0.017 (−0.020, −0.014) | −0.077 (−0.087, −0.066) | −0.084 (−0.095, −0.074) | −0.261 (−0.277, −0.245) |
| fCover −10 | −0.026 (−0.031, −0.021) | −0.059 (−0.068, −0.051) | −0.072 (−0.083, −0.062) | −0.235 (−0.251, −0.218) |
| Hmax +1 m | −0.001 (−0.003, 0.000) | 0.000 (−0.001, +0.001) | −0.001 (−0.002, +0.001) | −0.001 (−0.003, +0.001) |
| Hmax −1 m | 0.000 (−0.002, +0.002) | +0.001 (0.000, +0.003) | 0.000 (−0.002, +0.002) | −0.001 (−0.002, 0.000) |
| profil (réel−unif) | −0.033 (−0.041, −0.026) | −0.055 (−0.066, −0.044) | −0.043 (−0.059, −0.025) | +0.136 (+0.101, +0.178) |

## Leviers Fig 4 (par unité, médiane par archétype)
- Quantité (LAI/unité) : P1 +0.054 · P2 −0.105 · P3 −0.088 · **P4 −0.410**
- Couvert (par 10 pts) : P1 −0.146 · P2 −0.691 · P3 −0.675 · **P4 −2.650**
- Hauteur : négligeable partout (≤ 0.001 °C/m)
- Profil (réel−unif) : P1 −0.016 · P2 −0.052 · P3 −0.065 · P4 +0.104

## §3.2 texte — substitutions clés (ancien FR-Blo+vent → nouveau CHS41/sans-vent)
- P4 « adding 0.5 of leaf area cools by **0.26 °C** » → **0.203 °C** ; cover 0.27 → **0.261 °C**.
- P1 leaf-area (warming) 0.05 → **0.024 °C** ; asymétrie P4 (up −0.203 vs down −0.174).
- Hauteur « |effect| ≤ 0.004 °C/m » → **≤ 0.001 °C/m**.

## Fig 6 (test H2 contrôlé)
- Effet de forme du profil : **0.10 à 0.48 °C** across LAI (ancien 0.15–0.51). Amplifie à bas LAI,
  tamponne à haut LAI ; raising centroid cools à LAI ≥ 3.

## H1 réel-contre-réel (CHS41)
- |réel − forme moyenne| : P1 0.011 · P2 0.027 · P3 0.055 · P4 0.071 °C (tous < leviers quantité).
- |réel − uniforme| : P1 0.017 · P2 0.052 · P3 0.073 · P4 0.106.

## Forçage (Table H1 + §2.4 + Appendix H)
- Fichier : **`MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc`** (plus FR-Blo).
- Sources : Tair/RH/précip = CHS 41 ; radiation/vent/pression = ERA5 ; **h_sbl = MERRA-2, médiane 550 m**.
- **PAS de correction de vent** en baseline (le vent brut ERA5 est utilisé) ; la correction de vent
  log-profile devient une **annexe de sensibilité**.
- Retirer la claim « insensible au h_sbl » : vérifié que le h_sbl bas de FR-Blo (48.8 m) **aplatissait
  le levier LAI en P2** (−0.105 → −0.146) ; dense P3/P4 inchangé.

## Bouwen
- Accord sur : ordre quantité→couvert, hauteur faible predictor. Mais **désaccord sur la saturation LAI** :
  le levier LAI est **maximal en P4 dense (−0.41/unité, LAI ~5)** alors que Bouwen trouve le LAI insensible
  au-dessus de 4–6 m²/m². À formuler comme accord partiel + différence explicite, pas accord complet.
