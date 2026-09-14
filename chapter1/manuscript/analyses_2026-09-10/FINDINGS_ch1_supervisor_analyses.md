# Ch1 — supervisor analyses, findings (2026-09-10)

Six analyses requested by the encadrants. All on **CHS41-Rmerge / no wind correction**
(provenance verified). Figures in `out_files/Chapter1/figures/diag_2026-09-10/`. Scripts: `c1_diag_*.R`.

**INTEGRATED INTO THE MANUSCRIPT 2026-09-11 (Nathan's approved minimal set, voice pass):**
- **#5/#6** → **main Figure 8** "Sensitivity to the analysis footprint" + a paragraph in the
  Discussion §4.3 (`figures/Fig8_footprint_radius.png`). Promoted from appendix to main on Nathan's
  call: footprint robustness is part of the model's sensitivity analysis to structure variables.
- **#1** → one sentence in §2.5 (profiles converge near the top, so the 1 m readout sits where
  the canopy signal is largest).
- **#3 soil** → one sentence in §4.3 (the 1 m air offset understates the ground signal:
  soil surface spans 7.4 °C vs 0.5 °C at 1 m).
Manuscript now 41 pages, 62 keys, figures 1→7, 0 unresolved, 0 em-dash.
**HELD OUT (diagnostics only):** #2 (normalised, sign artefact) and #4 (add-block, different
experiment type). #3 leaf temperature = stated limitation.

Five of the six were **pure post-processing of data already on disk**; nothing heavy was
re-run. The sixth (#3 leaf temperature, #4 below-sensor add) is bounded below.

---

## #1 — Vertical T(z) by LAI: do the profiles converge at 1 m?  → NO (they converge at the top)
`c1_diag_vertical_by_lai.R` · `diag_vertical_Tprofile_by_LAI.png`
Gaussian LAD, peak μ = 0.5 Hmax, spread σ = 0.24, cover 0.87, Hmax 25 m, LAI 2→6.
- Air temperature at the **1 m readout** falls from 22.82 °C (LAI 2) to 22.10 °C (LAI 6):
  spread **0.71 °C** across LAI.
- Spread at the **true canopy top (z ≈ Hmax = 25 m)**: **0.031 °C** (and it stays ~0.03 °C up
  through the free air aloft; the air column extends to ~43 m).
- **The profiles are pinned together at the canopy top and fan out downward.** The 1 m
  sensor sits near the maximum LAI separation, not at a convergence point. This **supports
  the 1 m readout** as a place where the LAI signal is largest, and answers the encadrants'
  question directly: they rejoin at the top, not at 1 m.

## #2 — Normalising the levers by the plot's own buffering  → double-edged, touches a locked framing
`c1_diag_normalized_levers.R` · `diag_levers_normalised_by_buffering.png`
Each native-step effect / plot's own |ΔTmax|. Median fraction (all four archetypes; 29 plots
within |base| ≤ 0.05 °C of the amplify→buffer sign change dropped: 0 in P1, 2 P2, 14 P3, 13 P4):

| lever | P1 | P2 | P3 | P4 |
|---|---|---|---|---|
| Leaf area | +0.147 | −0.265 | −0.221 | −0.596 |
| Cover | −0.080 | −0.339 | −0.379 | −0.771 |
| Height | −0.001 | 0.000 | −0.008 | 0.002 |
| **Profile** | −0.089 | **−0.264** | **−0.449** | **+0.341** |

- **Sign structure matters:** the baseline ΔTmax *amplifies* (+) in P1–P3 (median +0.185, +0.194,
  +0.144) and *buffers* (−) in P4 (−0.295). So the fraction changes meaning across the transition;
  read within an archetype, not across it. (This is why P1's leaf-area fraction is positive: a
  warming step on an amplifying base.)
- The normalisation does "show more": the profile is a **larger fraction of local ΔTmax** than the
  absolute numbers suggest, up to **0.45 in P3**.
- **DECISION FLAG:** in P3 the profile fraction (0.45) **exceeds the cover fraction (0.38)** — a
  reordering the absolute view does not show. It is a fraction-of-a-small-baseline effect (P3 base
  is near the sign change) and does **not** change the absolute second-order conclusion (the
  profile's absolute effect is still ≤ the cover step). It runs against the locked "profile ≤ cover
  step" wording if read as a ranking. **Diagnostic only, not for the manuscript.**

## #3 — Full profile (canopy top, 1 m air, soil); leaf temperature  → soil is the big signal; leaf T not emitted
`c1_diag_full_profile.R` · `diag_full_profile_air_soil.png`, daytime JJAS mean per archetype (°C):

| archetype | canopy top (z≈Hmax) | free air aloft | 1 m air | soil surface |
|---|---|---|---|---|
| P1 (open) | 22.48 | 22.25 | 22.73 | **27.17** |
| P2 | 22.39 | 22.30 | 22.95 | 22.46 |
| P3 | 22.49 | 22.29 | 22.85 | 21.83 |
| P4 (dense) | 22.40 | 22.31 | 22.44 | **19.77** |

- Air is pinned near forcing at the canopy top (~22.4 °C everywhere, spread ~0.10 °C). The
  **soil surface** swings **7.4 °C**
  from open (27.2) to dense (19.8), far larger than the 1 m air spread. **The canopy-density
  signal is strongest at the ground, muted at 1 m.** This is the substance for the encadrants'
  "dommage de s'arrêter à la validation HOBO → sinon être prudent et en discuter": the 1 m air
  metric under-represents the surface response.
- **LIMITATION (not an analysis):** the v3.2.0 output has **no leaf/foliage temperature**
  variable (`Tair_z`, `T_soil`, `wind_z`, `wair_z` only; `nveg=10` but no `Tleaf`). Emitting it
  would require recompiling MuSICA. Report as a limitation, or schedule a model change (out of
  scope for a re-analysis).

## #4 — LAI peak position vs the 1 m sensor  → peak effect grows with density; "below sensor" not in the design
`h2_gaussian_grid_chs41.csv` (μ-sweep already on disk), σ = 0.24:

| LAI | μ-effect (°C) | best μ |
|---|---|---|
| 2 | 0.046 | 0.3 |
| 3 | 0.078 | 0.7 |
| 4 | 0.124 | 0.7 |
| 5 | 0.166 | 0.6 |
| 6 | 0.247 | 0.6 |

- Raising the peak matters more in denser canopy, with the optimum at mid-upper (μ ≈ 0.6),
  reproducing "mid-canopy overtakes top-heavy". This is already Fig 6.

### #4b — Add LAI at the sensor vs in the canopy (RAN, `c1_diag_addblock_dryrun.R`, `diag_addblock_position_by_LAI.png`, `h2_addblock_chs41.csv`)
Literal "below 1 m" is sub-grid: MuSICA re-bins LAD onto 10 layers (2.5 m each) and the 1 m
sensor sits in the lowest layer (air-grid centre 0.71 m). So the resolvable form: add +0.5 LAI
either LOW (0–2 m, at the sensor) or HIGH (~12.5 m, canopy peak) on a base Gaussian, compare
ΔTmax. Result (°C):

| base LAI | add-high | add-low | position (low−high) |
|---|---|---|---|
| 2 (sparse) | +0.144 | +0.215 | +0.071 |
| 4 (median) | −0.096 | −0.014 | +0.082 |
| 6 (dense edge) | −0.678 | −0.885 | **−0.207** |

- **The better placement reverses with density.** In sparse/median canopies, foliage added
  HIGH buffers more (understory foliage sits in the warm near-ground layer and worsens the
  daytime maximum). At the dense edge (LAI 6), it flips: understory foliage shades the sensor
  directly and buffers more, because the upper canopy already intercepts the radiation.
- **Caveat:** LAI 6 is the dense edge of the design (median 3.48, P4 ≈ LAI 5), the same edge
  flagged for the σ jump. The design-body result (LAI 2–4) is "upper placement buffers more".
- This extends H2 (top-heavy buffers more) to the "add" framing and reaches the near-sensor
  region the μ-sweep could not (μ min = 0.3 = 7.5 m). Diagnostic; not in the manuscript.

## #5 / #6 — Footprint radius, full sweep, everything as a function of radius  → 20 m robust; only Hmax moves, and it does not matter
`c1_diag_radius_dtmax.R` (ΔTmax only) and `c1_diag_radius_full.R` (all quantities) ·
`diag_radius_dtmax_clhs.png`, `diag_radius_full_vs_radius.png`.
**All 7 radii on disk are done: 5, 10, 12.5, 15, 20, 25, 50 m** (~400 cLHS plots each; 20 m = the
reference, so it is the 0 line, not a missing radius). The pipeline runs on the **400 cLHS** plots
(not the 53 HOBO) as the **paired difference** vs the 20 m reference (not a correlation).

Everything as a shift from the 20 m reference (pooled, 5 m → 50 m):

| quantity | at 5 m | at 50 m |
|---|---|---|
| ΔTmax (°C) | −0.028 | −0.009 |
| buffering slope (–) | −0.004 | +0.004 |
| one-sided LAI | −0.001 | +0.041 |
| fractional cover | 0.000 | 0.000 |
| **Hmax (m)** | **−2.10** | **+1.85** |
| LAD centroid (rel.) | +0.056 | −0.039 |

- Both **MuSICA responses barely move** with footprint (ΔTmax ≤ 0.028 °C pooled, ≤ 0.091 °C in the
  dense P4 at 5 m; slope ≤ 0.004). At 10–15 m the ΔTmax deviation is ≤ 0.02 °C for most archetypes.
- Of the four structural inputs, **only Hmax moves with footprint** (± ~2 m), because it is an
  extremal statistic: the tallest tree caught grows with the area sampled. LAI, cover and the LAD
  centroid are near-invariant.
- **Internal consistency:** the one variable the footprint perturbs (Hmax) is exactly the one that
  is thermally negligible (|effect| ≤ 0.001 °C m⁻¹), which is **why ΔTmax is robust to footprint**.
- **The 20 m footprint is robust**; footprint sensitivity is a second-order effect concentrated in
  dense stands.

---

## Decisions taken (Nathan: "je te laisse voir")
1. **#2 normalised view** stays a **diagnostic only, not in the manuscript.** It reorders profile
   above cover in P3 by a fraction-of-a-small-baseline effect and would run against the locked
   "profile ≤ cover step, second-order" wording. The absolute conclusion is unchanged and stands.
2. **#4 add LAI at the sensor vs canopy: RAN** (9 sims, resolvable form; literal below-1 m is
   sub-grid). The position effect reverses with density (see #4b). Diagnostic; not in the manuscript.

## One limitation to state
**#3 leaf temperature** is not emitted by MuSICA v3.2.0; reporting it needs a model recompile.
Air and soil are covered; the soil-surface response is the larger, under-reported signal.
