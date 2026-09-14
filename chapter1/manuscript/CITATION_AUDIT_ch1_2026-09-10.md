# Chapter 1 — full-text citation audit (2026-09-10)

Every substantive citation in `manuscript_chap1_FINAL_coherence_2026-09-09.md` was
verified against the cited PDF (not the abstract), on three axes: **fit** (does the
full article support this exact sentence), **direction** (does the source's sign and
conditioning of effect match what the sentence implies), **precision** (does the
sentence name the paper's finding, or only its topic). 43 substantive keys read across
4 agents; 18 routine named-object attributions (MuSICA, ERA5, MERRA-2, FPCA, cLHS, SA
methods, silhouette, extinction theory, ABL, VCI, foliage-profile method, ALS retrieval)
verified at entry level per the skill. No citation was recommended for removal.

## Verdict in one line

The citations are sound. Signs and conditioning are correct almost everywhere; the
load-bearing numbers are exact (de Frenne 4.1 ± 0.5 °C / 98 sites / five continents;
Gril 72.6 and 9.7 pts m⁻², 4.5 cm, 17 June 2021, *R*² = 0.91; Laslier 0.74 vs 0.92 at
Blois). Five items need a decision; the rest hold.

---

## Tier A — scientific decisions (do NOT edit silently)

### A1. "the vertical arrangement of foliage is often ranked first" [@parker; @ehbrecht] (L72–73) — BOTH cites fail to support it
- **Parker (2004)**: descriptive 3D old-growth Douglas-fir/hemlock study, foliage
  bottom-heavy, radiation absorbed low, shallow gradients. Does **no** ranking.
- **Ehbrecht (2019)**: ranks canopy **openness first**; vertical stratification (ENL)
  ranks **third** among heterogeneity metrics, and the heterogeneity effect is small.
  This **contradicts** "ranked first" and in fact **agrees with this chapter's own
  second-order conclusion**.
- This sentence sets up H2, so it is load-bearing.
- **Options:** (a) reframe to the defensible motivation, that the vertical position of
  foliage sets where radiation is intercepted and is *expected* to matter mechanistically
  (Parker documents exactly this), and move Ehbrecht to the Discussion as agreeing that
  the profile is second-order (a "but see" that strengthens §4.1); (b) keep "ranked first"
  and find sources that actually rank it first. Recommend (a).

### A2. maclean 2017 mis-cited on the concurvity clause (L81–82) — MISMATCH
- Sentence: "Their coefficients are unstable under concurvity, and individual
  contributions are not identifiable [@greiser; @maclean2017]." Maclean 2017 never
  discusses concurvity or identifiability; that critique is Greiser's. Maclean's adjacent
  point is that phenomenological models predict poorly outside their calibration range.
- **Clean fix (both cites kept, re-apportioned):** "Their coefficients are unstable under
  concurvity, so individual contributions are not identifiable [@greiser]. Phenomenological
  microclimate models also predict poorly outside their calibration range [@maclean2017]."
  Recommend applying.

### A3. park 2026 is a diversity experiment on VPD, not an observational temperature-offset gradient (L509)
- Sentence: "the canopy-density control on the offset reported from observational
  gradients [@kovacs; @park]." Kovács is a clean observational gradient on the temperature
  offset. Park is a **manipulative** tree-diversity experiment (FAB2) whose buffered
  variable is **daily VPD amplitude**, with density only a mediator of a diversity effect.
- **Fix:** let Kovács carry "observational gradients"; recharacterise or relocate Park
  (e.g. as diversity raising cover and damping VPD, alongside Schnabel at L449 which is
  the same kind of result). Recommend narrowing the sentence to Kovács and moving Park
  next to Schnabel.

### A4. Bouwen sentence (a) conflates two thesis chapters (L67–68)
- "A LiDAR-MuSICA coupling identified a canopy-density threshold across stand types
  [@bouwen]." The cross-stand-type threshold (PAIe ≈ 2.5) is a **MuSICA + measured-PAI**
  result; the **LiDAR**-MuSICA coupling is a *separate* maritime-pine vertical-complexity
  chapter. Direction is right, attribution is not.
- **Fix:** "MuSICA simulations across stand types identified a canopy-density threshold on
  buffering [@bouwen]." Recommend applying.

### A5. "CHS 41" attributed near a Gril cite (macroclimate-reference sentence) — minor
- The open-field grassland station at 1.5 m is confirmed in Gril 2023; the label "CHS 41"
  is your own MuSICA station name, absent from the paper. As written it is defensible if
  read as your label, but it should not be read as Gril's naming. Optional tidy.

---

## Bouwen thesis vs standalone article (your flag)

The standalone article `bouwenMicroclimateVariationsRoughness2026` (Agric. For. Meteorol.
389:111442) is a purely observational ICOS roughness-sublayer turbulence study: full-text
grep finds no MuSICA, no yoyo, no ERA5/Seidel, no h_ABL/blending height, no beta profile,
no ΔTmax. **Keep the thesis as the source for all five Bouwen sentences.** The article is
not a substitute for any of them. Its only added value is an **optional companion cite on
the boundary-layer-coupling sentence (L209)**, as the published, evaluated evidence that
roughness-sublayer flux corrections above forests are needed and work (the physical
foundation the yoyo rests on). Your call whether to add it.

Note on (e), h_sbl: what Bouwen took from ERA5 is h_ABL, then set blending height
z_b = 0.1·h_ABL. Keep the factor-of-ten explicit if you touch that sentence.

---

## Tier B — precision (topic vs finding): mostly DECLINE, by design

Agents flagged several intro/methods sentences as "topic-level" (baldocchi, lenoir,
kemppinen, meeussen, depauw, hes, kolstela, aklilutesfaye, greiser, lembrechts,
zellweger-seasonal, bonan). On review, most sit in the **Introduction's motivation list**
("microclimate controls communities, regeneration, microrefugia…") or in **Methods
mechanism statements** ("radiation attenuates following extinction theory"). In those
contexts a topical citation is correct; expanding each into a full finding sentence would
bloat the passage and break the intro's pace, against the voice budget. This matches the
skill's own nuance: routine and topical attributions are fine at entry level.

Two Tier-B items are worth considering because the sentence makes a specific claim:
- **zellweger-seasonal (L61)**: its buffering result is specifically **summer maximum**,
  driven by **cover** (not height); "taller" is only weakly supported. If you want, drop
  "taller" from "Denser, taller and more closed canopies", since diaz-calafat also does
  not test height.
- **aklilutesfaye (L629)**: could name that saturation is worst for NDVI and least for the
  red-edge TVI, but the current topical use is adequate for a forward-looking sentence.

Minor non-blocking caveats (no action): atkins predicts *soil* temperature in one
experiment ("across landscapes" leans on Jucker); depauw's warming treatment was a *null*
(its light×structure result carries the community claim, Haesen carries temperature);
meeussen's regeneration effect is species/latitude-dependent.

---

## Local / global coherence (whole-chapter view)

- **Neighbour positioning is complete.** Bouwen (§4.1 dedicated), Gril (§4.2 dedicated,
  same site), Kovács + Park (§4.1 observational counterpart), Schnabel (§4.1), Vinod
  (appendix). The two the chapter positions against, Bouwen and Gril, each get a paragraph.
- **No bare load-bearing claim.** The one novel, contested claim (leaf area warms P1,
  cools P4) is the chapter's own result, correctly stated uncited and corroborated where
  neighbours exist (Kovács, Park, Bouwen, Marsh).
- **No mis-stacked corroboration.** The `[@a; @b]` piles sit on the standard
  buffering-exists framing in the Introduction, where reinforcement belongs.
