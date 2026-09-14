# Cross-check: General Discussion (D.6–D.8) against Appendix G, the standalone perspectives file, the materials map, and the source tables

*Examiner pass, 2026-08-18. Read in full: `discussion_generale_these_EN.md`,
`annexe_GEDI_EN.md`, `disc_generale_perspectives_GEDI_EN.md`,
`disc_generale_materiaux_FR.md`, `discussion_generale_these_EN_harsh_review.md`,
`NC_Full/manuscripts/ch4/Chapter4_prototype_plan_EN.md`, and Tables 0, 5, 6, 10,
12–24 of `NC_Full/manuscripts/ch4/tables/`. Correlations were recomputed from the
CSVs (Tables 23c, 24b, 24c). Items the prior harsh review marks as fixed are not
repeated unless the fix itself introduces a new inconsistency. No file other than
this report was edited.*

---

## A. Number trace (Discussion D.6–D.8 → Appendix G → CSV)

| # | Discussion value (section) | Appendix G location | CSV source | Verdict |
|---|---|---|---|---|
| 1 | GEDI EOS 25–34 d later than S2 (D.7.1) | §G.3, Table G1 (304 vs 279 = 25) | Table19: EOS_trs50 304/279 → 25 d; the 34 d is EOS(S2) 270 vs 304 from the spline mid-fall metric (notebook §5.7), not in Table G1 | OK, but Table G1 supports only 25; the "34" metric is named nowhere in the appendix ("depending on metric") |
| 2 | SOS agreement (3 d, 122 vs 125) | §G.3, Table G1 | Table19 | OK in App G; **absent from D.7** (was in the standalone file) |
| 3 | FUSION_GEDIMAX summer RMSE 1.61, autumn *R*² 0.46 / RMSE 1.57 (D.7.1 "best forcing", numbers only in App G) | §G.5 | Table22: 1.609 / 0.464 / 1.568 | OK |
| 4 | Composition vs Ch3 fusion "significantly better in both windows" (D.7.1) | §G.5: summer ΔRMSE −0.066 [−0.084, −0.048], Δ*R*² +0.034 [+0.013, +0.056]; autumn −0.082 [−0.104, −0.058], +0.194 [+0.064, +0.406] | Table22b, identical | OK (paired *n* = 47 / 46 not stated in App G) |
| 5 | GEDIFALL autumn ΔRMSE −0.08, Δ*R*² +0.19 (D.7.1) | §G.4, Table G2: −0.083 [−0.105, −0.060]; +0.188 [+0.064, +0.404] | Table15b, identical | OK |
| 6 | Summer dip-filling gain (App G only): JJAS −0.066 [−0.085, −0.048], +0.035; August −0.158, +0.273 | §G.4, Table G2 | Table18d rows GEDIMAX − ANNUAL_FIX, identical | OK |
| 7 | GEDIFALL JJAS control −0.021, +0.022 | Table G2 | Table15b | OK |
| 8 | Table G3 (8 condensed rows, summer *R*²/RMSE/bias, autumn *R*²/RMSE) | §G.5 | Table20, all 40 cells match at 2 dp | OK |
| 9 | "every variant where GEDI supplies the per-plot magnitude fails in summer (*R*² ≤ 0.18)" (App G §G.5) | §G.5 | Table20: STATIC_GEDI_DIRECT 0.316, STATIC_GEDI_RATIO 0.49 (GEDI magnitude, not shown in Table G3) | **Mismatch**: holds only for the condensed rows, not the full 15-scenario matrix the sentence invokes |
| 10 | 2022 coupling *r* = −0.88 vs −0.92 (D.7.2) | §G.6, Table G4 | Table23: −0.881 (2022 hot p90) / −0.92 (2021 JJAS all days); 2021 hot p90 = −0.882 | OK numerically; the like-for-like pair (hot vs hot) is −0.882 vs −0.881, which is stronger and not quoted anywhere |
| 11 | Mean slope 0.93 vs 0.90 (D.7.2) | **absent from App G** (Table G4 has no slope row) | Table23: mean_slope 0.925 (2022 hot) / 0.901 (2021 all); 2021 hot p90 = 0.86 | **Discussion number not sourced in App G**; also mismatched pair (hot 2022 vs all-days 2021) |
| 12 | Offset −1.22 (2022 hot) vs −0.73 (2021 summer) (D.7.2) | §G.6, Table G4 | Table23: 2021 hot p90 = −1.13 | OK numerically; **mismatched pair**: 2021 hot days already at −1.13, so the 2022-attributable deepening is 0.09 °C, not 0.49 °C |
| 13 | Bias +1.16 (2022) vs +0.73 (2021) | §G.6, Table G4 | Table23b FUSION_GEDIMAX 1.16; Table22 summer bias 0.734 | OK |
| 14 | Fusion products lead 2022 (*R*² 0.40/0.39; RMSE 1.80–1.81); static ALS 0.08 | Table G4 | Table23b | OK |
| 15 | Vernal *r* = +0.48, partial +0.46, *n* = 31 (D.7.3) | §G.7 | Table23c: 0.482 / 0.464 / 31 rows | Numbers reproduce; **the SOS variable is defective** (see B.1); *n* = 31 not stated in D.7.3 |
| 16 | Winter → SOS partial *r* = 0.05 (D.7.3) | **absent from App G** (§G.7 and §G.9 present the winter test as future work) | Table24c: *r* = 0.11, partial 0.053, *n* = 33 | **Discussion number not sourced in App G**; App G is stale relative to Table24c |
| 17 | Understory greens up ~1 week before overstory (D.7.3) | §G.3 (116 vs 123) | Table19 | OK (Table17b gives 114 vs 119 by spline SOS, 107 vs 119 mid-greenup) |
| 18 | Night *r* = +0.83 vs day −0.89, *n* = 53; mean ΔTmin −0.66, SD 1.24 (D.6) | §G.8 | Table24b: 0.831 / −0.885 / −0.66 / 1.24 | OK numerically; interpretation "dense canopy warms the night" overstated (B.4) |
| 19 | 15–20 loggers suffice (D.7.4, D.8) | §G.8: *n* = 10 median −0.94 [−0.98, −0.79]; *n* = 20 [−0.96, −0.86]; slope ±25% vs ±10% with all 53 | Table24d: *n* = 10 −0.938 [−0.982, −0.789]; *n* = 20 [−0.962, −0.858]; b at 20: −0.095 [−0.118, −0.065] (−24/+32%); *n* = 53 row degenerate (subsampling without replacement) | OK except "±10% with all 53" is not in the table (±11% is the *n* = 40 row); "at a new site" overreaches a within-network subsample |
| 20 | Within ~150 m, 80–90% of ALS microclimate correlation (D.7.1) | §G.5: −0.73/−0.80 vs −0.86/−0.90 | Table10: those are the 100 m GEDI_k3 numbers (*n* = 11); at 150 m GEDI_k3 −0.70/−0.77 vs −0.86/−0.89 (*n* = 18) = 82–86%; single nearest footprint (GEDI_direct) at 150 m = 74–81% | OK for the 3-footprint mean; App G should say "mean of the three nearest power footprints", give *n* = 11–18, and quote the 150 m row if it says 150 m |
| 21 | Montane site: "the same pipeline yields noise" (D.7.1), "height gate and GEDI both weaken" (D.7.4) | **absent from App G** (Aigoual appears only in the footprint count) | Table0: Aigoual cor 0.13 (all beams), 0.12 hi-sens; 0.22 power beams (notebook §3.5, not in a table read); Table5 anchor over-correction bias −0.90 | **Discussion claim not sourced in App G**; also mis-described: the ~150 m microclimate pipeline was never run at Aigoual (no loggers); what fails is footprint-level r(PAI, ALS) and the anchor test |
| 22 | Footprint counts 9190 (516/2448/6226); Blois power beams *n* = 635 | §G.1 | Notebook §2; Table6 n = 635 | OK; Discussion carries no counts (standalone had "~9200") |
| 23 | Amplitude 2.3/3.1/3.3; winter floor 0.5–0.8; plateau 3.4–3.6 (App G §G.2) | §G.2 | Table13 2.31/3.08/3.26; Table12 Mormal Feb–Mar 0.52–0.56, Dec 0.75, Jul 3.59, Aug 3.37 | OK (dense-stratum floor is 0.89, outside "0.5–0.8"; Table21b gives a different amplitude set 2.67/3.21/3.86, provenance to reconcile) |
| 24 | June 6 2021 0.68 ± 0.41 below September (*z* = −3.2) (App G §G.3) | §G.3 | Table14 adjusted 3.07 vs 3.81 → 0.74 | Cannot reproduce 0.68 from Table14 (contrast may come from the lm date term); minor, verify |
| 25 | Cross-site shape *r* = 0.74, max diff 0.28 (App G §G.3) | §G.3 | Table21 (Blois DOY 101–318, Mormal 46–347, no shared DOY grid) | Not reproducible from Table21 as stored (different DOY grids); needs interpolation; provenance note |
| 26 | Density strata SOS 111 vs 123, EOS 305/305/307 (App G §G.3) | §G.3 | Table21b | OK |
| 27 | Scenario dispersion "below 0.28" (D.6) | Ch3, not App G | Table20 SDrec max = 0.280 (FUSION_H) | "≤ 0.28" is exact; "below" is not |

---

## B. Issues, ranked

### Severe

**B.1 The vernal-window correlation rests on a defective SOS metric (D.7.3, App G §G.7).**
`sos_s2` in Table23c is the first DOY at which the raw per-plot Sentinel-2 series
crosses min + half-range (`c4_perspectives_seeds.R` l. 86–89), not the
Whittaker-smoothed TRS50 used in §G.3. Sorted values: 26, 32, 33, 33, 34, 34, 35,
35, 37, 38, 38, 39, 40, 41, 56 | 70, 90, 94, 97, 99, 103, 109, 115, 117, 118,
131, 133, 133, 133, 133, 134. Fifteen of 31 plots "green up" between January 26
and February 25, which no oak stand does; these are the tall plots (cor with Hmax
= −0.40), whose winter S2 LAI floor or a single winter retrieval spike crosses the
half-range early. The *r* = +0.48 is therefore a two-cluster contrast (mean spring
Tmax 16.1 °C for the January cluster vs 18.5 °C for the April–May cluster); within
the physically plausible cluster (SOS ≥ 70, *n* = 16) *r* = 0.06. The winter test
(Table24c) uses the same SOS. Neither +0.48/+0.46 nor 0.05 can be quoted as they
stand.
Quote to retract in D.7.3: "later-leafing plots accumulate a warmer understory
spring (*r* = +0.48, partialled on summer leaf area +0.46). The reciprocal arrow
is absent in our data, since winter microclimate does not predict green-up across
plots (partial *r* = 0.05)."
Fix: recompute SOS on the smoothed per-plot series with the spring-window TRS50
(DOY 60–180, same metric as Table G1), rerun both correlations, and either
re-quote or drop the numbers. The layer-resolved GEDI vernal window (§G.3) is
unaffected.

**B.2 Two Discussion numbers and one claim have no source in Appendix G.**
(a) Slopes 0.93/0.90 (D.7.2) are not in Table G4. (b) The winter partial *r* =
0.05 (D.7.3) contradicts App G, which twice presents that test as not yet done
("sets up the reciprocal test... that the logger archive would support";
"which the reciprocal winter-temperature test would address"). (c) The montane
failure (D.7.1 "at the montane site the same pipeline yields noise"; D.7.4
"where the height gate and GEDI both weaken") appears nowhere in App G, and the
prior harsh review's "Appendix G says so" is false for the current text. The
appendix is billed as sourcing everything; these three break that contract.

**B.3 The 2022 comparison pairs unlike windows and understates its own best evidence (D.7.2, App G §G.6, Table G4).**
Table23 has four rows; the text uses "2021 JJAS all days" against "2022 Jun–Jul
hot p90". Like for like: hot 2021 vs hot 2022 gives ΔTmax −1.13 vs −1.22
(deepening 0.09 °C, not 0.49), *r*(slope, LAI) −0.882 vs −0.881 (identical), mean
slope 0.86 vs 0.925. The Discussion's own argument (a deeper offset is the
arithmetic of a sub-unity slope under hotter skies) is proven by the 2021 hot-day
row and would be stronger if quoted. App G §G.6 still reads "the buffering deepens
on hot days... Dense forest buffers most exactly when it matters most", the
framing the harsh review had the Discussion retract; appendix and discussion now
sit at different epistemic levels on the same table.

### Moderate

**B.4 "Dense canopy cools the day and warms the night" misreads the sign structure (D.6, App G §G.8).**
Table24b: dense plots (LAI > 4) have mean ΔTmin +0.09 °C (parity with the macro
reference); open plots (LAI < 1.5) −3.15 °C. Nothing is warmed; open stands cool
2–4 °C below the reference at night while closed stands sit at parity. The
correlation is robust (Spearman 0.80; *r* = 0.56 on LAI ≥ 2 only) but the
sentence should say "open stands radiate below the reference at night, closed
stands hold it", which is still a bilateral amplitude compression. The macro ΔTmin
comes from the pblh forcing file (`c4_disc_seeds2.R` l. 24), a caveat App G
does not state.

**B.5 App G §G.5 claim "*R*² ≤ 0.18" is false on the full matrix it cites.**
STATIC_GEDI_DIRECT (0.32) and STATIC_GEDI_RATIO (0.49) carry GEDI magnitude and
exceed 0.18 in Table20. Restrict the sentence to the temporal/regressed variants
or quote the exceptions.

**B.6 The 150 m footprint result is under-specified in both documents.**
"the direct GEDI measurement" is the mean PAI of the three nearest power
footprints (GEDI_k3); the quoted −0.73/−0.80 are the 100 m row (*n* = 11); at
150 m *n* = 18 and the recovery is 82–86%; the single nearest footprint gives
74–81%. Sample sizes of 11–18 plots must appear where the 80–90% is claimed.

**B.7 App G §G.8 "15–20 loggers suffice... at a new site" overreaches the design.**
Table24d is subsampling without replacement of one network at one site; the
*n* = 53 row is degenerate, so "±10% with all 53" is not from this table (±11% is
the *n* = 40 row). The result says how many of these 53 loggers reproduce this
site's coupling, not what a new site needs; site-to-site heterogeneity is untested.
D.7.4 and D.8 lean on it for portability ("cheap to test elsewhere").

**B.8 D.7.1 "replicating at site level the continental gap" is one site, one pooled year.**
Blois only (Mormal has no S2 annual series in Table G1), GEDI composites pooling
2021–2022 against a 2021 S2 series. Say "at one site".

**B.9 App G Table G1 vs text.** §G.3 says "nine dates", Table G1 says "8 dates",
Table19 says "adjusted, 9 dates" with n_pts = 8. Pick one formulation ("nine
acquisitions, eight composite points").

**B.10 The "34 d" upper bound of the EOS gap has no visible source.** Table G1
supports 25 d only; the 34 d (S2 270 vs GEDI 304, spline mid-fall metric,
notebook §5.7) is not in the appendix. Add the metric in a footnote to Table G1
or drop the range to "about 25 days (34 with the unsmoothed spline metric)".

**B.11 Table G3 formatting bug.** The last row runs straight into the paragraph
("| 0.30 / 1.60 | All low-bias scenarios carry..."); a blank line is missing, so
the table will not render and the paragraph loses its opening.

### Minor

- D.7.3 carries no *n*; App G has *n* = 31. Add it inline.
- D.7.3 has no Appendix cross-reference; add "(Appendix G, §G.3 and §G.7)".
- D.6 "scenario dispersion ratios below 0.28": FUSION_H SDrec = 0.280 exactly; write "at or below 0.28".
- App G §G.1: "9190 footprints, 2019–2023 acquisitions restricted here to 2021–2022" reads as if 9190 were the 2019–2023 count; the notebook says the rds holds 9190 footprints of 2021–2022.
- App G §G.2 "winter floor 0.5–0.8": the dense stratum is 0.89 (Table13); "0.5–0.9" or "0.5–0.8 by month".
- App G §G.3 "0.68 ± 0.41" is not the difference of the two adjusted means in Table14 (3.81 − 3.07 = 0.74); state which contrast is quoted.
- App G §G.3 cross-site *r* = 0.74 / max diff 0.28 cannot be recomputed from Table21 as stored (different DOY grids); note the interpolation.
- Two amplitude sets exist in the folder (Table13: 2.31/3.08/3.26; Table21b: 2.67/3.21/3.86); App G uses Table13. Reconcile or annotate provenance before defense.
- Fig10_c4_crosssite_shape.png exists and matches the §G.3 cross-site paragraph but is not referenced.
- App G §G.5 paired *n* (47 summer, 46 autumn) is stated in Table G2 but not for the composition contrasts.
- Style: em-dashes survive in App G headers and all eight figure/table captions ("Figure G1 — ", "Table G1 — ") and in the standalone title; house style forbids them (colon or period instead). The `[REF: Grulois 2026 — ...]` placeholder is the prescribed format and is fine, but it is still unresolved (materials map §3.9). Bold "**+0.83**" in body prose (§G.8) is off-style. US spelling checks clean; italic stats consistent; hedge stacking in the D.6 night paragraph (may / although / not as / cannot / applies) is one hedge too many, and the paragraph's cadence is uniformly long.
- Out of my remit but flagged: D.2 "sensor crossover at LAI ≈ 3.9" is, per project memory, a value whose derivation from Ch3 tables was not reproducible (`project_ch3_v320_v323_trap`); confirm before the defense copy.

---

## C. Suggested exact rewrites

**C.1 D.7.2, deepening sentence.**
Old: "The mean offset did deepen (−1.22 °C on 2022 hot days against −0.73 °C over summer 2021), but a deeper offset under hotter skies is the expected arithmetic of a sub-unity slope, not evidence of adaptation; the invariant coupling is the finding."
New: "The mean offset was deeper on hot days (−1.22 °C on the hottest 2022 days, against −1.13 °C on the hottest 2021 days and −0.73 °C over all of summer 2021), which is the expected arithmetic of a sub-unity slope under hotter skies rather than evidence of adaptation; the coupling itself did not move (*r* = −0.88 in both years' hot windows), and that invariance is the finding."

**C.2 D.7.2, slope clause.**
Old: "the mean buffering slope did not collapse (0.93 against 0.90)"
New: "the mean buffering slope did not collapse (0.93 on the hottest 2022 days against 0.86 on the hottest 2021 days; Appendix G, Table G4)"
and add to App G Table G4 a row: "Mean micro–macro slope | 0.90 (all days) / 0.86 (hot p90) | 0.93".

**C.3 App G §G.6, framing.**
Old: "and the buffering deepens on hot days (mean ΔTmax −0.73 °C over summer 2021, −1.22 °C on 2022 hot days), consistent with the GEDI observation that the 2022 canopy itself was intact (§G.3). Dense forest buffers most exactly when it matters most, the premise of the microrefugia literature"
New: "and the mean offset on hot days is deeper (−1.13 °C on 2021 hot days, −1.22 °C on 2022 hot days, against −0.73 °C over all of summer 2021), as a sub-unity slope under hotter skies implies; the GEDI composites show the 2022 canopy itself intact (§G.3). The invariant coupling, not the deeper offset, is the evidence relevant to the microrefugia premise"
Table G4: replace the two-row block by four rows (2021 all, 2021 hot p90, 2022 all, 2022 hot p90) for ΔTmax, *r*(slope, LAI), and mean slope.

**C.4 D.7.3, pending recomputation (B.1).** If the smoothed-TRS50 recomputation confirms a positive association:
Old: "later-leafing plots accumulate a warmer understory spring (*r* = +0.48, partialled on summer leaf area +0.46). The reciprocal arrow is absent in our data, since winter microclimate does not predict green-up across plots (partial *r* = 0.05)."
New: "later-leafing plots accumulate a warmer understory spring (*r* = [x], partialled on summer leaf area [y]; *n* = 31, green-up dated by the same smoothed half-amplitude threshold as Appendix G, Table G1). The reciprocal association is absent in our data: winter understory temperature does not predict green-up across plots (partial *r* = [z], *n* = 33; Appendix G, §G.7)."
If it does not: delete both numbers and keep only the GEDI layer lead sentence, ending "the thermal consequence of that window is the next test the logger archive supports."
App G §G.7: replace "sets up the reciprocal test... that the logger archive would support" by the winter result, and §G.9 "which the reciprocal winter-temperature test would address" by "the reciprocal winter test (§G.7) is null but shares the cross-sectional design".

**C.5 D.6, night sentence.**
Old: "dense canopy cools the day and warms the night, a bilateral amplitude compression under one structural control."
New: "closed canopy holds the night at the reference temperature while open stands fall 2–4 °C below it, the mirror of the daytime pattern: a bilateral amplitude compression under one structural control."
App G §G.8: same substitution, plus "the macro minimum is taken from the MuSICA forcing file, so cold-air pooling at the reference is not excluded."

**C.6 D.7.1, montane clause, and App G §G.9 addition.**
Old (D.7.1): "at the montane site the same pipeline yields noise, and that boundary is part of the perspective."
New: "at the montane site the footprint-level agreement with ALS is itself weak (*r* = 0.13–0.22, September only, conifer admixture; Appendix G, §G.9), so no such archive-based check is available there, and that boundary is part of the perspective."
App G §G.9, add: "Aigoual is a stated failure case: 516 footprints, September only, relief and conifer admixture; footprint-level *r*(PAI, ALS LAI) = 0.13 (all beams) and 0.22 (power beams), and a GEDI-anchored magnitude correction over-corrects (bias −0.90 m² m⁻²). All microclimate results are Blois-only by construction."

**C.7 App G §G.5, "≤ 0.18" sentence.**
Old: "every variant where GEDI supplies the per-plot magnitude fails in summer (*R*² ≤ 0.18)"
New: "every temporal variant where GEDI supplies the per-plot magnitude fails in summer (*R*² ≤ 0.18); the two static GEDI-magnitude scenarios score higher (0.32 direct at footprints, 0.49 ratio-scaled S2) but carry the largest warm biases of the matrix (+1.08, +1.23 °C)"

**C.8 App G §G.5 and D.7.1, footprint clause.**
Old (App G): "Within ~150 m of a power footprint, the direct GEDI measurement additionally recovers 80–90% of the ALS microclimate correlation (*r* = −0.73/−0.80 vs −0.86/−0.90 for plot ALS)"
New: "For the 11–18 plots within 100–150 m of a power footprint, the mean PAI of the three nearest footprints, used directly, recovers 80–90% of the ALS microclimate correlation (*r* = −0.70 to −0.73 / −0.77 to −0.80 against −0.86 / −0.89 to −0.90 for plot ALS)"
D.7.1: insert "(*n* = 11–18 plots)" after "80–90% of the ALS microclimate correlation".

**C.9 App G §G.8, loggers.**
Old: "and the regression coefficient is constrained to about ±25%, against ±10% with all 53. Fifteen to twenty well-spread loggers suffice to establish the coupling at a new site; the full network refines the slope."
New: "and the regression coefficient is constrained to about ±25–30%, against ±11% at *n* = 40. Within this network, fifteen to twenty well-spread loggers reproduce the coupling; whether that transfers to another site is untested, since the subsample inherits the site's own gradient."
D.7.4: "so the validation design of this thesis is portable at modest cost" → "so the validation design of this thesis is, at least within one network, reproducible at modest cost".

**C.10 D.7.1, replication clause.** "replicating at site level" → "replicating at one site".

**C.11 App G, captions and title.** "Figure G1 — GEDI monthly..." → "Figure G1. GEDI monthly..."; "Table G1 — TRS50..." → "Table G1. TRS50..."; "Appendix G — Exploratory analysis: GEDI..." → "Appendix G. Exploratory analysis: GEDI...". Insert a blank line after the last row of Table G3.

**C.12 D.8.** "its first two trials" → "its first two exploratory trials".

---

## D. The standalone perspectives file (`disc_generale_perspectives_GEDI_EN.md`)

**Recommendation: archive** (rename with a `_superseded` suffix or move to an
`archive/` folder; do not delete). D.7 supersedes it in structure and in caveats:
the pivot against the General Introduction (weaker requirement), the selection
caveat, the montane boundary, the sub-unity-slope reading of 2022, the winter
null, the Greiser nuance, and the logger-count note exist only in D.7. Keeping
both would duplicate ~70% of the content and would carry the older, less
guarded 2022 framing ("buffering deepens exactly when it matters most") into the
thesis twice.

Sentences in the standalone that D.7 currently lacks and should be ported before archiving:

1. SOS agreement: "Its start of season at Blois agrees with the Sentinel-2 series to three days, consistent with optical green-up being reliable [@soudaniEvaluationOnsetGreenup2008]" (D.7.1 cites Soudani for the shoulders but never states the leaf-out agreement, which is the half of the timing result that vindicates Sentinel-2).
2. Leaf/wood separation: "a winter floor of 0.5–0.8 m² m⁻² that is plausibly wood and branch area... A waveform sensor thus separates leaf from wood by season" (D.7.1 says only that the amplitude ranks density).
3. Marcescence caveat: "the late GEDI end of season mixes true leaves, marcescent leaves (typical of *Quercus* and *Fagus*), and exposed wood, although for radiation and microclimate purposes this structural signal is arguably the relevant one, which is precisely what the logger test supports" (absent from the whole Discussion; a jury member on phenology will ask).
4. Complementarity to wall-to-wall height products with the FORMS citation: "it complements rather than replaces the wall-to-wall products of Chapters 2 and 3 [@schwartzFORMSForestMultiple2023]" (D.4 cites FORMS for the gate; D.7.1 loses the link between GEDI's footprint-bound skill and the wall-to-wall products).
5. Optional closing image: "magnitude from LiDAR, timing from optics, and now a LiDAR check on the timing itself" (D.7.1's matrix sentence covers the content; keep only if D.7.1 wants a landing line).

Sentences that would be lost and should stay lost: the summer dip-smoothing
gain (a small, non-pre-specified effect; App G §G.9 says to treat it as
suggestive), the "~9200 footprints" count (belongs to App G), the "buffering
deepens exactly when it matters most" framing (retracted).

What D.7 says that the standalone does not (for the record, all keep): the
"weaker requirement" pivot; the selection-conditioned CI caveat; the montane
failure; the sub-unity-slope arithmetic; the slope 0.93/0.90 (to be re-sourced);
the winter null (to be recomputed); the Greiser nuance; D.7.4 in full.

---

## E. Materials-map red list: status in the Discussion

| # | Red-list item | Status | Note |
|---|---|---|---|
| 1 | No categorical "fusion value is spatial not temporal" | Respected | D.7.1 frames the gap as "no LiDAR truth outside summer", not as a Ch3 negative |
| 2 | Pivot against Intro §1.7 (GEDI = weaker requirement) | Respected | D.7.1 "The role it can hold is a different and weaker requirement" |
| 3 | No "systematic underestimation of LAI_S2" | Respected | D.1 "dynamic range rather than central tendency", means comparable at *k* = 0.65 |
| 4 | Confront Ch2 §4.5 "may serve as operational input" | Respected | D.3, turned as predicted anti-corollary |
| 5 | Height as classifier, not lever | Respected | D.4 |
| 6 | ΔTmin reconciled (sign, mechanism, Ch1 guard, two-sided ecology) | Partly | D.6 has guard and two-sidedness; sign reading overstated (B.4) |
| 7 | Collinearity trap stated once for the seeds | Respected | D.6, explicit; D.7.3 names confounds |
| 8 | 2022 as nuance, not refutation; absolute skill degrades | Respected in intent | Greiser nuance and +1.16 present; pairing of windows still inflates the deepening (B.3) |
| 9a | d_opt "6–10 m, 7 m combined" | Respected | D.1 |
| 9b | [REF: Grulois 2026] resolved | Not resolved | Placeholder still in D.7.4 |
| 9c | *r* = 0.50 (never 0.46) | Respected | D.1, D.3 |
| 9d | Quantity = LAI + cover together | Respected | D.1 |
| Étage 3 | No wall-to-wall map, no absolute accuracy, no attribution from correlations | Respected | D.8, D.5, D.6 |

---

## F. Verdict (5 lines)

1. All 40 cells of Table G3, all bootstrap contrasts (G2, G.5), the 2022 scenario table, the night and logger numbers reproduce from the CSVs; the appendix is numerically clean where it quotes.
2. Three Discussion items are unsourced by Appendix G (slopes 0.93/0.90, winter partial 0.05, montane failure); the appendix must be updated to carry them or the Discussion must drop them.
3. The vernal-window correlation (+0.48/+0.46, and the 0.05 null) is built on a defective SOS (15 of 31 plots "green up" in January); recompute with the smoothed TRS50 before it reaches the thesis, otherwise retract D.7.3's numbers.
4. The 2022 stress test compares hot 2022 days with all 2021 days; the like-for-like rows (−1.13 vs −1.22; *r* −0.882 vs −0.881) make the Discussion's own argument better and should replace them in both documents, with Appendix G §G.6 aligned to the Discussion's retracted-deepening framing.
5. Archive the standalone perspectives file after porting four sentences (SOS agreement, leaf/wood separation, marcescence caveat, FORMS complementarity); the D.7 structure is the right one and the red list is otherwise honored.
