# Cross-check: General Discussion vs Chapter 3

*Examiner pass, 2026-08-18. Target: `thesis_frame/discussion_generale_these_EN.md` (draft 2026-08-13,
post-harsh-review). Reference: `NC_Full/manuscripts/ch3/Chapter3_article_standalone_EN.md`
(manuscript text taken as truth; `tables/figK_hybrid_ch2.csv` and `figJ_stratified_cross386.csv`
consulted only to confirm). Items already marked fixed in `discussion_generale_these_EN_harsh_review.md`
(oracle numbers re-attributed, 2022 arithmetic, "accurate" → consistency, three-lines provenance,
spring causal verb, night adjudication, "not a skill ranking", footprint failure case, D.8 wording)
are NOT repeated; I re-verified that each fix is present in the current draft.*

Line numbers below refer to `sed -n` numbering of each file (Ch3 = "L", Discussion = "D").

---

## (A) Claim-by-claim table

| # | Discussion quote (§, D-line) | Ch3 quote (§, L-line) | Verdict |
|---|---|---|---|
| 1 | D.1 (D45–46) "dense-canopy ΔTmax *R*² = 0.60 against 0.02 for Sentinel-2" | 3.2 (L91) "ΔT~max~ *R²* = 0.60, 95 % CI 0.43–0.82, against 0.02, CI 0.00–0.21" | MATCH |
| 2 | D.1 (D46–47) "Sentinel-2 suffices where the canopy is open (0.71 against LiDAR near zero)" | 3.2 (L91) "in open canopy Sentinel-2 won (ΔT~max~ *R²* = 0.71 ... against near-zero for LiDAR, CI 0.00–0.09)" | MATCH |
| 3 | D.1 (D48–50) "In its ideal, LAI-keyed form the switch holds the lowest bias (+0.67 °C) and RMSE (1.60 °C) of the products tested" | 3.4 (L103) "lowest warm bias among the skill-competitive products ... warm bias +0.67 °C, RMSE 1.60 °C (against 1.75–1.76 °C for the single-sensor products)"; 3.4 (L117) oracle "RMSE of 1.57 °C" | VALUES MATCH; SCOPE drift: Ch3 qualifies "among the skill-competitive products" (naive open-stratum bias +0.08 °C, oracle RMSE 1.57 °C exist). Minor. |
| 4 | D.1 (D50–52) "gated on a canopy-height product at 17.8 m, reproduces the assignment for 85% of plots ... (*R*² = 0.54, bias +0.71 °C, RMSE 1.62 °C), stable under leave-one-out" | 3.4 (L103) "FORMS-H threshold of 17.8 m (which reproduces the LAI split for 85 % of plots ...) was nearly identical at ΔT~max~ *R²* = 0.536 ... stable under leave-one-out"; Table 2 (L127) FORMS-H-gated bias +0.71; 3.4 (L117) "1.62 °C for the operational FORMS-H switch" | MATCH (0.536 → 0.54 rounding). NB `figK_hybrid_switch.csv` gives 0.550/0.540 for a second switch implementation; the manuscript uses the `figK_hybrid_ch2` values, consistent with the Discussion. |
| 5 | D.1 (D52–54) "The pooled leaderboard that seemed to favor Sentinel-2 is a composition effect: the pooled lead is real, but it is carried entirely by the open half" | 3.2 (L91) "Pooled ... Sentinel-2 appeared ahead (0.56 ... against 0.24) but this pooled lead is a Simpson effect: it arises entirely from the open plots"; 4.2 (L146) "the pooled *R²* lead of Sentinel-2 is a genuine ranking result, but it is a Simpson effect" | MATCH in substance. Term drift: Ch3 says "Simpson effect" throughout; Discussion says "composition effect" and never names Simpson. |
| 6 | D.1 (D45) "LiDAR magnitude is indispensable exactly where buffering is deepest" | Abstract (L19) / 5 (L162) "LiDAR remains indispensable, among the products tested, where the canopy is dense"; 4.3 (L154) "the dense-canopy requirement is then one for some form of LiDAR, airborne or spaceborne, not for airborne LiDAR specifically" | SCOPE: Ch3 attaches "among the products tested" every time; Discussion drops it. Minor (D.4/D.7.1 partly recover it). |
| 7 | D.2 (D62–63) "the sensor crossover at LAI ≈ 3.9 (the stratified forcing comparison of Chapter 3)" | 2.3 (L51) "That threshold is a split of convenience, not a quantity we derive: it lies close to the median LiDAR LAI ... Nothing in what follows depends on its exact value"; 3.2 (L91) "held across thresholds from LAI 2.0 to 4.5 ... the 3.86 split is a convenience placed at the sample median" | OVER-READ: Ch3 does label it "sensor crossover" (Table 2, App. C) but explicitly refuses to treat 3.86 as an estimated crossover. In the "three lines" geometry the Discussion uses it as a located point. Moderate. |
| 8 | D.2 (D64–65) "the optical saturation threshold, placed at LAI 5–6 in the literature, with our own pooled Sentinel-2 response flattening near 4" | 1 (L27) "the retrieved Sentinel-2 LAI saturates in dense canopy, with a ceiling of 6–7 m²/m²"; 4.1 (L139) "only well above both does the Sentinel-2 retrieval reach its magnitude ceiling of 6–7 m²/m². The optical ranking therefore collapses well before the magnitude ceiling" | MISMATCH of value and currency: Ch3's own third line is the retrieval ceiling 6–7, not "5–6 (literature)" nor "near 4"; the "near 4" plateau is a mean-per-bin quantity at *k* = 0.5 (harsh review already flagged it as unverified against Ch2). Moderate. |
| 9 | D.2 (D69–71) "63% of the leaf area sits below the effective optical depth, and Chapter 3 showed that this hidden layer buffers as much per unit leaf area as the visible top" | 3.3 (L95) "in the dense plots 63 % of the leaf area sits below *d*~opt~ ... the sub-canopy that Sentinel-2 cannot see buffers as much as the top" | 63 % MATCH. "per unit leaf area" is NOT in Ch3 (Ch3's evidence is the +1.14 vs +0.73 °C bias contrast, not a per-unit slope). Moderate upgrade. |
| 10 | D.2 (D77–79) "the fusion rule of Chapter 3 is a genuine resolution rather than a compromise: each half of the domain is served by the instrument that resolves it" | 4.1 (L142) "The validity map is asymmetric in its support. The dense-canopy half ... robust ... model-free. The open-canopy half ... more contingent ... established only through the coupled model"; 4.1 (L138) Sentinel-2 in open stands "is biased in magnitude there and ... largely redundant with canopy height"; 3.5 (L131) fusion "was not significantly better than Sentinel-2 alone on the pooled ΔT~max~ ranking" | STRONGER than Ch3: symmetric "resolution" where Ch3 insists on an asymmetric, partly contingent verdict and a pooled tie. Moderate. |
| 11 | D.3 (D97–99) "Chapter 3 found that this same product [forest-tuned opt] is the worst microclimate forcing of the set: pooled skill indistinguishable from zero and ... sign-inverted (*r* = +0.19)" | 3.4 (L107) "'opt' retrieval, applied across all 53 plots, ranked near zero (pooled ΔT~max~ *R²* = 0.00); restricted to its calibration domain (39 plots) it was no better than the raw ATBD retrieval (*R²* = 0.07 against 0.01)"; 4.1 (L142) "this penalty is in part one of domain misuse"; Table 2 note (L119) opt "omitted from this table because ... its 53-plot value would not be comparable to the others"; 3.2 (L89) "+0.19, with the sign inverted" | *r* = +0.19 and pooled ≈ 0 MATCH. "worst ... of the set" is NOT a Ch3 statement: Ch3 says in-domain it equals ATBD (both saturated), calls the 53-plot value non-comparable, and attributes part of the penalty to domain misuse. Moderate. |
| 12 | D.3 (D99–101) "Truncating the LiDAR magnitude to the optical depth fails the same way, with the largest warm bias of all scenarios (+1.14 °C against +0.73 °C for the full column)" | 3.3 (L95) "carried the largest warm bias (+1.14 °C)"; "LiDAR (+0.73 °C)" | MATCH (Ch3's own wording; note Ch3 also reports stratum biases +1.86 and +2.29 °C for other products, so "of all scenarios" is pooled-only). |
| 13 | D.3 (D108–110) "the third established that the route, not the tuning, is the constraint" | App. B (L207) raw-reflectance RF "ranked the dense stratum at ΔT~max~ *R²* = 0.35 ... above the retrieved Sentinel-2 LAI (0.02) but below LiDAR (0.60) ... the deficit is a property of the retrieved leaf area rather than of the reflectance itself" | SLIGHTLY STRONGER: Ch3 locates part of the loss in the retrieval (reflectance retains dense information the PROSAIL retrieval discards). "Route not tuning" is defensible only with that nuance. Minor–moderate. |
| 14 | D.4 (D121–122) "the leaf quantity that optical indices stop resolving beyond LAI 5–6" | 4.1 (L139) "the optical between-plot ranking is already lost throughout the dense stratum defined at the near-median split (LAI ≥ 3.86, *R²* = 0.02) ... collapses well before the magnitude ceiling, so it is the ranking, not the ceiling, that sets the domain of validity" | MISMATCH: Ch3's operative loss is the ranking loss at ≥ 3.86 (and the reversal from 2.0), not a 5–6 ceiling. Also internally inconsistent with D.2's crossover at 3.9. Moderate. |
| 15 | D.4 (D124–125) canopy height "is an excellent classifier ... gates the fusion rule (17.8 m) because it proxies the density regime" | 3.2 (L89) "canopy height from FORMS-H moderately (*r* = −0.64)"; 3.4 (L101) height-only RF *R²* = 0.585, "the height carries the operational signal"; 3.4 (L103) 85 % agreement; 3.5 (L131) Aigoual only 31 % assigned | "excellent" is stronger than Ch3's "moderately"/85 %/site-dependent. Minor. |
| 16 | D.5 (D141–142) "functional arbitration (Chapter 3), in which 53 independent temperature loggers judge LAI products" ... (D147–148) "its non-independent loggers make the reported *R*² optimistic" | 2.5 (L79) "the plots are not fully independent spatially (pseudoreplication), so the effective *n* is below 53"; 4.3 (L152) "the loggers ... are not spatially independent" | *n* = 53 MATCH. "independent" contradicts Ch3 and the Discussion's own next sentence (intended sense: independent of the LAI products). Moderate wording. |
| 17 | D.5 (D145–147) "the consistency analysis carries a bounded prior circularity (≤ 17% ...)"; D.1 "regression slopes remain at 0.25–0.80" | Not in Ch3 (Ch2 numbers). | Out of Ch3 scope; previously verified in the harsh review. Not re-checked here. |
| 18 | D.6 (D155–157) "scenario dispersion ratios below 0.28 in Chapter 3", "nighttime warm bias of +1.5–1.7 °C" | 3.3 (L97) "SDrec ≤ ~0.28"; 3.5 (L131) "large warm bias (+1.5 to +1.7 °C)" | MATCH. Term drift: Ch3 = "amplitude recovery SDrec"; Discussion = "dispersion ratios". |
| 19 | D.6 (D169–171) "Chapter 3 treated the night as a sensor-agnostic regime that leaf area ranks only weakly; the correlation above suggests that part of that weakness may sit in the simulated field" | 4.2 (L148) "nighttime cooling is set by radiative loss and cold-air drainage that the leaf-area field ranks only weakly and the coupled model under-predicts"; but 3.5 (L131) "the two products ranking the plots comparably (pooled ΔT~min~ *R²* = 0.68 to 0.73) but with a large warm bias" | INTERNAL Ch3 tension inherited: Ch3's own night numbers (pooled ΔT~min~ *R²* 0.68–0.73, higher than daytime pooled) do not describe a field that "ranks only weakly"; the sentence "weakness may sit in the simulated field" then has no premise. Moderate. |
| 20 | D.7.1 (D190–193) "Every LiDAR quantity in the preceding chapters is a summer snapshot"; (D198–201) "Re-running the Chapter 3 machinery with the autumn limb ... improved the simulated microclimate significantly (autumn ΔRMSE −0.08 °C, Δ*R*² +0.19)" | 2.1 (L39) loggers "recorded hourly through the summer of 2021 (1 June to 30 September)"; App. A (L173) "Because the HOBO loggers cover the summer only, this test is ... a model-to-model scenario comparison, without ground truth, over the leaf-out and autumn windows ... The temporal skill of a resolved Sentinel-2 dynamic is therefore not testable beyond summer with these data" | CONTRADICTION: Ch3 tells the reader autumn cannot be validated against loggers; the Discussion reports a significant autumn logger-scored improvement (App. G) on the same site and loggers. Severe (see B1). Also *n*: App. G's autumn paired bootstrap uses 47 series plots (memory), Discussion says "the same 53 plots" (D204). |
| 21 | D.7.1 (D206–208) sensor matrix: "ALS supplies the magnitude, Sentinel-2 the plot-scale summer timing, and GEDI the seasonal arbitration" | App. A (L173) "the dynamic scenario converges to the static one on the foliar plateau ... The value of the fusion is spatial and regime-dependent, not temporal"; App. B (L199) "0.24 against 0.23 ... the deficit is structural, not temporal"; App. A "at leaf-out the parametric summer phenology matched or outperformed the observed Sentinel-2 timing" | MISMATCH: Ch3 finds Sentinel-2 timing worthless in summer (plateau) and no better than parametric at leaf-out; the Discussion assigns "plot-scale summer timing" to Sentinel-2 as its role. Moderate. |
| 22 | D.7.2 (D223–224) "*r* = −0.88 against −0.92 in 2021"; "mean buffering slope did not collapse (0.93 against 0.90)" | 3.1 (L87) "observed *r* = −0.92" (slope vs LAI, spatial) | MATCH for the 2021 anchor. (App. G values not checked.) |
| 23 | D.7.4 (D268–271) "Reprocessing the raw point clouds to retain sub-2 m returns would test the one mechanism we could not, the understory vegetation that may explain why Sentinel-2 ranks open plots"; "Sequential data assimilation ... the one untested fusion family"; "montane site ... where the height gate ... weaken[s]" | 4.3 (L152) "testing whether retaining sub-2 m returns would help requires reprocessing the raw point clouds"; 4.3 (L158) "the one fusion family our elimination did not test ... sequential data assimilation"; 3.5 (L131) Aigoual 31 % | MATCH |
| 24 | D.8 (D291–292) "a regime rule that says which sensor to trust on which half of the density gradient" | Abstract (L19) "an explicit domain of validity: a canopy-height-gated switch" | MATCH in substance; Discussion never uses Ch3's headline term "domain of validity" (see D-terminology). |
| 25 | Retired phrases: "co-lead", "artifact", "skill ranking" | grep of Discussion: none of the three appears | CLEAN. |
| 26 | VPD | Not in Ch3; not in Discussion | No cross-check possible; known gap (harsh review §9). |
| 27 | Plots *n* = 47 (series plots), *n* = 31 (spring) | Ch3 uses 53 throughout, never 47; 31 appears only as Aigoual "31 %" | Discussion cites neither 47 nor 31. D.7.3 should state *n* = 31; D.7.1 should state 47 if App. G's autumn/paired bootstrap is on 47. |

---

## (B) Contradictions, ranked

### Severe

**B1. Autumn validation exists in the Discussion but is declared impossible in Chapter 3.**
Discussion D.7.1: "Re-running the Chapter 3 machinery with the autumn limb of the seasonal series
constrained by GEDI improved the simulated microclimate significantly (autumn ΔRMSE −0.08 °C, Δ*R*²
+0.19)". Chapter 3 §2.1: loggers "recorded hourly through the summer of 2021 (1 June to 30 September)";
App. A: "Because the HOBO loggers cover the summer only, this test is ... a model-to-model scenario
comparison, without ground truth, over the leaf-out and autumn windows ... not testable beyond summer
with these data". A jury member reading Ch3 then D.7.1 will ask what the autumn ΔRMSE is computed
against. The logger archive in fact runs 2020-07 → 2023-08 (project memory), so Ch3's App. A sentence
is a scoping choice mis-stated as a data limit. Either Ch3 App. A must say "the summer window retained
for this chapter" and admit the autumn loggers exist, or the Discussion must say explicitly that
Chapter 3 restricted itself to summer and that Appendix G is the first use of the autumn logger
records (and give the *n*: 47 series plots, not 53, if that is what App. G used).

### Moderate

**B2. Night: "ranks only weakly" vs pooled ΔTmin *R*² = 0.68–0.73.**
Discussion D.6: "Chapter 3 treated the night as a sensor-agnostic regime that leaf area ranks only
weakly; the correlation above suggests that part of that weakness may sit in the simulated field".
Ch3 §3.5 reports "pooled ΔT~min~ *R²* = 0.68 to 0.73" (both sensors), i.e. the simulated night field
ranks the plots *better* than the daytime pooled field; only §4.2's gloss says "ranks only weakly".
The Discussion inherits Ch3's weaker gloss and then reasons from it. Recast around what Ch3 actually
measured: reversal dissolves, both sensors rank comparably, warm bias +1.5–1.7 °C.

**B3. Where does the optical signal fail: 3.9, 5–6, or 6–7?**
D.2 places the crossover at ≈ 3.9 and the saturation line at "5–6 (literature) / near 4 (ours)"; D.4
says optical indices "stop resolving beyond LAI 5–6". Ch3 §4.1: ranking "already lost throughout the
dense stratum (LAI ≥ 3.86)", magnitude ceiling "6–7 m²/m²", and "it is the ranking, not the ceiling,
that sets the domain of validity". Three different numbers for one loss in the Discussion, none of
them Ch3's 6–7. Use Ch3's ordering verbatim (breakpoint 2.89 → ranking lost from the near-median split,
stable 2.0–4.5 → magnitude ceiling 6–7).

**B4. Symmetric "resolution" vs Ch3's asymmetric, tied verdict.**
D.2: "a genuine resolution rather than a compromise: each half of the domain is served by the
instrument that resolves it". Ch3 §4.1: "The validity map is asymmetric in its support ... The
open-canopy half ... established only through the coupled model and resting on a mechanism we cannot
test"; §3.5: fusion "not significantly better than Sentinel-2 alone on the pooled ΔT~max~ ranking";
§4.1: Sentinel-2 in open stands "largely redundant with canopy height". The Discussion states the
open half in Ch3's dense-half register.

**B5. "Worst microclimate forcing of the set" (opt retrieval).**
D.3 vs Ch3 §3.4/§4.1: in-domain 0.07 vs ATBD 0.01, "in part one of domain misuse", 53-plot value
"would not be comparable to the others". Ch3 makes the point as a two-sided negative control (accuracy
against LiDAR ≠ microclimate skill), not as a leaderboard bottom.

**B6. Sentinel-2's role = "plot-scale summer timing".**
D.7.1 matrix vs Ch3 App. A/B: summer timing adds nothing (0.24 vs 0.23; 0.56 vs 0.54), the parametric
phenology matched or beat S2 timing at leaf-out, "the value of the fusion is spatial and
regime-dependent, not temporal". If the matrix keeps a temporal role for Sentinel-2 it must be
scoped to the shoulders and paired with Ch3's summer null.

**B7. "53 independent temperature loggers" (D.5).**
Ch3 §2.5/§4.3 and the Discussion's own next sentence say the loggers are not spatially independent.
The intended sense (independent of both LAI products) must be made explicit.

**B8. "hidden layer buffers as much per unit leaf area as the visible top" (D.2).**
Ch3 says only "buffers as much as the top", evidenced by the +1.14 vs +0.73 °C bias contrast. A
per-unit-LAI slope equality exists in project tables (β_below ≈ β_top) but is not in the Ch3
manuscript. Either add it to Ch3 §3.3 or drop "per unit leaf area".

**B9. "sensor crossover at LAI ≈ 3.9" as a located point (D.2).**
Ch3 §2.3/§3.2 repeatedly: "a split of convenience ... Nothing in what follows depends on its exact
value ... held across thresholds from LAI 2.0 to 4.5". The three-lines figure needs that clause or the
alignment is built on a median.

### Minor

**B10.** D.1 "lowest bias and RMSE of the products tested" → Ch3 "among the skill-competitive
products"; oracle RMSE 1.57 °C is the ceiling.
**B11.** D.1 "LiDAR magnitude is indispensable" → Ch3 always "among the products tested" and "some
form of LiDAR, airborne or spaceborne".
**B12.** D.4 "excellent classifier" → Ch3 "moderately (*r* = −0.64)", 85 %, Aigoual 31 %.
**B13.** D.3 "the route, not the tuning, is the constraint" → Ch3 App. B: raw reflectance recovers
dense *R*² = 0.35 (vs 0.02 retrieved, 0.60 LiDAR): the retrieval discards part of the signal.
**B14.** Table provenance: `figK_hybrid_switch.csv` carries a second switch implementation
(LAI-oracle 0.550 / +0.68 / 1.61; FORMS-H 0.540 / +0.70 / 1.615). The manuscript and Discussion use
`figK_hybrid_ch2.csv` (0.573 / +0.67 / 1.60; 0.536 / +0.71 / 1.62). Consistent with each other; make
sure the thesis figure set draws from the same file.

---

## (C) Suggested exact rewrites (Discussion only)

1. **D.1** OLD "LiDAR magnitude is indispensable exactly where buffering is deepest"
   → NEW "LiDAR magnitude is, among the products tested, indispensable exactly where buffering is
   deepest".

2. **D.1** OLD "In its ideal, LAI-keyed form the switch holds the lowest bias (+0.67 °C) and RMSE
   (1.60 °C) of the products tested;"
   → NEW "In its ideal, LAI-keyed form the switch holds the lowest bias (+0.67 °C) and RMSE (1.60 °C)
   among the skill-competitive products, within 0.03 °C of a per-plot oracle ceiling (1.57 °C), while
   remaining statistically tied with Sentinel-2 alone on the pooled ranking (Δ*R*² = 0.010, CI −0.11
   to +0.20);"

3. **D.1** (add after "carried entirely by the open half of the gradient.")
   → NEW "The two halves do not carry equal weight: the dense-canopy verdict is model-free (within-
   stratum *r* = −0.74) and survives every robustness check, whereas the open-canopy usability of
   Sentinel-2 is established only through the coupled model and rests on an untested understory
   mechanism."

4. **D.2** OLD "the sensor crossover at LAI ≈ 3.9 (the stratified forcing comparison of Chapter 3);
   and the optical saturation threshold, placed at LAI 5–6 in the literature, with our own pooled
   Sentinel-2 response flattening near 4."
   → NEW "the sensor reversal, stratified at a near-median split of LAI 3.86 in Chapter 3 but stable
   for any split from 2.0 to 4.5; and the retrieval ceiling of the Sentinel-2 leaf area at 6–7 m² m⁻²
   (Chapter 2, restated in Chapter 3), which the optical ranking loses well before reaching."

5. **D.2** OLD "this hidden layer buffers as much per unit leaf area as the visible top"
   → NEW "this hidden layer buffers as much as the visible top: discarding it costs +0.41 °C of warm
   bias (+1.14 against +0.73 °C)". (Or add the β_below ≈ β_top slope test to Ch3 §3.3 and keep the
   per-unit wording with a pointer.)

6. **D.2** OLD "the fusion rule of Chapter 3 is a genuine resolution rather than a compromise: each
   half of the domain is served by the instrument that resolves it."
   → NEW "the fusion rule of Chapter 3 is a division of labor rather than a compromise: the dense half
   is served by the only instrument that resolves it, and the open half by the cheaper one, where
   canopy height alone already carries most of the signal."

7. **D.3** OLD "Chapter 3 found that this same product is the worst microclimate forcing of the set:
   pooled skill indistinguishable from zero and,"
   → NEW "Chapter 3 found that this same product is no better a microclimate forcing than the untuned
   retrieval it was meant to improve on (in-domain *R*² = 0.07 against 0.01, pooled skill
   indistinguishable from zero) and,"

8. **D.3** OLD "the third established that the route, not the tuning, is the constraint."
   → NEW "the third established that the constraint lies upstream of the tuning, in the retrieved
   per-pixel leaf area: a learner on the raw reflectance recovers part of the dense ranking (*R*² =
   0.35) but still less than LiDAR (0.60)."

9. **D.4** OLD "the leaf quantity that optical indices stop resolving beyond LAI 5–6;"
   → NEW "the leaf quantity whose between-plot ranking the optical retrieval has already lost
   throughout the dense stratum, well before its 6–7 m² m⁻² magnitude ceiling;"

10. **D.4** OLD "it is an excellent classifier, cheap to obtain wall-to-wall"
    → NEW "it is a serviceable classifier (85% agreement with the leaf-area split at Blois), cheap to
    obtain wall-to-wall"

11. **D.5** OLD "in which 53 independent temperature loggers judge LAI products by what they are for."
    → NEW "in which 53 temperature loggers, independent of both LAI products, judge them by what they
    are for."

12. **D.6** OLD "Chapter 3 treated the night as a sensor-agnostic regime that leaf area ranks only
    weakly; the correlation above suggests that part of that weakness may sit in the simulated field
    rather than in the coupling,"
    → NEW "Chapter 3 treated the night as a sensor-agnostic regime: the daytime reversal dissolves,
    both products rank the plots comparably (pooled ΔTmin *R*² 0.68–0.73), and the model runs
    +1.5–1.7 °C warm. The correlation above suggests that the missing piece is amplitude and bias
    rather than ranking,"

13. **D.7.1** OLD "Re-running the Chapter 3 machinery with the autumn limb of the seasonal series
    constrained by GEDI improved the simulated microclimate significantly"
    → NEW "Chapter 3 confined its logger validation to summer; the archive extends through autumn, and
    re-running the Chapter 3 machinery over October–November with the autumn limb of the seasonal
    series constrained by GEDI improved the simulated microclimate against those loggers significantly
    (47 series plots; ...)". Then align Ch3 App. A ("cover the summer only" → "were used over the
    summer window only in this chapter").

14. **D.7.1** OLD "ALS supplies the magnitude, Sentinel-2 the plot-scale summer timing, and GEDI the
    seasonal arbitration at the shoulders"
    → NEW "ALS supplies the magnitude and the dense-canopy ranking, Sentinel-2 the open-canopy ranking
    and the only plot-scale timing available outside the LiDAR date (of no value on the summer
    plateau, Chapter 3), and GEDI the seasonal arbitration at the shoulders".

15. **D.7.3** add "*n* = 31" after "later-leafing plots accumulate a warmer understory spring".

16. **D.8** OLD "a regime rule that says which sensor to trust on which half of the density gradient"
    → NEW "a domain of validity, a regime rule that says which sensor to trust on which half of the
    density gradient" (adopts Ch3's headline term once).

---

## (D) Chapter 3 points the Discussion should pick up (or contradicts)

1. **The summer-scoped temporal null (DYN vs STATIC).** Ch3 App. A/B: dynamic Sentinel-2 series
   changes nothing (0.56 vs 0.54; 0.24 vs 0.23), the plateau leaves nothing to resolve, and the
   temporal claim is "not testable beyond summer with these data". The Discussion never states this
   verdict, yet D.7.1 builds the whole temporal perspective on top of it. One sentence in D.1 (or
   D.7.1's opening) is needed: "Chapter 3 found no summer value in a resolved Sentinel-2 phenology,
   and could not test the shoulders."
2. **The pooled tie.** Fusion vs Sentinel-2 alone: Δ*R*² = 0.010, CI −0.11 to +0.20; the fusion's
   advantage is bias, RMSE and within-stratum robustness. Missing from D.1.
3. **The height ablation.** Optical LAI adds nothing on top of FORMS-H height (Δ*R*² 0.044, CI spans
   zero); naive S2 + FORMS-H forcing ranks the open plots at 0.73. This is the direct evidence for
   D.4's "height as classifier" and should be cited there.
4. **Multi-site RF transfer failure.** Leave-Blois-out learner pooled 0.29, dense 0.04: training
   elsewhere does not substitute for local LiDAR. Relevant to D.8's "cheap to test elsewhere" and to
   D.7.4's portability claim; currently absent.
5. **H2 raw-reflectance result** (dense *R*² = 0.35): a real Ch3 finding about *where* the optical
   information is lost (in the retrieval); D.3 should carry it (rewrite 8).
6. **Robustness set** (*k* = 0.5 vs 0.65, floor 2/3/5 m, model-free invariance to *k*): D.5 mentions
   only single-site and non-independence; one clause on *k*-invariance of the model-free ranking
   would strengthen the "rankings not effect sizes" stance and links to D.3's 23 % argument (which is
   a *k*-rescale of the same product: the ranking is invariant, the coupled skill is not, and the
   Discussion should say both).
7. **ERA5 macro reference and Terando warm-logger caveat** are named in Ch3 §4.3; D.6 cites Terando
   but not the ERA5 dependence of absolute ΔTmax and slope. Minor.
8. **Forest-type dependence of the saturation level and CLMS/HR-VPP inheritance** (Ch3 §4.3): fits
   D.7.4's "conifer and mixed stands are the stated boundary"; one clause.
9. **"Some form of LiDAR, airborne or spaceborne"** (Ch3 §4.3) vs D.7.1 "GEDI ... cannot replace ALS
   as a per-plot structural source": not a contradiction once the 150 m local-truth clause is read,
   but D.7.1 should acknowledge that Ch3 already framed the requirement as LiDAR-generic.
10. **Ch3-internal issue to surface to the author (not a Discussion defect):** Ch3 App. A's leaf-out
    sentence ("the parametric summer phenology matched or outperformed the observed Sentinel-2
    timing") conflicts with project memory of the 2×LAI runs (dynamic beats CONST/STATIC at leaf-out,
    ΔRMSE +0.29, significant). Verify which is current before the Discussion says anything about
    spring timing (D.7.3 currently avoids it, correctly).

---

## (E) Verdict

The Discussion's core Ch3 numbers (0.60/0.02, 0.71, 0.573/0.536, +0.67/+0.71, 1.60/1.62, 17.8 m, 85 %,
63 %, +1.14/+0.73, SDrec ≤ 0.28, +1.5–1.7 °C, *n* = 53, *r* = −0.92) all match the manuscript; the
retired phrases are gone. One severe cross-document contradiction remains: D.7.1 reports a significant
autumn logger-scored gain while Ch3 App. A tells the reader autumn cannot be validated with these
loggers (and the *n* may be 47, not 53). Below that, the Discussion consistently states Ch3's verdict
one notch stronger than Ch3 does: symmetric "resolution" where Ch3 says asymmetric and pooled-tied,
"worst of the set" where Ch3 says "no better than ATBD in domain", "excellent classifier" where Ch3
says "moderately", a located crossover where Ch3 says a median split stable from 2.0 to 4.5, and a
saturation line (5–6 / near 4) that is not Ch3's own 6–7 ceiling. Fix B1–B9 with the rewrites in (C)
and add the summer temporal null and the pooled tie (D1–D2); the arc then reads as Ch3 wrote it.
