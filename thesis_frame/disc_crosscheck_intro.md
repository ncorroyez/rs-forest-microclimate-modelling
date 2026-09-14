# Cross-check: General Discussion vs General Introduction

*Examiner pass, 2026-08-18. Target: `discussion_generale_these_EN.md` (302 lines, draft
2026-08-13, post harsh review). Reference: `intro_generale_these_EN.md` (932 lines) and its
locked plan. Rules: `these-manuscript` SKILL.md (locked framings), `style-gril-bouwen`.
Items already marked "fixed" in `discussion_generale_these_EN_harsh_review.md` (§1-§8, §9
D.8) are not repeated. Where a chapter fact was needed to adjudicate, Ch1's live manuscript
(`chapter1/manuscript/manuscript_chap1_EN_native20.md`) and Appendix G (`annexe_GEDI_EN.md`)
were consulted. Line numbers: I = intro, D = discussion.*

---

## (A) Intro promises vs Discussion status

Legend: **Closed** = explicitly answered; **Implicit** = answered but not tied back to the
promise; **Dropped** = never returns; **Contradicted** = the Discussion says otherwise.

| # | Intro promise (quote, line) | Discussion status | Where |
|---|---|---|---|
| A1 | RQ1: "Which dimension of canopy structure governs sub-canopy ΔTmax buffering, leaf quantity or vertical arrangement" (I743-745) | **Closed.** "Leaf quantity ... governs the realized ΔTmax buffering at both ends of the density gradient; the vertical arrangement ... is the largest single lever in mid-density canopy, but only as an upper bound" (D15-19). Matches the locked framing. | D.1 |
| A2 | RQ1 second half: "can a mechanistic model attribute it despite the collinearity" (I744-745) | **Closed** via D.5 ("mechanistic attribution ... physical parameters not fitted to the loggers, on-manifold perturbations", D135-137) and the collinearity guard in D.6 (D181-183). | D.5, D.6 |
| A3 | RQ1 hypothesis 1: "vertical structure modulates radiative and turbulent transfer enough to make the LAD profile a key parameter for a faithful microclimate simulation" (I755-757) | **Implicit.** D.1 and D.4 say the profile's realized effect is "below what *n* = 53 loggers resolve" and "statistically indistinguishable from zero" (D19-20, D114-115), i.e. H1 not realized. Never named as a hypothesis outcome. | D.1, D.4 |
| A4 | RQ1 hypothesis 2: "at fixed leaf area and height, foliage concentrated high in the canopy attenuates the diurnal cycle more strongly than foliage spread through the column" (I757-758); §1.7 adds "against which the realistic profile buffers slightly less" (I814) | **Dropped.** The Discussion never states the direction of the profile effect. Ch1 (l. 288) reports "H2 holds in direction over most of the density range, but as a second-order lever ... Its direction is not universal either". One clause is owed. | none |
| A5 | RQ1 design promise: answer "stratified by structural archetype, so that the answer is allowed to differ between open and dense canopy" (I752-754); §1.7 "Leaf quantity leads throughout ... As the canopy closes and leaf quantity saturates, its lead erodes until the vertical profile becomes a comparable-magnitude lever ... it co-leads in dense canopy because the alternatives saturate" (I808-811) | **Contradicted** (by the Discussion AND by Ch1). Discussion: leaf-quantity "levers are strongest in dense canopy (−0.26 °C per 0.5 LAI in the densest archetype)" (D69-71); profile "largest single lever in mid-density canopy" (D17-18). Ch1 §3.2 title: "the profile lever peaks in intermediate canopy"; l. 225 "it strengthens toward the closed end"; l. 229 "in the dense P4 it is dwarfed by leaf quantity (0.08 versus 0.26 °C)". The intro has the geometry inverted and uses "co-leads", a red-listed word. See B1. | D.1, D.2 |
| A6 | Ch1 validation promise: "reproduces the between-plot buffering ranking well, while carrying a warm bias ... compressing the micro-macro slope toward unity, so that both buffering and amplification are under-expressed" (I817-819); amplification set up at length (I72-89, I153) | **Partly closed.** Ranking/bias/amplitude triplet is closed (D20-22, D157-160). Amplification, the phenomenon the intro insists must be "established on observational grounds alone" (I88), never appears in the Discussion (0 hits). Ch1's own text is ambivalent (44/9 vs 45/8 recovered vs "does not reproduce the amplifying plots", see factcheck memory). The Discussion should say whether the model reproduces the amplifying tail. | none |
| A7 | Ch1 residual: "the residual traces chiefly to sub-pixel canopy gaps that a one-dimensional column cannot represent" (I820, again I927) | **Implicit/incomplete.** D.6 lists "edge effects, lateral advection, and cold-air drainage" (D164) and protocol (D165-167) but not sub-pixel gaps, which the intro (and Ch1) name as the chief structured residual. | D.6 |
| A8 | VPD: "Where a vapor-pressure-deficit offset is computed, it remains model-internal for want of a matching observation, and is reported as a robustness check rather than as a validated finding" (I99-101); "temperature buffering and humidity buffering are two readings of one partitioning" (I241) | **Dropped.** 0 hits for VPD/vapor in the Discussion. Already flagged as open in harsh review §9 (not marked fixed). | none |
| A9 | Bouwen positioning: "Three things consequently remain open. The measured profile has not been substituted for the assumed one. The effect of the vertical arrangement itself has therefore not been isolated ... And the work was conducted in evergreen conifer stands" (I420-425); Bouwen's beta function "judged to have contributed to an underestimation of the simulated buffering" (I416-417) | **Dropped as an explicit closure.** The Discussion never cites Bouwen (0 hits) and never says: the substitution was made, the isolated effect is small, and the under-buffering persists with measured profiles (D157-160), so the beta function was not its cause. This is the single cleanest positioning statement available and it is missing. | none |
| A10 | Bouwen's second settled point: "the meteorological forcing must itself be corrected for canopy structure before scenarios ... can legitimately be compared" (I411-413); §1.7 "the dense-canopy dominance of LiDAR is robust to the model's boundary-layer configuration, whereas the open-canopy usability of Sentinel-2 is the configuration-dependent link" (I857-859) | **Dropped.** No mention of the boundary-layer configuration robustness in the Discussion. The intro promises an asymmetric result; the Discussion states the reversal symmetrically (D40-43). | none |
| A11 | RQ2: "Can Sentinel-2 optical LAI be reconciled with LiDAR-derived structural LAI, or does optical saturation decouple them?" (I760-761) | **Closed.** "the divergence is one of dynamic range rather than of central tendency: regression slopes remain at 0.25–0.80 after every correction, while the site means are comparable at *k* = 0.65" (D29-31). | D.1 |
| A12 | RQ2 driver 1, effective depth (I764-765) | **Closed** (D32-34, "6–10 m of canopy (7 m across sites)"). Consistent with I826 "upper 6 to 10 m". | D.1 |
| A13 | RQ2 driver 2, PROSAIL parameterization; "bounds the extent to which the improvement depends on LiDAR-derived priors" (I836-837) | **Closed** (D89-91 "forest-tuned PROSAIL configuration that cuts RMSE by 1.2 m² m⁻²"; D145-146 "bounded prior circularity (≤ 17%)"). | D.3, D.5 |
| A14 | RQ2 driver 3, "horizontal heterogeneity of the canopy surface within the satellite pixel" (I766); §1.7 "Residual discrepancies remain largest in stands whose canopy surface is horizontally heterogeneous ... and a structural overestimation of leaf area in sparse stands persists that no parameter tuning removes" (I837-840) | **Dropped.** Heterogeneity and the sparse-stand overestimation are absent from the Discussion (0 hits for "heterogene", "sparse" in the Ch2 sense). One of the three announced drivers has no landing. | none |
| A15 | RQ2 output "carried forward as a condition of use" (I772); domain caveat "silent about short, open stands" (I686-687); §1.7 "Chapter 2 is therefore consumed as a map of where each sensor is valid, not as a corrector" (I872-873) | **Implicit.** D.3 says "the route, not the tuning, is the constraint" (D107-108) but omits the intro's first explanation, that the optimized product "was calibrated on closed, tall canopy, so applying it across the full density gradient extrapolates it into the open, short stands its calibration domain excluded" (I868-869). D.3 gives only the flattening mechanism. | D.3 |
| A16 | Intro caveat 2 of §1.5.2: "better agreement in leaf area does not automatically propagate into better simulated microclimate ... Which one a correction actually improves is an empirical question, not a corollary" (I689-697) | **Closed**, and well: "the General Introduction promised the statement as an open, empirical question, and the chapters answered it three times" (D84-86). | D.3 |
| A17 | RQ3: "Which remotely sensed description ... should force a mechanistic microclimate model, in which canopy regime" (I774-775); hypothesis "no single sensor wins everywhere, and ... the operationally useful product is therefore a rule for choosing between them" (I792-794) | **Closed** ("No single sensor wins everywhere ... the operational product is therefore a rule for choosing between them", D39-44). | D.1 |
| A18 | RQ3 second half: "can a product be built that remains deployable beyond the footprints where LiDAR exists?" (I775-776); §1.7 answer: "the forcing itself still draws canopy height and profile shape from LiDAR in both regimes" (I875-876); "The forcing still requires LiDAR structure wherever the canopy is dense, so coverage is bounded by the flight footprint" (I894-895) | **Implicit and slightly over-sold.** D.1 calls the height-gated form "operational" (D45-46) without saying that only the *gate* is LiDAR-free. D.7.1 says GEDI "cannot replace ALS as a per-plot structural source" (D195-196), and D.8 says nothing about the footprint bound. A referee reading D.1 alone will think the deployability question is answered "yes". | D.1, D.8 |
| A19 | RQ3 evaluation "on four independent axes ... a model-free test in leaf-area space ..., the between-plot ranking, the warm bias, and the recovery of observed thermal amplitude" (I782-784); §1.3.1 "this thesis reports the three separately throughout" (I386-387) | **Closed** but the count wobbles: D.5 says "ranking, dispersion, and bias separately throughout" (D150-151); the model-free axis surfaces only in D.3 (D93-94). Harmless, but say "four" once. | D.3, D.5 |
| A20 | RQ3 "the adjudication is stratified by canopy density at the leaf-area value where the two sensors cross over" (I787-788), value deliberately unnamed in the intro | **Closed** with the number, correctly attributed to Ch3: "the sensor crossover at LAI ≈ 3.9 (the stratified forcing comparison of Chapter 3)" (D57-58). Consistent with the intro's non-attribution to Ch2. | D.2 |
| A21 | Open-canopy verdict: "where the canopy is open, Sentinel-2 leads, though by a narrower margin on the thermal-coupling slope than on the daytime offset" (I854-855) | **Contradicted in emphasis.** D.1: "Sentinel-2 suffices where the canopy is open (0.71 against LiDAR near zero)" (D42-43). "Near zero" holds for the ΔTmax offset axis; on the coupling slope LiDAR ranks open plots (Ch3: 0.57 vs 0.74). See B4. | D.1 |
| A22 | Scope: "instrumented with 53 loggers over one summer" (I900); "outside summer the loggers do not exist, so the question cannot be arbitrated with these data at all" (I908-909) | **Contradicted.** The Discussion arbitrates autumn (D206-208), spring (D254-259), winter (D258-259), and June–July 2022 (D232-237) against the same logger network. Appendix G: "The HOBO network stayed active through the record June and July 2022 heatwaves". See B2. | D.7.1-D.7.3, D.8 |
| A23 | Scope: "within summer a temporally resolved satellite series adds nothing ... there is no timing left to resolve" (I906-908) | **Contradicted in wording.** D.7.1: "ALS supplies the magnitude, Sentinel-2 the plot-scale summer timing" (D216-217). See B3. | D.7.1 |
| A24 | Scope: thresholds "calibrated in-sample and will require recalibration elsewhere; because the loggers are not spatially independent, the confidence intervals Chapter 3 derives from them are optimistic" (I900-903) | **Closed** (D147-148 "the arbitration is single-site, and its non-independent loggers make the reported *R*² optimistic"), but see B5 for the internal contradiction with D143. | D.5 |
| A25 | 1-D ceiling: "This sets a ceiling on absolute accuracy that is independent of how well the structural inputs are measured" (I466-467) | **Closed** and developed ("the margin of progress has shifted from the forcing to the model", D162-163). | D.6 |
| A26 | Perspective 1, GEDI: "Sampling geometry and geolocation uncertainty nonetheless keep it short, for now, of what plot-scale mechanistic microclimate modeling requires" (I921-922); Fig. 1.4 caption "GEDI is drawn hollow: it is a perspective taken up in the General Discussion" (I736) | **Closed** and quoted back faithfully (D191-193). | D.7.1 |
| A27 | Perspective 2, 3-D solver: "heterogeneity- and edge-aware formulations that couple a canopy scheme to a three-dimensional atmospheric solver ... relax exactly the assumption this thesis has had to accept" (I927-932) | **Closed**, near-verbatim (D270-273). See (B) redundancy R1. | D.7.4 |
| A28 | Boundary named first: "This thesis does not produce a wall-to-wall microclimate map. It establishes what a map would have to be built on" (I890-891) | **Closed**, near-verbatim (D291-293). Deliberate echo; acceptable once, but see R2. | D.8 |
| A29 | "Both halves of that question matter, and both are addressed in what follows": how hot on the worst days, and where (I39-44) | **Half-closed.** "How hot on the worst days" lands in D.6/D.7.2 (hottest days, heatwaves). "Where" is conceded as unmapped (D291). Fine, but the Discussion could say so in one clause. | D.6, D.8 |
| A30 | Microrefugia stake (I135-147) | **Closed** (D228-231). | D.7.2 |
| A31 | Roadmap Fig. 1.5: "each has to be closed before the next can be posed" (I798) | **Closed** ("The chapters bear the sequencing out", D11-13). | D.1 |

---

## (B) Contradictions and mismatches, ranked

### Severe

**B1. Where the profile lever peaks: intro says dense, Discussion (and Ch1) say intermediate; intro uses the red-listed "co-leads".**
- Intro I808-811: "As the canopy closes and leaf quantity saturates, its lead erodes until the vertical profile becomes a comparable-magnitude lever, having already overtaken cover and height. The vertical dimension is a small absolute lever everywhere; it co-leads in dense canopy because the alternatives saturate, not because it is large." Also I882-883: "with the vertical arrangement mattering only once quantity has saturated (Chapter 1)".
- Discussion D15-18: "Leaf quantity ... governs the realized ΔTmax buffering at both ends of the density gradient; the vertical arrangement of foliage is the largest single lever in mid-density canopy, but only as an upper bound". D69-71: "Its mechanistic levers are strongest in dense canopy (−0.26 °C per 0.5 LAI in the densest archetype)".
- Ch1 native20 l. 225 "Leaf quantity is the dominant lever at both ends of the density gradient, and it strengthens toward the closed end"; l. 229 "in the dense P4 it is dwarfed by leaf quantity (0.08 versus 0.26 °C) ... the swap is the largest single lever in P3".
- Skill lock: "largest lever in intermediate canopy (P2/P3) ... Never 'co-lead', 'small-but-leading', or 'leads in dense canopy'".
- Verdict: the Discussion is right, the intro is wrong on the geometry (leaf quantity does not saturate as a lever; it steepens) and violates the locked vocabulary. A jury member reading Ch1 §3.2 after the intro will catch it. Fix the INTRO (see C-intro 1).

**B2. Logger availability outside summer: the intro forbids what the Discussion does.**
- Intro I899-909: "The microclimate work rests on a single lowland oak site instrumented with 53 loggers over one summer ... outside summer the loggers do not exist, so the question cannot be arbitrated with these data at all."
- Discussion D206-208 (autumn ΔRMSE, ΔR² against loggers), D232-237 ("The logger network stayed active through the record June and July 2022 heatwaves"), D254-259 (spring temperatures, "winter microclimate does not predict green-up across plots"), D299-301 ("one out-of-year (the 2022 heatwaves)").
- Appendix G §G.6: "The HOBO network stayed active through the record June and July 2022 heatwaves ... n = 51 plots"; §G.4 autumn Oct–Nov scored on 53 loggers.
- Verdict: a flat factual contradiction. The intro sentence was written to justify the summer scoping of Chapter 3; it should scope the *chapter*, not deny the archive. Fix the INTRO (C-intro 2). The Discussion may add half a clause acknowledging that the chapters used the summer window and the perspectives use the rest of the archive.

### Moderate

**B3. "Sentinel-2 the plot-scale summer timing" vs "within summer a temporally resolved satellite series adds nothing".**
- Discussion D215-217: "ALS supplies the magnitude, Sentinel-2 the plot-scale summer timing, and GEDI the seasonal arbitration at the shoulders".
- Intro I906-908: "within summer a temporally resolved satellite series adds nothing, and does so for a mechanical reason, since leaf area is on its phenological plateau and there is no timing left to resolve".
- Appendix G Table G3: ALS+S2 (temporal, S2-timed) summer 0.24/1.76/0.73 = ALS static 0.24/1.76/0.73, i.e. the intro is right for the JJAS plateau; the S2 timing role Appendix G defends is that "a site-level GEDI shape cannot replace the plot-scale Sentinel-2 timing (summer R² 0.15)", a leaf-on-trajectory role, not a within-plateau gain.
- Verdict: the Discussion word "summer" is misleading; the role is per-plot leaf-on timing (spring/shoulders, and the summer *shape* against a site-level curve). Rephrase (C1).

**B4. Open-canopy LiDAR "near zero" vs intro's "narrower margin on the coupling slope".**
- Discussion D41-43: "Sentinel-2 suffices where the canopy is open (0.71 against LiDAR near zero)".
- Intro I853-855: "where the canopy is open, Sentinel-2 leads, though by a narrower margin on the thermal-coupling slope than on the daytime offset".
- Ch3 (factcheck memory): on the coupling slope in open canopy LiDAR 0.57 vs S2 0.74; "il ne cesse pas de classer". "Near zero" is the ΔTmax-offset axis only.
- Verdict: name the axis in D.1 (C2). Also reconcile with the intro's asymmetry claim (boundary-layer robustness, A10): the intro promises "the open-canopy usability of Sentinel-2 is the configuration-dependent link"; the Discussion should carry that qualifier or the intro should drop it.

**B5. Internal contradiction inside D.5 on logger independence (and against the intro).**
- D143: "53 independent temperature loggers judge LAI products".
- D147-148: "its non-independent loggers make the reported *R*² optimistic". Intro I902: "the loggers are not spatially independent".
- Verdict: "independent" at D143 must mean "independent of the model and of the LAI retrievals", not spatially independent. Say so (C3).

**B6. Bouwen positioning never closed (A9-A10).**
- Intro I409-426 builds Ch1 as the answer to Bouwen's closing recommendation and lists three open points; the Discussion never returns to Bouwen, never says the beta-function hypothesis for under-buffering is not supported (the warm bias survives measured profiles, D157-160), and never mentions the forcing-correction/boundary-layer robustness. Ch1 itself does the comparison (native20 l. 298 "The nearest comparable study sharpens this rather than contradicting it. @bouwenInteractionsEntreStructure imposed a generic beta distribution ..."). Add two sentences to D.4 or D.6 (C4).

**B7. RQ2 third driver dropped (A14).** Horizontal heterogeneity and the sparse-stand overestimation, both promised in RQ2 and reported in §1.7, do not land anywhere. One sentence in D.1's Chapter 2 paragraph (C5).

**B8. Deployability half of RQ3 (A18).** D.1's "operational form, gated on a canopy-height product" (D45-46) needs the intro's qualifier that height and profile shape still come from LiDAR in both regimes, otherwise D.1 promises more than I875-876 and I894-895 concede. Add one clause (C6).

**B9. "flattening near 4" in D.2 (D59-60) vs the intro's and Ch2's saturation numbers.**
- Intro I601-602: "largely exhausted somewhere above five". Ch2 alignment decision (factcheck memory, 2026-08-04): Ch2 quotes literature 5–6, an observed ceiling of 6–7, and a 4–5 "threshold"; "jamais « 4 »". The Discussion's "our own pooled Sentinel-2 response flattening near 4" comes from the Ch3 pixel pool (3.85/4.20 at *k* = 0.5), which harsh review §"could not verify" already flagged. Either attribute it to the Chapter 3 pool at *k* = 0.5 or drop the number (C7). Note that D.2 already mixes *k* conventions across its three lines and says so; adding a fourth unattributed number weakens the paragraph the harsh review just repaired.

**B10. Title of D.3 still says "Retrieval accuracy" (D77).** Harsh review §3 fixed the body but the heading keeps the word the intro spends I673-678 and I769 forbidding ("the target cannot be absolute accuracy"). Retitle (C8). The intro has the same slip once: I864-866 "accuracy in leaf area and skill in microclimate come apart. The forest-tuned optical retrieval ... is the better leaf-area product" (see C-intro 3).

**B11. Amplification (A6) and sub-pixel gaps (A7) not carried into D.6.** The intro invests three paragraphs (I72-89, I153, I819) in amplification as the phenomenon a model must not lose and names sub-pixel gaps as the chief residual; D.6 mentions neither word. Two clauses (C9).

### Minor

**B12. VPD robustness check (A8)** never mentioned; one clause in D.5 (C10) or accept knowingly.

**B13. Hypothesis 2 direction (A4)** never stated; one clause in D.4 (C11).

**B14. Terminology: four names for the micro-macro slope.** Intro: "slope and an equilibrium term" (I346), "micro-macro slope" (I818), "thermal-coupling slope" (I855). Discussion: "buffering slope" (D57, D236), "slope-LAI coupling" (D286). Skill: one term per concept. Pick "buffering slope (the micro-macro slope of Gril et al., 2023)" once and reuse. Also "sub-canopy" (intro-dominant) vs "understory" (Discussion-dominant): the skill locks "understory microclimate"; the intro's "sub-canopy" is defensible as the physical descriptor but the two texts should not each own one.

**B15. *k* = 0.65 appears undefined in the Discussion** (D31, D88). The intro never introduces the extinction coefficient by symbol (I524-525 "an extinction coefficient ... commonly fixed at a nominal value"). Gloss at first use in D.1: "the LiDAR extinction coefficient *k* = 0.65 adopted for inter-sensor work (0.5 in Chapter 1)". As it stands, "the site means are comparable at *k* = 0.65" is opaque to a reader coming from the intro.

**B16. Ch1 "field tests" (D19-20) not announced by the intro.** The intro (I815-817) says validation "establishes the model's credibility rather than testing the per-plot sensitivities directly"; the Discussion invokes "the field tests [that] place its realized effect below what *n* = 53 loggers resolve", i.e. the model-free forward-inclusion. Not a contradiction, but the intro could announce that observational check in one clause so the Discussion's "field tests" has an antecedent (C-intro 4).

**B17. Four axes vs three (A19).** Trivial; say "four axes (one model-free)" once in D.5.

**B18. Em-dash.** Only inside the shared placeholder "[REF: Grulois 2026 — coupled ARPS-MuSICA thesis]" (I929, D271-272). Will vanish when the key is resolved; make sure the replacement key does not keep the dash.

**B19. Skill/voice compliance otherwise clean.** US spelling verified (modeled, modeling, recenters, favor); no Shapley, no "co-lead" in the Discussion, italics on *r*, *R*², *n*; en-dashes for ranges; Bouwen problem-then-resolution architecture throughout. Two Gril-flavored flourishes survive ("the optical signal has already left the field", D67; "the clearest unwritten chapter of this thesis", D185); the harsh review kept them, and they are defensible, but D185 is the more poetic of the two and Nathan's stated preference is sober. Optional.

### Redundancy with the intro (near-verbatim)

- **R1.** I927-932 "heterogeneity- and edge-aware formulations that couple a canopy scheme to a three-dimensional atmospheric solver [REF: Grulois 2026; @vandewalleForEdgeClimV103D2026] ... Those formulations relax exactly the assumption this thesis has had to accept." ↔ D270-273 "Heterogeneity- and edge-aware formulations that couple a canopy scheme to a three-dimensional atmospheric solver ([REF: Grulois 2026 ...]; @vandewalleForEdgeClimV103D2026) relax exactly the assumption our one-dimensional column had to accept". Same sentence, same citations. Rewrite one side (C12).
- **R2.** I890-892 "This thesis does not produce a wall-to-wall microclimate map. It establishes what a map would have to be built on: which structural dimension the temperature responds to, what a satellite can and cannot see of that dimension, and which product should drive the model in which canopy regime." ↔ D291-296. Deliberate bookend; acceptable, but D.8's list already differs (adds seasonal arbitration and model ceiling), so lean on the difference rather than repeating the first two clauses.
- **R3.** I924-926 "A one-dimensional column omits edges and lateral advection, yet forest edges are systematically warmer than interiors [@reekForestEdgesAre2025]" ↔ D163-165 "The one-dimensional column omits edge effects, lateral advection, and cold-air drainage [@reekForestEdgesAre2025]". Mild.
- **R4.** I376-387 (three metrics reported separately) ↔ D149-151. Mild, and D.5 adds the reason ("the three axes disagree too often"), so keep.

---

## (C) Suggested exact rewrites

### C. Discussion (old → new)

**C1 (B3), D215-218.**
OLD: "ALS supplies the magnitude, Sentinel-2 the plot-scale summer timing, and GEDI the seasonal arbitration at the shoulders, where optical senescence runs early"
NEW: "ALS supplies the magnitude, Sentinel-2 the per-plot leaf-on timing (a role a site-level GEDI shape cannot take over, although within the summer plateau that timing adds little, as Chapter 3 found), and GEDI the seasonal arbitration at the shoulders, where optical senescence runs early"

**C2 (B4), D41-43.**
OLD: "LiDAR magnitude is indispensable exactly where buffering is deepest (dense-canopy ΔTmax *R*² = 0.60 against 0.02 for Sentinel-2), Sentinel-2 suffices where the canopy is open (0.71 against LiDAR near zero), and the operational product is therefore a rule for choosing between them."
NEW: "LiDAR magnitude is indispensable exactly where buffering is deepest (dense-canopy ΔTmax *R*² = 0.60 against 0.02 for Sentinel-2), Sentinel-2 suffices where the canopy is open (0.71 against LiDAR near zero on the daytime offset, by a narrower margin on the buffering slope), and the operational product is therefore a rule for choosing between them. The two halves of that reversal do not carry equal weight: the dense-canopy dominance of LiDAR survives the model's boundary-layer configuration, whereas the open-canopy usability of Sentinel-2 depends on it."

**C3 (B5), D143.**
OLD: "in which 53 independent temperature loggers judge LAI products by what they are for."
NEW: "in which 53 temperature loggers, independent of both the model and the LAI retrievals, judge LAI products by what they are for."

**C4 (B6), append to D.4 after D128 (or to D.6 after D162).**
NEW: "The distinction also closes the question the General Introduction inherited from @bouwenInteractionsEntreStructure. That study attributed part of its under-buffering to the parametric beta profile it had to assume and recommended substituting measured LiDAR profiles. We made the substitution, in a deciduous canopy, and isolated the profile at matched leaf area and height: its realized contribution is small, and the warm bias and compressed amplitude persist unchanged (§D.6). The assumed profile was not the cause of the under-buffering; the one-dimensional column is."

**C5 (B7), D.1 Chapter 2 paragraph, after D34 "against canopy tops of 17–32 m."**
NEW: "The third driver, horizontal heterogeneity of the canopy surface within the pixel, is where the residual discrepancy concentrates after depth and parameterization are accounted for, and a structural overestimation of leaf area in sparse stands survives every correction; both mark the edge of the domain over which consistency was established."

**C6 (B8), D45-48.**
OLD: "its operational form, gated on a canopy-height product at 17.8 m, reproduces the assignment for 85% of plots and concedes little"
NEW: "its operational form, in which only the gate is LiDAR-free (a canopy-height product at 17.8 m) while canopy height and profile shape still come from ALS in both regimes, reproduces the assignment for 85% of plots and concedes little"
And in D.8, after D294-295 "measured without saturation where buffering is deepest;": add "which still means airborne LiDAR wherever the canopy is dense, so coverage remains bounded by the flight footprint;".

**C7 (B9), D58-60.**
OLD: "placed at LAI 5–6 in the literature, with our own pooled Sentinel-2 response flattening near 4."
NEW: "placed at LAI 5–6 in the literature and at an observed ceiling of 6–7 in Chapter 2's pixel pool (at *k* = 0.5, the pooled Sentinel-2 response of the Chapter 3 site already flattens near 4)."
(or simply delete "with our own pooled Sentinel-2 response flattening near 4" if the 6–7 ceiling is not to be quoted here).

**C8 (B10), D77.**
OLD: "## D.3 Retrieval accuracy is not functional fidelity"
NEW: "## D.3 Inter-sensor consistency is not functional fidelity"
and D79-81 keep the Ch1 quotation as attributed, as already done.

**C9 (B11), D163-165.**
OLD: "The one-dimensional column omits edge effects, lateral advection, and cold-air drainage [@reekForestEdgesAre2025], and part of the measured gap is protocol rather than physics"
NEW: "The one-dimensional column omits the sub-pixel gaps that structure the Chapter 1 residual, together with edge effects, lateral advection, and cold-air drainage [@reekForestEdgesAre2025]; the amplifying tail of the gradient, the open plots that run warmer than the open air, is where this omission bites hardest, and part of the measured gap there is protocol rather than physics"

**C10 (B12), D.5 after D151.**
NEW: "The same discipline applies to the vapor-pressure-deficit offsets that Chapter 1 computed: they are model-internal, unvalidated by the loggers, and were carried only as a robustness check on the temperature result."

**C11 (B13), D.4 after D116 "(≤ 0.004 °C m⁻¹)."**
NEW: "The direction the General Introduction hypothesized, that foliage concentrated high in the canopy buffers more, holds over most of the density range but as a second-order lever, and reverses at the dense end."

**C12 (R1), D270-273.**
OLD: "Heterogeneity- and edge-aware formulations that couple a canopy scheme to a three-dimensional atmospheric solver ([REF: Grulois 2026 — coupled ARPS-MuSICA thesis]; @vandewalleForEdgeClimV103D2026) relax exactly the assumption our one-dimensional column had to accept, and §D.6 argues the next margin of progress lies there."
NEW: "The first is dimensional. §D.6 locates the next margin of progress in the column itself, and the coupled canopy-atmosphere formulations named in the General Introduction ([REF: Grulois 2026 coupled ARPS-MuSICA thesis]; @vandewalleForEdgeClimV103D2026) are where the regime rule of Chapter 3 should next be evaluated, in the edges and gaps a column cannot see."

**C13 (B15), D30-31.**
OLD: "while the site means are comparable at *k* = 0.65."
NEW: "while the site means are comparable once the LiDAR extinction coefficient is set to *k* = 0.65, the value adopted for inter-sensor work (Chapter 1 keeps *k* = 0.5)."

**C14 (B17), D150.** "ranking, dispersion, and bias separately throughout" → "ranking, dispersion, and bias separately throughout, with a model-free leaf-area test as a fourth axis in Chapter 3".

### C. INTRO (do not edit here; for the intro pass)

**C-intro 1 (B1), I808-811.**
OLD: "Leaf quantity leads throughout, by a wide margin in open canopy. As the canopy closes and leaf quantity saturates, its lead erodes until the vertical profile becomes a comparable-magnitude lever, having already overtaken cover and height. The vertical dimension is a small absolute lever everywhere; it co-leads in dense canopy because the alternatives saturate, not because it is large."
NEW: "Leaf quantity leads at both ends of the density gradient, and its lever steepens as the canopy closes. The vertical profile is a potent lever only when fully exercised: reported as a complete real-versus-uniform swap, it becomes the largest single lever in intermediate canopy, an upper bound that real canopies barely span, and in the densest canopy it is again dwarfed by leaf quantity."
And I882-883: "with the vertical arrangement mattering only once quantity has saturated (Chapter 1)" → "with the vertical arrangement a bounded, second-order lever that peaks in intermediate canopy (Chapter 1)".

**C-intro 2 (B2), I899-909.**
OLD: "instrumented with 53 loggers over one summer ... outside summer the loggers do not exist, so the question cannot be arbitrated with these data at all."
NEW: "instrumented with 53 loggers, of which the chapters use the summer 2021 window ... outside summer the chapters do not arbitrate the question; the logger archive extends through the following seasons and the 2022 heatwaves, and the General Discussion draws on it for perspectives only."

**C-intro 3 (B10), I864-866.**
OLD: "accuracy in leaf area and skill in microclimate come apart. The forest-tuned optical retrieval optimized in Chapter 2 is the better leaf-area product"
NEW: "consistency in leaf area and skill in microclimate come apart. The forest-tuned optical retrieval optimized in Chapter 2 is the product most consistent with LiDAR"

**C-intro 4 (B16), I815-817.** After "since the perturbation design points are not the logger plots." add: "A model-free check on the same loggers asks the observational form of the question, whether adding the measured profile to leaf quantity and cover improves the fit to the observed buffering."

**C-intro 5 (B14).** Choose one name for the slope in I346/I818/I855 ("buffering slope", introduced at I346 as "the slope-and-equilibrium formulation") so that D57/D236/D286 inherit it.

**C-intro 6 (B3).** If C1 is applied to the Discussion, I906-908 can stand; otherwise soften "adds nothing" to "adds nothing to the ranking of plots within the summer plateau".

---

## (D) Verdict

1. The Discussion closes every research question and both stated hypotheses in substance, and its D.2 gradient logic and D.5 epistemic stance are the intro's own (I723-729, I462-469, I673-697) carried through; the arc holds.
2. Two severe mismatches, both faults of the INTRO, not the Discussion: the intro places the profile's "co-lead" in dense canopy (Ch1 and the Discussion say intermediate, and "co-lead" is red-listed), and it states that loggers "do not exist" outside summer while the Discussion arbitrates autumn, spring, and 2022 against them.
3. Four moderate omissions in the Discussion: no return to Bouwen (the cleanest positioning statement in the thesis is left unsaid), no landing for RQ2's third driver (heterogeneity/sparse-stand overestimation), no boundary-layer-robustness asymmetry (promised at I857-859), and an over-sold "operational" fusion whose forcing still needs ALS.
4. Wording fixes: "summer timing" (D217) vs the intro's summer negative; "independent loggers" (D143) vs "non-independent" (D148); "Retrieval accuracy" heading (D77); undefined *k* = 0.65; "flattening near 4" needs its provenance or removal.
5. Style is compliant (US, no em-dash outside the shared REF placeholder, locked vocabulary respected in the Discussion); one near-verbatim repeat of the intro (3-D solver perspective) should be rewritten on one side.
