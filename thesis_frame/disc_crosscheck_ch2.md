# Cross-check — General Discussion vs Chapter 2 (RSE-D-25-04417, r2 + r3 inserts)

*Examiner pass, 2026-08-18. Target: `thesis_frame/discussion_generale_these_EN.md` (draft 2026-08-13,
post harsh-review fixes). Reference: Ch2 r2 plain text, r3 inserts (`Manuscript_Inserts_r3`),
`Ajouts_encadrants_r3` (retitle + §5c d_opt-circularity paragraph), `Tables_r3/*.csv`,
`reviewers_comments.md`. Items already marked fixed in the harsh review are not repeated
unless the fix left a residue.*

Two facts about the Ch2 end-state frame everything below:

- **Retitle (Ajouts_encadrants r3, block 0):** "Discrepancies Between S2- and LiDAR-derived LAI in
  Temperate Deciduous Forests linked to effective optical depth, PROSAIL LUT priors and
  intra-pixel canopy heterogeneity". Rationale given: move away from "Assessing and Reducing"
  (accuracy claim) to a **decomposition into three named components**. Jérôme's global
  instruction: chase residual "improvement / improved / gain" wording.
- **New circularity taxonomy (r3 A.16 + §5c):** (i) prior circularity, bounded ≤ 17 % of the RMSE
  reduction; (ii) scoring circularity, unbounded; (iii) **d_opt-selection circularity**, also
  reference-driven and unbounded ("not physical penetration depths").

---

## (A) Claim-by-claim table

| # | Discussion (section, quote) | Ch2 (location, quote) | Verdict |
|---|---|---|---|
| 1 | D.1 "Chapter 2 answered the consistency question." | Title r3: "Discrepancies … linked to effective optical depth, PROSAIL LUT priors and intra-pixel canopy heterogeneity"; §4.3 "reflects improved inter-sensor consistency rather than demonstrated accuracy gains"; §5 "the analysis quantifies inter-sensor consistency". | **OK on positioning**, but Ch2's final self-description is a *decomposition into three components*; the Discussion names only two (depth, LUT) and never the third (heterogeneity). See B-severe-2, D-1. |
| 2 | D.1 "Sentinel-2 LAI saturates against the LiDAR reference in dense canopy, and the divergence is one of dynamic range rather than of central tendency" | §3.1 "LAI_S2_ATBD saturated at a maximum value of 6–7 m²/m² … LAI_ALS reached up to 10 … and even 15"; §4.2 "LAI_S2 also tends to overestimate LAI_ALS below 4–5 m²/m²"; §5 "The persistence of slopes well below unity … is driven **primarily** by the amplified gap … at low LAI_ALS values: where LAI_ALS falls below ~2 m²/m², LAI_S2 remains substantially higher". | **Half of Ch2.** Ch2's Conclusion attributes the sub-unity slope *primarily* to low-LAI overshoot, not to dense-canopy saturation. "Dynamic range" is a fair gloss only if both ends are named. See B-moderate-1. |
| 3 | D.1 "regression slopes remain at 0.25–0.80 after every correction" | §4.3 "regression slopes remain well below unity (0.25–0.80 in Pareto-optimal configurations)". Table 6 \|Slope−1\| ranges imply slopes 0.15–0.47 (Aigoual), 0.49–0.97 (Blois), 0.42–0.71 (Mormal); Fig. 7 final slopes 0.49 / 0.51 / 0.62. | Quote is faithful to Ch2 text; **Ch2's own text and Table 6 disagree** (0.15–0.97 from the table). Not a Discussion error, but the number is fragile: quote the retained-configuration slopes (0.49–0.62) or fix Ch2. Minor. |
| 4 | D.1 "the site means are comparable at *k* = 0.65" | §3.1 "k was set to 0.65, providing a joint compromise across the three sites between RMSE minimisation and **bias re-centring**"; Fig. 7 baseline bias Aigoual −0.32, Blois +1.11, Mormal +0.03 m²/m²; §4.2 insert: Blois full pool LAI_S2_ATBD mean 4.65 vs LAI_ALS 3.44. | **Partly by construction, and false at Blois.** k = 0.65 was chosen in part to re-center the bias; and Blois shows a +1.1–1.2 m²/m² S2 excess. Say "at Aigoual and Mormal, with k chosen partly to re-center the bias". See B-moderate-2. |
| 5 | D.1 "Sentinel-2 reflectance is most informative about the upper 6–10 m of canopy (7 m across sites), a thickness we termed the effective optical depth d_opt, against canopy tops of 17–32 m" | Abstract "restricted to the upper 6–10 m"; Table 3 (r3) "Effective optical depth (Pareto selection): Aigoual 8 m; Blois 6 m; Mormal 10 m; common 7 m"; §4.2 "canopy thickness within which S2 reflectance is most informative"; "(17–32 m in dominant stands)". | **Numbers correct; term matches r3 title/Table 3.** But Ch2 body prose says "empirical proxy … not a direct measurement of optical penetration depth", the range holds only under CHM normalisation (DTM: 10 m at all sites, Table 5), and Blois 6 m is "a statistical compromise across weakly-differentiated criteria". The Discussion carries none of these three qualifiers. See B-moderate-3. |
| 6 | D.1 "The chapter improved inter-sensor consistency, never accuracy" | §4.3 verbatim spirit; r3 instruction: replace "improved/improvement" wording by "decomposition / reduced inter-sensor RMSE relative to ATBD". | Positioning correct; **verb "improved" is on the r3 red list**. Cosmetic but visible. |
| 7 | D.1 "it flagged the caveat that mattered most for what followed: better agreement in leaf area does not automatically propagate into better simulated microclimate" ; "Chapter 3 answered the forcing question, and the caveat became the result." | §4.5: "the optimised LAI_S2_opt product **may serve as an operational input** for downstream applications … Forest microclimate modelling … is one such application: coupling LAI_S2_opt with canopy-height products … (e.g., FORMS-H) could provide the structural variables needed to parameterise radiation transfer models". No sentence in Ch2 (r2, r3 inserts, Ajouts) states that agreement may not propagate to microclimate skill. | **Not in Ch2.** Ch2 *offers* the product for microclimate modelling; it does not caveat it. D.3 itself says so ("a sentence Chapter 2 offered in good faith … Chapter 3 tested that offer and declined it"), so D.1 and D.3 contradict each other about what Ch2 said. See B-severe-1. |
| 8 | D.1 "a uniform 23% error in leaf area erases the model's buffering skill" (Ch1) ; D.3 "a uniform 23% reduction of the LiDAR leaf area, the direction and order of the *k* = 0.65 harmonization adopted for inter-sensor work" | §2.5.3 "k … generally set to 0.5"; §4.1 "We retained k = 0.65 as the operational value, which falls close to the field-reported range … (k ≈ 0.66–0.81) … Because LAI_ALS scales multiplicatively with k, changing k globally rescales LAI_ALS without altering its linear relationship with LAI_S2". 0.5/0.65 = 0.77 → −23 %. | Arithmetic correct. Two wording issues: (a) "error" in D.1 implicitly makes k = 0.5 the truth, whereas Ch2 defends 0.65 from field-reported k; (b) "harmonization" is not Ch2's word (Ch2: "operational value", "joint compromise … RMSE minimisation and bias re-centring"). Ch3 also uses 0.65 (Ch3 §2). Minor, but say "rescaling" and name the k conventions once (Ch1 0.5; Ch2 and Ch3 0.65). |
| 9 | D.2 "the optical saturation threshold, placed at LAI 5–6 in the literature, with our own pooled Sentinel-2 response flattening near 4" | §1 "S2-derived LAI is known to saturate above 5–6 m²/m² (Tesfaye et al., 2020)"; §3.1 "saturated at a maximum value of **6–7** m²/m² across all three state forests"; §4.2/§4.4 divergence threshold "4–5". Memory note: "near 4" (3.85/4.20) is the Ch3 canonical pixel pool, at k = 0.5. | **Provenance mismatch, still open** (harsh review left it as "check before quoting"). Ch2's own observed ceiling is 6–7, not 4. If "near 4" is kept it must be tagged Ch3 (Blois plot pool) and its k-scale stated; otherwise use Ch2's 6–7 ceiling. See B-moderate-4. |
| 10 | D.2 "the sensor crossover at LAI ≈ 3.9 (the stratified forcing comparison of Chapter 3)" | Absent from Ch2 (memory: 3.86 not in RSE manuscript; Ch2's only LAI threshold is the qualitative 4–5). Ch3 standalone: "the 3.86 split is a convenience placed at the sample median". | **Provenance now correct** (Ch3, not Ch2). Residual: Ch3 itself calls it both "sensor crossover" and "convenience placed at the sample median"; D.2's "crossover" inherits the stronger reading. Minor; align with Ch3's final wording. |
| 11 | D.2 "In dense stands, 63% of the leaf area sits below the effective optical depth" | Absent from Ch2 (grep "63" → only unrelated hits). Source is Ch3 (Blois dense plots, PAD profile: 63 % vs 19 % open). | **Unattributed Ch3 number placed inside Ch2 material.** Add "(Chapter 3, Blois plots)". Minor-moderate. |
| 12 | D.3 "Chapter 2 produced the best retrieval of the thesis, a forest-tuned PROSAIL configuration that cuts RMSE by 1.2 m² m⁻² and recenters the bias" | Abstract "reduced the reported inter-sensor RMSE values in LAI by 1.2 m²/m² and the bias by 1.5 m²/m², on average, across sites and relative to the ATBD baseline"; Table A.16: RMSE 1.49→1.32 (Aig.), 3.23→0.83 (Blois), 1.65→0.71 (Mormal), all **against LAI_ALS_dopt**; Blois residual bias +0.78 (A.16) / +0.67 (Fig. 7); optimisation "run separately for each study site". | Numbers correct. Three drifts: "best retrieval" is accuracy language (retired); "forest-tuned" hides that the LUT is site-specific and that at Blois its LAI prior is LAI_ALS_dopt (Table A.16); "recenters" is true at Aigoual/Mormal, partial at Blois. Crucially, the 1.2 m² m⁻² is agreement with a **depth-truncated** LiDAR LAI. See B-severe-3. |
| 13 | D.3 "a sentence Chapter 2 offered in good faith, that the optimized product 'may serve as an operational input' for microclimate modeling" | §4.5 exact phrase. | **Correct.** Note Ch2 also proposes coupling with FORMS-H canopy height, i.e. it anticipates the ingredient Ch3 turns into the gate. See D-4. |
| 14 | D.3 "Truncating the LiDAR magnitude to the optical depth fails the same way (+1.14 °C against +0.73 °C)" | Ch3 number; Ch2 supplies the mechanism (§4.2 "deeper layers … contributing marginally to the S2 signal despite hosting substantial leaf area"). | Correct; could be sharpened by the A.16 fact that the Blois-optimal LUT prior *is* LAI_ALS_dopt, so the two failures share a cause by construction. |
| 15 | D.5 "the consistency analysis carries a bounded prior circularity (≤ 17% of the reported improvement) and an unbounded scoring circularity" | A.16 caption "at most ~17 % of the reported **RMSE reduction** … residual scoring circularity … not bounded here"; Ajouts §5c and A.16 phrase 2: "the d_opt selection … introduce[s] additional forms of circularity that this analysis does not address". | Two of three named. Add the **d_opt-selection circularity** and say "RMSE reduction" not "improvement". See B-moderate-5. |
| 16 | D.5 "No field LAI was measured in this thesis, at any site, in any chapter." | §2.1 "Complementary data collected at the same sites include hemispherical photographs … not directly used in the present analysis"; §4.1 "absence of reliable ground-based LAI measurements that are concomitant to the airborne LiDAR and Sentinel-2 acquisitions". | **Over-strong.** Ch2 says no *coincident* field LAI was used; DHP exists in the site network. A jury member who knows IMPRINT will ask. See B-moderate-6. |
| 17 | D.5 "the reference instruments themselves disagree substantially in closed forests [Woodgate 2015]" | §4.1 (r2 + r3 block 3): Woodgate "across diverse forest systems"; the closed-canopy inference is Ch2's own extension. | Slight over-specification; "across diverse forest systems" is the citable scope. Minor. |
| 18 | D.5 "the arbitration is single-site" (Ch3) | Ch2 §3.2/§4.2/§4.4: Blois = weakest inter-sensor agreement (r = 0.47), flat d_opt surface, largest residual bias (+3.1 after d_opt, +0.7–0.8 after LUT), narrowest S2 IQR (0.42). | Not contradicted, but the Discussion never says the single site of Ch3 is Ch2's *hardest* site. See D-3. |
| 19 | D.7.4 "Dedicated field LAI campaigns … would convert consistency statements into accuracy statements" ; "Conifer and mixed stands are the stated boundary of every chapter" | §4.5 direction 3 ("quantifying the absolute accuracy gain … separately from the inter-sensor consistency improvement"); §4.5 limitation (i). | **Correct.** |
| 20 | D.7.4 "Reprocessing the raw point clouds to retain sub-2 m returns would test … the understory vegetation that may explain why Sentinel-2 ranks open plots" | §4.4 "The h_min sensitivity analysis suggests that the LAI_S2 signal **at Blois** is partly supported by canopy returns in the 2–5 m stratum: removing this layer … degrades all metrics more steeply at Blois". | **Consistent and directly supported by Ch2 at the Ch3 site**, but Ch2 is not cited. See D-2. |
| 21 | D.7.1 "Every LiDAR quantity in the preceding chapters is a summer snapshot" | Table 1: Mormal ALS 2021-06-16 + 2021-07-19 (50 % reacquired one month later); §4.1 "foliage conditions remained near-maximal across both LiDAR flight dates". | Fine (Ch2 concludes marginal effect). One clause acknowledging the Mormal two-date acquisition would pre-empt the question. Minor. |
| 22 | D.1 / D.4 "optical indices saturate" / "optical indices stop resolving beyond LAI 5–6" | Ch2 uses hybrid PROSAIL inversion, not vegetation indices; §1 "Regardless of the method, S2-derived LAI is known to saturate above 5–6". | Terminology: "optical retrievals" (Ch2 is an RTM inversion). Minor. |
| 23 | Site names / pixel counts | Ch2: three sites, 5 000 stratified pixels per site, fCover > 90 %, h ≥ 10 m. Discussion never names Aigoual, Blois or Mormal, nor "three sites". | No error; D.5 would gain from stating "three sites in Ch2, one of them (Blois) in Ch3". |

---

## (B) Contradictions, ranked

### Severe

**B-severe-1 — D.1 attributes to Ch2 a caveat Ch2 does not contain, and D.3 says the opposite.**
D.1: "it flagged the caveat that mattered most for what followed: better agreement in leaf area
does not automatically propagate into better simulated microclimate." Ch2 §4.5 says instead that
LAI_S2_opt "may serve as an operational input … Forest microclimate modelling … is one such
application". D.3 quotes exactly that offer and says Ch3 "declined it". A reader who opens Ch2
finds the offer, not the caveat. Either (a) the thesis version of Ch2 gets a bridging sentence in
its conclusion (then D.1 is true), or (b) D.1 is rewritten so the caveat belongs to the General
Introduction/Ch3 design and Ch2 keeps its "offer". Option (b) is safer and keeps "the caveat
became the result" only if the caveat is re-homed. (Rewrite C-1.)

**B-severe-2 — Ch2's third factor (intra-pixel canopy heterogeneity, H3) is absent from the
Discussion's account of Ch2, although it is now in Ch2's title.** Ch2 §3.4/§4.4/§5: agreement
declines with DSM_SD; residual errors "remain largest in heterogeneous stands"; "a single fixed
d_opt is an effective averaging for closed, homogeneous canopies but breaks down for canopies with
discontinuous canopy tops"; the d_opt over-correction bias grows with heterogeneity. The
Discussion mentions heterogeneity only as an atmospheric-solver perspective (D.7.4). This is a
missing limb of the "consistency question" answer and it also feeds D.6 (edges, 1-D column) and
the Ch3 rule (a height gate applied to heterogeneous stands). (Rewrite C-2, and D-1.)

**B-severe-3 — D.3 "the best retrieval of the thesis … cuts RMSE by 1.2 m² m⁻²" hides that the
gain is measured against the depth-truncated LiDAR LAI and, at Blois, that the LUT prior is that
same truncated LAI.** Ch2 §2.5.4: metrics "computed between LAI_ALS_dopt and LAI_S2"; Table A.16
Blois full Pareto: "lai = LAIALS_dopt". So the product Ch3 uses at Blois was tuned to agree with a
6 m slice of the canopy. This is not a weakness of D.3's argument, it is its strongest support
("a retrieval tuned to match a saturating signal inherits the saturation" becomes literal), but
as written the sentence reads as accuracy language ("best retrieval") the thesis disowns in D.5
and Ch2 disowns in r3. (Rewrite C-3.)

### Moderate

**B-moderate-1 — Only the dense-end saturation is told; Ch2's Conclusion says the slope is driven
"primarily" by low-LAI overshoot.** D.1 "saturates … in dense canopy … divergence is one of dynamic
range" vs Ch2 §5 "driven primarily by the amplified gap … at low LAI_ALS values"; §4.2 "LAI_S2 also
tends to overestimate LAI_ALS below 4–5"; §4.5 "The optimisation was inherently overconstrained"
because of it. Both compress the slope. The Discussion should name both ends; the low end matters
because Ch3 lets S2 *rank* open plots even though Ch2 says S2 *overshoots* them (rank vs magnitude
is exactly D.5's distinction). (Rewrite C-4.)

**B-moderate-2 — "site means comparable at k = 0.65" is partly circular and untrue at Blois.**
Ch2 §3.1: k = 0.65 chosen for "RMSE minimisation and bias re-centring"; Blois S2 > ALS by ~1.1–1.2
m²/m². (Rewrite C-5.)

**B-moderate-3 — d_opt presented as a physical depth without Ch2's three qualifiers.** Ch2 §4.2
"empirical proxy … not a direct physical measurement"; Blois "statistical compromise … rather than
an inferred physical penetration depth"; Ajouts §5c "operationally useful descriptors of
inter-sensor agreement, but not … physical penetration depths"; DTM normalisation gives 10 m at
all sites. D.1's "Sentinel-2 reflectance is most informative about the upper 6–10 m" and D.2's
"the layers an optical sensor cannot see through" are physical readings. Ch3 is Blois-only and its
truncation scenario uses d_opt: the Blois caveat travels with it. (Rewrite C-6.)

**B-moderate-4 — "our own pooled Sentinel-2 response flattening near 4" contradicts Ch2 §3.1
(observed ceiling 6–7 m²/m²) unless attributed to Ch3.** Left open by the harsh review. (Rewrite C-7.)

**B-moderate-5 — D.5 lists two circularities; Ch2 r3 names three.** Missing: d_opt-selection
circularity (Ajouts §5c; A.16 phrase 2). Also "≤ 17 % of the reported improvement" → "of the RMSE
reduction". (Rewrite C-8.)

**B-moderate-6 — "No field LAI was measured in this thesis, at any site, in any chapter."** Ch2 §2.1
lists hemispherical photographs at the sites, "not directly used"; §4.1 frames the gap as
"concomitant" measurements. (Rewrite C-9.)

### Minor

- **B-minor-1** "improved inter-sensor consistency" (D.1), "improvement" (D.5): on Jérôme's r3
  red list; Ch2 now says "decomposition" and "reduced the inter-sensor RMSE relative to ATBD".
- **B-minor-2** "23% error" (D.1) → "23% rescaling"; "k = 0.65 harmonization" (D.3) → Ch2's
  "operational value … close to the field-reported range". State once that Ch1 uses k = 0.5 and
  Ch2/Ch3 use k = 0.65, so the D.2 thresholds (2–3, 3.9, 4, 5–6) are not all on one k-scale.
- **B-minor-3** "0.25–0.80" is Ch2's own text but Ch2 Table 6 implies 0.15–0.97; retained-config
  slopes are 0.49–0.62 (Fig. 7). Fix in Ch2 or quote the latter.
- **B-minor-4** "63% of the leaf area sits below the effective optical depth" is Ch3 (Blois);
  attribute.
- **B-minor-5** "optical indices" (D.1, D.4) → "optical retrievals" (Ch2 is an RTM inversion).
- **B-minor-6** "disagree substantially in closed forests" → "across diverse forest systems"
  (Woodgate's scope).
- **B-minor-7** "Three mutually independent supports" (D.5): the Ch2 consistency analysis at
  Blois and the Ch3 arbitration share the same ALS point cloud, k, h_min and d_opt; "independent
  in kind" is defensible, "mutually independent" is generous.
- **B-minor-8** D.3 heading "Retrieval accuracy is not functional fidelity" keeps the retired
  word in the one place a reader sees before the gloss; consider "Retrieval consistency is not
  functional fidelity" (the Ch1 quotation inside the section can stay attributed).

---

## (C) Suggested rewrites (old → new)

**C-1 (D.1, Ch2 paragraph, last sentence; addresses B-severe-1)**
OLD: "The chapter improved inter-sensor consistency, never accuracy, and it flagged the caveat
that mattered most for what followed: better agreement in leaf area does not automatically
propagate into better simulated microclimate."
NEW: "The chapter decomposed the divergence into three named components, the effective optical
depth, the LUT priors, and intra-pixel heterogeneity, and reduced the inter-sensor RMSE relative
to the ATBD baseline; it claimed consistency, never accuracy. It ended by proposing the optimized
product as an operational input for microclimate modeling. That proposal, not a caveat, is what
Chapter 3 inherited: whether better agreement in leaf area propagates into better simulated
microclimate was the open question the General Introduction had posed."
And in the next paragraph: OLD "Chapter 3 answered the forcing question, and the caveat became the
result." → NEW "Chapter 3 answered the forcing question, and the Introduction's caveat became the
result." (or "and the offer became the test").

**C-2 (D.1, insert one sentence after the d_opt sentence; addresses B-severe-2)**
NEW: "The third component was horizontal: agreement declined with within-pixel canopy-height
heterogeneity at all three sites, and the residual errors after both corrections remained largest
where canopy tops are discontinuous, so a single fixed d_opt is an averaging that holds for closed,
homogeneous canopies and breaks down elsewhere."

**C-3 (D.3, second currency; addresses B-severe-3)**
OLD: "Second, Chapter 2 produced the best retrieval of the thesis, a forest-tuned PROSAIL
configuration that cuts RMSE by 1.2 m² m⁻² and recenters the bias."
NEW: "Second, Chapter 2 produced the most inter-sensor-consistent retrieval of the thesis: a
site-specific PROSAIL LUT that reduces the RMSE against the depth-truncated LiDAR LAI by
1.2 m² m⁻² on average and re-centers the bias at two sites (a +0.7 m² m⁻² residual remains at
Blois). At Blois, the site of Chapter 3, the selected LAI prior is the truncated LiDAR
distribution itself, so the product was tuned, by construction, to the upper 6 m of the canopy."

**C-4 (D.1, first Ch2 sentence; addresses B-moderate-1)**
OLD: "Sentinel-2 LAI saturates against the LiDAR reference in dense canopy, and the divergence is
one of dynamic range rather than of central tendency: regression slopes remain at 0.25–0.80 after
every correction, while the site means are comparable at *k* = 0.65."
NEW: "Sentinel-2 LAI saturates against the LiDAR reference in dense canopy (ceiling 6–7 m² m⁻²
where LiDAR reaches 10–15) and overshoots it in sparse canopy below 4–5 m² m⁻²; the divergence is
therefore one of dynamic range rather than of central tendency, with regression slopes of about
0.5–0.6 in the retained configurations, while the site means are close at Aigoual and Mormal once
*k* = 0.65 is adopted, a value chosen in part to re-center the inter-sensor bias."

**C-5** is folded into C-4 (B-moderate-2).

**C-6 (D.1, d_opt sentence; addresses B-moderate-3)**
OLD: "Sentinel-2 reflectance is most informative about the upper 6–10 m of canopy (7 m across
sites), a thickness we termed the effective optical depth d_opt, against canopy tops of 17–32 m."
NEW: "The LiDAR integration depth at which the two products agree best, which we termed the
effective optical depth d_opt, is 6–10 m (7 m across sites, under top-of-canopy normalization)
against canopy tops of 17–32 m. It is an empirical proxy for the canopy thickness that dominates
the Sentinel-2 signal, not a measured penetration depth, and at Blois, the site of Chapter 3, the
agreement surface is flat enough that the 6 m value is a statistical compromise."

**C-7 (D.2; addresses B-moderate-4)**
OLD: "and the optical saturation threshold, placed at LAI 5–6 in the literature, with our own
pooled Sentinel-2 response flattening near 4."
NEW: "and the optical saturation threshold, placed at LAI 5–6 in the literature and observed as
a 6–7 m² m⁻² ceiling across the three sites of Chapter 2, with divergence from LiDAR growing
beyond 4–5 (the Blois plot pool of Chapter 3 flattens near 4)."
Also, one clause after "The three estimators differ": "and they do not share one extinction
coefficient (Chapter 1 at *k* = 0.5, Chapters 2 and 3 at 0.65)".

**C-8 (D.5; addresses B-moderate-5)**
OLD: "the consistency analysis carries a bounded prior circularity (≤ 17% of the reported
improvement) and an unbounded scoring circularity"
NEW: "the consistency analysis carries a bounded prior circularity (≤ 17% of the reported RMSE
reduction), an unbounded scoring circularity, and a d_opt selection that is itself reference-driven"

**C-9 (D.5, opening; addresses B-moderate-6)**
OLD: "No field LAI was measured in this thesis, at any site, in any chapter."
NEW: "No field LAI concurrent with the acquisitions was available, and none was used in any
chapter; the hemispherical photographs of the site network were not analyzed here."

**C-10 (D.2, 63 %)**
OLD: "In dense stands, 63% of the leaf area sits below the effective optical depth, and Chapter 3
showed…" → NEW: "In the dense Blois plots of Chapter 3, 63% of the leaf area sits below the
effective optical depth, and that chapter showed…"

**C-11 (D.1 and D.3, k wording)**
OLD "a uniform 23% error in leaf area" → NEW "a uniform 23% rescaling of leaf area (the ratio
between the *k* = 0.5 and *k* = 0.65 conventions)".
OLD "the *k* = 0.65 harmonization adopted for inter-sensor work" → NEW "the *k* = 0.65 operational
value adopted in Chapters 2 and 3".

**C-12 (D.7.4, sub-2 m sentence)**
OLD: "Reprocessing the raw point clouds to retain sub-2 m returns would test the one mechanism we
could not, the understory vegetation that may explain why Sentinel-2 ranks open plots."
NEW: "Reprocessing the raw point clouds to retain sub-2 m returns would test the one mechanism we
could not: Chapter 2 already found that, at Blois alone, removing the 2–5 m stratum degrades the
Sentinel-2/LiDAR agreement more than at the other sites, so a low-stature layer that Sentinel-2
sees and the buffering responds to may explain why it ranks open plots."

**C-13 (D.1 / D.4, "optical indices")** → "optical retrievals".

---

## (D) Ch2 points the Discussion should pick up but does not

1. **Heterogeneity (H3), now in the Ch2 title.** Agreement declines with DSM_SD; d_opt
   over-correction bias grows with heterogeneity at Aigoual/Mormal; "a single fixed d_opt …
   breaks down for canopies with discontinuous canopy tops"; adaptive/site-specific d_opt as a
   perspective (§4.4, §4.5). Natural hooks: D.1 (C-2), D.6 (the 1-D column and edges), D.7.4
   (heterogeneity-aware solvers are the same limitation seen from the atmosphere side).
2. **h_min sensitivity at Blois** (§4.4): the 2–5 m stratum supports the S2 signal specifically at
   Blois. Direct evidence for D.7.4's understory hypothesis and for why S2 works in open plots.
3. **Blois is Ch2's hardest site** (r = 0.47; flat d_opt surface; bias plateau +3.1 after
   truncation; narrowest S2 IQR 0.42; "intrinsic to the structural and optical specificity of the
   Blois forest"; multi-date S2 test rules out a date artifact). D.5's "single-site" caveat should
   say the single site is the one where inter-sensor consistency was weakest, which makes the Ch3
   verdict conservative for S2 in dense canopy but also means the transfer to Aigoual/Mormal is
   untested in the direction that would favor S2.
4. **Ch2 anticipated the FORMS-H coupling** (§4.5: "coupling LAI_S2_opt with canopy-height
   products … FORMS-H"). Ch3 realized it as a *gate*, not as a structural input. D.4 ("height as
   classifier, not lever") would gain from noting that Ch2 proposed height as an input and Ch3
   found its use as a switch; the two chapters converge on FORMS-H by different routes.
5. **d_opt robustness to (k, θ)** (§4.1: "d_opt remains within 5–8 m at Aigoual while it sticks
   to 6 m at Blois and 10 m at Mormal", A.12): this licenses transporting d_opt across the k
   conventions of the chapters; one clause in D.2 or D.5 would close the k question cleanly.
6. **Ch2's own perspective list** (§4.5): multi-date S2 for seasonal d_opt (fits D.7.1's
   phenology thread: does d_opt move through the season?); cross-site LUT transfer (which priors
   are generic forest priors); low-LAI overshoot as an unresolved mechanism. Only the field
   campaign and the conifer boundary are picked up.
7. **LiDAR is a proxy too** (§1, §4.1: "LAI_ALS should be interpreted as a structural proxy rather
   than as an absolute reference"; ~15 m²/m² tail; scan-angle 1.155× worst case; clumping
   underestimation). D.5 says "two imperfect instruments"; one clause on why the LiDAR magnitude
   is nonetheless the less saturating of the two would tie D.4's "unsaturated magnitude" to Ch2's
   own caveats.
8. **Mormal two-date acquisition** (Table 1, §4.1): marginal, but "every LiDAR quantity … is a
   summer snapshot" (D.7.1) can carry "(two summer dates at Mormal)".
9. **The retired-vocabulary sweep** requested for r3 ("decomposition, not improvement") applies to
   the Discussion's Ch2 sentences as much as to Ch2 itself.

## (E) Verdict

The Discussion respects Ch2's final positioning where it matters (consistency, no ground truth,
≤ 17 % bound, offer-then-decline in D.3), and every Ch2 number it quotes is traceable to r2/r3.
Three things need fixing before a jury reads both texts side by side: D.1 credits Ch2 with a
microclimate caveat Ch2 never states (and D.3 contradicts it); Ch2's third component,
heterogeneity, now in its title, is missing from the Ch2 account; and "best retrieval … cuts RMSE
by 1.2" must say the gain is against a depth-truncated LiDAR LAI whose distribution is the Blois
LUT prior. Then de-physicalize d_opt (Blois caveat, CHM-only), name the low-LAI overshoot,
re-home the "near 4" plateau to Ch3, and add the d_opt-selection circularity to D.5.
