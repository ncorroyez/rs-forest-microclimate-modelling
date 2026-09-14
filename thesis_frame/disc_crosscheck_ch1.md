# Cross-check: General Discussion vs Chapter 1

*Adversarial examiner pass, 2026-08-18. Target: `thesis_frame/discussion_generale_these_EN.md`
(draft 2026-08-13/14, 302 lines, "DISC" below). Chapter 1 reference: the canonical
`chapter1/manuscript/manuscript_chap1_EN_native20.md` (2026-08-14, 873 lines, "CH1" below).
Items already marked fixed in `discussion_generale_these_EN_harsh_review.md` (oracle numbers,
2022 arithmetic, "accurate" gloss, three-lines numerology, spring causal verb, night clause,
composition effect, montane GEDI failure, out-of-sample wording) are not repeated.*

## 0. A version trap that precedes every other finding

The file named as "CH1 manuscript" in the brief, `review/manuscript_chap1_EN.md`, is the
**stale 25 June copy** (76 KB, *r* = 0.93, +0.9 °C, "Shapley-style", trait steps ±1 LAI /
±5 m, no *k* = 0.65 test, no dense/sparse split, no Appendix F on the leaf-area scale). The
these-manuscript skill flags it explicitly as "do not source facts from it". Against that
file, **nearly every Chapter 1 number the Discussion quotes is wrong** (table 0 below).
Against the canonical native20 manuscript, the numbers are right. This is not a Discussion
error, but it is an assembly risk: if the thesis is assembled from `review/`, the jury will
read a Chapter 1 that contradicts its own General Discussion on the headline validation.

Table 0. Same claim, two Chapter 1 versions.

| DISC claim (line) | Stale `review/manuscript_chap1_EN.md` | Canonical `native20.md` |
|---|---|---|
| ranking *r* = 0.50 (l.20) | *r* = 0.93 (l.30, 413) | *r* = 0.50, CI [0.26, 0.68] (l.38, 249) |
| warm bias +0.98 °C (l.21) | "about +0.9 °C" (l.31, 414) | +0.98 °C (l.38, 249, 302) |
| ~15% of amplitude (l.21) | absent ("compresses the slope toward unity") | "about 15%" (l.39, 249) |
| 23% cut erases skill; *k* = 0.65 (l.25, 87) | absent | App. F l.612-631 |
| dense *r* = 0.69 / open +0.20 (l.72) | absent | l.249 |
| −0.26 °C per 0.5 LAI, densest (l.70) | −0.17 °C per 1 LAI in P4 (l.340) | −0.264 for LAI +0.5 in P4 (l.225, 670) |
| height ≤ 0.004 °C m⁻¹ (l.116) | "below 0.07 °C even at ±5 m" (l.346) | ≤ 0.004 °C m⁻¹ (l.225, 651) |
| profile: largest lever mid-density, upper bound (l.17-18) | "co-leads in dense P4" (l.328) | "largest single lever in P3", upper bound (l.229-231) |
| profile increment indistinguishable from 0 (l.114) | "lifts *r* 0.86 → 0.90 in P4" (l.425) | +0.026, CI [−0.007, +0.070] (l.253) |

**Action:** replace or clearly quarantine `review/manuscript_chap1_EN.md` before assembly
(do not delete without Nathan's go). Everything below is checked against native20.

## A. Claim-by-claim table

Verdicts: OK = matches; OK* = matches but a Chapter 1 caveat is dropped; DRIFT = wording or
scope differs; MISMATCH = number/fact differs; N/A = not a Chapter 1 claim (source elsewhere).

| # | DISC quote (line) | CH1 quote (line) | Verdict |
|---|---|---|---|
| A1 | "Leaf quantity, read as leaf area and cover together, governs the realized ΔTmax buffering at both ends of the density gradient" (l.15-17) | "Leaf area and cover together dominate the simulated ... ΔTmax at both ends of the density gradient" (l.33); "we read leaf area and cover together as leaf quantity throughout, and never rank one against the other" (l.296); title: "Horizontal canopy density ... governs simulated understory microclimate buffering" (l.2) | OK; DRIFT on the headline term: CH1's title term "horizontal canopy density" never appears in DISC |
| A2 | "the vertical arrangement of foliage is the largest single lever in mid-density canopy, but only as an upper bound that real canopies barely span" (l.17-19) | "On this metric and this footing the swap is the largest single lever in P3" (l.229); "that potency is an upper bound, not a realized effect" (l.231); "though not on the buffering slope" (l.36); "Its lead in intermediate canopy is therefore specific to ΔTmax" (l.235) | OK* (drops "on ΔTmax only, not on the slope") |
| A3 | "the field tests place its realized effect below what *n* = 53 loggers resolve" (l.19-20) | "With 53 loggers the field tests place the profile's realized effect below what they resolve, and a tendency in the H2 direction survives in the densest third (partial *r* = −0.41, *p* = 0.09, *n* = 18)" (l.342); "It is a bounded negative rather than an exclusion" (l.342) | OK* (drops "bounded, not excluded" and the dense-third tendency) |
| A4 | "validation triplet (ranking *r* = 0.50, warm bias +0.98 °C, roughly 15% of the observed amplitude)" (l.20-21) | "*r* = 0.50, 95% CI [0.26, 0.68], *n* = 53 ... warm bias of +0.98 °C and capturing about 15% of the observed ΔTmax amplitude" (l.249) | OK |
| A5 | "the mechanistic model ranks stands far better than it reproduces them" (l.22) | "reproduces the between-plot ranking of buffering **only moderately** (*r* = 0.50)" (l.249); "recovers the buffering *slope* ... much better than the amplitude (*r* = 0.78 ... versus 0.50)" (l.249); "the process model is not justified by predictive skill at this site" (l.308) | DRIFT: CH1 calls the ΔTmax ranking "moderate"; what it recovers "far better" is the slope (0.78), which DISC does not name |
| A6 | "a uniform 23% error in leaf area erases the model's buffering skill" (l.25-26) | "A uniform 23% **reduction** in retrieved leaf area, applied identically to every plot, erases the model's buffering skill" (l.338); only the reduction was run (App. F l.612-631) | DRIFT: "error" (two-sided) vs "reduction" (one direction tested) |
| A7 | "That is precisely the regime in which optical indices saturate" (l.26) | "which is precisely the regime in which optical indices saturate" (l.338) | OK (verbatim) |
| A8 | "Chapter 1 anticipated this geometry ... Its mechanistic levers are strongest in dense canopy (−0.26 °C per 0.5 LAI in the densest archetype)" (l.69-71) | "in the densest archetype P4, adding 0.5 of leaf area cools the understory by 0.26 °C" (l.225); but the profile lever "peaks in intermediate canopy" (l.223) and "in the dense P4 ... its sign flipping from cooling to warming" (l.229); in P1 "its leaf-area lever warms rather than cools (+0.05 °C)" (l.227) | OK for leaf quantity and cover; DRIFT for "its levers" (plural): the profile lever peaks in P3 and flips sign in P4; the P1 sign reversal is not mentioned |
| A9 | "validation skill concentrates in the same regime (dense *r* = 0.69; open *r* = +0.20, not significant)" (l.71-72) | "split at the median leaf-area index, the model tracks the dense ones (*r* = 0.69, CI [0.42, 0.85], *n* = 27) but has no significant skill on the sparse ones (*r* = +0.20, CI [−0.20, +0.55], *p* = 0.32, *n* = 26)" (l.249) | OK ("open" vs CH1 "sparse"; split = median LAI, n given only in CH1) |
| A10 | "The model and the optical sensor thus fail in opposite halves of the gradient" (l.72-73) | CH1: open-stand failure is "a structural property of the homogeneous one-dimensional column, which cannot generate the gap-driven, horizontally advected warming of an open canopy" (App. C l.483) | **TENSION with Ch3** (S2-forced MuSICA reaches ΔTmax *R*² = 0.71 in the open half): see B1 |
| A11 | quotation "the leaf area that is most accurate for retrieval is not the leaf area that best couples to the microclimate" (l.80-81) | "The leaf area that is most accurate for retrieval is therefore not the leaf area that best couples to the microclimate" (l.620) | OK (verbatim minus "therefore"); but CH1 l.613 asserts *k* = 0.65 "is more accurate for the leaf area", contradicting DISC §D.5 "No field LAI was measured" (see B4) |
| A12 | "a uniform 23% reduction of the LiDAR leaf area, the direction and order of the *k* = 0.65 harmonization adopted for inter-sensor work, degrades the ranking skill from *r* = 0.50 to 0.19" (l.87-89) | "*k* = 0.65 ... reduces it by a further 23% ... the ΔTmax correlation falls from *r* = 0.50 to 0.19, the amplitude captured from 15% to 3%, and the warm bias grows from +0.98 to +1.33 °C" (l.612-615) | OK on numbers; DRIFT on the name: CH1 "optimized value ... for leaf-area retrieval accuracy", DISC "harmonization" |
| A13 | "the profile increment is statistically indistinguishable from zero" (l.114-115) | "That increment is +0.026, with a 95% bootstrap interval of −0.007 to +0.070. It is not distinguishable from zero. The data exclude a profile increment larger than about 0.07" (l.253) | OK |
| A14 | "maximum height is a negligible lever (≤ 0.004 °C m⁻¹)" (l.115-116) | "Maximum height is negligible in every archetype, \|effect\| ≤ 0.004 °C per meter" (l.225); "The height step therefore reads the radiative and geometric pathway of height and not its aerodynamic one, which makes it a **lower bound** on the full height effect" (l.157) | OK* (drops the lower-bound caveat: the wind at canopy top was held fixed under the height step) |
| A15 | "the three-dimensional detail that distinguishes LiDAR ... is, for this application, mostly a means of measuring that quantity without saturating" (l.120-122) | CH1 4.4: "The vertical arrangement of foliage is the weaker realized lever, but it is not beside the point. Thinning from below and thinning from above shift the vertical balance in opposite directions. The controlled test moves ΔTmax by up to 0.51 °C" (l.328); App. B: "the vertical dimension matters ecologically even where it moves our 1 m response metric little" (l.389) | OK* (CH1's two reservations, silviculture and other strata, are dropped) |
| A16 | "mechanistic attribution (Chapter 1), whose honesty conditions were physical parameters not fitted to the loggers, on-manifold perturbations, and a validation kept separate from the design" (l.135-137) | "MuSICA is calibrated on nothing here" (l.308); "the cLHS design points are not co-located with the loggers, so validation is stand-level rather than point-by-point" (l.336) | OK |
| A17 | "the attribution has no direct empirical test of itself" (l.145) | "Nor does the density-dependence rest on the perturbation alone: the field-anchored forward inclusion (Section 3.3) and the model-free observational regression (Appendix D) independently place leaf quantity ahead of the profile" (l.316); "the model's split *between* leaf area and cover is a model-only attribution the field cannot confirm" (l.194); "Only the variable *ranking* ... is claimed to carry beyond this dataset" (l.630-631) | DRIFT: too strong on one side (the ranking is corroborated model-free), too weak on the other (what has no test is the magnitude, the LAI-vs-cover split, and the transfer of the absolute coupling) |
| A18 | "53 independent temperature loggers" (l.143) vs "its non-independent loggers make the reported *R*² optimistic" (l.148) | CH1: residual Moran's *I* = 0.049, *p* = 0.098, "this bounds the residual structure rather than excluding it" (l.489) | Internal DISC inconsistency in the word "independent" (independent of the LAI products, not spatially independent) |
| A19 | "a warm bias near +1 °C that survives wind and scan-angle corrections" (l.157) | "+0.92 °C with neither correction, +1.03 °C with the wind correction alone, and +0.98 °C with both" (l.302) | OK; but Ch3 reports +0.73 °C for the same LiDAR full-column forcing at the same 53 loggers (Ch3 §3, "+1.14 °C against 0.73 °C"), which DISC quotes at l.96 without reconciling the two biases (see B5) |
| A20 | "compressed between-plot amplitude (about 15% of the observed range in Chapter 1 ...)" (l.158-159) | "reach roughly 15% of the amplitude" (l.302) | OK |
| A21 | "a nighttime warm bias of +1.5–1.7 °C" listed among the "recurring pattern across Chapters 1 and 3" (l.155-160) | "Nighttime buffering ... is a distinct mechanism we do not examine" (l.330) | N/A: not a Chapter 1 result; attribute to Chapter 3 only |
| A22 | "a skill that degrades on the hottest days precisely where predictions matter most" (l.160-161) | hot days: ΔTmax *r* = 0.35, slope 0.63 (Fig. S3, l.815); "Every lever is larger under heat" (l.798); "the dense-canopy regime ... is most pronounced precisely when buffering is most consequential" | OK; CH1 adds that the levers strengthen on those same days (missing point, see D) |
| A23 | "part of the measured gap is protocol rather than physics, since **unshielded loggers** warm in open stands [@terando...]" (l.165-167) | "inside a white homemade PVC radiation shield (10 × 15 cm)" (l.121); "**non-aspirated shields** can read warm in open, sunlit, low-wind stands [@terando...]" (l.304); "the 1 m versus 1.5 m measurement height" (l.304) | **MISMATCH**: the loggers are shielded; the issue is non-aspirated shields (and the 1 m vs 1.5 m mismatch) |
| A24 | "The one-dimensional column omits edge effects, lateral advection, and cold-air drainage" (l.164) | "slope, aspect, cold-air pooling, lateral advection, edge effects and sub-pixel gaps are not inputs" (l.306) | OK |
| A25 | "the margin of progress has shifted from the forcing to the model" (l.162-163) | "The claim is that the bias is not an artifact of the two corrections applied, not that it is insensitive to the leaf-area scale" (l.302); "Reconciling [*k* = 0.5 with the ellipsoidal extinction inside MuSICA] is a worthwhile refinement rather than a detail" (l.334) | DRIFT: CH1 keeps the leaf-area *scale* open as a forcing lever (only a reduction was tested, and it worsened the bias); DISC closes it |
| A26 | "Reprocessing the raw point clouds to retain sub-2 m returns would test the one mechanism we could not" (l.276-278) | "Rebuilding the profiles at 0.5 m to capture the sub-1.5 m understory does not promote the profile. The recovered understory loads almost entirely onto leaf quantity ... and the ranking is preserved" (l.336); LAD "binned at 1 m resolution from 1.5 m above ground" (l.131) | DRIFT: CH1 already did a 0.5 m rebuild at the 53 loggers (ranking only, separate lineage); the untested part is its effect on the Ch2/Ch3 products and on the open-plot ranking, and the threshold is 1.5 m in CH1, 2 m in DISC |
| A27 | "Conifer and mixed stands are the stated boundary of every chapter" (l.281-282) | "In evergreen needleleaf canopies, plant area keeps mattering through needle clumping. The single-species scope is thus a testable mechanistic prediction, not only a limitation" (l.298) | OK* (CH1 turns the boundary into a prediction; DISC keeps it as a boundary) |
| A28 | "microclimate regime transition near LAI 2–3 (a segmented fit to the observed buffering slope)" (l.56-57) | no segmented fit in CH1 | N/A (Ch3: LAI ≈ 2.9) |
| A29 | "the mean buffering slope did not collapse (0.93 against 0.90)" (l.236) | "The micro-macro slope ... comes out at 0.91 simulated against 0.89 observed" (l.302) | MISMATCH minor: 2021 observed slope 0.90 (App. G, *n* = 51, mean) vs 0.89 (CH1, *n* = 53); harmonize or state the subset |
| A30 | "leaf quantity as the controlling structural variable, measured without saturation where buffering is deepest" (l.293-294) | "The leaf-quantity lever is strongest in closed canopy where buffering is deepest" (l.342) | OK |
| A31 | "LAI" throughout DISC | "LAI is the vertically integrated one-sided LAD, an effective plant area index ... a proxy for true leaf area"; "does not separate woody from foliar surface ... We neither quantify that fraction nor bound its consequence" (l.129, 336) | OK* (DISC never says the LiDAR "leaf area" is a plant-area proxy; relevant to §D.5) |
| A32 | "Shapley" | absent from both; CH1 uses "Shapley-style decomposition" only as the rejected alternative (l.314) | OK (no drift; the brief's "density-dependent Shapley" is the superseded framing) |

## B. Contradictions and gaps, ranked

### Severe

**B1. "The model fails in the open half" (Ch1) vs "Sentinel-2 suffices where the canopy is open" (Ch3), left unreconciled in §D.2 (l.72-75).**
CH1 diagnoses the open-canopy failure as structural: "The failure is therefore a structural
property of the homogeneous one-dimensional column, which cannot generate the gap-driven,
horizontally advected warming of an open canopy, not a clipping of the cover input" (App. C
l.483; also l.306, l.346 "fails to reproduce the amplification of the most open stands").
DISC §D.1 (l.42-43) reports that the same column, forced with Sentinel-2 LAI, ranks the open
plots at *R*² = 0.71. Both cannot be read literally: either the open-stand failure is a
*forcing* failure (LiDAR LAI does not carry the open-plot signal, but some other field does),
or the Sentinel-2 open-plot skill is a proxy correlation the model does not earn (Ch3's own
model-free test finds S2 LAI *sign-inverted* against observed buffering, *r* = +0.19, DISC
l.94; Ch3 suspects understory vegetation, DISC l.277-278). §D.2's sentence "The model and the
optical sensor thus fail in opposite halves" papers over this. A jury member who reads App. C
of Ch1 and Fig. J of Ch3 side by side will ask which it is. Add one paragraph that names the
tension and adopts a reading (the proxy reading is the one consistent with Ch1 l.227: "The
P1 sensitivities are the response of a configuration the model does not reproduce", and with
Ch1's P1 LAI lever being sign-reversed, +0.05 °C, so no LAI-pathway forcing can rank open
plots through the mechanism Ch1 attributes).

**B2. "unshielded loggers" (l.166) is factually wrong about the thesis's own network.**
CH1 l.121: loggers sit "inside a white homemade PVC radiation shield"; the Terando point is
about *non-aspirated* shields (l.304), plus the 1 m vs 1.5 m height mismatch and the
north-side trunk microhabitat. Correct the word; a jury will read it as not knowing one's
own protocol.

**B3. Stale Chapter 1 copy in `review/` (section 0).** Not a Discussion error, but the highest
practical risk: nine headline numbers differ between the two files.

### Moderate

**B4. Ch1 App. F asserts what §D.5 denies.** CH1 l.613: *k* = 0.65 "is more accurate for the
leaf area"; l.622: "optimized for a different objective, leaf-area retrieval accuracy against
a reference". DISC l.132: "No field LAI was measured in this thesis, at any site, in any
chapter", and l.82-83 glosses "accuracy" as inter-sensor consistency. The gloss is right, but
the quoted Ch1 sentence carries an accuracy claim into the Discussion by attribution. Either
Ch1 App. F is reworded before assembly (recommended: "more consistent with the Sentinel-2
retrieval", "optimized for inter-sensor consistency") or DISC l.80-84 must say that Ch1's
wording predates the consistency framing. Also terminology: DISC "harmonization" (l.88) vs
CH1 "optimized value" (l.612).

**B5. Two warm biases for one LiDAR forcing.** DISC quotes +0.98 °C (Ch1, l.21, 157) and
+0.73 °C for the "full column" (Ch3, l.96) and then folds both into "a warm bias near +1 °C"
(l.157). Same site, same 53 loggers, same nominal LiDAR forcing, 0.25 °C apart. If the two
chapters run different clocks, forcing lineages or cover floors, say so in one clause; if
not, one of the numbers is wrong. (App. G's observed 2021 mean ΔTmax is −0.73 °C, *n* = 51,
which suggests the Ch3 "bias" may be measured against a different baseline; check.)

**B6. §D.1 l.22 "ranks stands far better than it reproduces them" misreads Ch1's own
verdict.** Ch1 calls the ΔTmax ranking "only moderately" reproduced (*r* = 0.50); what it
recovers well is the buffering *slope* (*r* = 0.78, l.249) and absolute hourly temperature
(*R*² = 0.91, l.249). Ch1 also concedes that a fitted regression on the same loggers reaches
*R*² = 0.79 (l.308), so "the process model is not justified by predictive skill at this
site". The Discussion should say which quantity is ranked well (the slope) and carry the
model-free comparison, which strengthens §D.6 rather than weakening §D.5.

**B7. §D.6's "no forcing repairs" ceiling is asserted in the one direction Ch1 did not
test.** Ch1 App. F cut leaf area (bias +0.98 → +1.33 °C) and Ch3 truncated it (+0.73 →
+1.14 °C); nobody increased it (a *k* < 0.5 or a clumping correction would raise the LiDAR
LAI). Ch1 l.302 explicitly limits the invariance claim: "not that it is insensitive to the
leaf-area scale", and l.334 lists the *k* = 0.5 vs ellipsoidal-extinction reconciliation as
a "worthwhile refinement". §D.6 should keep the ceiling claim conditional on the tested
directions, and §D.7.4 could name the untested one (it is cheap: a rerun of the 53-logger
validation at *k* = 0.4).

**B8. §D.5 residue for Chapter 1 (l.145) is misassigned.** "no direct empirical test of
itself" contradicts Ch1 l.316 (ranking corroborated by forward inclusion and the model-free
regression). What has no test is: the magnitudes (model-internal), the LAI-vs-cover split
(l.194, "a model-only attribution the field cannot confirm"), and the transfer of the absolute
coupling ("Only the variable *ranking* ... is claimed to carry beyond this dataset", l.630).
The last point also qualifies §D.1 l.24-26: the "23% bar" is a site-conditional number
(single site, one summer, one species, unmeasured LAI scale), yet DISC uses it as the
transferable requirement on optical retrievals.

**B9. Height "negligible as a physical lever" (l.116, 123) drops Ch1's lower-bound caveat.**
CH1 l.157: the height step moves the LAD depth and the forcing level "but not the wind speed
delivered there ... a lower bound on the full height effect". §D.4's classifier/lever
distinction survives, but the lever statement should say "through its radiative pathway".

**B10. §D.7.4 l.276-278 presents as untested what Ch1 partly tested.** Ch1 rebuilt profiles
at 0.5 m at the 53 loggers (l.336): the recovered sub-1.5 m stratum "loads almost entirely
onto leaf quantity, its sensitivity roughly doubling", ranking preserved. What remains
untested is (i) the same rebuild in the Ch2/Ch3 pixel pipeline (h_min = 2 m there) and (ii)
whether that stratum explains the open-plot ranking that Sentinel-2 captures. Say that.
Also harmonize the threshold (Ch1: 1.5 m filter, cover above 2 m; DISC: "sub-2 m").

**B11. Ch1's P1 sign reversal is absent from §D.2.** In the open archetype the model's
leaf-area lever *warms* (+0.05 °C per 0.5 LAI, l.225, 227, 292) and the model "holds
ΔTmax near zero" where loggers amplify. This is the mechanistic reason the LiDAR-forced
model cannot rank open plots (A9) and it bears directly on B1. One clause in l.69-72.

### Minor

**B12.** l.17-18: add "on ΔTmax, not on the buffering slope" to the mid-density profile
lever (CH1 l.36, 235: on the slope the swap "never leads in the intermediate archetypes").
**B13.** l.19-20 / l.113-116: Ch1 frames the profile result as "a bounded negative rather
than an exclusion" with a dense-third tendency (partial *r* = −0.41, *p* = 0.09, *n* = 18,
l.342, 574); §D.4 "contributes little" should carry the bound.
**B14.** l.25 "23% error" → "23% reduction" (only the cut was run).
**B15.** l.143 "53 independent temperature loggers" vs l.148 "non-independent loggers":
say "independent of every LAI product".
**B16.** l.159-160 nighttime bias listed under "Chapters 1 and 3"; Ch1 does not examine
night (l.330). Attribute to Ch3.
**B17.** l.236: 2021 slope 0.90 (App. G, *n* = 51) vs Ch1 0.89 (*n* = 53): state the
subset or harmonize.
**B18.** Terminology: Ch1's headline term "horizontal canopy density" (title, l.284, 294,
342) never appears in DISC; DISC uses "leaf quantity" only. One bridging clause in §D.1.
Ch1 also says "variables" (not "traits") since the native20 revision; DISC is neutral. Ch1
"sparse"/"dense" halves vs DISC "open"/"dense": acceptable but say "split at the median
LAI" once.
**B19.** l.129 "leaf area" for the LiDAR product: Ch1 l.129/336 calls it an effective
plant area index (wood included, unquantified). §D.5 would be stronger for saying so.

## C. Suggested rewrites (old → new), Discussion only

C1 (l.20-22). OLD: "established at the same time that the mechanistic model ranks stands far
better than it reproduces them."
NEW: "established at the same time that the mechanistic model recovers each plot's position on
the buffering gradient (slope *r* = 0.78) far better than the amplitude of its offset
(*r* = 0.50, about 15%), and that a regression fitted on the same loggers ranks them better
still (*R*² = 0.79): the model earns its place by naming pathways, not by predictive skill."

C2 (l.24-26). OLD: "a uniform 23% error in leaf area erases the model's buffering skill."
NEW: "a uniform 23% reduction in leaf area, the size of the *k* = 0.5 → 0.65 rescaling, erases
the model's buffering skill (*r* 0.50 → 0.19); the bar is a single-site number and only the
ranking, not the coupling, was claimed to transfer."

C3 (l.17-19). OLD: "the vertical arrangement of foliage is the largest single lever in
mid-density canopy, but only as an upper bound that real canopies barely span, and the field
tests place its realized effect below what *n* = 53 loggers resolve."
NEW: "the vertical arrangement of foliage is the largest single lever in mid-density canopy,
on ΔTmax and not on the buffering slope, and only as a full real-versus-uniform swap that real
canopies barely span; the field tests place its realized effect below what *n* = 53 loggers
resolve, a bounded negative rather than an exclusion."

C4 (l.69-75, replace the paragraph's last two sentences). OLD: "The model and the optical
sensor thus fail in opposite halves of the gradient, which is why the fusion rule of
Chapter 3 is a genuine resolution rather than a compromise: each half of the domain is served
by the instrument that resolves it."
NEW: "In the open half the model's own leaf-area lever reverses sign (+0.05 °C per 0.5 LAI in
the open archetype) and the column holds the understory near neutral where the loggers
amplify; Chapter 1 reads that failure as structural. Chapter 3 nevertheless ranks the same
open plots at *R*² = 0.71 when the column is forced with Sentinel-2. We do not read this as
the column recovering a mechanism it lacks: the optical field is sign-inverted against
observed buffering in the model-free test (*r* = +0.19), so its open-canopy skill is a proxy
association, plausibly carried by understory or ground vegetation that the sub-2 m LiDAR
filter removes, and it should be used as a ranking device there, not as evidence that the
one-dimensional column resolves open stands. On that reading the fusion rule of Chapter 3 is
a resolution rather than a compromise: each half of the gradient is served by the field that
ranks it, for reasons that differ between halves."

C5 (l.113-116). OLD: "the profile increment is statistically indistinguishable from zero, and
maximum height is a negligible lever (≤ 0.004 °C m⁻¹)."
NEW: "the profile increment to the field fit is statistically indistinguishable from zero
(+0.026, 95% CI −0.007 to +0.070), and maximum height is a negligible lever through its
radiative pathway (≤ 0.004 °C m⁻¹, the aerodynamic pathway held fixed, so a lower bound)."

C6 (l.144-146). OLD: "the attribution has no direct empirical test of itself;"
NEW: "the attribution's ranking is corroborated model-free (nested regression, *R*² = 0.79)
but its magnitudes are model-internal, its split of leaf quantity into leaf area and cover
is one the collinear field cannot confirm, and only the ranking was claimed to transfer;"

C7 (l.155-160). OLD: "a nighttime warm bias of +1.5–1.7 °C,"
NEW: "a nighttime warm bias of +1.5–1.7 °C (Chapter 3; Chapter 1 did not examine the night),"

C8 (l.161-163). OLD: "These defects move together when the forcing changes, but they do not
close. We conclude that the margin of progress has shifted from the forcing to the model."
NEW: "These defects move together when the forcing changes, but they do not close, and every
forcing change tested so far moved leaf area down (a 23% cut in Chapter 1, an optical
truncation in Chapter 3), which worsened them; no chapter tested a higher leaf-area scale
(*k* < 0.5, or a clumping correction), the one direction that could reduce the warm bias
inside the present column. Subject to that untested direction, we conclude that the margin of
progress has shifted from the forcing to the model."

C9 (l.165-167). OLD: "since unshielded loggers warm in open stands [@terandoAdHocInstrumentation2017]."
NEW: "since non-aspirated shields read warm in open, sunlit, low-wind stands
[@terandoAdHocInstrumentation2017] and the loggers sit at 1 m against a 1.5 m reference."

C10 (l.276-278). OLD: "Reprocessing the raw point clouds to retain sub-2 m returns would test
the one mechanism we could not, the understory vegetation that may explain why Sentinel-2
ranks open plots."
NEW: "Chapter 1 rebuilt the logger profiles at 0.5 m and found the recovered sub-1.5 m stratum
loading onto leaf quantity without reordering the levers; carrying that rebuild into the
pixel-scale products of Chapters 2 and 3, which drop returns below 2 m, would test the one
mechanism we could not, the understory vegetation that may explain why Sentinel-2 ranks
open plots."

C11 (l.87-88). OLD: "the direction and order of the *k* = 0.65 harmonization adopted for
inter-sensor work,"
NEW: "the direction and order of the *k* = 0.65 rescaling that maximizes inter-sensor
consistency in Chapter 2 (Chapter 1's Appendix F calls it 'optimized for retrieval', a
wording that predates the consistency framing),"

C12 (l.15). OLD: "Leaf quantity, read as leaf area and cover together, governs"
NEW: "Horizontal canopy density, in Chapter 1's term, leaf quantity read as leaf area and
cover together, governs"

C13 (l.143). OLD: "53 independent temperature loggers judge"
NEW: "53 temperature loggers, independent of every LAI product, judge"

## D. Chapter 1 points the Discussion should pick up (or contradicts)

D1. **Silvicultural corollary** (CH1 l.328): thinning lowers both controlling levers; thinning
from above vs below shifts the vertical balance, a lever of up to 0.51 °C at the dense end;
"buffering lost to harvest returns only as leaf area and closure rebuild". The Discussion has
no management perspective at all, in a forest-microclimate thesis whose motivation (l.229-231
of DISC) is microrefugia. One paragraph in §D.7.4 or §D.8.

D2. **The process model loses to a fitted regression on its own validation set** (CH1 l.308,
*R*² = 0.79 vs *r* = 0.50): "not justified by predictive skill at this site. Its value is that
it names the pathway". This is the honest sentence §D.5 and §D.6 need; it also explains why
Chapter 3's arbitration leans on rankings.

D3. **Scale constraint on optical mapping** (CH1 §3.4, l.338): archetypes prefer different
footprints (closed 5-12.5 m; open never fits; pooled *r* 0.41-0.53 over 5-50 m); "a coarser
optical pixel inherits that compromise rather than escaping it". Directly relevant to the
10-20 m Sentinel-2 grid of Chapters 2-3 and absent from DISC.

D4. **Levers strengthen exactly where skill drops** (CH1 l.237, 798, 815): on the hottest 10%
of days every lever grows (P4 leaf area −0.36 °C) while ΔTmax *r* falls to 0.35. §D.6 quotes
only the skill drop; the pairing (highest sensitivity, lowest skill) is the sharper point for
§D.7.2's heatwave test.

D5. **k = 0.5 inversion vs ellipsoidal extinction inside MuSICA** (CH1 l.334): a named
inconsistency between retrieval and model, "worthwhile refinement". Fits §D.7.4 next to the
k-sensitivity of Chapter 2 and closes B7.

D6. **The Gril 2023 complementarity** (CH1 l.320): same site, statistical slope map at
*R*² = 0.91 from three collinear metrics; "the statistical model predicting where the forest
buffers and the mechanistic model explaining which variable makes it buffer". DISC §D.8
("not the map, but what the map must be made of") is the same idea and could cite it.

D7. **Deciduous vs needleleaf as a mechanistic prediction** (CH1 l.298), not just a boundary
(DISC l.281-282).

D8. **Plant area, not leaf area** (CH1 l.129, 336): the LiDAR "LAI" is an effective plant
area index with an unquantified woody fraction; §D.5's ground-truth paragraph should own it.

D9. **P1 sign reversal** (CH1 l.225-227, 292): adding foliage warms the open understory in the
model. Missing from §D.2 and load-bearing for B1.

D10. **Chapter 1's open-canopy limits also bound Chapter 3's "Sentinel-2 suffices where open"**:
no cLHS design point below cover 0.5 (CH1 l.328), several loggers below that floor, and the
model "least faithful" there. DISC l.42-43 states the open-canopy verdict without this bound.

## E. Verdict

1. Every Chapter 1 number in the Discussion matches the canonical native20 manuscript
   (0.50 / +0.98 / 15% / 0.69 / +0.20 / −0.26 / ≤0.004 / 0.026 / 23% / 0.19); against the
   stale `review/manuscript_chap1_EN.md` nine of them are wrong, so quarantine that file.
2. Two factual errors to fix now: "unshielded loggers" (they are shielded, non-aspirated) and
   the unreconciled +0.98 vs +0.73 °C warm bias for one LiDAR forcing.
3. One real cross-chapter tension is papered over in §D.2: Ch1 calls the open-canopy failure
   structural to the 1-D column, Ch3 ranks the same open plots at *R*² = 0.71 with S2; adopt
   the proxy reading explicitly (C4).
4. Chapter 1's caveats most often dropped: "on ΔTmax only", "bounded, not excluded", height as
   a radiative lower bound, ranking corroborated but magnitudes/split untested, only the ranking
   transfers, leaf-area scale still open (only cuts were tested).
5. Chapter 1 material the Discussion leaves on the table: silviculture, the fitted-regression
   beats-the-model admission, footprint/scale, levers-up-skill-down on hot days, k vs
   ellipsoidal extinction, plant area not leaf area.
