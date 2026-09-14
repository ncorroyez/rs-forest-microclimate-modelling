# Article plan — Chapter 1 (detailed, IMRaD)

**Angle:** ecological — the *relative* contribution of vertical leaf arrangement (LAD profile)
vs total leaf area (LAI) to summer understory buffering — with the methodological rigour
(attribution under collinearity, off-manifold control) given prominent weight.
**Scope:** single site (Blois, deciduous oak), one summer (2021), 53 HOBO loggers — stated up front.
**Headline:** **leaf quantity and cover (LAI & fCover, co-dominant) drive the simulated buffering — their
collinearity (r = 0.76) does NOT distort the attribution (conditional/off-manifold control, gap ≤ 0.02 °C)**;
the *vertical arrangement* adds almost nothing — and even where the model reproduces the between-plot
ranking (r ≈ 0.93) it hits a glass ceiling on the absolute level and on the amplifying plots.

ΔTmax = diurnal micro − macro offset (macro = open-field regional weather).

**Settled framing (no longer open):**
- *Method:* exact Shapley on MuSICA scenarios — the GAMM of the V1 is the **superseded** lock
  (collinearity / concurvity / black-box "Date"), cited in the Introduction as the bottleneck removed.
- *Story:* hypotheses expected the vertical profile to be decisive; the result is that **leaf quantity &
  cover (LAI, fCover) dominate and LAD ≈ 0**. This is reframed as a strong, publishable nuance of the
  "3-D is indispensable" dogma — not a failure. H2's mechanism is confirmed in silico but shown to be
  second-order and undetectable in the field.
- *Focus:* field-anchored. Sentinel-2 / H3 stays in the annex; the Chapter-2 bridge is discursive.

---

## Title (working)
> **Leaf quantity and cover, not vertical arrangement, govern understory microclimate buffering in
> LiDAR-driven MuSICA simulations.**
>
> *(Evolved from the conference title "Integrating LiDAR-Derived Vegetation Structural Attributes into the
> MuSICA Model to Map and Explain Fine-Scale Forest Microclimate Variability". Finding-forward, ~15 words,
> matching the declarative style of the target literature (e.g. "The (lack of) impact of canopy height and
> vertical complexity on understory microclimate in pine plantations"). Method woven in (LiDAR-driven MuSICA)
> rather than a colon tail; "Shapley" lives in the abstract/Methods, not the title; "map … variability"
> dropped — the glass ceiling means we capture the spatial ranking (r ≈ 0.93) but compress the amplitude, so
> the glass-ceiling caveat belongs in the abstract, not the title. US spelling "understory" — unify across MS.)*
> *(Framing: the contrast is **quantity/cover vs vertical profile**. LAI and fCover are **co-dominant** —
> their collinearity (r = 0.76) was tested with the same off-manifold conditional control used for LAD|Hmax
> (Annex A4b): the artefact is negligible (gap ≤ 0.02 °C, within CI), so each trait's buffering attribution
> is robust regardless of whether the partner is present. We report them as co-dominant but do NOT rank
> them (CIs fully overlap, φ_LAI ≈ φ_fCover ≈ −0.20 °C). The dominant contrast remains quantity/cover ≫ vertical.)*

## Research questions
1. How sensitive is the MuSICA-simulated understory temperature to the **vertical and horizontal**
   distribution of vegetation?
2. What is the **relative contribution of the vertical profile (LAD) vs total leaf area (LAI)** to ΔTmax?

## Hypotheses
- **H1** — Vertical structure modulates radiative/turbulent transfers, so the LAD profile is a key
  parameter for faithful microclimate simulation. *(partly refuted — see Results/Discussion.)*
- **H2** — At fixed LAI and height, **top-heavy** profiles (biomass concentrated high) cause stronger
  diurnal attenuation (radiation intercepted higher, heat sources shifted to the canopy top).
  *(mechanistically confirmed in silico, but second-order and undetectable in the field.)*
- **H3 (perspective, annex)** — The 2-D optical signal (Sentinel-2) gives an *illusion* of buffering
  (a structure proxy) → motivates Chapter 2.

---

## 1. Introduction
1. Forest microclimate decoupling and thermal refugia under warming (De Frenne et al.; Lembrechts
   et al.; Zellweger et al.).
2. Canopy structure as the driver — but **which dimension**: total leaf area vs the **3-D vertical
   profile**? Position relative to **Bouwen** (neighbouring LiDAR+MuSICA work; density threshold).
3. The V1 bottleneck: a monolithic GAMM with **trait collinearity / concurvity** and a black-box
   "Date" effect. This motivates an attribution that (i) explores **all trait combinations**,
   (ii) is robust under collinearity, and (iii) avoids off-manifold partial-dependence plots.
4. Our approach: **MuSICA scenarios + exact Shapley** on the 2⁴ trait lattice, validated against 53 HOBO.
5. Questions and H1–H2 (H3 as perspective).

## 2. Materials & Methods
> Reviewer guardrails baked in: scope paragraph stated early; **LAI ×2 = total (two-sided) leaf area**
> expected by MuSICA (never "clumping"); the vertical-gradient figure is a mechanistic illustration.

- **2.1 Site & sensors** — Blois oak forest, summer 2021, 53 HOBO at ~1 m. *Scope paragraph here.*
- **2.2 LiDAR & traits** — 0.5 m LAD profiles (resolving the understory); traits LAI, Hmax, fCover,
  **VCI**, and profile shape; LiDAR LAI is one-sided (projected) → **×2 = two-sided** area for MuSICA.
- **2.3 Typology** — **FPCA of the doubly-normalised LAD profile (normalised in height x AND density y
  → pure shape, scale-free of LAI/Hmax)**; **K-means on the FPCA scores** → **four archetypes P1–P4**
  (A1). The typology is therefore architectural-shape-based, independent of leaf amount. Clustering
  frozen; cLHS design.
- **2.4 MuSICA & simulation design** — **16 coalitions (2⁴)** per plot (each trait real vs baseline:
  LAI→mean, Hmax→mean, fCover→absent, LAD→uniform); fCover baseline = mean (0.869); cLHS LAI_b ≈ 6.73;
  legacy binary **v3.2.0**. (An in-silico constrained factorial and an in-situ stratified design are
  also possible; the implemented set is the coalition design.)
- **2.5 Extraction — two conventions (never merged)**: (a) absolute ΔTmax/ΔVPDmax → fixed 1 m interp,
  −2 h shift (HOBO-comparable); (b) micro–macro slope (buffering/amplification split) → nair==1, no shift.
- **2.6 Attribution** — **exact Shapley** φ per trait per plot (ΔTmax, ΔVPDmax), additivity verified,
  bootstrap CIs; this **isolates φ_LAD vs φ_LAI** (Question 2) in place of the GAMM. **Conditional Shapley**
  bounds the off-manifold artefact of correlated traits on **both** co-defined/collinear pairs:
  φ_LAD | Hmax (structural co-definition, A4) and φ_LAI | fCover & φ_fCover | LAI (statistical collinearity
  r = 0.76, A4b) — each trait's φ recomputed conditional on the partner being real vs baseline.
- **2.7 Validation** — REF vs 53 HOBO; daily ΔTmax; micro–macro slope (45/8 split). **ΔVPDmax is
  attributed but NOT validated** (HOBO measure air temperature only) → VPD results are model-internal,
  stated as such.

## 3. Results  (`claim → figure → number`)

### R1 — A four-archetype structural typology
→ **Fig 1** (`fig_archetypes_profiles` — mean LAD profile + ribbon + trait annotations, style MEB,
  LAI shown one-sided; + `tab_cluster_structure`)
→ four archetypes from open **P1** to dense **P4**.

### R2 — Leaf quantity & cover are co-dominant; vertical arrangement adds little *(core)*
*(Message: **{LAI, fCover} ≫ vertical**, and LAI & fCover are co-dominant. We pre-empt the collinearity
objection (r = 0.76) directly — see the off-manifold control below — rather than retreating to a vague
"combined factor".)*
→ **Fig 2** (`fig_shapley_ranking`) and **Fig 3** (`fig_shapley_by_cluster`); `tab_shapley_ranking.csv`
→ ranking **LAI ≈ fCover ≫ Hmax > LAD ≈ 0**; φ_LAI(ΔTmax) ≈ −0.2 °C, φ_fCover ≈ −0.2,
  φ_Hmax small, **φ_LAD ≈ 0 (ns)**; Σφ ≈ −0.38 °C (additivity holds); same order for VPD.
→ **Collinearity is not driving the split (Annex A4b, conditional Shapley LAI|fCover & fCover|LAI):** the
  same off-manifold *logic* used for LAD|Hmax — recompute each trait's φ conditional on its correlated
  partner being real vs stripped to baseline; agreement ⇒ no chimeric-coalition artefact — gives
  φ_LAI|fCover=real −0.204 vs |baseline −0.227 (**gap 0.022 °C**) and φ_fCover|LAI=real −0.205 vs |baseline
  −0.215 (**gap 0.010 °C**). Both gaps are an order of magnitude **below the effect size (≈ −0.20 °C)** and
  **inside the ±0.1 °C bootstrap CI**, so the co-dominance is not an artefact of correlated coalitions; the
  contributions are near-additive (if LAI and fCover were pure substitutes, |φ|partner=baseline would far
  exceed |φ|partner=real). *(Same conclusion as the LAD|Hmax control — the two conditionals agree — but a
  different mechanism: structural co-definition there, statistical collinearity here.)*
→ **Mean-level co-dominance ≠ spatial co-dominance — keep them distinct.** The pooled medians match
  (φ_LAI ≈ φ_fCover ≈ −0.20 °C) but their bootstrap CIs do not: φ_fCover is tight ([−0.216, −0.206], width
  0.01) while φ_LAI is wide ([−0.302, −0.103], width 0.20). The by-cluster decomposition (Fig 3) explains it:
  **φ_fCover is a near-constant buffering floor** (P2 −0.21, P3 −0.22, P4 −0.21; only sparse P1 +0.03),
  whereas **φ_LAI carries the between-plot variation** — a monotonic gradient with the sign-flip
  (P1 +0.15 amplify → P4 −0.56 buffer). So fCover and LAI are co-dominant *on average* but play different
  spatial roles: **fCover sets a uniform floor in closed canopy; LAI is the trait that differentiates plots**
  — which is precisely the "3-D distributes the between-plot variance" backbone (r ≈ 0.93). *Honest limit:*
  we report mean-level co-dominance but do **not rank** LAI vs fCover (CIs overlap), and they remain
  collinear in the field, so the precise 50/50 mean split is not itself the claim.
→ *At the ΔTmax scale, the vertical profile adds essentially nothing beyond leaf quantity and cover.*
  Lead sentence: *we attribute the buffering simulated by MuSICA, validated spatially against HOBO.*

### R3 — The 3-D distributes the variance, but a glass ceiling against the field
→ **Fig 4** (`fig_validation_perplot`; + `fig_validation_by_cluster`) and **Fig 5**
  (`fig_validation_forward_shapley`)
→ **between-plot ranking captured (r ≈ 0.93)**; warm **bias +0.9 °C**; field-HOBO split **45/8**, z05-model split **46/7**.
→ Per-plot slope, MuSICA vs HOBO, split buffering/amplifying (Annex A8b): the model **compresses the
  slope range toward 1** — under-buffers the bufferers (sim 0.94 vs obs 0.84) and under-amplifies the
  amplifiers (1.12 vs 1.27).
→ **Read A8b jointly with A8a — same phenomenon, not independent evidence:** the 8 amplifying plots ARE
  the gappy/steep/exposed ones (A8a), whose amplification is topo/heterogeneity-driven and **not a
  MuSICA input**, so the model "failing" to amplify them is a **scope limitation, not a structural
  error**. The only mismatch attributable to the model itself is the **slope compression on the
  closed-canopy buffering plots** (gap ≈ 0, topo negligible) — that is the genuine structural limit. The
  glass ceiling is therefore two things: (i) real structural under-buffering on closed canopy, and
  (ii) an out-of-scope inability to reproduce topo/gap-driven amplification.
→ **Night-time / ΔTmin** (Annex A5): the attribution and the buffering hierarchy hold for ΔTmin as well
  — stated in one sentence; the article focuses on daytime ΔTmax (the heat-buffering question), ΔTmin
  kept in annex.

### R4 — Mechanism: vertical gradient (illustration)
→ **Fig 6** (`fig_vertical_profiles_wind_rh_vpd_V9vsV11`) ⚠️ mechanistic illustration, no vertical field reference.

## 4. Discussion
- **D1 — Why leaf quantity & cover dominate and not the profile.** Near-saturated interception/extinction
  in dense oak (driven by total leaf area and cover, not its vertical placement); consistent with Bouwen's
  density threshold; **Canopy Ratio (Starck) adds nothing** beyond LAI (Annex). Nuances the "3-D
  indispensable" dogma.
- **D1bis — The theory ↔ field gap (closes the "does the vertical matter?" debate).**
  *Field result* (Annex A2): controlling LAI+Hmax+fCover, profile center-of-mass (COM) has **no effect**
  on observed buffering (partial r = −0.04, p = 0.79) — consistent with φ_LAD ≈ 0. The naive r = −0.43
  (controlling LAI only) was COM proxying fCover (r = 0.62) and Hmax (r = 0.47).
  *Mechanistic explanation* (Annex A3, in silico): H2's mechanism **exists** (top-heavy cools more) but
  is **small and saturates** with LAI. Signal/noise: model sensitivity ≈ 0.5–1.2 °C per COM unit, but the
  COM variation actually available at fixed structure is narrow (residual SD ≈ 0.09 vs the 0.50 span
  imposed in silico) → expected field effect ≈ **0.05–0.10 °C**, ~**10× below the between-plot noise
  (SD ≈ 0.65 °C)** → not identifiable. **No model↔reality disagreement**: the profile modulates ΔTmax in
  the right direction but, over the real architectural variability, an order of magnitude below the
  inter-plot signal → not an operative lever in real stands.
- **D2 — Glass ceiling.** Ranking good, absolute level not; the +0.9 to +2 °C warm bias is **most likely
  a forcing artefact (ERA5)** — the regional reanalysis above-canopy input, not a process error in the
  canopy scheme. (Secondary candidates: soil moisture, sensor radiative error.) To be revisited if the
  forcing is corrected.
- **D2bis — What MuSICA cannot see (topography & lateral heterogeneity).** MuSICA is a **plot-scale 1-D**
  model: slope, aspect, elevation, cold-air pooling, lateral advection, edge effects and soil-moisture
  heterogeneity are **not inputs and not simulated**. Key implication: these do **not** bias the
  attribution (they are absent from the 16 coalitions — φ stays valid *within the MuSICA world*); they
  are confined to the **validation residual**, i.e. they are prime suspects for the warm bias and the
  ~8 unreproduced amplifying plots (edges / locally exposed positions a 1-D model cannot generate).
  Scope mitigation: **Blois lies in the low-relief Loire plain**, so topographic forcing is minor.
  **Quantified (Annex A8a, `scripts/make_residual_vs_topo.R`):** regressing the per-plot validation
  residual (HOBO − model ΔTmax) on NC_Full topo/heterogeneity covariates (elevation, slope, TWI,
  northness, height-SD, rumple, gap fraction) explains **R² = 0.51 (p = 3e-4)** of the residual —
  i.e. **about half the glass ceiling is structured by what MuSICA cannot see**, not noise. The
  dominant term is **canopy gap fraction** (r = 0.58, p < 0.001; the only strong term in the multiple
  model), with secondary topographic exposure (elevation r = −0.29, slope +0.31, northness −0.31). The
  **8 amplifying plots are the gappy / steeper / more exposed ones** (gap 0.050 vs 0.002; slope 4.9 vs
  1.7; northness −1.0 vs −0.3). → the unreproduced amplification is **sub-pixel canopy openness +
  topographic exposure**, both absent from a homogeneous 1-D plot model.
  *Non-circularity (fCover vs gap fraction):* these are **not the same variable** — fCover is the
  **plot-mean cover** MuSICA actually ingests as a single scalar per plot (and which the Shapley
  attributes), whereas gap fraction here is the **within-pixel / sub-plot heterogeneity** that a single
  plot-mean fCover cannot represent. So gap fraction explaining the *residual* is not the attribution
  leaking back in: it is precisely the lateral openness the 1-D plot model is blind to, by construction.
  *Honest caveat:* gap fraction is zero-inflated (most plots closed-canopy), so the relationship is
  leveraged by the few gappy/amplifying plots; n ≈ 44 (9 plots outside raster coverage).
- **D3 — Methodological contribution.** Exact Shapley vs GAMM (resolves concurvity, replaces the
  black-box "Date" with physical ERA5 forcings); averages over **all** coalitions; robust under
  collinearity; conditional Shapley bounds the off-manifold artefact on **both** correlated pairs — φ_LAD ≈ 0
  holds on the coherent Hmax-real subset (+0.022 vs +0.023 off-manifold), and φ_LAI / φ_fCover are stable
  whether the collinear partner is real or baseline (gaps 0.022 / 0.010 °C, A4b), so neither the φ_LAD ≈ 0
  result nor the LAI/fCover co-dominance is a chimeric-coalition artefact. The **−0.43 → −0.04** field contrast is a direct
  demonstration of why full-control Shapley beats naive partial correlations / GAMM partial effects.
  Limit: attribution is **model-internal, spatially validated**.
- **D3bis — Why a mechanistic model (MuSICA) rather than a statistical one.** Not a horse-race against a
  regression: the rationale is **physical** (energy/radiative/turbulent balance, not a fitted slope or an
  equilibrium statistical model) and **generalisable** in principle (transferable parameters), which a
  site-fitted statistical model is not — even though we validate on a single site here. Statistical
  microclimate models (linear offset/sensitivity "slope", equilibrium models) describe; MuSICA lets us
  *attribute* to physical canopy traits and, downstream (Ch. 3), *spatialise* dynamically.
- **D4 — Robustness & scope.** Trait ranking **identical for v3.2.0 and v3.2.3** (Annex A6); v3.2.0 keeps
  the spatial structure (r 0.93 vs 0.63), v3.2.3 only trims the warm bias. Single-site scope restated.
  **Pseudoreplication / spatial structure (Moran's I):** the *observed* buffering slope is strongly
  spatially autocorrelated (I = 0.30, p < 0.001 kNN=5) — real forest structure is spatially organised;
  the *model* slope reproduces this spatial structure but **dampened** (I = 0.13, p = 0.02) — another
  facet of the compression/glass ceiling; and crucially the **validation residual is NOT autocorrelated**
  (I = 0.09, p = 0.07 kNN=5; I = 0.015, p = 0.29 inverse-distance) → the model captures the spatial
  pattern well enough that the residual is effectively independent, so the per-plot inference is not
  inflated by spatial clustering. cLHS sampling along the structural gradient further mitigates.
- **D5 — Perspective (Chapter 2).** Sentinel-2 optical illusion (Annex); the naive all-satellite run
  (S2 LAI + FORMS-H height → MuSICA) is a **Chapter-3 baseline** (cross-reference).

## 5. Conclusion
Leaf quantity and cover (LAI, fCover) are the co-dominant levers of simulated summer buffering — fCover a
near-constant floor, LAI the between-plot differentiator — while the **vertical profile adds little**
(a strong nuance to H1/H2). The 3-D structure distributes the between-plot variance well (r ≈ 0.93) but
the model plateaus on the absolute level and the amplifiers → an honest "glass ceiling", opening Chapter 2.

---

## Figures (numbered set: `outputs/figures_article_z05/`)
- **Fig 1** typology — LAD profiles of P1–P4
- **Fig 2** Shapley ranking *(key)* — LAI ≈ fCover ≫ Hmax > LAD ≈ 0
- **Fig 3** Shapley by cluster — redistribution along the density gradient
- **Fig 4** validation / glass ceiling — REF vs HOBO, r ≈ 0.93
- **Fig 5** forward validation in Shapley order
- **Fig 6** vertical gradient (illustration)

## Annexes
*(consolidated: A6+A7 merged → A6a/b/c; A10+A11 merged → A8a/b; old "LAD z0.5 vs z1" dropped)*
- **A1** FPCA of doubly-normalised LAD shape (→ clustering) + k-means elbow + cluster map (A1a–d)
- **A2** H2 field validation (COM vs buffering; null at full control)
- **A3** H2 controlled in-silico (top-heavy cools more; saturates) — supports D1bis
- **A4** conditional Shapley for LAD | Hmax (off-manifold bound, structural co-definition)
- **A4b** conditional Shapley for LAI | fCover & fCover | LAI (off-manifold bound, collinearity r=0.76; gaps ≤0.02 °C)
- **A5** extra metrics (ΔTmin / amplitude / stability) — backs the one-line Tmin statement
- **A6** version robustness v3.2.0 vs v3.2.3 (a scatter, b heatmap, c sim-vs-obs)
- **A7** Sentinel-2 / optical illusion = saturation (Chapter-2 bridge)
- **A8** glass ceiling decomposition (a topography R²=0.51, b slope buff/amp) — read jointly
- **A9** sampling-height robustness (nair==1 vs fixed 1 m) — justifies the two extraction conventions
- **A10** trait collinearity matrix (LAI/Hmax/fCover/COM/VCI) — underpins the Shapley rationale & COM↔fCover

## Reproducibility & data availability
**Pinned z05 expected values** (the article branch runs with `strict = FALSE`, so these asserts are NOT
machine-enforced — pin them here as the reference any rerun must reproduce):
- spatial validation **r = 0.931** (ΔTmax REF vs HOBO, between-plot)
- **Σφ = −0.378 °C** (additivity: metric(1111) − metric(0000))
- micro–macro slope split **46/7** (model z05) — distinct from the **45/8** field-HOBO split (cite each to
  its own source; never conflate the model split with the field split)
- φ-ranking **{LAI, fCover} ≫ Hmax > LAD ≈ 0** (φ_LAD ns)

**Sample sizes** — n varies by analysis (53 HOBO total; 49/44 where an analysis needs the NC_Full oak
raster): plots falling outside the NC_Full oak raster coverage are dropped from the topo/heterogeneity
analyses. State this once and cite the relevant n at each figure.

**Data availability** — the core pipeline (FPCA → clustering → MuSICA coalitions → Shapley → HOBO
validation) is self-contained in-repo. External dependencies to flag: **A8a (topo residual), A7 (S2), A10
(collinearity) and the topo covariates require the external `~/Documents/NC_Full` dataset**; the **frozen
cLHS design `.rds` has no in-repo generator** (provided as a frozen artefact). State both in the data
statement so a reader knows what is and isn't regenerable from the repo alone.

## Reviewer-2 pre-emption checklist
- [ ] Single-site/year/species scope stated early.
- [ ] Validation headline = **r ≈ 0.93 (spatial ranking)**; +0.9 °C bias + 8 amplifiers = the named ceiling (not r = 0.67).
- [ ] Attribution framed as **model-internal, spatially validated** (Results lead sentence).
- [ ] Method = exact Shapley throughout (GAMM = removed V1 lock).
- [ ] Combined **leaf-quantity/cover factor** dominates, LAD ≈ 0 — hypotheses reframed, not claimed confirmed.
- [ ] **LAI and fCover co-dominant, collinearity controlled** (A4b conditional Shapley: off-manifold gap
      ≤ 0.02 °C, within CI) — reported as co-dominant, NOT ranked (CIs overlap); the split is robust, not labile.
- [ ] **fCover (attributed, plot-mean) ≠ gap fraction (residual, sub-plot heterogeneity)** — stated to
      pre-empt the circularity objection in D2bis.
- [ ] **45/8 = field-HOBO split; 46/7 = z05-model split** — cited to source, never conflated.
- [ ] **n = 53 / 49 / 44** variation explained once (plots outside NC_Full oak raster coverage).
- [ ] Pinned z05 values (r = 0.931, Σφ = −0.378, ranking) recorded since branch runs strict = FALSE.
- [ ] LAI ×2 = two-sided leaf area (never "clumping"); display LAI shown one-sided (÷2).
- [ ] Vertical gradient = illustration; ET/forcing(ERA5) = one hypothesis among several; cluster↔LAI confound acknowledged.
- [ ] Slope mismatch split correctly: closed-canopy under-buffering = structural model limit; amplifier
      mismatch = topo/gap scope limit (A8b read with A8a, not as independent model failure).
- [ ] "autocorr. uncorrected" caveat kept in per-cluster significance prose; Moran obs-slope I = 0.30 sig.
- [ ] One-line caveats present: n ≈ 7–8 amplifiers (modest inference); −2 h shift assumption.
