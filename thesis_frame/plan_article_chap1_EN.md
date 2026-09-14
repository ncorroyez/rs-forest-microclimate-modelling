# Leaf quantity governs understory microclimate, but vertical arrangement takes over where the canopy saturates: a LiDAR–MuSICA attribution

## Abstract
Forest understories buffer macroclimatic extremes, yet which dimension of canopy structure drives this
buffering — total leaf quantity versus the vertical arrangement of foliage — remains debated. We combine
LiDAR-derived canopy traits (LAI, maximum height Hmax, fractional cover fCover, and the vertical
leaf-area-density profile) with the MuSICA biophysical canopy model to attribute simulated summer understory
microclimate (ΔTmax, the diurnal micro–macro temperature offset) at a deciduous oak forest (Blois, France;
summer 2021; 53 HOBO loggers). Across **400 LiDAR plots stratified into four structural archetypes** (100
each), we compute an **exact Shapley** attribution of ΔTmax over all 2⁴ trait combinations, toggling each
trait against its **own archetype's realistic canopy**. The dominant driver is **density-dependent**: leaf
quantity (LAI) governs microclimate variation in **open** canopy (|φ| up to 0.41 °C), but where the canopy is
**dense** and leaf quantity saturates, the **vertical profile (LAD) becomes the leading driver** (|φ| ≈
0.13 °C while LAI and cover collapse to ~0). A common-baseline landscape view shows leaf quantity and cover
dominating with the vertical profile ≈ 0 *on average* — but that average **hides the dense-canopy regime**
that the realistic per-archetype baseline reveals; the vertical profile is a small absolute lever everywhere,
leading in dense canopy only because quantity saturates. MuSICA reproduces the between-plot buffering ranking
well (r ≈ 0.93) but hits a glass ceiling on absolute level (+0.9 °C warm bias) and on topographic / sub-pixel
canopy-gap heterogeneity a 1-D plot model cannot represent. This density-dependent control connects to the
canopy-density threshold of neighbouring work and motivates the optical-saturation analysis of Chapter 2.

## 1. Introduction
- Forest microclimate is decoupled from the regional macroclimate and is the climate organisms actually
  experience: understory conditions can differ from open-field / regional weather by several degrees, with
  direct consequences for forest biota [REF De Frenne et al.; Lembrechts et al.; Zellweger et al.].
- This buffering shapes forest ecosystem functioning — it underpins thermal refugia, tree regeneration,
  understory biodiversity and species persistence under climate warming [REF].
- Canopy structure is the primary driver of buffering [REF]. Traditionally characterised from field
  inventories, it can now be derived from remote sensing — in particular airborne laser scanning (ALS/LiDAR),
  which resolves the full 3-D vegetation profile — and ALS-derived structure has been used in
  forest-microclimate studies [REF Zellweger; Jucker]. In a neighbouring LiDAR + MuSICA study, Bouwen [REF]
  reported a canopy-density threshold on buffering, leaving open whether the same holds in deciduous stands
  and, more fundamentally, **which structural dimension actually drives the effect**.
- But, more precisely, which structural variables contribute most to buffering remains unclear. LiDAR is
  valued *precisely* for resolving the vertical profile, and the vertical arrangement of foliage is widely
  assumed decisive — **yet whether the vertical dimension itself matters, or merely co-varies with the total
  quantity of leaves (LAI) and cover, has never been disentangled.** This is our core question: leaf quantity
  vs vertical arrangement.
- Statistical modelling of microclimate (linear offset / sensitivity "slope", equilibrium models, GAMMs) is
  data-cheap, predictive and widely used [REF] — but it cannot cleanly *attribute* buffering to individual
  structural traits: canopy traits are strongly collinear, so regression coefficients are unstable
  (concurvity) and the effect of any one trait cannot be isolated observationally.
- A physics-based proposition addresses exactly this: a mechanistic multilayer canopy model (MuSICA [REF])
  solving the energy / radiative / turbulent balance can simulate **controlled** trait combinations,
  attributing buffering to physical traits and disentangling correlated variables that a regression cannot
  separate — and being process-based, it is transferable in principle.
- In this study, we couple ALS-derived canopy traits with MuSICA and apply an exact-Shapley attribution over
  all trait combinations, validated against 53 HOBO loggers at Blois, to quantify the relative contribution
  of leaf quantity vs vertical arrangement to summer understory buffering.

We hypothesised that:
- **(H1)** vertical structure modulates radiative / turbulent transfers, so the LAD profile is expected to be
  a key parameter for faithful buffering simulation;
- **(H2)** at fixed LAI and height, top-heavy profiles (biomass concentrated high) cause stronger diurnal
  attenuation (radiation intercepted higher, heat sources shifted to the canopy top).

To achieve our objective and test these hypotheses, we: (1) derived a structural typology of the LAD profile
(FPCA + k-means → archetypes P1–P4); (2) simulated the 2⁴ trait coalitions per plot with MuSICA and computed
an exact-Shapley attribution of ΔTmax, with an off-manifold conditional control for correlated traits; and
(3) validated the simulated buffering against the HOBO network and decomposed the residual (the glass ceiling)
against topographic / heterogeneity covariates that MuSICA cannot see.

As a perspective toward Chapter 2, we further ask whether a 2-D optical proxy (Sentinel-2) can substitute for
LiDAR structure or merely gives an optical *illusion* of buffering (a structure proxy).

## 2. Material and Methods
### 2.1 Study site and in-situ microclimate
- Blois is a temperate lowland deciduous forest located in the centre of France, almost entirely composed of
  sessile oaks (*Quercus petraea*).
- 53 HOBO loggers, distributed along a canopy light / structure gradient, record hourly air-temperature
  measurements at 1 m above the ground (summer 2021).
- The macroclimate reference is open-field / regional weather; ΔTmax = diurnal micro − macro temperature
  offset (the buffering metric).

### 2.2 LiDAR canopy structure and traits
- 0.5 m vertical leaf-area-density (LAD) profiles from ALS (resolving the understory); traits: LAI, maximum
  height Hmax, fractional cover fCover, VCI, and profile shape.
- LiDAR LAI is one-sided (projected) → **×2 = two-sided leaf area** expected by MuSICA (not a clumping
  correction).

### 2.3 Structural typology
- FPCA of the **doubly-normalised LAD profile** (normalised in height *x* AND density *y* → pure shape,
  scale-free of LAI / Hmax); K-means on the FPCA scores → **four archetypes P1–P4** (architectural-shape-based,
  independent of leaf amount). Clustering frozen; cLHS design.

### 2.4 MuSICA model and simulation design
- MuSICA is a multilayer biophysical canopy model; it ingests **LAI as two-sided leaf area distributed per
  layer** and the **LAD profile as the vertical discretisation** of that leaf area.
- **Design = per-point cLHS, stratified by archetype**: a conditioned Latin-hypercube sample of **400 plots
  (100 per archetype P1–P4)** spanning the real trait space; **16 coalitions (2⁴) per plot** (each trait
  real vs baseline). 6 400 MuSICA runs; legacy binary **v3.2.0**.
- **Two baseline schemes (the choice IS the question):**
  - *Main — cluster-specific baseline:* each trait toggled against **its own archetype's mean** (LAD→uniform).
    Realistic / on-manifold per type; coalition 0000 = that archetype's canopy. → within-archetype drivers.
  - *Complement — common (global) baseline:* LAI→mean (LAI_b ≈ 6.73), Hmax→mean, fCover→mean (0.869),
    LAD→uniform. → landscape attribution (but off-manifold for atypical types — see R2).
- *Why per-point, not the 4 mean profiles:* MuSICA is nonlinear, so the Shapley of a cluster's **mean
  profile** is Jensen-biased (f(mean) ≠ mean(f)); per-point + aggregation is unbiased and yields the
  within-archetype distribution (the Jensen gap is itself reported: largest in the sparse P1).

### 2.5 Microclimate metric extraction (two conventions)
- (a) absolute ΔTmax / ΔVPDmax → fixed 1 m interpolation, −2 h shift (HOBO-comparable);
- (b) micro–macro slope (buffering / amplification split) → nair==1, no shift.

### 2.6 Shapley attribution (all combinations of traits)
- **Exact Shapley** φ per trait per plot (ΔTmax, ΔVPDmax), additivity verified (Σφ = full − base, to machine
  precision), bootstrap CIs; isolates φ_LAD vs φ_LAI in place of the GAMM.
- **Main readout (cluster-specific baseline):** mean φ ≈ 0 by design, so importance = mean **|φ|** per
  archetype → within-archetype trait importance.
- **Complement (common baseline):** signed φ → landscape ranking (report median / |φ|, **not** the
  design-pooled signed mean, which the broad cLHS tail skews).
- **Conditional Shapley** bounds the off-manifold artefact of correlated traits on the common-baseline run:
  φ_LAD | Hmax (structural co-definition) and φ_LAI | fCover & φ_fCover | LAI (statistical collinearity
  r = 0.76) — gaps ≤ 0.02 °C.

### 2.7 Validation and spatial-heterogeneity decomposition
- REF vs 53 HOBO; daily ΔTmax; micro–macro slope (45/8 split); Moran's I for spatial autocorrelation.
- ΔVPDmax is attributed but not validated (HOBO measure air temperature only); VPD results are model-internal.
- The validation residual (HOBO − model ΔTmax) is regressed on topographic / heterogeneity covariates from an
  external LiDAR dataset (elevation/DTM, slope, TWI, northness, canopy height-SD, rumple, gap fraction) —
  variables a 1-D plot model cannot represent.

## 3. Results

### R1 — A four-archetype structural typology
Four archetypes from open **P1** to dense **P4** (**Fig 1**; + `tab_cluster_structure`). Although the
typology is built on profile *shape* (FPCA of the doubly-normalised LAD), the four archetypes line up along a
**canopy-density / LAI gradient** (P1 sparsest → P4 densest). The typology is therefore a presentation device
that organises this gradient and lets us show *how the attribution redistributes along it* (R2, Fig 3) — it is
**not** the engine of the attribution, which is computed per plot, independently of the clustering.

### R2 — Which trait governs understory microclimate is density-dependent *(core)*

**Main analysis — within-archetype drivers (per-point cLHS, cluster-specific baseline; Fig 2).** For each
archetype we run an exact Shapley on **100 cLHS plots**, toggling each trait against **its own cluster's
realistic mean** (not a global mean, which would impose off-manifold canopies on atypical types — e.g. a
landscape-mean dense canopy on a sparse-tree cluster). Mean φ ≈ 0 by design, so importance = mean **|φ|**.
The leading driver of plot-to-plot ΔTmax variation **shifts with canopy density**:
- **Open canopy (P1):** LAI dominates (|φ| = 0.41 °C), with fCover (0.21), LAD (0.17) and Hmax (0.13) all
  contributing — the most variable type.
- **Intermediate (P2, P3):** LAI leads (0.20, 0.15); others minor.
- **Dense canopy (P4):** LAI and fCover **saturate** (|φ| = 0.10 and 0.01) → the **vertical profile (LAD)
  becomes the leading driver** (|φ| = 0.13 °C). Where leaf quantity can no longer vary the microclimate,
  *arrangement* governs the residual.
→ **Leaf quantity drives the open end of the gradient; the vertical profile governs the dense end where
quantity saturates** — a density-dependent control that connects directly to Bouwen's density threshold.

**Complement — landscape attribution against a common baseline (Fig 3, Annex).** With a *single* baseline
across the whole landscape, leaf quantity and cover dominate and φ_LAD ≈ 0 *on average* (φ_LAI ≈ φ_fCover
≈ −0.2 °C; Σφ ≈ −0.38; same order for VPD). This is what a single-baseline Shapley (or the V1 GAMM) reports,
and the conditional/off-manifold control (A4/A4b) confirms the collinearity (LAI↔fCover r = 0.76) does not
distort this global ranking (gaps ≤ 0.02 °C). But this common baseline is **off-manifold for atypical types**
and averages the density structure away: only the realistic per-archetype baseline reveals that
**"LAD ≈ 0" is a landscape average, not a universal — it hides the dense-canopy regime where the profile leads.**

**Reconciling with H2 (no contradiction).** The vertical profile is a **small absolute lever everywhere**
(|φ_LAD| ≤ 0.17 °C; H2's top-heavy mechanism is real but second-order and saturates — Annex A2/A3, field
partial r = −0.04 at full control). It "dominates" in dense canopy only **relatively**: LAI and fCover
collapse to ~0 leverage there, so the small vertical effect is what remains. **Small-but-leading, not large.**

### R3 — The 3-D distributes the variance, but a glass ceiling against the field
Between-plot ranking captured (**r ≈ 0.93**; **Fig 4**); warm **bias +0.9 °C**. With the slope now at the same
**1 m** as the loggers, the model reproduces the buffering/amplifying split (**44/9 vs observed 45/8**), but
compresses the slope **magnitude** toward 1 (median sim **0.95** vs obs **0.86**) — it under-expresses both
buffering and amplification rather than mis-counting plots. The glass ceiling is therefore (i) the slope
**amplitude compression** on the closed-canopy plots (gap ≈ 0, topo negligible) — the genuine model limit —
and (ii) the residual deviations, which are topo/gap-driven and out of scope for a 1-D model.

**Decomposing the residual (what the model misses; Fig 5, Annex A8).** Regressing the per-plot validation
residual (HOBO − model ΔTmax) on topographic / heterogeneity covariates explains **R² = 0.51 (p = 3e-4)** —
about half the glass ceiling is structured, not noise. The dominant term is **canopy gap fraction**
(r = 0.58, p < 0.001), with secondary topographic exposure (elevation/DTM r = −0.29, slope +0.31, northness
−0.31; all p ≈ 0.04). Crucially, **vertical canopy heterogeneity does not contribute** — CHM height-SD
(r = 0.04, ns) and rumple (r = 0.17, ns) are null. The 8 amplifying plots are the gappy/steeper/more exposed
ones (gap 0.050 vs 0.002; slope 4.9 vs 1.7). So the unreproduced amplification is **sub-pixel canopy openness
+ topographic exposure** — horizontal, not vertical. The attribution and buffering hierarchy also hold for
night-time ΔTmin (one sentence; daytime ΔTmax is the focus).

## 4. Discussion
### D1 — A density-dependent control: quantity in the open, arrangement in the dense
The leading driver shifts with canopy density (R2). In open canopy, near-unsaturated interception means
**total leaf quantity** sets the microclimate — adding leaves still changes radiative/turbulent transfer, so
LAI dominates. As the canopy closes, interception/extinction **saturate**: extra leaf area and cover no longer
move ΔTmax (|φ_LAI|, |φ_fCover| → ~0 in P4), so the remaining variation is governed by **where** the foliage
sits — the **vertical profile**. This is exactly Bouwen's density threshold expressed at the trait level:
below threshold, quantity; above it, arrangement. Crucially the vertical effect is **small in absolute terms
everywhere** (|φ_LAD| ≤ 0.17 °C; Canopy-Ratio adds nothing beyond LAI) — it *leads* in dense canopy by
saturation of the others, not by being large. The vertical gradient of wind / RH / VPD
(**Fig 6**, mechanistic illustration — no vertical field reference) shows the transfer mechanism operating but
saturating. H2 is the clean illustration of this: its mechanism is real (top-heavy cools more, in silico) yet
operationally negligible — the COM variation available at fixed structure is narrow (residual SD ≈ 0.09 vs the
0.50 span imposed in silico), so the expected field effect (≈ 0.05–0.10 °C) is ~10× below the between-plot
noise (SD ≈ 0.65 °C) and undetectable. The naive field signal (COM partial r = −0.43 controlling LAI only)
was COM proxying fCover (r = 0.62) and Hmax (r = 0.47); it vanishes at full control (−0.04). No model ↔ reality
disagreement: the profile acts in the right direction but is not an operative lever in real stands.

### D2 — The glass ceiling and what a 1-D model cannot see
MuSICA captures the between-plot ranking (r ≈ 0.93) but not the absolute level: the +0.9 to +2 °C warm bias is
most likely a **forcing artefact** (ERA5 above-canopy reanalysis input), not a process error in the canopy
scheme (secondary candidates: soil moisture, sensor radiative error); to be revisited if the forcing is
corrected. MuSICA is a plot-scale **1-D** model — slope, aspect, elevation, cold-air pooling, lateral
advection, edge effects and soil-moisture heterogeneity are not inputs. These do **not** bias the attribution
(absent from the 16 coalitions — φ stays valid within the MuSICA world); they are confined to the validation
residual, which is exactly why the residual decomposition (R3: R² = 0.51, gap-fraction-dominated, vertical
heterogeneity null) localises the missing amplification in **horizontal openness + topographic exposure**.
This is the spatial-heterogeneity mirror of φ_LAD ≈ 0: the vertical dimension is inoperative both inside the
model (attribution) and outside it (residual). Note fCover (the plot-mean cover MuSICA ingests and the Shapley
attributes) and gap fraction (sub-pixel heterogeneity a single plot-mean cannot represent) are different
variables, so gap fraction explaining the residual is not the attribution leaking back in. Spatial structure
is genuine but does not inflate the inference (Moran's I): the observed buffering slope is autocorrelated
(I = 0.30, p < 0.001), the model slope reproduces this but dampened (I = 0.13, p = 0.02 — another facet of the
ceiling), and crucially the **residual is not autocorrelated** (I = 0.09, p = 0.07), so per-plot inference is
not pseudoreplicated; cLHS sampling along the structural gradient further mitigates.

### D3 — Methodological contribution
Exact Shapley supersedes the GAMM (resolves concurvity, replaces the black-box "Date" with physical ERA5
forcings); it averages over all coalitions and is robust under collinearity, and conditional Shapley bounds
the off-manifold artefact on both correlated pairs — φ_LAD ≈ 0 holds on the coherent Hmax-real subset
(+0.022 vs +0.023 off-manifold), and φ_LAI / φ_fCover are stable whether the collinear partner is real or
baseline (gaps 0.022 / 0.010 °C), so neither result is a chimeric-coalition artefact. The −0.43 → −0.04 field
contrast demonstrates why full-control Shapley beats naive partial correlations / GAMM partial effects. The
deeper rationale for a mechanistic model is physical (energy/radiative/turbulent balance, not a fitted slope
or equilibrium model) and generalisable in principle (transferable parameters): statistical microclimate
models describe; MuSICA lets us **attribute** to physical traits and, downstream (Ch. 3), **spatialise**
dynamically. The attribution is model-internal but spatially validated.

### D4 — Limitations and perspectives
Single-site, single-year, single-species scope. Trait ranking is identical for v3.2.0 and v3.2.3 (v3.2.0
keeps the spatial structure, r 0.93 vs 0.63; v3.2.3 only trims the warm bias), so the conclusions are robust
to the model version. Toward Chapter 2: the Sentinel-2 optical signal saturates and gives an optical
illusion of buffering, and the naive all-satellite run (S2 LAI + FORMS-H height → MuSICA) is a Chapter-3
baseline.

## 5. Conclusion
Which canopy trait governs understory microclimate is **density-dependent**: leaf quantity (LAI) drives the
variation in open canopy, but where the canopy is dense and quantity saturates, the **vertical profile takes
over** as the leading driver. The vertical dimension is a small absolute lever everywhere — its prominence in
dense stands is a saturation effect, not a large effect — which both nuances the "3-D is indispensable" dogma
and explains why a landscape average reads "vertical ≈ 0". MuSICA reproduces the between-plot buffering
ranking well (r ≈ 0.93) but plateaus on absolute level and amplifiers → an honest glass ceiling. The
density-dependent control connects to the canopy-density threshold of neighbouring work and opens Chapter 2.

## Figures
- **Fig 1** typology — LAD profiles of P1–P4 (`fig_archetypes_profiles`, LAI one-sided)
- **Fig 2** *(core)* within-archetype drivers — |φ| per trait per archetype, cluster-specific baseline
  (`FigSh_shapley_clhs_cluster`); LAI leads open P1, LAD leads dense P4 → density-dependent
- **Fig 3** landscape attribution (common baseline) — leaf quantity & cover dominate, LAD ≈ 0 on average
  (`FigSh_shapley_clhs_global` / z05 ranking); the average that Fig 2 refines
- **Fig 4** validation / glass ceiling — REF vs HOBO, r ≈ 0.93
- **Fig 5** residual decomposition — residual vs gap fraction / topo (R² = 0.51)
- **Fig 6** vertical gradient — mechanistic illustration (cited in D1)

*(forward-validation-in-Shapley-order figure → Annex; two validation figures in the main text would
over-weight validation for an attribution paper.)*

## Annexes
- **A1** FPCA of doubly-normalised LAD shape (→ clustering) + k-means elbow + cluster map (A1a–d)
- **A2** H2 field validation (COM vs buffering; null at full control)
- **A3** H2 controlled in-silico (top-heavy cools more; saturates)
- **A4** conditional Shapley for LAD | Hmax (off-manifold bound, structural co-definition)
- **A4b** conditional Shapley for LAI | fCover & fCover | LAI (off-manifold bound, collinearity r = 0.76; gaps ≤ 0.02 °C)
- **A4c** Shapley vs leave-one-out vs forward (add-one) on ΔTmax (`FigA_shapley_vs_loo_forward`): Shapley is bracketed by the two; forward−LOO wedge widest for the collinear pair LAI/fCover, tight for Hmax/LAD; the modest wedge confirms near-additive co-dominance (ties to A4b) — justifies Shapley over LOO/selection
- **A5** extra metrics (ΔTmin / amplitude / stability)
- **A6** version robustness v3.2.0 vs v3.2.3 (a scatter, b heatmap, c sim-vs-obs)
- **A7** Sentinel-2 / optical illusion = saturation (Chapter-2 bridge)
- **A8** glass ceiling — slope buffering/amplifying panel (the topography R² = 0.51 panel is promoted to main Fig 5)
- **A9** sampling-height robustness (nair==1 vs fixed 1 m)
- **A10** trait collinearity matrix (LAI/Hmax/fCover/COM/VCI)
- **A11** forward validation in PER-CLUSTER addition order (+ All, pooled) — each archetype adds traits in its own |φ| order; LAD first in dense P4; n per cluster 8/16/15/14; |φ|-importance ≠ spatial skill (LAD modest alone in P4, LAI completes)
