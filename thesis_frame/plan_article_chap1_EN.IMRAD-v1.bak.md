# Leaf quantity and cover, not vertical arrangement, govern understory microclimate buffering in LiDAR-driven MuSICA simulations

## Abstract
Forest understories buffer macroclimatic extremes, yet which canopy structural traits drive this buffering
— total leaf quantity versus the vertical arrangement of foliage — remains debated. We combine LiDAR-derived
canopy traits (LAI, maximum height Hmax, fractional cover fCover, and the vertical leaf-area-density profile)
with the MuSICA biophysical canopy model to attribute simulated summer understory buffering (ΔTmax, the
diurnal micro–macro temperature offset) at a deciduous oak forest (Blois, France; summer 2021; 53 HOBO
loggers). For each plot we simulate the 2⁴ coalitions of traits (each real vs a baseline) and compute an
exact Shapley attribution of ΔTmax, validated spatially against the loggers. Leaf quantity and cover (LAI,
fCover) are co-dominant buffering traits (φ ≈ −0.2 °C each), whereas the vertical profile contributes
essentially nothing (φ_LAD ≈ 0). Their collinearity (r = 0.76) does not distort the attribution: an
off-manifold conditional Shapley shows each contribution is stable whether the correlated partner is present
or not (gap ≤ 0.02 °C). fCover acts as a near-constant buffering floor while LAI carries the between-plot
variation. MuSICA reproduces the between-plot ranking well (r ≈ 0.93) but compresses the amplitude and shows
a +0.9 °C warm bias — a glass ceiling: about half the validation residual is explained by topographic and
sub-pixel canopy-gap heterogeneity that a 1-D plot model cannot represent. Leaf quantity and cover, not
vertical arrangement, govern simulated summer buffering, nuancing the view that 3-D structure is
indispensable and motivating the optical-saturation analysis of Chapter 2.

## 1. Introduction
1. Forest microclimate decoupling and thermal refugia under warming (De Frenne et al.; Lembrechts et al.;
   Zellweger et al.).
2. Canopy structure as the driver — but **which dimension**: total leaf area vs the **3-D vertical
   profile**? Position relative to **Bouwen** (neighbouring LiDAR+MuSICA work; density threshold).
3. The limits of a monolithic GAMM with **trait collinearity / concurvity** and a black-box "Date" effect.
   This motivates an attribution that (i) explores **all trait combinations**, (ii) is robust under
   collinearity, and (iii) avoids off-manifold partial-dependence plots.
4. Our approach: **MuSICA scenarios + exact Shapley** on the 2⁴ trait lattice, validated against 53 HOBO.

### Research questions
1. How sensitive is the MuSICA-simulated understory temperature to the **vertical and horizontal**
   distribution of vegetation?
2. What is the **relative contribution of the vertical profile (LAD) vs total leaf area (LAI)** to ΔTmax?

### Hypotheses
- **H1** — Vertical structure modulates radiative/turbulent transfers, so the LAD profile is a key
  parameter for faithful microclimate simulation.
- **H2** — At fixed LAI and height, **top-heavy** profiles (biomass concentrated high) cause stronger
  diurnal attenuation (radiation intercepted higher, heat sources shifted to the canopy top).
- **H3** — The 2-D optical signal (Sentinel-2) gives an *illusion* of buffering (a structure proxy)
  → motivates Chapter 2.

## 2. Materials & Methods
- **2.1 Site & sensors** — Blois oak forest, summer 2021, 53 HOBO at ~1 m. Single site, single year,
  single species (scope).
- **2.2 LiDAR & traits** — 0.5 m LAD profiles (resolving the understory); traits LAI, Hmax, fCover, **VCI**,
  and profile shape; LiDAR LAI is one-sided (projected) → **×2 = two-sided** leaf area expected by MuSICA.
- **2.3 Typology** — FPCA of the **doubly-normalised LAD profile** (normalised in height *x* AND density *y*
  → pure shape, scale-free of LAI/Hmax); K-means on the FPCA scores → **four archetypes P1–P4**. The
  typology is architectural-shape-based, independent of leaf amount. Clustering frozen; cLHS design.
- **2.4 MuSICA & simulation design** — **16 coalitions (2⁴)** per plot (each trait real vs baseline:
  LAI→mean, Hmax→mean, fCover→mean, LAD→uniform); fCover baseline = mean (0.869); cLHS LAI_b ≈ 6.73;
  legacy binary **v3.2.0**.
- **2.5 Extraction — two conventions**: (a) absolute ΔTmax/ΔVPDmax → fixed 1 m interpolation, −2 h shift
  (HOBO-comparable); (b) micro–macro slope (buffering/amplification split) → nair==1, no shift.
- **2.6 Attribution** — **exact Shapley** φ per trait per plot (ΔTmax, ΔVPDmax), additivity verified,
  bootstrap CIs; this isolates φ_LAD vs φ_LAI in place of the GAMM. **Conditional Shapley** bounds the
  off-manifold artefact of correlated traits on both co-defined/collinear pairs: φ_LAD | Hmax (structural
  co-definition) and φ_LAI | fCover & φ_fCover | LAI (statistical collinearity r = 0.76) — each trait's φ
  recomputed conditional on the partner being real vs baseline.
- **2.7 Validation** — REF vs 53 HOBO; daily ΔTmax; micro–macro slope (45/8 split). ΔVPDmax is attributed
  but not validated (HOBO measure air temperature only); VPD results are model-internal.

## 3. Results

### R1 — A four-archetype structural typology
Four archetypes from open **P1** to dense **P4** (**Fig 1**; + `tab_cluster_structure`).

### R2 — Leaf quantity & cover are co-dominant; vertical arrangement adds little
Ranking **LAI ≈ fCover ≫ Hmax > LAD ≈ 0** (**Fig 2**, **Fig 3**); φ_LAI(ΔTmax) ≈ −0.2 °C, φ_fCover ≈ −0.2,
φ_Hmax small, **φ_LAD ≈ 0 (ns)**; Σφ ≈ −0.38 °C (additivity holds); same order for VPD.

Collinearity does not drive the split (conditional Shapley LAI|fCover & fCover|LAI): each trait's
attribution is stable whether the correlated partner is real or stripped to baseline — φ_LAI|fCover=real
−0.204 vs |baseline −0.227 (gap 0.022 °C); φ_fCover|LAI=real −0.205 vs |baseline −0.215 (gap 0.010 °C).
Both gaps are an order of magnitude below the effect size (≈ −0.20 °C) and inside the ±0.1 °C bootstrap CI;
contributions are near-additive (no substitution). Same conclusion as the LAD|Hmax control — the two
conditionals agree — but a different mechanism (structural co-definition there, statistical collinearity here).

Mean-level co-dominance ≠ spatial co-dominance. Pooled medians match (φ_LAI ≈ φ_fCover ≈ −0.20 °C) but their
CIs do not: φ_fCover is tight ([−0.216, −0.206]) while φ_LAI is wide ([−0.302, −0.103]). By cluster (Fig 3),
**φ_fCover is a near-constant buffering floor** (P2 −0.21, P3 −0.22, P4 −0.21; only sparse P1 +0.03) whereas
**φ_LAI carries the between-plot variation** — a monotonic gradient with the sign-flip (P1 +0.15 amplify →
P4 −0.56 buffer). fCover sets a uniform floor in closed canopy; LAI differentiates plots. We report
mean-level co-dominance but do not rank LAI vs fCover (CIs overlap). At the ΔTmax scale, the vertical profile
adds essentially nothing beyond leaf quantity and cover.

### R3 — The 3-D distributes the variance, but a glass ceiling against the field
Between-plot ranking captured (**r ≈ 0.93**; **Fig 4**, **Fig 5**); warm **bias +0.9 °C**; field-HOBO split
45/8, z05-model split 46/7. Per-plot slope, MuSICA vs HOBO, splits buffering/amplifying: the model compresses
the slope range toward 1 — under-buffers the bufferers (sim 0.94 vs obs 0.84) and under-amplifies the
amplifiers (1.12 vs 1.27). The 8 amplifying plots are the gappy/steep/exposed ones whose amplification is
topo/heterogeneity-driven and not a MuSICA input, so that mismatch is a scope limitation; the only mismatch
attributable to the model is the slope compression on the closed-canopy buffering plots (gap ≈ 0, topo
negligible) — the genuine structural limit. The glass ceiling is therefore (i) real structural under-buffering
on closed canopy and (ii) an out-of-scope inability to reproduce topo/gap-driven amplification. The
attribution and buffering hierarchy also hold for night-time ΔTmin (one sentence; daytime ΔTmax is the focus).

### R4 — Mechanism: vertical gradient
Vertical profiles of wind / RH / VPD illustrate the transfer mechanism (**Fig 6**; mechanistic illustration,
no vertical field reference).

## 4. Discussion
- **D1 — Why leaf quantity & cover dominate and not the profile.** Near-saturated interception/extinction in
  dense oak (driven by total leaf area and cover, not its vertical placement); consistent with Bouwen's
  density threshold; Canopy Ratio (Starck) adds nothing beyond LAI. Nuances the "3-D indispensable" dogma.
- **D1bis — The theory ↔ field gap.** Field: controlling LAI+Hmax+fCover, profile center-of-mass (COM) has
  no effect on observed buffering (partial r = −0.04, p = 0.79) — consistent with φ_LAD ≈ 0; the naive
  r = −0.43 (controlling LAI only) was COM proxying fCover (r = 0.62) and Hmax (r = 0.47). In silico: H2's
  mechanism exists (top-heavy cools more) but is small and saturates with LAI — model sensitivity
  ≈ 0.5–1.2 °C per COM unit, but the COM variation available at fixed structure is narrow (residual SD ≈ 0.09
  vs the 0.50 span imposed in silico) → expected field effect ≈ 0.05–0.10 °C, ~10× below the between-plot
  noise (SD ≈ 0.65 °C) → not identifiable. No model↔reality disagreement: the profile modulates ΔTmax in the
  right direction but, over the real architectural variability, an order of magnitude below the inter-plot
  signal → not an operative lever in real stands.
- **D2 — Glass ceiling.** Ranking good, absolute level not; the +0.9 to +2 °C warm bias is most likely a
  forcing artefact (ERA5 above-canopy reanalysis input), not a process error in the canopy scheme.
  (Secondary candidates: soil moisture, sensor radiative error.) To be revisited if the forcing is corrected.
- **D2bis — What MuSICA cannot see (topography & lateral heterogeneity).** MuSICA is a plot-scale 1-D model:
  slope, aspect, elevation, cold-air pooling, lateral advection, edge effects and soil-moisture heterogeneity
  are not inputs and not simulated. They do not bias the attribution (absent from the 16 coalitions — φ stays
  valid within the MuSICA world); they are confined to the validation residual, so they are prime suspects
  for the warm bias and the ~8 unreproduced amplifying plots. Blois lies in the low-relief Loire plain, so
  topographic forcing is minor. Regressing the per-plot residual (HOBO − model ΔTmax) on topo/heterogeneity
  covariates (elevation, slope, TWI, northness, height-SD, rumple, gap fraction) explains R² = 0.51
  (p = 3e-4) — about half the glass ceiling is structured by what MuSICA cannot see, not noise. The dominant
  term is canopy gap fraction (r = 0.58, p < 0.001), with secondary topographic exposure (elevation/DTM
  r = −0.29, slope +0.31, northness −0.31, all p ≈ 0.04). Crucially, **vertical canopy heterogeneity does
  NOT contribute** — CHM height-SD (r = 0.04, p = 0.82) and rumple (r = 0.17, ns) are null, as is TWI
  (r = −0.10, ns). So the residual is driven by **horizontal openness + topographic exposure, not vertical
  structural heterogeneity** — the spatial-heterogeneity mirror of φ_LAD ≈ 0: the vertical dimension is
  inoperative both inside the model (attribution) and outside it (residual). The 8 amplifying
  plots are the gappy/steeper/more exposed ones (gap 0.050 vs 0.002; slope 4.9 vs 1.7; northness −1.0 vs
  −0.3) → the unreproduced amplification is sub-pixel canopy openness + topographic exposure, both absent from
  a homogeneous 1-D plot model. fCover (the plot-mean cover MuSICA ingests and the Shapley attributes) and gap
  fraction (within-pixel heterogeneity a single plot-mean cannot represent) are not the same variable, so gap
  fraction explaining the residual is not the attribution leaking back in. Gap fraction is zero-inflated, so
  the relationship is leveraged by the few gappy/amplifying plots (n ≈ 44).
- **D3 — Methodological contribution.** Exact Shapley vs GAMM (resolves concurvity, replaces the black-box
  "Date" with physical ERA5 forcings); averages over all coalitions; robust under collinearity; conditional
  Shapley bounds the off-manifold artefact on both correlated pairs — φ_LAD ≈ 0 holds on the coherent
  Hmax-real subset (+0.022 vs +0.023 off-manifold), and φ_LAI / φ_fCover are stable whether the collinear
  partner is real or baseline (gaps 0.022 / 0.010 °C), so neither the φ_LAD ≈ 0 result nor the LAI/fCover
  co-dominance is a chimeric-coalition artefact. The −0.43 → −0.04 field contrast demonstrates why
  full-control Shapley beats naive partial correlations / GAMM partial effects. Attribution is model-internal,
  spatially validated.
- **D3bis — Why a mechanistic model (MuSICA) rather than a statistical one.** The rationale is physical
  (energy/radiative/turbulent balance, not a fitted slope or equilibrium statistical model) and generalisable
  in principle (transferable parameters), which a site-fitted statistical model is not. Statistical
  microclimate models describe; MuSICA lets us attribute to physical canopy traits and, downstream (Ch. 3),
  spatialise dynamically.
- **D4 — Robustness & scope.** Trait ranking identical for v3.2.0 and v3.2.3; v3.2.0 keeps the spatial
  structure (r 0.93 vs 0.63), v3.2.3 only trims the warm bias. Single-site scope restated. Moran's I: the
  observed buffering slope is strongly spatially autocorrelated (I = 0.30, p < 0.001 kNN=5); the model slope
  reproduces this but dampened (I = 0.13, p = 0.02) — another facet of the glass ceiling; the validation
  residual is not autocorrelated (I = 0.09, p = 0.07 kNN=5; I = 0.015, p = 0.29 inverse-distance) → the
  per-plot inference is not inflated by spatial clustering. cLHS sampling along the structural gradient
  further mitigates.
- **D5 — Perspective (Chapter 2).** Sentinel-2 optical illusion (saturation); the naive all-satellite run
  (S2 LAI + FORMS-H height → MuSICA) is a Chapter-3 baseline.

## 5. Conclusion
Leaf quantity and cover (LAI, fCover) are the co-dominant levers of simulated summer buffering — fCover a
near-constant floor, LAI the between-plot differentiator — while the vertical profile adds little (a strong
nuance to H1/H2). The 3-D structure distributes the between-plot variance well (r ≈ 0.93) but the model
plateaus on the absolute level and the amplifiers → an honest glass ceiling, opening Chapter 2.

## Figures
- **Fig 1** typology — LAD profiles of P1–P4 (`fig_archetypes_profiles`, LAI one-sided)
- **Fig 2** Shapley ranking — LAI ≈ fCover ≫ Hmax > LAD ≈ 0
- **Fig 3** Shapley by cluster — redistribution along the density gradient
- **Fig 4** validation / glass ceiling — REF vs HOBO, r ≈ 0.93
- **Fig 5** forward validation in Shapley order
- **Fig 6** vertical gradient (illustration)

## Annexes
- **A1** FPCA of doubly-normalised LAD shape (→ clustering) + k-means elbow + cluster map (A1a–d)
- **A2** H2 field validation (COM vs buffering; null at full control)
- **A3** H2 controlled in-silico (top-heavy cools more; saturates)
- **A4** conditional Shapley for LAD | Hmax (off-manifold bound, structural co-definition)
- **A4b** conditional Shapley for LAI | fCover & fCover | LAI (off-manifold bound, collinearity r = 0.76; gaps ≤ 0.02 °C)
- **A5** extra metrics (ΔTmin / amplitude / stability)
- **A6** version robustness v3.2.0 vs v3.2.3 (a scatter, b heatmap, c sim-vs-obs)
- **A7** Sentinel-2 / optical illusion = saturation (Chapter-2 bridge)
- **A8** glass ceiling decomposition (a topography R² = 0.51, b slope buff/amp)
- **A9** sampling-height robustness (nair==1 vs fixed 1 m)
- **A10** trait collinearity matrix (LAI/Hmax/fCover/COM/VCI)
