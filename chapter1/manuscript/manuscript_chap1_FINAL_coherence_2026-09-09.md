---
title: "Leaf quantity outweighs vertical arrangement in simulated understorey microclimate buffering: a LiDAR–MuSICA attribution in a temperate oak forest"
bibliography: ../../Bib_these.bib
---

*A preliminary version of this work was presented at the Microclimate Ecology and Biogeography (MEB)
conference (Montpellier, France, June 2026).*

**Keywords:** Microclimate buffering, Understorey temperature, Leaf area density profile, Airborne LiDAR, Multilayer canopy model, Sensitivity analysis, Temperate deciduous forest

# Key message

Airborne laser scanning coupled to a multilayer canopy model shows that the amount of foliage, more than its
vertical arrangement, controls the simulated temperature buffering of a temperate oak understorey. The vertical profile
has an effect in the model, but measured canopies span too narrow a range of shapes for it to change the
buffering between plots.

# Abstract

**Context** Forest understoreys buffer their occupants from macroclimatic extremes. Which dimension of canopy
structure controls that buffering is not settled. Airborne laser scanning (ALS) resolves the vertical
dimension. Canopy variables covary in real stands, so observational models cannot attribute the effect
to one of them.

**Aims** We ask which dimension controls the simulated buffering, leaf quantity or vertical arrangement, and
whether the leaf-area-density profile adds anything once leaf area, height and cover are fixed.

**Methods** We coupled four ALS-derived variables to the multilayer model MuSICA: leaf area index, maximum
height, fractional cover and the vertical leaf-area-density (LAD) profile. We perturbed each around the real
canopy of 400 plots in a temperate oak forest (Blois, France; summer 2021). We selected the plots by
conditioned Latin hypercube sampling.

**Results** Leaf area and cover together dominate the sensitivity of the simulated daytime maximum-temperature offset ΔT~max~,
and maximum height is negligible. In the most open stands every effect is below 0.03 °C, too small to rank, and maximum height is negligible (at most 0.001 °C per metre).
The vertical profile has a second-order effect. Outside the most open stands, the real-versus-uniform contrast
stays at or below the cover step. Its rank against the leaf-area step depends on the response
metric and the summary statistic. A controlled sweep identifies the vertical spread of the profile, not its
peak height, as the stronger axis over most of the leaf-area range. Real profiles vary too little in either
axis to exercise this mechanism, so the realised effect is about 0.05 °C.

**Conclusions** Along the Blois canopy gradient, leaf quantity dominates the between-plot sensitivity of the
simulated buffering, and the vertical LAD profile adds little to it. This result bounds what the vertical detail of ALS can
contribute to mapping forest microclimate.

# 1. Introduction

Beneath a forest canopy the climate experienced by organisms departs from the conditions recorded by
standard weather stations. Daytime maxima are attenuated and the diurnal range is compressed relative to the open [@defrenneGlobalBufferingTemperatures2019; @zellwegerForestMicroclimateDynamics2020]. This decoupling is termed microclimatic buffering. It controls the distribution and persistence of understorey plant
communities [@depauwForestUnderstoreyCommunities2022; @haesenMicroclimateRevealsTrue2023], the success of
tree regeneration [@meeussenInitialOakRegeneration2022], and the capacity of forests to act as thermal
microrefugia under macroclimatic warming [@lenoirClimaticMicrorefugiaAnthropogenic2017; @kemppinenMicroclimateImportantPart2024]. Accounting for vegetation-driven microclimatic variation reduces
the climate velocities species must track and alters their direction, at least in tropical forests
[@soiferMicroclimatesSlowAlter2026]. Forest maximum temperatures average 4.1 °C cooler than the open across
98 sites on five continents, and the offset grows as temperatures become more extreme
[@defrenneGlobalBufferingTemperatures2019]. Dedicated sensor networks and
gridded sub-canopy temperature products now quantify it at continental scale
[@lembrechtsMicroclimaticConditionsAnywhere2020; @haesenForestTempSubcanopyMicroclimate2021].

Canopy structure is the first-order control on buffering. Denser and more closed canopies intercept
more radiation and decouple the understorey further [@zellwegerSeasonalDriversUnderstorey2019; @diaz-calafatBroadleavesConifersEffect2023]. Dense old-growth or well-covered stands buffer thermal
extremes the most [@freySpatialModelsReveal2016; @delombaerdeMaintainingForestCover2022]. Airborne laser
scanning (ALS, airborne LiDAR) resolves the full vertical distribution of vegetation
[@bouvier7_GeneralizingPredictiveModels2015]. ALS-derived structural metrics predict sub-canopy temperature
across landscapes [@juckerCanopyStructureTopography2018; @atkinsEffectsForestStructural2023] and have been
used to map the summer temperature offset [@vandewieleMappingSpatialMicroclimate2023]. Simulations with the
multilayer canopy model MuSICA across stand types have identified a canopy-density threshold in buffering
[@bouwenInteractionsEntreStructure]. That study left open which variable drives the effect, and it imposed a
generic profile shape for lack of measured leaf-area-density profiles. We address both gaps here.

LAI is the established correlate of sub-canopy temperature [@hardwickRelationshipLeafArea2015]. The height
at which the foliage concentrates sets where radiation is intercepted and where the canopy's heat sources
and sinks lie [@parkerThreedimensionalStructureOldgrowth2004]. Two canopies of equal leaf quantity (leaf
area and cover) but different vertical arrangement could therefore buffer the understorey by different
amounts. This mechanism has rarely been tested under controlled conditions. Whether the vertical arrangement
itself controls the buffering, or only covaries with leaf area and cover, is not settled.
Observational models face two difficulties. Canopy variables are collinear in real stands, so the
coefficients of a site-fitted model are unstable under concurvity and the contribution of each variable is
not identifiable [@greiserMonthlyMicroclimateModels2018]. Such models cannot separate the role of leaf
quantity from that of arrangement. They also predict poorly outside their calibration range
[@macleanFinescaleClimateChange2017].

A mechanistic model addresses both difficulties. MuSICA is a one-dimensional multilayer
soil-vegetation-atmosphere model that solves the coupled energy, radiative and turbulent balance
[@ogeeMuSICACO2Water2003]. It can be driven with controlled combinations of canopy variables. Each variable's
contribution is then isolated by construction rather than inferred from correlated observations. Such models
have been used to predict and map forest microclimate [@zellwegerMicroclimateMappingUsing2024; @brusseMechanisticallyMappingNearsurface2024], but rarely to attribute it to individual canopy variables.
Because the model is mechanistic, its parameters are not fitted to one site and should transfer to
another. Other mechanistic microclimate models share this property [@macleanMicroclimcMechanisticModel2021; @kolstelaRevealingFinescaleVariability2024]. We coupled ALS-derived canopy variables to MuSICA to isolate
each variable's contribution to summer understorey buffering in a temperate deciduous oak forest. The
400-plot design samples the canopies the forest contains.

This study asks two questions. (1) Which dimension of canopy structure controls the simulated understorey
buffering: the quantity of foliage, or its vertical arrangement? (2) Does the vertical leaf-area-density
profile add anything to the buffering once leaf area, height and cover are fixed?

We hypothesise (H1) that vertical structure modulates radiative and turbulent transfer, so that the
leaf-area-density (LAD) profile is a controlling input of the simulated microclimate. We hypothesise (H2)
that at fixed leaf area, height and cover, top-heavy profiles attenuate the diurnal cycle more strongly,
because radiation is then intercepted higher and the canopy heat source shifts upward. We tested these
hypotheses in three steps. We derived a structural typology
of the LAD profile. We perturbed each variable around every plot's real canopy with MuSICA and computed the
local sensitivity of the daytime maximum-temperature offset ΔT~max~ in each structural type. For the
profile, we contrasted the real vertical profile with a uniform one. We prescribed a Gaussian profile and
swept its peak height and vertical spread at fixed leaf area, height and cover. We repeated the sweep over
the leaf-area range the forest occupies. That sweep tests the
direction H2 asserts and adds the spread of the foliage as a second profile axis.

The perturbation separates the dimensions, because leaf area and vertical arrangement enter the model by
different pathways even though they covary in real stands. Two
outcomes would count against the hypotheses. H1 fails if the real-versus-uniform contrast does not move the simulated buffering,
whatever the between-plot range of measured profiles. H2 fails if raising the peak height at fixed leaf
area does not cool the understorey further. We report
the size of each mechanism separately from the range real stands span in each variable, because a mechanism
can be large in the model while its variable varies little between stands.

# 2. Material and methods

## 2.1 Study site

Blois is a temperate lowland deciduous forest in the central Loire region of France (about 47.6° N,
1.3° E), almost entirely sessile oak (*Quercus petraea*). Elevation spans 94 to 141 m, so topographic
forcing of the microclimate is at most 0.31 °C by lapse rate (Section 2.5).

### 2.1.1 Airborne LiDAR acquisition

Canopy structure was characterised from a leaf-on ALS survey flown on 17 June
2021 [@grilUsingAirborneLiDAR2023]. The platform, scanner, flight parameters and return density are given in
Appendix D. The provider filtered and classified the cloud (vegetation, ground, other) and height-normalised
it against the ground returns. The LAD profiles and scalar variables used downstream were derived from this
normalised cloud (Section 2.2).

## 2.2 LiDAR canopy structure and variables

We reconstructed vertical
leaf-area-density (LAD) profiles by foliage-profile inversion of the vertical return distributions
[@macarthurFoliageProfileVertical1969; @bouvier7_GeneralizingPredictiveModels2015; @kamoskeLeafAreaDensity2019]. From them we derived the leaf area index (LAI), maximum canopy height
(*H*~max~), fractional cover (fCover) and the vertical complexity index (VCI). LAI is the vertically
integrated one-sided LAD, an effective plant area index. The MacArthur–Horn inversion assumes a spherical
leaf-angle distribution with the conventional extinction coefficient *k* = 0.5 and no foliage-clumping
correction, so the value is a proxy for true leaf area. We retained *k* = 0.5 and corrected every profile for
beam scan angle (Appendix D).

Fractional cover is the fraction of returns above a 2 m height threshold within the 20 m cell. All four
variables and the LAD profile were computed on the 20 m grid directly from the returns in each cell, not
aggregated from a finer product. LAD was binned at 1 m resolution from 1.5 m above ground. ALS
retrieves plant area in heterogeneous canopies at least as accurately as in-situ optical sensors
[@vincent6_MappingPlantArea2017], so we took the ALS retrieval as the structural reference.

VCI is the normalised Shannon entropy of the vertical distribution of returns above ground
[@vanewijkCharacterizingForestSuccession2011]. It summarises the evenness of that distribution rather than the
shape of the LAD profile. We report it only as a descriptor of structural complexity (Fig. 1) and never use
it as a model input. Appendix C gives its formula and Section 4.1 discusses what it tracks. We characterised
profile shape by a functional decomposition instead (Section 2.3).

## 2.3 Forest structural typology

To separate architecture from quantity, we applied a functional principal component analysis (FPCA) to the
LAD profiles, on a B-spline basis whose dimension was selected by generalised cross-validation
[@ramsayFunctionalDataAnalysis1997]. Each profile was first double-normalised, rescaled to z/*H*~max~ in
height and to unit integral in density. FPC1 to FPC3 capture about 92% of the profile-shape variance. They
serve only to define the typology and are not interpreted individually.

We clustered the plots by K-means on six z-standardised ALS metrics: LAI, *H*~max~, fCover and FPC1 to
FPC3. We fixed four archetypes for interpretability and tractability. The within-cluster sum of squares
flattens beyond four groups, and the average silhouette width is 0.27 at four clusters
[@rousseeuwSilhouettesGraphicalAid1987]. The four archetypes are labelled P1 to P4 by increasing canopy
density and were frozen before any simulation. Figure 1 summarises them by their mean LAD profiles. The typology only groups the plots; every perturbation is applied per plot.

Within each archetype we selected 100 plots by conditioned Latin hypercube sampling [cLHS,
@minasnyConditionedLatinHypercube2006], implemented in the *clhs* R package
[@roudierClhsConditionedLatin2012]. The sampler ran 10,000 iterations at a fixed seed. Fractional cover was
floored at 0.5, the floor the model applies to its cover input. cLHS selects among the existing 20 m
LiDAR grid cells and does not generate variable values. Each of the 400 plots perturbed in Section 2.4 is
an observed 20 m cell, with its own measured LAD profile and variable combination. The design
reproduces each archetype's variable marginals and correlations and covers the full range, including the
rare open, low-density stands a random draw would miss (Appendix C). The descriptive statistics of Table 1 and Figure 1 are means. Because the design
over-represents the open extremes, we report every effect size as a median over plots, with a bootstrap
interval (Table 2, Table E1). The
design also probes height independently of density, with the widest height spread in the open archetype P1
(Appendix C).

![**Figure 1** Structural typology of the Blois oak forest. **(a)** Forest map coloured by archetype (P1 open to P4 dense) with the 400 cLHS design plots as dots. **(b)** Per archetype, one representative ALS point cloud on a common 0 to 40 m scale, and the mean LAD profile of its 100 design plots with its interquartile band, cut at the median maximum height. The panel labels the profile PAD; annotated *H*~max~, one-sided LAI and VCI are archetype means. **(c)** Per-archetype half-violins of one-sided LAI, *H*~max~, fractional cover and VCI over the 400 design plots, median in white.](figures/article_v323/Fig1_typology_gril_native.png){width=98%}

The complete workflow is summarised in Figure 2.

![**Figure 2** Analysis pipeline. ALS yields four canopy variables at 20 m (LAI, fractional cover, maximum height and the LAD profile); a functional PCA of the profiles feeds a K-means typology of four archetypes and a cLHS design of 400 plots. MuSICA, forced by the CHS 41 station and ERA5 (Section 2.4), simulates ΔT~max~ and the micro-macro slope per plot; each scalar variable is perturbed in fixed native steps and the profile by a real-versus-uniform contrast.](figures/article_v323/Fig2_methodo.png){width=62%}

## 2.4 MuSICA model and simulation design

MuSICA is a multilayer, multileaf biophysical canopy model. It solves the coupled radiative, energy and
turbulent balance through a vertically discretised canopy [@ogeeMuSICACO2Water2003]. Radiation is
attenuated layer by layer following Beer–Lambert extinction theory [@baldocchiSolarRadiationOak1984; @campbellExtinctionCoefficientsRadiation1986]. Appendix F gives its full configuration. Every simulation
ran the whole of 2021 at an hourly time step, so the analysis window, 1 June to 30 September, is preceded by five months of
within-year spin-up (Appendix F). Each structural variable enters MuSICA by its own pathway. LAI
enters as two-sided leaf area distributed across the layers, and the LAD profile sets how that leaf area is
distributed with height. The LiDAR LAI is one-sided (projected) whereas MuSICA expects the total, two-sided leaf
area, so LAI is doubled before input (a units convention, not a clumping correction). Fractional
cover enters as the canopy clumping factor. It leaves the leaf area in the column unchanged and sets how
efficiently that leaf area intercepts radiation, so a closed canopy of a given LAI shades more than a
clumped one. Maximum height enters twice: it sets the depth over which the LAD profile is distributed
and places the forcing level 2 m above the canopy top.

Canopy structure also alters the meteorology just above the canopy [@bouwenMicroclimateVariationsRoughness2026],
so one fixed forcing applied to every structural scenario ignores that feedback. We therefore ran MuSICA version 3.2.3 (v3.2.3) with its atmospheric boundary-layer
coupling in iterative ("yoyo") mode [@bouwenInteractionsEntreStructure], so that canopy structure feeds
back on the meteorology above. The blending height of that coupling is a prescribed input, not derived
internally (Appendix F). The coupling
represents the canopy-atmosphere feedback that a fixed forcing omits
[@bonanModelingCanopyinducedTurbulence2018], at the cost of one extra forcing series, the boundary-layer height.

The forcing combines the local CHS 41 station (altitude 127 m) with ERA5 reanalysis
[@hersbachERA5GlobalReanalysis2020]. The station is the same open-field reference at 1.5 m used for the
offset (Section 2.5). Air temperature, humidity and precipitation come from the station; shortwave and
longwave radiation, wind and surface pressure come from ERA5. We took temperature from the station because
the reanalysis has a warm bias against it that varies with the time of day over the summer window (Appendix
A). The surface-boundary-layer height the yoyo coupling requires is the hourly MERRA-2 boundary-layer height for the
Blois grid cell [@gelaroModernEraRetrospectiveAnalysis2017], median 550 m. It is archived in the same
forcing file (Appendix F).

The ERA5 wind is the 10 m open-field value. All simulations reported here use it as delivered, without
correction. MuSICA needs the wind just above each canopy, and a neutral logarithmic profile applied at each
plot's canopy height would multiply the 10 m value by about 0.41 to 0.48. We treated that
correction as a sensitivity test (Appendix A). Because the same wind series drives every plot and every
step, perturbing *H*~max~ moves the LAD depth and the forcing level but not the wind speed delivered there.
The height step therefore tests the radiative and geometric pathway of height, not its aerodynamic one,
so it is a lower bound on the full height effect.

Every analysis ran on the 400-plot cLHS design (Section 2.3), around each plot's own measured canopy rather
than a synthetic or archetype-mean one. For each plot we ran a reference at its measured variables. We then
perturbed one variable at a time by a fixed step in its own native unit, with the other three at their
measured values. We also contrasted the real LAD profile against a uniform one. Both are specified in
Section 2.6, and all perturbations used the same v3.2.3
configuration. We perturbed per plot rather than the four mean profiles because MuSICA is nonlinear, so
*f*(mean) ≠ mean(*f*). Aggregating per-plot effects then recovers both the archetype mean and its
within-archetype spread, without the unobserved variable combinations a global baseline would impose.

## 2.5 Microclimate metric extraction

Every metric is read at 1 m, the standard sub-canopy measurement height, by vertical interpolation of the
simulated profile.

**(a) Absolute offsets (ΔT~max~, ΔVPD~max~).** ΔT~max~ = T~max,micro~ − T~max,macro~ is the offset between
the sub-canopy and the open-field reference at the daily maximum. Negative values denote buffering (a
cooler understorey) and positive values amplification. The two maxima are matched in time. For each day we located the hour
of the macroclimatic maximum and read the sub-canopy temperature at that same hour
[@bouwenInteractionsEntreStructure]. ΔVPD~max~ is defined identically for vapour-pressure deficit.

The macroclimatic reference is the local open-field station that also drives the model (CHS 41, air
temperature at 1.5 m; Section 2.4), not the free-air reanalysis aloft. The same station was the open-field
reference in the earlier LiDAR-microclimate mapping of this forest [@grilUsingAirborneLiDAR2023]. The
offset and the driving meteorology share one baseline. We did not elevation-correct it across
plots. Over the 94 to 141 m relief, a lapse-rate correction at 6.5 K km^-1^ would reach at most
0.31 °C. A per-plot offset of this kind cancels in the within-plot perturbation differences that define
every effect here.

**(b) Micro-macro slope.** The offset measures how much the understorey is cooled. The slope measures how
strongly it tracks the macroclimate. We obtained it by regressing the hourly 1 m sub-canopy temperature on
the station temperature at the same clock hour, over all hours of the analysis window. A slope below one
means the understorey damps the macroclimatic range and buffers; a slope above one means it amplifies.

Only the vertical-gradient illustrations depart from this 1 m readout and use the full set of 15 air layers (Fig. B1).

## 2.6 Variable-perturbation sensitivity analysis

From the perturbed simulations (Section 2.4) we computed each variable's effect on the two response
metrics of Section 2.5, ΔT~max~ and the micro-macro slope. Each was computed over the full summer and over
the hottest 10% of days, the 13 days of the window with the highest station daily maximum. For scalar variables
the sensitivity is the metric change a fixed native step produces at the plot's operating point; for the
profile it is the real-versus-uniform contrast. The analysis is local and one-at-a-time
[@saltelliGlobalSensitivityAnalysis2007; @pianosiSensitivityAnalysisEnvironmental2016].

We report the profile contrast raw, as the real-versus-uniform change in ΔT~max~. We did not rescale it by
VCI or per standard deviation, because a given VCI corresponds to infinitely many profiles. Uniform is the
neutral reference, with no preferred direction. The controlled sweep of the profile addresses the height at
which the foliage is concentrated (Hypothesis H2, below). The contrast is an upper bound on the profile's
effect, because it takes the profile to a uniform limit that real canopies reach only in part.

Each scalar variable was moved by a fixed step in its own native unit: leaf area by ±0.5
(one-sided), fractional cover by ±10 points and maximum height by ±1 m. The upward and the downward step are
reported separately, never averaged, so any curvature in the response is visible. Each variable was bounded
to the range the stand can physically occupy: leaf area at a floor of 0.1, fractional cover between 0.5 and
1, and maximum height between 2 and 40 m. Where a step would cross a bound it was truncated. The change
plotted is then the one simulated, without rescaling to the nominal step, so an effect near a bound is
understated, never inflated.

The leaf-area and height steps preserve profile shape by construction (Appendix F). The uniform
reference is a constant density from 1 m to the canopy top with the plot's own leaf area, maximum
height and cover, so it differs from the real profile in shape alone. The profile has no fixed step of this
kind. It is reported as the full contrast, which is not commensurable with the scalar steps. We keep every
effect's sign throughout, and report each for both metrics and both periods.

Table 1 gives each archetype's variable means and within-archetype standard deviations, so each fixed step
can be compared with the size of the variable itself.

**Table 1** Structural variables and fixed steps per archetype. Within-archetype mean and standard deviation (SD) of each scalar variable, in native units (one-sided LAI; fractional cover; maximum height in m). The fifth column gives each fixed step (LAI ±0.5, cover ±10 points, height ±1 m) as a percentage of the archetype mean. Leaf area and cover are stepped by 10 to 28% of their means, height by 3 to 8%. In P4 cover is 0.92 ± 0.03 against a ceiling of 1, so the cover step is truncated to about 8 points.

| Archetype | LAI (mean ± SD) | fCover (mean ± SD) | *H*~max~, m (mean ± SD) | fixed step, % of the archetype mean (LAI · cover · height) |
|:--|:--:|:--:|:--:|:--:|
| P1 (open) | 1.80 ± 0.86 | 0.56 ± 0.08 | 13.3 ± 7.8 | 28% · 18% · 8% |
| P2 | 3.15 ± 0.74 | 0.78 ± 0.09 | 26.0 ± 6.0 | 16% · 13% · 4% |
| P3 | 3.69 ± 0.81 | 0.82 ± 0.07 | 18.8 ± 7.1 | 14% · 12% · 5% |
| P4 (dense) | 5.03 ± 0.71 | 0.92 ± 0.03 | 32.6 ± 3.7 | 10% · 11% · 3% |

One qualification applies to the most collinear pair, leaf area and fractional cover (*r* = 0.94 across the
400 cLHS plots). Moving one while the other is held fixed leaves their observed joint distribution. The model's
split between leaf area and cover is therefore a model-only attribution the field cannot confirm
(Section 4.1). The profile contrast keeps leaf area, height and cover fixed by construction, so it is not exposed to
that collinearity. Vapour-pressure-deficit sensitivities are model-internal results, since no sub-canopy humidity
measurement is available for comparison.

**Controlled test of H2.** The real-versus-uniform contrast measures whether the profile departs from a
uniform one, but not how the foliage is arranged. To test Hypothesis H2 directly we prescribed a Gaussian
LAD profile and swept its two shape parameters, μ and σ. The same (μ, σ)
parameterisation describes vertical canopy profiles in Sentinel-2 sensitivity studies [@fujiwaraEvaluatingSensitivitySentinel22026]. The peak height μ sets the dominant canopy height.
The width σ sets the vertical spread of the foliage. Both are fractions of *H*~max~, and
each profile is renormalised to preserve leaf area. We crossed five values of μ from 0.30 to 0.70 with five of σ from 0.12 to 0.36,
at five one-sided leaf-area levels (LAI 2 to 6), a grid of 125 simulations. Cover was fixed at 0.87 and maximum height at 25 m. MuSICA
re-bins the profile onto ten vegetation layers, so σ below about 0.10 *H*~max~ is not resolved, and the swept σ values
stay above it. The sparse levels do not correspond to observed stands. A one-sided LAI of 2 under a cover of
0.87 has no counterpart in the open P1 (cover 0.56, height 13.3 m; Table 1).
The Gaussian is unimodal, so bimodal real profiles lie outside the sweep.

# 3. Results

## 3.1 Structural typology

The FPCA and K-means typology partitioned the plots into four architectural archetypes, from the open, sparse P1 to the dense, closed P4 (Fig. 1). The types order monotonically along a canopy-density and leaf-area gradient (Table 1). The clustering builds that ordering in (Section 2.3).

The typology did not depend on VCI. A principal component analysis of the four variables recovered the
same four groups along a single dominant density-closure axis with and without VCI (Fig. 3).

![**Figure 3** Principal component analysis of the canopy variables with VCI (LAI, fCover, *H*~max~, VCI) and without VCI (LAI, fCover, *H*~max~), points coloured by archetype. VCI's loading splits between the first two components instead of aligning with the density axis (Appendix C). The archetypes separate equally with and without it, so VCI is collinear with the density variables.](figures/article_v323/Fig3_pca_recover_clusters.png){width=98%}

## 3.2 Sensitivity of ΔT~max~ to each canopy variable, per archetype

Leaf quantity dominated from the intermediate to the closed archetypes and was largest in the closed P4
(Fig. 4, Table 2). In the densest archetype P4, adding 0.5 of leaf area cooled the understorey by 0.205 °C
and the cover step cooled it by 0.265 °C. In the open P1 the same steps gave +0.027 °C, a
warming, and −0.015 °C. The two steps are not commensurable, and the cover step is truncated at the cover ceiling in P4 (Table 1), so we do not rank leaf area against cover. Maximum height was
negligible in every archetype, |effect| ≤ 0.001 °C per metre.

We keep the open P1 in the attribution. Real profiles are most bottom-heavy there, with an LAD centroid at
0.40 of *H*~max~, the lowest of the four archetypes (Section 4.1). P1 is therefore the regime that most
tests the vertical dimension. Every effect was small there: 0.027 °C for the leaf-area step, 0.015 °C for
the cover step, 0.016 °C for the full profile contrast and 0.001 °C per metre or less for height. Including
P1 does not change the ordering found in P4.

The real-versus-uniform profile contrast, at fixed leaf area, cover and height, was a second-order effect
throughout. It moved ΔT~max~ by a median of 0.016 °C in P1, 0.052 °C in P2 and 0.065 °C in P3, cooling in
each. In P4 it moved by +0.104 °C, where the real profile warmed the understorey relative to a uniform one. From
P2 to P4 the contrast stayed at or below the cover step (0.052 against 0.069, 0.065 against 0.068, 0.104
against 0.265 °C). In the open P1 all four effects were below 0.03 °C, and we do not order them. The contrast's rank against the leaf-area step depended on the summary statistic. In P2 and P3 the contrast fell below that step under the mean and reached or exceeded it under the median, so we leave the two unordered. Per plot, the
contrast exceeded the leaf-area step in 47% of the 400 plots and the cover step in 42%.

The full contrast is an upper bound. It takes the canopy to the uniform limit, a range real canopies
span only in part. Because the clusters are defined partly on the profile components FPC1 to FPC3, most
profile-shape variance lies between archetypes. Within a type the profile varies little (Section 4.1). Replacing each plot's real profile
by its archetype's mean shape moved ΔT~max~ by a median of 0.011 °C in P1, 0.027 °C in P2, 0.055 °C in P3
and 0.071 °C in P4. These values are below the cover step in every archetype.

Table 2 gives each effect's median and 95% bootstrap interval over the plots within each archetype. Leaf
area and cover both strengthened towards the dense end: the cover step cooled by 0.015 °C in P1 and by
0.265 °C in P4. Maximum height stayed within 0.001 °C per metre of zero in every archetype and in both
directions, and its widest interval bound was 0.003 °C. The two directions of the leaf-area step were not
mirror images in P4 (0.205 against 0.172 °C), so the effect steepens as the canopy closes. Cover was
asymmetric in the same sense (0.265 against 0.224 °C). The leaf-area and profile intervals are disjoint in
P1 and P4, overlapping in P2 and nearly touching in P3. Leaf area exceeded the profile contrast in P4 (0.205
against 0.104 °C); in P1 every effect was below 0.03 °C.

**Table 2** Median change in ΔT~max~ (°C) per archetype for each fixed step (LAI one-sided), with 95% bootstrap interval over the plots within each archetype. Both directions are expressed as the change per increase of the variable, the downward step being sign-reversed, so two equal rows mean a symmetric response. The profile is the full real − uniform contrast. Steps are truncated at the variable bounds, so an effect near a bound is understated rather than inflated. Truncation also removes plots. Every row rests on all 100 plots of its archetype except cover −10 points, which uses 47 plots in P1 and 97 in P2, because cover there is already at or near the 0.5 floor.

| Effect | P1 | P2 | P3 | P4 |
|:--|:--:|:--:|:--:|:--:|
| LAI +0.5 | +0.027 (+0.024, +0.030) | −0.052 (−0.063, −0.044) | −0.044 (−0.054, −0.034) | −0.205 (−0.223, −0.191) |
| LAI −0.5 | +0.023 (+0.019, +0.028) | −0.033 (−0.040, −0.026) | −0.030 (−0.039, −0.021) | −0.172 (−0.196, −0.158) |
| fCover +10 pts | −0.015 (−0.018, −0.012) | −0.069 (−0.081, −0.057) | −0.068 (−0.077, −0.059) | −0.265 (−0.281, −0.243) |
| fCover −10 pts | −0.026 (−0.030, −0.017) | −0.054 (−0.063, −0.045) | −0.054 (−0.066, −0.050) | −0.224 (−0.264, −0.210) |
| *H*~max~ +1 m | −0.000 (−0.001, +0.000) | −0.000 (−0.000, +0.000) | −0.001 (−0.003, −0.000) | +0.000 (−0.001, +0.001) |
| *H*~max~ −1 m | +0.000 (−0.001, +0.001) | +0.001 (+0.001, +0.001) | −0.001 (−0.003, −0.000) | −0.001 (−0.002, −0.001) |
| profile (full contrast) | −0.016 (−0.033, −0.010) | −0.052 (−0.063, −0.035) | −0.065 (−0.075, −0.053) | +0.104 (+0.072, +0.143) |

On the buffering slope, maximum height stayed negligible everywhere, leaf area and cover remained the two
largest scalar sensitivities, and the leaf-area effect was again largest in the dense P4. The profile's rank against the leaf-area step differed on the slope. Its contrast was +0.002 in P1, −0.002 in P2, +0.008 in P3 and +0.027 in P4
(dimensionless). It matched the leaf-area step in P3 and exceeded it in the dense P4, with the opposite
sign, while staying at or below the cover step in every archetype (Appendix E, Table E1).

Over the hottest 10% of days the ordering was preserved and most effects strengthened. The profile
contrast grew in the intermediate archetypes, its median ΔT~max~ moving from −0.052 to −0.078 °C in P2 and
from −0.065 to −0.121 °C in P3. In P4 the leaf-area step moved from −0.205 to −0.268 °C, a factor of 1.31,
and the cover step from −0.265 to −0.284 °C (Fig. S1).

![**Figure 4** Sensitivity of ΔT~max~ to each canopy variable per archetype (P1 open to P4 dense) over the full summer: leaf area (±0.5, one-sided), fractional cover (±10 points) and maximum height (±1 m), with upward (dark) and downward (pale) steps, and the profile as the full real − uniform contrast at fixed leaf area, cover and height. Violins are the per-plot distributions, with median and interquartile range in white. Panels do not share a y-axis and the steps are not commensurable (Table 1), so panel heights are not a ranking.](figures/Fig4_attribution_chs41.png){width=95%}

Rescaled to one full one-sided unit of leaf area, the median leaf-area effect was +0.054 °C in the open P1,
a warming. It was −0.105 °C in P2, −0.088 °C in P3 and −0.410 °C in the dense P4.

The per-archetype sensitivities did not depend on the clustering. When all 400 plots are pooled without labels, each sensitivity is a smooth function of its own operating point. The archetypes are
bins along that function, and the sensitivities cross the boundaries without a step (Fig. 5). A variance partition agrees. The continuous operating point explained more of the
leaf-area and cover sensitivities than the label did: *R*² = 0.75 against 0.74 for leaf area and
0.74 against 0.72 for cover. For the profile the order reversed: the label explained 0.33 and the operating
point 0.28. Neither predictor explained the height sensitivity (*R*² = 0.00 for both).

![**Figure 5** Sensitivity of ΔT~max~ to each variable against its baseline operating point over the 400 cLHS plots, coloured by archetype with a robust-loess smooth: **(a)** leaf-area step, **(c)** cover step (plots at the 0.5 cover floor excluded) and **(d)** real − uniform profile contrast, all against baseline one-sided LAI, and **(b)** height step against baseline *H*~max~. Each panel gives the *R*² of the sensitivity against the continuous operating point and against the cluster label. The operating point explains more for leaf area and cover, the label more for the profile (0.33 against 0.28).](figures/Fig5_operating_point_chs41.png){width=98%}

## 3.3 Controlled test of the vertical arrangement (H2)

At fixed leaf area, height and cover, two axes of the profile moved ΔT~max~ (Fig. 6). The peak height μ tests
Hypothesis H2 directly. The vertical spread σ is a second axis the hypothesis does not name. The stand
amplified the open-field maximum at one-sided LAI 2 to 3, by up to +0.27 °C. It crossed zero near LAI 4 and
buffered the understorey at LAI 5 to 6, down to −0.78 °C.

Raising the peak cooled the understorey at every level, the direction H2 asserts. The effect grew with leaf
area, from 0.05 °C at LAI 2 to 0.25 °C at LAI 6. At the dense end it was not monotone: a mid-upper peak
(μ = 0.6) buffered most, and the top-heavy shape less.

The vertical spread moved ΔT~max~ more than the peak over most of the range. At fixed peak height it spanned
0.18 to 0.20 °C at LAI 2, 3 and 5, against 0.05 to 0.17 °C for the peak. A more diffuse profile warmed the sparse stands and cooled the dense
ones. Foliage spread lower feeds the
shallow warm layer when the canopy is open, and shades the sub-canopy when it is closed. At LAI 4, the level
nearest the median stand of 3.48, the peak was instead the larger axis (0.12 against 0.08 °C). At one-sided
LAI 6 the spread moved ΔT~max~ by 0.61 °C, three times its value at LAI 5. This is the dense end of the range, above the median stand.

![**Figure 6** Controlled test of H2: simulated ΔT~max~ over the peak height μ (x-axis) and vertical spread σ (y-axis) of a Gaussian LAD profile, both as fractions of *H*~max~, at five one-sided leaf-area levels (LAI 2 to 6; design median 3.48), with cover 0.87 and maximum height 25 m fixed. Red amplifies the open-field maximum and green buffers it; dashed lines mark the middle half of the spread of real profiles (σ 0.21 to 0.28). Raising the peak cools the understorey, and the spread moves ΔT~max~ more than the peak over most of the range; σ below 0.10 *H*~max~ is not resolved by the ten-layer grid.](figures/Fig6_h2_gaussian_chs41.png){width=98%}

# 4. Discussion

## 4.1 Leaf quantity as the dominant control and the realised range of the vertical profile

The quantity of foliage controls the simulated buffering, and its vertical arrangement adds little once
leaf area, cover and height are fixed. In tree-diversity experiments, the buffering gained with species
richness was likewise mediated by higher canopy density [@schnabelTreeDiversityIncreases2025]. In a young
plantation, diverse mixtures that raised canopy cover damped daily vapour-pressure-deficit swings
[@parkForestCompositionDiversity2026]. We
hypothesised that the vertical profile is a controlling input of the microclimate (H1) and that shifting
foliage upward buffers more (H2). Both mechanisms exist in the model, but over the range real stands span
their realised effect is second-order. As a real-versus-uniform contrast, the profile moved ΔT~max~ by a median of 0.02 to 0.07 °C in
the three sparser archetypes and by +0.10 °C in the dense P4. There the real profile warmed relative to a
uniform one. The contrast stayed at or below the cover step in P2 to P4, and its rank against the leaf-area
step depends on the response metric and the summary statistic. The contrast is an upper bound. Measured
profiles are nearly invariant within a structural type, and replacing each real profile by its archetype's
mean shape moved ΔT~max~ by 0.01 to 0.07 °C (Section 3.2). The between-stand variation of real profiles is too small to express the mechanism (H1).

The controlled test supports the direction of H2 and adds a second axis. It swept two axes of the profile at fixed leaf
area, height and cover, over the leaf-area range the forest occupies (Section 3.3, Fig. 6). Raising the
peak cooled the understorey at every level, so H2 holds in direction. Its effect
was the smaller of the two axes over most of the range, 0.05 to 0.17 °C. The larger axis was the vertical
spread, which H2 did not name. At fixed peak height it moved ΔT~max~ by 0.18 to 0.20 °C at LAI 2, 3 and 5. At
LAI 4, nearest the median stand, the peak was the larger axis instead. A more diffuse profile warmed the sparse stands and cooled the dense ones. Both axes steepened with leaf area and reached their widest
range at one-sided LAI 6, the dense end above the median stand of 3.48.

Real stands vary little along either axis. The middle half of measured spreads is σ 0.21 to 0.28, within a swept range of 0.12 to 0.36. Denser stands are both more diffuse and buffer more,
so spread and density cannot be separated. Over that realised spread the axis moves ΔT~max~ by about 0.05 °C, the same order as the contrast in the sparser archetypes (Table 2). The centroid range is as narrow: real profiles
span a middle half of 0.44 to 0.55 of *H*~max~, and the open P1 is lower, at 0.40. The vertical arrangement
is therefore a real control in the model but a second-order one in these stands. This controlled test does not support a controlling role of the vertical dimension beyond leaf quantity, the
question left open in Section 1.

The simulated within-canopy temperature profiles show why the profile contrast stays small at the readout
height (Fig. 7). For each archetype, the real LAD profile and a uniform profile of the same LAI, *H*~max~
and cover give temperature profiles that diverge inside the canopy and converge towards the ground. At 1 m
the gap is at most 0.09 °C. This is a 10:00 to 16:00 average, not the hour-of-maximum ΔT~max~ contrast of
Table 2. The two agree in size, both near 0.1 °C in the dense P4. The profile shape redistributes the heat
source within the canopy, and little of that redistribution reaches the 1 m air temperature.

![**Figure 7** Simulated vertical profile of air temperature per archetype (P1 open to P4 dense), for the real LAD profile (solid) and a uniform profile of the same LAI, *H*~max~ and cover (dashed), June to September, 10:00 to 16:00. The two profiles diverge inside the canopy and converge at the 1 m readout, where the gap is at most 0.09 °C.](figures/Fig7_vertical_Tprofile_chs41.png){width=98%}

The archetypes differ in VCI, and their simulated within-canopy temperature profiles differ by type
(Fig. 7). Both differences follow density. The four archetypes are recovered the same way with and without VCI in the input (Fig. 3). VCI co-varies with leaf area and cover (Appendix C). Because VCI tracks
density, a correlation between VCI and buffering does not show that arrangement matters.

The two ends of the gradient respond in opposite ways. In the open P1, added leaf area
warms the sub-canopy air. With a sparse canopy the ground and the foliage both absorb radiation into a
shallow, poorly mixed layer, and the warm surfaces heat the air faster than the extra foliage shades it. In
the dense P4 the same increment acts through shading and the redistribution of radiation down the profile,
so it cools. The sign flip is consistent with the mechanism, and the ±0.01 control rules out a solver artefact (Appendix E). One caveat qualifies
any comparison across archetypes. The nominal variable values differ between them by construction
(Table 1), so a per-unit sensitivity is read at a different operating point in each type. In the dense P4
cover is 0.92 ± 0.03 against a ceiling of one, so its step is truncated to about 8 points and the effect is read close to saturation. Truncation understates the P4 cover effect. In the open P1 the same fixed steps are large fractions of
small variables, 28% of the mean leaf area against 10% in P4. The open effect is read where the
relative perturbation is largest. Saturation and truncation compress the dense end and do not inflate it.

Leaf area and cover dominate the simulated buffering from the intermediate to the closed archetypes.
Their effect is largest in the closed canopy, where buffering is deepest. The profile remains a second-order
effect of inconsistent sign. In observational data, canopy openness dominates the diurnal temperature range and vertical stratification adds little [@ehbrechtEffectsStructuralHeterogeneity2019]. This ranking rests on the full
contrast, an upper bound. The leaf-quantity control is the variable-level counterpart of the canopy-density
control on the offset reported from observational gradients [@kovacsStandStructuralDrivers2017]. A sensor that measures leaf quantity accurately, without
three-dimensional detail, may therefore resolve the controlling structural information (Section 4.3).

One split remains model-internal: leaf area and fractional cover are too collinear to separate empirically
(*r* = 0.94 across the design). The perturbation supplies that split because the two variables enter MuSICA
by distinct radiative pathways (Section 2.4), but no observation in this study can confirm or refute it. We
therefore read leaf area and cover together as leaf quantity throughout, and never rank one against the
other.

The nearest comparable study agrees in part with this ordering [@bouwenInteractionsEntreStructure]. It imposed a generic beta distribution for the vertical profile, for lack of measured leaf-area-density profiles. The per-plot profiles here supply measured shapes, and the contrast they produce is 0.02 to 0.10 °C over the range real stands span. That study finds maximum height a poor predictor, as we do. Where it reports vertical
complexity correlating with buffering, our test suggests why: the index tracks leaf quantity (Appendix C). The two
studies differ on saturation. That study finds the deciduous offset insensitive to leaf area above about 4
to 6 m² m⁻², and set by cover beyond it. Here the leaf-area effect is largest in the dense P4, at −0.41 °C
per unit of one-sided LAI around a mean LAI of 5.0, and the cover effect is largest there too. The two studies agree on which variables matter, and not on the range over which leaf area acts. Whether the non-saturating
leaf-area effect seen here reflects the leaf-area scale and parameterisation of this canopy remains to
be tested in stands of other density and leaf habit.

## 4.2 Scope of the per-plot perturbation

The variable-perturbation sensitivity analysis resolves the two difficulties of observational microclimate models (Section 1). It measures a local partial response at the observed canopies, and it
attributes the effect to physical variables instead of concurvity-confounded coefficients. Two limits come
with it. Local sensitivities do not sum to the total as an additive budget, but they avoid the
baseline-choice and off-manifold ambiguities such a decomposition incurs under strong collinearity. Each
sensitivity is also conditional. It is read at the operating point of its archetype, so the ranking mixes a
variable's local effect with the archetype's position on the manifold. The density dependence of the
sensitivities is a result and not an artefact of that conditioning (Fig. 5).
The vertical profile is the one input that admits no scalar step, so its effect is a contrast against a
uniform reference shape that preserves the plot's own leaf area, height and cover. The per-plot design also avoids the bias of reading a nonlinear response at the mean profile (Jensen's inequality).

We chose a mechanistic over a statistical model for inference. A model that solves the energy balance
attributes buffering to physical variables. Because it is process-based and not fitted to a site, it
supports extrapolation where a fitted sensitivity or equilibrium model only describes
[@macleanMicroclimcMechanisticModel2021; @macleanFinescaleClimateChange2017]. Such models can be
applied outside the conditions they were built in, as when a land surface model is used to project whether
understorey buffering persists under future warming [@hesProjectingFutureForest2024].

We represent vertical structure by the LAD profile itself and add no further indices. The LAD profile is a physical input to MuSICA: it sets the leaf-area density of each of the ten vegetation layers over which the
radiative and turbulent balance is solved. VCI and any other index of the profile are summaries of that same profile. Adding them as perturbed variables would feed MuSICA nothing new and would reintroduce the
collinearity that rules out regression (Section 1).

At this same site, @grilUsingAirborneLiDAR2023 mapped the buffering slope statistically from three
collinear ALS metrics: maximum height, plant area index and VCI. That model reached *R*² = 0.91 at a 5 m
radius. Its authors note that it cannot say which dimension drives the slope. Coupling the same airborne
structure to MuSICA, we find that leaf area and cover drive the slope. The vertical arrangement,
resolved by the complexity index and the profile, is a second-order term that opposes the buffering at the
dense end (Table E1). The two approaches are complementary: the statistical model maps where the forest
buffers, and this one attributes the buffering to leaf quantity.

## 4.3 Scope and limits

The attribution is a sensitivity analysis internal to MuSICA. It rests on the model's physical
parameterisation (Appendix F) and on the 400-plot cLHS design. It has no independent field validation: no
understorey temperature measurement enters this chapter. The coupled model is evaluated against sub-canopy
measurements in the companion microclimate-forcing study (Corroyez et al., 2026a, in preparation). One model-internal check rules out numerical noise: a leaf-area control of ±0.01 shows the solver deterministic at the scale of the weakest effects (Appendix E). The attribution does not depend on the typology. A principal component analysis of the variables recovers the archetypes
at the 1 m vertical binning with and without VCI (Fig. 3). The ranking does not depend on the archetype
label (Fig. 5).

The four types could be read as developmental stages, but height does not order with density (Table 1). P1 is shortest on average and most variable in height (*H*~max~ 13.3 ± 7.8 m). The intermediate types invert, P2 tall
at 26.0 m and P3 shorter at 18.8 m. P4 is both tallest and most closed. P1 is therefore a heterogeneous open class, and its height spread suggests a mix of sparse plots and opened crowns of older stands. A testable reading orders the types from
open, light-rich canopies, whether young regeneration or thinned mature stands, to tall closed stands.
Linking structural type to developmental stage and management, including the slow recovery of buffering
after harvest [@starckSlowRecoveryMicroclimate2025; @aaltoQuantifyingImpactManagement2023], needs inventory
data this study lacks.

The dominance of leaf quantity has a silvicultural corollary. Thinning and shelterwood cuts lower both leaf area
and cover, the two controlling variables, and so weaken the buffer directly. The vertical arrangement of
foliage is the weaker control in the measured stands but stays relevant, because thinning from below and
thinning from above shift the vertical arrangement in opposite directions. The controlled test moves ΔT~max~ by up to 0.20 °C over the swept profile range at the leaf-area levels these stands occupy, and by more at the dense end
(Section 3.3, Fig. 6). Stands here do not differ enough in arrangement for the effect to appear, though a cut
could create that difference. Because the open-canopy
behaviour of a one-dimensional column is weakly constrained here (Section 3.2) and is tested in the companion
study (Corroyez et al., 2026a, in preparation), this corollary holds as a direction, without a quantitative
prescription. The design compounds the limit: the cLHS sample floors
fractional cover at 0.5, so the attribution has no design points in the regime a heavy cut would create.
The corollary restates the trade-off between opening the canopy for regeneration light and keeping a cool understorey. Buffering lost to harvest returns only as leaf area and cover rebuild, on the slow post-harvest timescale [@starckSlowRecoveryMicroclimate2025].

The study covers a single site, one leaf-on summer season and one species, and the analysis is internal to
MuSICA, so the conclusions are mechanistic and not a landscape prediction. Composition affects the offset:
mixing tree species can enlarge the offset in young plantations [@zhangTreeSpeciesMixing2022]. The leaf-quantity
control is nonetheless found beyond temperate oak. In a historically degraded tropical forest, loss of
canopy density likewise raises understorey temperature [@marshMeasuringModellingMicroclimatic2022].

The response metric is the daytime maximum offset ΔT~max~. This 1 m offset understates the between-type
contrast at the soil surface. Across the four archetype canopies the simulated daytime-mean soil-surface
temperature spans 7.4 °C from the open to the dense type, against 0.5 °C for the daytime-mean air
temperature at 1 m (not shown). Night-time buffering, controlled by longwave trapping, is a distinct
mechanism we do not examine. We claim the mechanistic ordering, not the exact magnitudes.

The result could depend on the forcing and on the leaf-area retrieval. Two corrections leave the variable ranking unchanged: the canopy-height wind correction (Appendix A) and the scan-angle correction (Appendix D).

The 20 m analysis footprint does not drive the result. We re-clipped each plot in circular footprints of 5
to 50 m radius and read every quantity as the shift from the 20 m footprint (Fig. 8). The two response
metrics stay within 0.03 °C for ΔT~max~ and 0.004 for the slope. Of the four structural inputs, only maximum
height moves with the footprint, by about 2 m, because the highest return in the footprint rises with the
area sampled. That variable is thermally negligible (Section 3.2).

![**Figure 8** Shift of the canopy inputs and the two response metrics from their 20 m value, for circular footprints of 5 to 50 m radius, pooled over the 400 cLHS plots by archetype. Top row: the two MuSICA responses with 95% intervals. Lower rows: the four structural inputs, of which only maximum height shifts, by about 2 m.](figures/Fig8_footprint_radius.png){width=98%}

The remaining dependences are properties of the leaf-area input. One assumption links the retrieval to the
model. Leaf area is inverted under a spherical leaf-angle distribution with *k* = 0.5 (Section 2.2), and MuSICA then attenuates radiation under an ellipsoidal one (Appendix F). The two distributions are not identical, and the leaf-area scale is not freely adjustable (Appendix D), so they remain to be reconciled.

A property of the leaf-area input itself is harder to bound. The ALS retrieval returns an
effective plant area (Section 2.2) and does not separate woody from foliar surface, so trunks and branches
are counted within it. We neither quantify that fraction nor bound its consequence. The scale is also
unvalidated. We have no independent ground measurement of leaf area at these plots, so we anchor the
retrieval on the published inversion alone. Because the absolute coupling is conditional on that unmeasured
scale (Appendix D), we report the variable ranking as the result, and not the absolute coupling.

Because leaf quantity controls the simulated buffering, the next step is to drive the same coupling with a
two-dimensional optical proxy of leaf area, for example from Sentinel-2, in place of ALS. The test is
whether it resolves the controlling structural information or saturates before it can
[@aklilutesfayeEvaluationSaturationProperty2021]. At this same site, optical and radar predictors already
recover part of the microclimate signal a LiDAR model captures
[@laslierMappingForestMicroclimates2023]. A uniform 23% reduction in retrieved leaf area follows from the optimised extinction coefficient of
Appendix D. That shift alone changes the simulated coupling, as the companion forcing study measures against
sub-canopy observations (Corroyez et al., 2026a, in preparation). An optical proxy should therefore keep its
leaf-area error below 23% wherever the absolute coupling matters, although the between-plot ranking survives a
shift of that size. Optical indices tend to saturate in dense canopy before reaching such accuracy.

# 5. Conclusion

Over the Blois canopy gradient, leaf quantity dominates the between-plot sensitivity of the simulated
buffering, and its effect is strongest in the closed canopy where buffering is deepest. As a full
real-versus-uniform contrast, the vertical LAD profile moves the offset by a median of 0.02 to 0.07 °C in
the three sparser archetypes, and by +0.10 °C in the densest. In P2 to P4 it stays at or below the cover step.
Replacing each real profile by its archetype mean moves ΔT~max~ by 0.01 to 0.07 °C, so over the design
gradient the profile is a real but second-order control (H1).

Shifting foliage upward buffers the understorey more, the direction H2 asserts. The controlled test finds the
vertical spread of the profile a stronger axis than its peak height over most of the range, up to 0.20 °C
at the leaf-area levels these stands occupy. Real stands vary too little in either axis to reach that effect,
so H2 holds in direction but the profile's realised effect is second-order.

Because leaf quantity controls the simulated buffering, forest management and remote sensing need to track leaf area and cover in these deciduous oak stands. A two-dimensional
measure of leaf area and cover may suffice to map buffering. Section 4.3 gives the accuracy such a
measure must reach. A companion study tests whether Sentinel-2 resolves that information or saturates
before it can (Corroyez et al., 2026b, in revision).

# Appendix A. Forcing: station temperature, and a canopy-height wind sensitivity test

The model is driven by the CHS 41 station for air temperature, humidity and precipitation, and by ERA5
reanalysis for the radiative and aerodynamic variables (Section 2.4).
Pooled over all hours, ERA5 is 0.96 °C too warm (RMSE 2.06 °C, *r* = 0.92; Fig. A1a). The bias depends on
the time of day: it peaks near +3.6 °C in the early morning, around 06:00 to 09:00 UTC, falls through the
day, and turns slightly negative in the evening (Fig. A1b). A forcing with this diurnal bias would
distort the simulated within-canopy diurnal cycle and the daytime offset ΔT~max~ the study targets. Only the
radiative and aerodynamic variables, less exposed to this near-surface bias, are taken from ERA5. Gaps in
the station record are filled from the reanalysis, so the driving series stays continuous and the station
series is used wherever it is available.

![**Figure A1** Bias of ERA5 air temperature against the CHS 41 open-field station, summer 2021. **(a)** Hourly ERA5 versus station temperature (bias +0.96 °C, RMSE 2.06 °C, *r* = 0.92). **(b)** Mean diurnal cycle of the bias, ERA5 minus station, with the interquartile range; it peaks near +3.6 °C in the early morning and turns negative in the evening.](figures/article_v323/FigA1_era5_station_bias.png){width=92%}

The baseline drives every simulation with the ERA5 10 m open-field wind as delivered (Section 2.4). MuSICA
needs the wind just above each canopy, so we tested what a height correction would change. Under a neutral
surface layer, a logarithmic wind profile maps the open-field 10 m wind to the wind at *h* + 2 m above a canopy of
height *h*. The correction factor is *U*(*h*+2m) / *U*(10m) = [ln(2 + *h* − *d*) − ln(*z*~0~)] / [ln(10) −
ln(*z*~0,ERA~)]. The zero-plane displacement is *d* = 0.7*h* and the canopy roughness length
*z*~0~ = 0.1*h*. The ERA5 surface
roughness for the Blois grid cell, taken from the reanalysis forecast-surface-roughness field, is 0.44 m and
nearly constant over the summer window (0.438 to 0.443 m from June to September). Open grassland has a roughness of about 0.01 m and closed forest about 1 m, so the Blois value lies
between them. At that roughness the correction factor
is about 0.41 to 0.48 across the 13 to 33 m archetype canopy heights (Fig. A2). The baseline ERA5 10 m wind
is then roughly 2.1 to 2.4 times the corrected estimate just above the canopy. As a sensitivity test,
we reran the attribution with the corrected wind applied offline per plot at its canopy height, on a subset
of eight plots of the CHS 41 configuration, two per archetype, each rerun with and without the correction. Per unit of leaf area, the leaf-area effect is −0.42 °C in P4
with the correction against −0.41 °C without. It is −0.15 against −0.09 °C in P3, −0.12 against −0.11 °C in P2,
and about +0.09 against +0.05 °C in P1. P4 stays the largest and P1 the only warming. P2 and P3 differ by
at most 0.03 °C in both runs and swap order. The variable ranking is unchanged (Section 4.3).

![**Figure A2** Neutral log-profile correction factor *U*(*h*+2m) / *U*(10m) from the ERA5 10 m open-field wind to the wind just above a canopy of height *h*, with displacement *d* = 0.7*h* and canopy roughness *z*~0~ = 0.1*h*. The four curves are reanalysis roughness lengths *z*~0,ERA~ from 0.01 m (grassland) to 1 m (closed forest), with the Blois value of 0.44 m in bold; at that value the factor is 0.41 to 0.48 over the 13 to 33 m canopy heights of this study. The correction is used only in the sensitivity test of this appendix, not in the baseline simulations.](figures/article_v323/FigA2_wind_profile_correction.png){width=72%}

# Appendix B. The within-canopy gradient, per variable

Microclimate, leaf temperature and gas exchange all vary across forest strata
[@vinodThermalSensitivityForest2023], so the vertical dimension matters ecologically even where it moves our
1 m response metric little. Figure 7 shows the within-canopy temperature profile for the real and the
uniform LAD profile per archetype. Figure B1 extends the same analysis to each variable separately, on
each archetype's median canopy, for temperature, wind, humidity and vapour-pressure deficit. Each
archetype's median canopy was perturbed by the fixed native steps of Section 2.6. Profiles were interpolated
onto a common relative-height grid before averaging, and the axis stops at 1.4 because not all per-plot grids
extend above it. Because the model's layers scale with canopy height, each run was normalised by its own
*H*~max~ (the perturbed run by *H*~max~ + 1 m). Every difference is taken at a matched fraction of canopy
height, so the height step is read at constant relative height. Figure 4 reads it at a fixed 1 m and finds
|effect| ≤ 0.001 °C per metre. These within-canopy profiles are model constructs, with no multi-height measurement to
anchor them.

Among the three scalar steps the ordering is that of the main attribution: leaf area and cover exceed
maximum height in every archetype. On the 10:00 to 16:00 mean, the leaf-area step warms the sub-canopy air
by 0.055 °C in P1 and cools it by 0.227 °C in P4. This mean is a different measure from ΔT~max~. The wind response to the leaf-area step is negative
throughout, as expected when foliage is added. On wind the profile contrast exceeds the scalar steps, but a full real-versus-uniform contrast and a
fixed step are not commensurable. These
representative-canopy panels are consistent with the full-contrast effect of Section 3.2 but do not
establish it independently. Figure 4 remains the per-plot result.

![**Figure B1** Mean within-canopy response of air temperature, wind, relative humidity and vapour-pressure deficit (rows) to each perturbation (coloured lines) per archetype, averaged over June to September and 10:00 to 16:00 on the forcing clock. The dashed line marks no effect, and ΔLAD is the real − uniform profile contrast, not a step.](figures/FigB1_pervariable_gradient_chs41.png){width=98%}

# Appendix C. Variable correlation structure, per archetype and overall

We report the pairwise Pearson
correlations among the canopy variables over all 400 cLHS plots and within each archetype (Fig. C1). The
variables are LAI, fractional cover and maximum height as scalars, and the vertical complexity index VCI as
a scalar descriptor of structural complexity. Over the cLHS design the variables are intercorrelated
(LAI-fCover *r* = 0.94, LAI-VCI *r* = 0.73). A single regression therefore cannot attribute the thermal
effect to one variable. We used a per-archetype variable-perturbation sensitivity analysis instead, around
each plot's own measured canopy (Section 2.4).

This collinearity persists within each archetype. The within-archetype LAI-fCover correlation is
weakest in the open P1, *r* = 0.79. It is 0.97 in P2, 0.93 in P3 and 0.95 in P4.

The per-archetype perturbation does not rely on the variables being decorrelated within an archetype.
It isolates each variable's partial response by moving that variable a small step while the others are
held at their real values. Because LAI and fractional cover enter the canopy scheme as separate inputs with
distinct mechanisms, the model assigns each a distinct sensitivity (Table 2). Their relative size depends on
the choice of steps and is not interpreted as a ranking (Section 3.2). The split between leaf area and cover
is a model construct the collinear field data cannot confirm (Section 4.1). Where they are most collinear,
the single-variable step also departs most from the joint distribution the canopy follows. The
leaf-quantity-versus-arrangement contrast is different in kind: it is a contrast of profile shape at the
plot's own leaf area, height and cover. Leaf area is held by construction, not by statistical
control.

![**Figure C1** Pearson correlation matrices among LAI, fCover, *H*~max~ and VCI (structural-complexity scalar), within each archetype P1 to P4 and over all 400 cLHS points. LAI-fCover reaches *r* = 0.97 in P2, 0.93 in P3 and 0.95 in P4, and 0.79 at its weakest, in the open P1.](figures/article_v323/FigC1_trait_collinearity.png){width=98%}

The design's dispersion reproduces real variability. One-sided LAI disperses comparably across
archetypes and fractional cover saturates towards one in the dense types. Height behaves differently. In
the open P1 height is decoupled from cover (within-archetype Pearson *r* = 0.11). It spans 4 to 35 m there,
almost the whole forest height range. In the dense archetypes height and cover are correlated. In P1 the
design therefore probes height independently of density, and *H*~max~ dispersion is largest there.

VCI is defined as the Shannon entropy of the vertical distribution of returns above ground, normalised to
range from zero to one [@vanewijkCharacterizingForestSuccession2011],

$$\mathrm{VCI} = -\frac{\sum_{i=1}^{HB} p_i \ln p_i}{\ln(HB)},$$

where $p_i$ is the proportion of returns in height bin $i$ and $HB$ the number of 1 m height bins spanning
the canopy. The bins are in absolute metres above ground, so VCI is computed on the metric profile and not
rescaled by maximum height. It is distinct from the height-normalised (z/*H*~max~) LAD profile
used to build the typology (Section 2.3).

To make the redundancy of VCI explicit, we ran a principal component analysis on the four scalar variables
with and without VCI and asked whether the four archetypes are recovered either way (Fig. 3). Including VCI
changes little. Leaf area and fractional cover load together on the first component, the density axis
(loadings −0.55 and −0.56). PC1 explains about 75% of the variance. VCI loads on both components (PC1
−0.47, PC2 −0.55). The case for the redundancy of VCI rests on two other facts. VCI correlates with leaf area and cover (*r* = 0.73 and 0.74; Fig. C1), and the plots
separate into the four archetypes the same way with and without it. VCI adds no independent structural
dimension to the typology: it is a scalar summary of the density gradient the other variables already
describe. The vertical profile therefore enters the model as the full LAD profile (Section 2.4) and not as VCI.

# Appendix D. LiDAR scan-angle correction of leaf area

Correcting the LiDAR leaf area for beam scan angle lowers it by about 5% and leaves the variable ranking
unchanged. The foliage-profile inversion (Section 2.2) counts returns per height layer and converts them to
leaf-area density assuming the beams travel vertically, so the path length through a layer equals the layer
thickness. Off-nadir beams travel a longer path, *l* = *dz* · ⟨sec θ⟩, where θ is the scan angle. The path-corrected
leaf-area density is then the vertical estimate divided by the mean secant of the scan angles over the
returns. Scan geometry biases area-based LiDAR metrics
[@dayalEnhancingForestAttribute2023], and consistent estimates of light attenuation across airborne sensors
require intercalibration [@vincentMultisensorAirborneLidar2023].

The survey was flown by ALTOA with a RIEGL VQ-780i near-infrared scanner (Partenavia P68 at 900 m altitude,
2000 kHz pulse frequency, ±30° nominal scan angle). The mean point density was 72.6 pts m^-2^ (9.7 ground
returns m^-2^) at 4.5 cm vertical precision [@grilUsingAirborneLiDAR2023].

At Blois the return scan angles run from −32° to +33°, just beyond the ±30° nominal field of view. The
excess comes from aircraft roll. The median angle magnitude is 16° and the 95th percentile 29°. This gives
⟨sec θ⟩ ≈ 1.05 and a leaf-area reduction of about 5%, nearly uniform across plots (standard deviation
1.5%; Fig. D1). We applied this correction to all leaf areas reported here. Like the wind sensitivity test
(Appendix A), it leaves the variable ranking unchanged. We verified this by re-running the attribution. A near-uniform 5% rescaling moves every plot's leaf area
by the same fraction, so it shifts the absolute magnitudes without reordering the sensitivities.

![**Figure D1** Scan-angle correction of leaf-area density, characterised on a 53-plot ALS sample at Blois and applied to all 400 design plots. **(a)** Mean leaf-area-density profile from the standard lidR inversion (grey) and after the scan-angle correction (red dashed); the corrected profile lies about 5% below at every height. **(b)** Per-plot ⟨sec θ⟩ factor: all values lie between 1.02 and 1.09.](figures/article_v323/FigD1_scanangle_control.png){width=92%}

We also considered the extinction coefficient itself. An optimised value *k* = 0.65 was obtained for this ALS
product in our companion Sentinel-2/LiDAR retrieval study (Corroyez et al., 2026b, in revision). It
maximises inter-sensor consistency of the leaf area but reduces the leaf area by a further 23%. The leaf area most
consistent between sensors is therefore not guaranteed to be the one that best couples to the microclimate.
We kept the conventional *k* = 0.5 for the coupling rather than import a coefficient optimised for a
different objective. The consequence of
the 23% rescaling for the simulated coupling is assessed in the companion
microclimate-forcing study (Corroyez et al., 2026a, in preparation). There a shift of that size changes the
simulated coupling, so the leaf-area scale is not freely adjustable. A 23% shift in retrieved leaf area is within
the range of clumping and extinction-coefficient uncertainty, so the absolute coupling is not demonstrated
to transfer. At a site whose leaf area departs from the *k* = 0.5 proxy, the absolute effects would not hold. Only the
variable ranking, which the near-uniform corrections leave intact, transfers beyond this
dataset.

# Appendix E. Uncertainty of the variable effects

Each plot contributes both directions of every scalar step, so the data give both the sign and the curvature of each
effect. Table 2 of the main text gives the ΔT~max~ effects with their
intervals. Table E1 gives the same design on the buffering slope and shows which of those patterns depend
on the response metric. Maximum height stays at or below 0.001 per metre. Leaf area and cover keep
both their sign and their strengthening towards the dense end. The three scalar sensitivities therefore
transfer to the slope. The profile contrast is small in the open and intermediate archetypes (+0.002 in P1, −0.002
in P2, +0.008 in P3) and largest in the dense P4 (+0.027). In magnitude it matches the leaf-area step in P3 and
exceeds it in P4 (0.027 against 0.022), with the opposite sign. It stays at or below the cover step in both
(0.012 and 0.030 in magnitude). Its positive sign in P3 and P4 means the real profile weakens buffering relative to the
uniform one.

A separate control rules out numerical noise behind the small effects. Perturbing leaf area on the CHS 41
forcing by a physically negligible ±0.01, a step of 0.02 in total, moved ΔT~max~ by a median of 0.003 °C
over 40 plots. That response is consistent with the linear extrapolation of the ±0.5 step and keeps the
archetype-specific sign, warming by 0.001 °C in the open P1 and cooling by 0.007 °C in the dense P4. In P4 the 0.007 °C matches the
0.008 °C predicted by proportion from the 0.5-step effect (0.205 °C, Section 3.2). MuSICA is therefore
deterministic at this scale, and no numerical floor masks the weak effects. The same
control calibrates the height effect. At 0.001 °C per metre or less, moving maximum height by a whole metre
shifts ΔT~max~ less than moving leaf area by two hundredths of a unit.

**Table E1.** Median change in the micro-macro buffering slope (dimensionless) per archetype for each fixed step, with 95%
bootstrap confidence interval, on the same design and the same plots as Table 2. Unlike Table 2, each step direction is reported with its own raw sign, so the up and down rows of a monotonic response carry opposite signs.

| Effect | P1 | P2 | P3 | P4 |
|:--|:--:|:--:|:--:|:--:|
| LAI +0.5 | +0.002 (+0.001, +0.002) | −0.012 (−0.013, −0.011) | −0.008 (−0.009, −0.007) | −0.022 (−0.024, −0.021) |
| LAI −0.5 | −0.002 (−0.002, −0.002) | +0.011 (+0.010, +0.012) | +0.007 (+0.006, +0.008) | +0.020 (+0.019, +0.021) |
| fCover +10 pts | −0.002 (−0.003, −0.002) | −0.016 (−0.017, −0.014) | −0.012 (−0.013, −0.011) | −0.030 (−0.032, −0.029) |
| fCover −10 pts | +0.003 (+0.003, +0.004) | +0.015 (+0.014, +0.016) | +0.012 (+0.010, +0.013) | +0.029 (+0.027, +0.030) |
| *H*~max~ +1 m | −0.000 (−0.000, −0.000) | −0.000 (−0.000, −0.000) | −0.001 (−0.001, −0.001) | −0.000 (−0.000, +0.000) |
| *H*~max~ −1 m | +0.000 (+0.000, +0.000) | +0.001 (+0.001, +0.001) | +0.001 (+0.001, +0.001) | +0.001 (+0.000, +0.001) |
| profile (full contrast) | +0.002 (+0.002, +0.003) | −0.002 (−0.005, +0.000) | +0.008 (+0.004, +0.015) | +0.027 (+0.023, +0.031) |

# Appendix F. MuSICA model configuration

We simulated stand microclimate with MuSICA v3.2.3, a multilayer soil-vegetation-atmosphere transfer model
[@ogeeMuSICACO2Water2003], using the official binary released under the EUPL. The model ran for a
single broadleaf deciduous species, sessile oak (*Q. petraea*), over the whole of 2021 at an hourly time
step. The forcing file is `MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc`. The
four-month analysis window, 1 June to 30 September, is preceded by five months of spin-up. The
atmosphere is coupled in iterative surface-boundary-layer mode ("yoyo"). The vertical domain comprises 15
air layers and 10 vegetation layers, both fixed by air_resolution_level = 2 and not chosen here. The LiDAR
LAD profile is mapped onto the 10 vegetation layers. Within each time step the model reconciles the canopy-atmosphere fluxes
with the atmosphere up to a blending height. Above that height the flow no longer responds to landscape
heterogeneity [@garrattReviewAtmosphericBoundary1994]. The model takes that blending height as a prescribed input and does not derive it
internally. The surface-boundary-layer height the yoyo requires is an hourly series archived in
the site forcing file, which the model reads directly as the blending height. We took that
series from the MERRA-2 reanalysis for the Blois grid cell [@gelaroModernEraRetrospectiveAnalysis2017]; its
median is 550 m. @bouwenInteractionsEntreStructure instead derived that height from ERA5, following the reanalysis
boundary-layer-height climatology of @seidelClimatologyPlanetaryBoundary2012.

The soil is resolved on an irregular grid to 1.6 m depth with van Genuchten retention. Phenology
follows a parametric deciduous scheme with a single leaf cohort scaled to each plot's leaf area. The leaf
energy balance is solved separately for sunlit and shaded and for wet and dry leaf fractions. The vertical
leaf-area-density profile is the per-plot LiDAR profile, or its uniform counterpart for the profile
contrast; no parametric default is used.

The scalar perturbations of Section 2.6 preserve that profile's shape by construction. A change in leaf
area rescales the whole density profile by a single factor, so every layer keeps its share. A change in
maximum height stretches the profile in relative height, z/*H*~max~, before it is re-interpolated onto the
1 m grid. Neither redistributes foliage vertically.

The main settings are given in Table F1 and the soil-hydraulic and leaf parameters in Table F2; the
complete namelist files are archived with the code.

**Table F1.** MuSICA v3.2.3 configuration, identical for every simulation reported here, for the 400 cLHS design plots of the variable perturbation at Blois (central Loire, France) in 2021. Only the leaf-area-density profile, leaf area index, fractional cover and maximum height vary between plots; they are the per-plot LiDAR inputs. Soil-hydraulic and leaf parameters are in Table F2.

| Component | Setting |
|---|---|
| Model and version | MuSICA v3.2.3 [@ogeeMuSICACO2Water2003], INRAE, EUPL-licensed, official binary |
| Atmosphere coupling | Iterative surface-boundary-layer mode ("yoyo", ABL_flag = 'iter'): near-surface scalars computed up to the surface-boundary-layer height, then recomputed at the reference height |
| Boundary-layer height | Hourly surface-boundary-layer height (h_sbl) from MERRA-2 [@gelaroModernEraRetrospectiveAnalysis2017], archived in the forcing file (median 550 m), the input the yoyo requires |
| Forcing file | `MuSICA_in_CHS41-Blois_2021-station-hsbl-Rmerge.nc`, hourly |
| Forcing sources | Air temperature, relative humidity and precipitation from the CHS 41 station (Blois); incoming shortwave and longwave radiation, 10 m wind speed (as delivered, no height correction; Appendix A), and surface pressure from ERA5; surface-boundary-layer height from MERRA-2, archived in the same file |
| Simulation period | Full year 2021 (analysis window 1 June to 30 September) |
| Time step | Hourly (3600 s) |
| Vertical grid | 15 air layers and 10 vegetation layers, both set by air_resolution_level = 2 and not tuned (the LiDAR leaf-area-density profile is mapped onto the 10 vegetation layers); one leaf-age class |
| Leaf energy balance | Solved separately for sunlit and shaded, and wet and dry leaf fractions |
| Radiation, extinction | Leaf-inclination (ellipsoidal) distribution, leaf inclination index 0.63; LiDAR leaf area entered two-sided (one-sided LAI doubled before input) |
| Canopy clumping | CLUMPING_FACTOR, set per plot to the LiDAR fractional cover (Section 2.4); the namelist default is not used |
| Radiation, soil | Soil albedo 0.15 (visible) and 0.28 (near-infrared); soil emissivity 0.97 |
| Soil, grid | Irregular vertical grid (15 layers) to 1.6 m depth; depth-dependent parameters at 4 nodes (0.04, 0.14, 0.28, 0.8 m) |
| Soil, hydraulics | van Genuchten retention; *θ*~sat~ 0.54 m^3^ m^-3^ at the surface grading to 0.47 m^3^ m^-3^ at depth; residual water content 0.02 m^3^ m^-3^ |
| Soil, initialisation | Initialised on 1 January at 0.40 m^3^ m^-3^ and 283 K at the surface grading to 278 K at depth, then integrated forward with no multi-year spin-up, so the four-month analysis window is preceded by five months of model time |
| Phenology | Parametric deciduous scheme; budburst on day-of-year 115; single leaf cohort with maximum leaf area set to each plot's two-sided leaf area (the vertical distribution is set by the LiDAR profile, not the parametric default) |
| Stand | Single broadleaf deciduous species, sessile oak (*Q. petraea*) |
| LiDAR | Extinction coefficient *k* = 0.5 (Appendix D) |

**Table F2.** Additional soil-hydraulic and leaf parameters as set in the run namelists, identical across plots. Where a parameter is depth-dependent, the four values are the soil nodes of Table F1 in order (0.04, 0.14, 0.28 and 0.8 m depth). Units follow the namelist convention (UNITS FOR *h*~s~ AND *K*~sat~ TO BE SUPPLIED).

| Parameter | Value |
|---|---|
| Soil retention model | van Genuchten (RETENTION_CURVE_MODEL_FLAG = 3), with van Genuchten hydraulic conductivity (HYDRAULIC_COND_MODEL_FLAG = 2) |
| Air-entry parameter *h*~s~ (by depth) | 1.48, 2.17, 2.81, 2.81 |
| Retention shape parameter *n* | 1.05 (all depths) |
| Exponent *m* | 0.5 (all depths); not used under this flag combination, since the retention model fixes its exponent internally to 1 and *m* enters only the Brooks and Corey conductivity branch |
| Saturated hydraulic conductivity *K*~sat~ (by depth) | 0.106, 0.041, 0.013, 0.001 |
| Leaf characteristic dimension | 0.05 m |
| Leaf mass per area (canopy top) | 0.1 kg m^-2^ |
| Leaf inclination index (radiation scheme) | 0.63 |
| Stomatal conductance | slope *g*~1~ = 10, intercept *g*~0~ = 0.001 mol m^-2^ s^-1^, hypostomatous |
| Stomatal water-potential limitation | half-response at leaf xylem potential −1.3 MPa, shape 2.6 |
| *J*~max~ temperature response | optimum 38 °C, curvature *θ* = 0.7 |
| Photosynthetic capacity (*V*~cmax25~, *J*~max25~) | MuSICA v3.2.3 default broadleaf parameterisation |

# Supplementary figures

The main figures are computed over all summer days. We repeated the analyses over the hottest 10% of days (13 days, ranked by the station daily maximum),
because buffering matters most during heat extremes. The result is unchanged: leaf quantity controls the buffering on hot days as
over the full summer. Most effects are larger under heat, and so is their spread. The per-plot dispersion of
the leaf-area effect widens by a factor of about 1.4 to 1.7 relative to the full summer, so we interpret the hot-day
medians as a direction only.

![**Figure S1** Variable attribution as in Figure 4 over the hottest 10% of days. Leaf quantity remains dominant in the dense P4 (leaf area +0.5 giving −0.268 °C, cover +10 points −0.284 °C). The full-contrast profile effect grows in magnitude in the intermediate archetypes P2 and P3 (−0.078 and −0.121 °C, against −0.052 and −0.065 °C over the full summer), and maximum height stays negligible.](figures/FigS1_attribution_hot_chs41.png){width=90%}

# Declarations

## Ethics approval and consent to participate

Not applicable. This study involves no human participants, human data, human tissue, or animals.

## Consent for publication

Not applicable. This manuscript contains no data from any individual person.

## Availability of data and materials

The analysis is implemented as a staged pipeline (`run_chapter1.R`) that regenerates every figure, table and
number reported here from the archived simulation outputs, in a few minutes and without re-running the model.
Regenerating the archived outputs themselves requires MuSICA and the raw airborne laser scanning data, and
takes far longer. The pipeline, the MuSICA run namelists of Appendix F, and the derived variable, design and
simulation products are archived at (REPOSITORY AND DOI TO BE SUPPLIED). MuSICA v3.2.3
[@ogeeMuSICACO2Water2003] is distributed by INRAE under the EUPL (DISTRIBUTION SOURCE TO BE SUPPLIED). The
ERA5 reanalysis [@hersbachERA5GlobalReanalysis2020] is publicly available from the Copernicus Climate Data
Store. The airborne laser scanning point cloud and the CHS 41 station record are available on request, subject
to the agreements under
which they were obtained (PROVENANCE AND CONTACT TO BE CONFIRMED). Software versions used for the retrieval,
typology, design and statistical analysis are to be listed here (R AND PACKAGE VERSIONS TO BE SUPPLIED).

## Competing interests

(TO BE COMPLETED BY THE AUTHORS. If none: "The authors declare that they have no competing interests.")

## Funding

(TO BE COMPLETED BY THE AUTHORS: all sources of funding for the reported research, and whether any funder had
a role in conceptualisation, design, data collection, analysis, the decision to publish, or manuscript
preparation.)

## Authors' contributions

(TO BE COMPLETED BY THE AUTHORS, using author initials, e.g. "NC designed the study, ran the simulations and
wrote the manuscript. All authors read and approved the final manuscript.")

## Acknowledgements

(TO BE COMPLETED BY THE AUTHORS: anyone who contributed but does not meet the authorship criteria. Permission
to acknowledge must be obtained from those named. If none, write "Not applicable".)

# References

::: {#refs}
:::
