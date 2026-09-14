---
title: "General Introduction"
bibliography: ../Bib_these.bib
header-includes:
  - \usepackage{caption}
  # figures carry their own "Figure 1.n." label in the caption text, so
  # pandoc's automatic "Figure n:" prefix would number them twice
  - \captionsetup{labelformat=empty}
---

# 1. General Introduction

## 1.1 The ecological importance of forest microclimates in a changing climate

### 1.1.1 Climate change and the increasing frequency of extreme phenomena

Climate change raises mean temperature, but its impact on living systems operates mostly through
the frequency and severity of extremes. The Earth's energy imbalance, an integrative measure of the
pace of warming, has more than doubled relative to the late twentieth century, and the state of the
climate system is tracked annually against IPCC Sixth Assessment methods
[@forsterIndicatorsGlobalClimate2026]. Hot days, heatwaves, and compound heat and drought episodes
push organisms past physiological thresholds that a slowly warming annual mean never reaches. Under
two degrees of mean warming, a forest acquires more days on which leaf temperature approaches the
threshold for photosynthetic damage, more nights that fail to relieve the accumulated water deficit,
and longer runs of both in succession.

These extremes damage forests directly. Extreme drought reduces the vitality of individual trees,
with effects modulated by tree size and neighborhood [@riederTreeSizeNeighbourhood2026], and
drought-induced dieback and severe wildfire have intensified at stand scale over recent decades
[@cassellWidespreadSevereWildfires2019], with European disturbance regimes projected to intensify
further under continued warming [@grunigClimateChangeWill2026]. Disturbance also removes the canopy
that attenuates the extremes beneath it, so structure lost to one heat and drought episode is
unavailable to buffer the next.

For forest ecology, the relevant questions are how hot it gets on the worst days, and where.
Near-ground temperature responds to stand structure most strongly at the warm extremes, whereas
average and minimum temperatures differ far less between cut and uncut stands
[@potterImpactForestStructure2001], so the daily maximum carries the clearest structural signature.
The second question remains largely unanswered: macroclimatic projections are produced on grids of
one kilometer or coarser, whereas the conditions that determine whether a seedling survives a
heatwave vary over meters.

### 1.1.2 The forest canopy as a physical modifier of the sub-canopy climate

The forest canopy intercepts radiation, slows turbulent exchange, and transpires, so the air beneath
it is partly decoupled from the free atmosphere above [@geigerClimateGround1995;
@defrenneGlobalBufferingTemperatures2019]. The resulting sub-canopy microclimate is quantified
throughout this thesis by a single operational metric, the canopy temperature offset
ΔTmax = Tmicro − Tmacro, the difference between the sub-canopy and above-canopy daily maximum
temperature. The daily maximum is the component most relevant to the heat extremes of Section 1.1.1.
A negative ΔTmax indicates microclimatic buffering, a positive ΔTmax amplification.

Across biomes, closed canopies cool hot daytime maxima and warm cold nights, compressing the
temperature range experienced below [@defrenneGlobalBufferingTemperatures2019;
@zellwegerForestMicroclimateDynamics2020]. The daytime and nighttime effects have distinct causes.
During the day, foliage intercepts shortwave radiation before it heats the forest floor; at night,
the same foliage restricts the escape of longwave radiation to a cold sky. Only the daytime effect
attenuates heat extremes, which motivates ΔTmax rather than the mean or minimum offset as the metric
of interest. Old-growth stands with deep, continuous canopies buffer most strongly
[@freySpatialModelsReveal2016], and maintaining canopy cover is among the few management levers that
preserve this capacity under future climate [@delombaerdeMaintainingForestCover2022].

The sign and strength of the offset depend on canopy structure. Where the canopy is sparse or open,
the offset can reverse and the understory overheats relative to the open air
[@kovacsStandStructuralDrivers2017; @zellwegerSeasonalDriversUnderstorey2019]. Structural
heterogeneity widens the diurnal range in stands whose canopy is uneven rather than uniformly closed
[@ehbrechtEffectsStructuralHeterogeneity2019], and severe disturbance can push a stand from
buffering to amplification [@atkinsEffectsForestStructural2023]. Species composition and stand
density modulate the effect across latitudes, with conifer and broadleaf canopies buffering
differently [@diaz-calafatBroadleavesConifersEffect2023], and management leaves its own signature on
the heterogeneity of buffering within and between stands [@mengeImpactsForestManagement2023]. At
forest edges, sub-canopy conditions converge toward the open within a few tens of meters
[@meeussenMicroclimaticEdgetointeriorGradients2021].

Amplification is documented in field observations, independently of any model, although radiative
overheating of the sensors themselves in open plots can complicate the measurement
[@grilUsingAirborneLiDAR2023]. This observational basis matters when amplification later appears in
a simulation, where it could otherwise be mistaken for an artifact. Buffering also varies through
time: it weakens as soils dry [@greiserHigherSoilMoisture2024] and tracks the phenological cycle of
leaf area through the season [@zellwegerSeasonalDriversUnderstorey2019].

Two qualifications bound the scope of this thesis. First, temperature is not the only buffered
variable. The canopy modifies atmospheric humidity, and with it vapour pressure deficit, at least as
strongly as temperature, and vapour pressure deficit is more proximate to plant water stress but
harder to measure and to downscale [@burtonDownscalingVaporPressure2024]. We focus on temperature
because it is what the logger networks record, and therefore the only quantity against which
simulations can be validated here; where a vapor-pressure-deficit offset is computed, it remains
model-internal and is reported as a robustness check. Second, canopy structure is not the only
driver. Topography redistributes radiation and cold air, soil moisture modulates the evaporative
term, and both interact with structure [@juckerCanopyStructureTopography2018;
@davisMicroclimaticBufferingForests2019]. The microclimate work of this thesis is conducted at a
lowland site of low relief, which limits topographic forcing at the landscape scale and allows
structure to be read with little topographic confounding, without implying that topography is
unimportant in general.

### 1.1.3 Ecological consequences of fine-scale temperature variability

Fine-scale variability is the climate understory organisms experience. Sub-canopy temperature
governs the composition and dynamics of understory plant communities, and buffered forests can slow
the thermophilization of their flora under macroclimate warming
[@defrenneMicroclimateModeratesPlant2013; @depauwForestUnderstoreyCommunities2022]. When species
distributions are re-examined against microclimatic rather than macroclimatic temperatures, the
thermal niches inferred for forest plants shift by several degrees, implying that a large part of
the published climate-change vulnerability of forest understory species has been estimated against
the wrong temperature [@haesenMicroclimateRevealsTrue2023]. Vascular plants and bryophytes differ in
their affinity for buffered conditions, so the microclimate filters communities rather than shifting
them wholesale [@grilAffinityVascularPlants2024].

Regeneration is the second process at stake. Seedling survival and early growth respond to the
temperature and moisture of the first meter above the ground, not to the conditions logged at a
weather station kilometers away [@meeussenInitialOakRegeneration2022;
@vonarxMicroclimateForestsVarying2013]. Juveniles and adults of the same species do not experience
the same climate, and the thermal gap between them has widened over recent decades
[@caronThermalDifferencesJuveniles2021]. Since the juvenile stage occupies the sub-canopy layer that
the overstory buffers, the demographic bottleneck of temperate forests plays out inside the microclimate.
Fauna exploit the same structure: large mammals measurably use canopy-derived thermal shelters
during summer heat [@melinMooseLcesAlces2014].

Under a warming climate, structurally buffered stands can act as thermal microrefugia, where
cool-adapted species persist while the surrounding region warms
[@lenoirClimaticMicrorefugiaAnthropogenic2017; @lenoirUnveilUnseenUsing2022]. Accounting for them
changes the velocity at which species must move to track a shifting climate
[@soiferMicroclimatesSlowAlter2026], and microclimatic heterogeneity within a stand is associated
with structural complexity and with biodiversity [@ehbrechtQuantifyingStandStructural2017]. The
microclimate also modulates ecosystem functioning, mediating part of the effect of plant diversity
on process rates [@beugnonMicroclimateModulationOverlooked2024], and ignoring the sub-canopy layer
biases expectations of biodiversity change [@lembrechtsIncorporatingMicroclimateSpecies2019;
@kemppinenMicroclimateImportantPart2024]. Locating and mapping such refugia requires predicting the
microclimate rather than describing it in a few places [@zellwegerForestMicroclimateDynamics2020].

A microclimate prediction is only as good as the canopy description it is given, and the errors that
matter are structured. A structural input that under-reports foliage in dense canopy makes buffered
stands appear less buffered than they are; one that misses the openness of a gap hides the
amplification that kills a seedling. Both errors compress the simulated range of conditions, and
this range is what determines where a microrefugium lies and whether regeneration survives a
heatwave. A model given such an input can fit the observations on average while misplacing the
extremes. The quality of the structural description is therefore part of the ecological question
itself.

### 1.1.4 Observing microclimate: logger networks and their spatial limitation

The empirical foundation of microclimate ecology is the in-situ temperature logger. Networks of
small sensors, typified by HOBO loggers and coordinated efforts such as SoilTemp, supply large,
standardized archives of sub-canopy conditions [@lembrechtsMicroclimaticConditionsAnywhere2020;
@lembrechtsGlobalMapsSoil2022]. These measurements are accurate, resolve the daily cycle, and
constitute the reference against which any microclimate prediction must be judged
[@macleanMeasurementMicroclimate2021]. Every simulation reported in this thesis is validated against
them.

Logger deployment requires care. A sensor exposed to direct sunlight reads its own radiative load
rather than the air temperature, so shielding, orientation, and mounting height are part of the
measurement, and the bias is worst in the open plots where amplification is expected
[@macleanMeasurementMicroclimate2021; @grilUsingAirborneLiDAR2023]. Network design matters equally:
the placement of loggers determines which part of the structural gradient is sampled, and a network
that oversamples closed canopy underestimates the range of conditions present
[@lembrechtsDesigningCountrywideRegional2021]. Dense local networks can be interpolated into
high-resolution grids where the instrumentation effort is sufficient
[@brunaHighresolutionMicroclimaticGrids2026], but the effort scales with area.

A logger nevertheless measures one point. It cannot map how ΔTmax varies continuously across a
heterogeneous forest landscape, because the sensor network is always far sparser than the structural
variation it samples [@zellwegerAdvancesMicroclimateEcology2019]. Fifty loggers in a forest of
several thousand hectares cover a vanishing fraction of the structural configurations present, and
the configurations they miss are not random. Wall-to-wall knowledge of the microclimate therefore
requires predicting temperature between the points, from variables that can be observed everywhere.

Predicting fine-scale microclimate over a landscape thus requires combining a plot-scale mechanistic
understanding of how canopy structure controls temperature with continuous, remotely sensed
descriptions of that structure. Each sensor imposes a trade-off: airborne LiDAR resolves
three-dimensional structure but only once and locally, whereas satellite optical imagery is repeated
but two-dimensional and prone to saturation. The three research chapters address this tension.

## 1.2 Drivers of forest microclimate and the role of canopy structure

### 1.2.1 Linking three-dimensional vegetation structure to sub-canopy microclimate

Sub-canopy buffering results from four coupled energy-balance processes, each of which depends on
the vertical position of the foliage.

The first is shortwave interception. Solar radiation is attenuated approximately exponentially with
cumulative leaf area, so the energy available to heat the understory decreases as the canopy
thickens. The extinction coefficient governing this attenuation is not a universal constant: it
depends on leaf inclination, solar zenith angle, and foliage clumping
[@baldocchiSolarRadiationOak1984; @campbellExtinctionCoefficientsRadiation1986], so two canopies
with identical leaf area but different leaf angle distributions transmit different amounts of light.
The second is longwave exchange, which partly reverses the first at night: foliage emits downward
and screens the understory from a cold sky, which keeps closed stands warmer than the open after
dark. The buffering of minima and of maxima are therefore physically distinct phenomena.

The third is turbulent transfer. Foliage exerts aerodynamic drag, so wind speed decays with depth
into the canopy and exchange with the atmosphere above weakens [@geigerClimateGround1995]. The
consequence is not uniformly cooling: weak mixing insulates the understory from warm air aloft but
also traps heat generated below, so a sparse canopy that intercepts little radiation while still
suppressing ventilation can amplify rather than buffer. Representing this exchange within and just
above the canopy remains one of the harder problems in land surface modeling
[@bonanModelingCanopyinducedTurbulence2018].

The fourth is latent heat. Energy consumed by evaporation does not warm the air, and canopy
transpiration is only part of this flux: the understory vegetation and the soil surface evaporate as
well, and where soil water remains available this is a substantial daytime energy sink.
Transpiration in temperate oak stands tracks leaf area closely and varies strongly between years
with water availability [@bredaIntraInterannualVariations1996], and understory leaf area is itself a
non-negligible and poorly measured quantity [@georgeMethodComparisonIndirect2021]. Buffering
consequently weakens as the local water balance dries, so the same canopy buffers less in a drought
year than in a wet one [@davisMicroclimaticBufferingForests2019; @greiserHigherSoilMoisture2024].

Two smaller terms close the balance. Heat conducted into and out of the soil stores energy during
the day and releases it at night, so soil and air offsets are correlated but not interchangeable,
and logger networks increasingly record both [@lembrechtsGlobalMapsSoil2022]. The energy
partitioning is also expressed in humidity, since the latent flux that does not warm the air raises
its vapor content, making temperature and humidity buffering two readings of the same partitioning
[@burtonDownscalingVaporPressure2024].

Canopy models represent these evaporative terms least evenly, usually collapsing the understory
layer and the soil surface into a lower boundary condition rather than resolving them with the
detail afforded to the overstory. This asymmetry matters for interpreting simulated sub-canopy
temperature, because energy a model fails to route into evaporation is routed into warming the air
instead.

Each of these processes is height-dependent, so canopy structure is three-dimensional rather than
scalar. Radiation is intercepted where the leaves are, drag is exerted where the foliage is, and the
canopy heat source forms at the height where the foliage is concentrated. Detailed three-dimensional
descriptions of old-growth canopies were built to reason about this coupling between architecture,
radiation balance, and gas exchange [@parkerThreedimensionalStructureOldgrowth2004], and modeling
studies indicate that vertical canopy architecture affects transpiration and leaf thermoregulation
at fixed total leaf area [@banerjeeEffectVerticalCanopy2018]. Two stands with the same total leaf
area, one with foliage packed near the top and one with it spread down the profile, present
different vertical arrangements of sources and sinks to the air beneath. Whether this difference
produces a measurable difference in ΔTmax, and how large it is relative to changing leaf quantity,
is the question of the next section.

### 1.2.2 Which dimension of structure? Leaf quantity and vertical arrangement

Canopy structure is usually summarized by the leaf area index (LAI), the one-sided leaf area per
unit ground area [@chenDefiningLeafArea1992]. LAI is a single number, it enters radiative transfer
directly, and it correlates strongly with observed buffering across gradients from closed forest to
open plantation [@hardwickRelationshipLeafArea2015; @vonarxMicroclimateForestsVarying2013]. The term
itself requires care. Most instruments and remote-sensing retrievals return an effective LAI that
assumes a random foliage distribution and underestimates true leaf area wherever the canopy is
clumped [@chenDefiningLeafArea1992; @chenEvaluationVegetationIndices1996]; optical instruments and
gap-fraction inversions return plant area, which includes woody elements
[@brownNearinfraredDigitalHemispherical2024]; and radiative-transfer schemes differ in whether they
expect a one-sided or a two-sided leaf area. A large part of the disagreement between LAI products
in the literature is definitional before it is physical.

Leaf area is not a constant of a stand. It rises through leaf-out, plateaus through summer, and
declines through senescence, and in deciduous forest the amplitude of this cycle exceeds most of the
between-stand variation observed at any single date [@bredaIntraInterannualVariations1996;
@cotrina-sanchezPhenologyEuropeanForests2026]. A single-date structural description is therefore a
snapshot of a moving quantity, which is why the temporal signature of each sensor (Section 1.4)
matters.

LAI is also an integral, and integrals discard information. The same LAI can be realized by a dense
canopy covering part of the ground or a thinner canopy covering all of it, and by foliage
concentrated in a narrow upper crown layer or distributed through a deep, multi-layered profile.
Complementary descriptors capture what LAI discards: fractional cover (fCover) the horizontal
dimension, maximum height (*H*max) the depth of the column available for attenuation, and the
vertical profile of leaf area density (LAD) the arrangement within that column. Characterizing a
canopy by its foliage profile is an old idea [@macarthurFoliageProfileVertical1969] that became
routinely measurable only with laser scanning [@kamoskeLeafAreaDensity2019;
@arnqvistRobustProcessingAirborne2020]. A parallel literature compresses the three-dimensional point
cloud into scalar indices of structural complexity [@mcelhinnyForestWoodlandStand2005;
@beckschaferEnhancedStructuralComplexity2013; @kaneComparisonsFieldLiDARbased2010], and a further
strand distinguishes horizontal from vertical heterogeneity, showing that the two can act in
opposite directions on ecological responses [@carrascoMetricsLidarDerived3D2019].

![**Figure 1.1.** Three schematic canopies with the same leaf area index per unit ground area. In
(a) and (b) the shaded areas are equal and only the vertical arrangement of the foliage differs; in
(c) the same leaf area is gathered over part of the ground, so the profile is locally denser. A
single number cannot separate the three, although they present different surfaces to radiation and
to the air beneath; axes are schematic and carry no measured
values.](figures/Fig_intro_lai_integral.pdf){width=100%}

Which of these dimensions governs the thermal effect is not settled. Some studies report canopy
density or cover as the dominant control [@kovacsStandStructuralDrivers2017;
@zellwegerSeasonalDriversUnderstorey2019], others find that structural complexity and vertical
heterogeneity carry independent explanatory power for the diurnal temperature range
[@ehbrechtQuantifyingStandStructural2017; @ehbrechtEffectsStructuralHeterogeneity2019], and others
show the two changing together along disturbance and land-use gradients
[@atkinsEffectsForestStructural2023; @juckerCanopyStructureTopography2018]. The same ambiguity
recurs outside the microclimate literature, where canopy structural complexity has been proposed as
a driver of productivity in ways equally hard to separate from leaf quantity
[@hardimanRoleCanopyStructural2011; @faheyDefiningSpectrumIntegrative2019].

These variables covary in real stands: tall stands tend to be dense, dense stands tend to be closed,
and closed stands tend to have deep, multi-layered profiles. A study finding LAI dominant and a
study finding vertical complexity dominant may be reading the same underlying gradient through
different predictors, and additional field data collected along that gradient cannot separate them.
We therefore treat the question as open and state it in the terms this thesis can test: does the
total quantity of foliage govern sub-canopy buffering, or does its vertical arrangement? Answering
it requires observations in which the two are moved independently, which no forest provides and only
a model can supply.

## 1.3 Modeling the forest microclimate

### 1.3.1 Statistical versus mechanistic approaches

Two families of models predict sub-canopy temperature, with different limitations.

Statistical models regress the observed offset on structural and topographic predictors, using
linear models, generalized additive models, or machine learning trained on logger networks. Monthly
microclimate models have been built for managed boreal landscapes
[@greiserMonthlyMicroclimateModels2018], fine-scale warming rates have been mapped from terrain and
canopy predictors [@macleanFinescaleClimateChange2017], generalized additive models have been used
to interpolate microclimate within stands [@burnettUsingGeneralizedAdditive2019], and continental
gridded products such as ForestTemp and ForestClim provide sub-canopy bioclimatic variables for
European forests [@haesenForestTempSubcanopyMicroclimate2021;
@haesenForestClimBioclimaticVariables2023], alongside global maps of soil temperature
[@lembrechtsGlobalMapsSoil2022]. Parsimonious formulations summarizing the micro-macro relationship
by a slope (hereafter the buffering slope) and an equilibrium term have proved robust for comparing
sites [@grilSlopeEquilibriumParsimonious2023]. All of these interpolate well within the range of
conditions they were fitted on.

Their limitation is the collinearity of Section 1.2.2. When predictors are strongly correlated,
regression coefficients become unstable and their standard errors inflate; the additive analogue is
concurvity, in which one smooth term can be partly reproduced by the others, so the decomposition
into individual contributions is not identifiable [@dominiciUseGeneralizedAdditive2002]. Machine
learning is subject to the same limitation: permutation-based variable importance is unstable under
predictor correlation, and correlated predictors share importance in ways that depend on the
algorithm rather than on the data [@nicodemusBehaviourRandomForest2010;
@gregoruttiCorrelationVariableImportance2017]. A model of this kind can rank plots correctly while
attributing the ranking to the wrong variable. It is also site-fitted, so it extrapolates poorly to
structural combinations absent from its training data, and it offers no mechanism to separate the
radiative role of leaf quantity from the turbulent role of leaf arrangement.

Mechanistic models solve the coupled energy, water, and radiative balance that produces the offset,
so their parameters are in principle transferable rather than calibrated to one site
[@macleanMicroclimcMechanisticModel2021]. The family includes dedicated microclimate models
[@macleanMicroclimcMechanisticModel2021], validated against observations in temperate understories
[@brusseMechanisticallyMappingNearsurface2024] and applied to fine-scale variability in boreal
forests [@kolstelaRevealingFinescaleVariability2024], radiative-transfer approaches to microclimate
mapping [@zellwegerMicroclimateMappingUsing2024], and land surface models applied to projections of
future forest microclimate [@hesProjectingFutureForest2024]. Mechanistic and empirical approaches
have been compared directly in degraded tropical forest, with neither dominating on accuracy alone
[@marshMeasuringModellingMicroclimatic2022].

Model accuracy itself requires definition. A model can be judged on how close its simulated offset
is to the observed one, on whether it orders plots from most to least buffered, and on whether it
reproduces the observed spread between plots. A model carrying a constant warm bias may rank plots
perfectly while being wrong everywhere in absolute terms, and a model that compresses the simulated
range may look accurate on average while being unable to distinguish a buffered plot from an
amplifying one. Since the ecological use of a microclimate map is usually to locate the cool places
rather than to predict an exact temperature, ranking and spread must be reported alongside bias.
This thesis reports the three separately throughout.

This thesis requires a mechanistic model for controllability rather than for accuracy. Because the
canopy description enters as an input rather than as a fitted covariate, leaf quantity can be
changed while arrangement is held fixed, and arrangement changed while quantity is held fixed. The
collinearity that prevents attribution in observational data is absent from a designed set of
simulations, so attribution follows by construction.

### 1.3.2 The MuSICA model

MuSICA is a multilayer soil-vegetation-atmosphere transfer model that solves energy, water, and
carbon exchange through a vertically discretized canopy [@ogeeMuSICACO2Water2003]. Radiative
transfer is computed layer by layer, separately for sunlit and shaded foliage; the leaf energy
balance and stomatal conductance are solved at each level; turbulent transfer links the layers to
the air above; and soil water and heat are treated in a coupled multilayer column. The model was
developed and evaluated at European forest flux sites across timescales from hourly to yearly
[@ogeeMuSICACO2Water2003; @ogeePartitioningNetEcosystem2003], and has been coupled to LiDAR-derived
structure to study sub-canopy microclimate in conifer stands [@bouwenInteractionsEntreStructure].

That conifer study established two results and one constraint. A multilayer canopy model can
reproduce the ordering of sub-canopy conditions across stands of contrasting density, and the
meteorological forcing must itself be corrected for canopy structure before scenarios of differing
structure can be compared. The constraint concerns the canopy input: for lack of measured vertical
profiles, the distribution of foliage was described by a parametric beta function, a unimodal shape
that cannot represent the stratified canopies produced by partial harvesting, and whose use was
judged to have contributed to an underestimation of the simulated buffering. Its closing
recommendation was to feed measured LiDAR-derived profiles of plant area density into the model in
place of that function.

Three questions consequently remain open. The measured profile has not been substituted for the
assumed one; the effect of the vertical arrangement has therefore not been isolated, since a profile
held to a single parametric family cannot be contrasted against an alternative arrangement at
matched leaf area and height; and the work was conducted in evergreen conifer stands, where the
seasonal cycle of foliage and the geometry of interception differ from those of a temperate
deciduous canopy. How the structural input might be obtained beyond the LiDAR footprint is a
separate problem, taken up in Section 1.3.3.

Three properties make MuSICA the appropriate tool for the question of Section 1.2.2.

First, MuSICA is multilayer: a vertical profile of leaf area density is a native input, whereas a
big-leaf scheme collapses the canopy into a single effective surface and cannot represent the
vertical dimension. Testing the arrangement of foliage against its quantity requires a model that
can accept both. MuSICA belongs to a small family of canopy schemes with this property, alongside
multilayer biophysical models developed for temperate deciduous stands
[@baldocchiHowEnvironmentCanopy2002], three-dimensional radiation and gas-exchange models
[@sinoquetRATPModelSimulating2001], integrated soil-canopy radiative and energy-balance models
[@vandertolIntegratedModelSoilcanopy2009], and the multilayer canopy parameterizations introduced
into land surface models [@bonanModelingCanopyinducedTurbulence2018;
@bonanModelingStomatalConductance2014]. Empirical microclimate downscalers take canopy cover or a
height metric as a scalar predictor and could not accept a profile.

![**Figure 1.2.** Two canopy representations and what each can accept as input. In (a), a multilayer
scheme discretizes the canopy into layers, each holding sunlit and shaded foliage with its own
energy balance, so the vertical arrangement of leaf area can be altered while its total is held
fixed. In (b), a big-leaf scheme collapses the column into one effective surface that can only be
given a scalar, so the profile is lost and the question of Section 1.2.2 cannot be posed; layer
count and profile shape are illustrative and carry no measured
values.](figures/Fig_intro_musica.pdf){width=100%}

Second, the coupling between structure and temperature is physical rather than empirical. No
transfer function is fitted between the two, so a sensitivity computed by perturbing an input is a
property of the physics rather than an estimate contaminated by the covariance structure of a
training set.

Third, the canopy inputs are explicit and separable. Leaf area, canopy height, cover, and the shape
of the vertical profile enter as distinct quantities, which makes controlled perturbation possible
and allows the measured profile to be contrasted against a leaf-area- and height-preserving uniform
alternative.

One limitation bounds all subsequent analyses. MuSICA is a one-dimensional column, so lateral
advection, edge effects, and sub-pixel canopy gaps lie outside its representation, and a plot whose
surroundings differ from the column above it, for instance through a gap the mean structural input
does not resolve, cannot be reproduced. This sets a ceiling on absolute accuracy that is independent
of how well the structural inputs are measured, so the argument developed here rests on relative
comparisons between plots rather than on absolute temperatures.

### 1.3.3 From plot to landscape

Mechanistic canopy models are plot-scale by construction: a one-dimensional column requires, for
every simulated location, a complete structural description of the canopy above that point together
with the meteorological forcing above it. Running one column is computationally cheap; the obstacle
to running one everywhere is the structural input.

Measuring LAI in the field, whether by hemispherical photography, optical plant-canopy analyzers, or
destructive sampling, is slow, weather-dependent, and subject to its own methodological uncertainty
[@weissReviewMethodsSitu2004; @chenEvaluationHemisphericalPhotography1991;
@holstMeasuringModellingPlant2004]. Different indirect methods disagree even when applied at the
same sites, most strongly in the understory layer [@georgeMethodComparisonIndirect2021], and
automated processing has improved reproducibility without changing the effort involved
[@brownHemiPyPythonModule2023]. Measuring a vertical LAD profile in the field is harder still, and a
campaign of that kind yields tens of plots, occasionally hundreds, whereas a landscape contains
millions of model pixels. Bridging plot-scale mechanism and landscape-scale microclimate therefore
requires a structural dataset acquired remotely, continuously, and at fine spatial resolution.

The working resolution is itself a modeling decision. A one-dimensional column assumes that the
canopy above the simulated point is horizontally uniform over the area it represents, so the pixel
must be small enough for this assumption to hold and large enough for the structural retrieval to be
reliable. Too coarse, and a single column stands for a mixture of gaps and closed canopy whose
average behaves like neither [@garriguesInfluenceLandscapeSpatial2006]; too fine, and retrieval
noise dominates the structural signal. No single resolution serves every purpose: attributing a
thermal effect to canopy structure calls for footprints large enough to define a stable vertical
profile, whereas comparing satellite against airborne retrievals is constrained to the 10-m grid the
satellite imposes. Both sit within the range over which sub-canopy temperature varies in temperate
stands [@zellwegerAdvancesMicroclimateEcology2019; @mengeImpactsForestManagement2023]. The next
section examines which remote-sensing instrument can populate that grid, and at what cost in
structural fidelity or temporal coverage.

## 1.4 Remote sensing of forest structure: bridging the spatial scale

### 1.4.1 Airborne LiDAR

Airborne laser scanning (ALS) is the only operational technique that resolves the interior of a
forest canopy over an area. Laser pulses penetrate gaps in the foliage and return a
three-dimensional point cloud; once the ground surface has been identified and the cloud normalized
to height above ground, the vertical gap fraction can be computed layer by layer and inverted,
through a Beer-Lambert formulation, into a profile of plant area density and an integrated LAI
[@bouvier7_GeneralizingPredictiveModels2015; @richardsonModelingApproachesEstimate2009;
@zhengRetrievalEffectiveLeaf2013]. The approach has been established for two decades and validated
across forest types and sensors [@lefskyLidarBlackwellScienceLtdRemote2002;
@morsdorfEstimationLAIFractional2006; @chenUsingLidarEffective2004;
@kwakEstimationEffectivePlant2010; @peduzziEstimatingLeafArea2012], and the reconstruction of full
LAD profiles rather than integrated totals is now routine [@linRetrievalEffectiveLeaf2016;
@kamoskeLeafAreaDensity2019; @arnqvistRobustProcessingAirborne2020].

The retrieval is not parameter-free. The conversion from gap fraction to leaf area rests on an
extinction coefficient that encodes leaf angle distribution and is commonly fixed at a nominal
value, and on an assumption of random foliage distribution that real canopies violate through
clumping [@huUsingAirborneLaser2018]. Flying height, scan angle, and pulse density all affect the
signal's penetration depth, uncertainty in voxel-scale plant area density is substantial and should
be propagated [@pimont5_EstimatorsConfidenceIntervals2018], and canopy height models derived from
the same cloud carry their own processing choices [@khosravipourGeneratingPitfreeCanopy2014]. A
LiDAR LAI is therefore a model-based estimate resting on stated assumptions, not a direct
measurement, and this matters as soon as such an estimate is set against an optical one: neither
side of the comparison is a reference truth.

The same point cloud also yields a digital terrain model, a canopy height model, fractional cover,
and the height and density metrics on which forest inventory modeling rests
[@khosravipourGeneratingPitfreeCanopy2014; @bouvier7_GeneralizingPredictiveModels2015]. An ALS
survey therefore provides a structural description of a stand rather than a leaf-area measurement
alone, and canopy height in particular is both easier to retrieve robustly than leaf area and
correlated with it, making it a natural fallback wherever a full profile is unavailable. Terrestrial
and portable laser scanning resolve structure in still finer detail, reconstructing leaf area
density at the scale of individual crowns [@hosoiVoxelBased3DModeling2006;
@liEstimatingLeafArea2017], but they cover plots rather than landscapes. Airborne scanning covers an
area while still seeing inside the canopy.

ALS-derived structure has become a standard input to microclimate studies, whether to map buffering
directly [@grilUsingAirborneLiDAR2023; @vandewieleMappingSpatialMicroclimate2023], to constrain it
jointly with topography [@juckerCanopyStructureTopography2018], or to explain the thermal shelters
used by wildlife [@melinMooseLcesAlces2014]. For the question of Section 1.2.2, ALS is the sensor
that delivers the vertical arrangement itself rather than a proxy for it.

This structural detail comes with a restrictive spatio-temporal signature. ALS is acquired by
aircraft campaign, which makes it expensive, spatially bounded by the flight plan, and in practice
single-date: a survey describes the canopy as it was on one summer afternoon, says nothing about
budburst, senescence, or the following year, and does not extend beyond its footprint.

### 1.4.2 Optical satellite imagery

Sentinel-2 presents the opposite trade-off. The two-satellite constellation images every land
surface every five days in thirteen spectral bands, four delivered at 10 m and six at 20 m, among
which three red-edge bands sensitive to canopy chlorophyll and leaf area
[@druschSentinel2ESAsOptical2012; @delegidoEvaluationSentinel2RedEdge2011]. The mission was designed
for operational vegetation monitoring, and leaf area index is one of its standard biophysical
products. The mixed resolution matters for forest work: restricting an inversion to the 10 m bands
buys spatial detail at the price of leaving the red-edge and shortwave infrared information unused,
while including the 20 m bands does the reverse.

LAI retrieval proceeds by inverting a canopy radiative-transfer model against the observed
reflectance, with PROSAIL, the coupling of a leaf optical model with a canopy scattering scheme, as
the standard model [@jacquemoudPROSPECT+SAILModelsReview2009; @verhoefLightScatteringLeaf1984;
@feretPROSPECTDModelingLeaf2017; @feretProsailPROSAILLeaf2024]. The inversion is ill-posed, since
different parameter combinations produce nearly identical spectra, and is regularized by prior
information on the parameter ranges [@combalRetrievalCanopyBiophysical2002]. Two operational
strategies exist: the neural network distributed with the Sentinel toolbox, trained on a simulated
database built from prescribed parameter distributions [@weissS2ToolBoxLevel22020], and hybrid
inversion, which trains a machine-learning regressor on a comparable look-up table
[@verrelstMachineLearningRegression2012; @verrelstExperimentalSentinel2LAI2015;
@verrelstQuantifyingVegetationBiophysical2019]. Both have been applied to forests
[@chrysafisRetrievalLeafArea2020; @fernandesNotJustPretty2024] and validated against field data
[@brownValidationBaselineModified2021]. The result is wall-to-wall, repeated, and free.

Producing a usable time series requires a substantial operational chain. Top-of-atmosphere
reflectance must be atmospherically corrected, clouds and their shadows masked, and observation and
illumination geometry recorded, since the inversion depends on them
[@feretPreprocS2PreprocessingSentinel22024]. The surviving observations are irregularly spaced,
because cloud cover in temperate Europe removes a large share of the nominal five-day revisit, so a
continuous trajectory has to be reconstructed by gap-filling and smoothing
[@atzbergerTimeSeriesMonitoring2011; @liangUsingEnhancedGapFilling2023]. The satellite therefore
delivers a smoothed estimate of a seasonal trajectory rather than a set of measurements, and its
principal advantage over an airborne survey lies in that trajectory: the timing of leaf-out and
senescence is observed per pixel rather than assumed [@cotrina-sanchezPhenologyEuropeanForests2026].

Four constraints qualify the product and define the gap set out in Section 1.5.2.
The signal is a two-dimensional reflectance measured from above, so it carries information about the
upper canopy and only indirectly about what lies beneath. Saturation follows: sensitivity declines
as leaf area accumulates and is exhausted somewhere above five, the exact point depending on the
retrieval, the index, and the canopy [@gaoEvaluatingSaturationEffect2023;
@aklilutesfayeEvaluationSaturationProperty2021], and closed temperate forest canopies sit at or
beyond that limit for much of the growing season. The consequence is a compression rather than an
abrupt ceiling: differences in leaf area between dense stands are progressively under-expressed in
the reflectance before they disappear from it. The prior information regularizing the inversion is
tuned primarily for crops, and transferring it to closed forest canopies biases the retrieval, while
correcting the bias trades against variance [@fernandesEvidenceBiasvarianceTrade2024]. Finally,
because the reflectance-to-LAI relationship is non-linear, a pixel containing a mixture of canopy
heights or densities does not return the average of its parts
[@garriguesInfluenceLandscapeSpatial2006; @maImpactSpatialLAI2008]. Continuous temporal coverage
therefore comes at the cost of structural depth, and the cost is largest where the canopy is densest
and most heterogeneous.

## 1.5 Identifying the scientific gap

### 1.5.1 The attribution gap

The first gap concerns causal attribution. Canopy structural variables covary strongly in real
stands, so correlative microclimate models cannot uniquely assign the buffering effect to any one of
them: coefficients are unstable, smooth terms are concurve, and machine-learning importance measures
redistribute themselves among correlated predictors [@dominiciUseGeneralizedAdditive2002;
@nicodemusBehaviourRandomForest2010; @gregoruttiCorrelationVariableImportance2017]. The question of
Section 1.2.2, whether the quantity of foliage or its vertical arrangement governs sub-canopy
buffering, therefore remains open despite an extensive literature on both sides
[@kovacsStandStructuralDrivers2017; @ehbrechtEffectsStructuralHeterogeneity2019].

The gap is methodological: additional field data collected along the same natural gradient cannot
close it, because the gradient itself confounds the predictors. Closing it requires observations in
which the structural variables are moved independently of one another, and since no forest offers
such observations, they must be generated with a mechanistic model driven by controlled combinations
of canopy inputs. Applying the perturbations around each plot's own realistic canopy keeps the
sensitivities physically meaningful rather than describing configurations that could not exist.

Substituting simulation for observation is defensible only under conditions. First, the model's
parameters are physical and are not fitted to the response being explained, so the sensitivity of
simulated temperature to leaf area is a consequence of radiative and turbulent transfer, not of a
coefficient estimated from the gradient at issue. Second, the perturbations are small and centered
on each plot's measured canopy, so every simulated configuration stays close to one that occurs in
the field and the physics is not extrapolated into untested regimes. Third, the model's ability to
reproduce observed between-plot buffering is established separately against the logger network, so
the attribution is made with a model whose behavior on real stands has been characterized
independently. What the design cannot provide is a direct empirical test of the attribution itself:
the attribution is as good as the physics, and the physics is testable in ways a fitted coefficient
is not.

### 1.5.2 Optical versus structural LAI

The second gap concerns the measurement of leaf area itself. Airborne LiDAR and Sentinel-2 both
deliver a quantity called LAI, but they measure different physical things: one integrates
intercepted returns through the whole canopy column, the other inverts a reflectance signal that
originates mostly near the top and stops responding once the canopy is dense. The two are known to
disagree, and the disagreement is expected to grow with canopy density, but its magnitude and
drivers have not been characterized jointly for temperate deciduous forest.

Three candidate drivers plausibly contribute, and they are rarely assessed together. The attenuation
of the optical signal limits the canopy depth contributing to satellite reflectance, so the
satellite may describe only the upper part of the column that the LiDAR integrates in full. The
parameterization of the radiative-transfer inversion supplies prior ranges tuned for crops rather
than forest canopies, and those priors propagate directly into the retrieved leaf area
[@combalRetrievalCanopyBiophysical2002; @fernandesEvidenceBiasvarianceTrade2024]. Horizontal
heterogeneity of the canopy surface within a satellite pixel, combining variation in canopy height
with variation in the underlying terrain, strains the one-dimensional assumption on which the
inversion rests, and does so non-linearly [@garriguesInfluenceLandscapeSpatial2006].

The question is one of inter-sensor consistency rather than accuracy. No ground-based LAI
measurement coincident with the airborne and satellite acquisitions exists at these sites, so
nothing external can arbitrate between the two retrievals. What can and must be established, before
a satellite LAI is used as a structural input to a microclimate model, is what part of the canopy
that number describes.

Two caveats condition everything downstream. The first concerns domain. Consistency established
between two sensors holds over the structural domain on which it was established, and a depth-based
comparison must restrict itself to canopy deep and closed enough to have an interior: pixels below a
minimum canopy height are excluded, since a search over successive canopy depths cannot converge on
stands shallower than the depths it scans, and a high fractional-cover threshold keeps gaps and
edges out of the sample. Whatever agreement is achieved on that restricted domain is silent about
short, open stands, which are the stands where an optical retrieval behaves differently.

The second caveat is that better agreement in leaf area need not propagate into better simulated
microclimate. A mechanistic canopy model responds non-linearly to leaf area, steeply in open canopy
and with strong saturation once the canopy closes, so an identical reduction in leaf-area error is
worth different amounts depending on where along that curve it applies. A retrieval tuned to
minimize leaf-area discrepancy over one part of the range can be a worse microclimate driver over
another, and can even reverse the ordering of plots if the tuning compresses the gradient the model
relies on. Leaf-area accuracy and microclimate skill are therefore related but distinct objectives,
and which one a correction improves is an empirical question.

![**Figure 1.3.** Schematic of the canopy column read by each sensor. A laser pulse travels to the
ground and back, so gap-fraction inversion integrates the whole profile, whereas reflected sunlight
stops carrying information once enough leaf area lies above, so an optical retrieval responds only
to the leaf area shallower than an effective optical depth. From (a) to (b), leaf area increases
while the optical depth stays fixed, so the satellite reads a smaller share of the column, which is
the structural origin of optical saturation; axes are schematic and carry no measured
values.](figures/Fig_intro_sensor_column.pdf){width=100%}

### 1.5.3 The spatio-temporal trade-off

The third gap follows from the first two. ALS resolves three-dimensional structure but is acquired
once and over a limited footprint; Sentinel-2 is repeated and spatially complete but
two-dimensional and saturating. No sensor currently delivers fine three-dimensional structure and
temporal dynamics together.

Driving a mechanistic model beyond the plots where LiDAR exists therefore forces a compromise
between structural fidelity and coverage. This compromise has usually been made implicitly, by using
whichever product was available, and evaluated by pooling all plots into a single accuracy
statistic. Neither the choice nor the evaluation has been examined against the regime dependence
that Sections 1.4.2 and 1.5.2 predict: if the optical signal fails specifically where the canopy is
dense, a single pooled verdict on which sensor is better may conceal two opposite verdicts in the
two halves of the density gradient. Which compromise is right, and whether it is the same everywhere
in a landscape, has not been established.

The three gaps must be addressed in order: the structural dimension the microclimate responds to
determines what a remote sensor needs to measure, and what the optical sensor sees of the canopy
determines whether it measures that dimension. Any judgment about which product should force a
microclimate model over a landscape depends on both. The order of the chapters follows the order of
the gaps.

![**Figure 1.4.** The spatio-temporal trade-off between structural sensors, each placed by how much
of the three-dimensional canopy it recovers and how often it revisits; spatial coverage is stated in
words rather than encoded in symbol size. The corner a mechanistic microclimate model would need,
full structure repeated through the season, is empty: this emptiness is the third gap. Airborne
LiDAR and Sentinel-2 are the two sensors used in this thesis; GEDI is drawn hollow because it is
used in no chapter and is taken up in the General Discussion, and positions are qualitative with no
axis tick values.](figures/Fig_intro_tradeoff.pdf){width=100%}

## 1.6 Research questions and objectives

**RQ1 (Chapter 1).** *Which dimension of canopy structure governs sub-canopy ΔTmax buffering, leaf
quantity or vertical arrangement, and can a mechanistic model attribute it despite the collinearity
of structural variables?*

The objective is to couple ALS-derived canopy structure (LAI, fCover, *H*max, and the LAD profile)
to MuSICA and to run a controlled trait-perturbation sensitivity analysis, in which each structural
variable is moved by a small, native step around every plot's own canopy and the measured vertical
profile is contrasted against a leaf-area- and height-preserving uniform profile. Because each
perturbation is applied per plot around its own operating point, the sensitivities describe the
local response of a real canopy rather than an extrapolation. The effects are stratified by
structural archetype, so the answer is allowed to differ between open and dense canopy, and the
simulated buffering is validated against an in-situ network of 53 temperature loggers. We
hypothesize that vertical structure modulates radiative and turbulent transfer enough to make the
LAD profile a key parameter for a faithful microclimate simulation, and that, at fixed leaf area and
height, foliage concentrated high in the canopy attenuates the diurnal cycle more strongly than
foliage spread through the column.

**RQ2 (Chapter 2).** *Can Sentinel-2 optical LAI be reconciled with LiDAR-derived structural LAI, or
does optical saturation decouple them?*

The objective is to intercompare the two retrievals across three temperate deciduous forests and to
assess the three candidate drivers of their discrepancy identified in Section 1.5.2: the effective
canopy depth contributing to the satellite reflectance, the parameterization of the PROSAIL
inversion, and the horizontal heterogeneity of the canopy surface within the satellite pixel. The
three are addressed in sequence, each correction applied on top of the previous one, so that their
individual contributions can be read. The target is inter-sensor consistency rather than absolute
accuracy, since no coincident field measurement of LAI is available to arbitrate between the
sensors. The expected output is a quantified statement of how well, and under what canopy
conditions, the two retrievals can be brought into agreement, to be carried forward as a condition
of use.

**RQ3 (Chapter 3).** *Which remotely sensed description of canopy structure should force a
mechanistic microclimate model, in which canopy regime, and can a product be built that remains
deployable beyond the footprints where LiDAR exists?*

The objective is to force MuSICA with a set of competing LAI scenarios in which the leaf-area
magnitude is drawn from LiDAR, from Sentinel-2, or from their combination, while canopy height and
the shape of the vertical profile are held fixed at the LiDAR values, a design that isolates the
effect of the magnitude source. The scenarios are adjudicated against the logger network on four
independent axes: a model-free test in leaf-area space that does not involve MuSICA, the
between-plot ranking, the warm bias, and the recovery of the observed thermal amplitude. Two
properties of the configuration govern how these axes are read: every scenario under-disperses the
between-plot amplitude of the simulated buffering, so the coefficient of determination measures how
faithfully a product ranks plots rather than absolute skill, and the adjudication is stratified by
canopy density at the leaf-area value where the two sensors cross over, so that a regime-dependent
answer can be detected. A specific test follows from the second caveat of Section 1.5.2: whether an
optical retrieval tuned to agree with LiDAR over closed, tall canopy remains the better microclimate
driver once applied across the whole density gradient, including the short, open stands its tuning
excluded. We hypothesize that no single sensor wins everywhere, and that the operationally useful
product is a rule for choosing between them rather than a single best retrieval. Any subsequent
landscape-scale mapping depends on this step.

![**Figure 1.5.** Structure of the thesis. Each of the three gaps identified in Section 1.5 gives
rise to one research question and one research chapter. The vertical arrows on the left indicate
that each gap has to be closed before the next can be posed, which fixes the order of the
chapters.](figures/Fig_intro_roadmap.pdf){width=100%}

## 1.7 Overview of the thesis

**Chapter 1** answers RQ1. ALS-derived structure from a lowland sessile oak forest drives MuSICA at
several hundred plots stratified into four structural archetypes, from open to dense, and each
structural variable is perturbed around every plot's own canopy so that its contribution to
simulated ΔTmax is isolated by construction. The
answer is density-dependent, which reconciles the conflicting literature of Section 1.2.2. Leaf
quantity leads at both ends of the density gradient, and its effect steepens as the canopy closes.
The vertical profile, reported as a complete real-versus-uniform contrast, becomes the largest single
factor in intermediate canopy, but this is an upper bound that real canopies barely span, and in the
densest canopy it is again dwarfed by leaf quantity. Because the profile is a function rather than a
scalar, it is isolated by a categorical contrast against a leaf-area- and height-preserving uniform
profile, against which the realistic profile buffers slightly less. A stand-level validation against
the 53 loggers establishes the model's credibility rather than testing the per-plot sensitivities
directly, since the perturbation design points are not the logger plots, and a model-free check on
the same loggers asks whether adding the measured profile to leaf quantity and cover improves the
fit to the observed buffering. The validation shows that the model reproduces the between-plot
buffering ranking well, while carrying a warm bias on the absolute offset and compressing the
micro-macro slope toward unity, so that both buffering and amplification are under-expressed in
magnitude; the residual traces mainly to sub-pixel canopy gaps that a one-dimensional column cannot
represent. Since leaf quantity is what the microclimate responds to, the next question is whether a
satellite can measure it.

**Chapter 2** answers RQ2. Sentinel-2 and LiDAR LAI are intercompared across three temperate
deciduous forests, a montane beech stand and two lowland oak-dominated ones, and the three candidate
drivers of their discrepancy are addressed in sequence. Satellite LAI agrees best with the LiDAR
estimate when the latter is restricted to the upper 6 to 10 m of the canopy, depending on the site.
This effective canopy depth is an empirical proxy, selected by a multi-criteria compromise between
weakly differentiated criteria: it names the depth at which the two sensors are most consistent, not
a measured optical penetration depth, and at one site the individual criteria disagree enough that
the selected value is a statistical compromise.

Adapting the PROSAIL parameter ranges to forest canopies through a multi-criteria Pareto selection
then improves agreement, though unevenly: the gain is large at one site and marginal at another, and
part of what it recovers is agreement lost at the preceding depth-restriction step. The analysis
separates the parameters that can keep their default values from those that need forest-adapted
priors, and shows that most of the improvement survives when the LiDAR-derived priors are replaced
by LiDAR-free ones. Residual discrepancies remain largest in stands whose canopy surface is
horizontally heterogeneous within the satellite pixel, so vertical thickness and horizontal
heterogeneity have to be accounted for jointly, and a structural overestimation of leaf area in
sparse stands persists that no parameter tuning removes. What improves is inter-sensor consistency,
not demonstrated accuracy against a field truth.

The domain of validity is part of the result. The depth-based comparison is posed only where the
canopy has an interior, so the analysis retains pixels of high fractional cover and of canopy height
above the depth range being probed, and excludes the short, open stands below them. The optimized
parameterization is calibrated on closed, tall canopy and carries no guarantee outside it. Chapter 3
uses these results as a statement of where the optical signal can be trusted, not as a correction to
be applied everywhere.

**Chapter 3** answers RQ3. MuSICA is forced at the logger plots, not over a landscape, with
competing LAI products, and the sensor verdict depends on canopy regime. Where the canopy is dense and the optical signal has saturated, LiDAR carries
the between-plot ranking and Sentinel-2 carries none of it; where the canopy is open, Sentinel-2
leads, though by a narrower margin on the buffering slope than on the daytime offset. The apparent
pooled advantage of one sensor over the other dissolves once the plots are stratified at a
near-median leaf-area split. The two halves of that reversal do not carry equal weight: the
dense-canopy dominance of LiDAR is robust to the model's boundary-layer configuration, whereas the
open-canopy usability of Sentinel-2 is the configuration-dependent link. The mechanism is the one
diagnosed in Chapter 2: most of the leaf area of a dense canopy lies below the effective optical
depth and is invisible from above, which also explains why truncating the LiDAR column to that depth
degrades the simulated buffering rather than improving it.

Consistency in leaf area and skill in microclimate come apart, as anticipated in Section 1.5.2. The
forest-tuned optical retrieval optimized in Chapter 2 is the product most consistent with LiDAR, yet
it is the worst microclimate driver tested, inverting the sign of the between-plot relationship. Two
factors account for this, and neither is a defect of the optimization: it was calibrated on closed,
tall canopy, so applying it across the full density gradient extrapolates it into the open, short
stands its calibration domain excluded; and because a mechanistic canopy model responds non-linearly
to leaf area, a tuning that reduces mean discrepancy can flatten the between-plot contrast the model
needs, which costs ranking even as it buys accuracy. Chapter 2 is therefore used as a map of where
each sensor is valid, not as a corrector to be applied. The practical outcome is a regime-dependent
fusion, switched on canopy height so that no per-pixel LiDAR is needed to decide which product to
use, although the forcing itself still draws canopy height and profile shape from LiDAR in both
regimes. Correcting the optical retrieval where it does not saturate makes matters worse, because
the switch has already placed it where its saturation does not bite. The chapter delivers
inter-sensor consistency and an explicit domain of validity, not the elimination of LiDAR: it maps
the trade-off of Section 1.5.3 rather than resolving it, and each sensor is used where it is valid.

Taken together, the three chapters trace one argument. Leaf quantity is what the microclimate
responds to, with the vertical arrangement a bounded, second-order factor that peaks in intermediate
canopy (Chapter 1); the optical sensor that could supply leaf quantity continuously does not see all
of it, and the part it misses is the part that dense canopies hold (Chapter 2); and the consequence
for mechanistic microclimate modeling is a regime-dependent choice of structural input rather than a
single best sensor (Chapter 3). Each chapter narrows the question the next one asks.

One boundary of scope should be stated first. This thesis does not produce a wall-to-wall
microclimate map; it establishes what such a map would have to be built on: which structural
dimension the temperature responds to, what a satellite can and cannot see of that dimension, and
which product should drive the model in which canopy regime. Mapping remains out of reach for two
reasons established in the chapters: the forcing still requires LiDAR structure
wherever the canopy is dense, so coverage is bounded by the flight footprint, and a one-dimensional
column cannot represent the gaps, edges, and lateral exchanges of a real landscape. This boundary
motivates the perspectives of the General Discussion.

The remaining limits are stated here once. The microclimate work rests on a single lowland oak site
instrumented with 53 loggers, of which the chapters use the summer 2021 window, so the operational
thresholds proposed in Chapter 3, in particular the canopy height at which the fusion switches, are
calibrated in-sample and will require recalibration elsewhere; because the loggers are not spatially
independent, the confidence intervals Chapter 3 derives from them are optimistic. The analysis is
confined to the summer window, when the coupling between leaf area and the sub-canopy thermal regime
is expressed and the logger series are comparable across plots. The seasonal question is addressed
but the answer is bounded: within summer a temporally resolved satellite series adds nothing, since
leaf area is on its phenological plateau and there is no timing left to resolve, and outside summer
the chapters do not arbitrate. The logger archive extends through the following seasons and the 2022
heatwaves, and the General Discussion draws on it for perspectives only. No ground LAI coincident
with the acquisitions is available, so Chapter 2 establishes consistency between sensors rather than
accuracy against a truth. And the mechanistic model is one-dimensional throughout, leaving edges,
lateral advection, and sub-pixel gaps outside its representation. Each of these limits marks where
the next study should begin.

Two avenues are treated in the General Discussion rather than in any chapter. The first is
spaceborne LiDAR. GEDI and the global canopy-height products derived from it promise
three-dimensional structure with satellite coverage [@dubayahGlobalEcosystemDynamics2020;
@langHighresolutionCanopyHeight2022; @schwartzFORMSForestMultiple2023], plant area index products
are being validated from the same instrument [@brownStage1Validation2023], and fusion frameworks
with Sentinel-2 are emerging for wall-to-wall retrieval
[@jiaGEDISentinel2IntegrationFramework2026]. Sampling geometry and geolocation uncertainty
nonetheless keep it short, for now, of what plot-scale mechanistic microclimate modeling requires
[@schleichImprovingGEDIFootprint2023].

The second is dimensional. A one-dimensional column omits edges and lateral advection, yet forest
edges are warmer than interiors [@reekForestEdgesAre2025], and the residuals of Chapter 1 are
structured by the sub-pixel gaps such a column cannot see. The continuation is toward
heterogeneity- and edge-aware formulations that couple a canopy scheme to a three-dimensional
atmospheric solver [REF: Grulois 2026, coupled ARPS-MuSICA thesis;
@vandewalleForEdgeClimV103D2026], in which the regime-dependent forcing developed in Chapter 3 can
be evaluated where lateral fluxes matter. Those formulations relax the assumption this thesis has
had to accept.
