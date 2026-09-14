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

Climate change raises mean temperatures and, with them, the frequency and severity of extremes. The
Earth's energy imbalance, an integrative measure of the pace of heating tracked annually against the
methods of the IPCC's Sixth Assessment, has more than doubled relative to the late twentieth century
[@forsterIndicatorsGlobalClimate2026]. Living systems, however, respond less to the shift in the
annual mean than to the hot days, heatwaves, and compound heat and drought episodes that ride on top
of it, because these short, severe episodes push organisms past physiological thresholds that a
slowly warming mean never reaches. A forest under two degrees of mean warming does not warm by two
degrees everywhere and at all times; it acquires more days on which leaf temperature approaches the
threshold for photosynthetic damage, more nights that fail to relieve the accumulated water deficit,
and longer runs of both in succession.

Extremes also act on forests through damage. Extreme drought leaves a measurable imprint on the
vitality of individual trees, modulated by their size and neighborhood
[@riederTreeSizeNeighbourhood2026]; at stand scale, drought-induced dieback and severe wildfire have
intensified over recent decades [@cassellWidespreadSevereWildfires2019], and European disturbance
regimes are projected to intensify further under continued warming
[@grunigClimateChangeWill2026]. Disturbance removes canopy, and canopy attenuates the extremes
beneath it, so structure lost to one heat and drought episode is unavailable to buffer the next.

For forest ecology, the operative question is how hot it gets on the worst days, and where.
Near-ground temperature responds to stand structure most strongly in the warm extremes, whereas
average and minimum temperatures differ far less between cut and uncut stands
[@potterImpactForestStructure2001], so the daily maximum carries the clearest thermal signature of
canopy structure. The "where" is answered less often: macroclimatic projections are produced on
grids of one kilometer or coarser, whereas the conditions that determine whether a seedling survives
a heatwave vary over meters.

### 1.1.2 The forest canopy as a physical modifier of the sub-canopy climate

The canopy intercepts radiation, slows turbulent exchange, and transpires, so the air beneath it is
partly decoupled from the free atmosphere above [@geigerClimateGround1995;
@defrenneGlobalBufferingTemperatures2019]. The result is a distinct sub-canopy microclimate. We
quantify it throughout this thesis with a single metric, the canopy temperature offset
ΔTmax = Tmicro − Tmacro, the difference between the sub-canopy and the above-canopy daily maximum
temperature. The daily maximum is the component most relevant to the heat extremes of
Section 1.1.1. A negative ΔTmax indicates buffering, the understory being cooler than the open air;
a positive ΔTmax indicates amplification.

Across biomes, closed canopies buffer: they cool hot daytime maxima and warm cold nights,
compressing the temperature range experienced below [@defrenneGlobalBufferingTemperatures2019;
@zellwegerForestMicroclimateDynamics2020]. The daytime and nighttime causes differ. During the day,
foliage intercepts shortwave radiation before it heats the forest floor, so the understory stays
cooler than the open. At night, the same foliage restricts the escape of longwave radiation to a
cold sky, so the understory stays warmer. Both effects narrow the diurnal range, but only the
daytime one attenuates heat extremes, hence the focus on ΔTmax rather than on the mean or minimum
offset. Old-growth stands, with deep and continuous canopies, buffer most strongly
[@freySpatialModelsReveal2016], and maintaining canopy cover is among the few management levers
shown to preserve this capacity under future climate [@delombaerdeMaintainingForestCover2022].

The sign and strength of the offset depend on canopy structure. Where the canopy is sparse or open,
the offset can reverse and the understory overheats relative to the open air
[@kovacsStandStructuralDrivers2017; @zellwegerSeasonalDriversUnderstorey2019]. Structural
heterogeneity acts in the same direction, widening the diurnal range in stands whose canopy is
uneven rather than uniformly closed [@ehbrechtEffectsStructuralHeterogeneity2019], and severe
disturbance can push a stand from buffering to amplification
[@atkinsEffectsForestStructural2023]. Species composition and stand density modulate the effect
across latitudes, with conifer and broadleaf canopies buffering differently
[@diaz-calafatBroadleavesConifersEffect2023], and management leaves its own signature on the
heterogeneity of buffering within and between stands [@mengeImpactsForestManagement2023]. Forest
edges add a further gradient, over which sub-canopy conditions converge toward the open within a few
tens of meters [@meeussenMicroclimaticEdgetointeriorGradients2021].

Amplification is documented in field measurements independently of any model, although radiative
overheating of the sensors themselves in open plots can complicate the measurement
[@grilUsingAirborneLiDAR2023]. This observational grounding matters: without it, amplification could
be mistaken for an artifact wherever it appears in a simulation. Buffering also varies through time,
weakening as soils dry [@greiserHigherSoilMoisture2024] and tracking the phenological cycle of leaf
area through the season [@zellwegerSeasonalDriversUnderstorey2019].

Two qualifications bound what this thesis claims. First, temperature is not the only buffered
variable. The canopy modifies atmospheric humidity, and with it vapour pressure deficit, at least as
strongly as temperature; vapour pressure deficit is more proximate to plant water stress but harder
to measure and to downscale [@burtonDownscalingVaporPressure2024]. We lead on temperature because it
is what the logger networks record, and therefore the only quantity against which simulations can be
validated here. Where a vapor-pressure-deficit offset is computed, it remains model-internal for
want of a matching observation and is reported as a robustness check. Second, canopy structure is
not the only driver. Topography redistributes radiation and cold air, soil moisture modulates the
evaporative term, and both interact with structure [@juckerCanopyStructureTopography2018;
@davisMicroclimaticBufferingForests2019]. The microclimate work of this thesis is conducted at a
lowland site of low relief, which limits topographic forcing at the landscape scale. Structure can
be read there with little topographic confounding, but topography remains important in general.

### 1.1.3 Ecological consequences of fine-scale temperature variability

Fine-scale variability is the climate understory organisms experience. Sub-canopy temperature
governs the composition and dynamics of understory plant communities, and buffered forests can slow
the thermophilization of their flora under macroclimate warming
[@defrenneMicroclimateModeratesPlant2013; @depauwForestUnderstoreyCommunities2022]. When species
distributions are re-examined against microclimatic rather than macroclimatic temperatures, the
thermal niches inferred for forest plants shift by several degrees, so part of the published
climate-change vulnerability of forest understory species has been estimated against the wrong
temperature [@haesenMicroclimateRevealsTrue2023]. Vascular plants and bryophytes differ in their
affinity for buffered conditions, so the microclimate filters communities rather than shifting them
wholesale [@grilAffinityVascularPlants2024].

Regeneration is the process most directly tied to canopy structure. Seedling survival and early
growth respond to the temperature and moisture of the first meter above the ground, not to
conditions logged at a weather station kilometers away [@meeussenInitialOakRegeneration2022;
@vonarxMicroclimateForestsVarying2013]. Juveniles and adults of the same species do not experience
the same climate, and the thermal gap between them has widened over recent decades
[@caronThermalDifferencesJuveniles2021]. The juvenile stage occupies the sub-canopy layer that the
overstory buffers, so the demographic bottleneck of temperate forests sits inside the microclimate.
Fauna also exploit the structure: large mammals use canopy-derived thermal shelters during summer
heat [@melinMooseLcesAlces2014].

Under a warming climate, structurally buffered stands can act as thermal microrefugia, pockets
where cool-adapted species persist while the surrounding region warms
[@lenoirClimaticMicrorefugiaAnthropogenic2017; @lenoirUnveilUnseenUsing2022]. Accounting for them
changes the velocity at which species must move to track a shifting climate
[@soiferMicroclimatesSlowAlter2026], and microclimatic heterogeneity within a stand is associated
with structural complexity and with biodiversity [@ehbrechtQuantifyingStandStructural2017]. The
microclimate also modulates ecosystem functioning, mediating part of the effect of plant diversity
on process rates [@beugnonMicroclimateModulationOverlooked2024]. Ignoring the sub-canopy layer
therefore biases expectations of biodiversity change
[@lembrechtsIncorporatingMicroclimateSpecies2019; @kemppinenMicroclimateImportantPart2024]. Locating
and mapping such refugia depends on the ability to predict the microclimate, not only to describe it
in a few places [@zellwegerForestMicroclimateDynamics2020].

A microclimate prediction is only as good as the canopy description it is given, and the errors
that matter are structured. An input that under-reports foliage in dense canopy makes buffered
stands look less buffered than they are; one that misses the openness of a gap hides the
amplification that kills a seedling. Both errors compress the range of conditions a forest
contains, and that range decides where a microrefugium lies and whether regeneration survives a
heatwave. A model given such an input can fit the observations on average while misplacing the
extremes, so the quality of the structural description conditions the ecological conclusion.

### 1.1.4 Logger networks and their spatial limitation

The empirical foundation of microclimate ecology is the in-situ temperature logger. Networks of
small sensors, typified by HOBO loggers and coordinated efforts such as SoilTemp, now supply large,
standardized archives of sub-canopy conditions [@lembrechtsMicroclimaticConditionsAnywhere2020;
@lembrechtsGlobalMapsSoil2022]. These measurements are accurate, resolve the daily cycle, and are
the reference against which any microclimate prediction must be judged
[@macleanMeasurementMicroclimate2021]. Every simulation reported in this thesis is answerable to
them.

Deployment requires care. A sensor exposed to direct sunlight reads its own radiative load rather
than the air temperature, so shielding, orientation, and mounting height are part of the
measurement, and the bias is worst in the open plots where amplification is expected
[@macleanMeasurementMicroclimate2021; @grilUsingAirborneLiDAR2023]. Network design carries the same
weight: the placement of loggers determines which part of the structural gradient is sampled, and a
network that oversamples closed canopy underestimates the range of conditions present
[@lembrechtsDesigningCountrywideRegional2021]. Dense local networks can be interpolated into
high-resolution grids where the instrumentation effort is sufficient
[@brunaHighresolutionMicroclimaticGrids2026], but the effort scales with area.

The core limitation is that a logger measures one point. It cannot map how ΔTmax varies
continuously across a heterogeneous forest landscape, because the sensor network is always far
sparser than the structural variation it samples [@zellwegerAdvancesMicroclimateEcology2019]. Fifty
loggers in a forest of several thousand hectares cover a vanishing fraction of the structural
configurations present, and the configurations they miss are not random. Wall-to-wall knowledge of
the microclimate therefore requires a way to predict temperature between the points, from variables
that can be observed everywhere.

Predicting fine-scale microclimate over a landscape means combining a plot-scale mechanistic
understanding of how canopy structure controls temperature with continuous, remotely sensed
descriptions of that structure. Each sensor imposes a trade-off: airborne LiDAR resolves
three-dimensional structure but only once and locally, whereas satellite optical imagery is
repeated but two-dimensional and prone to saturation. The three research chapters are three
responses to this tension.

## 1.2 Drivers of forest microclimate and the role of canopy structure

### 1.2.1 Energy-balance coupling between canopy structure and sub-canopy climate

Buffering is an energy-balance phenomenon. What reaches the forest floor, and what the sub-canopy
air does with it, is set by four coupled processes, each depending on where the foliage sits.

The first is shortwave interception. Solar radiation is attenuated approximately exponentially with
cumulative leaf area as it passes through the canopy, so the energy available to heat the understory
falls as the canopy thickens. The rate of attenuation is governed by an extinction coefficient that
is not a universal constant: it depends on leaf inclination, solar zenith angle, and foliage
clumping [@baldocchiSolarRadiationOak1984; @campbellExtinctionCoefficientsRadiation1986]. Two
canopies with identical leaf area but different leaf angle distributions therefore transmit
different amounts of light. The second is the longwave exchange that partly reverses the first at
night: foliage emits downward and screens the understory from a cold sky, so closed stands are
warmer than the open after dark, and the buffering of minima and of maxima are physically distinct
phenomena.

The third is turbulent transfer. Foliage exerts aerodynamic drag, so wind speed decays with depth
into the canopy and exchange with the atmosphere above weakens, leaving the sub-canopy air partly
insulated from the free atmosphere [@geigerClimateGround1995]. The consequence is not uniformly
cooling: weak mixing insulates the understory from warm air aloft but also traps heat generated
below, which is one reason a sparse canopy that intercepts little radiation while still suppressing
ventilation can amplify rather than buffer. Representing this exchange within and just above the
canopy remains one of the harder problems in land surface modeling
[@bonanModelingCanopyinducedTurbulence2018].

The fourth term is latent heat. Every unit of energy consumed by evaporation is a unit that does
not warm the air. Canopy transpiration is only part of it: the understory vegetation and the soil
surface evaporate as well, and where soil water remains available that flux is a substantial sink
for daytime energy. Transpiration in temperate oak stands tracks leaf area closely and varies
strongly between years with water availability [@bredaIntraInterannualVariations1996], and
understory leaf area is itself a non-negligible and poorly measured quantity
[@georgeMethodComparisonIndirect2021]. Buffering consequently weakens as the local water balance
dries, so the same canopy buffers less in a drought year than in a wet one
[@davisMicroclimaticBufferingForests2019; @greiserHigherSoilMoisture2024].

Two smaller terms close the balance. Heat conducted into and out of the soil stores energy during
the day and releases it at night, so soil and air offsets are correlated but not interchangeable,
and logger networks increasingly record both [@lembrechtsGlobalMapsSoil2022]. The energy
partitioning is also expressed in humidity: the latent flux that fails to warm the air raises its
vapor content, so temperature buffering and humidity buffering are two readings of one partitioning
[@burtonDownscalingVaporPressure2024].

Canopy models represent the evaporative terms least evenly, usually collapsing the understory layer
and the soil surface into a lower boundary condition rather than resolving them with the detail
afforded to the overstory. This matters whenever a simulated sub-canopy temperature is interpreted,
because energy a model fails to route into evaporation warms the air instead.

Each of these processes is height-dependent, which is why structure is three-dimensional rather
than scalar. Radiation is intercepted where the leaves are, drag is exerted where the foliage is,
and transpiration occurs where the stomata are, so the canopy heat source sits at whatever height
the foliage is concentrated. Detailed three-dimensional descriptions of old-growth canopies were
built to reason about this coupling between architecture, radiation balance, and gas exchange
[@parkerThreedimensionalStructureOldgrowth2004], and modeling studies indicate that vertical canopy
architecture affects transpiration and leaf thermoregulation at fixed total leaf area
[@banerjeeEffectVerticalCanopy2018]. Two stands with the same total leaf area, one with foliage
packed near the top and one with it spread down the profile, do not present the same vertical
arrangement of sources and sinks to the air beneath. Whether that difference translates into a
measurable difference in ΔTmax, and how large it is relative to having more or fewer leaves, is the
question the next section poses.

### 1.2.2 Leaf quantity and vertical arrangement

Canopy structure is usually summarized by the leaf area index (LAI), the one-sided leaf area per
unit ground area [@chenDefiningLeafArea1992]. LAI is attractive because it is a single number,
enters radiative transfer directly, and correlates strongly with observed buffering across gradients
from closed forest to open plantation [@hardwickRelationshipLeafArea2015;
@vonarxMicroclimateForestsVarying2013]. The term itself requires care. Most instruments and most
remote-sensing retrievals return an effective LAI that assumes a random distribution of foliage and
therefore underestimates true leaf area wherever the canopy is clumped [@chenDefiningLeafArea1992;
@chenEvaluationVegetationIndices1996]. Optical instruments and gap-fraction inversions return plant
area, which includes woody elements, rather than leaf area alone
[@brownNearinfraredDigitalHemispherical2024]; radiative-transfer schemes differ in whether they
expect a one-sided or a two-sided leaf area; and comparisons between sensors are only meaningful
once these conventions are aligned. Much of the disagreement between LAI products in the literature
is definitional.

Leaf area is not a constant of a stand. It follows a phenological cycle, rising through leaf-out,
plateauing through summer, and declining through senescence, and in deciduous forest the amplitude
of that cycle exceeds most of the between-stand variation observed at any single date
[@bredaIntraInterannualVariations1996; @cotrina-sanchezPhenologyEuropeanForests2026]. A single-date
structural description is a snapshot of a moving quantity, which is one reason the temporal
signature of each sensor (Section 1.4) matters.

LAI is also an integral, and integrals discard information. The same LAI can be realized by a dense
canopy covering part of the ground or a thinner canopy covering all of it, by foliage concentrated
in a narrow upper crown layer or distributed through a deep, multi-layered profile. Several
complementary descriptors address this. Fractional cover (fCover) captures the horizontal
dimension, maximum height (*H*max) the depth of the column available for attenuation, and the
vertical profile of leaf area density (LAD) the arrangement within that column. The idea of
characterizing a canopy by its foliage profile is old [@macarthurFoliageProfileVertical1969], but it
became routinely measurable only with laser scanning [@kamoskeLeafAreaDensity2019;
@arnqvistRobustProcessingAirborne2020]. A parallel literature compresses the three-dimensional point
cloud into indices of structural complexity [@mcelhinnyForestWoodlandStand2005;
@beckschaferEnhancedStructuralComplexity2013; @kaneComparisonsFieldLiDARbased2010], and a further
strand distinguishes the horizontal from the vertical component of heterogeneity, showing that the
two can act in opposite directions on ecological responses [@carrascoMetricsLidarDerived3D2019].

![**Figure 1.1.** Schematic of three canopies with the same leaf area index per unit ground area. In
(a) and (b) the shaded areas are equal and only the vertical arrangement of the foliage differs; in
(c) the same leaf area is gathered over part of the ground, so the profile is locally denser. No
axis carries measured values.](figures/Fig_intro_lai_integral.pdf){width=100%}

The literature has not settled which of these dimensions governs the thermal effect. Some studies
report canopy density or cover as the dominant control [@kovacsStandStructuralDrivers2017;
@zellwegerSeasonalDriversUnderstorey2019]. Others find that structural complexity and vertical
heterogeneity carry independent explanatory power for the diurnal temperature range
[@ehbrechtQuantifyingStandStructural2017; @ehbrechtEffectsStructuralHeterogeneity2019]. Others show
the two changing together along disturbance and land-use gradients
[@atkinsEffectsForestStructural2023; @juckerCanopyStructureTopography2018]. The same ambiguity
recurs outside the microclimate literature, where canopy structural complexity has been proposed as
a driver of productivity in ways equally hard to separate from leaf quantity
[@hardimanRoleCanopyStructural2011; @faheyDefiningSpectrumIntegrative2019].

The disagreement has a structural cause: in real stands these variables covary. Tall stands tend to
be dense, dense stands tend to be closed, and closed stands tend to have deep multi-layered
profiles. A study that finds LAI dominant and a study that finds vertical complexity dominant may be
reading the same gradient through different predictors, and additional field data collected along
that gradient cannot separate them. We therefore treat the question as open and state it in the
terms this thesis can test: does the total quantity of foliage govern sub-canopy buffering, or does
its vertical arrangement? Answering it requires observations in which the two can be moved
independently, which no forest provides and only a model can supply.

## 1.3 Modeling the forest microclimate

### 1.3.1 Statistical versus mechanistic approaches

Two families of models predict sub-canopy temperature. Statistical, or correlative, models regress
the observed offset on structural and topographic predictors, using linear models, generalised
additive models, or machine learning trained on logger networks. Monthly microclimate models have
been built for managed boreal landscapes [@greiserMonthlyMicroclimateModels2018], fine-scale
warming rates have been mapped from terrain and canopy predictors
[@macleanFinescaleClimateChange2017], generalized additive models have been used to interpolate
microclimate within stands [@burnettUsingGeneralizedAdditive2019], and continental gridded products
such as ForestTemp and ForestClim now provide sub-canopy bioclimatic variables for European forests
[@haesenForestTempSubcanopyMicroclimate2021; @haesenForestClimBioclimaticVariables2023], alongside
global maps of soil temperature [@lembrechtsGlobalMapsSoil2022]. Parsimonious formulations that
summarize the micro-macro relationship by a slope (the buffering slope, as we call it throughout)
and an equilibrium term have proved robust for comparing sites
[@grilSlopeEquilibriumParsimonious2023]. Within the range of conditions they were fitted on, all of
these interpolate well.

Their limitation is the question of Section 1.2.2. When predictors are strongly collinear, a fitted
coefficient is not an effect. Regression coefficients become unstable and their standard errors
inflate; the additive analogue is concurvity, in which one smooth term can be partly reproduced by
the others, so the decomposition into individual contributions is not identifiable
[@dominiciUseGeneralizedAdditive2002]. Machine learning relocates the problem: permutation-based
variable importance is unstable under predictor correlation, and correlated predictors share
importance in ways that depend on the algorithm rather than on the data
[@nicodemusBehaviourRandomForest2010; @gregoruttiCorrelationVariableImportance2017]. A model of this
kind can rank plots correctly while attributing the ranking to the wrong variable. It is also
site-fitted and descriptive, so it extrapolates poorly to structural combinations absent from its
training data, and it offers no mechanism to separate the radiative role of leaf quantity from the
turbulent role of leaf arrangement.

Mechanistic models solve the coupled energy, water, and radiative balance that produces the offset
rather than fitting it, so their parameters are in principle transferable
[@macleanMicroclimcMechanisticModel2021]. The family includes dedicated microclimate models
[@macleanMicroclimcMechanisticModel2021], validated against empirical observations in temperate
understories [@brusseMechanisticallyMappingNearsurface2024] and applied to fine-scale variability in
boreal forests [@kolstelaRevealingFinescaleVariability2024]; radiative-transfer approaches to
microclimate mapping [@zellwegerMicroclimateMappingUsing2024]; and land surface models applied to
projections of future forest microclimate [@hesProjectingFutureForest2024]. Mechanistic and
empirical approaches have been compared directly in degraded tropical forest, with neither
dominating on accuracy alone [@marshMeasuringModellingMicroclimatic2022].

What "accurate" should mean for a microclimate model governs how any of them can be compared. A
model can be judged on how close its simulated offset is to the observed one, on whether it orders
plots from most to least buffered, or on whether it reproduces the observed spread between plots. A
model carrying a constant warm bias may rank plots perfectly while being wrong everywhere in
absolute terms; a model that compresses the simulated range may look accurate on average while
failing to distinguish a buffered plot from an amplifying one. Since the ecological use of a
microclimate map is usually to find the cool places rather than to predict an exact temperature,
ranking and spread deserve to be reported alongside bias. This thesis reports the three separately
throughout.

This thesis needs a mechanistic model for controllability. Because the canopy description enters as
an input rather than as a fitted covariate, we can change leaf quantity while holding arrangement
fixed, and change arrangement while holding quantity fixed. The collinearity that prevents
attribution in observational data is absent from a designed set of simulations, because the
experimenter, not the forest, decides which combinations are run. Attribution then follows by
construction.

### 1.3.2 The MuSICA model

MuSICA is a multilayer soil-vegetation-atmosphere transfer model that solves energy, water, and
carbon exchange through a vertically discretized canopy [@ogeeMuSICACO2Water2003]. Radiative
transfer is computed layer by layer and separately for sunlit and shaded foliage, the leaf energy
balance and stomatal conductance are solved at each level, and turbulent transfer links the layers
to the air above, while soil water and heat are treated in a coupled multilayer column. The model
was developed and evaluated at European forest flux sites across timescales from hourly to yearly
[@ogeeMuSICACO2Water2003; @ogeePartitioningNetEcosystem2003], and it has been coupled to
LiDAR-derived structure to study sub-canopy microclimate in conifer stands
[@bouwenInteractionsEntreStructure].

That study is the closest antecedent of this work. It established that a multilayer canopy model
can reproduce the ordering of sub-canopy conditions across stands of contrasting density, and that
the meteorological forcing must itself be corrected for canopy structure before scenarios of
differing structure can be compared. It also identified the constraint that most limits such
simulations: for lack of measured vertical profiles, the distribution of foliage was described by a
parametric beta function, a unimodal shape that cannot represent the stratified canopies produced by
partial harvesting, and whose use was judged to have contributed to an underestimation of the
simulated buffering. Its closing recommendation was to feed measured LiDAR-derived profiles of plant
area density into the model in place of that function.

Three things remain open. The measured profile has not been substituted for the assumed one. The
effect of the vertical arrangement has therefore not been isolated, since a profile held to a single
parametric family cannot be contrasted against an alternative arrangement at matched leaf area and
height. And the work was conducted in evergreen conifer stands, where the seasonal cycle of foliage
and the geometry of interception differ from those of a temperate deciduous canopy. It also does not
address how the structural input might be obtained beyond the LiDAR footprint, the problem taken up
in Section 1.3.3.

Three properties make MuSICA the appropriate tool for the question of Section 1.2.2.

The first is that it is multilayer. A vertical profile of leaf area density is a native input,
whereas a big-leaf scheme collapses the canopy into a single effective surface and cannot represent
the vertical dimension at all. If the arrangement of foliage is to be tested against its quantity,
the model must see both. MuSICA belongs to a small family of canopy schemes with this property,
alongside multilayer biophysical models developed for temperate deciduous stands
[@baldocchiHowEnvironmentCanopy2002], three-dimensional radiation and gas-exchange models
[@sinoquetRATPModelSimulating2001], integrated soil-canopy radiative and energy-balance models
[@vandertolIntegratedModelSoilcanopy2009], and the multilayer canopy parameterizations now being
introduced into land surface models [@bonanModelingCanopyinducedTurbulence2018;
@bonanModelingStomatalConductance2014]. Empirical microclimate downscalers, by contrast, take canopy
cover or a height metric as a scalar predictor and could not accept a profile.

![**Figure 1.2.** Schematic of what a canopy scheme can accept as input. In (a) a vertical profile
of leaf area density is a native input: the canopy is discretized into layers, each holding sunlit
and shaded foliage with its own energy balance, so the arrangement of leaf area can be altered while
its total is held fixed. In (b) the column is collapsed into one effective surface that can only be
given a scalar, so the question of Section 1.2.2 cannot be posed; layer count and profile shape are
illustrative and carry no measured values.](figures/Fig_intro_musica.pdf){width=100%}

The second is that the coupling between structure and temperature is physical. No transfer function
is fitted between the two, which removes the attribution problem of Section 1.3.1: a sensitivity
computed by perturbing an input of a physical model is a property of the physics, not an estimate
contaminated by the covariance structure of a training set.

The third is that the canopy inputs are explicit and separable. Leaf area, canopy height, cover,
and the shape of the vertical profile enter as distinct quantities, which makes controlled
perturbation possible and allows the vertical profile to be contrasted against a leaf-area- and
height-preserving uniform alternative.

One limitation bounds everything that follows. MuSICA is a one-dimensional column, so lateral
advection, edge effects, and sub-pixel canopy gaps lie outside its representation. A column knows
only the canopy directly above it, and a plot whose surroundings differ from that column, for
instance one containing a gap the mean structural input does not resolve, cannot be reproduced by
such a model. This sets a ceiling on absolute accuracy independent of how well the structural
inputs are measured, and it is why the argument developed here rests on relative comparisons
between plots rather than on absolute temperatures.

### 1.3.3 From plot to landscape

Mechanistic canopy models are plot-scale by construction. They require, for every location
simulated, a complete structural description of the canopy above that point together with the
meteorological forcing above it. Running one column is cheap; running one everywhere is limited by
inputs, not by computation.

The structural inputs are the bottleneck. Measuring LAI in the field, whether by hemispherical
photography, optical plant-canopy analyzers, or destructive sampling, is slow, weather-dependent,
and subject to its own methodological uncertainty [@weissReviewMethodsSitu2004;
@chenEvaluationHemisphericalPhotography1991; @holstMeasuringModellingPlant2004]. Different indirect
methods disagree even when applied at the same sites, and the disagreement is largest in the
understory layer [@georgeMethodComparisonIndirect2021]. Automated processing has improved
reproducibility [@brownHemiPyPythonModule2023] but not the arithmetic: measuring a vertical LAD
profile in the field is harder still, and a campaign of that kind yields tens of plots, occasionally
hundreds. A landscape contains millions of model pixels, so any route from plot-scale mechanism to
landscape-scale microclimate map passes through a structural dataset acquired remotely,
continuously, and at fine spatial resolution.

The working resolution follows from the model assumptions. A one-dimensional column assumes that
the canopy above the simulated point is horizontally uniform over the area it represents, so the
pixel must be small enough for that assumption to hold and large enough for the structural
retrieval to be reliable. Too coarse, and a single column stands for a mixture of gaps and closed
canopy whose average behaves like neither [@garriguesInfluenceLandscapeSpatial2006]; too fine, and
retrieval noise dominates the structural signal. No single resolution serves every purpose:
attributing a thermal effect to canopy structure calls for footprints large enough to define a
stable vertical profile, whereas comparing satellite against airborne retrievals is constrained to
the 10 m grid the satellite imposes. Both sit within the range over which sub-canopy temperature
varies in temperate stands [@zellwegerAdvancesMicroclimateEcology2019;
@mengeImpactsForestManagement2023]. The next section examines which remote-sensing instrument can
populate that grid, and at what cost in structural fidelity or temporal coverage.

## 1.4 Remote sensing of forest structure

### 1.4.1 Airborne LiDAR

Airborne laser scanning (ALS) is the only operational technique that resolves the interior of a
forest canopy over an area. Laser pulses penetrate gaps in the foliage and return a
three-dimensional point cloud; after the ground surface has been identified and the cloud
normalized to height above ground, the vertical gap fraction can be computed layer by layer and
inverted, through a Beer-Lambert formulation, into a profile of plant area density and an
integrated LAI [@bouvier7_GeneralizingPredictiveModels2015;
@richardsonModelingApproachesEstimate2009; @zhengRetrievalEffectiveLeaf2013]. The approach has been
established for two decades and validated across forest types and sensors
[@lefskyLidarBlackwellScienceLtdRemote2002; @morsdorfEstimationLAIFractional2006;
@chenUsingLidarEffective2004; @kwakEstimationEffectivePlant2010; @peduzziEstimatingLeafArea2012],
and the reconstruction of full LAD profiles rather than integrated totals is now routine
[@linRetrievalEffectiveLeaf2016; @kamoskeLeafAreaDensity2019; @arnqvistRobustProcessingAirborne2020].

The retrieval is not parameter-free. The conversion from gap fraction to leaf area rests on an
extinction coefficient that encodes leaf angle distribution and is commonly fixed at a nominal
value, and on an assumption of random foliage distribution that real canopies violate through
clumping [@huUsingAirborneLaser2018]. Acquisition geometry leaves its own imprint: flying height,
scan angle, and pulse density all affect penetration depth. Uncertainty in voxel-scale plant area
density is therefore substantial and should be propagated rather than ignored
[@pimont5_EstimatorsConfidenceIntervals2018], and canopy height models derived from the same cloud
carry their own processing choices [@khosravipourGeneratingPitfreeCanopy2014]. These limitations do
not invalidate the method. They mean that a LiDAR LAI is a model-based estimate resting on stated
assumptions, not a direct measurement, so neither side of a LiDAR-optical comparison is a reference
truth.

ALS also delivers more than leaf area. The same point cloud yields a digital terrain model, a
canopy height model, fractional cover, and the height and density metrics on which forest inventory
modeling rests [@khosravipourGeneratingPitfreeCanopy2014;
@bouvier7_GeneralizingPredictiveModels2015]. An ALS survey is best thought of as a structural
description of a stand rather than a leaf-area sensor, and canopy height in particular is easier to
retrieve robustly than leaf area and correlated with it, which makes it a natural fallback where a
full profile is unavailable. Terrestrial and portable scanning resolve structure in finer detail,
reconstructing leaf area density at the scale of individual crowns [@hosoiVoxelBased3DModeling2006;
@liEstimatingLeafArea2017], but they cover plots rather than landscapes. Airborne scanning covers
an area while still seeing inside the canopy.

ALS-derived structure has become a standard input to microclimate studies, whether to map buffering
directly [@grilUsingAirborneLiDAR2023; @vandewieleMappingSpatialMicroclimate2023], to constrain it
jointly with topography [@juckerCanopyStructureTopography2018], or to explain the thermal shelters
used by wildlife [@melinMooseLcesAlces2014]. For the question of Section 1.2.2 its value is
specific: it is the sensor that delivers the vertical arrangement, not a proxy for it.

The cost is spatio-temporal. ALS is acquired by aircraft campaign, which makes it expensive,
spatially bounded by the flight plan, and, in practice, single-date. A survey describes the canopy
as it was on one summer afternoon; it says nothing about budburst, senescence, or the following
year, and it does not extend beyond its footprint.

### 1.4.2 Optical satellite imagery

Sentinel-2 presents the opposite trade-off. The two-satellite constellation images every land
surface every five days in thirteen spectral bands, four delivered at 10 m and six at 20 m, among
which three red-edge bands sensitive to canopy chlorophyll and leaf area
[@druschSentinel2ESAsOptical2012; @delegidoEvaluationSentinel2RedEdge2011]. The mission was designed
for operational vegetation monitoring, and leaf area index is one of its standard biophysical
products. The mixed resolution matters for forest work: restricting an inversion to the 10 m bands
buys spatial detail at the price of leaving the red-edge and shortwave infrared information unused,
while including the 20 m bands does the reverse.

Retrieval proceeds by inverting a canopy radiative-transfer model against the observed reflectance.
The standard model couples a leaf optical model with a canopy scattering scheme, giving PROSAIL
[@jacquemoudPROSPECT+SAILModelsReview2009; @verhoefLightScatteringLeaf1984;
@feretPROSPECTDModelingLeaf2017; @feretProsailPROSAILLeaf2024]. Inverting it directly is ill-posed,
since different parameter combinations produce nearly identical spectra, so the inversion is
regularized by prior information on the parameter ranges [@combalRetrievalCanopyBiophysical2002].
Two operational strategies follow: the neural network distributed with the Sentinel toolbox is
trained on a simulated database built from prescribed parameter distributions
[@weissS2ToolBoxLevel22020], and hybrid inversion trains a machine-learning regressor on a
comparable look-up table [@verrelstMachineLearningRegression2012;
@verrelstExperimentalSentinel2LAI2015; @verrelstQuantifyingVegetationBiophysical2019]. Both have
been applied to forests [@chrysafisRetrievalLeafArea2020; @fernandesNotJustPretty2024] and
validated against field data [@brownValidationBaselineModified2021]. The result is wall-to-wall,
repeated, and free.

Producing a usable time series requires an operational chain. Top-of-atmosphere reflectance must be
atmospherically corrected, clouds and their shadows masked, and observation and illumination
geometry recorded, since the inversion depends on them [@feretPreprocS2PreprocessingSentinel22024].
The surviving observations are irregularly spaced, because cloud cover in temperate Europe removes
a large share of the nominal five-day revisit, so a continuous trajectory has to be reconstructed
by gap-filling and smoothing before it can force anything [@atzbergerTimeSeriesMonitoring2011;
@liangUsingEnhancedGapFilling2023]. The satellite therefore delivers a smoothed estimate of a
seasonal trajectory rather than a set of raw measurements. Its principal advantage over an airborne
survey lies in that trajectory: the timing of leaf-out and senescence is observed per pixel rather
than assumed [@cotrina-sanchezPhenologyEuropeanForests2026].

Four constraints qualify it, and together they define the gap set out in Section 1.5.2. The first
is geometric: the signal is a two-dimensional reflectance measured from above, so it carries
information about the upper canopy and only indirectly about what lies beneath. The second is the
saturation that follows: sensitivity declines as leaf area accumulates and is exhausted somewhere
above five, the exact point depending on the retrieval, the index, and the canopy
[@gaoEvaluatingSaturationEffect2023; @aklilutesfayeEvaluationSaturationProperty2021]. Closed
temperate forest canopies sit at or beyond that limit for much of the growing season, and the
consequence is a compression rather than an abrupt ceiling: differences in leaf area between dense
stands are progressively under-expressed in the reflectance before they disappear from it. The
third is the prior information that regularizes the inversion: operational parameter ranges are
tuned primarily for crops, transferring them to closed forest canopies biases the retrieval, and
correcting the bias trades against variance [@fernandesEvidenceBiasvarianceTrade2024]. The fourth
arises from spatial heterogeneity: because the reflectance-to-LAI relationship is non-linear, a
pixel containing a mixture of canopy heights or densities does not return the average of its parts
[@garriguesInfluenceLandscapeSpatial2006; @maImpactSpatialLAI2008]. Continuous temporal coverage
therefore comes at the cost of structural depth, and the cost is largest where the canopy is
densest and most heterogeneous.

## 1.5 Scientific gaps

### 1.5.1 The attribution gap

The first gap is one of causal attribution. Canopy structural variables covary strongly in real
stands, so correlative microclimate models cannot uniquely assign the buffering effect to any one of
them: coefficients are unstable, smooth terms are concurve, and machine-learning importance
measures redistribute themselves among correlated predictors in ways that depend on the algorithm
[@dominiciUseGeneralizedAdditive2002; @nicodemusBehaviourRandomForest2010;
@gregoruttiCorrelationVariableImportance2017]. The question posed in Section 1.2.2, whether the
quantity of foliage or its vertical arrangement governs sub-canopy buffering, therefore remains open
despite an extensive literature on both sides [@kovacsStandStructuralDrivers2017;
@ehbrechtEffectsStructuralHeterogeneity2019].

The gap is methodological. Additional field data collected along the same natural gradient will not
close it, because the gradient itself confounds the predictors. Resolving it requires observations
in which the structural variables can be moved independently. Since no forest offers such
observations, they must be generated: a mechanistic model driven by controlled combinations of
canopy inputs, with the perturbations applied around each plot's own realistic canopy so that the
sensitivities remain physically meaningful.

Substituting simulation for observation must not trade a regression whose coefficients are
unidentifiable for a model whose output cannot be checked. Three properties guard against this.
First, the model's parameters are physical and are not fitted to the response being explained, so no
calibration step can inject the answer; the sensitivity of simulated temperature to leaf area is a
consequence of radiative and turbulent transfer, not of a coefficient estimated from the gradient at
issue. Second, the perturbations are small and centered on each plot's measured canopy, so every
simulated configuration stays close to one that occurs in the field and the physics is never
extrapolated into regimes where it has not been evaluated. Third, the model's ability to reproduce
observed between-plot buffering is established separately against the logger network, so the
attribution statement rests on a model whose behavior on real stands has been characterized
independently. The design provides no direct empirical test of the attribution itself; no design
can. The attribution is as good as the physics, and the physics is testable in ways a fitted
coefficient is not.

### 1.5.2 Optical versus structural LAI

The second gap concerns the measurement of leaf area itself. Airborne LiDAR and Sentinel-2 both
deliver a quantity called LAI, but they measure different physical things: one integrates
intercepted returns through the whole canopy column, the other inverts a reflectance signal that
originates mostly near the top and stops responding once the canopy is dense. The two are known to
disagree, and the disagreement is expected to grow with canopy density, but its magnitude and
drivers have not been characterized jointly for temperate deciduous forest.

At least three candidates contribute, and they are rarely assessed together. The attenuation of the
optical signal limits the canopy depth that contributes to satellite reflectance, so the satellite
may describe only the upper part of the column that the LiDAR integrates in full. The
parameterization of the radiative-transfer inversion supplies prior ranges tuned for crops rather
than forest canopies, and those priors propagate directly into the retrieved leaf area
[@combalRetrievalCanopyBiophysical2002; @fernandesEvidenceBiasvarianceTrade2024]. Horizontal
heterogeneity of the canopy surface within a satellite pixel, which combines variation in canopy
height with variation in the underlying terrain, strains the one-dimensional assumption on which
the inversion rests, and does so non-linearly [@garriguesInfluenceLandscapeSpatial2006].

The question is one of inter-sensor consistency, a prerequisite rather than a microclimate claim.
No ground-based LAI measurement coincident with the airborne and satellite acquisitions exists at
these sites, so nothing external can arbitrate between the two retrievals, and the target cannot be
absolute accuracy. What can and must be established, before a satellite LAI is used as a structural
input to a microclimate model, is what part of the canopy that number describes.

Two caveats condition everything downstream. The first concerns domain. Consistency established
between two sensors holds over the structural domain on which it was established, and a depth-based
comparison is restricted to canopy deep and closed enough to have an interior. A search over
successive canopy depths cannot converge on stands shallower than the depths it scans, so pixels
below a minimum canopy height are excluded, and a high fractional-cover threshold retains dense
canopy and keeps gaps and edges out of the pixel. Whatever agreement is achieved on that restricted
domain is silent about short, open stands, which are the stands where an optical retrieval behaves
differently.

The second caveat is that better agreement in leaf area does not propagate into better simulated
microclimate. A mechanistic canopy model responds non-linearly to leaf area, steeply in open canopy
and with strong saturation once the canopy closes, so an identical reduction in leaf-area error is
worth different amounts depending on where along that curve it applies. A retrieval tuned to
minimize leaf-area discrepancy over one part of the range can be a worse driver of the microclimate
over another, and can even reverse the ordering of plots if the tuning compresses the gradient the
model relies on. Leaf-area accuracy and microclimate skill are related but distinct objectives, and
which one a correction improves is an empirical question.

![**Figure 1.3.** Schematic of the part of the canopy column each sensor reads. A laser pulse
travels to the ground and back, so gap-fraction inversion integrates the whole profile, whereas
reflected sunlight stops carrying information once enough leaf area lies above, so an optical
retrieval responds only to the leaf area shallower than an effective optical depth. Holding that
optical path fixed while leaf area increases, from (a) to (b), leaves the satellite reading a
smaller share of the column; no axis carries measured
values.](figures/Fig_intro_sensor_column.pdf){width=100%}

### 1.5.3 The spatio-temporal trade-off

The third gap follows from the first two. ALS resolves three-dimensional structure but is acquired
once and over a limited footprint. Sentinel-2 is repeated and spatially complete but optical,
two-dimensional, and saturating. No single sensor currently delivers fine three-dimensional
structure and temporal dynamics together.

Driving a mechanistic model beyond the plots where LiDAR exists therefore forces a compromise
between structural fidelity and coverage. That compromise has usually been made implicitly, by
using whichever product was available, and evaluated by pooling all plots into a single accuracy
statistic. Neither the choice nor the evaluation has been examined against the regime dependence
that Sections 1.4.2 and 1.5.2 predict: if the optical signal fails where the canopy is dense, a
single pooled verdict on which sensor is better may conceal two opposite verdicts in the two halves
of the density gradient. Which compromise is right, and whether it is the same everywhere in a
landscape, has not been established.

The three gaps have to be addressed in order. The structural dimension the microclimate responds to
determines what a sensor needs to measure; the part of the canopy the optical sensor sees
determines whether it measures that dimension; and only once both are settled can the choice of
forcing product be adjudicated. The order of the chapters follows the order of the gaps.

![**Figure 1.4.** Schematic of the spatio-temporal trade-off between structural sensors, each placed
by how much of the three-dimensional canopy it recovers and how often it revisits; spatial coverage
is stated in words rather than encoded in symbol size. The corner a mechanistic microclimate model
would need, full structure repeated through the season, is empty. GEDI is drawn hollow: it is used
in none of the three chapters and is taken up in the General Discussion; positions are qualitative
and no axis carries tick values.](figures/Fig_intro_tradeoff.pdf){width=100%}

## 1.6 Research questions and objectives

**RQ1 (Chapter 1).** *Which dimension of canopy structure governs sub-canopy ΔTmax buffering, leaf
quantity or vertical arrangement, and can a mechanistic model attribute it despite the collinearity
of structural variables?*

The objective is to couple ALS-derived canopy structure (LAI, fCover, *H*max, and the LAD profile)
to MuSICA and to run a controlled trait-perturbation sensitivity analysis, in which each structural
variable is moved by a small, native step around every plot's own realistic canopy, and the
measured vertical profile is contrasted against a leaf-area- and height-preserving uniform profile.
Because each perturbation is applied per plot around its own operating point, the sensitivities
describe the local response of a real canopy rather than an extrapolation. Effects are stratified
by structural archetype, so the answer is allowed to differ between open and dense canopy, and the
simulated buffering is validated against an in-situ network of 53 temperature loggers. Two
hypotheses are posed: that vertical structure modulates radiative and turbulent transfer enough to
make the LAD profile a key parameter for a faithful microclimate simulation, and that, at fixed
leaf area and height, foliage concentrated high in the canopy attenuates the diurnal cycle more
strongly than foliage spread through the column.

**RQ2 (Chapter 2).** *Can Sentinel-2 optical LAI be reconciled with LiDAR-derived structural LAI, or
does optical saturation decouple them?*

The objective is to intercompare the two retrievals across three temperate deciduous forests and to
assess the three candidate drivers of their discrepancy identified in Section 1.5.2: the effective
canopy depth contributing to the satellite reflectance, the parameterization of the PROSAIL
inversion, and horizontal heterogeneity of the canopy surface within the satellite pixel. The three
are addressed in sequence, each correction applied on top of the previous one, so their individual
contributions can be read. The target is inter-sensor consistency, not absolute accuracy, since no
coincident field measurement of LAI is available to arbitrate between the two sensors. The expected
output is a quantified statement of how well, and under what canopy conditions, the two retrievals
can be brought into agreement, so that it can be carried forward as a condition of use.

**RQ3 (Chapter 3).** *Which remotely sensed description of canopy structure should force a
mechanistic microclimate model, in which canopy regime, and can a product be built that remains
deployable beyond the footprints where LiDAR exists?*

The objective is to force MuSICA with a set of competing LAI scenarios in which the leaf-area
magnitude is drawn from LiDAR, from Sentinel-2, or from their combination, while canopy height and
the shape of the vertical profile are held fixed at the LiDAR values. This design isolates the
effect of the magnitude source rather than confounding it with a change of structure. The scenarios
are adjudicated against the logger network on four independent axes: a model-free test in leaf-area
space that does not involve MuSICA, the between-plot ranking, the warm bias, and the recovery of
observed thermal amplitude. Two properties of the configuration govern how those axes are read.
Every scenario under-disperses the between-plot amplitude of the simulated buffering field, so the
coefficient of determination measures how faithfully a product ranks plots, not absolute skill; and
the adjudication is stratified by canopy density at the leaf-area value where the two sensors cross
over, so a regime-dependent answer can be detected if one exists. A specific test follows from the
second caveat of Section 1.5.2: whether an optical retrieval tuned to agree with LiDAR over closed,
tall canopy is still the better microclimate driver once applied across the whole density gradient,
including the short, open stands its tuning excluded. We hypothesize that no single sensor wins
everywhere, and that the operationally useful product is a rule for choosing between them rather
than a single best retrieval. This step conditions any subsequent landscape-scale mapping.

![**Figure 1.5.** Structure of the thesis. Each of the three gaps identified in Section 1.5 gives
rise to one research question and one research chapter. The vertical arrows on the left indicate
that each gap has to be closed before the next can be posed, which fixes the order of the
chapters.](figures/Fig_intro_roadmap.pdf){width=100%}

## 1.7 Overview of the thesis

**Chapter 1** answers RQ1. ALS-derived structure from a lowland sessile oak forest drives MuSICA at
several hundred plots stratified into four structural archetypes, from open to dense, and each
structural variable is perturbed around every plot's own canopy so that its contribution to
simulated ΔTmax is isolated by construction. The answer is density-dependent, which reconciles the
conflicting literature of Section 1.2.2 rather than adjudicating it. Leaf quantity leads at both
ends of the density gradient, and its effect steepens as the canopy closes. The vertical profile
matters most in intermediate canopy: reported as a complete real-versus-uniform swap, it becomes
the largest single effect there, an upper bound that real canopies barely span, and in the densest
canopy it is again dwarfed by leaf quantity. Because the profile is a function rather than a
scalar, it is isolated by a categorical contrast against a leaf-area- and height-preserving uniform
profile, against which the realistic profile buffers slightly less. A separate stand-level
validation against the 53 loggers establishes the model's credibility rather than testing the
per-plot sensitivities directly, since the perturbation design points are not the logger plots. A
model-free check on the same loggers asks the observational form of the question: whether adding
the measured profile to leaf quantity and cover improves the fit to the observed buffering. The
validation shows that the model reproduces the between-plot buffering ranking well, while carrying
a warm bias on the absolute offset and compressing the micro-macro slope toward unity, so that both
buffering and amplification are under-expressed in magnitude; the residual traces mainly to
sub-pixel canopy gaps that a one-dimensional column cannot represent. Because leaf quantity is what
the microclimate responds to, the next question is whether a satellite can measure it.

**Chapter 2** answers RQ2. Sentinel-2 and LiDAR LAI are intercompared across three temperate
deciduous forests, a montane beech stand and two lowland oak-dominated ones, and the three
candidate drivers of their discrepancy are addressed in sequence. Satellite LAI agrees best with
the LiDAR estimate when the latter is restricted to the upper 6 to 10 m of the canopy, depending on
the site. This effective canopy depth is an empirical proxy, selected by a multi-criteria
compromise between weakly differentiated criteria: it names the depth at which the two sensors are
most consistent, not a measured optical penetration depth, and at one site the individual criteria
disagree enough that the selected value is a statistical compromise and no more.

Adapting the PROSAIL parameter ranges to forest canopies through a multi-criteria Pareto selection
then improves agreement, though unevenly: the gain is large at one site and marginal at another,
and part of what it recovers is agreement lost at the preceding depth-restriction step. The
analysis separates the parameters that can keep their default values from those that need
forest-adapted priors, and it bounds the extent to which the improvement depends on LiDAR-derived
priors, showing that most of it survives when those priors are replaced by LiDAR-free ones.
Residual discrepancies remain largest in stands whose canopy surface is horizontally heterogeneous
within the satellite pixel, so vertical thickness and horizontal heterogeneity have to be accounted
for jointly, and a structural overestimation of leaf area in sparse stands persists that no
parameter tuning removes. What improves is inter-sensor consistency, not demonstrated accuracy
against a field truth.

The domain of validity bounds the result. The depth-based comparison is posed only where the canopy
has an interior, so the analysis retains pixels of high fractional cover and of canopy height above
the depth range being probed, and excludes the short, open stands below them. The optimized
parameterization is calibrated on closed, tall canopy and carries no guarantee outside it.
Chapter 3 reads these results as a statement of where the optical signal can be trusted, and uses
them in that role rather than as a correction to be applied everywhere.

**Chapter 3** answers RQ3. MuSICA is forced at the logger plots, not over a landscape, with
competing LAI products, and the sensor verdict depends on canopy regime. Where the canopy is dense
and the optical signal has saturated, LiDAR carries the between-plot ranking and Sentinel-2 carries
none of it; where the canopy is open, Sentinel-2 leads, though by a narrower margin on the
buffering slope than on the daytime offset. The apparent pooled advantage of one sensor over the
other dissolves once the plots are stratified at a near-median leaf-area split. The two halves of
that reversal do not carry equal weight: the dense-canopy dominance of LiDAR is robust to the
model's boundary-layer configuration, whereas the open-canopy usability of Sentinel-2 is the
configuration-dependent link. The mechanism is the one diagnosed in Chapter 2: most of the leaf
area of a dense canopy lies below the effective optical depth and is invisible from above, which
also explains why truncating the LiDAR column to that depth degrades the simulated buffering rather
than improving it.

Consistency in leaf area and skill in microclimate come apart, as Section 1.5.2 anticipated. The
forest-tuned optical retrieval optimised in Chapter 2 is the product most consistent with LiDAR,
yet it is the worst microclimate driver tested, inverting the sign of the between-plot
relationship. Two things account for this, and neither is a defect of the optimisation. It was
calibrated on closed, tall canopy, so applying it across the full density gradient extrapolates it
into the open, short stands its calibration domain excluded. And because a mechanistic canopy model
responds non-linearly to leaf area, a tuning that reduces mean discrepancy can flatten the
between-plot contrast the model needs, which costs ranking even as it buys accuracy. Chapter 2 is
therefore consumed as a map of where each sensor is valid, not as a corrector to be applied. The
practical outcome is a regime-dependent fusion, switched on canopy height so that no per-pixel
LiDAR is needed to decide which product to use, although the forcing itself still draws canopy
height and profile shape from LiDAR in both regimes. Correcting the optical retrieval where it does
not saturate makes matters worse, because the switch has already placed it where its saturation
does not bite. The chapter delivers inter-sensor consistency and an explicit domain of validity,
not the elimination of LiDAR: it maps the trade-off of Section 1.5.3 rather than resolving it, and
each sensor is used where it is valid.

The three chapters form one argument. Leaf quantity is what the microclimate responds to, with the
vertical arrangement a bounded, second-order effect that peaks in intermediate canopy (Chapter 1).
The optical sensor that could supply leaf quantity continuously does not see all of it, and the
part it misses is the part that dense canopies hold (Chapter 2). The consequence for mechanistic
microclimate modeling is a regime-dependent choice of structural input rather than a single best
sensor (Chapter 3). Each chapter narrows the question the next one asks.

One boundary of scope comes first: this thesis does not produce a wall-to-wall microclimate map. It
establishes what such a map would have to be built on: which structural dimension the temperature
responds to, what a satellite can and cannot see of that dimension, and which product should drive
the model in which canopy regime. Mapping remains out of reach for two reasons established in the
chapters rather than assumed. The forcing still requires LiDAR structure wherever the canopy is
dense, so coverage is bounded by the flight footprint; and a one-dimensional column cannot
represent the gaps, edges, and lateral exchanges that a landscape is made of. The General
Discussion starts from this boundary.

Four further limits apply. The microclimate work rests on a single lowland oak site instrumented
with 53 loggers, of which the chapters use the summer 2021 window, so the operational thresholds
proposed in Chapter 3, in particular the canopy height at which the fusion switches, are calibrated
in-sample and will require recalibration elsewhere; because the loggers are not spatially
independent, the confidence intervals Chapter 3 derives from them are optimistic. The microclimate
analysis is confined to the summer window, when the coupling between leaf area and the sub-canopy
thermal regime is expressed and the logger series are comparable across plots. The seasonal
question is addressed but the answer is bounded: within summer a temporally resolved satellite
series adds nothing, for a mechanical reason, since leaf area is on its phenological plateau and
there is no timing left to resolve; outside summer the chapters do not arbitrate the question. The
logger archive extends through the following seasons and the 2022 heatwaves, and the General
Discussion draws on it for perspectives only. On the optical side, no ground LAI coincident with
the acquisitions is available, so Chapter 2 establishes consistency between sensors rather than
accuracy against a truth. And the mechanistic model is one-dimensional throughout, which leaves
edges, lateral advection, and sub-pixel gaps outside its representation. Each limit marks where the
next study should begin.

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
