# Appendix G. Exploratory analysis: GEDI as the temporal LiDAR dimension

*(Evidence base for the General Discussion perspectives. Exploratory status:
single summer of ALS truth, microclimate tests at Blois only. All scripts and
tables in `NC_Full/manuscripts/ch4/` and the simulation repository; ~30 min of
compute end to end.)*

## G.1 Data and processing

We used the GEDI L2A/L2B footprints already extracted over the three sites
(9190 footprints, 2019–2023 acquisitions restricted here to 2021–2022: Aigoual
516, Blois 2448, Mormal 6226), each carrying relative-height metrics, total
plant area index (PAI), the vertical PAI profile in 5 m bins
[@dubayahGlobalEcosystemDynamics2020], geolocation-corrected coordinates, and
ALS LAI/LAD extracted at the corrected footprint position (k = 0.5, rescaled by
0.769 to the study convention k = 0.65). Quality filtering follows the
literature and our own sweep: full-power beams only, sensitivity ≥ 0.90, PAI > 0,
ALS LAI ≥ 2 m² m⁻². The beam filter is the single most effective lever: at Blois
the footprint-level correlation with ALS LAI rises from *r* = 0.49 (all beams)
to *r* = 0.67 (power beams, *n* = 635). Because successive GEDI dates sample
different footprints, all temporal composites are adjusted to a constant
structure by regressing PAI on acquisition date and ALS LAI, and reported at the
common mean ALS LAI; a common-support control gives the same curves.
Phenological metrics use the 50% seasonal-amplitude threshold (TRS50) on
weekly-interpolated, Whittaker-smoothed series [@atzbergerTimeSeriesMonitoring2011],
the metric of @cotrina-sanchezPhenologyEuropeanForests2026; the Elmore
double-logistic misbehaves on one pooled season of sparse points and is kept
only as a control.

Microclimate tests reuse the Chapter 3 machinery unchanged: MuSICA
[@ogeeMuSICACO2Water2003] at the 53 HOBO plots of Blois, v3.2.3-iter binary,
per-plot wind correction, total (two-sided) leaf area, scenarios validated
against observed ΔTmax and the hourly micro–macro slope.

## G.2 A waveform seasonal cycle that ranks canopy density

Monthly composites at Mormal (10 months sampled) trace a full deciduous cycle:
winter floor 0.5–0.8 m² m⁻² (December–March, wood and branch area), leaf-out by
May, summer plateau 3.4–3.6, decline through October (Fig. G1). The seasonal
amplitude, a leaf-only proxy since the wood floor subtracts out, increases
monotonically across ALS-LAI strata (2.3 / 3.1 / 3.3 m² m⁻² for strata 2–3 /
3–4 / >4), and the summer plateau ranks density at both sites (Fig. G2). This is
the temporal counterpart of the Chapter 2 result: where the optical signal
saturates, the waveform amplitude still discriminates.

![Figure G1. GEDI monthly composites (power beams, structure-adjusted; ribbon =
IQR; dashed = adjusted). Orange: the Sentinel-2 ATBD annual series at Blois.
Winter floor = wood/branch area.](figs_gedi/Fig6_c4_gedi_seasonal.png)

![Figure G2. Winter floor, summer plateau, and seasonal amplitude by ALS-LAI
stratum. The amplitude ranks canopy density where the summer optical signal is
flat.](figs_gedi/Fig7_c4_seasonal_amplitude.png)

## G.3 Timing: leaf-out confirmed, senescence arbitrated

At Blois, per-acquisition composites adjusted to constant structure (nine acquisitions
forming eight composite points, 2021 and 2022; Fig. G3) yield three timing results. (i) The 2021 summer plateau
was not flat: June 6 sat 0.68 ± 0.41 m² m⁻² below September (*z* = −3.2), and
GEDI lands on the Sentinel-2 June values, confirming that the S2 dip reflected
the late, cold spring of 2021 rather than inversion noise. (ii) The 2022 drought
summer left no canopy signature (September 2022 slightly above September 2021 at
matched structure), while May 14, 2022 was already at seasonal maximum: GEDI
resolves interannual phenology that a single-summer ALS snapshot cannot.
(iii) TRS50 phenometrics (Table G1) put the GEDI start of season within three
days of Sentinel-2 (122 vs 125), but the end of season about 25 days later
(34 days with the unsmoothed spline metric of the lab notebook), replicating at site level the ~35-day continental gap of
@cotrina-sanchezPhenologyEuropeanForests2026. Layer-resolved composites
(Fig. G4), transposing the method of @oliveiraUpperCanopyUnderstory2025, show
the understory (0–10 m) greening up about seven days before the overstory
(TRS50 116 vs 123), a measured vernal window.

![Figure G3. Blois per-acquisition GEDI composites, adjusted to a common ALS
LAI (95% CI), against the Sentinel-2 2021 series.](figs_gedi/Fig8_c4_blois_summer.png)

![Figure G4. Layer-resolved GEDI phenology (understory 0–10 m vs upper canopy),
structure-adjusted monthly composites normalized by each layer's seasonal
maximum.](figs_gedi/Fig9_c4_layer_phenology.png)

Table G1. TRS50 phenometrics (day of year; weekly-interpolated,
Whittaker-smoothed series).

| Series | SOS | EOS |
|---|---|---|
| Blois, Sentinel-2 ATBD series (19 dates, 2021) | 125 | 279 |
| Blois, GEDI adjusted composites (eight composite points) | 122 | 304 |
| Mormal, GEDI total (10 months) | 123 | 305 |
| Mormal, GEDI understory (0–10 m) | 116 | 298 |
| Mormal, GEDI upper canopy (>10 m) | 123 | 312 |

Two robustness observations complete the timing picture. The phenometrics
transfer across sites: SOS 122 vs 123 and EOS 304 vs 305 between Blois and
Mormal (~250 km apart), although the full leaf-fraction curves overlap only
moderately (*r* = 0.74, maximum difference 0.28), because the Blois curve
carries the real 2021 spring anomaly that the Mormal months do not sample; the
shoulders transfer, the mid-season carries local interannual signal. And at
Mormal the TRS50 end of season is invariant across canopy-density strata
(305 / 305 / 307 for ALS-LAI 2–3 / 3–4 / >4) while the start of season comes
about 12 days earlier in the densest stratum (111 vs 123), an observation to
read with species-composition caveats but directly relevant to the
structure–phenology–microclimate coupling of @wuCanopyStructureRegulates2024.

## G.4 The arbitration matters for microclimate

We cloned the Chapter 3 seasonal scenario (LiDAR magnitude × Sentinel-2 timing)
and replaced only the autumn limb (day of year ≥ 250) with the GEDI-measured
leaf fraction. Against the 53 loggers, the GEDI-timed autumn significantly
improved both metrics over the S2-timed series (paired bootstrap, 47 plots:
ΔRMSE −0.083 °C [−0.105, −0.060]; Δ*R*² +0.188 [+0.064, +0.404]), turning the
series' worst window (autumn *R*² 0.10, below a constant-LAI control) into its
best RMSE. A second, pre-registered-in-spirit correction, filling the summer
dips of the S2 series up to the plateau while keeping the peak, added a smaller
but significant summer gain (JJAS ΔRMSE −0.066 [−0.085, −0.048] vs the raw
series). Month-by-month decomposition attributes the summer effect to September
and August, confirms June 2021 as real signal, and flags the mid-summer S2
wobble as noise, exactly as the GEDI plateau suggested (Table G2).

Table G2. Paired-bootstrap contrasts against the S2-timed seasonal series
(47 plots with a full series, *B* = 3000; ΔRMSE in °C on ΔTmax; brackets =
95% CI; source tables 15–18 of the working folder).

| Contrast | Window | ΔRMSE | Δ*R*² |
|---|---|---|---|
| GEDI autumn limb | Oct–Nov | −0.083 [−0.105, −0.060] | +0.188 [+0.064, +0.404] |
| GEDI autumn limb | JJAS (control) | −0.021 [−0.032, −0.012] | +0.022 [+0.006, +0.044] |
| + summer dip-filling | JJAS | −0.066 [−0.085, −0.048] | +0.035 [+0.014, +0.060] |
| + summer dip-filling | August alone | −0.158 [−0.211, −0.109] | +0.273 [+0.135, +0.488] |

## G.5 One irreplaceable role per sensor

A 15-scenario matrix ({ALS, S2, GEDI} × {static, temporal}, same machinery,
Table G3) condenses the exploration.

Table G3. Sensor-combination matrix, summer (JJAS) and autumn, condensed to
the informative rows (full 15-scenario table in the working folder, Table 20).
ΔTmax metrics vs 53 loggers; pooled summer *R*² inherits the open-canopy
composition effect of Chapter 3 and is read jointly with RMSE and bias.

| Combo (type) | Summer *R*² / RMSE / bias | Autumn *R*² / RMSE |
|---|---|---|
| S2 (static) | 0.56 / 1.75 / 0.79 | 0.31 / 1.63 |
| ALS (static) | 0.24 / 1.76 / 0.73 | 0.35 / 1.63 |
| ALS+S2, Ch3 layered fusion | 0.50 / **1.66** / 0.79 | 0.34 / 1.63 |
| ALS+S2 (temporal, S2-timed) | 0.24 / 1.76 / 0.73 | 0.10 / 1.66 |
| ALS+GEDI (temporal, GEDI shape) | 0.15 / 1.91 / 0.89 | 0.30 / **1.58** |
| S2+GEDI, no ALS (temporal) | 0.18 / 1.84 / 0.84 | 0.33 / 1.63 |
| GEDI only (temporal) | 0.01 / 2.13 / 1.16 | 0.29 / 1.64 |
| ALS+S2+GEDI (temporal, corrected) | 0.27 / 1.72 / **0.68** | 0.30 / 1.60 |

All low-bias scenarios carry the ALS magnitude; every temporal variant where
GEDI supplies the per-plot magnitude fails in summer (*R*² ≤ 0.18), and the two
static GEDI-magnitude scenarios score higher (0.32 direct at footprints, 0.49
ratio-scaled Sentinel-2) but carry the largest warm biases of the matrix
(+1.08, +1.23 °C), because GEDI's dense-canopy skill,
strong at its footprints (*r* = 0.67; bottom-layer agreement *r* = 0.84 where
Sentinel-2 is blind), does not regress beyond them through wall-to-wall
features. A site-level GEDI shape cannot replace the plot-scale Sentinel-2
timing either (summer *R*² 0.15). What GEDI adds, uniquely, is the seasonal
arbitration: both GEDI-timed scenarios take the best autumn RMSE, and the
triple combination (ALS magnitude + S2 timing + GEDI correction) is the best
balanced package (best summer bias 0.68 °C, second RMSE in both windows).
For the 11–18 plots within 100–150 m of a power footprint, the mean PAI of
the three nearest footprints, used directly, recovers 80–90% of the ALS
microclimate correlation (*r* = −0.70 to −0.73 / −0.77 to −0.80 against
−0.86 / −0.89 to −0.90 for plot ALS), decaying beyond 200 m at the pace of the forest's
own spatial decorrelation.

The composition closes the exploration. Applying the GEDI temporal correction
inside the Chapter 3 height-gated fusion (tall plots: ALS magnitude ×
GEDI-corrected series; short plots: unchanged) produced the best product tested
in either window: summer RMSE 1.61 °C (previous best 1.66), autumn *R*² 0.46
and RMSE 1.57, with significant paired gains over the Chapter 3 fusion in both
windows (summer ΔRMSE −0.066 [−0.084, −0.048], Δ*R*² +0.034 [+0.013, +0.056];
autumn ΔRMSE −0.082 [−0.104, −0.058], Δ*R*² +0.194 [+0.064, +0.406]). The
Chapter 3 spatial rule and the GEDI temporal arbitration are complementary and
compose without interference.

## G.6 The 2022 heat test: the forcing rule transfers to an extreme year

The HOBO network stayed active through the record June and July 2022 heatwaves,
and the simulations of the scenario machinery run to late July 2022, so the
whole validation transfers to an out-of-calibration extreme year at no cost
(Table G4). Two results. First, the observed LAI–buffering coupling holds under
extreme heat: *r*(slope, LAI) = −0.88 on the p90 hottest days of June–July 2022
(−0.88 on the p90 hottest days of 2021 as well; −0.92 over all of summer
2021), the mean micro-macro slope does not collapse (0.93 against 0.86 on the
hottest 2021 days), and the mean offset on hot days is deeper (−1.13 °C on
2021 hot days, −1.22 °C on 2022 hot days, against −0.73 °C over all of summer
2021), as a sub-unity slope under hotter skies implies; the GEDI composites
show the 2022 canopy itself intact (§G.3). The invariant coupling, not the
deeper offset, is the evidence relevant to the microrefugia premise
[@lenoirClimaticMicrorefugiaAnthropogenic2017; @defrenneGlobalBufferingTemperatures2019].
Second, the scenario ranking transfers: the two fusion products dominate the
2022 window (ΔTmax *R*² 0.40 and 0.39, RMSE 1.80–1.81 °C) while the static ALS
scenario collapses (*R*² 0.08). One honest nuance: the warm bias grows in the
extreme year (+1.16 °C vs +0.73 °C in 2021), so the model under-buffers more
when heat is most severe, a clear target for model improvement.

Table G4. Observed coupling and scenario validation in the 2022 heatwave
window (June 1 – July 29, 2022; *n* = 51 plots; 2021 rows *n* = 53).

| Quantity | 2021 JJAS, all days | 2021 JJAS, hot p90 | 2022 Jun–Jul, all days | 2022 Jun–Jul, hot p90 |
|---|---|---|---|---|
| Mean observed ΔTmax (°C) | −0.73 | −1.13 | −1.05 | −1.22 |
| *r*(slope, ALS LAI) | −0.92 | −0.88 | −0.89 | −0.88 |
| Mean micro-macro slope | 0.90 | 0.86 | 0.90 | 0.93 |

| Scenario (2022 window) | ΔTmax *R*² | RMSE (°C) | Bias (°C) |
|---|---|---|---|
| Ch3 fusion (FUSION_H) | 0.40 | 1.81 | +1.17 |
| + GEDI correction (FUSION_GEDIMAX) | 0.39 | 1.80 | +1.16 |
| Sentinel-2 static | 0.27 | 2.12 | +1.41 |
| ALS static | 0.08 | 2.11 | +1.32 |

## G.7 The vernal window in degrees: which way the association points

Per-plot green-up dates were taken from the full-year Sentinel-2 series of
each plot (the annual ATBD series of Chapter 3, 47 plots) with the smoothed
half-amplitude threshold of Table G1 (weekly interpolation, Whittaker
smoothing, last upward crossing before the seasonal peak); they fall between
day of year 106 and 134. An earlier version of this test used the summer-only
series and a raw first crossing, which dated half the plots in January; those
numbers are withdrawn. Crossing the corrected dates with the plot's own
observed spring temperature (mean daily maximum, March 15 – April 30, 2021;
*n* = 33 plots with both) gives *r* = +0.55, and +0.40 after partialling out
ALS LAI (+0.33 after LAI and maximum height). Plots that leaf out later stay
open longer and accumulate a warmer understory spring. The reciprocal test
(pre-season winter temperature, January 1 – March 14, against green-up, in the
spirit of @wuCanopyStructureRegulates2024; *n* = 35) gives *r* = +0.40 but
+0.19 once ALS LAI is partialled out, so most of the winter association runs
through leaf area. Green-up date is itself correlated with ALS LAI (*r* =
−0.67, denser plots dated earlier), which is why the partial coefficients are
the ones to read. Both tests are cross-sectional and share their confounds
(April canopy openness, composition, topography); together they say that the
spring understory temperature follows phenology more closely than winter
temperature anticipates it, and they connect the GEDI understory lead of §G.3
to the energy balance it implies.

## G.8 Ancillary logger analyses (night, and how many loggers)

Two observation-only analyses of the 2021 logger archive complete the
evidence base of the General Discussion. First, the night. Computing ΔTmin
(daily-minimum micro minus macro) alongside ΔTmax over summer 2021 gives a
near-mirror pair of couplings to LiDAR leaf area: *r* = −0.89 for ΔTmax and
+0.83 for ΔTmin (*n* = 53). Closed canopy (LAI > 4) holds the night at the
reference temperature (mean ΔTmin +0.09 °C) while open stands (LAI < 1.5) fall
2–4 °C below it (−3.15 °C), the mirror of the daytime pattern and a bilateral
amplitude compression under a single structural control; the mean ΔTmin is
−0.66 °C with a between-plot SD of 1.24 °C, so the night field carries real
structure. The macro minimum is taken from the MuSICA forcing file, so
cold-air pooling at the reference is not excluded. This is a consistency check against the
Chapter 1 attribution, not an attribution itself.

Second, network size. Bootstrap subsampling of the 53 loggers shows the
observed slope-LAI coupling is establishable from small networks: at
*n* = 10 the median *r* is −0.94 (95% CI −0.98 to −0.79), at *n* = 20 the CI
tightens to −0.96/−0.86 and the regression coefficient is constrained to
about ±25–30%, against ±11% at *n* = 40. Within this network, fifteen to
twenty well-spread loggers reproduce the coupling; whether that transfers to
another site is untested, since the subsample inherits the site's own gradient.

## G.9 Limits

Seven caveats bound these perspectives. The microclimate tests are single-site
(Blois) and the composites pool 2021–2022, so interannual phenology is both
signal and confound. Aigoual is a stated failure case for the footprint route:
516 footprints (203 matched), September only, relief and conifer admixture;
footprint-level *r*(PAI, ALS LAI) = 0.13 (all beams) and 0.22 (power beams),
and a GEDI-anchored magnitude correction over-corrects there (bias
−0.90 m² m⁻²). All microclimate results are Blois-only by construction. The late GEDI end of season mixes true leaves, marcescent
leaves (typical of *Quercus* and *Fagus*), and exposed wood
[@cotrina-sanchezPhenologyEuropeanForests2026]; the logger test suggests this
structural signal is nevertheless the functionally relevant one for radiation
and microclimate. Successive acquisitions sample different footprints, handled
by structure adjustment, not by repeated measures. Several scenario
variants were evaluated on the same logger network (47 series plots in the
paired contrasts); the large effects (August,
autumn) were pre-specified by the GEDI diagnosis, the small ones (July,
September) were not and should be treated as suggestive. The 2022 test recycles
the 2021 seasonal curves (the summer plateau makes this benign, and the GEDI
composites support an unchanged 2022 canopy, but the leaf-out of 2022 was
earlier and is not represented). And the vernal-window correlations (§G.7, *n* = 33 and 35)
are cross-sectional: they order two associations, spring temperature after
green-up and winter temperature before it, and do not identify a feedback.
