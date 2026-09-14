# General Introduction — LOCKED outline

**The three research chapters** (the General Introduction builds toward these):
- **Chapter 1** — ALS-derived canopy structure → MuSICA → *mechanistic attribution* of sub-canopy microclimate
  buffering (leaf quantity vs vertical arrangement; density-dependent), validated against in-situ loggers.
- **Chapter 2** — Sentinel-2 vs LiDAR LAI: *optical saturation / the "optical illusion"* — does the optical
  signal carry the structural information that governs microclimate?
- **Chapter 3** — *Spatialisation* of fine-scale microclimate from satellite remote sensing (Sentinel-2
  time-series), navigating the spatio-temporal trade-off; plot → landscape upscaling.

**Red thread (stated early — end of 1.1, head of 1.4):** *can mechanistic, plot-scale microclimate be mapped
over landscapes from remote sensing, when every sensor imposes a trade-off — ALS resolves 3-D structure but
only once and locally; Sentinel-2 is repeated but optical and saturating?* The three chapters are three
responses to this single tension. *(Spaceborne LiDAR is reserved for the thesis Discussion as a future avenue,
not treated in any chapter.)*

---

## 1.1 The ecological importance of forest microclimates in a changing climate
- **1.1.1** Climate change and the increasing frequency of temperature/weather extremes.
- **1.1.2** The forest canopy as a physical modifier of the sub-canopy climate: **buffering *and* amplifying**
  effects (amplifying retained — it sets up the amplifying plots / glass ceiling of Chapter 1).
- **1.1.3** Why fine-scale temperature variability matters for understory biodiversity, tree regeneration and
  thermal microrefugia.
- **1.1.4** Observing microclimate: in-situ logger networks (HOBO; SoilTemp-type efforts) and their intrinsic
  **spatial limitation** — point measurements cannot map a landscape → motivates modelling and remote sensing,
  and provides the validation backbone of Chapter 1.

## 1.2 Drivers of forest microclimate and the role of canopy structure
- **1.2.1** Linking 3-D vegetation structure to sub-canopy microclimate (radiative and turbulent transfer).
- **1.2.2** *Which dimension of structure?* — **total leaf quantity (LAI) vs the vertical arrangement of
  foliage (LAD profile)**, framed as an **open question** (not "LAI is the key metric"; Chapter 1's result is
  density-dependent and must not be pre-judged). The vertical profile (LAD) is introduced here alongside LAI.

## 1.3 Modelling the forest microclimate
- **1.3.1** Empirical/statistical versus mechanistic approaches; the limits of correlative models under trait
  collinearity (sets up the attribution gap, 1.5.1).
- **1.3.2** The MuSICA model: mechanistic energy, water and radiative transfer in a multilayer canopy.
- **1.3.3** From plot to landscape: mechanistic models are plot-scale; spatialising them requires **continuous
  structural inputs** → motivates remote sensing (1.4) and the upscaling problem of Chapter 3.

## 1.4 Remote sensing of forest structure: bridging the spatial scale
*(Each subsection closes on the sensor's spatio-temporal signature, so the trade-off gap 1.5.3 lands.)*
- **1.4.1** Airborne LiDAR (ALS): high-resolution 3-D structural mapping — but local and acquired once.
  *(Input to Chapter 1.)*
- **1.4.2** Optical satellite imagery (Sentinel-2): continuous, repeated time-series — but a 2-D optical signal
  prone to **saturation**. *(Input to Chapters 2 and 3.)*

## 1.5 Identifying the scientific gap
- **1.5.1 [→ Chapter 1] The attribution gap.** Structural traits are strongly collinear, so statistical /
  correlative microclimate models cannot uniquely *attribute* buffering to a trait (unstable coefficients,
  concurvity); which structural dimension — quantity vs vertical arrangement — drives the effect is unresolved
  → calls for a mechanistic, all-trait-combinations attribution.
- **1.5.2 [→ Chapter 2] Optical versus structural LAI.** Discrepancies between optical (Sentinel-2) and
  structural (LiDAR) estimates of LAI, and optical saturation in dense canopy → does the optical signal carry
  the structural information that governs microclimate, or only an illusion of it?
- **1.5.3 [→ Chapter 3] The spatio-temporal trade-off.** ALS gives fine 3-D structure but is not repeated;
  Sentinel-2 is repeated but optical and saturating; no single sensor delivers 3-D structure *and* temporal
  dynamics → mapping mechanistic microclimate forces a compromise.

## 1.6 Research questions and objectives
- **RQ1 / Chapter 1:** *Which dimension of canopy structure governs sub-canopy microclimate buffering — leaf
  quantity or vertical arrangement — and can a mechanistic model attribute it under trait collinearity?*
  → Objective: attribute MuSICA-simulated ΔT_max to LiDAR-derived traits with an exact-Shapley decomposition
  over all trait combinations, validated against in-situ loggers.
- **RQ2 / Chapter 2:** *Can an optical satellite signal (Sentinel-2 LAI) substitute for LiDAR-derived
  structure in representing the microclimate-relevant canopy, or does optical saturation decouple them?*
  → Objective: intercompare Sentinel-2 and LiDAR LAI and their respective microclimatic signal.
- **RQ3 / Chapter 3:** *Can fine-scale forest microclimate be spatialised over the landscape from satellite
  remote sensing, and at what accuracy cost given the spatio-temporal trade-off?*
  → Objective: upscale the mechanistic microclimate using remote-sensing structural inputs (Sentinel-2
  time-series).

## 1.7 Overview of the thesis
- One paragraph per research chapter (1, 2, 3), each tied to its research question and to the red thread.
- Closing perspective: spaceborne LiDAR (GEDI) as a future avenue for spaceborne 3-D structure — currently
  limited by spatial discontinuity and geolocation uncertainty — **reserved for the thesis Discussion, not
  treated in any chapter.**

---

### Critique decisions baked in (for the record)
1. **1.5.1 attribution gap** added — the intro previously set up Chapters 2/3 but not Chapter 1.
2. **1.2.2 reframed** from "LAI is the key metric" to the open *quantity vs vertical arrangement* question;
   the LAD profile is introduced in 1.2.
3. **1.1.4** in-situ observation backbone added (validation + spatial-limitation motivation).
4. **1.3.3** plot → landscape scaling added (red thread to Chapter 3).
5. **GEDI removed from the outline** — no false symmetry with the two sensors actually used; kept as a
   Discussion-only perspective.
6. Minor: "ecosystem engineer" → "physical modifier of the sub-canopy climate".
