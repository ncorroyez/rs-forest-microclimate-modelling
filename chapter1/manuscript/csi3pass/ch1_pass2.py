# -*- coding: utf-8 -*-
"""Second batch: split every sentence of 40 words or more, and the remaining
mannerisms. Numbers, citation keys and cross-references untouched."""
import io
R = "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/redaction/"
F = R + "manuscript_chap1_EN_native20_csi3pass_2026-09-02.md"
E = []
def e(fault, old, new): E.append((fault, old, new))

e("47 mots", "We coupled four ALS-derived variables, leaf area index, maximum height, fractional cover and the\nvertical leaf-area-density profile, to the multilayer model MuSICA, perturbing each around the real canopy of\n400 plots selected by conditioned Latin hypercube sampling in a temperate oak forest (Blois, France; summer 2021).",
  "We coupled four ALS-derived variables to the multilayer model MuSICA: leaf area index, maximum\nheight, fractional cover and the vertical leaf-area-density profile. We perturbed each around the real canopy\nof 400 plots, selected by conditioned Latin hypercube sampling in a temperate oak forest (Blois, France;\nsummer 2021).")

e("64 mots", "We hypothesize that (H1) vertical structure modulates radiative and turbulent transfer, so that the leaf-area-density (LAD)\nprofile is a key parameter for a faithful microclimate simulation; and (H2) at fixed leaf area, height and cover,\ntop-heavy profiles, with leaf area concentrated high in the canopy, attenuate the diurnal cycle more strongly,\nbecause radiation is intercepted higher and the canopy heat source is shifted upward.",
  "We hypothesize (H1) that vertical structure modulates radiative and turbulent transfer, so that the\nleaf-area-density (LAD) profile is a key parameter for a faithful microclimate simulation. We hypothesize (H2)\nthat at fixed leaf area, height and cover, top-heavy profiles attenuate the diurnal cycle more strongly, the leaf\narea being concentrated high in the canopy. Radiation is then intercepted higher and the canopy heat source is\nshifted upward.")

e("79 mots", "To test these\nhypotheses, we (1) derive a structural typology of the LAD profile, (2) perturb each variable around every plot's realistic canopy with MuSICA and measure the local sensitivity of ΔT~max~ within each structural type,\ncontrasting the real vertical profile against a uniform one, and (3) test the direction H2 asserts by sweeping the vertical balance of the profile from bottom- to top-heavy at fixed leaf area, height and cover, repeated across the leaf-area range the forest occupies.",
  "We test these hypotheses in\nthree steps. First, we derive a structural typology of the LAD profile. Second, we perturb each variable around\nevery plot's realistic canopy with MuSICA and measure the local sensitivity of ΔT~max~ within each structural\ntype, contrasting the real vertical profile against a uniform one. Third, we sweep the vertical balance of the\nprofile from bottom- to top-heavy at fixed leaf area, height and cover, repeated across the leaf-area range the\nforest occupies. That sweep tests the direction H2 asserts.")

e("42 mots", "LAI is the vertically integrated one-sided LAD, an effective plant area index: the MacArthur-Horn inversion assumes a spherical leaf-angle distribution with the conventional extinction coefficient *k* = 0.5 and no foliage-clumping correction, so the value is a proxy for true leaf area.",
  "LAI is the vertically integrated one-sided LAD, an effective plant area index. The MacArthur-Horn inversion assumes a spherical leaf-angle distribution with the conventional extinction coefficient *k* = 0.5 and no foliage-clumping correction, so the value is a proxy for true leaf area.")

e("43 mots", "All four variables and the LAD profile are computed natively on a 20 m grid, directly from the returns in each 20 m cell (MacArthur-Horn inversion; LAD binned at 1 m resolution from 1.5 m above ground), not aggregated from a finer product.",
  "All four variables and the LAD profile are computed natively on a 20 m grid, directly from the returns in each 20 m cell, and not aggregated from a finer product. The inversion is MacArthur-Horn, with LAD binned at 1 m resolution from 1.5 m above ground.")

e("45 mots", "To separate architecture from quantity, we apply a functional principal component analysis (FPCA) to the LAD profiles after double normalization (each profile rescaled to z/*H*~max~ in height and to unit integral in density), on a B-spline basis whose dimension is selected by generalized cross-validation [@ramsayFunctionalDataAnalysis1997].",
  "To separate architecture from quantity, we apply a functional principal component analysis (FPCA) to the LAD profiles, on a B-spline basis whose dimension is selected by generalized cross-validation [@ramsayFunctionalDataAnalysis1997]. Each profile is first double-normalized, rescaled to z/*H*~max~ in height and to unit integral in density.")

e("40 mots", "Within each archetype we select 100 plots (10,000 iterations, fixed seed, fractional cover floored at 0.5, the same floor the model applies to its cover input) by conditioned Latin hypercube sampling [cLHS, @minasnyConditionedLatinHypercube2006], implemented in the *clhs* R package [@roudierClhsConditionedLatin2012].",
  "Within each archetype we select 100 plots by conditioned Latin hypercube sampling [cLHS, @minasnyConditionedLatinHypercube2006], implemented in the *clhs* R package [@roudierClhsConditionedLatin2012]. The sampler runs 10,000 iterations at a fixed seed, with fractional cover floored at 0.5, the same floor the model applies to its cover input.")

e("41 mots", "MuSICA is a multilayer, multileaf biophysical canopy model solving the coupled radiative, energy and turbulent balance through a vertically discretized canopy [@ogeeMuSICACO2Water2003], with radiation attenuated layer by layer following canopy-extinction theory [@baldocchiSolarRadiationOak1984; @campbellExtinctionCoefficientsRadiation1986]; its full configuration is given in Appendix H.",
  "MuSICA is a multilayer, multileaf biophysical canopy model. It solves the coupled radiative, energy and turbulent balance through a vertically discretized canopy [@ogeeMuSICACO2Water2003]. Radiation is attenuated layer by layer following canopy-extinction theory [@baldocchiSolarRadiationOak1984; @campbellExtinctionCoefficientsRadiation1986]. Appendix H gives its full configuration.")

e("43 mots", "We map one to the other with a neutral logarithmic profile, applied offline per plot at each plot's canopy height, and use the corrected wind in all simulations reported here; the profile, its roughness parameters and its magnitude are given in Appendix A.",
  "We map one to the other with a neutral logarithmic profile, applied offline per plot at each plot's canopy height. The corrected wind is used in all simulations reported here. Appendix A gives the profile, its roughness parameters and its magnitude.")

e("44 mots", "For each plot we run a reference at its measured variables, then perturb one variable at a time by a fixed step in its own native unit, holding the other three real, and contrast the real LAD profile against a uniform one (Section 2.6).",
  "For each plot we run a reference at its measured variables. We then perturb one variable at a time by a fixed step in its own native unit, holding the other three real, and contrast the real LAD profile against a uniform one (Section 2.6).")

e("45 mots", "From the perturbed simulations (Section 2.4) we read each variable's effect on two response metrics, ΔT~max~ and the micro-macro buffering slope, each over the full summer and over the hottest 10% of days, the 13 days of the window with the highest station daily maximum.",
  "From the perturbed simulations (Section 2.4) we read each variable's effect on two response metrics, ΔT~max~ and the micro-macro buffering slope. Each is read over the full summer and over the hottest 10% of days, the 13 days of the window with the highest station daily maximum.")

e("53 mots", "One qualification applies to the most collinear pair, leaf area and fractional cover (*r* = 0.94 across the 400 cLHS plots): moving one while the other is held fixed steps off their joint distribution, so the model's split *between* leaf area and cover is a model-only attribution the field cannot confirm (Section 4.1).",
  "One qualification applies to the most collinear pair, leaf area and fractional cover (*r* = 0.94 across the 400 cLHS plots). Moving one while the other is held fixed steps off their joint distribution. The model's split *between* leaf area and cover is therefore a model-only attribution the field cannot confirm (Section 4.1).")

e("44 mots", "Stepped in their native units (Fig. 4), leaf area and canopy closure both act: in the densest archetype P4, adding 0.5 of leaf area cools the understory by 0.26 °C and adding 10 points of cover cools it by 0.27 °C, against 0.05 and 0.02 °C in the open P1.",
  "Stepped in their native units (Fig. 4), leaf area and canopy closure both act. In the densest archetype P4, adding 0.5 of leaf area cools the understory by 0.26 °C and adding 10 points of cover cools it by 0.27 °C. The same steps give 0.05 and 0.02 °C in the open P1.")

e("42 mots", "There the model holds ΔT~max~ near zero and its leaf-area lever warms rather than cools (+0.05 °C for a 0.5 step), whereas observed open stands at this site depart strongly from neutral, the most open ones amplifying (Corroyez et al., in preparation).",
  "There the model holds ΔT~max~ near zero, and its leaf-area lever warms rather than cools (+0.05 °C for a 0.5 step). Observed open stands at this site depart strongly from neutral, the most open ones amplifying (Corroyez et al., in preparation).")

e("51 mots", "Defining the clusters partly on the profile components FPC1 to FPC3 places most profile-shape variance between archetypes, leaving the profile nearly invariant within a type, and the measured centroids occupy a narrow band of the swept axis (Section 4.1), so the design gradient spans only a small part of the contrast.",
  "Defining the clusters partly on the profile components FPC1 to FPC3 places most profile-shape variance between archetypes. The profile is therefore nearly invariant within a type, and the measured centroids occupy a narrow band of the swept axis (Section 4.1). The design gradient spans only a small part of the contrast.")

e("41 mots", "Those two ranges are not commensurable: a 0.51 °C span is a complete bottom-to-top rearrangement\nof the profile, whereas 0.26 °C is a ±0.5 step in leaf area, so the profile is the larger number\nonly under the more generous definition.",
  "Those two ranges are not commensurable. A 0.51 °C span is a complete bottom-to-top rearrangement\nof the profile, whereas 0.26 °C is a ±0.5 step in leaf area. The profile is the larger number only\nunder the more generous definition.")

e("48 mots, maniéré: carries", "Its result rests on the model's\nphysical parameterization (Appendix H) and on the 400-plot cLHS design, and carries no\nindependent field validation: no understory temperature measurement enters this chapter, and\nthe coupled model is evaluated against sub-canopy measurements in the companion\nmicroclimate-forcing study (Corroyez et al., in preparation).",
  "Its result rests on the model's\nphysical parameterization (Appendix H) and on the 400-plot cLHS design, and has no independent\nfield validation. No understory temperature measurement enters this chapter. The coupled model\nis evaluated against sub-canopy measurements in the companion microclimate-forcing study\n(Corroyez et al., in preparation).")

e("41 mots", "Two model-internal checks\nsupport the configuration: the reported offsets are insensitive to the prescribed\nboundary-layer height across a fiftyfold span (Appendix H), and a leaf-area control of ±0.01\nshows the solver deterministic at the scale of the weakest levers (Appendix G).",
  "Two model-internal checks\nsupport the configuration. The reported offsets are insensitive to the prescribed boundary-layer\nheight across a fiftyfold span (Appendix H). A leaf-area control of ±0.01 shows the solver\ndeterministic at the scale of the weakest levers (Appendix G).")

# ── remaining mannerisms ─────────────────────────────────────────────────────
e("maniéré: carries", "The mechanistic priority carries a silvicultural corollary.",
  "The mechanistic priority has a silvicultural corollary.")
e("maniéré: carries", "so the attribution carries no design points in the regime",
  "so the attribution has no design points in the regime")
e("maniéré: carries", "whether Sentinel-2 carries that information or saturates",
  "whether Sentinel-2 resolves that information or saturates")
e("métaphore positionnelle", "A lever is read at the operating point where its archetype sits, so the ranking mixes",
  "A lever is read at the operating point of its archetype, so the ranking mixes")

s = io.open(F, encoding="utf-8").read()
miss = []
for fault, old, new in E:
    if old not in s: miss.append((fault, old[:65])); continue
    s = s.replace(old, new, 1)
io.open(F, "w", encoding="utf-8").write(s)
print(f"{len(E)-len(miss)}/{len(E)} corrections appliquées")
for f, o in miss: print("  NON TROUVE |", f, "|", o.replace("\n"," "))
