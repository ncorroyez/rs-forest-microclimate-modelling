"""
Corroyez_MEB2026_v9.pptx → v10.pptx.

Push the full set of new speaker notes (style of the v2 supervisor notes,
adapted to the v9 narrative : LOO only, profiles P1-P4, forward selection,
buffering convention ΔTmax = Tmicro − Tmacro).
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Pt

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v9.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v10.pptx")

prs = Presentation(SRC)


def set_notes(slide, text):
    nt = slide.notes_slide.notes_text_frame; nt.clear()
    for i, line in enumerate(text.split("\n")):
        p = nt.paragraphs[0] if i == 0 else nt.add_paragraph()
        r = p.add_run(); r.text = line; r.font.size = Pt(11)


NOTES = [
# Slide 1 — Title
"Hello, my name is Nathan Corroyez, I am a French PhD student in my final year "
"at INRAE here in Montpellier. My work focuses on using remote sensing data to "
"improve forest microclimate models.\n"
"Today I will present a talk on \"Integrating LiDAR-Derived Vegetation "
"Structural Attributes into the MuSICA Model to Map and Explain Fine-Scale "
"Forest Microclimate Variability\".",

# Slide 2 — Forest buffers microclimate
"Because of the canopy, the microclimate under the trees is buffered against "
"extremes outside.\n"
"This understory microclimate is a key component of many forest services. "
"It matters for microrefugia, regeneration, biogeochemical feedbacks.\n"
"And critically, understory biodiversity responds more strongly to this "
"microclimate than to the macroclimate. Forest management practices alter "
"canopy structure, affecting understory microclimate through mechanisms "
"not yet fully understood, potentially impacting forest resilience to "
"climate change.\n"
"In this work we focus on the understory microclimate at air temperature "
"at 1 m above the ground.",

# Slide 3 — Stats hint → mechanistic
"Understory microclimate is generally modelled with statistical approaches. "
"A widely used model is the slope-and-equilibrium approach, popularized by "
"the work of Eva Gril, who I assume is in the room!\n"
"She modelled the microclimate in a French temperate deciduous forest, "
"deriving three structural metrics from airborne LiDAR — maximum height, "
"Vertical Complexity Index, and Plant Area Index — and explained the "
"buffering slope statistically with these, reaching an R² of 91% against "
"field measurements.\n"
"Although widely used, this approach has two key limitations. First, it "
"assumes a linear relationship between microclimate and macroclimate, "
"which is not necessarily realistic because many parameters come into play. "
"Second, it does not take into account the laws of physics. So it tells us "
"that canopy structure relates to buffering, not which variable drives it "
"through which mechanism.\n"
"So we have statistical approaches that work, but we want a mechanistic "
"model that needs canopy structure as input — to see whether LiDAR-derived "
"structure alone is enough to predict the buffering, and ultimately to "
"spatialize the prediction.",

# Slide 4 — MuSICA model
"Now I will present the physical model used in this work, MuSICA — Multilayer "
"Simulator of the Interactions between a vegetation Canopy and the "
"Atmosphere. It mechanistically simulates exchanges of energy, water and "
"carbon between canopy and atmosphere.\n"
"The key advantage of this model for forest microclimate is that it divides "
"the canopy into multiple vertical layers, each with its own meteorological "
"state, which is ideal for vertical stratification of the microclimate.\n"
"It takes forest structure data as input — maximum height, vertical LAI "
"profile, fractional cover — plus soil data and meteorological forcing "
"above the canopy, here taken from ERA5 reanalysis. It provides hourly "
"predictions; for our case we look at air temperature at 1 m above the "
"ground.",

# Slide 5 — Study site
"This work was done on one study site: the Blois State Forest, a temperate "
"deciduous lowland forest in central France, almost entirely composed of "
"sessile oaks.\n"
"About the available data: a HOBO network of 53 micro-sensors deployed "
"inside the forest along a light gradient, recording hourly air temperature "
"at 1 m above the ground for 122 summer days, plus a weather station "
"outside the forest for macroclimate.\n"
"We also have a leaf-on airborne LiDAR acquisition from June 2021. The "
"point clouds were normalized then we computed a 1-m canopy height model. "
"From there we derived canopy variables at 25-m resolution: maximum height, "
"the vertical Leaf Area Density (LAD) profile via Beer-Lambert's law, the "
"total Leaf Area Index (LAI) integrated from the profile, fractional cover "
"(fCover), and we captured the dominant LAD vertical-shape mode with a "
"joint functional PCA (FPC1).\n"
"Using these LiDAR metrics, we ran a k-means clustering and drew 100 plots "
"per cluster using conditioned Latin Hypercube Sampling, so that the 400 "
"plots optimally cover the variable space at Blois. We see their location "
"on the map alongside the 53 HOBO sensors.",

# Slide 6 — Four structural profiles
"Here are the four typical structural profiles obtained from the k-means "
"clustering.\n"
"We can group them by structure: P1 is sparse, regenerating canopies with "
"low LAI and short trees. P2 is the dense, closed canopy with the highest "
"LAI and tall trees. P3 is mixed-height stands with intermediate density. "
"P4 is tall canopy with moderate LAI.\n"
"These four profiles will become the aggregation scale for the LOO "
"attribution we'll show in two slides — instead of looking plot by plot, "
"we summarize the attribution within each typical structure.",

# Slide 7 — Research questions
"Two research questions follow.\n"
"First, how sensitive is summer Tmax buffering to 3D canopy architecture? "
"Our hypothesis H1: aggregate variables — total leaf area, fractional "
"cover — dominate, while vertical structure plays a secondary role.\n"
"Then, does the vertical foliage profile, what we call LAD profile, add "
"information beyond total LAI? Hypothesis H2: LAD shape contributes "
"marginally once LAI, Hmax and fCover are accounted for.",

# Slide 8 — Methods 2-step (factorial + baselines + LOO + forward-selection + why-LOO)
"We use a two-step procedure.\n"
"The experimental design is a 2⁴ factorial: each of the four LiDAR-derived "
"traits is independently set either to its real plot value, or to a "
"baseline. Four traits gives 16 canopy configurations, which we simulate "
"with MuSICA for each plot, so 6400 simulations total.\n"
"For the baselines, shown bottom-left: LAI and Hmax are replaced by the "
"cLHS spatial mean — 3.1 and 21 m. fCover is replaced by the spatial mean, "
"floored at 0.5 — three of four profiles already have fCover above 0.98, "
"so the floor only changes sparse-canopy plots. And LAD is replaced by a "
"uniform vertical profile that preserves the plot's own LAI and Hmax.\n"
"Step 1 is the theoretical attribution with Leave-One-Out. Starting from "
"the full canopy with all four traits real, we remove each trait one at a "
"time and look at the resulting ΔTmax. Δv equals Tmax of LOO_v minus Tmax "
"of REF. A positive Δv means removing v warmed the canopy — so v "
"contributes to buffering.\n"
"Step 2 is the HOBO validation by forward selection. We start from the "
"baseline — all traits uniform — and add canopy traits one at a time in "
"the LOO order, then quantify how well the simulated ΔTmax matches what "
"each of our 53 HOBO sensors measured. We use Spearman ρ for rank "
"agreement and linear regression for magnitude.\n"
"A word on why we chose LOO and not LVA or Shapley: LVA tests each trait "
"alone vs. the null, which is unrealistic and misses interactions. "
"Shapley adds desirable axioms but behaves similarly to LOO here. LOO is "
"the simplest, most transparent, baseline-anchored attribution — and it "
"matches the forward-selection narrative.\n"
"About the convention: ΔTmax equals Tmicro minus Tmacro, so negative "
"means the canopy buffers.",

# Slide 9 — Heatmap LOO × profiles
"This heatmap shows the LOO attribution scores Δv for the four LiDAR-"
"derived traits across the four structural profiles. Red means the trait "
"contributes positively to buffering (removing it warms the canopy); blue "
"means the opposite. The colour intensity is the magnitude in degrees.\n"
"LAI clearly dominates the dense profile P2 — Δ_LAI is around +0.47 °C, "
"meaning removing LAI alone warms the canopy by half a degree on the "
"season-mean daily Tmax. In sparser or mixed-height profiles like P1 or "
"P3, Hmax and fCover gain weight. LAD vertical shape stays marginal "
"across all four profiles — never above 0.1 °C in magnitude.\n"
"So the take-home is: the attribution depends on canopy structure — there "
"is no single ranking that holds for every forest — but aggregate traits "
"dominate everywhere, supporting H1.",

# Slide 10 — HOBO forward selection
"Now we move from theory to data. For each of the 53 HOBO sensors we "
"computed the observed ΔTmax — Tmicro_obs minus Tmacro — and the simulated "
"ΔTmax from MuSICA at each forward-selection step. The convention is the "
"same on both axes: negative means buffering.\n"
"We start from the baseline — uniform canopy — in the first panel, then "
"we add traits one at a time in the LOO ranking order: fCover, LAI, Hmax, "
"and finally LAD which brings us to the full REF simulation.\n"
"Reading the panels: Spearman ρ improves from −0.2 with the baseline up "
"to +0.66 with the full structure restored. R² also improves but stays "
"modest — the simulated buffering spans about 2 °C while the observed "
"range is around 7 °C. This is the \"glass-ceiling\" effect — MuSICA "
"underestimates the absolute magnitude.\n"
"Note also: ρ is negative for the first two steps — only fCover, then "
"fCover plus LAI — and flips to positive only once Hmax is added. So "
"Hmax is the critical trait for recovering the correct plot ordering, "
"even though it ranks third in the LOO MAE. The MAE-rank and the "
"ordering-rank are not the same question.",

# Slide 11 — Takeaways
"A few key takeaways.\n"
"First, aggregate canopy traits — LAI, fCover, Hmax — drive the simulated "
"ΔTmax buffering. LAD vertical shape adds little, supporting H1 and H2.\n"
"Second, MuSICA reproduces the *ordering* of buffering across HOBO "
"sensors, but underestimates the absolute *magnitude*. A "
"\"glass-ceiling\" effect to investigate, possibly related to the "
"forcing bias.\n"
"Third, the LOO ranking is robust to method (LVA and Shapley converge), "
"baseline, metric (MAE, RMSE), and aggregation scale (plot vs profile). "
"LAI is the universal driver; Hmax is the order-flipping one on HOBO.\n"
"Future work: spatialize the buffering predictions across the whole "
"forest using the LiDAR rasters; focus on hottest summer days to see if "
"the ranking shifts under extremes; correct the macroclimate forcing to "
"close the glass-ceiling gap; and extend the approach to other temperate "
"forests.",

# Slide 12 — Conclusion
"To conclude. What we did: a mechanism-based attribution of LiDAR-derived "
"canopy variables to summer ΔTmax buffering at Blois forest, via a "
"16-coalition factorial with MuSICA, using Leave-One-Out attribution, "
"and validated on 53 HOBO sensors with a forward-selection scheme.\n"
"What we found: aggregate canopy variables drive understory buffering — "
"LAI as the universal driver. The vertical foliage profile is marginal "
"once total leaf area, height, and cover are accounted for. Dense "
"canopies buffer the macroclimate, sparse stands amplify it — consistent "
"with what statistical approaches found, but now with a mechanism-aware "
"explanation.\n"
"Why it matters: we can now predict and map fine-scale microclimate "
"using airborne LiDAR coupled with MuSICA. Key in a climate-change "
"context, for biodiversity conservation and forest regeneration. And for "
"MEB specifically: the factorial-attribution framework is transferable "
"to other microclimate-rich systems — tropical canopies, soil "
"microclimates, water bodies, urban canopies.\n"
"Thanks for your attention.",
]

for i, note in enumerate(NOTES):
    set_notes(prs.slides[i], note)

prs.save(DST)
print(f"Saved {DST}")
print()

# Sanity : word counts
prs2 = Presentation(DST)
total = 0
for i, slide in enumerate(prs2.slides, 1):
    wc = len(slide.notes_slide.notes_text_frame.text.split())
    total += wc
    print(f"  Slide {i:2d} : {wc} words")
print(f"  Total : {total} words")
print(f"  Estimated speech @ 150 wpm : {total/150:.1f} min")
print(f"  Estimated speech @ 130 wpm : {total/130:.1f} min")
