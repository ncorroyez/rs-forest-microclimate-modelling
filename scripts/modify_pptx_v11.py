"""
Corroyez_MEB2026_v10.pptx → v11.pptx.

1. Slide 8 title → "Two-step framework: attribution + validation"
2. Rewrite every slide's speaker notes:
   - apply ALL meeting feedback (simplification of Eva Gril intro, no
     bit/2^4 notation, mechanistic transition + spatialisation,
     'aggregate' wording, 2-step methods, MAE-rank vs ordering, glass-
     ceiling, etc.)
   - target ~1300 words total → ≈10 min @ 130 wpm
   - simpler vocabulary (no 'axioms', no '2^4')
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Pt

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v10.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v11.pptx")

prs = Presentation(SRC)


def set_paragraph_text(shape, new_lines, font_size=None, bold=False):
    if not shape.has_text_frame: return
    tf = shape.text_frame; tf.clear()
    for i, line in enumerate(new_lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        run = p.add_run(); run.text = line
        if font_size: run.font.size = Pt(font_size)
        if bold:      run.font.bold = True


def set_notes(slide, text):
    nt = slide.notes_slide.notes_text_frame; nt.clear()
    for i, line in enumerate(text.split("\n")):
        p = nt.paragraphs[0] if i == 0 else nt.add_paragraph()
        r = p.add_run(); r.text = line; r.font.size = Pt(11)


# ============================================================================
# Slide 8 title change : F = "Two-step framework: attribution + validation"
# ============================================================================
for shape in prs.slides[7].shapes:
    if shape.has_text_frame and "Methods" in shape.text_frame.text \
       and "validation" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Two-step framework: attribution + validation",
        ], font_size=22, bold=True)
        break

# ============================================================================
# Rewrite all speaker notes — meeting-feedback compliant, ~1300 words total
# ============================================================================
NOTES = [
# Slide 1 — Title
"Hello, my name is Nathan Corroyez, French PhD student at INRAE Montpellier. "
"I work on remote sensing to improve forest microclimate models. Today I "
"present \"Integrating LiDAR-Derived Vegetation Structural Attributes into "
"the MuSICA Model to Map and Explain Fine-Scale Forest Microclimate "
"Variability\" — or more simply, improving auxiliary-data integration and "
"model efficiency in process-based microclimate simulation.",

# Slide 2 — Forest buffers microclimate
"Because of the canopy, the microclimate under trees is buffered against "
"extremes. This understory microclimate matters for microrefugia, "
"regeneration, biogeochemical feedbacks. Critically, understory "
"biodiversity responds more strongly to it than to macroclimate. Forest "
"management alters canopy structure and thus buffering, through mechanisms "
"not yet fully understood — a key uncertainty for forest resilience to "
"climate change. We focus here on air temperature at 1 m above the ground.",

# Slide 3 — Stats hint → mechanistic (compressed, meeting-feedback flow)
"Understory microclimate is usually modelled statistically — the slope-and-"
"equilibrium approach popularised by Eva Gril (probably in the room!): she "
"derived three LiDAR metrics and reached R² around 91 % on field data.\n"
"\n"
"Statistical models work, but they're linear and physics-free. We want a "
"mechanistic complement, where we feed LiDAR-derived canopy structure into "
"a process-based model, ultimately to spatialise the buffering prediction "
"across whole forests — and ask: is LiDAR-derived structure enough?",

# Slide 4 — MuSICA model
"MuSICA — Multilayer Simulator of the Interactions between a Canopy and the "
"Atmosphere. A mechanistic model of energy / water / carbon exchange, that "
"divides the canopy into vertical layers, ideal for microclimate "
"stratification. Inputs: canopy structure (Hmax, LAI / LAD profile, "
"fCover), soil, and ERA5 macroclimate forcing. Output: hourly air "
"temperature at every layer — including 1 m above the ground.",

# Slide 5 — Study site (compressed)
"Blois State Forest — a temperate deciduous lowland forest in central "
"France, almost entirely sessile oaks. Two data streams. 53 HOBO "
"micro-sensors deployed along a light gradient, recording hourly air "
"temperature at 1 m for 122 summer days, plus a weather station for "
"macroclimate.\n"
"\n"
"A leaf-on airborne LiDAR acquisition from June 2021. We derived canopy "
"variables at 25-m resolution: Hmax, vertical LAD profile via Beer-Lambert, "
"total LAI from the profile, fCover, plus the dominant LAD-shape mode from "
"a joint functional PCA (FPC1). We then k-means clustered the plots and "
"sampled 100 per cluster using conditioned Latin Hypercube — 400 plots "
"optimally covering the variable space, shown on the map alongside the "
"53 HOBO sensors.",

# Slide 6 — Four structural profiles
"Here are the four typical profiles from the k-means clustering. P1 = "
"sparse, regenerating; P2 = dense, closed canopy with the highest LAI and "
"tall trees; P3 = mixed-height stands with intermediate density; P4 = tall "
"canopy with moderate LAI. These four profiles become the aggregation "
"scale for the LOO attribution in two slides — instead of looking plot by "
"plot, we summarise the attribution within each typical structure.",

# Slide 7 — Research questions
"Two questions follow. H1: how sensitive is summer Tmax buffering to 3D "
"canopy architecture? We expect aggregate variables — total leaf area, "
"fractional cover, height — to dominate, vertical structure to play a "
"secondary role. H2: does the vertical foliage profile, the LAD shape, add "
"information beyond bulk LAI? We expect LAD shape to contribute marginally "
"once LAI, Hmax, fCover are accounted for.",

# Slide 8 — Two-step framework (compressed, no 2⁴, no jargon)
"Two-step framework.\n"
"\n"
"Step 1 — Attribution. We use a 16-configuration factorial design: each of "
"the four LiDAR traits is independently set to its real plot value or to a "
"baseline. Baselines, shown bottom-left: LAI and Hmax replaced by the cLHS "
"spatial mean; fCover by the same mean floored at 0.5; LAD by a uniform "
"vertical profile. From the full real canopy we then remove each trait one "
"at a time — Leave-One-Out. Δ_v = Tmax(LOO_v) − Tmax(REF). Positive Δ_v "
"means removing v warmed the canopy, so v contributes to buffering.\n"
"\n"
"Step 2 — HOBO validation. We start from the uniform baseline and add "
"traits one at a time in the LOO order, then quantify how well the "
"simulated ΔTmax matches the observed one at each step. Spearman ρ for "
"rank agreement, linear regression for magnitude.\n"
"\n"
"Why LOO and not LVA / Shapley? LVA tests each trait alone vs. the null — "
"unrealistic, misses interactions. Shapley uses all 16 configurations but "
"its ranking ends up close to LOO here. LOO is the simplest, most "
"transparent attribution — and matches the forward-selection narrative.\n"
"\n"
"Convention: ΔTmax = Tmicro − Tmacro, so negative means the canopy buffers.",

# Slide 9 — Heatmap LOO × profiles
"This heatmap shows the LOO attribution Δ_v of the four LiDAR traits "
"across the four structural profiles. Red = trait contributes positively "
"to buffering (removing it warms the canopy); blue = the opposite; "
"intensity is the magnitude in degrees.\n"
"\n"
"How does the score evolve from homogeneous to complex canopies? In the "
"dense P2, LAI dominates — Δ_LAI ≈ +0.47 °C. In sparser or mixed-height "
"profiles (P1, P3), Hmax and fCover gain weight. LAD shape stays marginal "
"everywhere — never above 0.1 °C. So the attribution depends on canopy "
"structure — no single ranking holds for every forest — but aggregate "
"traits dominate everywhere, supporting H1.",

# Slide 10 — HOBO forward selection
"From theory to data. For each of the 53 HOBO sensors we computed observed "
"ΔTmax and the simulated ΔTmax at each forward-selection step — same sign "
"convention on both axes: negative = canopy buffers. We start from the "
"uniform baseline, then add traits in the LOO order: fCover, LAI, Hmax, "
"finally LAD = REF.\n"
"\n"
"Reading: Spearman ρ improves from −0.2 with the baseline to +0.66 with "
"full structure restored. R² also improves but stays modest — simulated "
"buffering spans ~2 °C, observed range ~7 °C. This is the glass-ceiling: "
"MuSICA gets the ordering, not the magnitude.\n"
"\n"
"Notable: ρ is negative for the first two steps (fCover, fCover+LAI), then "
"flips to positive once Hmax is added. So Hmax is critical for recovering "
"the plot ordering — even though it only ranks third in the LOO MAE. The "
"MAE-rank and the ordering-rank are not the same question.",

# Slide 11 — Takeaways
"Three takeaways.\n"
"First, aggregate canopy traits — LAI, fCover, Hmax — drive the simulated "
"ΔTmax buffering. LAD vertical shape adds little, supporting H1 and H2.\n"
"Second, MuSICA reproduces the ordering of buffering across HOBO sensors "
"but underestimates the absolute magnitude — a glass-ceiling to "
"investigate, likely forcing-bias related.\n"
"Third, the LOO ranking is robust to method (LVA and Shapley converge), "
"baseline, metric and aggregation scale. LAI is the universal driver; "
"Hmax is the order-flipping one on HOBO.\n"
"Future work: spatialise the predictions across the whole forest using "
"LiDAR rasters; focus on hottest summer days for extremes; correct the "
"macroclimate forcing; extend the approach to other temperate forests.",

# Slide 12 — Conclusion
"What we did: a mechanism-based attribution of LiDAR-derived canopy "
"variables to summer ΔTmax buffering at Blois forest, via a 16-"
"configuration factorial with MuSICA, using Leave-One-Out, validated on "
"53 HOBO sensors with a forward-selection scheme.\n"
"What we found: aggregate canopy variables drive understory buffering — "
"LAI as the universal driver. The vertical foliage profile is marginal "
"once total leaf area, height and cover are accounted for. Dense canopies "
"buffer the macroclimate, sparse stands amplify it — consistent with "
"statistical approaches, now with a mechanism-aware explanation.\n"
"Why it matters: we can now predict and map fine-scale microclimate using "
"airborne LiDAR coupled with MuSICA — relevant in a climate-change "
"context for biodiversity conservation and forest regeneration. The "
"factorial-attribution framework is transferable to other microclimate-"
"rich systems — tropical canopies, soil microclimates, water bodies, "
"urban canopies.\n"
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
