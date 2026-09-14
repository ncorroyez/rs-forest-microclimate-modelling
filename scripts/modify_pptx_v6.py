"""
Corroyez_MEB2026_v5.pptx → v6.pptx.

Pulls the speaker notes from Corroyez_MEB2026_notes.pptx (supervisor-feedback
version) and adapts them to the v5 narrative (12 slides, LOO only, forward-
selection HOBO, profiles slide inserted, 10-min target = ~60 s / slide).
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Pt

SRC   = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v5.pptx")
NOTES = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_notes.pptx")
DST   = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v6.pptx")

prs = Presentation(SRC)


def set_notes(slide, text):
    nt = slide.notes_slide.notes_text_frame
    nt.clear()
    lines = text.split("\n")
    for i, line in enumerate(lines):
        p = nt.paragraphs[0] if i == 0 else nt.add_paragraph()
        run = p.add_run(); run.text = line
        run.font.size = Pt(11)


# ============================================================================
# Slide 1 — Title
# ============================================================================
set_notes(prs.slides[0], (
"Hello, my name is Nathan Corroyez. I'm a French PhD student in my final year "
"at INRAE Montpellier, working on remote sensing to improve forest microclimate "
"models. Today I'll present \"Integrating LiDAR-Derived Vegetation Structural "
"Attributes into the MuSICA Model to Map and Explain Fine-Scale Forest "
"Microclimate Variability\".\n"
"\n"
"(Tagline alt for slides/abstract: improving auxiliary-data integration and "
"model efficiency in process-based microclimate simulation.)"
))

# ============================================================================
# Slide 2 — Forest buffers microclimate, and that matters
# ============================================================================
set_notes(prs.slides[1], (
"We focus on understory microclimate — air temperature at 1 m above the "
"ground. Because of the canopy, the microclimate under trees is buffered "
"against extremes outside. It matters for microrefugia, regeneration, "
"biogeochemical feedbacks, and understory biodiversity responds more strongly "
"to it than to macroclimate. Forest management alters canopy structure and "
"thus this buffering — through mechanisms not yet fully understood."
))

# ============================================================================
# Slide 3 — This buffering depends on canopy structure
# ============================================================================
set_notes(prs.slides[2], (
"Understory microclimate is usually modelled statistically — the slope-and-"
"equilibrium approach popularised by Eva Gril, who modelled buffering in a "
"French temperate deciduous forest using three LiDAR-derived metrics (height, "
"complexity, plant area index) and reached R² ≈ 91 % against in-situ data.\n"
"\n"
"Two limits: this is a linear, statistical link — and it doesn't tell us "
"which canopy trait drives buffering through which physical mechanism. "
"Hence the need for a process-based, mechanistic attribution."
))

# ============================================================================
# Slide 4 — MuSICA model
# ============================================================================
set_notes(prs.slides[3], (
"MuSICA — Multilayer Simulator of the Interactions between a Canopy and the "
"Atmosphere. A mechanistic model of energy/water/carbon exchange between "
"canopy and atmosphere, that divides the canopy into vertical layers (ideal "
"for microclimate stratification). Inputs: canopy structure (Hmax, LAI/LAD "
"profile, fCover), soil, and ERA5 macroclimate forcing. Output: hourly air "
"temperature at each layer, including 1 m above ground.\n"
"\n"
"Design: a factorial of 4 canopy traits each toggled real vs. baseline → 16 "
"canopy configurations. Baselines (histograms below): LAI / H_max are "
"replaced by the cLHS spatial mean, fCover by the same mean floored at 0.5, "
"LAD by a uniform vertical profile."
))

# ============================================================================
# Slide 5 — Study site
# ============================================================================
set_notes(prs.slides[4], (
"Blois state forest — a temperate deciduous lowland forest in central "
"France, almost entirely sessile oak. \n"
"\n"
"Two data streams. First, 53 HOBO micro-sensors deployed along a light "
"gradient, recording hourly air temperature at 1 m above the ground. "
"Second, a leaf-on full-waveform LiDAR acquisition from June 2021, "
"25–50 pts/m². We derived canopy traits at 25-m resolution: Hmax from CHM, "
"the LAD vertical profile via Beer–Lambert, total LAI integrated from the "
"profile (FPC1 from a joint FPCA captures the dominant vertical-shape mode), "
"and fCover. From these we ran a k-means clustering with balanced cLHS sampling "
"of 100 plots per cluster, yielding 400 plots and 4 typical structural profiles."
))

# ============================================================================
# Slide 6 — Four structural profiles (NEW slide)
# ============================================================================
set_notes(prs.slides[5], (
"Here are the four typical profiles obtained from the k-means clustering on "
"the LiDAR features. We can group them by structure: P1 = sparse, "
"regenerating canopies with low LAI and low Hmax. P2 = dense, closed "
"canopies with high LAI and tall trees. P3 = mixed-height stands with "
"intermediate density. P4 = tall canopy with moderate LAI.\n"
"\n"
"These profiles will be used as the aggregation scale for our attribution: "
"the heatmap later shows how each LiDAR trait contributes within each "
"profile."
))

# ============================================================================
# Slide 7 — Research questions
# ============================================================================
set_notes(prs.slides[6], (
"Two research questions follow.\n"
"\n"
"First, which LiDAR-derived canopy traits most contribute to summer Tmax "
"buffering? Our H1: aggregate traits (LAI, H_max, fCover) dominate; vertical "
"structure (LAD shape) plays a secondary role.\n"
"\n"
"Second, does the vertical LAD profile add information beyond bulk LAI? "
"H2: LAD shape contributes marginally once LAI, H_max, fCover are accounted "
"for.\n"
"\n"
"We test these by attribution within MuSICA and validate the ranking on the "
"HOBO sensors."
))

# ============================================================================
# Slide 8 — Methods: from theory (LOO) to validation (HOBO)
# ============================================================================
set_notes(prs.slides[7], (
"Two-step procedure.\n"
"\n"
"Step 1 — Theoretical attribution with Leave-One-Out. Starting from the full "
"canopy with all four traits real, we remove each trait one at a time (set "
"to its baseline) and look at the resulting ΔTmax. Δ_v = T_max(LOO_v) − "
"T_max(REF). A positive Δ_v means removing v warmed the canopy — so v "
"contributes to buffering.\n"
"\n"
"Step 2 — Validation by forward selection on HOBO. We start from the "
"baseline (all uniform) and add canopy traits one at a time, in the order "
"given by the LOO ranking, then quantify how well the simulated ΔTmax "
"matches what each of our 53 HOBO sensors measured. We use Spearman ρ for "
"rank agreement and linear regression for magnitude.\n"
"\n"
"Why LOO over LVA / Shapley? LVA tests each trait alone vs the null and "
"misses interactions. Shapley adds desirable axioms but behaves similarly "
"to LOO here. LOO is the simplest, most transparent, baseline-anchored "
"attribution — and matches the forward-selection narrative best."
))

# ============================================================================
# Slide 9 — Heatmap LOO × profiles
# ============================================================================
set_notes(prs.slides[8], (
"Heatmap of LOO attribution scores (Δ_v) for the four LiDAR traits across "
"the four structural profiles. Red = trait contributes positively to "
"buffering (removing it warms the canopy); blue = the opposite. The "
"colour intensity is the magnitude of the contribution in °C.\n"
"\n"
"Reading: LAI dominates in dense profiles (P2: +0.47 °C). In sparser or "
"mixed-height profiles (P1, P3), H_max and fCover gain weight. LAD vertical "
"shape stays marginal across all four profiles — never above 0.1 °C in "
"magnitude. So the attribution depends on canopy structure: there is no "
"single ranking that holds for all forests, but aggregate traits dominate "
"everywhere."
))

# ============================================================================
# Slide 10 — HOBO forward selection
# ============================================================================
set_notes(prs.slides[9], (
"For each of the 53 HOBO sensors we computed observed ΔTmax = "
"Tmicro_obs − Tmacro (negative when the canopy buffers, on the y-axis here) "
"and the simulated ΔTmax from MuSICA at each forward-selection step "
"(x-axis, same convention).\n"
"\n"
"We start from the baseline (uniform canopy, panel 1), then add traits in "
"the LOO order: fCover, LAI, H_max, then LAD = REF.\n"
"\n"
"Key reading: Spearman ρ improves from −0.2 with the baseline to +0.66 "
"once the full structure is back in MuSICA. R² also improves but stays "
"modest — the simulated buffering range is narrower than the observed one, "
"a 'glass-ceiling' that we discuss in the takeaways.\n"
"\n"
"Notable finding: ρ is negative for the first two steps (only fCover or "
"fCover+LAI), then flips to positive once H_max is added. So H_max is the "
"critical trait for getting the plot ordering right — even though it ranks "
"3rd in the LOO MAE."
))

# ============================================================================
# Slide 11 — Takeaways
# ============================================================================
set_notes(prs.slides[10], (
"Three takeaways.\n"
"\n"
"First, aggregate canopy traits (LAI, fCover, H_max) drive the simulated "
"ΔTmax buffering — LAD vertical shape is marginal, in line with our H1/H2.\n"
"\n"
"Second, MuSICA reproduces the *ordering* of buffering across HOBO sensors "
"but underestimates the absolute *magnitude* — the buffering range is "
"~7 °C in the data and ~2 °C in the simulations. A 'glass-ceiling' to "
"investigate in follow-up work, possibly forcing-bias related.\n"
"\n"
"Third, the LOO ranking is robust to method (LVA, Shapley converge), "
"baseline, metric (MAE, RMSE), and aggregation scale (cLHS plots vs. "
"profiles). LAI is the universal driver; H_max is the order-flipping one."
))

# ============================================================================
# Slide 12 — Conclusion
# ============================================================================
set_notes(prs.slides[11], (
"To conclude: we performed a mechanism-based attribution of LiDAR-derived "
"canopy traits to summer ΔTmax buffering at Blois forest, via a 16-coalition "
"factorial with MuSICA, using Leave-One-Out attribution, and validated on "
"53 HOBO sensors with a forward-selection scheme.\n"
"\n"
"Aggregate canopy traits dominate; LAD shape is marginal. Dense canopies "
"buffer the macroclimate, sparse stands amplify it — consistent with what "
"statistical models found, but now mechanism-aware. We highlight the "
"potential of incorporating airborne LiDAR-derived metrics into a "
"mechanistic model to map fine-scale microclimate fluctuations, which is "
"key in a climate-change context for biodiversity conservation and forest "
"management.\n"
"\n"
"Thanks for listening!"
))

prs.save(DST)
print(f"Saved {DST}")

# Quick sanity check — confirm notes are populated on every slide
print()
print("Notes word count per slide:")
prs2 = Presentation(DST)
for i, slide in enumerate(prs2.slides, 1):
    notes = slide.notes_slide.notes_text_frame.text.strip()
    wc = len(notes.split())
    print(f"  Slide {i:2d} : {wc} words")
total_words = sum(len(s.notes_slide.notes_text_frame.text.split()) for s in prs2.slides)
print(f"  Total : {total_words} words")
print(f"  Estimated speech time @ 150 wpm : {total_words / 150:.1f} min")
print(f"  Estimated speech time @ 130 wpm : {total_words / 130:.1f} min")
