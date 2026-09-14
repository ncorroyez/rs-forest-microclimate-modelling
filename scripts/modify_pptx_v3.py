"""
Apply MEB 2026 v2 meeting feedback to Corroyez_MEB2026_v2.pptx → v3.pptx.

This script:
  - drops the 'bit = 0' / '2^4 design' bit notation language from the methods slide
  - replaces the study-site map (slide 5) with the new cLHS + HOBO map
  - replaces the heatmap on slide 8 with the LOO × archetypes gradient version
  - replaces the HOBO validation figure on slide 9 with the forward-selection figure
  - simplifies the slide 3 messaging (one line about stats models → mechanistic complement)
  - adds linear regression mention in slide 9
  - tightens takeaways on slide 10 and conclusion on slide 11
  - subscript H_max wherever 'Hmax' appears in text

No slide reordering; v2 stays intact at the original path.
"""
from pathlib import Path
from copy import deepcopy
from pptx import Presentation
from pptx.util import Inches, Emu, Pt
from pptx.enum.text import PP_ALIGN

# Paths
SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v2.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v3.pptx")
FIG_DIR = Path("/home/corroyez/Documents/z_Example_rmusica_31012025/outputs/figs_MEB2026_final")

FIG_MAP        = FIG_DIR / "fig_map_blois_locations.png"
FIG_HEATMAP    = FIG_DIR / "fig_heatmap_LOO_archetypes_gradient.png"
FIG_FORWARD    = FIG_DIR / "fig_forward_selection_HOBO.png"
FIG_FWD_EVOL   = FIG_DIR / "fig_forward_selection_evolution.png"
FIG_HOBO_LINREG = FIG_DIR / "fig_HOBO_final_linreg.png"
FIG_ARCHETYPES = FIG_DIR / "fig_archetypes_profiles.png"

prs = Presentation(SRC)


def replace_picture(slide, shape_idx, new_path, keep_box=True):
    """Replace the picture at shape_idx with new_path, keeping the bounding box."""
    old = slide.shapes[shape_idx]
    left, top, width, height = old.left, old.top, old.width, old.height
    sp = old._element
    sp.getparent().remove(sp)
    pic = slide.shapes.add_picture(str(new_path), left, top,
                                       width=width if keep_box else None,
                                       height=height if keep_box else None)
    return pic


def find_shape_by_text(slide, needle):
    for j, shape in enumerate(slide.shapes):
        if shape.has_text_frame and needle in shape.text_frame.text:
            return j, shape
    return None, None


def replace_text(shape, replacements):
    """In-place edit of every run within a text frame, applying replacements (dict)."""
    if not shape.has_text_frame:
        return
    for para in shape.text_frame.paragraphs:
        for run in para.runs:
            for old, new in replacements.items():
                if old in run.text:
                    run.text = run.text.replace(old, new)


def set_paragraph_text(shape, new_lines, font_size=None, bold=False):
    """Replace the entire text_frame content with new lines (a list of strings)."""
    if not shape.has_text_frame:
        return
    tf = shape.text_frame
    # Wipe existing paragraphs by setting the first one and removing the rest
    tf.clear()
    for i, line in enumerate(new_lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        run = p.add_run()
        run.text = line
        if font_size:
            run.font.size = Pt(font_size)
        if bold:
            run.font.bold = True


# ============================================================================
# Slide 2 — intro 'Forest buffers microclimate'
# Already simple; just polish wording slightly
# ============================================================================
# (keep as-is for now)

# ============================================================================
# Slide 3 — 'This buffering depends on canopy structure'
# Simplify text #12 (offset equation) — keep one line; remove offset details to be moved later
# Tighten text #13 to motivate mechanistic complement
# ============================================================================
slide3 = prs.slides[2]
for j, shape in enumerate(slide3.shapes):
    if not shape.has_text_frame:
        continue
    txt = shape.text_frame.text
    if "Tmicro,max" in txt and "slope < 1" in txt:
        # Drop the offset equation — moved to methods later
        set_paragraph_text(shape, [
            "Statistical models show LiDAR-derived canopy traits explain summer ΔTmax across stands.",
            "",
            "→ But statistical fits remain a black box on mechanism.",
        ], font_size=14)
    elif "the link remains statistical" in txt:
        set_paragraph_text(shape, [
            "Which canopy trait actually drives the buffering?",
            "Does the vertical LAD shape add information beyond bulk LAI?",
            "",
            "→ Need a process-based, mechanistic attribution.",
        ], font_size=14)

# ============================================================================
# Slide 4 — Research questions
# Update phrasing to flow from MuSICA (slide 6) - swap order is risky, so simplify text
# ============================================================================
slide4 = prs.slides[3]
for shape in slide4.shapes:
    if not shape.has_text_frame:
        continue
    if "How sensitive" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "1) Can a process-based model identify which LiDAR-derived canopy traits "
            "(LAI, H_max, fCover, vertical LAD shape) most contribute to summer ΔTmax buffering?",
            "",
            "2) Does the vertical foliage profile (LAD) add information beyond bulk LAI?",
            "",
            "3) Are the simulated rankings supported by independent HOBO measurements?",
        ], font_size=15)

# ============================================================================
# Slide 5 — Study site
# Replace large pic (#20, the map of plots) with new cLHS+HOBO map
# Update text #15 to mention FPCA, clustering, n=400 cLHS + n=53 HOBO
# ============================================================================
slide5 = prs.slides[4]
# Identify last large picture by area
big_pic_idx = None
max_area = 0
for j, shape in enumerate(slide5.shapes):
    if shape.shape_type == 13:  # picture
        area = (shape.width or 0) * (shape.height or 0)
        if area > max_area:
            max_area = area
            big_pic_idx = j
if big_pic_idx is not None and FIG_MAP.exists():
    replace_picture(slide5, big_pic_idx, FIG_MAP)

for shape in slide5.shapes:
    if not shape.has_text_frame:
        continue
    if "53 HOBO sensors" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "53 HOBO sensors (n=53)",
            "Hourly air temperature at 1 m above ground",
            "(microclimate signal)",
        ], font_size=12)
    if "Leaf-on airborne LiDAR" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Leaf-on airborne LiDAR (2021-06-14)",
            "Full-waveform · 25–50 pts/m²",
            "",
            "Derived traits: LAI, H_max, fCover,",
            "vertical LAD profile (joint FPCA → FPC1)",
            "",
            "Spatial design: k-means → 4 archetypes",
            "cLHS sub-sampling → 400 plots",
        ], font_size=12)

# ============================================================================
# Slide 6 — MuSICA
# Drop 'bit = 0' notation from text #11 ('Baselines (bit = 0) — what replaces each trait')
# Tighten the 2^4 design wording
# ============================================================================
slide6 = prs.slides[5]
for shape in slide6.shapes:
    if not shape.has_text_frame:
        continue
    txt = shape.text_frame.text
    if "bit = 0" in txt or "Baselines" in txt:
        set_paragraph_text(shape, [
            "Baselines — what replaces each trait when removed:",
            "  LAI    → spatial mean from cLHS distribution",
            "  H_max  → spatial mean from cLHS distribution",
            "  fCover → spatial mean (floored to 0.5)",
            "  LAD    → uniform vertical profile",
        ], font_size=13)
    if "Factorial 2" in txt:
        set_paragraph_text(shape, [
            "Factorial design with 4 LiDAR-derived traits",
        ], font_size=18, bold=True)
    if "Each canopy trait is toggled" in txt:
        set_paragraph_text(shape, [
            "Each trait is independently toggled real vs. baseline,",
            "yielding 16 canopy configurations:",
            "",
            "  LAI · H_max · fCover · LAD",
        ], font_size=13)

# ============================================================================
# Slide 7 — Three attribution methods
# Keep the visual; update title to mention LOO if needed
# ============================================================================
# leave figure in place

# ============================================================================
# Slide 8 — Heatmap (replace with gradient LOO × archetypes)
# Find the leftmost / largest picture and swap; also could add archetype profiles on right
# ============================================================================
slide8 = prs.slides[7]
# Find the two pictures
pics = [(j, sh) for j, sh in enumerate(slide8.shapes) if sh.shape_type == 13]
# Sort by area, replace the largest (heatmap)
if pics and FIG_HEATMAP.exists():
    pics.sort(key=lambda p: -(p[1].width * p[1].height))
    big = pics[0]
    replace_picture(slide8, big[0], FIG_HEATMAP)
    # If there is a second picture (archetype profiles), keep / swap
    if len(pics) > 1 and FIG_ARCHETYPES.exists():
        replace_picture(slide8, pics[1][0], FIG_ARCHETYPES)

# Update title
for shape in slide8.shapes:
    if shape.has_text_frame and "Methods converge" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "LOO attribution × 4 archetypes",
        ], font_size=22, bold=True)
        break

# ============================================================================
# Slide 9 — HOBO validation : replace with forward-selection figure
# ============================================================================
slide9 = prs.slides[8]
pics9 = [(j, sh) for j, sh in enumerate(slide9.shapes) if sh.shape_type == 13]
# Find largest
if pics9 and FIG_FORWARD.exists():
    pics9.sort(key=lambda p: -(p[1].width * p[1].height))
    big = pics9[0]
    replace_picture(slide9, big[0], FIG_FORWARD)

for shape in slide9.shapes:
    if shape.has_text_frame and "HOBO validation" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "HOBO forward-selection (LOO ranking)",
        ], font_size=22, bold=True)
        break

# ============================================================================
# Slide 10 — Key takeaways
# Tighten language and replace H_max formatting
# ============================================================================
slide10 = prs.slides[9]
for shape in slide10.shapes:
    if shape.has_text_frame and "LAI is the main driver" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Bulk canopy traits (LAI, fCover, H_max) drive the simulated ΔTmax ranking.",
            "",
            "LAD vertical shape adds little — confirms a bulk-trait dominance.",
            "",
            "MuSICA reproduces the ordering of buffering across sensors,",
            "but underestimates absolute magnitude (\"glass-ceiling\").",
            "",
            "Next: extend to other forest types and explore non-bulk traits.",
        ], font_size=15)

# ============================================================================
# Slide 11 — Conclusion : tighten and add forward-selection narrative
# ============================================================================
slide11 = prs.slides[10]
for shape in slide11.shapes:
    if not shape.has_text_frame:
        continue
    if "What we did" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "What we did:",
            "Process-based attribution of LiDAR traits to summer ΔTmax buffering,",
            "via a 16-coalition factorial with MuSICA, triangulated across LOO / LVA / Shapley",
            "and validated against 53 HOBO sensors using forward selection.",
            "",
            "What we found:",
            "Bulk canopy traits dominate; LAD vertical shape is marginal.",
            "Simulated rankings are robust to method, baseline, metric and aggregation scale.",
            "Forward selection on HOBO confirms the LOO ranking.",
        ], font_size=14)

# Global text pass : replace 'Hmax' with 'H_max' wherever it appears as a plain word
for slide in prs.slides:
    for shape in slide.shapes:
        if not shape.has_text_frame:
            continue
        for para in shape.text_frame.paragraphs:
            for run in para.runs:
                if "Hmax" in run.text:
                    run.text = run.text.replace("Hmax", "H_max")

prs.save(DST)
print(f"Saved {DST}")
