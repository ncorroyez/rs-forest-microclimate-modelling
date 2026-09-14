"""
Apply post-meeting feedback to Corroyez_MEB2026_v3.pptx → v4.pptx.

Changes vs v3:
  1. Swap slides 4 (Research questions) and 6 (MuSICA model) — MuSICA before QR
  2. Add Methods 2-step slide (theory + validation), replacing 3-methods illustration
  3. Add baselines distribution figure to MuSICA slide
  4. Rename "Archetype" / "Archetypes" → "Profile" / "Profiles" globally
  5. Replace "bulk" / "Bulk" → "aggregate" / "Aggregate"
  6. Add equation offset block on the new methods slide
  7. Add title speaker notes (auxiliary data integration, model efficiency)
  8. Add structural profiles visual to methods slide
"""
from pathlib import Path
from copy import deepcopy
from pptx import Presentation
from pptx.util import Inches, Emu, Pt
from lxml import etree

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v3.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v4.pptx")
FIG_DIR = Path("/home/corroyez/Documents/z_Example_rmusica_31012025/outputs/figs_MEB2026_final")

FIG_HEATMAP_PROFILES = FIG_DIR / "fig_heatmap_LOO_profiles_gradient.png"
FIG_PROFILES         = FIG_DIR / "fig_profiles_structural.png"
FIG_BASELINES        = FIG_DIR / "fig_baselines_distributions.png"
FIG_LOO_SCHEMA       = FIG_DIR / "fig_LOO_schema.png"
FIG_FORWARD          = FIG_DIR / "fig_forward_selection_HOBO.png"
FIG_MAP              = FIG_DIR / "fig_map_blois_locations.png"

prs = Presentation(SRC)


def replace_picture(slide, shape_idx, new_path, keep_box=True):
    old = slide.shapes[shape_idx]
    left, top, width, height = old.left, old.top, old.width, old.height
    sp = old._element
    sp.getparent().remove(sp)
    return slide.shapes.add_picture(str(new_path), left, top,
                                       width=width if keep_box else None,
                                       height=height if keep_box else None)


def set_paragraph_text(shape, new_lines, font_size=None, bold=False):
    if not shape.has_text_frame:
        return
    tf = shape.text_frame; tf.clear()
    for i, line in enumerate(new_lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        run = p.add_run(); run.text = line
        if font_size: run.font.size = Pt(font_size)
        if bold:      run.font.bold = True


def global_text_replace(prs, replacements):
    for slide in prs.slides:
        for shape in slide.shapes:
            if not shape.has_text_frame: continue
            for para in shape.text_frame.paragraphs:
                for run in para.runs:
                    for old, new in replacements.items():
                        if old in run.text:
                            run.text = run.text.replace(old, new)


# ============================================================================
# Step 0 : global text substitutions (Archetype→Profile, bulk→aggregate)
# ============================================================================
global_text_replace(prs, {
    "Archetypes": "Profiles",
    "archetypes": "profiles",
    "Archetype": "Profile",
    "archetype":  "profile",
    "Bulk":       "Aggregate",
    "bulk":       "aggregate",
})

# ============================================================================
# Step 1 : swap slides 4 (QR) and 6 (MuSICA)
# Manipulate sldIdLst : slides 4 and 6 become 6 and 4 (positions 3 and 5 in 0-idx)
# ============================================================================
sldIdLst = prs.slides._sldIdLst
slide_ids = list(sldIdLst)
# Python list indexing : positions 3 (slide 4 = QR) and 5 (slide 6 = MuSICA)
slide_ids[3], slide_ids[5] = slide_ids[5], slide_ids[3]
# Reinsert in new order
for s in list(sldIdLst): sldIdLst.remove(s)
for s in slide_ids:      sldIdLst.append(s)

# After reordering :
#   v4 slide 4 = previously slide 6 (MuSICA)
#   v4 slide 5 = previously slide 5 (Study site) — unchanged
#   v4 slide 6 = previously slide 4 (QR)

# Re-read in new order
prs_slides = list(prs.slides)

# ============================================================================
# Step 2 : update page numbers on text boxes (originally "4", "5", "6")
# Find page-number text frames and renumber
# ============================================================================
for i, slide in enumerate(prs_slides, 1):
    for shape in slide.shapes:
        if not shape.has_text_frame: continue
        # heuristic : small text box at bottom right with just a number
        txt = shape.text_frame.text.strip()
        if txt.isdigit() and len(txt) <= 2:
            # Test if positioned near bottom-right (within 1.5" of slide edge)
            try:
                if (shape.left > Inches(8) and shape.top > Inches(4.5)):
                    for para in shape.text_frame.paragraphs:
                        for run in para.runs:
                            run.text = str(i)
            except:
                pass

# ============================================================================
# Step 3 : MuSICA slide (new slide 4) — add baselines figure
# ============================================================================
slide_musica = prs_slides[3]   # new slide 4
# Check if baselines figure fits ; insert at bottom
if FIG_BASELINES.exists():
    # Find an empty area at bottom — add picture
    pic = slide_musica.shapes.add_picture(
        str(FIG_BASELINES),
        left=Inches(0.3), top=Inches(3.4),
        width=Inches(9.4), height=Inches(1.6))

# Update baselines text to be more concise + mention rasters/histograms
for shape in slide_musica.shapes:
    if not shape.has_text_frame: continue
    txt = shape.text_frame.text
    if "Baselines" in txt and "replaces each trait" in txt:
        set_paragraph_text(shape, [
            "Baselines — what replaces each trait when removed:",
            "  LAI / H_max → cLHS spatial mean",
            "  fCover      → cLHS mean (floored to 0.5)",
            "  LAD         → uniform vertical profile",
        ], font_size=11)

# ============================================================================
# Step 4 : Methods slide (new slide 7) — replace 3-methods illustration
# with LOO 2-step explanation + LOO schema figure
# ============================================================================
slide_methods = prs_slides[6]   # new slide 7

# Update title
for shape in slide_methods.shapes:
    if shape.has_text_frame and "Three attribution methods" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Methods: from theory (LOO) to validation (HOBO)",
        ], font_size=22, bold=True)
        break
    if shape.has_text_frame and "attribution methods" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Methods: from theory (LOO) to validation (HOBO)",
        ], font_size=22, bold=True)
        break

# Replace large picture with LOO schema
pics = [(j, sh) for j, sh in enumerate(slide_methods.shapes) if sh.shape_type == 13]
if pics and FIG_LOO_SCHEMA.exists():
    pics.sort(key=lambda p: -(p[1].width * p[1].height))
    # Big illustration → LOO schema (resized to fit)
    big_j = pics[0][0]
    old = slide_methods.shapes[big_j]
    left = Inches(0.3); top = Inches(0.9)
    width = Inches(5.8); height = Inches(3.3)
    sp = old._element; sp.getparent().remove(sp)
    slide_methods.shapes.add_picture(str(FIG_LOO_SCHEMA), left, top, width, height)

# Add 2-step text explanation on the right side
right_textbox = slide_methods.shapes.add_textbox(
    Inches(6.3), Inches(0.9), Inches(3.5), Inches(3.5))
tf = right_textbox.text_frame; tf.clear()
p1 = tf.paragraphs[0]; p1.add_run().text = "Two-step procedure:"
p1.runs[0].font.bold = True; p1.runs[0].font.size = Pt(14)
p2 = tf.add_paragraph(); p2.add_run().text = ""
p3 = tf.add_paragraph(); r3 = p3.add_run()
r3.text = "Step 1 — Attribution (theory)"; r3.font.bold = True; r3.font.size = Pt(13)
p4 = tf.add_paragraph(); r4 = p4.add_run()
r4.text = "Remove each trait one at a time from MuSICA. Δ_v = ΔTmax due to that trait."
r4.font.size = Pt(11)
p5 = tf.add_paragraph(); p5.add_run().text = ""
p6 = tf.add_paragraph(); r6 = p6.add_run()
r6.text = "Step 2 — Validation (data)"; r6.font.bold = True; r6.font.size = Pt(13)
p7 = tf.add_paragraph(); r7 = p7.add_run()
r7.text = ("Forward-selection: add traits in the LOO order and check the "
           "Spearman + linear-regression fit against 53 HOBO sensors.")
r7.font.size = Pt(11)
p8 = tf.add_paragraph(); p8.add_run().text = ""
p9 = tf.add_paragraph(); r9 = p9.add_run()
r9.text = "Offset model:"; r9.font.bold = True; r9.font.size = Pt(12)
p10 = tf.add_paragraph(); r10 = p10.add_run()
r10.text = "Tmicro,max ≈ slope × Tmacro,max + offset"; r10.font.size = Pt(11)
p11 = tf.add_paragraph(); r11 = p11.add_run()
r11.text = "ΔTmax = Tmicro,max − Tmacro,max < 0 → canopy buffers"
r11.font.size = Pt(11); r11.font.italic = True

# ============================================================================
# Step 5 : Heatmap slide (slide 8) — swap heatmap to LOO×profiles + profiles
# ============================================================================
slide_results = prs_slides[7]
pics = [(j, sh) for j, sh in enumerate(slide_results.shapes) if sh.shape_type == 13]
pics.sort(key=lambda p: -(p[1].width * p[1].height))

if len(pics) >= 1 and FIG_HEATMAP_PROFILES.exists():
    old = slide_results.shapes[pics[0][0]]
    left, top = old.left, old.top
    sp = old._element; sp.getparent().remove(sp)
    slide_results.shapes.add_picture(
        str(FIG_HEATMAP_PROFILES), left, top,
        width=Inches(5.5), height=Inches(2.3))

if len(pics) >= 2 and FIG_PROFILES.exists():
    old = slide_results.shapes[pics[1][0]]
    left, top = old.left, old.top
    sp = old._element; sp.getparent().remove(sp)
    slide_results.shapes.add_picture(
        str(FIG_PROFILES), left, top,
        width=Inches(4.0), height=Inches(4.0))

for shape in slide_results.shapes:
    if shape.has_text_frame and ("LOO attribution" in shape.text_frame.text
                                    or "Profiles" in shape.text_frame.text
                                    or "Profile" in shape.text_frame.text):
        if "×" in shape.text_frame.text or "attribution" in shape.text_frame.text:
            set_paragraph_text(shape, [
                "LOO attribution × 4 structural profiles",
            ], font_size=22, bold=True)
            break

# ============================================================================
# Step 6 : HOBO forward-selection slide (slide 9) — already done in v3, no-op
# ============================================================================

# ============================================================================
# Step 7 : Update takeaways (slide 10) — use 'aggregate' wording
# ============================================================================
slide_takeaway = prs_slides[9]
for shape in slide_takeaway.shapes:
    if shape.has_text_frame and "Aggregate canopy traits" in shape.text_frame.text:
        set_paragraph_text(shape, [
            "Aggregate canopy traits (LAI, fCover, H_max) drive the simulated ΔTmax ranking.",
            "",
            "LAD vertical shape adds little — confirms aggregate-trait dominance.",
            "",
            "MuSICA reproduces the ordering of buffering across sensors;",
            "absolute magnitude is underestimated (\"glass-ceiling\").",
            "",
            "Forward-selection on HOBO: H_max addition flips the ranking from",
            "anti- to positive correlation → critical for predicting plot ordering.",
        ], font_size=13)

# ============================================================================
# Step 8 : Conclusion (slide 11) — already touched in v3, ensure phrasing OK
# ============================================================================
slide_concl = prs_slides[10]
for shape in slide_concl.shapes:
    if shape.has_text_frame and ("What we did" in shape.text_frame.text):
        set_paragraph_text(shape, [
            "What we did:",
            "Process-based attribution of LiDAR-derived canopy traits to summer ΔTmax buffering,",
            "via a 16-coalition factorial with MuSICA and Leave-One-Out (LOO),",
            "validated against 53 HOBO sensors with a forward-selection scheme.",
            "",
            "What we found:",
            "Aggregate canopy traits dominate; LAD vertical shape is marginal.",
            "Simulated rankings are robust to baseline, metric, and aggregation scale.",
            "H_max addition flips the HOBO ordering — a feature beyond MAE-based rank.",
        ], font_size=13)

# ============================================================================
# Step 9 : speaker notes on title slide — auxiliary data integration tagline
# ============================================================================
slide_title = prs_slides[0]
notes_tf = slide_title.notes_slide.notes_text_frame
notes_tf.clear()
notes_tf.paragraphs[0].add_run().text = (
    "Title alt for press / tagline: \"Improving auxiliary data integration "
    "and model efficiency in process-based microclimate simulation.\""
)

# Notes for slide 4 (MuSICA) — explain factorial design and baselines
notes_musica = prs_slides[3].notes_slide.notes_text_frame
notes_musica.clear()
notes_musica.paragraphs[0].add_run().text = (
    "MuSICA = multilayer mechanistic model (Ogée et al., 2003). We feed it ERA5 "
    "macroclimate forcing plus LiDAR canopy structure. We toggle each of the 4 "
    "LiDAR-derived traits ON (real) or OFF (replaced by the baseline distribution "
    "mean — see histograms below), yielding a 2^4 = 16 factorial design."
)

# Notes for slide 7 (methods) — explain offset equation context
notes_methods = prs_slides[6].notes_slide.notes_text_frame
notes_methods.clear()
notes_methods.paragraphs[0].add_run().text = (
    "Step 1 = theoretical attribution: each trait is removed and we compute "
    "Δ_v = T_max(LOO_v) - T_max(REF). Δ_v > 0 means v contributes to buffering. "
    "Step 2 = HOBO validation via forward selection: we add traits one at a time "
    "in the LOO ranking order and quantify the predictive skill (Spearman + linear "
    "regression on observed ΔTmax = Tmicro,max - Tmacro,max). The offset equation "
    "links micro and macro daily Tmax."
)

prs.save(DST)
print(f"Saved {DST}")
