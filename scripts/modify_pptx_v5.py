"""
Corroyez_MEB2026_v4.pptx → v5.pptx.

Three remaining items from the meeting notes:
  A) Add a "why LOO" 1-2 line block on the methods slide (now slide 7)
  B) Add a narrative annotation on the heatmap slide about how attribution
     scores evolve between homogeneous (P1) and dense (P2) profiles
  C) Insert a dedicated methods slide showing the 4 structural profiles,
     positioned between Study Site (slide 5) and Research Questions (slide 6)
     ; and slim down the heatmap slide so the heatmap is bigger
"""
from copy import deepcopy
from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from lxml import etree

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v4.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v5.pptx")
FIG_DIR = Path("/home/corroyez/Documents/z_Example_rmusica_31012025/outputs/figs_MEB2026_final")

FIG_PROFILES         = FIG_DIR / "fig_profiles_structural.png"
FIG_HEATMAP_PROFILES = FIG_DIR / "fig_heatmap_LOO_profiles_gradient.png"

prs = Presentation(SRC)


def set_paragraph_text(shape, new_lines, font_size=None, bold=False, italic=False):
    if not shape.has_text_frame:
        return
    tf = shape.text_frame; tf.clear()
    for i, line in enumerate(new_lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        run = p.add_run(); run.text = line
        if font_size: run.font.size = Pt(font_size)
        if bold:      run.font.bold = True
        if italic:    run.font.italic = True


def duplicate_slide(prs, src_idx):
    """Copy slide at src_idx and append it to the deck. Returns the new slide."""
    src_slide = prs.slides[src_idx]
    new_slide = prs.slides.add_slide(src_slide.slide_layout)
    # Copy all shapes from the source slide
    for shp in src_slide.shapes:
        el = shp.element
        new_el = deepcopy(el)
        new_slide.shapes._spTree.insert_element_before(new_el, 'p:extLst')
    return new_slide


def move_slide(prs, from_idx, to_idx):
    """Move slide from from_idx to to_idx in the slide order."""
    sld_id_lst = prs.slides._sldIdLst
    ids = list(sld_id_lst)
    moving = ids[from_idx]
    sld_id_lst.remove(moving)
    # Re-fetch after removal
    remaining = list(sld_id_lst)
    # Insert at to_idx in the (new) remaining list
    if to_idx >= len(remaining):
        sld_id_lst.append(moving)
    else:
        # Insert before the slide currently at to_idx
        sld_id_lst.insert(sld_id_lst.index(remaining[to_idx]), moving)


# ============================================================================
# A) "Why LOO" block on slide 7 (Methods)
# ============================================================================
slide_methods = prs.slides[6]   # slide 7 (0-idx = 6)
# Add a small block at the bottom-right summarising why LOO was chosen
why_box = slide_methods.shapes.add_textbox(
    Inches(6.3), Inches(4.4), Inches(3.5), Inches(0.9))
tf = why_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run(); r.text = "Why LOO over LVA / Shapley?"
r.font.bold = True; r.font.size = Pt(11); r.font.italic = True
p2 = tf.add_paragraph()
r2 = p2.add_run()
r2.text = ("LVA isolates each trait alone (unrealistic). Shapley adds axioms but "
           "behaves similarly to LOO here. LOO = simplest, transparent, baseline-"
           "anchored attribution.")
r2.font.size = Pt(9)

# ============================================================================
# B) Narrative on heatmap slide (slide 8)
# ============================================================================
slide_heatmap = prs.slides[7]   # slide 8
# Add a single-line takeaway under the heatmap area
note_box = slide_heatmap.shapes.add_textbox(
    Inches(0.0), Inches(4.3), Inches(10), Inches(0.8))
tf = note_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run()
r.text = ("From homogeneous (P1) → dense (P2) profiles, LAI dominates the "
          "attribution; H_max and fCover gain weight on mixed-height stands "
          "(P3, P4). LAD shape stays marginal across all profiles.")
r.font.size = Pt(12); r.font.italic = True

# ============================================================================
# C) Insert a dedicated profiles slide in the methods section
# Steps :
#   1. Duplicate slide 5 (study site) → new slide at end of deck
#   2. Clear its content and rebuild as "4 structural profiles"
#   3. Move it from end of deck to position 5 (between study site and QR)
#   4. Remove the profiles figure from slide 8 (heatmap+profiles) and enlarge heatmap
# ============================================================================
src_idx = 4   # slide 5 (study site) — for layout reuse
new_slide = duplicate_slide(prs, src_idx)
# Clear all shapes
for shp in list(new_slide.shapes):
    sp = shp.element
    sp.getparent().remove(sp)

# Title
title_box = new_slide.shapes.add_textbox(Inches(0), Inches(0), Inches(10), Inches(0.6))
tf = title_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run(); r.text = "Four structural profiles from k-means clustering"
r.font.bold = True; r.font.size = Pt(22)

# Profiles figure
new_slide.shapes.add_picture(
    str(FIG_PROFILES),
    left=Inches(0.2), top=Inches(0.8),
    width=Inches(9.6), height=Inches(3.7))

# Caption / explanation
cap_box = new_slide.shapes.add_textbox(Inches(0.2), Inches(4.6), Inches(9.6), Inches(0.7))
tf = cap_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run()
r.text = ("400 cLHS plots clustered into 4 representative LAD profiles (P1–P4). "
          "P2 = dense / closed canopy ; P1 = sparse / regenerating ; "
          "P3 = mixed-height ; P4 = tall canopy with moderate LAI.")
r.font.size = Pt(12); r.font.italic = True

# Footer (page number)
foot_box = new_slide.shapes.add_textbox(Inches(9.6), Inches(4.9), Inches(0.4), Inches(0.4))
tf = foot_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run(); r.text = "6"
r.font.size = Pt(11); r.font.italic = True

# Move new slide from end of deck (currently position 11 in 0-idx) to position 5
# Old order positions 0-10 (11 slides) + 1 new slide at the end (idx 11).
# Target : insert between slide 5 (idx 4 = study site) and slide 6 (idx 5 = QR)
move_slide(prs, from_idx=11, to_idx=5)

# Renumber page numbers (text boxes containing a small integer in bottom-right)
prs_slides = list(prs.slides)
for i, slide in enumerate(prs_slides, 1):
    for shape in slide.shapes:
        if not shape.has_text_frame: continue
        txt = shape.text_frame.text.strip()
        if txt.isdigit() and len(txt) <= 2:
            try:
                if (shape.left > Inches(8) and shape.top > Inches(4.5)):
                    for para in shape.text_frame.paragraphs:
                        for run in para.runs:
                            run.text = str(i)
            except:
                pass

# C-bis : remove profiles figure from heatmap slide (now slide 9, formerly 8)
# Find the heatmap slide by title
for slide in prs_slides:
    title_txt = ""
    for shape in slide.shapes:
        if shape.has_text_frame and "LOO attribution" in shape.text_frame.text:
            title_txt = shape.text_frame.text
            # Also retitle while we're here (profiles moved out)
            set_paragraph_text(shape, [
                "LOO attribution on the 4 structural profiles",
            ], font_size=22, bold=True)
            break
    if not title_txt:
        continue
    # This is the heatmap slide
    pics = list(slide.shapes)
    pics = [sh for sh in pics if sh.shape_type == 13]
    ratio = lambda sh: max(sh.width, sh.height) / min(sh.width, sh.height)
    def is_square_profile(sh):
        # square-ish AND area > 4 in² → the profiles facet plot
        return ratio(sh) < 1.4 and (sh.width * sh.height) > Inches(2) * Inches(2)
    def is_heatmap(sh):
        # wide, positioned in the body area (top > 0.5"), substantial size
        return (sh.width > sh.height
                  and sh.top > Inches(0.4)
                  and (sh.width * sh.height) > Inches(3) * Inches(1.5))
    # Step 1 : remove profile picture
    for sh in [p for p in pics if is_square_profile(p)]:
        sp = sh.element; sp.getparent().remove(sp)
    # Step 2 : resize the heatmap (after profile removal)
    pics = [sh for sh in slide.shapes if sh.shape_type == 13]
    heatmaps = [sh for sh in pics if is_heatmap(sh)]
    if heatmaps:
        hp = heatmaps[0]
        hp.left   = Inches(0.2)
        hp.top    = Inches(0.9)
        hp.width  = Inches(9.6)
        hp.height = Inches(3.3)
    break

prs.save(DST)
print(f"Saved {DST}")
print()
print("Slide order check :")
for i, slide in enumerate(prs.slides, 1):
    title_str = ""
    for shape in slide.shapes:
        if shape.has_text_frame:
            txt = shape.text_frame.text.strip()
            if txt and len(txt) < 120 and "Integrating LiDAR" not in txt and txt[:2] not in ("0","1","2","3","4","5","6","7","8","9"):
                title_str = txt.split("\n")[0][:80]
                break
    print(f"  Slide {i:2d} : {title_str}")
