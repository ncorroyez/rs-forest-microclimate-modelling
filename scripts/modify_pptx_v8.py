"""
Corroyez_MEB2026_v7.pptx → v8.pptx.

Move the factorial-design text + baselines text + baselines histogram figure
from slide 4 (MuSICA model) to slide 8 (Methods). Rationale: the factorial
design is methodological, MuSICA slide stays purely about the tool.
"""
from pathlib import Path
from copy import deepcopy
from pptx import Presentation
from pptx.util import Pt, Inches

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v7.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v8.pptx")
FIG_BASELINES = Path("/home/corroyez/Documents/z_Example_rmusica_31012025/"
                      "outputs/figs_MEB2026_final/fig_baselines_distributions.png")

prs = Presentation(SRC)


def set_paragraph_text(shape, new_lines, font_size=None, bold=False, italic=False):
    if not shape.has_text_frame: return
    tf = shape.text_frame; tf.clear()
    for i, line in enumerate(new_lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        run = p.add_run(); run.text = line
        if font_size: run.font.size = Pt(font_size)
        if bold:      run.font.bold = True
        if italic:    run.font.italic = True


# ============================================================================
# Slide 4 (MuSICA) : remove factorial + baselines text — keep MuSICA intro only
# ============================================================================
slide4 = prs.slides[3]
to_remove = []
for shape in slide4.shapes:
    if not shape.has_text_frame: continue
    txt = shape.text_frame.text
    if ("Factorial design" in txt
        or "Each trait is independently toggled" in txt
        or ("Baselines" in txt and "replaces each trait" in txt)):
        to_remove.append(shape)
for shape in to_remove:
    sp = shape.element; sp.getparent().remove(sp)
print(f"Slide 4 : removed {len(to_remove)} text blocks (factorial / baselines)")

# ============================================================================
# Slide 8 (Methods) : restructure into a 2-column layout with
# factorial + baselines + baselines fig + LOO schema on the left,
# 2-step explanation + why-LOO on the right
# ============================================================================
slide8 = prs.slides[7]
# Find the existing 2-step text box and the why-LOO box, and the LOO schema
# Remove them — we'll rebuild
keep_shapes_idx = set()
for i, shape in enumerate(slide8.shapes):
    if shape.has_text_frame:
        txt = shape.text_frame.text
        if ("Methods:" in txt or "page" in txt.lower() or txt.strip().isdigit()
            or "Integrating LiDAR" in txt):
            keep_shapes_idx.add(i)
            continue
    elif shape.shape_type == 13:
        # Keep small logos (corner) — width < 2 inches; remove large schema
        if shape.width < Inches(2):
            keep_shapes_idx.add(i)
            continue
# Remove everything not in keep
for i, shape in enumerate(list(slide8.shapes)):
    if i not in keep_shapes_idx:
        sp = shape.element; sp.getparent().remove(sp)

# Now add the new content
# Left column : factorial + baselines text + baselines histogram fig
left_text = slide8.shapes.add_textbox(Inches(0.2), Inches(0.85),
                                          Inches(4.7), Inches(2.0))
tf = left_text.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run()
r.text = "Factorial design: 4 LiDAR traits × {real, baseline} = 16 configurations"
r.font.size = Pt(13); r.font.bold = True
p2 = tf.add_paragraph(); p2.add_run().text = ""
p3 = tf.add_paragraph()
r3 = p3.add_run(); r3.text = "Baselines — what replaces each trait when removed:"
r3.font.size = Pt(12); r3.font.bold = True
for line in [
    "  LAI / H_max → cLHS spatial mean",
    "  fCover      → cLHS mean (floored to 0.5)",
    "  LAD         → uniform vertical profile",
]:
    pp = tf.add_paragraph(); rr = pp.add_run(); rr.text = line; rr.font.size = Pt(11)

# Baselines fig left
slide8.shapes.add_picture(str(FIG_BASELINES),
                              left=Inches(0.2), top=Inches(2.9),
                              width=Inches(4.7), height=Inches(1.3))

# LOO schema in left bottom
from pathlib import Path as _P
FIG_LOO = _P("/home/corroyez/Documents/z_Example_rmusica_31012025/"
              "outputs/figs_MEB2026_final/fig_LOO_schema.png")
slide8.shapes.add_picture(str(FIG_LOO),
                              left=Inches(0.2), top=Inches(4.3),
                              width=Inches(4.7), height=Inches(0.95))

# Right column : 2-step + why-LOO
right_text = slide8.shapes.add_textbox(Inches(5.1), Inches(0.85),
                                           Inches(4.7), Inches(4.0))
tf = right_text.text_frame; tf.clear()
def add(line, size=11, bold=False, italic=False):
    p = tf.add_paragraph() if tf.paragraphs[0].text else tf.paragraphs[0]
    r = p.add_run(); r.text = line
    r.font.size = Pt(size); r.font.bold = bold; r.font.italic = italic

# First paragraph
p0 = tf.paragraphs[0]; r0 = p0.add_run()
r0.text = "Two-step procedure"; r0.font.size = Pt(13); r0.font.bold = True
tf.add_paragraph()  # blank
p1 = tf.add_paragraph(); r1 = p1.add_run()
r1.text = "Step 1 — Attribution (theory)"; r1.font.size = Pt(12); r1.font.bold = True
p2 = tf.add_paragraph(); r2 = p2.add_run()
r2.text = ("Remove each trait one at a time from MuSICA. "
           "Δ_v = T_max(LOO_v) − T_max(REF). Positive Δ_v ⇒ v buffers.")
r2.font.size = Pt(10)
tf.add_paragraph()
p3 = tf.add_paragraph(); r3 = p3.add_run()
r3.text = "Step 2 — Validation (data)"; r3.font.size = Pt(12); r3.font.bold = True
p4 = tf.add_paragraph(); r4 = p4.add_run()
r4.text = ("Forward-selection: add traits in the LOO order, check fit to "
           "53 HOBO sensors (Spearman + linear regression).")
r4.font.size = Pt(10)
tf.add_paragraph()
p5 = tf.add_paragraph(); r5 = p5.add_run()
r5.text = "Offset model:"; r5.font.size = Pt(11); r5.font.bold = True
p6 = tf.add_paragraph(); r6 = p6.add_run()
r6.text = "Tmicro,max ≈ slope × Tmacro,max + c"; r6.font.size = Pt(10)
p7 = tf.add_paragraph(); r7 = p7.add_run()
r7.text = "ΔTmax = Tmicro − Tmacro < 0 ⇒ canopy buffers"
r7.font.size = Pt(10); r7.font.italic = True

# Why-LOO box bottom right
why_box = slide8.shapes.add_textbox(Inches(5.1), Inches(4.3),
                                         Inches(4.7), Inches(0.95))
tf = why_box.text_frame; tf.clear()
p = tf.paragraphs[0]
r = p.add_run(); r.text = "Why LOO over LVA / Shapley?"
r.font.size = Pt(10); r.font.bold = True; r.font.italic = True
p2 = tf.add_paragraph()
r2 = p2.add_run()
r2.text = ("LVA isolates each trait alone (unrealistic). Shapley behaves "
           "similarly to LOO here. LOO = simplest, baseline-anchored.")
r2.font.size = Pt(9)

prs.save(DST)
print(f"Saved {DST}")

# Sanity check
print()
print("=== Slide 4 shapes (after) ===")
for shape in prs.slides[3].shapes:
    if shape.has_text_frame:
        t = shape.text_frame.text[:80].replace("\n", " | ")
        if t.strip(): print(f"  TXT: {t}")
    elif shape.shape_type == 13:
        print(f"  PIC: {shape.width/914400:.2f}x{shape.height/914400:.2f}in "
              f"@ ({shape.left/914400:.2f},{shape.top/914400:.2f})")
print()
print("=== Slide 8 shapes (after) ===")
for shape in prs.slides[7].shapes:
    if shape.has_text_frame:
        t = shape.text_frame.text[:80].replace("\n", " | ")
        if t.strip(): print(f"  TXT: {t}")
    elif shape.shape_type == 13:
        print(f"  PIC: {shape.width/914400:.2f}x{shape.height/914400:.2f}in "
              f"@ ({shape.left/914400:.2f},{shape.top/914400:.2f})")
