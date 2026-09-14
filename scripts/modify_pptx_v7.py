"""
Corroyez_MEB2026_v6.pptx → v7.pptx.

Final touches :
  1. Remove the baselines histogram figure from slide 4 (MuSICA) — user will
     place it themselves elsewhere
  2. Add 'spatialisation' to the next-steps line on slide 11 (takeaways) + notes
  3. Sanity sweep : no residual 'bit = 0' or '2^4 bits' language anywhere
"""
from pathlib import Path
from pptx import Presentation
from pptx.util import Pt, Inches

SRC = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v6.pptx")
DST = Path("/home/corroyez/Downloads/1Papers/files/Corroyez_MEB2026_v7.pptx")

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
# 1) Remove the baselines histogram figure from slide 4 (MuSICA)
# ============================================================================
slide4 = prs.slides[3]
# Find the baselines picture : large wide picture in the bottom half of the slide
# It was inserted at (0.3, 3.4) with size 9.4 x 1.6 in v4
removed = False
for shape in list(slide4.shapes):
    if shape.shape_type != 13: continue
    # Identify the baselines fig : positioned at y > 3" and width > 8"
    if shape.top > Inches(3) and shape.width > Inches(7):
        sp = shape.element; sp.getparent().remove(sp)
        removed = True
        print("Removed baselines picture from slide 4")
        break
if not removed:
    print("No baselines picture found on slide 4 (already removed?)")

# ============================================================================
# 2) Add spatialisation to slide 11 takeaways (visible text + notes)
# ============================================================================
slide_takeaway = prs.slides[10]
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
            "Next: spatialisation of predictions across the forest, "
            "focus on hottest days, climate-forcing correction, extension to other forests.",
        ], font_size=13)
        break

# Also update notes for slide 11
set_notes(slide_takeaway, (
"Three takeaways.\n"
"\n"
"First, aggregate canopy traits (LAI, fCover, H_max) drive the simulated "
"ΔTmax buffering — LAD vertical shape is marginal, in line with our H1/H2.\n"
"\n"
"Second, MuSICA reproduces the *ordering* of buffering across HOBO sensors "
"but underestimates the absolute *magnitude* — the buffering range is "
"~7 °C in the data and ~2 °C in the simulations. A 'glass-ceiling' to "
"investigate.\n"
"\n"
"Third, the LOO ranking is robust to method, baseline, metric, and "
"aggregation scale.\n"
"\n"
"Next steps : spatialise the simulated buffering predictions across the "
"whole forest using the LiDAR rasters, focus on hottest summer days, "
"correct the macro-climate forcing to address the glass-ceiling bias, "
"and extend the approach to other temperate forests."
))

# ============================================================================
# 3) Sanity sweep : verify no 'bit', '2^4', 'binary' language
# (defensive — already clean in v6 per the previous check)
# ============================================================================
suspects = ["bit = 0", "bit=0", "2⁴", "2^4", "binary", "ON / OFF"]
for i, slide in enumerate(prs.slides, 1):
    for shape in slide.shapes:
        if not shape.has_text_frame: continue
        tf = shape.text_frame
        for para in tf.paragraphs:
            for run in para.runs:
                for s in suspects:
                    if s in run.text:
                        # replace with neutral language
                        if s == "ON / OFF":
                            run.text = run.text.replace(s, "real / baseline")
                        elif s in ("2⁴", "2^4"):
                            run.text = run.text.replace(s, "16")
                        elif s == "binary":
                            run.text = run.text.replace(s, "")
                        else:
                            run.text = run.text.replace(s, "")
                        print(f"Slide {i}: scrubbed '{s}'")
    if slide.has_notes_slide:
        nt = slide.notes_slide.notes_text_frame
        for para in nt.paragraphs:
            for run in para.runs:
                for s in suspects:
                    if s in run.text:
                        if s == "ON / OFF":
                            run.text = run.text.replace(s, "real / baseline")
                        elif s in ("2⁴", "2^4"):
                            run.text = run.text.replace(s, "16")
                        elif s == "binary":
                            run.text = run.text.replace(s, "")
                        else:
                            run.text = run.text.replace(s, "")
                        print(f"Slide {i} (notes): scrubbed '{s}'")

prs.save(DST)
print(f"\nSaved {DST}")
