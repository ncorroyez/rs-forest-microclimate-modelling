# -*- coding: utf-8 -*-
"""CSI 3A pass on Chapter 1. Sentence-level only: no number, no citation key and no
section cross-reference is touched. Framing decision: `lever` is kept as a noun,
`exercise` and `potent` are removed."""
import io, sys
R = "/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/redaction/"
SRC = R + "manuscript_chap1_EN_native20_cites_2026-08-27.md"
DST = R + "manuscript_chap1_EN_native20_csi3pass_2026-09-02.md"

E = []   # (fault, old, new)
def e(fault, old, new): E.append((fault, old, new))

# ── Key message ──────────────────────────────────────────────────────────────
e("maniéré: carries, exercise",
  "the amount of foliage a\ncanopy carries, more than its vertical arrangement, controls the temperature buffering that\ntemperate oak stands realize in the understory. The vertical profile is a mechanistic lever real\ncanopies barely exercise.",
  "the amount of foliage, more than its\nvertical arrangement, controls the temperature buffering of temperate oak understories. The\nvertical profile has a mechanistic effect, and real canopies vary too little in shape to realize it.")

# ── Abstract ─────────────────────────────────────────────────────────────────
e("maniéré: potent, exercises",
  "The vertical profile is potent only when fully\nexercised: as a real-versus-uniform contrast it is the largest single lever in mid-density canopy, though not on\nthe buffering slope, and a controlled bottom- to top-heavy reshaping spans only 0.15 to 0.51 °C. Both are upper\nbounds: measured profiles are nearly invariant within a structural type, so the design gradient exercises\nonly a small part of either range.",
  "The vertical profile acts only at its full range. As\na real-versus-uniform contrast it is the largest single lever in mid-density canopy, though not on the\nbuffering slope, and a controlled bottom- to top-heavy reshaping spans 0.15 to 0.51 °C. Both are upper\nbounds. Measured profiles are nearly invariant within a structural type, so the design gradient spans\nonly a small part of either range.")

# ── Introduction ─────────────────────────────────────────────────────────────
e("maniéré: credited with primacy", "has often been credited with primacy", "is often ranked first")
e("métaphore positionnelle", "Two difficulties\nstand in the way.", "Two difficulties follow.")
e("maniéré: go unexercised",
  "We report both the size of each mechanism and the range real stands span,\nbecause a mechanism can exist and still go unexercised, and we keep the two separate.",
  "We report the size of each mechanism and the range real stands span, and keep the two separate. A\nmechanism can be large in the model and still vary little between real stands.")

# ── Methods ──────────────────────────────────────────────────────────────────
e("maniéré: buys", "The coupling buys a more complete", "The coupling gives a more complete")
e("pivot de négation (défensif)",
  "so a lever near a bound is understated rather than inflated",
  "so a lever near a bound is understated (never inflated)")

# ── Results ──────────────────────────────────────────────────────────────────
e("pivot de négation (rhétorique)",
  "That ordering is built in rather than recovered (Section 2.3).",
  "That ordering is built in (Section 2.3).")
e("énumération compressée, 84 mots",
  "The vertical arrangement of foliage is a potent lever, but only when fully exercised. Isolated as the real-versus-uniform contrast at fixed leaf area and height, it acts unevenly along the density gradient: it is small in the open P1 (0.01 °C); across the intermediate P2 and P3 it moves ΔT~max~ by up to 0.14 °C (P3, 0.14 °C; P2, 0.09 °C), exceeding the leaf-area step in 79% of P3 plots but only 58% of P2 plots; and in the dense P4 it is dwarfed by leaf quantity (0.08 versus 0.26 °C), its sign flipping from cooling to warming. On ΔT~max~, read as the full contrast, it is the largest single lever in P3.",
  "The vertical arrangement of foliage is a large lever at its full range. Isolated as the real-versus-uniform contrast at fixed leaf area and height, it acts unevenly along the density gradient. It is small in the open P1 (0.01 °C). Across the intermediate archetypes it moves ΔT~max~ by 0.09 °C in P2 and 0.14 °C in P3, exceeding the leaf-area step in 58% of P2 plots and 79% of P3 plots. In the dense P4 leaf quantity is three times larger (0.26 against 0.08 °C), and the profile contrast changes sign from cooling to warming. Read as the full contrast on ΔT~max~, the profile is the largest single lever in P3.")
e("maniéré: potency, exercises",
  "That potency is an upper bound: the full contrast carries the canopy all the way to the uniform limit, a range real canopies span only in part.",
  "That value is an upper bound. The full contrast takes the canopy all the way to the uniform limit, a range real canopies span only in part.")
e("maniéré: exercises", "so the design gradient exercises only a small part of the contrast",
  "so the design gradient spans only a small part of the contrast")
e("maniéré: firms up", "The profile contrast firms up in the intermediate archetypes",
  "The profile contrast grows in the intermediate archetypes")
e("maniéré: room to move", "Leaf quantity also has the most room to move within a structural type",
  "Leaf quantity also varies most within a structural type")
e("maniéré: fall as bins", "The archetypes fall as bins along it and cross the boundaries without a step (Fig. 5).",
  "The archetypes are bins along that function, and they cross the boundaries without a step (Fig. 5).")

# ── Discussion 4.1 ───────────────────────────────────────────────────────────
e("maniéré: potent",
  "The mechanism itself\nis potent: isolated as a real-versus-uniform contrast, the profile moves ΔT~max~, exceeding\nthe leaf-quantity step in the intermediate archetypes and staying small only at the density\nextremes.",
  "The mechanism itself is\nmeasurable. Isolated as a real-versus-uniform contrast, the profile moves ΔT~max~ by up to\n0.14 °C, exceeding the leaf-quantity step in the intermediate archetypes and staying small only at\nthe density extremes.")
e("maniéré: barely span",
  "The contrast is in any\ncase an upper bound real canopies barely span",
  "The contrast is an upper\nbound that real canopies span only in part")
e("maniéré: exercise",
  "H1 is not\nrealized: the mechanism exists, but the between-stand variation of real profiles does not\nexercise it.",
  "H1 is not\nrealized. The mechanism exists in the model, and the between-stand variation of real profiles is\ntoo small to express it.")
e("métaphore positionnelle: sits", "The open P1 sits lower,\nat a median of 0.40.",
  "The open P1 has a lower median,\n0.40.")
e("verbe grandiloquent: survive",
  "The primacy\noften attributed to the vertical dimension (Section 1) does not survive a controlled, on-manifold\ntest.",
  "A controlled, on-manifold test does not support the primacy often attributed to\nthe vertical dimension (Section 1).")
e("pivot de négation (défensif)",
  "The sign flip corroborates the mechanism rather than signaling an artifact.",
  "The sign flip corroborates the mechanism (it is not a solver artifact).")
e("métaphore positionnelle: sits", "In the dense P4 cover sits at 0.92 ± 0.03",
  "In the dense P4 cover is 0.92 ± 0.03")
e("pivot de négation (défensif)",
  "Saturation and truncation compress the dense end rather than inflate it.",
  "Saturation and truncation compress the dense end (they do not inflate it).")
e("maniéré: carry", "may carry the controlling structural information",
  "may resolve the controlling structural information")
e("maniéré: potent",
  "yet the profile proves potent in principle but barely exercised",
  "yet the profile is large as a contrast and small over the range real stands span")
e("maniéré: keeps mattering", "plant area keeps mattering through needle clumping",
  "plant area still matters, through needle clumping")

# ── Discussion 4.2 ───────────────────────────────────────────────────────────
e("méta-narration: that is what",
  "That is what allows such models to be pushed outside the conditions they were built in, as when",
  "Such models can therefore be applied outside the conditions they were built in, as when")

# ── Discussion 4.3 ───────────────────────────────────────────────────────────
e("pivot de négation (rhétorique)",
  "P1 is therefore not a young stand but a heterogeneous open class mixing sparse plots with the opened crowns of older stands.",
  "P1 is therefore a heterogeneous open class, mixing sparse plots with the opened crowns of older stands.")
e("maniéré: survives",
  "Some of those dependences are corrected here, and the ordering survives each correction.",
  "We correct some of those dependences here, and the ordering is unchanged by each correction.")
e("métaphore: margin to clear",
  "The margin such a proxy must clear is quantitative: a uniform 23% reduction in retrieved leaf area, the shift produced by the optimized extinction coefficient of Appendix F, removes the simulated coupling to the observed buffering (Corroyez et al., in preparation). An optical estimate must therefore resolve leaf quantity well inside that margin, which is the regime in which optical indices saturate.",
  "The requirement on such a proxy is quantitative. A uniform 23% reduction in retrieved leaf area, the shift produced by the optimized extinction coefficient of Appendix F, removes the simulated coupling to the observed buffering (Corroyez et al., in preparation). An optical estimate must therefore resolve leaf quantity to better than 23%, and that is the regime in which optical indices saturate.")

# ── Conclusion ───────────────────────────────────────────────────────────────
e("maniéré: potent, exercised",
  "The vertical profile is potent when fully exercised, largest in\nmid-density canopy on the offset and not on the buffering slope.",
  "The vertical profile acts at its full range, and is largest\nin mid-density canopy on the offset, not on the buffering slope.")
e("pivot de négation",
  "That lever grows rather than saturates as the canopy closes, reaching 0.51 °C at the dense end",
  "That lever grows as the canopy closes, without saturating, and reaches 0.51 °C at the dense end")
e("métaphore: margin to clear",
  "The requirement on such a measure is\nquantitative, and Section 4.3 sets the margin it must clear.",
  "The requirement on such a measure is\nquantitative, and Section 4.3 gives the value it must reach.")

s = io.open(SRC, encoding="utf-8").read()
miss = []
for fault, old, new in E:
    if old not in s:
        miss.append((fault, old[:70])); continue
    s = s.replace(old, new, 1)
io.open(DST, "w", encoding="utf-8").write(s)
print(f"{len(E)-len(miss)}/{len(E)} corrections appliquées")
for f, o in miss: print("  NON TROUVE |", f, "|", o)
