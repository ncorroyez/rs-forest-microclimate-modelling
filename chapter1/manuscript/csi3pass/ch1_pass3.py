# -*- coding: utf-8 -*-
"""Third batch: sentences of 38 to 41 words, matched across line breaks."""
import io, re
F = ("/home/corroyez/Documents/z_Example_rmusica_31012025/Chapitre1/redaction/"
     "manuscript_chap1_EN_native20_csi3pass_2026-09-02.md")
E = []
def e(fault, old, new): E.append((fault, old, new))

e("41 mots",
 "The attribution rests on the typology only weakly: a principal component analysis of the variables recovers the archetypes at the 1 m vertical binning, with and without VCI (Figure 3), and the ranking does not depend on the archetype label (Fig. 5).",
 "The attribution rests on the typology only weakly. A principal component analysis of the variables recovers the archetypes at the 1 m vertical binning, with and without VCI (Figure 3). The ranking does not depend on the archetype label (Fig. 5).")

e("40 mots",
 "VCI, the normalized Shannon entropy of the vertical distribution of returns above ground [@vanewijkCharacterizingForestSuccession2011], summarizes the evenness of that distribution rather than the shape of the leaf-area-density profile, so we report it only as a descriptor of structural complexity (Fig. 1) and never as a model input; its definition, and what it tracks, are given in Appendix E and Section 4.1.",
 "VCI is the normalized Shannon entropy of the vertical distribution of returns above ground [@vanewijkCharacterizingForestSuccession2011]. It summarizes the evenness of that distribution, not the shape of the leaf-area-density profile, so we report it only as a descriptor of structural complexity (Fig. 1) and never as a model input. Appendix E and Section 4.1 give its definition and what it tracks.")

e("39 mots, métaphore positionnelle",
 "Where the foliage sits sets the height at which radiation is intercepted and where the canopy's heat sources and sinks lie, so two canopies of equal leaf quantity but different vertical arrangement could buffer the understory to different depths.",
 "The height of the foliage sets where radiation is intercepted and where the canopy's heat sources and sinks lie. Two canopies of equal leaf quantity but different vertical arrangement could therefore buffer the understory to different depths.")

e("39 mots",
 "cLHS *selects* among the existing 20 m LiDAR grid cells rather than generating variable values, so each of the 400 plots perturbed in Section 2.4 is a real, observed pixel with its own measured LAD profile and variable combination.",
 "cLHS *selects* among the existing 20 m LiDAR grid cells, and does not generate variable values. Each of the 400 plots perturbed in Section 2.4 is therefore a real, observed pixel, with its own measured LAD profile and variable combination.")

e("39 mots",
 "The macroclimatic reference is the local open-field station that also drives the model (CHS 41, air temperature at 1.5 m [@grilUsingAirborneLiDAR2023]; Section 2.4), rather than the free-air reanalysis aloft, so the offset and the driving meteorology share one baseline.",
 "The macroclimatic reference is the local open-field station that also drives the model (CHS 41, air temperature at 1.5 m [@grilUsingAirborneLiDAR2023]; Section 2.4), and not the free-air reanalysis aloft. The offset and the driving meteorology therefore share one baseline.")

e("38 mots",
 "Structure that once had to be read from field inventories can now be retrieved over continuous landscapes by remote sensing, and airborne laser scanning (ALS/LiDAR) resolves the full vertical distribution of vegetation rather than a two-dimensional projection [@bouvier7_GeneralizingPredictiveModels2015].",
 "Structure that once had to be read from field inventories can now be retrieved over continuous landscapes by remote sensing. Airborne laser scanning (ALS/LiDAR) resolves the full vertical distribution of vegetation, not a two-dimensional projection [@bouvier7_GeneralizingPredictiveModels2015].")

e("38 mots",
 "Because an index that tracks density can be mistaken for evidence that arrangement matters, VCI is a descriptor here and never a model input (Section 2.2); profile shape is characterized by the measured LAD profile itself (Section 2.3).",
 "An index that tracks density can be mistaken for evidence that arrangement matters, so VCI is a descriptor here and never a model input (Section 2.2). Profile shape is characterized by the measured LAD profile itself (Section 2.3).")

s = io.open(F, encoding="utf-8").read()
miss = []
for fault, old, new in E:
    pat = re.compile(r"\s+".join(re.escape(w) for w in old.split()))
    m = pat.search(s)
    if not m: miss.append((fault, old[:60])); continue
    s = s[:m.start()] + new + s[m.end():]
io.open(F, "w", encoding="utf-8").write(s)
print(f"{len(E)-len(miss)}/{len(E)} corrections appliquées")
for f, o in miss: print("  NON TROUVE |", f, "|", o)
