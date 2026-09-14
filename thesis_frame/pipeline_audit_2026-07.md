# Reviewer-2 audit — Chapter 1 (current state, 2026-07)

Fresh cold read of `chapter1/manuscript/manuscript_chap1_EN.md` and its pipeline, after this
session's changes (sensitivity not Shapley; v3.2.3 iter + station forcing; Fig 1↔2 swap;
appendices A–G with new B = ERA5 bias and G = wind; Appendix D switched gap-fraction → fCover).
Supersedes most of `review/pipeline_critique.md` (that was the Shapley / r=0.93 era).

Tiers: **T1 science at risk** · **T2 missing analysis** · **T3 presentation/provenance**.
Each item: severity · evidence · fix. **[NEW]** = introduced or changed this session.

---

## T1 — Science at risk

### 1.1 [HIGH] [NEW] The open-canopy failure is described with the wrong sign
The manuscript repeatedly says the model "reverses sign in the sparsest stands, holding the
understory near neutral **where the loggers cool**" (§3.2, §3.3, §4.2, Abstract), and Appendix D
now says the model "**under-cools** the open, low-cover stands."

**Evidence (data).** For the low-cover plots that dominate the residual (fCover < 0.6):
obs ΔT~max~ mean = **+3.81 °C** (strong **warming / amplification**), model = **+0.22 °C** (neutral);
per plot 41_18 obs +4.80/sim +0.21, 41_39 +4.53/+0.23, 41_49 +4.22/+0.18. High-cover plots:
obs −1.44, model −0.25. So the open failure is the model **failing to reproduce the observed
warming (amplification)**, not "loggers cool / model reverses sign," and it is **not** "under-cooling"
(the model is *cooler* than obs there, not warmer).

**Why it matters.** This is the headline description of the glass ceiling; it is currently backwards
for exactly the plots invoked. The floor mechanism (below) is coherent with the *correct* direction,
so only the prose needs fixing — but it must be fixed everywhere (Abstract, §3.2/3.3, §4.2, App D).

**Fix.** Reword to: "in the open, low-cover stands the loggers **amplify** (understory hotter than
the open), while the model holds them near neutral — it fails to reproduce the amplification." Replace
"under-cools the open stands" with "fails to warm / under-amplifies."

### 1.2 [HIGH] [NEW] The fCover-floor residual (Appendix D) is a weaker, partly self-inflicted argument than the gap fraction it replaced
Switching Fig D1 from an external gap-fraction metric to the model's **own fractional-cover input,
floored at 0.5**, regresses the residual on a *model input* mixed with external covariates. Part of
the correlation is mechanical: the 4 open plots are pinned at the floor and cannot respond by
construction, so "residual is largest at the floor" is partly "we clipped our own input." The old
gap-fraction version regressed on something the model genuinely **cannot see** — a cleaner
"1-D limitation" claim.

**Why it matters.** It invites the obvious reviewer move: "then lower the floor and re-run the open
plots." It is the single sharpest actionable thread in the paper and it connects straight to 1.1.

**Fix (choose one).**
- (a) Revert Appendix D to the external openness metric (gap fraction, or raw uncapped cover ≈ 1−gap),
  keeping the honest "structure the model can't see" framing; or
- (b) Keep fCover but **test lowering the floor**: re-run the 4 floored open plots at their real cover
  and report whether the amplification failure shrinks. If it does, that is a real, fixable model
  limitation (strong result); if it does not, the floor is not the cause and the wording in 1.1 must
  drop the floor mechanism.
Decide before submission; do not ship the current half-way version.

### 1.3 [HIGH] The "uniform offset cancels, so the ranking survives" defence is in tension with a metric-dependent model
The model captures ~13 % of the amplitude and reverses behaviour by regime, yet §4.2 defends the
ranking with "a fixed offset cancels in the up-minus-down difference." A **sign reversal in the open
regime** and a model that gets the **slope right (0.90≈0.90, r=0.78) but the amplitude wrong (13 %)**
are, by definition, **non-uniform / metric-dependent** compression — exactly what the offset argument
does not cover. The manuscript half-admits this (the "trait-dependent compression" caveat) but leaves
the tension unresolved.

**Mitigation that already exists — state it.** The ranking **is** shown in the faithful (closed) regime:
Fig 3 P4 has LAI −0.28 °C/SD vs profile ≈ 0. So "LAI ≫ profile" holds on the plots the model gets right,
independently of the open-canopy failure. Say this explicitly as the rebuttal, and lean the ranking on
the two **model-free** lines (Fig 5 forward, Fig 6 nested) rather than on MuSICA's radiative gain.

**Fix.** Reconcile "slope-right/amplitude-wrong" with "compression ≈ uniform" in one sentence; state
the closed-canopy-only ranking as the robustness check; make explicit that the ranking's load is borne
by the model-free evidence, not the under-buffering model.

---

## T2 — Missing analyses

### 2.1 [MEDIUM] No uncertainty on the central result (Fig 3)
The attribution violins are distributions, but there is **no CI/bootstrap** on the per-archetype median
sensitivity per trait, and no test that LAI > profile is significant. "LAI dominates, profile weakest"
rests on point estimates — a standard reviewer ask.
**Fix.** Bootstrap plots within each archetype, report a CI on each trait's median sensitivity, and
confirm the LAI interval does not overlap the profile interval (per archetype). One extra panel or a table.

### 2.2 [MEDIUM] Forward inclusion (Fig 5) is order-dependent
Fixed order [LAI, H~max~, fCover, LAD] guarantees LAD-adds-least. The **order-free** answer is Fig 6
(nested regression: structure adds ΔR² = 0.00). Right now Fig 5 carries an order-dependent claim and
Fig 6 rescues it only implicitly.
**Fix.** One sentence: name Fig 6 as the order-free confirmation of the forward-inclusion order; optionally
show LAD-first as a robustness inset.

### 2.3 [MEDIUM] Per-archetype forward *r* on n≈8 is near-noise
Fig 5's per-cluster lines and the "−0.09 to +0.22" per-archetype LAD increments sit on n = 8–15.
Already caveated as "small per-archetype subsets," but the coloured per-cluster lines invite
over-reading.
**Fix.** Annotate n per cluster, or grey the per-cluster lines and keep the pooled "All" as the claim.

### 2.4 [LOW] H2 controlled test (Fig C1) is at 4 centroids, not per plot
Everywhere else the design is per-plot to avoid Jensen bias (f(mean) ≠ mean(f)); the H2 test is at the
4 archetype centroids. Defensible (it is a counterfactual sweep), but the asymmetry should be named.
**Fix.** One clause acknowledging Fig C1 is centroid-based.

### 2.5 [LOW] Night-time buffering (ΔT~min~) not examined
Longwave trapping at night may weight height/profile differently than daytime shading; refugia care
about both extremes. Acknowledged as a limitation — keep, but consider flagging that the profile could
matter more at night (an honest bound on the "profile is negligible" headline).

---

## T3 — Presentation / provenance

### 3.1 [CLEARED, keep visible] Provenance is now consistent
Attribution (`metrics6_frblo`), validation and forward (`tab_hobo_frblo_validation.csv`,
`musica_hobo_z05_v323station_iter`, forcing `FR-Blo_2021_v2.nc`) all use the **same station forcing and
v3.2.3**. The prior critique's "three simulation branches" concern is resolved. The one **intentional**
difference is LAD resolution — 1 m for the attribution, 0.5 m for the validation/robustness — which §4.4
justifies as conservative. Keep that sentence; it is the one provenance seam a reviewer could pull.

### 3.2 [LOW] Number-consistency after the session's edits — OK, but re-sweep once
fCover r = −0.80 / adj R² = 0.71, wind z0 = 0.44 / factor ~0.4 / ~2.3×, Moran 0.07 (resid) / 0.15 (field)
are consistent across Abstract, §3.3, §4.2, App D, App G. The Abstract now reflects the fCover/wind
edits. Do a final numeric grep before submission.

### 3.3 [LOW] Internal filename/label crossing
Figure 1 → `Fig2_site_typology_composite.png`, Figure 2 → `Fig1_methodo.png`, and the appendix files are
letter-shifted (Fig C1 → `B1_h2_...png`, etc.). Invisible to readers, but a maintenance hazard.
**Fix (optional).** Rename files to match labels in one pass, or add a mapping note in the repo.

---

## Discussion points worth adding (the user's "pts de disc")

- **Open stands as heat traps, not refugia.** The amplification the model misses (1.1) is ecologically
  meaningful: open/gappy stands can be *hotter* than the open field. The management message ("keep leaf
  area and closure") is reinforced by this, and the model's blindness to it is a stated bound on using
  the model for management — currently framed only as error, could be framed ecologically.
- **Species/leaf-habit scope.** The "where foliage sits is redundant once quantity is known" result is
  for a *closed deciduous oak* canopy. Evergreen/conifer stands with deeper, more differentiated crowns
  are exactly where the vertical profile might still matter — a one-line honest hedge strengthens the
  claim rather than weakens it.
- **The cover floor as a concrete model-improvement target** (ties to 1.2): naming it turns a weakness
  into a perspective ("relaxing the sub-canopy cover floor is the first step to closing the open-canopy
  gap").

---

## Priority order before the jury / submission
1. **1.1** fix the open-canopy sign wording everywhere (quick, and currently wrong).
2. **1.2** decide gap-fraction-revert vs fCover-with-floor-test for Appendix D (the sharpest thread).
3. **1.3** reconcile slope-right/amplitude-wrong vs "uniform offset"; state the closed-canopy ranking.
4. **2.1** add a CI/bootstrap to Fig 3.
5. **2.2/2.3** name Fig 6 as the order-free answer; tame the n≈8 per-cluster lines.
