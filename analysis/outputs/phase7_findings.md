# Phase 7 — does the 3D/QM layer actually explain the serum gap?

**Question.** Phases 1–5 established that flat 2D descriptors (rotatable bonds,
clogP, QED, the composite design-rule score) only weakly track the observed
serum MIC of the fusacandin C-6′ ester series (Spearman ρ = 0.32, p ≈ 0.20,
n.s.). Phase 6 proposed that the discriminating signal is 3D/conformational —
specifically the Boltzmann-weighted exposed **hydrophobic SASA** and its
ensemble **spread** (a rigidity metric). Phase 7 tests that proposal directly.

**Design.** Same compounds and same endpoint as the 2D baseline: the Yeung-1996
fusacandin analogs 6a–6u plus the reference natural products (n = 24), all with
matched serum-free / serum whole-cell MIC vs *C. albicans*. Each known compound
was given a self-consistent RDKit ETKDGv3 + MMFF94 conformer ensemble (40 confs,
pruned), written in CREST `crest_conformers.xyz` format, and pushed through the
**identical** Phase-6 SASA/shape pipeline used on the 12 novel candidates.
Rank correlations are internally valid because every compound is treated by one
method on one shared scaffold. (This is an MMFF *proxy* for CREST — see Limits.)

## Result

| 3D descriptor | Spearman ρ vs log₁₀(serum MIC) | p | vs 2D baseline |
|---|---|---|---|
| **hydrophobic_sasa_mean** | **−0.45** | **0.029** | beats 0.32, now significant |
| rg_mean | −0.39 | 0.062 | beats 0.32 (borderline) |
| hydrophobic_sasa_std (rigidity) | −0.16 | 0.453 | no |
| polar_sasa_mean | −0.16 | 0.442 | no |
| hydrophobic_fraction_mean | −0.12 | 0.564 | no |

**A 3D descriptor does beat the 2D baseline** — hydrophobic SASA reaches
significance (p = 0.029) where the 2D score never did. But two results overturn
the original mechanistic story:

1. **The sign is backwards from the albumin-sequestration hypothesis.** ρ = −0.45
   means *more* exposed hydrophobic surface goes with a *lower* (better) serum
   MIC. The five most serum-tolerant analogs (6j/6k/6n/6o/6p, serum MIC 12.5)
   carry the **largest** extended-aromatic acyl groups and have the highest
   hydrophobic SASA (874–967 Å²) and Rg (6.3–6.9 Å). The dominant axis is
   aromatic **extent / molecular size**, which co-varies with intrinsic FKS1
   potency (these biphenyl/terphenyl esters are the most potent serum-free, MIC
   0.78), so they retain activity even under serum sequestration. Hydrophobic
   SASA is reading potency, not an albumin-avoidance handle.

2. **The rigidity/flexibility hypothesis does not validate.** Ensemble SASA
   spread is flat against serum MIC (ρ = −0.16, n.s.), and the serum-tolerant
   vs serum-killed groups are not separated on any descriptor (Mann–Whitney
   p ≥ 0.21). The "rigid → serum-tolerant" intuition from Phase 1 is not
   supported by the 3D MMFF proxy on this set.

## What this means for the project

- The serum gap in this series is **not** cleanly explained by exposed
  hydrophobic surface or conformational flexibility. To first order it tracks
  the same thing intrinsic potency does: how much rigid aromatic is hung off
  C-6′. "Minimize hydrophobic SASA to dodge albumin" would be the **wrong**
  design rule here — the data point the other way.
- Practically this *agrees* with the Phase-5 multi-objective ranking, whose top
  novel candidates are exactly the extended rigid aromatics (quinoline,
  naphthoyl, pyridylphenyl). So finalist selection is unchanged; what changes is
  the *rationale* and the descriptor we carry forward (hydrophobic SASA / size,
  not SASA-spread rigidity).
- The cleanest unconfounded test — does any 3D descriptor predict serum
  tolerance *after* controlling for intrinsic (serum-free) potency? — needs the
  serum-free MIC as a covariate and ideally true CREST ensembles. That is the
  recommended next experiment (see Limits).

## Limits (why this is a screen, not the verdict)

- **MMFF, not CREST.** Known-compound ensembles are MMFF94, not GFN-FF/GFN2.
  Absolute SASA values are not comparable to the candidates' CREST values; only
  within-known-set ranks are used. A positive signal here is exactly the
  justification to spend cluster time CREST-ing the known set for confirmation.
- **Censored, tied endpoint.** 11/24 serum MICs are censored at 100 (`=100`/`>100`)
  and values cluster at 12.5/25/50/100, so Spearman rests on heavy ties; ρ is a
  trend, not a clean separator (note 6q/6a have high SASA yet serum MIC 100).
- **Potency confound, uncontrolled.** Serum MIC mixes intrinsic potency and serum
  sequestration. The next analysis should regress serum MIC on serum-free MIC +
  3D descriptors to isolate a genuine serum-tolerance term.

## Artifacts
- `phase7_retrospective_qm.py` — pipeline
- `phase7_known_qm_descriptors.csv` — 3D descriptors per known compound
- `phase7_retrospective_qm.csv` — descriptors + serum MIC merged
- `phase7_retrospective_qm.png` — best descriptor vs serum MIC + |ρ| bar vs 2D baseline
