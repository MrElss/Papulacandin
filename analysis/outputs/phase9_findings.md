# Phase 9 — electronic / solvation descriptors: convergence, not a breakthrough

**Setup.** Phase 8 closed the shape hypothesis: ensemble SASA/shape does not
predict serum tolerance once intrinsic potency is controlled; the one coherent
(weak) lead was "expose polar, not hydrophobic, surface" (polar SASA vs serum
shift ρ = −0.33, n.s.). Phase 9 tested the axis shape misses — **electronics /
solvation** — with cheap GFN2-xTB single points (water + octanol) on each known
compound's Boltzmann-populated conformers (24 compounds, 102 conformers, 204
single points), Boltzmann-averaged and run through the **same** Phase-8 statistics.

## Results (n = 24)

| Descriptor | ρ vs serum MIC | partial ρ \| potency | ρ vs serum SHIFT |
|---|---|---|---|
| **polarizability α(0)** | **−0.54 (p=0.01)** | −0.37 (p=0.07) | −0.02 (p=0.94) |
| **QM logP** (octanol/water) | −0.31 (p=0.14) | +0.09 (p=0.67) | **+0.30 (p=0.15)** |
| dipole moment | −0.14 (p=0.51) | −0.25 (p=0.24) | −0.15 (p=0.49) |
| aqueous G_solv | +0.05 (p=0.80) | +0.19 (p=0.37) | +0.13 (p=0.55) |
| HOMO–LUMO gap | −0.13 (p=0.54) | −0.02 (p=0.92) | +0.14 (p=0.51) |

Two things stand out, and neither is a clean win:

### 1. Polarizability is the strongest serum-MIC correlate in the whole project — but it's molecular size
α(0) vs serum MIC is ρ = −0.54 (p = 0.01), the single biggest correlation we've
seen (beating 2D 0.32, MMFF hydrophobic SASA −0.45, CREST −0.31). But α(0) is a
**bulk/size descriptor**: it correlates with radius of gyration (+0.52),
hydrophobic SASA (+0.53), and — crucially — with **serum-free** potency
(−0.42, p = 0.04). Against the potency-independent serum **shift** it is
ρ = −0.02 (nothing). So α(0) just sharpens the Phase-8 conclusion — bigger, more
polarizable, more potent molecules have lower serum MIC — rather than revealing a
serum-tolerance handle. (Its partial ρ = −0.37, p = 0.07 is borderline but is not
corroborated by the shift endpoint, so we do not claim it.)

### 2. QM logP independently corroborates the polar-surface lead
The octanol/water QM logP tracks hydrophobic SASA at **ρ = +0.76** — the
electronic and shape pictures are measuring the same physical thing — and against
the serum shift it gives **ρ = +0.30** (more hydrophobic → larger shift → worse
serum tolerance). That mirrors Phase-8's polar SASA (−0.33: more polar → better
tolerance) in both direction and magnitude. Dipole points the same way weakly
(−0.15). None reach significance at n = 24, but **two independent descriptor
families (3D shape SASA, QM solvation logP) now converge on one qualitative
rule**, which is stronger evidence for that rule than either alone.

## Bottom line

- **Electronics did not break the significance barrier either** — no descriptor
  predicts the serum shift at p < 0.05 on this n = 24 / 11-censored dataset.
- **But the project's signal is now consistent across methods**: serum *potency*
  is dominated by molecular size/lipophilicity (α(0), Rg, hydrophobic SASA all
  track serum-free and serum MIC), while serum *tolerance* (the shift) is weakly
  but repeatably associated — across both shape and electronic descriptors — with
  exposing **polar rather than hydrophobic** surface (|ρ| ≈ 0.30–0.33, p ≈ 0.12–0.15).
- **This is the actionable design hypothesis the project can defend**: among
  equipotent analogs, bias the exposed surface toward polar/H-bonding character
  (lower QM logP, higher polar SASA) to reduce serum loss — with the explicit
  caveat that it is a directional lead, not a validated predictor, and needs more
  uncensored serum data to confirm.

## Why no single descriptor reaches significance — and what would
The ceiling is the data, not the descriptors: n = 24 with 11/24 serum MICs
censored at 100 and values tied on {12.5, 25, 50, 100}. Three convergent
descriptors at |ρ| ≈ 0.30 with p ≈ 0.13 is what a real-but-modest effect looks
like under that much censoring. Confirmation needs (a) more analogs with
uncensored serum MICs, ideally a fresh chemotype, and/or (b) the direct
mechanistic test — explicit binding to human serum albumin — rather than more
single-molecule descriptors.

## Artifacts
- `phase9_electronic_descriptors.csv` — per-compound dipole, gap, α(0), G_solv(water), QM logP
- `phase9_electronic_stats.csv` — ρ vs serum MIC, partial (potency-controlled), ρ vs shift
- `phase9_electronic.png` — |ρ vs serum shift| per descriptor vs the Phase-8 shape lead
