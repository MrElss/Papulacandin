# Phase 8 — the 3D/QM hypothesis at real CREST quality: a potency artifact

**Setup.** Phase 7 used fast RDKit-MMFF ensembles for the 24 known serum-gap
compounds and found Boltzmann-weighted hydrophobic SASA tracked serum MIC at
Spearman ρ = −0.45 (p = 0.029), beating the 2D baseline (ρ = 0.32, n.s.). That
was promising enough to justify spending cluster time on real CREST/GFN-FF
ensembles. Phase 8 ran them — all **24** compounds with healthy multi-hundred
conformer ensembles (132–597 each; PAPU-0078/6h required a `--notopo` rerun to
get past a biphenylene topology-filter artifact, then gave 310 conformers) — and
re-tested the same correlations with the identical phase6 descriptor engine.

## Result 1 — the signal shrinks and loses significance at QM quality

| Descriptor vs log₁₀(serum MIC), n=24 | MMFF proxy (Ph7) | real CREST (Ph8) |
|---|---|---|
| hydrophobic_sasa_mean | −0.45 (p=0.029) | **−0.31 (p=0.135)** |
| rg_mean | −0.39 (p=0.062) | −0.35 (p=0.096) |
| hydrophobic_sasa_std (rigidity) | −0.16 (n.s.) | +0.17 (n.s.) |
| polar_sasa_mean | −0.16 (n.s.) | −0.06 (n.s.) |

The *direction* reproduces, but the better ensembles give a weaker effect that
**no longer clears p = 0.05**. The MMFF-proxy ρ = −0.45 (p = 0.029) was an
over-estimate from cruder conformers. Group separation (serum-tolerant vs
serum-killed) is non-significant for every descriptor (Mann–Whitney p ≥ 0.12).

## Result 2 (decisive) — it's intrinsic potency, not a 3D serum effect

The known compounds are **not** equipotent: their serum-free MIC ranges
0.78 → 100 µg/mL and by itself tracks serum MIC at **ρ = +0.79 (p < 0.001)**.
Intrinsic potency dominates the serum MIC. Controlling for it removes the 3D
signal entirely:

| Descriptor | partial ρ vs serum MIC \| serum-free | vs serum SHIFT (potency-free) |
|---|---|---|
| hydrophobic_sasa_mean | **+0.02 (p=0.92)** | +0.24 (p=0.25) |
| hydrophobic_fraction_mean | +0.19 (p=0.38) | +0.28 (p=0.19) |
| hydrophobic_sasa_std | −0.01 (p=0.98) | −0.11 (p=0.61) |
| polar_sasa_mean | −0.31 (p=0.13) | −0.33 (p=0.12) |
| rg_mean | −0.23 (p=0.28) | −0.05 (p=0.80) |

Once serum-free potency is partialled out, **hydrophobic SASA's correlation with
serum MIC vanishes (ρ = 0.02)**. The apparent Phase-7 signal was compounds that
are bigger / more aromatic being more potent to begin with, and that potency
carrying through into serum — not a 3D serum-tolerance property.

## What *is* (weakly) left, and in which direction

Tested against the **serum SHIFT fold** — the pure, baseline-potency-independent
serum-tolerance endpoint — nothing is significant at n = 24, but the signs flip
to the **mechanistically expected** albumin direction:
- more exposed **hydrophobic** surface → **larger** shift (worse tolerance): ρ = +0.24
- more exposed **polar** surface → **smaller** shift (better tolerance): ρ = −0.33

**polar_sasa_mean** is the single most consistent and coherent lead (partial
−0.31, shift −0.33, both p ≈ 0.12–0.13): exposing polar rather than hydrophobic
surface is the only effect pointing the way the albumin-sequestration model
predicts. It is *not* significant here, but unlike hydrophobic SASA it is not a
potency artifact, so it is the descriptor worth testing with more data.

## Bottom line for the project

- **The 3D/QM funnel, done rigorously, does not deliver a serum-tolerance
  predictor on this dataset.** Ensemble-averaged SASA/shape and the rigidity
  (SASA-spread) metric do not predict the serum shift once potency is controlled.
  The Phase-1 rigidity hypothesis and the Phase-7 hydrophobic-SASA result both
  fail confirmation at CREST quality.
- **Raw serum MIC is the wrong endpoint** — it is dominated by intrinsic
  potency (ρ = 0.79). Future serum-tolerance modelling must use the serum SHIFT
  (or serum-free MIC as a covariate), as done here.
- **The honest state of the evidence:** the only coherent, non-confounded trend
  is "expose polar, not hydrophobic, surface" (polar SASA vs shift ρ = −0.33),
  and it needs more uncensored serum data to confirm — 11/24 serum MICs are
  censored at 100, which caps the achievable correlation.
- **Where the QM funnel still pays off:** not for shape, but potentially for
  *electronics*. The finalist DFT step (B3LYP ESP/dipole, Phase 6/8b) gives
  descriptors this analysis never had — surface electrostatics / H-bonding
  capacity — which are the natural next variables to test against the serum
  shift once the finalist single points complete.

## Methods / caveats
- Real CREST/GFN-FF ensembles, `--gfnff --alpb water --quick -ewin 6`; phase6
  in-house Shrake–Rupley SASA, Boltzmann weights at 298.15 K.
- n = 24, all compounds with proper ensembles. PAPU-0078/6h initially failed
  (its strained **biphenylene** ring made CREST's CREGEN topology filter discard
  all 25,760 sampled conformers); a `--notopo` rerun recovered a normal 310-conformer
  ensemble. Folding the real 6h in (vs the earlier n=23) shifts no number
  materially and changes no conclusion.
- Endpoint heavily censored/tied (serum MIC ∈ {12.5, 25, 50, 100}); Spearman
  rests on ties, so all ρ are trends, not clean separations.

## Artifacts
- `phase8_known_crest_descriptors.csv` — real-CREST 3D descriptors (24 compounds)
- `phase8_retrospective_crest.csv` — descriptors + serum data merged
- `phase8_confound_analysis.csv` — partial correlations + shift correlations
- `phase8_retrospective_crest.png` — best descriptor vs serum MIC; |ρ| vs MMFF proxy & 2D
