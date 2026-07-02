# Phase 12 — serum-tolerance-biased generative design (findings)

## What was generated
Off the serum-active lead **PAPU-0080**, three branches:
- **ester** — re-esterify the validated aromatic C-6' handle (Phase 5 chemistry).
- **polaraxis** — designed acyls deliberately spanning hydrophobic -> polar.
- **notail** — the fatty acyl tail removed (deacylation) before re-esterifying (ibrexafungerp-inspired tail-free analogs).

30 unique analogs (28 novel), scored by the
reward below. 15 tail-free analogs generated.

## The reward (and its honest status)
Reward = **mean exposed POLAR surface fraction** over a 2-conformer ETKDG
ensemble (RDKit rdFreeSASA) — the operational form of the Phase 8/9/11 lead
("expose polar not hydrophobic surface"). It is a HYPOTHESIS, not a validated
oracle, and a crude single-few-conformer proxy that Phase 7->8 showed will
OVERSTATE effects. Bulk QED/Ro5/clogP terms were deliberately dropped (Phases 5,
11 showed them uninformative here).

Retrospective check on the 24 knowns: exposed polar fraction vs serum shift
Spearman **rho=-0.33 (p=0.11, n=24)** — the expected direction (more exposed polar surface -> smaller serum shift). Treat as a
sanity check on the reward's sign, not proof.

## Top novel candidates (full table: phase12_generated_library.csv)
- **notail::polaraxis_ax8_seryl** — exposed polar frac 0.48, novelty 0.49, MW 714
- **notail::polaraxis_ax9_polyhydroxy** — exposed polar frac 0.47, novelty 0.49, MW 775
- **notail::polaraxis_ax7_succinamoyl** — exposed polar frac 0.46, novelty 0.47, MW 726
- **notail::polaraxis_ax4_glycolyl** — exposed polar frac 0.46, novelty 0.47, MW 685
- **notail::polaraxis_ax6_carboxyethanoyl** — exposed polar frac 0.45, novelty 0.47, MW 727

## The deliverable that matters: a DISCRIMINATING SERIES
`phase12_discriminating_series.csv` is a small novel set that SPANS the
exposed-polar axis (low / mid / high) at comparable size on one scaffold. Its
purpose is to let a serum assay **falsify or confirm** the polar-surface lead:
if serum tolerance rises monotonically with exposed polar fraction across this
matched series, the lead is real and prospective; if not, the hypothesis is
rejected on-scaffold. This is worth more than any single "best" molecule.

## How to use this
1. Push `phase12_top_candidates.sdf` through the existing CREST -> Phase 6/8/9 QM
   funnel to confirm the exposed-polar property at QM-ensemble quality (the reward
   here is only a fast proxy).
2. Synthesize the **discriminating series** (SA set aside per instruction) and run
   the Phase-11 assay playbook — protein-adjusted MIC with albumin titration +
   equilibrium-dialysis fraction-unbound — reporting the serum SHIFT, not raw MIC.
3. Feed the measured shifts back as labels to train **Track B** (a REINVENT-style
   generative network on the external FKS pretraining set) with this same reward —
   at which point the reward stops being a hypothesis and becomes data-anchored.

## Caveats
Single-scaffold, proxy reward, no validated oracle; the tail-free branch produces
chemically aggressive structures by design (synthetic accessibility set aside).
The value is a testable, hypothesis-spanning series plus a reusable reward — not a
finished drug candidate.
