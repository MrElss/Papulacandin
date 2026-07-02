# Phase 13 — fatty-tail optimization, round 1 (findings + protocol)

## Scope
Core (sugar + spiroketal + aromatic C-6' acyl) FROZEN; only the long-chain fatty
acid is varied. Template **PAPU-0080**; native tail cleaved:
`CCCCC/C=C/C=C/[C@@H](O)C/C=C/C=C/C(=O)[18*:1]` (a C16 polyene). 15 tail
variants scored (native = control), all sharing one identical core.

## Reward (reused from Phase 12, validated rho=-0.33 on 24 knowns)
Mean exposed POLAR surface fraction over an ETKDG ensemble (rdFreeSASA). Because
the core is identical across all variants, differences in this readout — and in
the absolute exposed HYDROPHOBIC SASA — are attributable to the tail. Goal:
raise exposed polar fraction / lower exposed hydrophobic SASA vs the native tail
(the Phase-8 direction), WITHOUT touching the pharmacophore.

## Result (fast proxy — to be confirmed by CREST on your cluster)
Native baseline: exposed polar fraction **0.27**,
exposed hydrophobic SASA **932 A^2**.
11 novel tails beat the native baseline on the reward. Top:

- **C8_omega_sulfonate** — exposed polar frac 0.33 (native 0.27), exposed hydrophobic SASA 776, C8, novelty 0.38
- **oxa_PEG3** — exposed polar frac 0.32 (native 0.27), exposed hydrophobic SASA 786, C8, novelty 0.37
- **mid_chain_diol_C12** — exposed polar frac 0.31 (native 0.27), exposed hydrophobic SASA 830, C11, novelty 0.31
- **C12_omega_CONH2** — exposed polar frac 0.30 (native 0.27), exposed hydrophobic SASA 842, C12, novelty 0.35
- **C12_omega_COOH** — exposed polar frac 0.30 (native 0.27), exposed hydrophobic SASA 843, C12, novelty 0.35
- **amide_split_C12** — exposed polar frac 0.29 (native 0.27), exposed hydrophobic SASA 868, C12, novelty 0.34

Patterns to expect and read from the QM confirmation: shorter / terminally-capped
/ heteroatom-broken tails reduce exposed hydrophobic area; charged heads
(sulfonate, phosphocholine-like) raise polar exposure most but are the most
aggressive chemically (synthetic accessibility deliberately set aside).

## Deliverable: the discriminating series
`phase13_discriminating_series.csv` — one core, tails spanning LOW->HIGH exposed
polarity. Synthesize/assay this set so the serum SHIFT can be regressed on the
tail's exposed polarity DIRECTLY, on a fixed scaffold. That regression is the
round-1 test of the polar-surface lead.

## EXACT computational steps (run on your platform, upload results)
CREST/xTB are CPU codes (semiempirical GFN-FF/GFN2) — parallelize over cores, not
GPU. Full protocol in `phase13_qm_runs/SUBMIT.md`; summary:
1. **CREST ensembles** (GFN-FF, ALPB water, `--quick -ewin 6`, 52 cores) for each
   `phase13_qm_runs/<cand>/<cand>.xyz` -> `crest_conformers.xyz`. Screening tier;
   full-GFN2 search is intractable for this size.
2. **Upload** each `crest_conformers.xyz` back to the repo.
3. **Parse** at QM quality with the existing Phase-6 engine (command in SUBMIT.md)
   -> `phase13_qm_descriptors.csv`; compare polar/hydrophobic SASA vs native.
4. **Finalists**: GFN2 re-rank (`crest --screen ... --gfn2`) + Phase-9 xTB
   electronics (QM logP in water/octanol) before committing to synthesis.

## After round 1
Feed the measured serum shifts back as labels; that turns the reward from a
hypothesis into data and seeds Track B (a generative network over tail space).
The GPU on your platform is for THAT step, not for the CREST search here.
