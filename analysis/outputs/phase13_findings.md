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


---

# Step 3 — QM confirmation (real CREST/GFN-FF ensembles)

Parsed all 12 candidates' real CREST ensembles (75–528 conformers each) with the
Phase-6 engine (Boltzmann-weighted at GFN-FF energies), and compared to the
native tail **PAPU-0080** taken from the Phase-8 known set at identical quality
and settings. Same core across all → whole-molecule SASA differences are the tail.

**Native baseline:** exposed hydrophobic SASA **622 Å²**, polar
**451 Å²**, hydrophobic fraction **0.58**.

## Result — only two tails actually improve at QM quality
- **C8_omega_sulfonate** (C8): hydrophobic SASA 530 (-93 vs native), polar 462 (+11), hydrophobic fraction 0.53 (-0.046)
- **oxa_PEG3** (C8): hydrophobic SASA 543 (-80 vs native), polar 420 (-32), hydrophobic fraction 0.56 (-0.016)

Every other tail is equal-or-WORSE than native on hydrophobic fraction once the
full conformer ensemble is accounted for. The fast ETKDG proxy actually RANKED the
tails well (Spearman ρ=0.82, p=0.00 vs the CREST hydrophobic fraction) — it
is a valid cheap screen. What it got wrong was the NATIVE BASELINE: it scored
native from a single ETKDG conformer (hydrophobic fraction ~0.73), which overstates
exposed hydrophobe, so many tails looked like wins. Against the properly ensembled
597-conformer CREST native (0.58), only two survive. Lesson (echoing Phase 7→8):
screen with the proxy, but judge "beats native" only against a same-fidelity native
— which is what this step does. Mechanistically, **chain SHORTENING + heteroatom /
charged content wins; an ω-polar cap on a still-long C12 chain buries the cap and
re-exposes hydrophobe** (the C12 ω-OH/NH2/COOH tails are all ≥ native).

## Finalists to promote to GFN2 (+ xTB electronics)
1. **t01_C8_omega_sulfonate** — the only tail that both cuts hydrophobic SASA
   (~−90 Å²) and RAISES polar SASA above native; largest hydrophobic-fraction drop.
2. **t02_oxa_PEG3** — second-best; C8-length ether backbone, lowers hydrophobic
   SASA ~−80 Å² without a formal charge (a more conservative amphiphile than the
   sulfonate).
3. **t07_C8_saturated** — promote as a MECHANISTIC CONTROL: it isolates the pure
   chain-shortening effect (C8, no polar head) from the polar-head contribution in
   t01/t02. If t01 ≫ t07 the polar head matters; if similar, length is doing the work.

## Exact next steps on your platform (Step 4)
For each finalist directory, GFN2 re-rank the existing ensemble (no re-search),
then GFN2-xTB electronics in water + octanol for QM logP:
```
cd analysis/outputs/phase13_qm_runs/t01_C8_omega_sulfonate
crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 6 --T 52   # -> crest_ensemble.xyz
# then Phase-9 electronics (see analysis/gen_known_xtb_inputs.py + phase9_electronic.py)
```
Upload `crest_ensemble.xyz` per finalist; re-run this ranking on the GFN2 set to
confirm the ordering holds before committing to synthesis. Note the sulfonate/
phosphocholine carry a formal charge — set `--chrg` accordingly in GFN2/xTB
(the fast SASA proxy and GFN-FF ran them neutral).

## Honest caveat
These are exposed-surface descriptors — the *hypothesis* for serum tolerance, not
the endpoint. Step 3 narrows 12 tails to 2 (+1 control) worth carrying forward; the
serum SHIFT of that small set, measured in vitro, is what actually tests the lead.
