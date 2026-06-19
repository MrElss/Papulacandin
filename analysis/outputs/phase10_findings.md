# Phase 10 — explicit HSA docking does not explain the serum gap

**Goal.** Test the mechanistic hypothesis directly: if serum loss is albumin
sequestration, compounds that bind human serum albumin (HSA) more strongly should
show a larger serum SHIFT. Phases 8–9 only ever gave a weak, indirect lead
(expose polar not hydrophobic surface); Phase 10 dockss the compounds to HSA and
correlates binding with the shift.

## Method (and why this one)
These ligands are 1000–1200 Da glycolipids with **~38 rotatable bonds** — outside
the regime where docking is reliable. We confirmed this empirically:
- **rigid into a Sudlow drug pocket** → +42 kcal/mol (the ligand is far too big to
  fit a pocket sized for ~400 Da drugs; pure clash);
- **flexible (38 torsions)** → >3 min/ligand and badly under-sampled (~2.5× Vina's
  reliable torsion ceiling) — intractable and unreliable for 24×N.

So we used **rigid ensemble surface docking**: HSA chain A (PDB 1AO6, defatted)
as a rigid receptor; the top Boltzmann CREST/xTB conformers (real QM geometries)
kept rigid; large (40 Å) boxes over the two Sudlow regions (subdomains IIA/IIIA)
so the rigid amphiphile settles into the nearest surface groove. Score = best
(most negative) Vina affinity over {≤3 conformers × 2 sites} = strongest
accessible HSA association. 24 compounds, 67 conformers, 134 docks.

## Result — no relationship (n = 24)

| HSA best-affinity vs | Spearman ρ | p |
|---|---|---|
| **serum SHIFT** | **+0.22** | 0.30 |
| serum MIC | −0.02 | 0.94 |
| serum-free MIC | −0.16 | 0.47 |
| partial(shift \| serum-free potency) | +0.16 | 0.46 |
| partial(shift \| polarizability/size) | +0.22 | 0.30 |

**HSA binding does not predict the serum shift.** It is not significant, and the
sign is even *opposite* to the sequestration prediction (we expected stronger
binding → larger shift, i.e. negative ρ; we got +0.22). The relationship is flat
in every slicing — raw, potency-controlled, and size-controlled.

One methodological positive: the docking score is **not** merely a size proxy
(HSA affinity vs polarizability ρ = −0.15, n.s.), so the large-box rigid surface
docking captured something beyond contact area — that something just doesn't track
serum tolerance.

## Interpretation
Two non-exclusive readings, both honest:
1. **The model is too crude for this chemotype.** Vina's scoring function is
   parameterized for drug-sized ligands, not 1200 Da amphiphilic glycolipids;
   rigid surface docking ignores induced fit and the fatty-acid binding sites
   (FA1–7) that these lipid-tailed molecules might actually use, and ignores
   other serum carriers (α1-acid glycoprotein, lipoproteins). A null here cannot
   exclude albumin sequestration — it can only say *this* computable proxy for it
   carries no signal.
2. **HSA drug-site binding may simply not be the discriminating variable.** The
   serum-tolerant and serum-killed analogs bind the modeled HSA surface with
   indistinguishable strength (−11 to −15 kcal/mol, fully overlapping).

## Where this leaves the project
Across every computable angle now tried — 2D descriptors (ρ=0.32), 3D shape/SASA
(potency artifact), QM electronics/solvation (converge weakly on polar-surface,
|ρ|≈0.30 n.s.), and explicit HSA docking (null) — **no approach yields a
significant, mechanistically-clean serum-tolerance predictor on this dataset.**
The consistent, defensible takeaways are:
- serum *potency* is governed by intrinsic potency × molecular size/lipophilicity;
- serum *tolerance* (the shift) has, at most, a weak association with exposing
  polar rather than hydrophobic surface, seen in two independent descriptor
  families but never crossing significance.

The binding constraint is the **data**, not the methods: n = 24 with 11/24 serum
MICs censored at 100 and four tied levels. The highest-value next step is no
longer another descriptor or docking variant but **experimental** — more analogs
with *uncensored* serum MICs (ideally a fresh chemotype), or a direct
HSA-binding assay (e.g. equilibrium dialysis / fluorescence displacement) to
ground-truth the sequestration hypothesis that docking could not resolve.

## Reproduce
Receptor prep (one-time): download `1AO6.pdb`, keep chain A protein atoms →
`hsa_chainA.pdb`, then `mk_prepare_receptor.py --read_pdb hsa_chainA.pdb -o
hsa_receptor -p`. Box centers are the Sudlow site I/II residue centroids
(in `phase10_dock_hsa.py`). Then `python3 phase10_dock_hsa.py`.

## Artifacts
- `phase10_dock_hsa.py` — pipeline (ligand prep + rigid ensemble docking + stats)
- `phase10_docking/phase10_hsa_scores.csv` — per conformer/site affinities
- `phase10_docking/phase10_hsa_per_compound.csv` — best/mean HSA affinity + serum data
- `phase10_hsa.png` — HSA affinity vs serum shift
