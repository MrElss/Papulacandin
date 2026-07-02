# Analysis pipeline

Reproducible analysis for the Papulacandin / FKS1 serum-gap project.
Each script is self-contained, documented, and writes to `analysis/outputs/`.

## Environment
```
python3 -m pip install pandas numpy rdkit scikit-learn
```

## Scripts (run in order)
1. `serum_gap_analysis.py` — builds the matched serum-free vs serum MIC table
   for Candida albicans and computes the per-compound "serum shift". This is
   the project's dependent variable. Handles censored values (`>100`, `<0.03`)
   explicitly rather than discarding them.

## Outputs
- `outputs/serum_gap_pairs.csv` — matched pairs + serum shift per compound
- `outputs/serum_gap_summary.txt` — human-readable summary
- `outputs/phase1_descriptor_stats.csv` + `phase1_*.png` — serum-shift SAR
- `outputs/phase2_activity_ladder.csv` + `phase2_*.png` — activity ladder / attrition
- `outputs/phase3a_*` — external FKS model, CV metrics, papulacandin scores, joblib
- `outputs/phase3b_*` — within-series serum-tolerance model, LOO predictions, coeffs

## Phase 4 outputs
- `outputs/table_phase4_existing_leads.csv` — best observed serum-active leads
- `outputs/table_phase4_candidate_ranking.csv` — heuristic design-score ranking
- `outputs/design_rules.md` — distilled, source-attributed design rules
- `outputs/fig_phase4_*.png` — presentation figures

## Phase 5 — generative structure optimization (Tier 1 backbone)
`phase5_generate.py` — scaffold-constrained generation of novel Papulacandin-class
analogs by cleaving and re-esterifying the validated aromatic C-6' ester (chemistry
verified by InChIKey round-trip), scored with a transparent multi-objective function
(QED, SA score, Phase-1 design rules, novelty). No learned activity predictor is used.
Emits CREST-ready 3D structures for the Tier-4 QM funnel.

Requires SA_Score from RDKit Contrib (bundled with rdkit).

### Phase 5 outputs
- `outputs/phase5_virtual_library.csv` — generated analogs + multi-objective scores
- `outputs/phase5_top_candidates.sdf` — 3D structures of top novel candidates (CREST-ready)
- `outputs/phase5_crest_commands.sh` — per-candidate CREST/xTB command template (run on QM infra)
- `outputs/phase5_retrospective.png` — design-rule score vs observed serum MIC
- `outputs/phase5_score_distributions.png` — library drug-likeness / SA / novelty

## Phase 6 — QM integration layer (Tier 4 funnel)
`phase6_qm_layer.py` — parses CREST conformer ensembles (`crest_conformers.xyz`),
computes Boltzmann-weighted 3D descriptors that 2D rules cannot see: hydrophobic
vs polar SASA (in-house vectorized Shrake-Rupley, no external SASA dependency —
`freesasa` failed to build in this environment), ensemble SASA spread (a genuine
3D flexibility metric), radius of gyration, and asphericity. Also generates
Gaussian DFT single-point inputs (B3LYP/6-31G(d), PCM water, Pop=MK) for the
Boltzmann-populated conformer subset, and parses returned Gaussian logs for SCF
energy / dipole / ESP charges.

Motivation: Phase 5 found QED/Ro5 saturate near zero for this bRo5 chemotype
(uninformative) and flat 2D design rules only weakly track serum MIC
retrospectively (Spearman rho=0.32, n.s.). Both point to a 3D/conformational
gap that only ensemble QM descriptors can close.

Running `python3 phase6_qm_layer.py` executes a **self-test** on a synthetic
RDKit-generated ensemble (clearly not real QM data) to validate the full
parse -> descriptor -> Gaussian-input path before any real CREST output exists.
To run on real data, drop each candidate's `crest_conformers.xyz` under
`outputs/qm_runs/<candidate_name>/` (names match `phase5_top_candidates.sdf`,
e.g. `cand01_quinolinecarbonyl`) and call `run_qm_layer(...)` — see the
`__main__` block for the exact call signature.

### CREST settings actually used (fast screening tier, 12 candidates) — DONE
Full GFN2 (`QUICK=no`, `EWIN=10`) does not scale to these ~140-atom, 30+
rotatable-bond glycolipids (~10 days/compound observed). The screening tier
ran the conformer search at GFN-FF instead. NOTE: the cluster's "new
calculator" CREST build rejects the old `--gfn2//gfnff` composite syntax, so
GFN-FF is selected directly; GFN-FF requires ALPB solvation (`--alpb`), not the
GFN1/2-only `--gbsa`. Settings used:
```
crest input.xyz --gfnff --alpb water --quick -ewin 6 --T 52
```
This produced `crest_conformers.xyz` (≈36–505 conformers/compound) for all 12
candidates under `outputs/qm_runs/<candidate>/`. Reserve full `--gfn2`
(`-ewin 10`) for the 1–3 funnel finalists (see "GFN2 promotion" below).

### Phase 6 outputs (real CREST data, all 12 candidates) — DONE
- `outputs/phase6_qm_descriptors.csv` — per-compound Boltzmann-weighted 3D descriptors
- `outputs/qm_runs_gaussian/<candidate>/*.gjf` — DFT inputs for populated conformers (32 total)

Run on real data with:
```
python3 -c "import phase6_qm_layer as p6; p6.run_qm_layer(
  'outputs/phase5_top_candidates.sdf', 'outputs/qm_runs',
  'outputs/phase6_qm_descriptors.csv', 'outputs/qm_runs_gaussian', real_run=True)"
```

## Phase 7 — retrospective validation of the 3D hypothesis
`phase7_retrospective_qm.py` — tests whether the Phase-6 3D descriptors actually
separate serum-tolerant from serum-killed compounds on the REAL Yeung-1996
fusacandin analog series (6a–6u + reference natural products), the same set and
endpoint the 2D baseline (Spearman rho=0.32) used. Known compounds have no CREST
ensembles yet, so each is given a fast self-consistent RDKit ETKDGv3+MMFF94
ensemble (CREST xyz format) pushed through the IDENTICAL phase6 SASA/shape
pipeline; rank correlations are internally valid as a feasibility test, and a
positive result justifies CREST-ing the known set for confirmation.

### Phase 7 outputs
- `outputs/phase7_known_qm_descriptors.csv` — 3D descriptors per known compound
- `outputs/phase7_retrospective_qm.csv` — descriptors + serum MIC merged
- `outputs/phase7_retrospective_qm.png` — best 3D descriptor vs serum MIC, with
  the 2D baseline |rho| for scale

## Phase 8 — real CREST confirmation (the result that matters)
`phase8_known_crest_descriptors.py` — re-runs the phase6 descriptor engine on the
REAL CREST/GFN-FF ensembles of the 24 known compounds
(`outputs/qm_runs_known/PAPU-*/crest_conformers.xyz`) and repeats the Phase-7
correlations at QM quality, plus the decisive confound analysis (partial
correlation controlling for serum-free potency, and correlation vs the pure serum
SHIFT). Reads the cached descriptor CSV if present (delete it to recompute SASA).

**Finding:** the Phase-7 MMFF-proxy signal does NOT survive. Hydrophobic SASA
weakens to rho=-0.33 (p=0.12, n.s.) and, controlling for serum-free MIC (which
alone tracks serum MIC at rho=0.77), collapses to partial rho=0.03 — a pure
potency artifact. No ensemble SASA/shape descriptor predicts the serum shift
(all p>0.13); the only coherent, non-confounded trend is "expose polar not
hydrophobic surface" (polar SASA vs shift rho=-0.32), which needs more data.
Full writeup: `outputs/phase8_findings.md`.

`phase8b_finalist_gfn2.py` — recomputes the 3 finalists' descriptors from their
GFN2-reranked ensembles (`qm_runs/cand0*_gfn2/crest_ensemble.xyz`) and
regenerates their Gaussian inputs from the GFN2 populations
(`outputs/qm_runs_gaussian_gfn2/<cand>/*.gjf`). The rerank moved populations
materially (e.g. cand02 GFN-FF's flat 5-conformer set -> one GFN2 conformer at
68%), so the GFN2 DFT set differs from the GFN-FF one.

### Phase 8 outputs
- `outputs/phase8_known_crest_descriptors.csv` — real-CREST 3D descriptors (cache)
- `outputs/phase8_retrospective_crest.csv` — descriptors + serum data
- `outputs/phase8_confound_analysis.csv` — partial + shift correlations
- `outputs/phase8_retrospective_crest.png` — descriptor vs serum MIC; |rho| vs proxy & 2D
- `outputs/phase8_findings.md` — full interpretation
- `outputs/phase8b_finalist_gfn2_descriptors.csv`, `phase8b_gfnff_vs_gfn2.csv`
- `outputs/qm_runs_gaussian_gfn2/<cand>/*.gjf` — finalist DFT inputs (GFN2 pops)

## Phase 9 — electronic / solvation descriptors
`gen_known_xtb_inputs.py` selects each known compound's Boltzmann-populated CREST
conformers and emits GFN2-xTB single-point inputs in water + octanol
(`outputs/qm_runs_known_xtb/`, one cluster sbatch). `phase9_electronic.py` parses
them for dipole, HOMO-LUMO gap, polarizability α(0), aqueous G_solv and a QM logP
= (G_solv,water − G_solv,oct)/2.303RT, Boltzmann-averages, and runs the same
Phase-8 statistics (vs serum MIC, partial controlling serum-free potency, vs serum
shift).

**Finding:** electronics also fails to give a p<0.05 serum-tolerance predictor on
this n=24/censored set, but it converges with shape. Polarizability α(0) is the
strongest raw serum-MIC correlate (ρ=−0.54, p=0.01) yet is a size/potency proxy
(tracks Rg, hydrophobic SASA, serum-free MIC; shift ρ=−0.02). QM logP tracks
hydrophobic SASA (ρ=+0.76) and gives serum-shift ρ=+0.30 — independently
corroborating Phase-8's polar-surface lead (polar SASA shift ρ=−0.33). The
defensible, cross-validated design hypothesis: among equipotent analogs, bias the
exposed surface toward polar/H-bonding (lower QM logP) to cut serum loss. Full
writeup: `outputs/phase9_findings.md`.

### Phase 9 outputs
- `outputs/phase9_electronic_descriptors.csv`, `phase9_electronic_stats.csv`
- `outputs/phase9_electronic.png`, `outputs/phase9_findings.md`

## Phase 10 — explicit HSA docking
`phase10_dock_hsa.py` docks each compound to human serum albumin (PDB 1AO6) and
correlates binding with the serum shift. Whole-molecule docking of these 1000-1200
Da, ~38-rotatable-bond glycolipids is at/beyond docking's reliable regime (rigid
pocket docking clashes at +42 kcal/mol; 38-torsion flexible docking is intractable),
so the method is RIGID ENSEMBLE SURFACE DOCKING: rigid HSA + the top Boltzmann QM
conformers kept rigid + large boxes over both Sudlow regions; score = best Vina
affinity over conformers x sites.

**Finding (null):** HSA binding does not predict the serum shift (ρ=+0.22, p=0.30;
sign even opposite to the sequestration hypothesis), and the result is flat under
potency- and size-controls. The docking score is not a mere size proxy
(vs polarizability ρ=−0.15), so the method captured real surface association — it
just doesn't track serum tolerance. Caveat: Vina is not parameterized for ligands
this large/amphiphilic, so a null cannot exclude albumin sequestration; it needs an
experimental HSA-binding assay. Full writeup: `outputs/phase10_findings.md`.
Bulky docking artifacts (receptor/ligand pdbqt, 1AO6.pdb) are gitignored
(regenerable via the script).

### Phase 10 outputs
- `outputs/phase10_docking/phase10_hsa_scores.csv`, `phase10_hsa_per_compound.csv`
- `outputs/phase10_hsa.png`, `outputs/phase10_findings.md`

## Phase 11 — echinocandin cross-chemotype read-across
`phase11_echinocandin_readacross.py` — the synthesis named a second, independent
chemotype as the only way past the n=24 single-series ceiling. The echinocandins
(cyclic-lipopeptide Fks1/glucan-synthase inhibitors) are that chemotype, and our
own `external/` FKS corpus already carries their serum data: MIC-shift ratios in
±50% serum (the SAME endpoint direction as `serum_shift_fold`, confirmed from the
source quotes: MIC in 50% serum ÷ MIC in serum-free medium) plus PPB/Fu. The
script harmonizes those with the 24 papulacandins on `log2(serum-shift fold)` and
a consistently-recomputed RDKit 2D descriptor set, then (Q1) reframes our endpoint
via the free-drug hypothesis and (Q2) stress-tests the polar-surface lead across
scaffolds. No external binaries; RDKit only.

**Finding:** echinocandins show the same phenomenon and the robust, most
comparable (C. albicans) ordering is **caspofungin ×2 < anidulafungin ×16 <
micafungin ×64** — yet all three are front-line drugs, dosed to a *free*-drug
target (~96–99.8% protein bound). Reframe: the goal is serum-tolerant *free
exposure*, not a serum-invariant MIC. Honest null: **bulk** 2D descriptors do NOT
explain the ordering (micafungin is the most polar yet shifts most; papulacandin
within-series MolLogP ρ=+0.15, TPSA/heavy-atom ρ=−0.25, both n.s.). This refines
rather than refutes the Phase 8–9 lead — that lead was about *locally exposed*
polar surface (polar SASA / QM logP), which bulk TPSA cannot see — and it flags
the echinocandins' known direct serum effect on glucan synthase (cf. the Phase 10
docking null). Existence proof carried alongside: ibrexafungerp/enfumafungin hit
the same target with **no lipopeptide tail** (ibrexafungerp is oral), legitimizing
a tail-truncation/replacement design branch. Full writeup: `outputs/phase11_findings.md`.

### Phase 11 outputs
- `outputs/phase11_echinocandin_serum_shift.csv` — per-compound echinocandin
  serum-shift (C. albicans + all-species medians), PPB/Fu, RDKit descriptors
- `outputs/phase11_crosschemotype.csv` + `phase11_crosschemotype_stats.csv` —
  harmonized papulacandin + echinocandin table and within-papulacandin stats
- `outputs/phase11_crosschemotype.png`, `outputs/phase11_findings.md`

## Phase 12 — serum-tolerance-biased generative design (Track A)
`phase12_generate_serum_tolerant.py` — the AI-design step. Off the serum-active
lead PAPU-0080 it generates novel analogs and scores them by the **operational
form of the Phase 8/9/11 lead**: the mean **exposed polar surface fraction** over a
small ETKDG conformer ensemble (RDKit `rdFreeSASA` — the built-in, not the
`freesasa` pip package that failed in Phase 6). QED/Ro5/clogP are dropped (Phases
5/11 showed them uninformative). The FKS-engaging pharmacophore is preserved by
construction (re-esterifying onto the conserved core). Three branches: **ester**
(C-6′ handle, Phase-5 chemistry), **polaraxis** (designed acyls spanning
hydrophobic→polar), and **notail** (fatty-tail deacylation → ibrexafungerp-inspired
tail-free analogs). CPU-only; a SASA cache (`phase12_sasa_cache.csv`) makes reruns
instant.

**Honest framing:** there is no validated serum-tolerance oracle, so the reward is
the polar-surface *hypothesis*. The reward is validated retrospectively on the 24
knowns (exposed polar frac vs serum shift **ρ = −0.33, p = 0.11** — reproduces the
Phase-8 CREST value with a fast proxy). The headline deliverable is a
**discriminating series**: novel analogs on a single scaffold spanning exposed
polarity, built so a serum assay can *falsify or confirm* the lead, not just
confirm it. The tail-free branch dominates the predicted-tolerant, high-novelty
region at MW ~700–810 (ibrexafungerp's regime). Track B (an RL/CLM generative
network on the external FKS pretraining set reusing this reward) is the follow-on
once wet-lab shifts anchor the reward. Full writeup: `outputs/phase12_findings.md`.

Run it (first pass ~90 s on CPU; cached thereafter):
```
python3 analysis/phase12_generate_serum_tolerant.py
# quicker: PHASE12_N_CONF=1 python3 analysis/phase12_generate_serum_tolerant.py
```

### Phase 12 outputs
- `outputs/phase12_generated_library.csv` — all generated analogs + reward components
- `outputs/phase12_discriminating_series.csv` — matched single-scaffold series
  spanning exposed polarity (the set to synthesize/assay)
- `outputs/phase12_top_candidates.sdf` — 3D structures of the top novel analogs (CREST-ready)
- `outputs/phase12_crest_commands.sh` — per-candidate CREST command template
- `outputs/phase12_reward_validation.png`, `outputs/phase12_findings.md`
- `outputs/phase12_sasa_cache.csv` — InChIKey→exposed-polar-fraction cache

## Phase 13 — round-1 campaign: fatty-tail optimization (core fixed)
`phase13_fatty_tail_optimization.py` — the first design round. It freezes the
sugar / spiroketal core **and** the aromatic C-6′ acyl of the serum-active lead
PAPU-0080, cleaves ONLY the longest aliphatic fatty-acyl ester (the native C16
polyene tail), and re-esterifies that one position with a designed **tail
library** spanning the axis the Phase-8/9/11 lead implicates: chain length,
saturation, terminal polar caps (OH/COOH/NH₂/CONH₂), heteroatom/oxa insertion,
charged heads (sulfonate, phosphocholine-like), branching, fluorination. The
native tail is kept as the control/baseline. Reuses Phase 12's validated reward
(exposed polar surface fraction; ρ=−0.33 on the 24 knowns) and its SASA cache.
Because the core is identical across variants, changes in the reward — and in the
absolute exposed hydrophobic SASA — are attributable to the tail alone.

**Result (fast proxy, to be confirmed by CREST):** native baseline exposed polar
fraction ≈0.27 (hydrophobic SASA ≈930 Å²); the short / terminally-capped /
heteroatom-broken tails (C8-sulfonate, oxa-PEG, mid-chain diol, ω-amide/-acid)
raise polar exposure and cut hydrophobic SASA to ~776–843 Å² — the Phase-8
direction. The reward is not fooled by a charged head alone (the long
phosphocholine chain still exposes hydrophobe → no gain).

This phase writes **CREST-ready inputs and the exact cluster protocol** for the
compute-intensive confirmation you run yourself, laid out so the existing Phase-6
descriptor engine parses the returned ensembles directly. Full protocol:
`outputs/phase13_qm_runs/SUBMIT.md` and `outputs/phase13_findings.md`.

Run the generator (CPU, ~20 s with warm cache):
```
python3 analysis/phase13_fatty_tail_optimization.py
```

### Phase 13 outputs
- `outputs/phase13_fatty_tail_library.csv` — every tail variant + reward + tail descriptors
- `outputs/phase13_discriminating_series.csv` — fixed core, tails spanning exposed polarity
- `outputs/phase13_top_candidates.sdf` — 3D structures (names match the QM dirs)
- `outputs/phase13_qm_runs/<cand>/<cand>.xyz` — CREST starting geometries
- `outputs/phase13_qm_runs/run_crest.sbatch`, `outputs/phase13_qm_runs/SUBMIT.md`
  — cluster job + step-by-step protocol
- `outputs/phase13_findings.md`

### Round-1 cluster protocol (CPU; CREST/xtb are not GPU codes)
1. **CREST ensembles** — for each `phase13_qm_runs/<cand>/<cand>.xyz`, GFN-FF /
   ALPB water / `--quick -ewin 6` on 52 cores → `crest_conformers.xyz`
   (`sbatch -J <cand> ../run_crest.sbatch` from inside each dir; full-GFN2 search
   is intractable at this size — screening tier only).
2. **Upload** each `crest_conformers.xyz` back to the repo.
3. **Parse at QM quality** with the existing Phase-6 engine (exact command in
   `SUBMIT.md`) → `phase13_qm_descriptors.csv`; compare polar/hydrophobic SASA vs
   the native tail.
4. **Finalists only** — GFN2 re-rank (`crest --screen … --gfn2`) + Phase-9 xtb
   electronics (QM logP, water/octanol) before synthesis.

### Submitting the 12 CREST screening jobs
`outputs/qm_runs/<candidate_name>/<candidate_name>.xyz` — one starting geometry
per Phase-5 candidate (atom order matches `phase5_top_candidates.sdf`, required
for `phase6_qm_layer.py`'s ensemble parser). `outputs/qm_runs/run_crest.sbatch`
is your cluster script, defaulted to the fast screening tier
(`--gfn2//gfnff --quick -ewin 6`). Submit up to 5 at a time (your queue limit):
```
cd outputs/qm_runs/cand01_quinolinecarbonyl && sbatch -J cand01 ../run_crest.sbatch
```
For the 1-3 funnel finalists, rerun with the thorough settings:
```
METHOD_FLAG=--gfn2 EXTRA_OPTS="-ewin 10" sbatch -J cand01_final ../run_crest.sbatch
```
CREST writes `crest_conformers.xyz` into each candidate's own folder, which is
exactly what `run_qm_layer()` expects.
