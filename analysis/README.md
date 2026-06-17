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
