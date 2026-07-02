# Phase 13 — cluster protocol (CPU; CREST/xTB are not GPU codes)

Candidates: 12. One directory each, with a starting geometry
`<name>.xyz`. CREST does the conformer search (the compute-intensive step);
the returned ensembles are parsed by the EXISTING Phase-6 engine.

## Step 1 — CREST conformer ensembles (GFN-FF screening tier)
Submit up to your queue limit at a time (each writes `crest_conformers.xyz`):
```
cd phase13_qm_runs
cd t01_C8_omega_sulfonate && sbatch -J t01_C8_omega_sulfonate ../run_crest.sbatch && cd ..
cd t02_oxa_PEG3 && sbatch -J t02_oxa_PEG3 ../run_crest.sbatch && cd ..
cd t03_mid_chain_diol_C12 && sbatch -J t03_mid_chain_diol_C12 ../run_crest.sbatch && cd ..
cd t04_C12_omega_CONH2 && sbatch -J t04_C12_omega_CONH2 ../run_crest.sbatch && cd ..
cd t05_C12_omega_COOH && sbatch -J t05_C12_omega_COOH ../run_crest.sbatch && cd ..
# ...repeat for the rest
```
Cost lever: these ~140-atom, 30+ rotatable-bond glycolipids are cheap at
GFN-FF, intractable at full GFN2 search (~10 days/compound). Keep the search at
GFN-FF; refine only finalists.

## Step 2 — upload results
Commit each `phase13_qm_runs/<name>/crest_conformers.xyz` back to the repo
(that file is the only thing the parser needs).

## Step 3 — QM-quality exposed-surface descriptors (parse; CPU-cheap, local)
Re-scores the tails at real-ensemble quality with the SAME descriptor engine
the project already validated (Phase 6/8):
```
python3 -c "import sys; sys.path.insert(0,'analysis'); import phase6_qm_layer as p6; \
  p6.run_qm_layer('analysis/outputs/phase13_top_candidates.sdf', \
                  'analysis/outputs/phase13_qm_runs', \
                  'analysis/outputs/phase13_qm_descriptors.csv', \
                  'analysis/outputs/phase13_qm_gaussian', real_run=True)"
```
Compare `hydrophobic_sasa`/`polar_sasa` across tails vs the native baseline.

## Step 4 (finalists only) — GFN2 re-rank + electronics
```
# GFN2 re-rank of an existing ensemble (no re-search):
crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 6 --T 52
# QM logP / dipole / polarizability (Phase 9 machinery), water + octanol:
#   see analysis/gen_known_xtb_inputs.py + analysis/phase9_electronic.py
```
