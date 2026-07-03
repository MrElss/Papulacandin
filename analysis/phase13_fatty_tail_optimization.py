#!/usr/bin/env python3
"""
phase13_fatty_tail_optimization.py
==================================
PHASE 13 — Round-1 campaign: optimize the LONG-CHAIN FATTY ACID moiety only,
holding the sugar / spiroketal core (and the aromatic C-6' handle) fixed.

Rationale
---------
The fatty tail is the serum-liability driver: Phase 8 found more EXPOSED
HYDROPHOBIC surface -> larger serum shift (rho +0.24) and more EXPOSED POLAR
surface -> smaller shift (rho -0.33); the tail is the dominant hydrophobic-surface
contributor. Phase 11 (echinocandins) independently shows the acyl tail governs
protein binding. Freezing the core (a) preserves the FKS-engaging pharmacophore so
intrinsic potency is roughly held constant (Phase 8: serum MIC = potency x
tolerance), and (b) removes the core as a confounder, so any change in the reward
is attributable to the tail. This is the cleanest possible one-variable first round.

What it does
------------
1. On a serum-active lead (PAPU-0080), locate the longest genuine ALIPHATIC
   fatty-acyl ester (the C16 polyene tail) and cleave ONLY that ester — the sugar,
   spiroketal, and the aromatic C-6' acyl are untouched.
2. Re-esterify that single position with a designed TAIL LIBRARY spanning the
   axis the Phase-8/9/11 lead implicates: chain length, saturation, terminal polar
   caps (OH/COOH/NH2/CONH2), heteroatom (oxa/PEG) insertion, charged heads
   (sulfonate, phosphocholine-like), branching, fluorination. The native tail is
   kept as a control.
3. Score each analog by the fast 3D reward from Phase 12 (mean exposed polar
   surface fraction over an ETKDG ensemble, RDKit rdFreeSASA) — reused, so the
   reward and its retrospective validation (rho=-0.33 on the 24 knowns) carry over.
4. Emit a DISCRIMINATING SERIES (one core, tail spanning the axis) and, for the
   compute-intensive confirmation you will run on the cluster, CREST-ready inputs
   laid out so the EXISTING phase6/phase8 descriptor engine parses the returned
   ensembles directly. The exact cluster protocol is written to
   phase13_findings.md and phase13_qm_runs/SUBMIT.md.

Honest status: the reward is the polar-surface HYPOTHESIS (no validated oracle).
Round 1's job is to produce a tail series whose serum shift, measured in vitro,
tests it. QM ensembles only confirm the exposed-surface property the fast proxy
predicts; they are not themselves the serum endpoint.

Outputs (analysis/outputs/)
---------------------------
* phase13_fatty_tail_library.csv       — every tail variant + reward + tail descriptors
* phase13_discriminating_series.csv    — matched core, tail spanning exposed polarity
* phase13_top_candidates.sdf           — 3D structures (names match the QM dirs)
* phase13_qm_runs/<cand>/<cand>.xyz    — CREST starting geometries (+ run_crest.sbatch, SUBMIT.md)
* phase13_findings.md                  — interpretation + the EXACT cluster protocol
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phase12_generate_serum_tolerant as p12  # noqa: E402  (guarded __main__)

ROOT = p12.ROOT
CORE = p12.CORE
OUT = p12.OUT
QM_DIR = os.path.join(OUT, "phase13_qm_runs")
TEMPLATE = p12.SERUM_TOLERANT_TEMPLATE  # PAPU-0080, a serum-active lead
TOP_N_FOR_QM = 12

# The tail library: designed fatty-acyl replacements (attach as the ester acyl
# via [*:1]C(=O)R) spanning the lipophilic -> polar axis. Native kept as control.
TAIL_LIBRARY = {
    # --- length / saturation (modulate exposed hydrophobic AREA) ---
    "native_C16_polyene": None,  # filled at runtime from the template
    "C16_saturated": "[*:1]C(=O)CCCCCCCCCCCCCCC",
    "C12_saturated": "[*:1]C(=O)CCCCCCCCCCC",
    "C8_saturated": "[*:1]C(=O)CCCCCCC",
    # --- terminal polar cap on a mid-length chain (bury less hydrophobe) ---
    "C12_omega_OH": "[*:1]C(=O)CCCCCCCCCCCO",
    "C12_omega_COOH": "[*:1]C(=O)CCCCCCCCCCC(=O)O",
    "C12_omega_NH2": "[*:1]C(=O)CCCCCCCCCCCN",
    "C12_omega_CONH2": "[*:1]C(=O)CCCCCCCCCCC(N)=O",
    # --- heteroatom / oxa insertion (break up the hydrophobic stretch) ---
    "oxa_PEG3": "[*:1]C(=O)COCCOCCOCC",
    "amide_split_C12": "[*:1]C(=O)CCCCC(=O)NCCCCCC",
    # --- charged / lipid-mimetic heads (aggressive; SA set aside) ---
    "C8_omega_sulfonate": "[*:1]C(=O)CCCCCCCS(=O)(=O)O",
    "phosphocholine_like": "[*:1]C(=O)CCCCCCCCOP(=O)([O-])OCC[N+](C)(C)C",
    # --- shape modulation at ~constant lipophilicity ---
    "C12_branched": "[*:1]C(=O)CCCCCC(CC)CCC",
    "C12_omega_CF3": "[*:1]C(=O)CCCCCCCCC(F)(F)F",
    "mid_chain_diol_C12": "[*:1]C(=O)CCCC(O)CC(O)CCCC",
}


def cleave_longest_fatty_tail(mol):
    """Cleave ONLY the longest genuine aliphatic fatty-acyl ester.

    Returns (core_with_dummy, native_acyl_with_dummy) or (None, None). The acyl
    fragment is SANITIZED before the aromatic-ring test so an aromatic acyl
    (the C-6' handle) can never be mistaken for the aliphatic fatty tail.
    """
    best = None  # (n_carbons, oxy, carb)
    for match in mol.GetSubstructMatches(p12.ALIPH_ESTER):
        oxy, carb = match[0], match[1]
        rw = Chem.RWMol(mol)
        b = rw.GetBondBetweenAtoms(oxy, carb)
        if b is None:
            continue
        frag = Chem.FragmentOnBonds(rw, [b.GetIdx()], addDummies=True)
        parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
        if len(parts) != 2:
            continue
        acyl = min(parts, key=lambda p: p.GetNumAtoms())
        n_c = sum(1 for a in acyl.GetAtoms() if a.GetAtomicNum() == 6)
        if rdMolDescriptors.CalcNumAromaticRings(acyl) == 0 and n_c >= 8:
            if best is None or n_c > best[0]:
                best = (n_c, oxy, carb)
    if best is None:
        return None, None
    _, oxy, carb = best
    rw = Chem.RWMol(mol)
    frag = Chem.FragmentOnBonds(rw, [rw.GetBondBetweenAtoms(oxy, carb).GetIdx()],
                                addDummies=True)
    parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
    core = max(parts, key=lambda p: p.GetNumAtoms())
    acyl = min(parts, key=lambda p: p.GetNumAtoms())
    p12._tag_dummies(core)
    p12._tag_dummies(acyl)
    return core, acyl


def tail_descriptors(acyl_mol):
    """Descriptors of the tail fragment itself (dummy ignored)."""
    n_c = sum(1 for a in acyl_mol.GetAtoms() if a.GetAtomicNum() == 6)
    n_o_n = sum(1 for a in acyl_mol.GetAtoms() if a.GetAtomicNum() in (7, 8))
    n_db = sum(1 for b in acyl_mol.GetBonds()
               if b.GetBondType() == Chem.BondType.DOUBLE
               and b.GetBeginAtom().GetAtomicNum() == 6
               and b.GetEndAtom().GetAtomicNum() == 6)
    return dict(tail_n_carbon=n_c, tail_n_O_N=n_o_n, tail_CC_double_bonds=n_db,
                tail_clogp=round(Crippen.MolLogP(acyl_mol), 2),
                tail_tpsa=round(Descriptors.TPSA(acyl_mol), 1))


def build_tail_acyls(native_acyl):
    acyls = {}
    for name, smi in TAIL_LIBRARY.items():
        if name == "native_C16_polyene":
            acyls[name] = native_acyl
            continue
        am = Chem.MolFromSmiles(smi)
        if am:
            p12._tag_dummies(am)
            acyls[name] = am
    return acyls


def _fatty_tail_series(novel, n_bins=3, per_bin=2):
    """Pick novel tails spanning LOW/MID/HIGH exposed polarity on the fixed core,
    so the serum assay can regress serum shift on tail polarity directly."""
    if novel.empty:
        return novel
    n = novel.copy()
    try:
        n["polar_bin"] = pd.qcut(n["exposed_polar_fraction"], q=n_bins,
                                 labels=["low", "mid", "high"], duplicates="drop")
    except ValueError:
        n["polar_bin"] = "all"
    picks = [grp.sort_values("novelty", ascending=False).head(per_bin)
             for _, grp in n.groupby("polar_bin", observed=True)]
    out = pd.concat(picks).sort_values("exposed_polar_fraction")
    return out[["polar_bin", "tail_name", "exposed_polar_fraction",
                "exposed_hydrophobic_sasa", "tail_n_carbon", "tail_clogp",
                "novelty", "mw", "inchikey", "smiles"]]


def main():
    cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
    p12._load_sasa_cache()

    known_fps, known_ik = [], {}
    for smi, cid in zip(cm["smiles_canonical"], cm["compound_id"]):
        if isinstance(smi, str):
            m = Chem.MolFromSmiles(smi)
            if m:
                known_fps.append(p12.ecfp(m))
                known_ik[Chem.MolToInchiKey(m)] = cid

    tmpl = Chem.MolFromSmiles(
        cm.set_index("compound_id").loc[TEMPLATE, "smiles_canonical"])
    core, native_acyl = cleave_longest_fatty_tail(tmpl)
    assert core is not None, "no aliphatic fatty tail found on the template"
    native_smiles = Chem.MolToSmiles(native_acyl)

    acyls = build_tail_acyls(native_acyl)
    rows, mol_by_key, seen = [], {}, set()
    for tname, acyl in acyls.items():
        full = p12.attach(core, acyl)
        if full is None:
            continue
        ik = Chem.MolToInchiKey(full)
        if ik in seen:
            continue
        seen.add(ik)
        frac, total_sasa, n_ok = p12.cached_exposed_polar_fraction(full)
        if n_ok == 0:
            continue
        from rdkit.Chem import DataStructs
        sims = DataStructs.BulkTanimotoSimilarity(p12.ecfp(full), known_fps)
        rows.append(dict(
            branch="fatty_tail",
            tail_name=tname,
            is_native=(tname == "native_C16_polyene"),
            exposed_polar_fraction=round(frac, 4),
            exposed_hydrophobic_fraction=round(1 - frac, 4),
            exposed_hydrophobic_sasa=round(total_sasa * (1 - frac), 1),
            total_sasa=round(total_sasa, 1),
            n_conf_ok=n_ok,
            novelty=round(1 - max(sims), 3) if sims else 1.0,
            mw=round(Descriptors.MolWt(full), 1),
            clogp=round(Descriptors.MolLogP(full), 2),
            inchikey=ik,
            matches_known=known_ik.get(ik, ""),
            smiles=Chem.MolToSmiles(full),
            **tail_descriptors(acyl),
        ))
        mol_by_key[ik] = full
    p12._save_sasa_cache()

    lib = pd.DataFrame(rows)
    lib["serum_tolerance_reward"] = lib["exposed_polar_fraction"]
    lib["is_novel"] = lib["matches_known"] == ""
    # Rank by reward; the native tail is the baseline to beat.
    lib = lib.sort_values("serum_tolerance_reward", ascending=False).reset_index(drop=True)
    lib.to_csv(os.path.join(OUT, "phase13_fatty_tail_library.csv"), index=False)

    native = lib[lib["is_native"]].iloc[0]
    novel = lib[lib["is_novel"]].copy()
    series = _fatty_tail_series(novel)
    series.to_csv(os.path.join(OUT, "phase13_discriminating_series.csv"), index=False)

    top = novel.head(TOP_N_FOR_QM)
    _write_qm_inputs(top, mol_by_key)
    _write_findings(lib, native, series, native_smiles)

    print("PHASE 13 — fatty-tail optimization (core fixed)")
    print("=" * 58)
    print(f"Template: {TEMPLATE}   native tail: {native_smiles}")
    print(f"Tail variants scored: {len(lib)}  (novel: {int(lib['is_novel'].sum())})")
    print(f"Native tail exposed polar fraction (baseline): "
          f"{native['exposed_polar_fraction']:.3f}  "
          f"hydrophobic SASA {native['exposed_hydrophobic_sasa']:.0f}")
    print("\nTop tails by reward (higher exposed polar fraction = predicted more tolerant):")
    cols = ["tail_name", "exposed_polar_fraction", "exposed_hydrophobic_sasa",
            "tail_n_carbon", "novelty", "mw"]
    print(top[cols].to_string(index=False))
    print(f"\nCREST inputs for {len(top)} candidates -> {os.path.relpath(QM_DIR, ROOT)}/")
    print("Wrote: phase13_fatty_tail_library.csv, phase13_discriminating_series.csv,")
    print("       phase13_top_candidates.sdf, phase13_qm_runs/, phase13_findings.md")


RUN_CREST_SBATCH = """#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=52
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J p13_crest
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Phase 13 fatty-tail candidates — GFN-FF screening tier (same settings as the
# validated Phase 6 funnel). Run from inside each candidate directory:
#   cd phase13_qm_runs/t01_<tail> && sbatch -J t01 ../run_crest.sbatch
# Finalists only, re-rank at GFN2 (cheap, no re-search):
#   crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 6 --T 52
set -euo pipefail
ulimit -s unlimited
APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-52}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
cd "${SLURM_SUBMIT_DIR:-$PWD}"
INP=$(ls -1 *.xyz 2>/dev/null | head -n1)
[[ -n "$INP" ]] || { echo "[ERR] no .xyz here"; exit 2; }
METHOD_FLAG=${METHOD_FLAG:---gfnff}
SOLV_FLAG=${SOLV_FLAG:---alpb water}
EXTRA_OPTS=${EXTRA_OPTS:---quick -ewin 6}
set -x
crest "$INP" $METHOD_FLAG $SOLV_FLAG $EXTRA_OPTS --T "$OMP_NUM_THREADS" \\
  --chrg 0 --uhf 0 2>&1 | tee run.log
"""


def _write_qm_inputs(top, mol_by_key):
    os.makedirs(QM_DIR, exist_ok=True)
    writer = Chem.SDWriter(os.path.join(OUT, "phase13_top_candidates.sdf"))
    names = []
    for i, (_, r) in enumerate(top.iterrows(), 1):
        name = f"t{i:02d}_{r['tail_name']}"
        m = p12._embed(mol_by_key[r["inchikey"]])
        if m is None:
            continue
        m.SetProp("_Name", name)
        m.SetProp("serum_tolerance_reward", str(r["serum_tolerance_reward"]))
        m.SetProp("exposed_polar_fraction", str(r["exposed_polar_fraction"]))
        writer.write(m)
        d = os.path.join(QM_DIR, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f"{name}.xyz"), "w") as fh:
            fh.write(Chem.MolToXYZBlock(m))
        names.append(name)
    writer.close()
    with open(os.path.join(QM_DIR, "run_crest.sbatch"), "w") as fh:
        fh.write(RUN_CREST_SBATCH)
    with open(os.path.join(QM_DIR, "SUBMIT.md"), "w") as fh:
        fh.write(_submit_md(names))


def _submit_md(names):
    submit = "\n".join(
        f"cd {n} && sbatch -J {n} ../run_crest.sbatch && cd .." for n in names[:5])
    return f"""# Phase 13 — cluster protocol (CPU; CREST/xTB are not GPU codes)

Candidates: {len(names)}. One directory each, with a starting geometry
`<name>.xyz`. CREST does the conformer search (the compute-intensive step);
the returned ensembles are parsed by the EXISTING Phase-6 engine.

## Step 1 — CREST conformer ensembles (GFN-FF screening tier)
Submit up to your queue limit at a time (each writes `crest_conformers.xyz`):
```
cd phase13_qm_runs
{submit}
# ...repeat for the rest
```
Cost lever: these ~140-atom, 30+ rotatable-bond glycolipids are cheap at
GFN-FF, intractable at full GFN2 search (~10 days/compound). Keep the search at
GFN-FF; refine only finalists.

## Step 2 — upload results (which CREST files to keep)
The Phase-6 parser reads exactly ONE file per candidate: `crest_conformers.xyz`
(conformer energies are in its comment lines). Keep, per directory:
- **`crest_conformers.xyz`** — REQUIRED (the only parser input).
- `crest_best.xyz`, `run.log`, `crest.energies` — small; provenance / QC / finalist re-rank.
Discard the bulky regenerable scratch (already in `.gitignore`):
`crest_dynamics.trj`, `confcross.xyz`, `crest_rotamers.xyz`, `gfnff_topo`,
`crest.restart`, `crestopt.log`, `*.xtbrestart`, `wbo`, the slurm `*.<jobid>.out`.
One-liner to delete only the known bulky scratch (explicit denylist — safe):
```
find phase13_qm_runs -type f \\( -name 'crest_dynamics.trj' -o -name 'confcross.xyz' \\
  -o -name 'crest_rotamers.xyz' -o -name 'gfnff_topo' -o -name 'crest.restart' \\
  -o -name 'crestopt.log' -o -name '*.xtbrestart' -o -name '*.xtbtopo.mol' \\
  -o -name 'wbo' -o -name '*.[0-9]*.out' -o -name '*.[0-9]*.err' \\) -delete
```
Then commit each `phase13_qm_runs/<name>/crest_conformers.xyz` back to the repo.

## Step 3 — QM-quality exposed-surface descriptors (parse; CPU-cheap, local)
Re-scores the tails at real-ensemble quality with the SAME descriptor engine
the project already validated (Phase 6/8):
```
python3 -c "import sys; sys.path.insert(0,'analysis'); import phase6_qm_layer as p6; \\
  p6.run_qm_layer('analysis/outputs/phase13_top_candidates.sdf', \\
                  'analysis/outputs/phase13_qm_runs', \\
                  'analysis/outputs/phase13_qm_descriptors.csv', \\
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
"""


def _write_findings(lib, native, series, native_smiles):
    top = lib[lib["is_novel"]].head(6)
    beats = lib[(lib["is_novel"]) &
                (lib["exposed_polar_fraction"] > native["exposed_polar_fraction"])]
    top_lines = "\n".join(
        f"- **{r['tail_name']}** — exposed polar frac {r['exposed_polar_fraction']:.2f} "
        f"(native {native['exposed_polar_fraction']:.2f}), exposed hydrophobic SASA "
        f"{r['exposed_hydrophobic_sasa']:.0f}, C{r['tail_n_carbon']}, novelty {r['novelty']:.2f}"
        for _, r in top.iterrows())
    md = f"""# Phase 13 — fatty-tail optimization, round 1 (findings + protocol)

## Scope
Core (sugar + spiroketal + aromatic C-6' acyl) FROZEN; only the long-chain fatty
acid is varied. Template **{TEMPLATE}**; native tail cleaved:
`{native_smiles}` (a C{int(native['tail_n_carbon'])} polyene). {len(lib)} tail
variants scored (native = control), all sharing one identical core.

## Reward (reused from Phase 12, validated rho=-0.33 on 24 knowns)
Mean exposed POLAR surface fraction over an ETKDG ensemble (rdFreeSASA). Because
the core is identical across all variants, differences in this readout — and in
the absolute exposed HYDROPHOBIC SASA — are attributable to the tail. Goal:
raise exposed polar fraction / lower exposed hydrophobic SASA vs the native tail
(the Phase-8 direction), WITHOUT touching the pharmacophore.

## Result (fast proxy — to be confirmed by CREST on your cluster)
Native baseline: exposed polar fraction **{native['exposed_polar_fraction']:.2f}**,
exposed hydrophobic SASA **{native['exposed_hydrophobic_sasa']:.0f} A^2**.
{len(beats)} novel tails beat the native baseline on the reward. Top:

{top_lines}

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
"""
    with open(os.path.join(OUT, "phase13_findings.md"), "w", encoding="utf-8") as fh:
        fh.write(md)


if __name__ == "__main__":
    main()
