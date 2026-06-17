#!/usr/bin/env python3
"""
gen_known_crest_inputs.py
=========================
Generate cluster-ready CREST input folders for the 24 KNOWN serum-gap compounds
(Yeung-1996 fusacandin analogs 6a-6u + Fusacandin A/B + Papulacandin B), so the
Phase-7 retrospective (currently an MMFF proxy) can be confirmed with real
CREST/GFN-FF ensembles on the same footing as the 12 novel candidates.

Output mirrors the reference format of analysis/outputs/qm_runs/cand*/:
  outputs/qm_runs_known/<PAPU-ID>_<name>/
      <PAPU-ID>_<name>.xyz   -- single starting geometry (XYZ; atomcount /
                                "<folder> starting geometry for CREST" / atoms)
      run_crest.sbatch       -- per-folder copy of the screening-tier script

The starting geometry is the lowest-energy MMFF94 conformer from a short ETKDGv3
embedding (CREST does the real conformer search; it only needs one sane 3D
input). Atom order is irrelevant for CREST here (these knowns are not parsed by
phase6 against a fixed reference topology), but Hs are added so xtb sees a
complete structure.

Submit on the cluster exactly like the candidates (queue limit 5 at a time):
  cd outputs/qm_runs_known/<folder> && sbatch -J <shortjob> run_crest.sbatch
The sbatch auto-detects the lone *.xyz and runs the fast screening tier
(--gfnff --alpb water --quick -ewin 6).
"""
from __future__ import annotations
import os
import re
import shutil
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
SDF_DIR = os.path.join(HERE, "..", "curated", "structures", "sdf")
REF_SBATCH = os.path.join(OUT, "qm_runs", "run_crest.sbatch")
DEST = os.path.join(OUT, "qm_runs_known")
N_EMBED = 15
RANDOM_SEED = 0xC0FFEE


def sanitize(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")
    return re.sub(r"_+", "_", s)


def find_sdf(compound_id: str):
    import glob
    hits = glob.glob(os.path.join(SDF_DIR, f"{compound_id}_*.sdf"))
    return hits[0] if hits else None


def lowest_energy_geometry(sdf_path: str):
    """Return (mol_with_Hs, conf_id) of the lowest-MMFF conformer, or (None, None)."""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        return None, None
    mol = Chem.AddHs(mol, addCoords=False)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    params.useRandomCoords = True
    params.numThreads = 0
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_EMBED, params=params))
    if not cids:
        params.maxIterations = 2000
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_EMBED, params=params))
    if not cids:
        return None, None
    props = AllChem.MMFFGetMoleculeProperties(mol)
    best_e, best_cid = None, None
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid) if props else None
        if ff is None:
            continue
        ff.Minimize(maxIts=1000)
        e = ff.CalcEnergy()
        if best_e is None or e < best_e:
            best_e, best_cid = e, cid
    return (mol, best_cid) if best_cid is not None else (None, None)


def write_xyz(mol, cid, folder_name, path):
    conf = mol.GetConformer(cid)
    with open(path, "w") as fh:
        fh.write(f"{mol.GetNumAtoms()}\n")
        fh.write(f"{folder_name} starting geometry for CREST\n")
        for atom in mol.GetAtoms():
            p = conf.GetAtomPosition(atom.GetIdx())
            fh.write(f"{atom.GetSymbol():<2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}\n")


def main():
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    os.makedirs(DEST, exist_ok=True)
    made, skipped = 0, []
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        sdf = find_sdf(cid)
        if sdf is None:
            skipped.append((cid, "no SDF"))
            continue
        folder_name = f"{cid}_{sanitize(r['name'])}"
        folder = os.path.join(DEST, folder_name)
        os.makedirs(folder, exist_ok=True)
        mol, conf = lowest_energy_geometry(sdf)
        if mol is None:
            skipped.append((cid, "embed/MMFF failed"))
            continue
        write_xyz(mol, conf, folder_name, os.path.join(folder, f"{folder_name}.xyz"))
        shutil.copy2(REF_SBATCH, os.path.join(folder, "run_crest.sbatch"))
        made += 1
        print(f"[ok] {folder_name}  ({mol.GetNumAtoms()} atoms)")
    print(f"\nWrote {made} CREST input folders under {os.path.relpath(DEST, HERE)}")
    if skipped:
        print("Skipped:")
        for cid, why in skipped:
            print(f"  {cid}: {why}")


if __name__ == "__main__":
    main()
