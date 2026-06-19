#!/usr/bin/env python3
"""
gen_known_xtb_inputs.py  (Phase 9 — electronic descriptors, input stage)
========================================================================
Phase 8 showed ensemble SHAPE/SASA does not predict serum tolerance once
intrinsic potency is controlled; the one coherent (if weak) lead was "expose
polar, not hydrophobic, surface" (polar SASA vs serum shift rho=-0.33). That
points at ELECTRONICS / solvation, the one axis the shape descriptors miss.

Phase 9 tests it with cheap GFN2-xTB single points (semi-empirical QM; we already
run GFN2 on the cluster) on each known compound's Boltzmann-populated conformers,
in TWO implicit solvents:
  * water    -> dipole, G_solv(water), HOMO-LUMO gap, atomic charges
  * octanol  -> G_solv(octanol)
The pair gives a principled, QM-level hydrophobicity:
    logP_xtb  ~  (G_solv(octanol) - G_solv(water)) / (2.303 R T)
i.e. preference for the nonpolar phase -- a far better hydrophobicity measure
than the 2D clogP that failed in Phase 1, and the natural variable for an
albumin-sequestration mechanism.

This script picks, per compound, the top Boltzmann conformers covering >=50%
population (cap MAX_CONFS) from the real CREST ensemble, writes each as a plain
xyz under outputs/qm_runs_known_xtb/<cid>/, plus a per-compound weights.csv, and
emits ONE cluster sbatch (run_xtb_electronic.sbatch) that loops every compound x
conformer x {water, octanol} and runs the GFN2 single points.

Phase 9 parse/analysis (after the runs return) lives in phase9_electronic.py.
"""
from __future__ import annotations
import os
import glob
import numpy as np
import pandas as pd
import phase6_qm_layer as p6

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
KNOWN_ROOT = os.path.join(OUT, "qm_runs_known")
DEST = os.path.join(OUT, "qm_runs_known_xtb")
POP_CUTOFF = 0.50
MAX_CONFS = 5
SOLVENTS = ["water", "octanol"]


def select_conformers(frames, weights):
    order = np.argsort(weights)[::-1]
    cum, chosen = 0.0, []
    for idx in order:
        chosen.append(int(idx))
        cum += weights[idx]
        if cum >= POP_CUTOFF or len(chosen) >= MAX_CONFS:
            break
    return chosen


def write_xyz(frame, path, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(frame['elements'])}\n{comment}\n")
        for el, xyz in zip(frame["elements"], frame["coords"]):
            fh.write(f"{el:2s} {xyz[0]:14.8f} {xyz[1]:14.8f} {xyz[2]:14.8f}\n")


def main():
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    os.makedirs(DEST, exist_ok=True)
    made = 0
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        folders = glob.glob(os.path.join(KNOWN_ROOT, f"{cid}_*"))
        if not folders:
            print(f"[skip] no CREST folder for {cid}")
            continue
        xyz = os.path.join(folders[0], "crest_conformers.xyz")
        if not os.path.exists(xyz):
            print(f"[skip] no ensemble for {cid}")
            continue
        frames = p6.parse_crest_ensemble(xyz)
        weights, _ = p6.boltzmann_weights([f["energy_hartree"] for f in frames])
        chosen = select_conformers(frames, weights)
        cdir = os.path.join(DEST, cid)
        os.makedirs(cdir, exist_ok=True)
        rows = []
        for rank, idx in enumerate(chosen, 1):
            fn = f"conf{rank:02d}.xyz"
            write_xyz(frames[idx], os.path.join(cdir, fn),
                      f"{cid} conf{rank} pop={weights[idx]:.4f}")
            rows.append(dict(conformer=f"conf{rank:02d}", pop=round(float(weights[idx]), 4)))
        pd.DataFrame(rows).to_csv(os.path.join(cdir, "weights.csv"), index=False)
        made += 1
        print(f"[ok] {cid}: {len(chosen)} conformers (cum pop "
              f"{sum(weights[i] for i in chosen):.2f})")
    write_sbatch()
    print(f"\nWrote xtb inputs for {made} compounds under "
          f"{os.path.relpath(DEST, HERE)} + run_xtb_electronic.sbatch")


def write_sbatch():
    sb = os.path.join(DEST, "run_xtb_electronic.sbatch")
    with open(sb, "w") as fh:
        fh.write(r"""#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH -J xtb_elec
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#
# Phase 9 electronic descriptors: GFN2-xTB single points on the Boltzmann-
# populated conformers of every known compound, in water and octanol.
# Each run writes <conf>_<solvent>.out (parsed by phase9_electronic.py for
# dipole, G_solv, HOMO-LUMO gap, charges; logP from the solvent pair).
#
# Submit ONCE from this directory:  sbatch run_xtb_electronic.sbatch
set -euo pipefail
ulimit -s unlimited
APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"          # provides xtb on PATH (same module as CREST)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OMP_STACKSIZE=4G
which xtb

ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$ROOT"
for d in */ ; do
  d="${d%/}"
  [[ -f "$d/conf01.xyz" ]] || continue
  for xyz in "$d"/conf*.xyz ; do
    base="${xyz%.xyz}"
    for solv in water octanol ; do
      out="${base}_${solv}.out"
      [[ -s "$out" ]] && { echo "skip $out"; continue; }   # resume-friendly
      echo "[$(date +%H:%M:%S)] xtb $xyz  --gfn2 --alpb $solv"
      # --gfn 2 single point; --alpb implicit solvent; charges/dipole/gap in stdout
      xtb "$xyz" --gfn 2 --alpb "$solv" --chrg 0 --uhf 0 \
          --namespace "${base}_${solv}" > "$out" 2>&1 || echo "  [warn] nonzero exit for $out"
    done
  done
done
echo "[$(date)] Phase 9 xtb electronic single points done."
""")


if __name__ == "__main__":
    main()
