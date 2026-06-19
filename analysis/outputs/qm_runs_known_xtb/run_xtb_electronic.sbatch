#!/bin/bash
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
