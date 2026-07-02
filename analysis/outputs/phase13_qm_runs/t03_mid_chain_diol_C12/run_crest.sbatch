#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=25
#SBATCH --mem=64G
#SBATCH -J t03
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
# Phase 13 fatty-tail candidates — GFN-FF screening tier (same settings as the
# validated Phase 6 funnel). Run from inside each candidate directory:
#   cd phase13_qm_runs/t01_<tail> && sbatch -J t01 ../run_crest.sbatch
# Finalists only, re-rank at GFN2 (cheap, no re-search):
#   crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 6 --T 25
set -euo pipefail
ulimit -s unlimited
APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-25}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
cd "${SLURM_SUBMIT_DIR:-$PWD}"
INP=$(ls -1 *.xyz 2>/dev/null | head -n1)
[[ -n "$INP" ]] || { echo "[ERR] no .xyz here"; exit 2; }
METHOD_FLAG=${METHOD_FLAG:---gfnff}
SOLV_FLAG=${SOLV_FLAG:---alpb water}
EXTRA_OPTS=${EXTRA_OPTS:---quick -ewin 6}
set -x
crest "$INP" $METHOD_FLAG $SOLV_FLAG $EXTRA_OPTS --T "$OMP_NUM_THREADS" \
  --chrg 0 --uhf 0 2>&1 | tee run.log
