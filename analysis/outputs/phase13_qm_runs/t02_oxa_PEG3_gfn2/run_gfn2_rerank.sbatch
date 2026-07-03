#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=35
#SBATCH --mem=64G
#SBATCH -J t02_oxa_PEG3_gfn2
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#
# Phase 13 finalist GFN2 RE-RANK (t02_oxa_PEG3)
# =======================================
# Re-rank the GFN-FF screening-tier conformer POOL at GFN2-xTB WITHOUT a new
# conformer search ("--screen"). GFN2 ranks these polyol glycolipids' many
# intramolecular H-bonds far better than the GFN-FF force field used for the
# search, so the resulting Boltzmann populations -- which drive the SASA
# descriptors that decided the Step-3 ranking -- are more trustworthy. Cheap
# relative to the original search.
#
# CHARGE: this molecule is NET NEUTRAL as built (t02_oxa_PEG3 formal charge 0), so
# --chrg 0. NOTE for t01: it is the neutral sulfonic ACID (-SO3H). To instead
# model the physiological sulfonATE anion, rebuild the geometry with [O-] and
# set --chrg -1 -- that is a SEPARATE calculation, not this re-rank of the
# neutral ensemble.
#
# Input : pool_gfnff.xyz  (the GFN-FF ensemble, already in this folder)
# Output: crest_ensemble.xyz  (GFN2-reranked; the file to return), tee'd to rerank.log
# Submit from inside this folder:
#   sbatch run_gfn2_rerank.sbatch
set -euo pipefail
ulimit -s unlimited

APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS
echo "[$(date)] node=${SLURM_NODELIST:-NA} T=$OMP_NUM_THREADS jobid=${SLURM_JOB_ID:-NA}"
which crest; which xtb

cd "${SLURM_SUBMIT_DIR:-$PWD}"
POOL=${POOL:-pool_gfnff.xyz}
[[ -f "$POOL" ]] || { echo "[ERR] no conformer pool '$POOL' in $PWD"; exit 2; }

METHOD_FLAG=${METHOD_FLAG:---gfn2}
SOLV_FLAG=${SOLV_FLAG:---alpb water}
EWIN=${EWIN:--ewin 10}          # wide window: keep GFN-FF-misranked conformers for GFN2 to re-judge
CHRG=${CHRG:-0}
UHF=${UHF:-0}

set -x
crest --screen "$POOL" $METHOD_FLAG $SOLV_FLAG $EWIN --chrg "$CHRG" --uhf "$UHF" \
  --T "$OMP_NUM_THREADS" 2>&1 | tee rerank.log
rc=${PIPESTATUS[0]}
set +x
echo "[$(date)] Exit code = $rc"
echo "Return to the repo: crest_ensemble.xyz  (GFN2-reranked ensemble)"
exit $rc
