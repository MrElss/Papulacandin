#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=52
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J jobname_gfn2rerank
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#
# GFN2 RE-RANK (finalist tier)
# ============================
# Re-evaluate and re-rank the GFN-FF screening-tier conformer POOL at GFN2-xTB,
# WITHOUT a new conformer search ("--screen" / multilevel refinement). GFN2 is a
# semi-empirical QM method and ranks these polyol glycolipids' many intramolecular
# H-bonds far better than the GFN-FF force field used for the search, so the
# resulting Boltzmann populations -- which drive the SASA descriptors and decide
# which conformers go to DFT -- are more trustworthy. Cheap relative to the
# original search, and cheap insurance before the expensive Gaussian step.
#
# Input : crest_conformers.xyz  (the GFN-FF pool, already in this folder)
# Output: crest_ensemble.xyz / crest.out  (GFN2-reranked), tee'd to rerank.log
#
# Submit from inside this folder:
#   sbatch -J cand0X_gfn2 run_gfn2_rerank.sbatch

set -euo pipefail
ulimit -s unlimited

# ========== 1) environment (mirrors run_crest.sbatch) ==========
APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "[$(date)] node=${SLURM_NODELIST:-NA} T=$OMP_NUM_THREADS jobid=${SLURM_JOB_ID:-NA}"
which crest; which xtb

# ========== 2) input ==========
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"

ENS=${ENS:-crest_conformers.xyz}
[[ -f "$ENS" ]] || { echo "[ERR] no conformer pool '$ENS' in $WORKDIR"; exit 2; }
# work on a copy under a distinct name so CREST never clobbers the input pool
cp -f "$ENS" pool_gfnff.xyz

# ========== 3) parameters ==========
# GFN2 supports ALPB solvation (same water model as the GFN-FF search). The
# wider -ewin 10 keeps everything within 10 kcal/mol of the NEW (GFN2) minimum,
# so conformers that GFN-FF mis-ranked are not lost before re-evaluation.
METHOD_FLAG=${METHOD_FLAG:---gfn2}
SOLV_FLAG=${SOLV_FLAG:---alpb water}
EWIN=${EWIN:--ewin 10}
OUT=${OUT:-rerank.log}

opts=(--screen pool_gfnff.xyz)
[[ -n "$METHOD_FLAG" ]] && opts+=($METHOD_FLAG)
[[ -n "$SOLV_FLAG"   ]] && opts+=($SOLV_FLAG)
[[ -n "$EWIN"        ]] && opts+=($EWIN)

# ========== 4) run ==========
echo "[$(date)] crest ${opts[*]} --T $OMP_NUM_THREADS"
set -x
crest "${opts[@]}" --T "$OMP_NUM_THREADS" 2>&1 | tee "$OUT"
rc=${PIPESTATUS[0]}
set +x

# ========== 5) finish ==========
echo "[$(date)] Exit code = $rc"
echo "GFN2-reranked ensemble -> crest_ensemble.xyz (and/or crest_conformers.xyz)"
exit $rc
