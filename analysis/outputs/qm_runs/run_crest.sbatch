#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=52
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -J jobname_crest
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -euo pipefail
ulimit -s unlimited

# ========== 1) 环境 ==========
APPS=/home/duxiaonan/share/duxiaonan/apps
source "$APPS/bin/use-crest"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS
export MKL_NUM_THREADS=$OMP_NUM_THREADS

echo "[$(date)] node=${SLURM_NODELIST:-NA} T=$OMP_NUM_THREADS jobid=${SLURM_JOB_ID:-NA}"
which crest; which xtb

# ========== 2) 目录与输入 ==========
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$WORKDIR"

# 显式指定或自动寻找输入（优先 .pdb）
INP=${INP:-}
if [[ -z "$INP" ]]; then
  INP=$(ls -1 *.pdb 2>/dev/null | head -n1 || true)
  [[ -n "$INP" ]] || INP=$(ls -1 *.xyz 2>/dev/null | head -n1 || true)
fi
[[ -n "$INP" ]] || { echo "[ERR] No *.pdb or *.xyz in $WORKDIR"; exit 2; }

OUT=${OUT:-run.log}

# ========== 3) 计算参数 ==========
# FAST SCREENING TIER (12 funnel candidates): GFN-FF search, GFN2 refine on
# survivors only ("gfn2//gfnff" composite) -- the prior pure-GFN2/QUICK=no/
# EWIN=10 setup was projected at ~10 days/compound for these ~140-atom,
# 30+ rotatable-bond glycolipids, i.e. ~120 days serial for the funnel.
# Reserve pure --gfn2 (QUICK=no, EWIN=10) for the 1-3 finalists that survive
# this screening tier, by overriding METHOD_FLAG/EXTRA_OPTS at submission time:
#   METHOD_FLAG=--gfn2 EXTRA_OPTS="-ewin 10" sbatch -J <name> run_crest.sbatch
METHOD_FLAG=${METHOD_FLAG:---gfn2//gfnff}      # screening default; finalists: --gfn2
SOLV_FLAG=${SOLV_FLAG:---gbsa water}           # 不用溶剂：SOLV_FLAG=""
CHARGE=${CHARGE:-0}
UHF=${UHF:-0}
EXTRA_OPTS=${EXTRA_OPTS:---quick -ewin 6}      # screening default; finalists: "-ewin 10"

# 把可选参数安全地累加到一个数组（即使为空也不会报错）
opts=()
[[ -n "$METHOD_FLAG" ]] && opts+=($METHOD_FLAG)
[[ -n "$SOLV_FLAG"   ]] && opts+=($SOLV_FLAG)
[[ -n "$EXTRA_OPTS"  ]] && opts+=($EXTRA_OPTS)

# ========== 4) 运行 ==========
echo "[$(date)] Run CREST: $INP  T=$OMP_NUM_THREADS  ${METHOD_FLAG}  solvent='${SOLV_FLAG:-none}'  q=$CHARGE uhf=$UHF  extra='${EXTRA_OPTS:-none}'"
set -x
crest "$INP" "${opts[@]}" --T "$OMP_NUM_THREADS" --chrg "$CHARGE" --uhf "$UHF" \
  2>&1 | tee "$OUT"
rc=${PIPESTATUS[0]}
set +x

# ========== 5) 收尾 ==========
echo "[$(date)] Exit code = $rc"
echo "Results: crest.out, crest_best.xyz, crest_ensemble.xyz, crest_conformers/"
exit $rc
