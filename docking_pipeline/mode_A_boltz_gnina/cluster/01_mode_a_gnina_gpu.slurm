#!/usr/bin/env bash
# Mode A only: minimize and CNN-rescore all aligned Boltz poses on one A100.
#SBATCH --job-name=modeA_gnina
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713/logs/modeA_gnina_%j.out

set -euo pipefail
APP_ROOT=/home/share/zhangz/Apps/gnina-1.3.2
GNINA="$APP_ROOT/gnina"
DEPS="$APP_ROOT/deps"
JOB_ROOT=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713

CUDNN_LIB=$(find "$DEPS" -type d -path '*/nvidia/cudnn/lib' -print -quit)
CUBLAS_LIB=$(find "$DEPS" -type d -path '*/nvidia/cublas/lib' -print -quit)
[[ -n "$CUDNN_LIB" ]] || { echo "cuDNN library directory not found" >&2; exit 1; }
export LD_LIBRARY_PATH="$CUDNN_LIB${CUBLAS_LIB:+:$CUBLAS_LIB}:/home/share/zhangz/Apps/cuda-compat-12-4/usr/local/cuda-12.4/compat:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export TMPDIR="$JOB_ROOT/tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$JOB_ROOT/out_gnina"

echo "START=$(date --iso-8601=seconds)"
echo "HOST=$(hostname)"
echo "CPUS=$SLURM_CPUS_PER_TASK"
nvidia-smi
"$GNINA" --version

count=0
for ligand in "$JOB_ROOT"/cofold_poses/*.sdf; do
  [[ -e "$ligand" ]] || continue
  count=$((count + 1))
  tag=$(basename "$ligand" .sdf)
  template=${tag%%__*}
  receptor="$JOB_ROOT/receptors/${template}_H.pdb"
  output="$JOB_ROOT/out_gnina/${tag}_min.sdf"
  log="$JOB_ROOT/out_gnina/${tag}_min.log"
  [[ -f "$receptor" ]] || { echo "Missing receptor: $receptor" >&2; exit 1; }

  if [[ -s "$output" && -s "$log" ]]; then
    echo "[skip] $tag"
    continue
  fi

  echo "### $count/15 Mode A minimize+rescore: $tag"
  "$GNINA" -r "$receptor" -l "$ligand" \
    --minimize --cnn_scoring rescore -o "$output" 2>&1 | tee "$log"
  [[ -s "$output" ]] || { echo "Missing output for $tag" >&2; exit 1; }
done

[[ "$count" -eq 15 ]] || { echo "Expected 15 poses, found $count" >&2; exit 1; }
echo "OUTPUT_SDF_COUNT=$(find "$JOB_ROOT/out_gnina" -maxdepth 1 -name '*_min.sdf' -size +0c | wc -l)"
echo "END=$(date --iso-8601=seconds)"
echo "MODE_A_GNINA_OK"
