#!/usr/bin/env bash
# Mode A only: minimize and CNN-rescore 15 aligned Boltz poses with GNINA v1.1.
#SBATCH --job-name=modeA_gnina11
#SBATCH --partition=gpu
#SBATCH --gres=gpu:A100:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713/logs/modeA_gnina11_%j.out

set -euo pipefail
GNINA=/home/share/zhangz/Apps/gnina-1.1/gnina
JOB_ROOT=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713
export LD_LIBRARY_PATH=/home/share/zhangz/Apps/envs/diffdock/lib:/home/share/zhangz/Apps/cuda-compat-11-8/usr/local/cuda-11.8/compat:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export TMPDIR="$JOB_ROOT/tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$JOB_ROOT/out_gnina"

echo "START=$(date --iso-8601=seconds)"
echo "HOST=$(hostname)"
echo "CPUS=$SLURM_CPUS_PER_TASK"
if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "[warn] nvidia-smi is not in PATH; GNINA --device 0 remains the functional GPU check"
fi
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
  rm -f "$output" "$log"

  echo "### $count/15 Mode A minimize+rescore: $tag"
  "$GNINA" -r "$receptor" -l "$ligand" --device 0 \
    --minimize --cnn_scoring rescore -o "$output" 2>&1 | tee "$log"
  [[ -s "$output" ]] || { echo "Missing output for $tag" >&2; exit 1; }
done

[[ "$count" -eq 15 ]] || { echo "Expected 15 poses, found $count" >&2; exit 1; }
echo "OUTPUT_SDF_COUNT=$(find "$JOB_ROOT/out_gnina" -maxdepth 1 -name '*_min.sdf' -size +0c | wc -l)"
echo "END=$(date --iso-8601=seconds)"
echo "MODE_A_GNINA_V11_OK"
