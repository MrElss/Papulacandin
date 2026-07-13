#!/usr/bin/env bash
# Isolated Mode A only: minimize and CNN-rescore cofolded Boltz/AF3 ligand poses.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
MODE_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
PIPELINE_ROOT=$(cd -- "$MODE_ROOT/.." && pwd)

CONDA_SH=${CONDA_SH:-"$HOME/miniforge3/etc/profile.d/conda.sh"}
if [[ -f "$CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
else
  eval "$(conda shell.bash hook)"
fi
conda activate dock

GNINA=${GNINA:-"$PIPELINE_ROOT/bin/gnina"}
[[ -x "$GNINA" ]] || { echo "[ERR] GNINA is not executable: $GNINA" >&2; exit 1; }
mkdir -p "$MODE_ROOT/out_gnina" "$MODE_ROOT/logs"

found=0
for ligand in "$MODE_ROOT"/cofold_poses/*.sdf; do
  [[ -e "$ligand" ]] || continue
  found=1
  tag=$(basename "$ligand" .sdf)
  template=${tag%%__*}
  receptor="$PIPELINE_ROOT/receptors/${template}_H.pdb"
  output="$MODE_ROOT/out_gnina/${tag}_min.sdf"
  log="$MODE_ROOT/out_gnina/${tag}_min.log"

  [[ -f "$receptor" ]] || { echo "[ERR] missing receptor: $receptor" >&2; exit 1; }
  if [[ -s "$output" && -s "$log" ]]; then
    echo "[skip] existing result: $output"
    continue
  fi

  echo "### Mode A minimize+rescore: $tag"
  "$GNINA" -r "$receptor" -l "$ligand" \
    --minimize --cnn_scoring rescore -o "$output" 2>&1 | tee "$log"
done

[[ "$found" -eq 1 ]] || { echo "[ERR] no SDF files in $MODE_ROOT/cofold_poses" >&2; exit 1; }
echo "Mode A GNINA complete: $MODE_ROOT/out_gnina"
