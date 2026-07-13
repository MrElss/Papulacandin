#!/usr/bin/env bash
#SBATCH --job-name=cache_esm2_650M
#SBATCH --partition=cpu-256G
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/home/share/zhangz/Jobs/cache_esm2_650M_%j.out

set -euo pipefail

APP_ROOT=/home/share/zhangz/Apps
ENV_ROOT="$APP_ROOT/envs/diffdock"
export TORCH_HOME="$APP_ROOT/torch_cache"
CACHE_DIR="$TORCH_HOME/hub/checkpoints"
MODEL=esm2_t33_650M_UR50D

mkdir -p "$CACHE_DIR"
wget -c --tries=20 --timeout=60 \
  -O "$CACHE_DIR/$MODEL.pt" \
  "https://dl.fbaipublicfiles.com/fair-esm/models/$MODEL.pt"
wget -c --tries=20 --timeout=60 \
  -O "$CACHE_DIR/$MODEL-contact-regression.pt" \
  "https://dl.fbaipublicfiles.com/fair-esm/regression/$MODEL-contact-regression.pt"

ls -lh "$CACHE_DIR/$MODEL.pt" "$CACHE_DIR/$MODEL-contact-regression.pt"
"$ENV_ROOT/bin/python" -c \
  'from esm import pretrained; pretrained.load_model_and_alphabet("esm2_t33_650M_UR50D"); print("ESM-2 cache validation passed")'
