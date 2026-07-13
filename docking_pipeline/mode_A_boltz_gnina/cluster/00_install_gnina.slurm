#!/usr/bin/env bash
# Install the official GNINA compatibility binary and cuDNN 9 on a compute node.
#SBATCH --job-name=install_gnina
#SBATCH --partition=cpu-256G
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713/logs/install_gnina_%j.out

set -euo pipefail
APP_ROOT=/home/share/zhangz/Apps/gnina-1.3.2
GNINA="$APP_ROOT/gnina"
DEPS="$APP_ROOT/deps"
URL=https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2
EXPECTED_SIZE=1426790536

mkdir -p "$APP_ROOT" "$DEPS" /home/zhangzhuan/tmp
export TMPDIR=/home/zhangzhuan/tmp

current_size=0
[[ -f "$GNINA" ]] && current_size=$(stat -c %s "$GNINA")
if [[ "$current_size" -ne "$EXPECTED_SIZE" ]]; then
  echo "Downloading GNINA v1.3.2 compatibility binary"
  rm -f "$GNINA"
  tmp="$GNINA.download"
  rm -f "$tmp"
  curl --fail --location --retry 5 \
    --connect-timeout 60 --output "$tmp" "$URL"
  mv "$tmp" "$GNINA"
fi
[[ $(stat -c %s "$GNINA") -eq "$EXPECTED_SIZE" ]] || { echo "GNINA size check failed" >&2; exit 1; }
chmod 755 "$GNINA"

echo "--- GNINA linkage before cuDNN setup ---"
ldd "$GNINA" || true

if ldd "$GNINA" | grep -q 'libcudart.so.11'; then
  CUDNN_PACKAGE=nvidia-cudnn-cu11
elif ldd "$GNINA" | grep -q 'libcudart.so.12'; then
  CUDNN_PACKAGE=nvidia-cudnn-cu12
else
  echo "Cannot determine GNINA CUDA major version from ldd" >&2
  exit 1
fi

if ! find "$DEPS" -name 'libcudnn.so.9' -print -quit | grep -q .; then
  echo "Installing $CUDNN_PACKAGE under $DEPS"
  /opt/app/anaconda3/bin/python -m pip install --upgrade --target "$DEPS" "$CUDNN_PACKAGE"
fi

CUDNN_LIB=$(find "$DEPS" -type d -path '*/nvidia/cudnn/lib' -print -quit)
CUBLAS_LIB=$(find "$DEPS" -type d -path '*/nvidia/cublas/lib' -print -quit)
[[ -n "$CUDNN_LIB" ]] || { echo "cuDNN library directory not found" >&2; exit 1; }
export LD_LIBRARY_PATH="$CUDNN_LIB${CUBLAS_LIB:+:$CUBLAS_LIB}:${LD_LIBRARY_PATH:-}"

echo "--- GNINA linkage after cuDNN setup ---"
ldd "$GNINA"
if ldd "$GNINA" | grep -q 'not found'; then
  echo "Unresolved GNINA shared libraries remain" >&2
  exit 1
fi

"$GNINA" --version
echo "GNINA_INSTALL_OK"
