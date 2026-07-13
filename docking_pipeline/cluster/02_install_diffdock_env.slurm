#!/usr/bin/env bash
#SBATCH --job-name=install_diffdock_env
#SBATCH --partition=cpu-256G
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/home/share/zhangz/Jobs/install_diffdock_env_%j.out

set -euo pipefail

APP_ROOT=/home/share/zhangz/Apps
MM_DIR="$APP_ROOT/micromamba"
MAMBA_ROOT_PREFIX="$APP_ROOT/micromamba_root"
ENV_DIR="$APP_ROOT/envs/diffdock"
SOURCE_DIR="$APP_ROOT/DiffDock"
SOURCE_ARCHIVE="$APP_ROOT/DiffDock_source_85c49b6.tar.gz"
TORCH_WHEEL="$APP_ROOT/wheels/torch-1.13.1+cu117-cp39-cp39-linux_x86_64.whl"
PUBLIC_CUDA=/opt/app/cuda/11.8

export MAMBA_ROOT_PREFIX
mkdir -p "$MM_DIR" "$MAMBA_ROOT_PREFIX" "$APP_ROOT/envs" "$SOURCE_DIR"

if [[ ! -x "$MM_DIR/micromamba" ]]; then
  curl -LsS --retry 10 https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$MM_DIR" --strip-components=1 bin/micromamba
fi

test -s "$SOURCE_ARCHIVE"
test -s "$TORCH_WHEEL"
test -x "$PUBLIC_CUDA/bin/nvcc"
test -s "$PUBLIC_CUDA/include/cusparse.h"
tar -xzf "$SOURCE_ARCHIVE" -C "$SOURCE_DIR"

if [[ ! -x "$ENV_DIR/bin/python" ]]; then
  "$MM_DIR/micromamba" create -y -p "$ENV_DIR" -c conda-forge \
    python=3.9.18 pip python-devtools setuptools=69.5.1 prody=2.2.0 scipy=1.12.0
fi

PY="$ENV_DIR/bin/python"
TORCH_LIB="$ENV_DIR/lib/python3.9/site-packages/torch/lib"
export LD_LIBRARY_PATH="$ENV_DIR/targets/x86_64-linux/lib:$ENV_DIR/lib:$TORCH_LIB:${LD_LIBRARY_PATH:-}"
"$PY" -m pip install "$TORCH_WHEEL"
"$PY" -m pip install --timeout 180 --retries 10 \
  --find-links https://data.pyg.org/whl/torch-1.13.1+cu117.html \
  e3nn==0.5.1 'fair-esm[esmfold]==2.0.0' networkx==2.8.4 pandas==1.5.1 \
  pybind11==2.11.1 pytorch-lightning==1.9.5 rdkit==2022.03.3 \
  scikit-learn==1.1.0 torch-cluster==1.6.0+pt113cu117 \
  torch-geometric==2.2.0 torch-scatter==2.1.0+pt113cu117 \
  torch-sparse==0.6.16+pt113cu117 torch-spline-conv==1.2.1+pt113cu117 \
  torchmetrics==0.11.0 'gradio==3.50.*' requests

"$PY" -m pip install "$SOURCE_DIR/external/dllogger"

if [[ ! -f "$ENV_DIR/.openfold_sm80_built" ]]; then
  "$MM_DIR/micromamba" install -y -p "$ENV_DIR" \
    -c conda-forge gxx_linux-64=11
  "$PY" -m pip uninstall -y openfold || true
  rm -rf "$SOURCE_DIR/external/openfold/build" "$SOURCE_DIR/external/openfold/openfold.egg-info"
  export CUDA_HOME="$PUBLIC_CUDA"
  export PATH="$CUDA_HOME/bin:$PATH"
  export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$LD_LIBRARY_PATH"
  export CC="$ENV_DIR/bin/x86_64-conda-linux-gnu-cc"
  export CXX="$ENV_DIR/bin/x86_64-conda-linux-gnu-g++"
  export CPLUS_INCLUDE_PATH="$CUDA_HOME/include:${CPLUS_INCLUDE_PATH:-}"
  export LIBRARY_PATH="$CUDA_HOME/lib64:$ENV_DIR/lib:${LIBRARY_PATH:-}"
  export TORCH_CUDA_ARCH_LIST=8.0
  export MAX_JOBS="${SLURM_CPUS_PER_TASK:-8}"
  "$PY" -m pip install --no-build-isolation "$SOURCE_DIR/external/openfold"
  touch "$ENV_DIR/.openfold_sm80_built"
fi

cd "$SOURCE_DIR"
"$PY" utils/precompute_series.py
"$PY" -c 'import torch, torch_geometric, torch_scatter, torch_sparse, torch_cluster, e3nn, esm, openfold; print(torch.__version__, torch.version.cuda, torch_geometric.__version__)'
