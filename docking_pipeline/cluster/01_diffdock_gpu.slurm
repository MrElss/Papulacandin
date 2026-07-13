#!/usr/bin/env bash
#SBATCH --job-name=papB_diffdock
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/share/zhangz/Jobs/docking_pipeline_20260713/logs/diffdock_%j.out

set -euo pipefail

JOB_ROOT=/home/share/zhangz/Jobs/docking_pipeline_20260713
APP_ROOT=/home/share/zhangz/Apps
ENV_ROOT="$APP_ROOT/envs/diffdock"
SOURCE_ROOT="$APP_ROOT/DiffDock"
MODEL_ROOT="$APP_ROOT/DiffDock_models/v1.1"
PUBLIC_CUDA=/opt/app/cuda/11.8
CUDA_COMPAT="$APP_ROOT/cuda-compat-12-4/usr/local/cuda-12.4/compat"

mkdir -p "$JOB_ROOT"/{logs,out_diffdock,cache,tmp}
test -x "$ENV_ROOT/bin/python"
test -s "$MODEL_ROOT/score_model/best_ema_inference_epoch_model.pt"
test -s "$MODEL_ROOT/confidence_model/best_model_epoch75.pt"
test -s "$JOB_ROOT/diffdock_input.csv"

export PYTHONNOUSERSITE=1
export PYTHONFAULTHANDLER=1
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export TORCH_HOME="$APP_ROOT/torch_cache"
TORCH_LIB="$ENV_ROOT/lib/python3.9/site-packages/torch/lib"
test -s "$CUDA_COMPAT/libcuda.so.1"
export LD_LIBRARY_PATH="$CUDA_COMPAT:$PUBLIC_CUDA/lib64:$ENV_ROOT/lib:$TORCH_LIB:${LD_LIBRARY_PATH:-}"
"$ENV_ROOT/bin/python" -c 'import torch; assert torch.cuda.is_available(); print("GPU:", torch.cuda.get_device_name(0))'

cd "$SOURCE_ROOT"
"$ENV_ROOT/bin/python" -X faulthandler -m inference \
  --config "$JOB_ROOT/diffdock_inference.yaml"

find "$JOB_ROOT/out_diffdock" -maxdepth 2 -type f -printf '%P\t%s bytes\n' | sort
