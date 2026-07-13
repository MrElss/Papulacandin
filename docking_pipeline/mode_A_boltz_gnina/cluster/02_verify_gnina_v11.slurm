#!/usr/bin/env bash
# Verify GNINA v1.1 against the cluster's CentOS 7-era runtime on a compute node.
#SBATCH --job-name=verify_gnina11
#SBATCH --partition=cpu-256G
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713/logs/verify_gnina11_%j.out

set -euo pipefail
APP_ROOT=/home/share/zhangz/Apps/gnina-1.1
GNINA="$APP_ROOT/gnina"
EXPECTED_SIZE=306470832

[[ -f "$GNINA" ]] || { echo "GNINA v1.1 binary is missing" >&2; exit 1; }
[[ $(stat -c %s "$GNINA") -eq "$EXPECTED_SIZE" ]] || { echo "GNINA v1.1 size check failed" >&2; exit 1; }
chmod 755 "$GNINA"

# Reuse the newer libstdc++ already present in the shared DiffDock environment.
export LD_LIBRARY_PATH=/home/share/zhangz/Apps/envs/diffdock/lib:${LD_LIBRARY_PATH:-}
echo "--- linkage ---"
ldd "$GNINA"
if ldd "$GNINA" | grep -q 'not found'; then
  echo "Unresolved GNINA shared libraries remain" >&2
  exit 1
fi
"$GNINA" --version
"$GNINA" --help 2>&1 | grep -E -- '--device|--no_gpu'
echo "GNINA_V11_VERIFY_OK"
