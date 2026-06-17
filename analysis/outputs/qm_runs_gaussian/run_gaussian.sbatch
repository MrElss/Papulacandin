#!/bin/bash
#SBATCH -p cpu-256G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH -t 48:00:00
#SBATCH -J g16
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err
#
# run_gaussian.sbatch
# ===================
# Final funnel step: B3LYP/6-31G(d) / PCM(water) / Pop=MK single points on the
# Boltzmann-populated conformers of the funnel finalists, to get DFT dipoles and
# ESP (Merz-Kollman) charges -> the electrostatic complement to the SASA shape
# descriptors. The .gjf inputs were written by phase6_qm_layer.py
# (write_gaussian_inputs); each requests %nprocshared=8 / %mem=16GB, matched here.
#
# Usage — submit every populated conformer of one finalist:
#   for f in cand01_quinolinecarbonyl/*.gjf; do
#     sbatch -J $(basename "$f" .gjf) run_gaussian.sbatch "$f"
#   done
# or a single input:
#   sbatch -J cand03_conf01 run_gaussian.sbatch cand03_pyridylphenyl/cand03_pyridylphenyl_conf01_pop0.85.gjf
#
# Finalists (top-3 novel by Phase-5 combined score; all extended rigid aromatics,
# the direction Phase-7 found tracks serum tolerance):
#   cand01_quinolinecarbonyl, cand02_naphthoyl_6OH, cand03_pyridylphenyl
#
# NB: if you first GFN2-re-rank the ensembles (see ../qm_runs/promote_finalists.sh),
# regenerate the .gjf inputs from the GFN2 populations before submitting here, so
# the DFT set matches the refined Boltzmann weights.

set -euo pipefail
GJF=${1:?usage: run_gaussian.sbatch <input.gjf>}
[[ -f "$GJF" ]] || { echo "[ERR] no such input: $GJF"; exit 2; }

# ---- environment (edit to your site's Gaussian module) ----
# module load gaussian/g16
export GAUSS_SCRDIR=${GAUSS_SCRDIR:-/tmp/g16scr.$SLURM_JOB_ID}
mkdir -p "$GAUSS_SCRDIR"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

base="${GJF%.gjf}"
echo "[$(date)] g16 < $GJF  (nproc=$OMP_NUM_THREADS)"
which g16 || { echo "[ERR] g16 not on PATH — load the Gaussian module above"; exit 3; }

g16 < "$GJF" > "${base}.log" 2>&1
rc=$?
echo "[$(date)] exit=$rc  ->  ${base}.log"

# quick parse-back with the phase6 reader (Boltzmann-averages dipole/ESP up to
# the compound level using the CREST-derived weights):
#   python3 -c "import phase6_qm_layer as p6; print(p6.parse_gaussian_log('${base}.log'))"
rm -rf "$GAUSS_SCRDIR"
exit $rc
