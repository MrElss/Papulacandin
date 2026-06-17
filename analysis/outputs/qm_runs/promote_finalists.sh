#!/bin/bash
# promote_finalists.sh
# ====================
# Promote the 1-3 Phase-7-validated funnel finalists from the GFN-FF screening
# tier to a GFN2 refinement, WITHOUT redoing the expensive conformer search.
#
# Two promotion modes (pick per finalist; --screen is the cheap default):
#
#   (A) GFN2 RE-RANK of the existing GFN-FF ensemble  [cheap, recommended first]
#       Reranks crest_conformers.xyz at GFN2//ALPB(water) and re-prunes within a
#       wider 10 kcal/mol window. No new sampling; minutes-to-hours, not days.
#         crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 10 --T 52
#
#   (B) Full GFN2 conformer SEARCH from the GFN-FF best geometry  [expensive]
#       Only if (A) reshuffles the populated set enough to matter. Submits the
#       finalist's crest_best.xyz back through run_crest.sbatch with GFN2 flags.
#
# Usage:
#   ./promote_finalists.sh cand07_biphenyl_4CONH2 cand03_pyridylphenyl ...
# Edit FINALISTS below, or pass candidate folder names as arguments.

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

FINALISTS=("$@")
if [[ ${#FINALISTS[@]} -eq 0 ]]; then
  echo "No finalists given. Pass candidate folder names, e.g.:"
  echo "  ./promote_finalists.sh cand07_biphenyl_4CONH2 cand03_pyridylphenyl"
  exit 1
fi

for cand in "${FINALISTS[@]}"; do
  src="$HERE/$cand"
  [[ -d "$src" ]] || { echo "[skip] no folder $src"; continue; }
  ens="$src/crest_conformers.xyz"
  [[ -f "$ens" ]] || { echo "[skip] no crest_conformers.xyz in $src"; continue; }

  fin="$HERE/${cand}_gfn2"
  mkdir -p "$fin"
  cp "$ens" "$fin/crest_conformers.xyz"
  # carry the best geometry too, for optional mode (B)
  [[ -f "$src/crest_best.xyz" ]] && cp "$src/crest_best.xyz" "$fin/"

  echo "[ready] $cand -> $fin"
  echo "  (A) cheap GFN2 re-rank (recommended):"
  echo "      cd '$fin' && crest --screen crest_conformers.xyz --gfn2 --alpb water -ewin 10 --T 52 2>&1 | tee rerank.log"
  echo "  (B) full GFN2 search from best geometry (expensive):"
  echo "      cd '$fin' && INP=crest_best.xyz METHOD_FLAG=--gfn2 SOLV_FLAG='--alpb water' EXTRA_OPTS='-ewin 10' sbatch -J ${cand}_final ../run_crest.sbatch"
done

echo
echo "After GFN2 refinement, re-run phase6 on the *_gfn2 ensembles to refresh"
echo "the Boltzmann-weighted descriptors with GFN2-quality populations:"
echo "  python3 -c \"import phase6_qm_layer as p6; p6.run_qm_layer("
echo "    'outputs/phase5_top_candidates.sdf', 'outputs/qm_runs',"
echo "    'outputs/phase6_qm_descriptors_gfn2.csv', 'outputs/qm_runs_gaussian_gfn2', real_run=True)\""
echo "(point crest_root at a tree whose subfolders are the *_gfn2 ensembles, or"
echo " symlink them to the candidate names phase6 expects.)"
