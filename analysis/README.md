# Analysis pipeline

Reproducible analysis for the Papulacandin / FKS1 serum-gap project.
Each script is self-contained, documented, and writes to `analysis/outputs/`.

## Environment
```
python3 -m pip install pandas numpy rdkit scikit-learn
```

## Scripts (run in order)
1. `serum_gap_analysis.py` — builds the matched serum-free vs serum MIC table
   for Candida albicans and computes the per-compound "serum shift". This is
   the project's dependent variable. Handles censored values (`>100`, `<0.03`)
   explicitly rather than discarding them.

## Outputs
- `outputs/serum_gap_pairs.csv` — matched pairs + serum shift per compound
- `outputs/serum_gap_summary.txt` — human-readable summary
- `outputs/phase1_descriptor_stats.csv` + `phase1_*.png` — serum-shift SAR
- `outputs/phase2_activity_ladder.csv` + `phase2_*.png` — activity ladder / attrition
