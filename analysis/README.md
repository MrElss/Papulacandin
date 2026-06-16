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
- `outputs/phase3a_*` — external FKS model, CV metrics, papulacandin scores, joblib
- `outputs/phase3b_*` — within-series serum-tolerance model, LOO predictions, coeffs

## Phase 4 outputs
- `outputs/table_phase4_existing_leads.csv` — best observed serum-active leads
- `outputs/table_phase4_candidate_ranking.csv` — heuristic design-score ranking
- `outputs/design_rules.md` — distilled, source-attributed design rules
- `outputs/fig_phase4_*.png` — presentation figures

## Phase 5 — generative structure optimization (Tier 1 backbone)
`phase5_generate.py` — scaffold-constrained generation of novel Papulacandin-class
analogs by cleaving and re-esterifying the validated aromatic C-6' ester (chemistry
verified by InChIKey round-trip), scored with a transparent multi-objective function
(QED, SA score, Phase-1 design rules, novelty). No learned activity predictor is used.
Emits CREST-ready 3D structures for the Tier-4 QM funnel.

Requires SA_Score from RDKit Contrib (bundled with rdkit).

### Phase 5 outputs
- `outputs/phase5_virtual_library.csv` — generated analogs + multi-objective scores
- `outputs/phase5_top_candidates.sdf` — 3D structures of top novel candidates (CREST-ready)
- `outputs/phase5_crest_commands.sh` — per-candidate CREST/xTB command template (run on QM infra)
- `outputs/phase5_retrospective.png` — design-rule score vs observed serum MIC
- `outputs/phase5_score_distributions.png` — library drug-likeness / SA / novelty
