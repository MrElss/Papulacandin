# Stage 1a — intrinsic-potency oracle (serum-free MIC)

Potency is the learnable, dominant axis and the confound serum-tolerance work must beat (guardrail #1). In the design-make-test loop this oracle is the **filter**: only test compounds predicted active enough serum-free to give an informative serum-tolerance result.

## Active/inactive classifier (the loop filter) — USABLE

- 132 compounds (93 active, 39 inactive serum-free by the >100 ug/mL ceiling rule).
- **Scaffold-split ROC-AUC = 0.966**, balanced acc = 0.912 (majority baseline 0.705).
- **Verdict:** strongly beats baseline and generalizes across scaffolds — whether a compound is active serum-free is genuinely learnable. Use this as the loop filter.

## pMIC regressor (graded ranking) — DOES NOT GENERALIZE

- 93 active compounds, pMIC 3.79–7.0 (a narrow ~3-log window).
- **Scaffold-split R² = -0.297** (WORSE than the mean-predictor baseline R² = 0.0), MAE = 0.597 log-units.
- **Verdict:** fine potency ranking among active compounds does NOT transfer across scaffolds — the data is dominated by one core, and within-window potency differences are scaffold-specific. Do **not** rely on `pred_pmic` to rank across novel scaffolds; treat it as in-domain, low-confidence only. This is guardrail #2 in action: honest scaffold validation catches a model that a random split would have flattered.

## Applicability domain

- 129/132 library compounds in domain (the potency oracle is trained on this chemotype, so coverage is far better than the drug-trained Stage-1b PPB oracle — but still check the flag per design).

## Use in the loop

`potency_predictions.csv` gives `pred_active_prob` and `pred_pmic` (+SD, +AD flag) for every library compound. Panel selection keeps only candidates with high `pred_active_prob`, then chooses among them for diversity and serum-tolerance informativeness. A candidate predicted inactive serum-free is dropped — it cannot teach us about serum tolerance.
