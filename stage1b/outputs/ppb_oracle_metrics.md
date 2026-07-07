# Stage 1b — free-fraction / PPB oracle

Endpoint: **log10(% unbound)** (fu = 10^pred / 100). Trained on the public Fang et al. 2023 Computational-ADME human-PPB set.

## Training data

- 194 human-PPB labelled drugs (MW 151-666).

## Validation (guardrail #2: validate the hard way)

- **Scaffold-split hold-out (honest):** R² = 0.388, MAE = 0.596 log-units, Spearman = 0.659 (n = 39).
- Random 5-fold (optimistic, for contrast): R² = 0.451, MAE = 0.484. The gap is the usual scaffold penalty — trust the scaffold number.

## Anchor check — anidulafungin

- Predicted 96.02% bound (3.979% unbound); in-repo observed ~87.7% bound. In domain: False. (Literature anidulafungin PPB ≈ 99%; a large echinocandin, so an out-of-domain query — read with the AD flag.)

## Applicability domain — the decisive result

- **Papulacandins in domain: 2/137 (1%).**
- Echinocandins/FKS in domain: 91/321 (28%).

The papulacandin glycolipids (median MW ~930) sit almost entirely outside the training drugs' MW range (≤666), so the oracle must **extrapolate** for them and its predictions there are not trustworthy — exactly guardrail #4 (off-the-shelf QSAR does not transfer to this bRo5 chemotype). This operationalizes guardrails #5 and #8: for the real design targets, free fraction must be **measured** (equilibrium dialysis), not predicted. Use the oracle to *rank in-domain* analogs and to *flag* which designs even fall in a modellable region — not to assign an fu to a papulacandin.

## Files

- `free_fraction_predictions.csv` — per-compound fu + uncertainty + AD flag.
- `applicability_domain_report.csv` — in-domain fractions by chemotype.
