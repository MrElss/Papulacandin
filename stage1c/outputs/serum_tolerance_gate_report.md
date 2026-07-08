# Stage 1c — serum-tolerance oracle & generation gate

**Endpoint:** retained vs lost (Stage-0 categorical serum-tolerance call). **Model:** L2 logistic regression, class-balanced. The point of this stage is the GATE decision, validated adversarially against a potency baseline.

## The label set is one congeneric series

- 23 compounds (13 retained, 10 lost); **22 are fusacandin C-6' ester analogs** from one study. Only 1 sits outside that series.
- The envisioned echinocandin->papulacandin transfer is unsupported: no echinocandin serum-shift labels exist (Stage 0). So this is a single series, not a transfer problem.

## Validation (adversarial)

| model | LOO ROC-AUC | balanced acc |
|---|---|---|
| majority baseline | – | 0.565 |
| **potency only** | 0.838 | 0.735 |
| structure only | 0.731 | – |
| potency + structure | 0.885 | 0.823 |

- Structure adds **+0.046** AUC over potency alone (guardrail #1 — a serum model must beat the potency confound).
- Scaffold-grouped AUC (the gate's hard test): **0.885** — but across only 12 generic scaffolds that all share the fusacandin core, so this measures across-ester interpolation, not new-chemotype generalization.

## Within-series SAR hypothesis (constructive)

- Retained analogs carry larger extended-aromatic C-6' esters (median aromatic rings 3.0 vs 1.5 for lost; fsp3 0.49 vs 0.56) and are more potent (pMIC 6.13 vs 5.37).
- **Hypothesis (testable):** biphenyl / terphenyl / naphthoyl C-6' esters favour serum retention in fusacandins. **Caveat:** confounded with potency and unproven beyond this one series — a lead to TEST, not an objective to optimize.

## GATE DECISION

- beats_majority_baseline: **True**
- structure_beats_potency_confound: **False**
- chemotype_diverse_labels: **False**
- scaffold_split_validated: **True**

### Gate: **NOT PASSED** — do NOT proceed to unconstrained generation on the basis of this oracle.

The oracle is, at best, a within-fusacandin-series interpolator whose signal is entangled with potency; with essentially one scaffold it cannot be validated for new chemistry (guardrail #8: the binding constraint is DATA). **Recommendation:** do not gate generation on a 'validated' serum oracle. Either (a) run the Stage-4 design-make-test loop first to generate a scaffold-diverse serum dataset, or (b) if generating now, keep it scaffold-constrained to the fusacandin series, treat the aromatic-ester hypothesis as an experiment, and confirm serum tolerance in the wet lab rather than trusting the score.
