# Stage 1c — serum-tolerance oracle & the generation gate

START_HERE (§3) puts a **gate** here: *"do not proceed to generation until the
combined serum-tolerance oracle is validated on a scaffold-split hold-out."*
This stage builds the oracle and renders that gate decision honestly.

## Run

```bash
python stage1c/build_serum_tolerance_oracle.py
python -m pytest stage1c/test_stage1c.py -q      # 5 smoke tests
```

## What the data forced us to reframe

- **No transfer to pretrain.** Stage 0 showed there are no echinocandin
  serum-shift pairs, so the envisioned echinocandin→papulacandin transfer isn't
  supported. This is a single-dataset problem, not a transfer problem.
- **The labels are one congeneric series.** Of 23 modellable compounds, **22 are
  fusacandin C-6′ ester analogs from one study**; only papulacandin B sits
  outside. Out-of-chemotype generalization is essentially untestable (n=1).
- **The endpoint is categorical** (retained/lost) after the Stage-0 ceiling rule.

## Adversarial validation

Logistic regression, class-balanced, leave-one-out CV:

| model | LOO ROC-AUC | balanced acc |
|---|---|---|
| majority baseline | – | 0.57 |
| **potency only** | **0.84** | 0.74 |
| structure only | 0.73 | – |
| potency + structure | 0.89 | 0.82 |

**Potency alone already predicts serum tolerance at AUC 0.84.** Adding all
structural features lifts AUC by only **+0.05** — the "serum oracle" is ~90% a
potency model (guardrail #1: control for potency; a serum signal must beat it).
The scaffold-grouped AUC (0.89) looks strong but spans only generic scaffolds
that all share the fusacandin core, so it measures across-ester interpolation,
not new-chemotype generalization.

## Constructive SAR hypothesis (a lead, not an objective)

Retained analogs carry **larger extended-aromatic C-6′ esters** (median 3 aromatic
rings vs 1.5 for lost; lower fsp3) and are more potent (pMIC 6.1 vs 5.4).
**Hypothesis:** biphenyl / terphenyl / naphthoyl C-6′ esters favour serum
retention in fusacandins — testable, but confounded with potency and unproven
beyond this one series. Treat it as an experiment to run, not a reward to optimize.

## GATE DECISION — **NOT PASSED**

| criterion | result |
|---|---|
| beats majority baseline | ✅ |
| structure beats the potency confound | ❌ (+0.05 AUC) |
| chemotype-diverse labels | ❌ (22/23 fusacandin) |
| scaffold-split validated for new chemistry | ❌ (single core) |

**Do NOT gate generation on this oracle.** It is at best a within-series
interpolator entangled with potency; with one scaffold it cannot be validated for
new chemistry (guardrail #8: the binding constraint is DATA, not methods).

**Recommendation.** Either (a) run the Stage-4 design–make–test loop first to
generate a **scaffold-diverse serum dataset**, or (b) if generating now, keep it
**scaffold-constrained to the fusacandin series**, treat the aromatic-ester
hypothesis as an experiment, and **confirm serum tolerance in the wet lab** rather
than trusting the score. This is the gate doing its job: it stops the prior
project's mistake of generating first and validating the scoring later.

## Outputs (`stage1c/outputs/`)

| file | what it is |
|------|------------|
| `serum_tolerance_gate_report.md` | full validation + gate decision |
| `serum_tolerance_predictions.csv` | per-compound LOO retained-probability + features |

## Where this leaves Stage 1

- **1a** (potency oracle) is still worth building — potency is the dominant,
  learnable axis and the confound everything else must beat.
- **1b** (free fraction) is built, and its AD flag says fu must be measured for
  full-size molecules.
- **1c** (this stage) says the serum-tolerance endpoint is **not yet a validated,
  generalizable objective** — the honest gate outcome is to generate data before
  trusting a serum oracle.
