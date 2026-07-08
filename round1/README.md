# Round 1 — informative serum-test panel

The Stage-1c gate did not pass: the serum labels are one fusacandin series, so no
serum-tolerance oracle can be validated (guardrail #8 — the binding constraint is
DATA). This selector turns the oracles into the concrete first experiment that
*generates* the missing data, following `planning/round1_synthesis_panel_plan.md`.

## Run

```bash
python round1/build_panel.py
python -m pytest round1/test_round1.py -q       # 5 smoke tests
```

## What it selects (and why)

Objective = **maximum information about serum tolerance per cost**, not predicted
potency (that would re-pick potent compounds that die in serum — the class's
50-year failure).

- **Filter (Stage 1a):** only compounds predicted **active serum-free** — a
  broth-dead compound can't demonstrate serum tolerance, so its serum result is
  uninformative.
- **Novelty:** only compounds with **no serum label yet**.
- **Diversity:** greedy coverage of **scaffold × clogP** cells, with a bonus for
  scaffolds **absent from the current serum set** — the mechanism-key axis is
  lipophilicity / amphiphile sequestration (`serum_mechanism_evidence.md`).
- **Cheap first:** obtainable-by-fermentation naturals + cheap semisynthetic
  analogs — **round 1 needs no new chemistry**, just the serum assay on a diverse
  obtainable set.
- **Controls:** known serum-retained and serum-lost papulacandins to calibrate
  the assay.

## The round-1 panel

- **12 candidates + 4 controls.** 9 candidates sit on **scaffolds absent from the
  fusacandin-dominated serum set** — pestiorosins, pestiocandin, saricandin,
  furanocandin, corynecandin, F-10748, plus 2 cheap papulacandin B semisynthetics.
- Candidates span **6 scaffolds, clogP ~1.7–3.1**.
- Controls: fusacandin 6l/6j (retained) and 6u/6i (lost).

This is deliberately a **scaffold-diversity survey of obtainable papulacandins**:
its job is to test whether serum tolerance generalises beyond the one series we
have data for — the exact blind spot that failed the gate — at minimal cost.

## Assay → loop

Measure directly: protein-adjusted MIC ± albumin (binary retained/lost + shift)
and equilibrium-dialysis fu. Feed results back into the Stage-0 binary labels,
re-validate the Stage-1c gate, then **round 2** uses generative scaffold-decoration
to probe the clogP / aromatic-ester lever with **potency-matched pairs** (to
decouple serum tolerance from potency).

## Outputs (`round1/outputs/`)

| file | what it is |
|------|------------|
| `round1_synthesis_panel.csv` | ranked panel: role, compound, scaffold, clogP, MW, pred_active_prob, fu, reason |
| `round1_panel_report.md` | the design rationale |

## Caveat

`fu_pred` is shown but most candidates are **out of the Stage-1b oracle's domain**
(large glycolipids) — treat fu as measure-not-predict (guardrails #4/#5). The
round-1 equilibrium-dialysis fu is what calibrates it.
