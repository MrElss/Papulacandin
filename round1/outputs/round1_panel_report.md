# Round 1 — informative serum-test panel

Turns the oracles into a concrete first experiment that *generates* the scaffold-diverse serum data the Stage-1c gate lacked. Objective: maximum information about serum tolerance per cost — not predicted potency.

## The panel

- **12 test candidates** + **4 controls** (known serum-retained / serum-lost, to calibrate the assay).
- **9 candidates sit on scaffolds ABSENT from the current serum set** — these directly test whether serum tolerance generalises beyond the fusacandin series (the gate's blind spot).
- Candidates span 6 scaffolds and clogP 1.7–3.1 (the mechanism-key lipophilicity/sequestration axis).
- All are obtainable without new synthesis (fermentation naturals + cheap semisynthetic), so round 1 is just the serum assay on a diverse obtainable set — the cheapest way to break the single-series limit.

## Why this design (and not 'make the predicted-best')

Ranking by predicted potency would re-pick potent compounds that may die in serum — the class's 50-year failure. Instead the potency oracle is a FILTER (all candidates are predicted active serum-free, so a serum result is interpretable), and selection maximises scaffold + clogP diversity so the resulting labels can finally *separate serum tolerance from potency* and *generalise past fusacandins*.

## Assay & loop

Measure serum-relevant endpoints directly: protein-adjusted MIC ± albumin (binary retained/lost + shift) and equilibrium-dialysis fu. Feed results back to Stage 0 (binary labels) and re-validate the Stage-1c gate. Round 2 then uses generative decoration to probe the clogP / aromatic-ester lever with potency-matched pairs.

## Note on the free-fraction column

`fu_pred` is shown for reference but most candidates are OUT of the Stage-1b oracle's domain (large glycolipids) — treat fu as measure-not-predict (guardrails #4/#5). The equilibrium-dialysis fu in round 1 is what calibrates it.
