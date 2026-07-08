# Stage 1a — intrinsic-potency oracle (serum-free MIC)

Potency is the dominant, **learnable** axis and the confound every serum-tolerance
analysis must beat (guardrail #1). In the design–make–test loop
(`planning/round1_synthesis_panel_plan.md`) this oracle is the **filter**: only
synthesize/test compounds predicted active enough serum-free, because a compound
dead in broth gives an uninformative serum-tolerance label.

## Run

```bash
python stage1a/build_potency_oracle.py
python -m pytest stage1a/test_stage1a.py -q      # 5 smoke tests
```

## Two models, honestly validated (scaffold split, guardrail #2)

| model | endpoint | scaffold-split result | verdict |
|---|---|---|---|
| **classifier** | active vs inactive serum-free (MIC<100 µg/mL) | **ROC-AUC 0.97**, bal-acc 0.91 (majority 0.71) | ✅ **use it** — the loop filter |
| regressor | pMIC among active compounds | **R² −0.30** (worse than mean), MAE 0.60 | ❌ does not generalize across scaffolds |

**Why the split.** Whether a molecule is active serum-free is driven by having the
right pharmacophore/size — very learnable, and it transfers across scaffolds.
Fine potency *ranking* among already-active analogs is subtle and scaffold-specific;
with the data dominated by one core and a narrow ~3-log potency window, it doesn't
transfer. This is exactly what honest scaffold validation is for — a random split
would have flattered the regressor.

## Endpoint & data

- 132 compounds with serum-free MIC (µg/mL → pMIC via `mw_exact`); 93 active, 39
  inactive/censored-only. The Stage-0 ceiling rule (MIC > 100 µg/mL = inactive) is
  applied; censored values are excluded from the regression, not imputed.
- Features: the 14 shared physicochemical descriptors (`stage1b/featurize.py`).
- Uncertainty (tree spread) and the Stage-1b applicability-domain flag ship with
  every prediction; because this oracle is trained on the target chemotype, its
  in-domain coverage is far better than the drug-trained Stage-1b PPB oracle.

## Outputs (`stage1a/outputs/`)

| file | what it is |
|------|------------|
| `potency_predictions.csv` | per-compound `pred_active_prob`, `pred_pmic` (+SD), `in_domain` for the whole library |
| `potency_oracle_report.md` | validation + per-model verdicts |

## Use in the loop

Panel selection keeps candidates with high `pred_active_prob`, then chooses among
them for diversity and serum-tolerance informativeness. Drop anything predicted
inactive serum-free — it cannot teach us about serum tolerance. Do **not** rank
across novel scaffolds by `pred_pmic` (it failed validation); use it only
in-domain and low-confidence.

## Note on the structure-based half of Stage 1a

START_HERE also envisioned structure-based scoring against the FKS1 cryo-EM
structure. That needs a docking engine and a verified target structure not
available in this environment; it is deferred. The QSAR filter above is the piece
the loop actually needs now, and it is validated.
