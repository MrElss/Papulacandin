# Stage 0 — unify data & fix the endpoint

First step of the methodology in [`../START_HERE.md`](../START_HERE.md) (§3).
Its only job is to turn the raw activity records into **one potency-independent
endpoint both chemotypes can share** — the *serum shift* — and to carry the
confounders every downstream oracle needs. **No modeling happens here.**

```
serum_shift = MIC(serum-containing) / MIC(serum-free)      # same compound,
                                                           # same organism
```

## Run

```bash
python stage0/build_serum_shift_table.py     # pure stdlib, no deps
python -m pytest stage0/test_stage0.py -q     # 6 smoke tests
```

## Outputs (`stage0/outputs/`)

| file | what it is |
|------|------------|
| `papulacandin_serum_pairs.csv` | **24 anchor compounds**, 143 matched serum-free/serum MIC pairs — the legitimate serum-shift labels |
| `unified_serum_shift_table.csv` | the Stage-0 endpoint table (papulacandin pairs today; schema ready for more) |
| `echinocandin_serum_context.csv` | teacher corpus: serum-context MIC rows, **UNPAIRED** (see note) |
| `echinocandin_free_fraction_seed.csv` | Fu / PPB rows — seed for the Stage-1b free-fraction oracle |
| `applicability_domain.csv` | per-dataset coverage (compounds, organisms, shift range, censoring, refs) |
| `stage0_summary.md` | human-readable report |

## Key columns in the endpoint table

- `serum_shift`, `log10_serum_shift` — the endpoint.
- `shift_relation` — `=` exact, `>=` right-censored (serum MIC hit the assay
  ceiling, so the shift is a **lower bound**), `<=` / `~` ambiguous. **62 of 143
  pairs are `>=`** and must be modelled as censored, not as point values.
- `intrinsic_potency_free_mic_uM`, `intrinsic_potency_pMIC` — potency **covariate**
  (guardrail #1: raw serum MIC is dominated by potency, so the shift is always
  analysed controlling for it). In this data Spearman(potency, shift) ≈ 0.72,
  largely because censored serum MICs pin at the ceiling — direct confirmation
  that both the covariate and the censoring handling are necessary.
- `match_level` / `pairing_confidence` — all 143 pairs are strain-matched within
  a single study (`high`).

## Data-reality finding: no echinocandin serum pairs

START_HERE anticipated pooling papulacandin **and** echinocandin serum shifts.
Building this stage surfaced that **the external echinocandin table has no
within-study matched serum pairs**: serum-containing and serum-free MICs never
share a source study, and pooling across studies mixes wild-type with
FKS-resistant mutants (one drug's free MIC spans 0.0001–512 µg/mL). Any ratio
from that is a cross-study, resistance-confounded artifact.

Per the project's own guardrails (model the shift honestly; *DATA* is the
binding constraint), we **do not fabricate echinocandin serum shifts**. The
echinocandin serum data is instead staged **unpaired**, in the right shape for a
Stage-1 hierarchical / transfer model that carries `study` and `serum_state` as
effects — which is the statistically correct way to borrow that signal.

## Next (Stage 1)

- **1a** potency oracle (structure-based + pooled-MIC QSAR).
- **1b** free-fraction / PPB oracle — start from `echinocandin_free_fraction_seed.csv`
  plus public ADMET protein-binding data (the most data-rich, highest-value piece).
- **1c** transfer-learning serum-shift oracle: pretrain on the echinocandin
  serum-context corpus, fine-tune on the 24 papulacandin pairs; **validate on a
  scaffold-split hold-out before any generation (the gate).**
