# Stage 0 — unify data & fix the endpoint

First step of the methodology in [`../START_HERE.md`](../START_HERE.md) (§3).
Its only job is to turn the raw activity records into **one well-defined
serum-tolerance endpoint both chemotypes can share**, carrying the confounders
every downstream oracle needs. **No modeling happens here.**

The endpoint has two faces:

```
serum_tolerance = retained | lost      # categorical headline label
serum_shift     = MIC(serum) / MIC(serum-free)   # continuous, where meaningful
```

## Run

```bash
python stage0/build_serum_shift_table.py     # pure stdlib, no deps
python -m pytest stage0/test_stage0.py -q     # 8 smoke tests
```

## Two curation rules

- **Activity ceiling (point 1).** A MIC **> 100 µg/mL = no detectable activity** —
  a categorical *inactive* call, not a censored number to chase or model
  numerically. So a compound active serum-free but above the ceiling in serum has
  simply **lost** serum activity. Implemented in `activity_call()`.
- **In-vivo efficacy proxy (point 2).** In-vivo efficacy (curing/reducing an
  animal infection) requires activity in the bloodstream, where serum proteins
  are present → **presumed serum-tolerant** positive, flagged and lower-confidence
  than a direct serum MIC. **Whole-cell / cellular activity alone is NOT counted**
  — losing serum activity despite a good broth MIC is this class's failure mode.

## Outputs (`stage0/outputs/`)

| file | what it is |
|------|------------|
| `papulacandin_serum_pairs.csv` | 24 compounds, 143 matched serum-free/serum MIC pairs with `serum_tolerance` calls |
| `serum_tolerance_labels.csv` | **unified label set**: 135 direct pair labels + 52 in-vivo proxies (187 rows) |
| `invivo_serum_tolerance_proxies.csv` | presumed-positive labels from in-vivo efficacy (both chemotypes) |
| `echinocandin_serum_context.csv` | teacher corpus: serum-context MIC rows, **unpaired** (see note) |
| `echinocandin_free_fraction_seed.csv` | Fu / PPB rows — seed for the Stage-1b free-fraction oracle |
| `applicability_domain.csv` | per-dataset coverage |
| `stage0_summary.md` | human-readable report |

## The labels

**Direct (papulacandin matched pairs):** 24 compounds → **61 retained, 74 lost**,
plus a few ambiguous/uninformative. All strain-matched within a single study.
`serum_shift` is kept only where the serum MIC is a real active number
(`serum_shift_meaningful = True`); potency is carried as a covariate
(`intrinsic_potency_pMIC`) per guardrail #1.

**In-vivo proxies:** 52 presumed positives, including the 5 canonical
echinocandin anchors (caspofungin, anidulafungin, micafungin, cilofungin…) —
exactly the serum-tolerant drugs we want as positives. Each proxy is
cross-referenced to the higher-confidence direct pairs via `direct_pair_evidence`;
when a compound also has a direct `lost`/`mixed` call the proxy is down-weighted
to `confidence = low`. **One such conflict exists — Papulacandin B** — which is
genuinely mixed (retained in some organisms, dead in 50% mouse serum) yet carries
an in-vivo ED50: the textbook caution that an in-vivo readout does not guarantee
serum tolerance. Prefer the direct label where both exist.

## Data-reality finding: no echinocandin serum pairs

The external echinocandin table has **no within-study matched serum pairs**:
serum-containing and serum-free MICs never share a source study, and pooling
across studies mixes wild-type with FKS-resistant mutants (one drug's free MIC
spans 0.0001–512 µg/mL). We therefore do **not** fabricate echinocandin serum
shifts. That data is staged **unpaired** for a Stage-1 hierarchical/transfer
model, and the echinocandins' serum signal is instead borrowed through their
in-vivo efficacy proxies (point 2).

## Binary reframe (serum-active / inactive)

`build_binary_serum_activity.py` reframes the endpoint to **binary
`serum_active`** (retains activity in serum: yes/no) — the discovery-relevant
form, and robust to the finding that serum inactivation isn't fully computable
(`planning/serum_mechanism_evidence.md`). It pools papulacandin serum MICs,
echinocandin serum MICs, in-vivo efficacy, and clinical use.

- **Trustworthy training set: 31 compounds** (20 active / 11 inactive) — direct
  serum-MIC evidence + echinocandin in-vivo/clinical positives.
- **43 papulacandin in-vivo-only compounds are held out as _test candidates_,
  not training labels** — an in-vivo ED50 does not prove serum activity for this
  class (papulacandin B has one yet is serum-lost), so pooling them would inject
  likely-false positives.
- Only 6 echinocandins have direct serum evidence in-repo → the **ChEMBL pull is
  the real expansion path**. Drop a CSV at `stage0/data/chembl_echinocandin_serum.csv`
  (`compound_name, serum_active, source_ref`) and re-run to merge it automatically.

Outputs: `binary_serum_activity_labels.csv`, `binary_serum_activity_observations.csv`,
`binary_serum_activity_summary.md`.

## Next (Stage 1)

- **1a** potency oracle (structure-based + pooled-MIC QSAR).
- **1b** free-fraction / PPB oracle — start from `echinocandin_free_fraction_seed.csv`
  plus public ADMET protein-binding data.
- **1c** transfer serum-tolerance oracle: pretrain on the echinocandin serum-context
  corpus + in-vivo proxies, fine-tune on the 24 papulacandin direct labels;
  **validate on a scaffold-split hold-out before any generation (the gate).**
