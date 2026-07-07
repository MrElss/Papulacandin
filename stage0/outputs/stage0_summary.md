# Stage 0 — unify data & fix the serum-tolerance endpoint

Endpoint = whether a compound **keeps activity in serum**. Two faces: a categorical call (`retained` / `lost`) plus a continuous serum shift where meaningful. Intrinsic potency (serum-free MIC) is kept as a covariate.

## Curation rules applied

- **Activity ceiling (point 1):** MIC > 100 ug/mL = no detectable activity (a categorical INACTIVE call, not a censored number to model). So a compound active serum-free but above the ceiling in serum has simply **lost** serum activity.
- **In-vivo efficacy proxy (point 2):** in-vivo efficacy implies activity in the bloodstream (serum present) -> a **presumed serum-tolerant** positive. Whole-cell / cellular activity alone is NOT counted (losing serum activity despite a good broth MIC is this class's failure mode).

## Papulacandin direct labels (matched pairs)

- 24 compounds, 143 matched pairs. Categorical calls: **61 retained, 74 lost**, 0 ambiguous, 8 uninformative (inactive serum-free).
- The 'lost' pairs are the class's signature failure (potent in broth, dead in serum). Continuous serum shift is retained only where the serum MIC is a real active number (61 pairs).

## In-vivo efficacy proxies (presumed positives)

- 52 compounds flagged presumed serum-tolerant from in-vivo efficacy: 47 papulacandin, 5 canonical echinocandin (the caspofungin / anidulafungin / micafungin anchors). Labelled `retained_presumed`, confidence medium.

## Unified label set: 187 rows -> serum_tolerance_labels.csv

## Applicability domain

- **papulacandin_curated_matched_pairs** (direct_serum_tolerance_labels): 143 rows / 24 compounds; 3 refs. within-study matched serum-free/serum MIC; ceiling rule applied.
- **invivo_efficacy_proxies** (presumed_positive_labels): 52 rows / 49 compounds; 9 refs. in-vivo efficacy -> presumed serum-tolerant; cellular-only excluded.
- **echinocandin_external_serum_context** (teacher_corpus_unpaired): 6238 rows / 185 compounds; 2944 refs. 121 serum-present, 1126 canonical-FKS; UNPAIRED.
- **echinocandin_external_free_fraction** (free_fraction_oracle_seed): 9 rows / 1 compounds; 9 refs. Fu + PPB rows; seed for Stage-1b free-fraction oracle.

## Guardrails honored
1. Endpoint is serum TOLERANCE (retained/lost), not raw serum MIC.
2. Intrinsic potency retained as a covariate.
3. Ceiling rule: >100 ug/mL = inactive (point 1); no chasing numbers above it.
4. In-vivo efficacy is a flagged proxy; cellular-only is not a serum label (point 2).
5. No fabricated cross-study echinocandin pairs — DATA is the binding constraint.
