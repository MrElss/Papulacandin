# Stage 0 — unify data and fix the endpoint

Endpoint: **serum shift = MIC(serum-containing) / MIC(serum-free)** for the *same compound in the same organism* — a potency-independent measure. Intrinsic potency (the serum-free MIC) is carried as a covariate.

## Papulacandin (curated) — the serum-shift labels

- **24 compounds** carry matched serum-free + serum-containing MIC — the ~24 anchor labels referenced in START_HERE.
- **143 matched pairs** total (a compound recurs across organisms / strains); 143 are strain-matched within a single study (high confidence).
- **70/143 pairs are right-censored** (serum MIC reported as ">"), so their shift is a LOWER BOUND, flagged in `shift_relation` (`>=`). They must be modelled as censored, never as point values.

## Echinocandin (external) — teacher corpus, NOT fabricated pairs

**Key data-reality finding:** the external echinocandin table contains **no within-study matched serum pairs**. Serum-containing and serum-free MICs never share a source study, and pooling across studies mixes wild-type with FKS-resistant mutants (one drug's free MIC spans 0.0001–512 ug/mL). A ratio built from that is a cross-study, resistance-confounded artifact. Honoring the project's guardrails, we do **not** manufacture echinocandin serum shifts. Instead:

- `echinocandin_serum_context.csv`: 6238 whole-cell MIC rows (121 explicitly serum-present, 1125 canonical-FKS), kept UNPAIRED with their source study — the corpus for a Stage-1 hierarchical / transfer model with `study` and `serum_state` as effects.
- `echinocandin_free_fraction_seed.csv`: 9 Fu / PPB rows — the seed for the Stage-1b free-fraction oracle.

## Applicability domain

- **papulacandin_curated_matched_pairs** (serum_shift_labels): 143 rows / 24 compounds; shift 1.0–500.0x (median 16.03x); 3 refs. legitimate within-study matched serum-free/serum MIC.
- **echinocandin_external_serum_context** (teacher_corpus_unpaired): 6238 rows / 185 compounds; 2944 refs. 121 serum-present rows; 1125 canonical-FKS rows; UNPAIRED (no within-study serum pairs exist).
- **echinocandin_external_free_fraction** (free_fraction_oracle_seed): 9 rows / 1 compounds; 9 refs. Fu + PPB rows; seed for Stage-1b free-fraction oracle.

## Guardrails honored
1. Endpoint is the serum SHIFT, not the raw serum MIC (guardrail #1).
2. Intrinsic potency retained as a covariate (dominant confound).
3. Censoring tracked explicitly (`shift_relation`, `serum_censored`).
4. Applicability domain recorded per dataset.
5. No fabricated cross-study pairs — DATA is the binding constraint (guardrail #8).
