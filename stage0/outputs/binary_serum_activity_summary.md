# Stage 0 (binary reframe) — serum-active / inactive label set

Endpoint reframed to BINARY `serum_active` (retains activity in serum: yes/no) to pool more, and more diverse, data than the quantitative serum-shift endpoint could — the right call for a discovery goal, and robust to the fact that serum inactivation is not fully computable (planning/serum_mechanism_evidence.md).

## Trustworthy training set (use this)

- **31 compounds with usable labels** (20 serum-active, 11 inactive) — direct serum-MIC evidence plus echinocandin in-vivo/clinical positives.
- Direct serum-MIC compounds: 30 (24 papulacandin, 6 echinocandin).
- vs the **24** matched-pair compounds the quantitative endpoint could use — a modest but real gain, and now the endpoint is the discovery-relevant binary form.

## Caveated (NOT training labels)

- **43 papulacandin compounds known only from an in-vivo proxy**: serum activity is UNKNOWN for this class (in-vivo ED50s coexist with serum loss — papulacandin B is the archetype). These are **test candidates** for the wet-lab panel, not positive training labels — pooling them would inject ~43 likely-false positives.
- 10 compounds are mixed (active in some serum assays, lost in others).

## Where the real expansion comes from: ChEMBL

- Only 6 echinocandins have direct serum-MIC evidence in-repo. The binary endpoint can absorb many more echinocandin/FKS serum-shift records — that is the high-value ChEMBL pull (see drop-in below), the path to a genuinely larger serum-labelled set.
- 313 total observations. By source: papulacandin_serum_mic 143, echinocandin_serum_mic 115, in_vivo_proxy 52, clinical_approved 3.

## ChEMBL drop-in

- ChEMBL echinocandin serum records: not present yet. Place a CSV at `stage0/data/chembl_echinocandin_serum.csv` (columns: compound_name, serum_active, source_ref) and re-run to merge them automatically. (Pending the ChEMBL connector.)

## Use

`binary_serum_activity_labels.csv` is the pooled training/target set for a serum-active classifier and for panel selection. Keep the potency oracle (Stage 1a) as the filter and control potency as a covariate — the binary label is still potency-correlated within a series. `evidence_level` and `confidence` let you weight direct serum MICs above in-vivo/clinical proxies.
