# External FKS data-only clean handoff

This package contains only non-empty data files and data documentation for the external FKS / beta-1,3-glucan synthase inhibitor side of the Papulacandin/FKS1 project.

No trained models, modeling result reports, analysis conclusions, or pipeline scripts are included.

## Why this clean package exists

The project folder contains several scaffold CSV templates such as:

- `external_sources.csv`
- `external_structure_qc.csv`
- `external_compounds.csv`
- `external_activity_table.csv`
- `external_targets.csv`

In the current repository snapshot, those top-level scaffold files contain headers only and have 0 data rows. They are useful as schema templates, but they are distracting for an independent re-analysis. Therefore they are excluded from this clean handoff.

The actual external FKS data are in:

- `data/external/fks_inhibitors/source_exports/`
- `data/external/fks_inhibitors/source_notes/`
- `data/processed/external_fks_model_ready_v0_1/`
- `data/processed/pretraining_v0_1/`

## Included

### Documentation

- `data/external/fks_inhibitors/README.md`
- `data/external/fks_inhibitors/external_dataset_manifest.md`
- `data/external/fks_inhibitors/data_dictionary.md`
- `data/external/fks_inhibitors/curation_protocol_zh.md`
- `data/external/fks_inhibitors/search_strategy_zh.md`
- `data/external/fks_inhibitors/processed/README.md`
- `data/external/fks_inhibitors/structures/README.md`
- `data/external/fks_inhibitors/source_exports/README.md`
- `data/external/fks_inhibitors/source_notes/README.md`
- `data/external/fks_inhibitors/source_notes/*.md`

### Non-empty source-export CSV data

Included from `data/external/fks_inhibitors/source_exports/`.

These are the candidate exports, strict FKS/glucan-synthase reviewed sets, background/decoy data, source tables, activity tables, and inclusion decisions.

### Processed external data

Included from `data/processed/external_fks_model_ready_v0_1/`.

### Pretraining dataset data

Included from `data/processed/pretraining_v0_1/`.

## Excluded

- all 0-row scaffold CSV files
- `results/`
- `results/models/`
- trained `.joblib` model files
- modeling reports and metrics
- `scripts/`
- `data/external/fks_inhibitors/raw_downloads/`
- temporary outputs
- Python virtual environments
- `node_modules`

## Suggested use

Ask the receiving tool to independently:

1. inspect the external source-export data,
2. reconstruct inclusion/exclusion decisions,
3. rebuild its own external FKS-positive/background dataset,
4. recalculate descriptors/fingerprints if desired,
5. compare the recalculated matrices with the included processed matrices.
