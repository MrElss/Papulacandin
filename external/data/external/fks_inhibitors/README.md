# External FKS / Glucan-Synthase Inhibitor Database

Purpose: build an auxiliary database of non-Papulacandin FKS1/FKS2 or
beta-1,3-glucan-synthase inhibitors for transfer learning, mechanism comparison,
and chemical-space context.

Important boundary:

- Do not merge these rows into `data/curated/compounds_master.csv`.
- Do not train Papulacandin exact MIC regression by simply pooling external rows.
- Use external rows for pretraining, structure-space comparison, broad classification,
  and target/mechanism context.

Folder layout:

```text
raw_downloads/      # original downloaded files from databases or papers
source_exports/     # exported CSV/SDF/TSV files from ChEMBL, PubChem, BindingDB, etc.
source_notes/       # one note per source search or paper
structures/         # imported SDF/MOL/SMILES files for external compounds
processed/          # generated standardized external tables, not manually edited
```

Template tables:

```text
external_compounds.csv
external_activity_table.csv
external_targets.csv
external_sources.csv
external_structure_qc.csv
external_inclusion_decisions.csv
```

Guides:

```text
external_dataset_manifest.md
data_dictionary.md
search_strategy_zh.md
curation_protocol_zh.md
```

Automated candidate fetch:

```powershell
# Conservative first pass with all seed queries.
.venv-rdkit\Scripts\python.exe scripts\fetch_external_fks_inhibitors.py `
  --chembl-activity-limit 200 `
  --chembl-assay-search-limit 3 `
  --chembl-assay-activity-limit 150 `
  --pubchem-max-cids 1 `
  --pubchem-assay-row-limit 75 `
  --bindingdb-max-compounds 40
```

Generated candidate files are written to `source_exports/*candidate_v0_1.csv`
and raw JSON responses are cached under `raw_downloads/auto_fetch_external_fks_v0_1/`.
These candidate files are review inputs; do not treat them as curated external
tables until assay, structure, source, and endpoint fields have been checked.

Expanded v0.2 review set:

```powershell
# Rebuild the 500-compound review set from cached ChEMBL JSON.
.venv-rdkit\Scripts\python.exe scripts\build_external_fks_candidates_from_cache_v0_2.py `
  --output-version v0_2 `
  --max-output-compounds 500

# Build the human-review priority queue.
.venv-rdkit\Scripts\python.exe scripts\build_external_fks_review_queue_v0_2.py
```

Start manual review from
`source_exports/external_review_queue_candidate_v0_2.csv`, not directly from the
full activity table. The v0.2 files include broad antifungal/background rows as
candidate pretraining material; only rows promoted after review should enter the
formal external curated tables.

Strict FKS/glucan-synthase subset:

```powershell
.venv-rdkit\Scripts\python.exe scripts\build_strict_fks_glucan_synthase_set_v0_1.py
```

This creates a narrower review set under
`source_exports/strict_fks_glucan_synthase_*_candidate_v0_1.csv`. Use this set
when the current round should include only FKS1/FKS or beta-1,3-glucan-synthase
inhibitor candidates. Broad antifungal background compounds are excluded from
this strict set unless they have direct FKS/glucan-synthase inhibition evidence.

Database-field automatic review:

```powershell
python scripts\review_strict_fks_glucan_synthase_candidates_v0_1.py
```

This creates `source_exports/strict_fks_glucan_synthase_auto_*_v0_1.csv`.
The review is based on exported database fields such as target, assay method,
endpoint, source reference, and compound family. It is useful for first-pass
candidate acceptance, but it is not a replacement for source-paper confirmation
when a compound is used in a final claim or benchmark.

External model-ready package:

```powershell
python scripts\build_external_fks_model_ready_v0_1.py
```

This creates `data/processed/external_fks_model_ready_v0_1/`, including a
compound-level model-ready table, RDKit 2D descriptor matrix, and ECFP4/Morgan
fingerprint matrix. These outputs are for external pretraining, mechanism
classification, and chemical-space comparison. They should not be pooled into
the curated Papulacandin exact activity regression target.

Papulacandin versus external FKS chemical-space comparison:

```powershell
python scripts\compare_papulacandin_external_fks_space_v0_1.py
```

This creates a combined descriptor-space table, Papulacandin-to-external nearest
neighbor table, and a Chinese report under `results/external_fks_inhibitors/`.
Use this as the bridge between the curated Papulacandin core project and the
external FKS pretraining/background set.

FKS-like pretraining dataset and diagnostic model:

```powershell
python scripts\build_fks_pretraining_dataset_v0_1.py
python scripts\train_fks_pretraining_model_v0_1.py
```

This creates `data/processed/pretraining_v0_1/`, ChEMBL background-decoy exports,
and `results/modeling/fks_pretraining_*_v0_1.*`. The negative/background label
in this package is a weak background/decoy label, not an experimentally proven
inactive label.
