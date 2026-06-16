# Papulacandin main-library data-only handoff

This package contains the curated Papulacandin main-library data and structure files only.

It intentionally excludes processed datasets, modeling results, trained models, scripts, reports, and candidate-ranking outputs. The purpose is to let another tool independently inspect and calculate from the curated source tables and structures.

## Suggested reading order

1. `core_tables/compounds_master.csv`
2. `core_tables/activity_table.csv`
3. `core_tables/enzyme_assays.csv`
4. `core_tables/structure_quality_notes.tsv`
5. `core_tables/structure_qc_manual_review.csv`
6. `core_tables/synthesis_feasibility.csv`
7. `structures/sdf/`
8. `structures/mol/`
9. `structures/cdx/`
10. `after_analysis_reference/sar_annotations.csv`

## Core tables

- `compounds_master.csv`: one row per Papulacandin-family compound.
- `activity_table.csv`: curated whole-cell / activity records.
- `enzyme_assays.csv`: enzyme-level assay records, important for enzyme-cell and exposure-translation analysis.
- `structure_quality_notes.tsv`: structure-source and stereochemistry notes.
- `structure_qc_manual_review.csv`: manual structure-QC checklist.
- `synthesis_feasibility.csv`: synthesis and follow-up feasibility information.

## Structure files

- `structures/cdx/`: ChemDraw source files.
- `structures/mol/`: MOL files exported from ChemDraw/RDKit workflow.
- `structures/sdf/`: SDF files for cheminformatics parsing.

Only `cdx`, `mol`, and `sdf` are included to avoid confusing the receiving tool with multiple alternative structure exports such as V2000 folders or combined files.

## SAR annotation boundary

`after_analysis_reference/sar_annotations.csv` is included as an after-analysis reference only.

Recommended instruction to the receiving tool:

First analyze SAR from `activity_table.csv`, `compounds_master.csv`, and structure files independently. Then compare its conclusions against `after_analysis_reference/sar_annotations.csv`.

This avoids using existing human/literature SAR annotations as the starting assumption.

## Excluded

- `data/processed/`
- `results/`
- trained models
- scripts
- reports and PPT files
- external FKS datasets
- raw literature PDFs
- backup folders
- alternative structure export folders such as `combined`, `mol_v2000`, and `sdf_v2000`
