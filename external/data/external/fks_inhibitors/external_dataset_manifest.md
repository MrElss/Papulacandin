# External FKS Inhibitor Dataset Manifest

Version: v0.1 scaffold

## Scope

This dataset collects non-Papulacandin compounds relevant to FKS1/FKS2,
beta-1,3-glucan synthase, fungal cell-wall inhibition, or broad antifungal
pretraining.

The dataset is auxiliary. It supports pretraining, comparison, and design
context. It does not redefine the curated Papulacandin SAR database.

## Tiers

| Tier | Meaning | Typical examples | Recommended use |
| --- | --- | --- | --- |
| Tier 1 | Direct FKS / beta-1,3-glucan synthase inhibitors | echinocandins, pneumocandins, enfumafungin / ibrexafungerp, clearly annotated glucan-synthase inhibitors | mechanism comparison, transfer-learning fine-tuning candidate pool |
| Tier 2 | Fungal cell-wall related but not direct FKS | chitin synthase inhibitors, mannan/cell-wall stress related compounds | weakly related pretraining, contrastive biology |
| Tier 3 | Broad antifungal activity without direct target | general Candida/Aspergillus active molecules | broad pretraining only |

## Required Tables

| File | Purpose | Manual edit? |
| --- | --- | --- |
| `external_compounds.csv` | one row per external compound | yes |
| `external_activity_table.csv` | one row per activity measurement | yes |
| `external_targets.csv` | target metadata | yes |
| `external_sources.csv` | database/paper provenance | yes |
| `external_structure_qc.csv` | structure QC status | yes at v0.1, later generated |
| `external_inclusion_decisions.csv` | modeling-use decisions | yes |

## Inclusion Boundary

Include a compound if at least one is true:

1. It is reported as a beta-1,3-glucan synthase, FKS1, FKS2, or glucan-synthesis inhibitor.
2. It is a clinically or historically recognized echinocandin/pneumocandin/aculeacin/mulundocandin/enfumafungin-like compound.
3. It has curated antifungal activity and is explicitly cell-wall-related.
4. It is selected as a broad antifungal pretraining row and is clearly marked Tier 3.

Exclude or mark as low priority if:

1. Only cytotoxicity is reported.
2. Only antibacterial activity is reported.
3. The compound is a Papulacandin already curated in the PAPU database.
4. Target/mechanism is unknown and antifungal data are sparse.

## Modeling Boundary

External rows may be used for:

- pretraining,
- chemical-space comparison,
- broad active/inactive classification,
- target/mechanism contrast,
- external validation of descriptors.

External rows should not be used for:

- direct Papulacandin exact MIC regression,
- claiming Papulacandin-specific SAR,
- replacing source-level Papulacandin curation.

