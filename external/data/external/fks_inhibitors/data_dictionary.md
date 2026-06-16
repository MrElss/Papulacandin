# External FKS Inhibitor Data Dictionary

## `external_compounds.csv`

| Field | Meaning |
| --- | --- |
| `external_compound_id` | Stable ID, recommended format `EXT-FKS-0001`. |
| `preferred_name` | Main compound name, e.g. caspofungin, micafungin, ibrexafungerp. |
| `synonyms` | Other names separated by semicolon. |
| `compound_family` | e.g. echinocandin, pneumocandin, triterpenoid, chitin_synthase_inhibitor, broad_antifungal. |
| `chemical_class` | More structural description, e.g. lipopeptide, triterpenoid glycoside. |
| `source_type` | approved_drug, natural_product, semisynthetic, synthetic, database_only. |
| `mechanism_class` | fks_glucan_synthase, cell_wall_related, broad_antifungal_unknown_target. |
| `target_relevance_tier` | Tier 1, Tier 2, or Tier 3. |
| `smiles_raw` | Source SMILES before standardization. |
| `smiles_canonical` | RDKit canonical isomeric SMILES after standardization. |
| `inchikey` | InChIKey for deduplication. |
| descriptor fields | RDKit descriptors after standardization. |
| `structure_source` | ChEMBL, PubChem, BindingDB, paper_SI, supplier, manual. |
| `structure_confidence` | high, medium, low. |
| `source_ids` | External identifiers, e.g. ChEMBL ID, PubChem CID. |
| `curation_status` | planned, imported, curated, needs_review, excluded. |

## `external_activity_table.csv`

| Field | Meaning |
| --- | --- |
| `external_activity_id` | Stable row ID, recommended format `EXT-ACT-000001`. |
| `external_compound_id` | Link to `external_compounds.csv`. |
| `endpoint_type` | MIC, MIC50, IC50, EC50, Ki, inhibition_percent, in_vivo_ED50, etc. |
| `endpoint_value`, `endpoint_relation`, `unit` | Source value exactly as curated. Preserve `>`, `<`, `>=`, `<=`. |
| `converted_value_uM` | Only when reliable molecular weight and unit conversion are available. |
| `p_activity` | Optional `-log10(M)` potency. |
| `organism`, `strain` | Whole-cell assay organism and strain. |
| `target_name`, `target_gene`, `target_organism` | Enzyme/target context. |
| `assay_type` | biochemical, whole_cell, cell_wall_readout, in_vivo, binding, other. |
| `assay_context` | no_serum, serum_present, unknown_serum, enzyme_assay, in_vivo, etc. |
| `target_confidence` | direct, inferred, broad_cell_wall, unknown. |
| `mechanism_confidence` | direct_assay, literature_mechanism, database_annotation, unknown. |
| `activity_relevance_tier` | direct_fks, cell_wall_related, broad_antifungal, context_only. |

## `external_inclusion_decisions.csv`

| Field | Meaning |
| --- | --- |
| `include_in_tier1_fks_set` | yes/no. |
| `include_in_tier2_cell_wall_set` | yes/no. |
| `include_in_tier3_broad_antifungal_set` | yes/no. |
| `use_for_pretraining` | yes/no. |
| `use_for_finetuning` | yes/no; only high-confidence close-task rows. |
| `use_for_papu_exact_regression` | almost always no. |
| `do_not_mix_with_papu_reason` | Required when exact regression is no. |

