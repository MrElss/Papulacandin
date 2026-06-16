# 外部 FKS 抑制剂库整理流程

## 1. 每个化合物的最小录入标准

一个外部化合物至少要有：

1. `preferred_name`
2. `compound_family`
3. `target_relevance_tier`
4. `smiles_raw` 或结构文件来源
5. `structure_source`
6. 至少一个 `source_id` 或文献来源
7. `curation_status`

如果没有活性值，也可以先作为 structure/context row 放入 `external_compounds.csv`，
但不能进入 activity modeling。

## 2. 每条活性记录的最小录入标准

一条外部活性至少要有：

1. `external_compound_id`
2. `endpoint_type`
3. `endpoint_value`
4. `endpoint_relation`
5. `unit`
6. `assay_type`
7. `target_confidence`
8. `mechanism_confidence`
9. `source_id`

如果是 whole-cell MIC，必须写 `organism` 和尽量写 `strain`。

如果是 enzyme/binding assay，必须写 `target_name`、`target_gene`、`target_organism`。

## 3. 证据等级

### target_relevance_tier

| 值 | 含义 |
| --- | --- |
| Tier1_direct_fks | 直接 FKS / beta-1,3-glucan synthase 抑制剂 |
| Tier2_cell_wall | 真菌细胞壁相关，但不一定直接 FKS |
| Tier3_broad_antifungal | 广义抗真菌活性 |
| exclude_or_context_only | 暂不用于建模，只作背景 |

### target_confidence

| 值 | 含义 |
| --- | --- |
| direct | 有直接 enzyme/binding/genetic target 证据 |
| literature_supported | 文献机制明确，但本行不是直接靶点 assay |
| database_annotation | 数据库标注，未人工复核机制 |
| inferred | 从家族/药物类别推断 |
| unknown | 机制未知 |

### activity_relevance_tier

| 值 | 含义 |
| --- | --- |
| direct_fks | 直接靶点活性 |
| whole_cell_fks_drug | 已知 FKS 药物的 whole-cell 活性 |
| cell_wall_related | 细胞壁相关活性 |
| broad_antifungal | 广义抗真菌活性 |
| context_only | 不进入训练 |

## 4. 建模用途判断

### 可用于 pretraining

满足以下任一：

1. 结构可靠，有抗真菌活性。
2. 结构可靠，有 FKS/细胞壁相关机制。
3. 是高质量数据库或文献来源的 broad antifungal row。

### 可用于 fine-tuning

必须更严格：

1. Tier 1 direct FKS 或高度接近的 glucan-synthase inhibitor。
2. 活性 endpoint 和 assay context 清楚。
3. 结构可靠。
4. 不与 Papulacandin exact regression 直接混合，除非明确作为 external validation。

### 不用于 Papulacandin exact regression

默认所有外部数据都是 `use_for_papu_exact_regression = no`。

原因通常写：

```text
external_non_papulacandin_scaffold;different_assay_context;pretraining_only
```

## 5. 推荐工作流

1. 在 ChEMBL/PubChem/BindingDB/PubMed 找一个化合物或一个家族。
2. 把原始下载放到 `raw_downloads/`。
3. 在 `source_notes/` 写一条检索记录。
4. 把候选化合物填入 `external_compounds.csv`。
5. 把可靠活性填入 `external_activity_table.csv`。
6. 填 `external_inclusion_decisions.csv`，先决定用途，再建模。
7. 后续再用 RDKit 标准化结构，生成 processed 数据。

