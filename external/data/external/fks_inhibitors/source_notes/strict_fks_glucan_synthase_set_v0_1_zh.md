# Strict FKS / beta-1,3-glucan-synthase candidate set v0.1（中文）

生成时间：2026-06-11 14:12

## 1. 本轮边界

本轮不再混入三类候选。只保留以下对象：

1. 有直接 `FKS`、`1,3-beta-glucan synthase`、`beta-1,3-glucan synthase` 抑制证据的化合物。
2. 属于明确 FKS/glucan-synthase inhibitor 家族的候选化合物，例如 echinocandin、pneumocandin、enfumafungin、aculeacin、mulundocandin。

不纳入：

- azole/polyene/allylamine/antimetabolite 等广谱抗真菌背景药，除非有直接 FKS/glucan-synthase 抑制证据。
- FKS 突变菌株 MIC 背景但没有直接酶/靶点抑制证据的化合物。
- UDP-glucose 等底物、辅因子或 assay reagent。

## 2. 数量

| 项目 | 数量 |
| --- | ---: |
| 原始 v0.2 候选化合物 | 500 |
| strict 纳入化合物 | 321 |
| strict 排除化合物 | 179 |
| strict 纳入化合物的全部候选活性行 | 3113 |
| 其中 direct FKS/glucan-synthase 抑制证据行 | 468 |
| 使用到的 source/cache 记录 | 215 |
| strict 集合有 SMILES 的化合物 | 321 |
| strict 集合有 InChIKey 的化合物 | 321 |

## 3. 纳入类别

| strict_inclusion_class | n_compounds |
| --- | --- |
| direct_FKS_or_beta_1_3_glucan_synthase_inhibitor_evidence | 305 |
| known_FKS_or_glucan_synthase_inhibitor_family_candidate | 16 |

## 4. 排除原因

| strict_exclusion_reason | n_compounds |
| --- | --- |
| no_direct_FKS_or_glucan_synthase_inhibition_evidence_and_not_known_FKS_family | 153 |
| broad_antifungal_background_not_FKS_or_glucan_synthase_inhibitor | 25 |
| substrate_or_assay_reagent_not_inhibitor | 1 |

## 5. direct evidence target 预览

| target_name | n_rows |
| --- | --- |
| Beta-1,3-glucan synthase | 249 |
| Unchecked | 62 |
| 1,3-beta-glucan synthase | 52 |
| Glucan synthase | 43 |
| 1,3-beta-glucan synthase component GSC2 | 29 |
| 1,3-beta-D-glucan synthase subunit 1 | 15 |
| 1,3-beta-glucan synthase component FKS1 | 12 |
| Fksp | 3 |
| 1,3-beta-D-glucan synthase subunit | 2 |
| Candida albicans | 1 |

## 6. 人工复核建议

先打开：

`data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_review_queue_candidate_v0_1.csv`

优先复核 `P1_verify_direct_inhibition_evidence`。复核时不要只看 compound 名称，要打开对应 direct evidence 行：

`data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_direct_evidence_candidate_v0_1.csv`

确认 `assay_method` 或 `target_name` 是否真的是 FKS / beta-1,3-glucan synthase inhibition。只有文献或数据库明确写 FKS1 时，才标成 FKS1-specific；否则标成 `FKS_or_beta_1_3_glucan_synthase_inhibitor`。

## 7. 输出文件

| 文件 | 用途 |
| --- | --- |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_compounds_candidate_v0_1.csv` | strict 化合物候选 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_activity_table_candidate_v0_1.csv` | strict 化合物的全部候选活性行 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_direct_evidence_candidate_v0_1.csv` | direct FKS/glucan-synthase 抑制证据行 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_exclusions_candidate_v0_1.csv` | 从 500 个候选中排除的化合物和原因 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_review_queue_candidate_v0_1.csv` | 人工复核队列 |
