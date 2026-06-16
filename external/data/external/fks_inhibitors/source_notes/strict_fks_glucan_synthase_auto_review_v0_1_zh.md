# Strict FKS/glucan-synthase automatic review v0.1（中文）

生成时间：2026-06-11 14:23

## 1. 这次复核做了什么

本轮没有逐篇阅读原文，而是基于 ChEMBL 已导出的字段做数据库字段级复核：

- `target_name`
- `assay_method`
- `endpoint_type`
- `source_reference`
- `compound_family`

判断目标是：把 strict 候选化合物分成自动接受、家族候选但需回源、以及需要进一步确认几类。它不是最终文献级审定，但已经可以替你完成第一轮大批量复核。

## 2. 化合物级复核结果

| auto_review_decision | n_compounds |
| --- | --- |
| auto_accept_direct_FKS_or_glucan_synthase | 305 |
| auto_accept_known_FKS_family_need_source_confirmation | 16 |

## 3. 置信度分布

| auto_review_confidence | n_compounds |
| --- | --- |
| high | 265 |
| medium_high | 40 |
| medium | 16 |

## 4. FKS1 specificity 自动判断

| auto_fks1_specificity | n_compounds |
| --- | --- |
| FKS_or_beta_1_3_glucan_synthase_inhibitor | 318 |
| FKS1_specific_candidate | 3 |

说明：只有 `assay_method` 或 `target_name` 明确出现 FKS1/subunit 1/component FKS1 时，才标为 `FKS1_specific_candidate`。其余统一为 `FKS_or_beta_1_3_glucan_synthase_inhibitor`。

## 5. direct evidence 行级置信度

| direct_row_confidence | n_rows |
| --- | --- |
| high | 366 |
| medium_high | 102 |

## 6. 仍建议抽查的对象

### 6.1 target 字段为 Unchecked 但 assay_method 支持 direct inhibition

这些不是错误；很多 ChEMBL 行 target 字段较粗，但 assay_method 写得很清楚。建议作为高优先级候选时抽查原文。

| external_compound_id | preferred_name | accepted_direct_evidence_rows | target_unchecked_direct_rows | auto_review_confidence |
| --- | --- | --- | --- | --- |
| EXT-FKS-0001 | CASPOFUNGIN | 53 | 16 | high |
| EXT-FKS-0002 | ANIDULAFUNGIN | 52 | 18 | high |
| EXT-FKS-0003 | MICAFUNGIN | 47 | 17 | high |
| EXT-FKS-0017 | ACULEACIN A | 2 | 1 | high |
| EXT-FKS-0022 | CHEMBL3037784 | 2 | 1 | high |
| EXT-FKS-0023 | ENFUMAFUNGIN | 3 | 3 | medium_high |
| EXT-FKS-0194 | CHALCONE | 1 | 1 | medium_high |
| EXT-FKS-0409 | CHEMBL4751360 | 1 | 1 | medium_high |
| EXT-FKS-0410 | CHEMBL4758391 | 1 | 1 | medium_high |
| EXT-FKS-0411 | CHEMBL4759855 | 1 | 1 | medium_high |
| EXT-FKS-0412 | CHEMBL4800505 | 1 | 1 | medium_high |
| EXT-FKS-0028 | IBREXAFUNGERP | 1 | 1 | medium_high |

### 6.2 已知 FKS 家族但本表没有 direct evidence 的候选

这些可以保留为 family candidate，但用于最终训练前最好回源确认。

| external_compound_id | preferred_name | compound_family | recommended_action |
| --- | --- | --- | --- |
| EXT-FKS-0493 | CHEMBL1160279 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0494 | CHEMBL1160280 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0214 | CHEMBL1790210 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0413 | CHEMBL1790214 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0454 | CHEMBL2369955 | mulundocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0497 | CHEMBL266844 | pneumocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0498 | CHEMBL267318 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0016 | CHEMBL3760049 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0500 | CHEMBL405465 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0218 | CHEMBL412839 | pneumocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0219 | CHEMBL412976 | pneumocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0416 | CHEMBL452101 | pneumocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0417 | CHEMBL459922 | echinocandin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0020 | CHEMBL4764856 | enfumafungin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0191 | CHEMBL5414434 | enfumafungin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |
| EXT-FKS-0192 | CHEMBL5426159 | enfumafungin | keep_as_family_candidate; source_check_needed_for_direct_assay_or mechanism statement |

## 7. 输出文件

| 文件 | 用途 |
| --- | --- |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_auto_review_v0_1.csv` | 化合物级自动复核总表 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_auto_accepted_compounds_v0_1.csv` | 自动接受候选 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_auto_hold_compounds_v0_1.csv` | 暂缓/需进一步确认候选 |
| `data\external\fks_inhibitors\source_exports\strict_fks_glucan_synthase_direct_evidence_reviewed_v0_1.csv` | direct evidence 行级复核 |

## 8. 使用建议

你下一步可以先接受 `auto_accept_direct_FKS_or_glucan_synthase` 作为 strict 外部机制库的第一版；对 `medium_high` 和 `medium` 置信度对象做抽查即可。真正需要读原文的是家族候选、target unchecked 但非常关键的药物、以及后续准备进入训练/论文汇报的代表性行。
