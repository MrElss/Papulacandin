# External FKS candidate review queue v0.2（中文）

生成时间：2026-06-11 13:57

## 1. 用途

这个表用于人工复核 `external_*_candidate_v0_2.csv`。它不是新数据源，而是把 500 个候选化合物按复核优先级重新排队。

## 2. 复核优先级

| review_priority | n_compounds |
| --- | --- |
| P1_direct_fks_verify_first | 306 |
| P2_broad_antifungal_pretraining_review | 178 |
| P1_known_fks_family_verify_first | 16 |

## 3. 建议复核顺序

1. 先复核 `P1_direct_fks_verify_first`：确认是否真的是 FKS / beta-1,3-glucan synthase 直接抑制 assay。
2. 再复核 `P1_known_fks_family_verify_first`：确认结构和家族，例如 echinocandin、pneumocandin、enfumafungin 等。
3. `P2_broad_antifungal_pretraining_review` 主要作为预训练或背景抗真菌数据，不能默认当成 FKS 抑制剂。
4. `P3_context_or_structure_review` 暂不建议建模使用。

## 4. 前 30 个复核对象

| external_compound_id | preferred_name | compound_family | direct_fks | candidate_activity_rows | review_priority |
| --- | --- | --- | --- | --- | --- |
| EXT-FKS-0001 | CASPOFUNGIN | echinocandin | 53 | 817 | P1_direct_fks_verify_first |
| EXT-FKS-0002 | ANIDULAFUNGIN | echinocandin | 52 | 783 | P1_direct_fks_verify_first |
| EXT-FKS-0003 | MICAFUNGIN | echinocandin | 47 | 707 | P1_direct_fks_verify_first |
| EXT-FKS-0023 | ENFUMAFUNGIN | enfumafungin | 3 | 26 | P1_direct_fks_verify_first |
| EXT-FKS-0017 | ACULEACIN A | aculeacin | 2 | 32 | P1_direct_fks_verify_first |
| EXT-FKS-0022 | CHEMBL3037784 | echinocandin | 2 | 27 | P1_direct_fks_verify_first |
| EXT-FKS-0026 | URIDINE DIPHOSPHATE GLUCOSE | mixed_fks | 2 | 28 | P1_direct_fks_verify_first |
| EXT-FKS-0121 | CHEMBL1684737 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0122 | CHEMBL1684744 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0123 | CHEMBL1770510 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0131 | CHEMBL3693938 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0132 | CHEMBL3693939 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0133 | CHEMBL3693940 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0134 | CHEMBL3693941 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0135 | CHEMBL3693942 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0136 | CHEMBL3693943 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0137 | CHEMBL3693944 | mixed_fks | 2 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0028 | IBREXAFUNGERP | enfumafungin | 1 | 35 | P1_direct_fks_verify_first |
| EXT-FKS-0018 | CILOFUNGIN | echinocandin | 1 | 29 | P1_direct_fks_verify_first |
| EXT-FKS-0027 | CHEMBL2387203 | pneumocandin | 1 | 26 | P1_direct_fks_verify_first |
| EXT-FKS-0044 | Echinocandin B | echinocandin | 1 | 13 | P1_direct_fks_verify_first |
| EXT-FKS-0035 | CHEMBL437803 | pneumocandin | 1 | 12 | P1_direct_fks_verify_first |
| EXT-FKS-0050 | CHEMBL267988 | pneumocandin | 1 | 3 | P1_direct_fks_verify_first |
| EXT-FKS-0053 | CHEMBL407395 | pneumocandin | 1 | 3 | P1_direct_fks_verify_first |
| EXT-FKS-0054 | CHEMBL412195 | pneumocandin | 1 | 3 | P1_direct_fks_verify_first |
| EXT-FKS-0055 | CHEMBL414690 | pneumocandin | 1 | 3 | P1_direct_fks_verify_first |
| EXT-FKS-0072 | CHEMBL261821 | pneumocandin | 1 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0073 | CHEMBL262998 | pneumocandin | 1 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0074 | CHEMBL263952 | pneumocandin | 1 | 2 | P1_direct_fks_verify_first |
| EXT-FKS-0075 | CHEMBL264601 | pneumocandin | 1 | 2 | P1_direct_fks_verify_first |

仅显示前 30 行；完整表见 CSV。

## 5. 输出文件

`data\external\fks_inhibitors\source_exports\external_review_queue_candidate_v0_2.csv`
