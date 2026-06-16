# 外部 FKS / 葡聚糖合酶抑制剂自动检索 v0_2

生成时间：2026-06-11 13:56

## 1. 本轮做了什么

本脚本根据 `config/external_fks_seed_queries_expanded_v0_2.csv` 中的种子关键词，从已缓存的 ChEMBL JSON 重建候选结构、活性、assay 描述和来源编号；PubChem/BindingDB 未纳入本次 v0.2 主扩容表。

重要边界：

- 本轮输出是 `candidate_v0_2`，还不是人工复核后的正式外部库。
- 没有覆盖 `external_compounds.csv`、`external_activity_table.csv` 等人工模板。
- 外部数据默认不进入 Papulacandin exact MIC 回归。
- 本轮 v0.2 只使用 ChEMBL 缓存；PubChem/BindingDB 可在后续作为补结构和补证据来源。
- 本轮设置 `--max-output-compounds 500`，若原始命中超过上限，则按 direct FKS 证据、活性行数和结构完整性优先保留。

## 2. 输出文件

| 文件 | 说明 |
| --- | --- |
| `data/external/fks_inhibitors/source_exports/external_compounds_candidate_v0_2.csv` | 候选外部化合物结构与来源 ID |
| `data/external/fks_inhibitors/source_exports/external_activity_table_candidate_v0_2.csv` | 候选活性、assay、target 和来源编号 |
| `data/external/fks_inhibitors/source_exports/external_sources_candidate_v0_2.csv` | 每一次数据库请求和原始 JSON 文件位置 |
| `data/external/fks_inhibitors/source_exports/external_inclusion_decisions_candidate_v0_2.csv` | 候选纳入层级和建模用途初判 |
| `data/external/fks_inhibitors/raw_downloads/auto_fetch_external_fks_v0_2` | 原始 JSON 缓存目录 |

## 3. 数量概览

| 项目 | 数量 |
| --- | ---: |
| 种子关键词 | 37 |
| 候选化合物 | 500 |
| 候选活性行 | 10613 |
| 直接 FKS/glucan-synthase 候选活性行 | 470 |
| Tier 1 候选化合物 | 322 |
| source/API 请求记录 | 1174 |
| RDKit 描述符可用 | yes |

## 4. 活性来源分布

| source_database | n_activity_rows |
| --- | --- |
| ChEMBL | 10613 |

## 5. Tier 1 候选化合物预览

| external_compound_id | preferred_name | compound_family | activity_relevance_tier |
| --- | --- | --- | --- |
| EXT-FKS-0001 | CASPOFUNGIN | echinocandin | broad_antifungal;cell_wall_related;direct_fks |
| EXT-FKS-0002 | ANIDULAFUNGIN | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0003 | MICAFUNGIN | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0016 | CHEMBL3760049 | echinocandin | broad_antifungal |
| EXT-FKS-0017 | ACULEACIN A | aculeacin | broad_antifungal;direct_fks |
| EXT-FKS-0018 | CILOFUNGIN | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0020 | CHEMBL4764856 | enfumafungin | broad_antifungal |
| EXT-FKS-0022 | CHEMBL3037784 | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0023 | ENFUMAFUNGIN | enfumafungin | broad_antifungal;direct_fks |
| EXT-FKS-0026 | URIDINE DIPHOSPHATE GLUCOSE | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0027 | CHEMBL2387203 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0028 | IBREXAFUNGERP | enfumafungin | broad_antifungal;direct_fks |
| EXT-FKS-0035 | CHEMBL437803 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0043 | FR-901469 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0044 | Echinocandin B | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0050 | CHEMBL267988 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0051 | CHEMBL385455 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0052 | CHEMBL405670 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0053 | CHEMBL407395 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0054 | CHEMBL412195 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0055 | CHEMBL414690 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0069 | CHEMBL1793852 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0070 | CHEMBL2371681 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0071 | CHEMBL2371765 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0072 | CHEMBL261821 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0073 | CHEMBL262998 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0074 | CHEMBL263952 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0075 | CHEMBL264601 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0076 | CHEMBL267407 | mixed_fks | broad_antifungal;direct_fks |
| EXT-FKS-0077 | CHEMBL268666 | pneumocandin | broad_antifungal;direct_fks |

仅显示前 30 行；完整表见 CSV。

## 6. 需要人工复核的重点

1. 先看 `activity_relevance_tier = direct_fks` 的行，确认 assay 描述是否真的是 FKS / beta-1,3-glucan synthase。
2. 对 `mixed_fks` 或广谱抗真菌行，要确认它们是直接 FKS 证据、FKS 突变菌株背景，还是仅可作为预训练背景数据。
3. ChEMBL 的 `assay_type = B/F` 只是数据库分类，不等同于你课题里的生化/全细胞层级，关键行仍需回源文献。
4. PubChem/BindingDB 未进入本轮 v0.2 主表；后续如补充这些来源，需要单独复核 AID、PMID、相似结构和 endpoint。
5. 复核后再把可靠行复制或合并进正式 `external_*.csv`。

## 7. 警告与失败请求

无。