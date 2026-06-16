# 外部 FKS / 葡聚糖合酶抑制剂自动检索 v0.1

生成时间：2026-06-11 12:21

## 1. 本轮做了什么

本脚本根据 `config/external_fks_seed_queries.csv` 中的种子关键词，自动访问 ChEMBL、PubChem 和 BindingDB，拉取候选结构、同义名、活性、assay 描述和来源编号。

重要边界：

- 本轮输出是 `candidate_v0_1`，还不是人工复核后的正式外部库。
- 没有覆盖 `external_compounds.csv`、`external_activity_table.csv` 等人工模板。
- 外部数据默认不进入 Papulacandin exact MIC 回归。
- PubChem BioAssay 和 BindingDB 相似结构命中需要人工复核后才能用于建模。

## 2. 输出文件

| 文件 | 说明 |
| --- | --- |
| `data/external/fks_inhibitors/source_exports/external_compounds_candidate_v0_1.csv` | 候选外部化合物结构与来源 ID |
| `data/external/fks_inhibitors/source_exports/external_activity_table_candidate_v0_1.csv` | 候选活性、assay、target 和来源编号 |
| `data/external/fks_inhibitors/source_exports/external_sources_candidate_v0_1.csv` | 每一次数据库请求和原始 JSON 文件位置 |
| `data/external/fks_inhibitors/source_exports/external_inclusion_decisions_candidate_v0_1.csv` | 候选纳入层级和建模用途初判 |
| `data/external/fks_inhibitors/raw_downloads/auto_fetch_external_fks_v0_1` | 原始 JSON 缓存目录 |

## 3. 数量概览

| 项目 | 数量 |
| --- | ---: |
| 种子关键词 | 17 |
| 候选化合物 | 65 |
| 候选活性行 | 1011 |
| 直接 FKS/glucan-synthase 候选活性行 | 348 |
| Tier 1 候选化合物 | 65 |
| source/API 请求记录 | 144 |
| RDKit 描述符可用 | yes |

## 4. 活性来源分布

| source_database | n_activity_rows |
| --- | --- |
| BindingDB | 39 |
| ChEMBL | 644 |
| PubChem BioAssay | 328 |

## 5. Tier 1 候选化合物预览

| external_compound_id | preferred_name | compound_family | activity_relevance_tier |
| --- | --- | --- | --- |
| EXT-FKS-0001 | 14C-ANIDULAFUNGIN | echinocandin | no_activity_imported |
| EXT-FKS-0002 | ACULEACIN A | aculeacin | broad_antifungal;direct_fks |
| EXT-FKS-0003 | ANIDULAFUNGIN | echinocandin | broad_antifungal;context_only;direct_fks |
| EXT-FKS-0004 | BAY-1093884 | aculeacin | no_activity_imported |
| EXT-FKS-0005 | CASPOFUNGIN | echinocandin | broad_antifungal;direct_fks |
| EXT-FKS-0006 | CASPOFUNGIN ACETATE | echinocandin | no_activity_imported |
| EXT-FKS-0007 | CHEMBL217727 | echinocandin | broad_antifungal |
| EXT-FKS-0008 | CHEMBL217728 | echinocandin | broad_antifungal |
| EXT-FKS-0009 | CHEMBL2369955 | mulundocandin | broad_antifungal |
| EXT-FKS-0010 | CHEMBL260117 | echinocandin | broad_antifungal |
| EXT-FKS-0011 | CHEMBL265537 | mixed_fks | direct_fks |
| EXT-FKS-0012 | CHEMBL265871 | mixed_fks | direct_fks |
| EXT-FKS-0013 | CHEMBL266048 | mixed_fks | direct_fks |
| EXT-FKS-0014 | CHEMBL269275 | mixed_fks | direct_fks |
| EXT-FKS-0015 | CHEMBL3037784 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0016 | CHEMBL3679922 | enfumafungin | context_only |
| EXT-FKS-0017 | CHEMBL3760049 | echinocandin | broad_antifungal |
| EXT-FKS-0018 | CHEMBL385234 | mixed_fks | direct_fks |
| EXT-FKS-0019 | CHEMBL385455 | mixed_fks | direct_fks |
| EXT-FKS-0020 | CHEMBL387213 | mixed_fks | direct_fks |
| EXT-FKS-0021 | CHEMBL405670 | mixed_fks | direct_fks |
| EXT-FKS-0022 | CHEMBL408597 | echinocandin | broad_antifungal |
| EXT-FKS-0023 | CHEMBL410953 | mixed_fks | direct_fks |
| EXT-FKS-0024 | CHEMBL412195 | pneumocandin | broad_antifungal;direct_fks |
| EXT-FKS-0025 | CHEMBL412790 | mixed_fks | direct_fks |
| EXT-FKS-0026 | CHEMBL412839 | pneumocandin | broad_antifungal |
| EXT-FKS-0027 | CHEMBL413521 | mixed_fks | direct_fks |
| EXT-FKS-0028 | CHEMBL414288 | mixed_fks | direct_fks |
| EXT-FKS-0029 | CHEMBL414343 | echinocandin | broad_antifungal |
| EXT-FKS-0030 | CHEMBL414362 | mixed_fks | direct_fks |

仅显示前 30 行；完整表见 CSV。

## 6. 需要人工复核的重点

1. 先看 `activity_relevance_tier = direct_fks` 的行，确认 assay 描述是否真的是 FKS / beta-1,3-glucan synthase。
2. 对 BindingDB 行，要检查相似结构是否等同于目标化合物；相似命中不能自动当成同一化合物。
3. 对 PubChem BioAssay 行，要打开 AID 或 PMID，确认 endpoint、单位、物种和 assay 语境。
4. ChEMBL 的 `assay_type = B/F` 只是数据库分类，不等同于你课题里的生化/全细胞层级，关键行仍需回源文献。
5. 复核后再把可靠行复制或合并进正式 `external_*.csv`。

## 7. 警告与失败请求

无。