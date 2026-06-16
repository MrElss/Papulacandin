# 外部 FKS / 葡聚糖合酶抑制剂数据检索策略

## 1. 总原则

先找“最相关、最可信”的小集合，不要一开始下载几万条广义抗真菌数据。

推荐顺序：

1. Tier 1: 明确 FKS1/FKS2 或 beta-1,3-glucan synthase 抑制剂。
2. Tier 2: 真菌细胞壁相关，但不是直接 FKS。
3. Tier 3: 广义抗真菌活性小分子，只用于大规模预训练。

## 2. Tier 1 种子化合物/家族

优先从这些名字开始检索：

```text
caspofungin
micafungin
anidulafungin
rezafungin
cilofungin
echinocandin B
pneumocandin B0
aculeacin A
mulundocandin
enfumafungin
ibrexafungerp
SCY-078
beta-1,3-glucan synthase inhibitor
1,3-beta-glucan synthase inhibitor
FKS1 inhibitor
FKS2 inhibitor
glucan synthase inhibitor
```

## 3. 数据源优先级

### 3.1 ChEMBL

用途：

- 查 compound ID、canonical SMILES、InChIKey。
- 查 target、assay、activity。
- 找文献来源和机制标注。

推荐检索方式：

1. 先按 compound name 找 molecule。
2. 再按 molecule ChEMBL ID 拉 activity。
3. 同时按 target keyword 搜索 `glucan synthase`, `FKS1`, `FKS2`。
4. 导出 activity 时保留 assay description、target name、standard type、standard value、standard units、document。

注意：

- ChEMBL 的 target assignment 可能是 inferred 或 broad target，不能自动当作直接 FKS 抑制。
- 对大环肽/脂肽类，结构可能复杂，要优先保留 source SMILES/SDF 和 InChIKey。

### 3.2 PubChem

用途：

- 查结构、CID、同义名、IUPAC、供应商/专利来源。
- 补结构而不是作为唯一活性来源。

推荐检索方式：

1. 用 compound name 搜 CID。
2. 导出 Canonical SMILES、Isomeric SMILES、InChIKey、SDF。
3. 查 synonyms，补充 `synonyms` 字段。
4. 对有 BioAssay 的化合物，记录 AID 和 assay 描述，但先不要直接混入训练。

### 3.3 BindingDB

用途：

- 查蛋白-小分子结合数据。
- 如果有 glucan synthase/FKS 相关 target，作为高价值补充。

注意：

- BindingDB 更偏 binding affinity，不一定覆盖复杂天然产物和 whole-cell MIC。
- 若没有命中，也要在 `source_notes/` 记录“已检索但无可靠数据”。

### 3.4 UniProt

用途：

- 统一 target metadata：FKS1/FKS2、物种、UniProt ID、同义名。
- 不直接提供小分子活性，但对 `external_targets.csv` 很有用。

建议先建目标表：

```text
Saccharomyces cerevisiae FKS1
Candida albicans FKS1 / GSC1
Candida albicans FKS2
Candida auris FKS1
Aspergillus fumigatus FKS1-like glucan synthase
```

### 3.5 PubMed / 文献

用途：

- 查数据库没有结构化收录的历史天然产物、半合成物、专利前药物。
- 复核机制、体内/血清/PK 线索。

推荐关键词组合：

```text
"beta-1,3-glucan synthase inhibitor" antifungal
"1,3-beta-glucan synthase" echinocandin MIC
"FKS1" inhibitor antifungal compound
"enfumafungin" glucan synthase IC50
"ibrexafungerp" glucan synthase MIC
"pneumocandin" glucan synthase
"aculeacin" glucan synthase
"mulundocandin" antifungal glucan synthase
```

## 4. 第一轮人工检索目标

第一轮只做 30-80 个高相关化合物，不追求大而全。

建议优先列表：

1. 已上市/临床 FKS 靶向药：caspofungin、micafungin、anidulafungin、rezafungin、ibrexafungerp。
2. 历史 FKS 天然产物/半合成物：echinocandin B、pneumocandin B0、cilofungin、aculeacin A、mulundocandin。
3. enfumafungin 系列：enfumafungin、SCY-078/ibrexafungerp 相关公开类似物。
4. 文献中明确写 beta-1,3-glucan synthase inhibition 的其他结构。

## 5. 不要做的事

1. 不要把外部 MIC 直接追加到 `data/curated/activity_table.csv`。
2. 不要把外部化合物重新编号为 `PAPU-xxxx`。
3. 不要逐个手画所有外部化合物。
4. 不要把 broad antifungal rows 当成 FKS inhibitor。
5. 不要把 ChEMBL/PubChem 的所有 assay 盲目导入训练。

