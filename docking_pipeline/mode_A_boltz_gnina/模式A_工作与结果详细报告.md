# 模式 A：Boltz 共折叠构象—GNINA 最小化与重打分工作及结果报告

## 1. 最终结论

模式 A 已完成：使用全长酿酒酵母 FKS1（UniProt P38631）与 Papulacandin B 生成 5 个 Boltz-2.1 共折叠复合物，将每个复合物中的配体分别转移到 3 个实验模板坐标系，并对所得 15 个构象运行 GNINA `--minimize --cnn_scoring rescore`。

本轮计算**没有得到满足预设可靠性标准的共折叠构象**：

- 15/15 个构象均完成最小化和 CNN 重打分，但 **0/15 的 CNNscore 达到 0.7**；最高值为 0.5564。
- 仅 2/15 个构象的 GNINA minimizedRMSD 小于 2.5 Å，而且都只出现在 9WZU 模板；没有一个 Boltz 样本在 3 个模板上都保持低 RMSD。
- 各 Boltz 样本的跨模板平均 CNNscore 仅为 0.1741–0.3548；没有出现“跨模板一致高分且最小化后保持原位”的构象。

因此，按本项目在 `README_DOCKING.md` 中预先设定的判据，本轮模式 A **不能作为可靠的正向计算证据**。这不否定 Papulacandin B 位点：该位点由 S643/R1357 突变、靶点结合及结构数据（B3–B5）实验确定；本轮共折叠和对接只检验候选姿势在这个已知口袋中的计算自洽性，**没有用于发现或确定结合位点**。

## 2. 任务范围与科学定位

Papulacandin B 具有 64 个重原子、较大分子量和大量可旋转键，常规从头对接不易稳定采样。模式 A 的设计是：

1. 在实验已知口袋中生成蛋白–配体共折叠构象；
2. 将同一共折叠构象映射到多个实验模板；
3. 用 GNINA 做局部能量最小化和 CNN 重打分；
4. 检查姿势是否在最小化后保持稳定，以及分数能否跨模板复现。

Boltz 输入使用 `force: false` 的软口袋约束。它把采样引导到实验位点，但不会把某一具体配体构象固定在口袋中。因此，本轮属于**实验位点内的定向姿势检验**，不是盲位点搜索。

## 3. 输入、参数与费用

| 项目 | 内容 |
|---|---|
| 蛋白 | Saccharomyces cerevisiae FKS1，UniProt P38631 |
| 蛋白长度 | 1876 aa，已校验 |
| 配体 | Papulacandin B，SMILES 来自本项目 `03_prep_ligand.py` |
| Boltz 模型 | Boltz-2.1 |
| Boltz 样本数 | 5 |
| Boltz 任务 ID | `sab_pred_emWQwlk77rweNMp4v9CV` |
| Boltz 运行名 | `papb-fks1-mode-a-20260713` |
| 估算费用 | USD 0.8000，用户确认后提交 |
| 口袋约束 | P38631 残基 635、638、639、643、695、1357、1360、1361、1364、1365、1368 的 0-based 索引；`max_distance_angstrom: 6.0`；`force: false` |
| 实验模板 | `9WZU_apo`、`apo_8WL6`、`caspo_T2` |
| 每个 Boltz 样本的模板数 | 3 |
| GNINA 构象总数 | 5 × 3 = 15 |

Papulacandin B 有 64 个重原子，超过 Boltz API ligand–protein binding 指标要求的少于 50 个原子的范围。因此本任务只请求结构预测，没有请求不适用的 Boltz 结合指标。

## 4. 完成的工作流程

### 4.1 建立独立目录

所有模式 A 文件均放在：

`docking_pipeline/mode_A_boltz_gnina/`

与父目录中已完成的模式 B 输出分开保存。

### 4.2 获取并校验 FKS1 序列

- 从 UniProt 获取 P38631 FASTA。
- 去除 FASTA 标题和换行后长度为 1876 aa。
- 保存为 `inputs/P38631.fasta`。
- 清理了首次下载时因 Windows 中文路径编码造成的误建目录；最终文件仅保留在正确的模式 A 目录。

### 4.3 安装和使用 Boltz API CLI

- 安装 Boltz API CLI v0.37.1：`C:/Users/dxn-ud/AppData/Local/Programs/Boltz/bin/boltz-api.exe`。
- 完成 OAuth 设备登录；认证信息保存在用户凭据存储中，不在项目目录内，也不应上传 GitHub。
- 先运行费用估算，用户确认 USD 0.8000 后提交任务。
- 后台下载完成状态：`succeeded / ready`。
- 下载归档包含 `metrics.json`、5 个 CIF 和 5 个 PAE NPZ。

### 4.4 提取和对齐 Boltz 构象

- 将归档解压到 `predictions/`。
- 修复了 `05b_align_pose.pml` 中两个非 ASCII 破折号，解决 Windows PyMOL 按 GBK 读取时的解码错误；计算逻辑没有改变。
- 使用 `C:/ProgramData/pymol/Scripts/pymol.exe`，将每个预测复合物的蛋白与各模板叠合。
- 从每次叠合后的复合物中提取非聚合物配体，写入 `cofold_poses/`。
- 验证得到 15 个 SDF，每个输入 SDF 均有 64 个配体原子。

### 4.5 集群软件与作业环境

集群检查结果：PATH、模块和常见公共软件目录中未发现可直接使用的 GNINA，但存在 Singularity、Slurm 和 A100 GPU 分区。

远程目录：

- 软件：`/home/share/zhangz/Apps/gnina-1.1/gnina`
- 作业：`/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713`
- 作业输入：`cofold_poses/`、`receptors/`、`provenance/`
- 作业输出：`out_gnina/`
- 作业日志：`logs/`

最初尝试 GNINA v1.3.2，但集群宿主 glibc 2.17 低于该二进制要求的 glibc 2.28/2.29，同时缺少 cuDNN 9 和多项 CUDA 12 运行库。最终改用 GNINA v1.1 官方二进制：

- 文件大小：306,470,832 字节；
- SHA-256：`cd5f12d88f63bef2f637ef183beb24f1620ba97944d88173b76fbd373e46d84f`；
- 要求的最高 GLIBC 为 2.14，可在集群 glibc 2.17 上运行；
- 复用 `/home/share/zhangz/Apps/envs/diffdock/lib` 中较新的 `libstdc++`/`libgcc_s`；
- `--version`、动态链接和 GPU 参数检查均在计算节点通过。

已删除远程不可用的 `/home/share/zhangz/Apps/gnina-1.3.2`，释放约 1.4 GB；也删除了本地用于传输的两个大型二进制，但保留 SHA-256 记录和全部安装/作业脚本。

### 4.6 GNINA 最小化和重打分

最终 Slurm 作业申请：

- 分区：`gpu`
- CPU：10 核
- GPU：1 × A100
- 内存：64 GB
- 命令核心：`gnina -r receptor -l pose --device 0 --minimize --cnn_scoring rescore -o output`

作业获得了 A100 配额，但 GNINA v1.1 日志报告 `No GPU detected`，因此 CNN 实际在 CPU 上执行。15 个构象仍在 2 分 52 秒内全部完成；算法步骤、输出字段和构象数均完整。该资源利用限制不应被误写为“GPU 加速成功”。

## 5. Slurm 作业历史与故障处理

| 作业 ID | 状态 | 作用/结果 |
|---|---|---|
| 524823 | FAILED | v1.3.2 安装；计算节点 `wget` 无法建立 GitHub SSL 连接 |
| 524824 | CANCELLED | 依赖 524823，未运行 |
| 524825 | FAILED | 改用 `curl`；旧版 curl 不支持 `--retry-all-errors` |
| 524826 | CANCELLED | 依赖 524825，未运行 |
| 524827 | CANCELLED | 修复 curl 后下载速度约 40 KB/s，预计超过作业时限；主动停止 |
| 524828 | CANCELLED | 依赖 524827，未运行 |
| 524829 | CANCELLED | 本机下载并校验 v1.3.2 后检查出 glibc/CUDA/cuDNN 不兼容；停止该路线 |
| 524830 | CANCELLED | 依赖 524829，未运行 |
| 524831 | COMPLETED | GNINA v1.1 动态链接、版本和参数验证通过 |
| 524832 | FAILED | GPU 作业在 GNINA 启动前因 `nvidia-smi` 不在 PATH 而退出 |
| 524833 | COMPLETED | 将 `nvidia-smi` 改为非致命诊断后，15/15 个模式 A 构象完成 |

完整 `sacct` 记录保存为 `logs/cluster/sacct_job_history.tsv`，各次 Slurm 日志保存在 `logs/cluster/`。

## 6. Boltz 结果

`metrics.json` 没有给每项显式样本编号；下表按归档中的 `all_sample_results` 顺序映射为 sample 0–4。sample 0 的指标与 `best_sample.metrics` 完全相同。

| 样本 | structure confidence | pTM | ligand ipTM | complex pLDDT | complex ipLDDT | complex PDE | complex iPDE |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.7778 | 0.8774 | 0.7585 | 0.7826 | 0.6964 | 1.1298 | 5.6836 |
| 1 | 0.7772 | 0.8815 | 0.7703 | 0.7789 | 0.6950 | 1.1040 | 5.5737 |
| 2 | 0.7767 | 0.8843 | 0.7744 | 0.7773 | 0.6809 | 1.0861 | 5.2226 |
| 3 | 0.7734 | 0.8768 | 0.7482 | 0.7797 | 0.6846 | 1.1303 | 5.8659 |
| 4 | 0.7734 | 0.8780 | 0.7481 | 0.7797 | 0.6888 | 1.1246 | 5.8611 |

这些指标显示 5 个样本的 Boltz 内部置信度相近，ligand ipTM 约为 0.748–0.774。但它们只是模型内部置信度，不能证明姿势能在不同实验模板中稳定，也不能替代结合实验、自由能计算或膜环境 MD。后续 GNINA 检验正是对这种跨模板稳定性的独立检查。

## 7. GNINA 15 个构象的详细结果

`minimizedRMSD` 为 GNINA 输出 SDF 中的字段。GNINA 将输入的 64 原子 SDF 规范化为 67 原子输出，因此没有用简单逐原子直接差值替代该程序报告的 RMSD。

| 模板 | Boltz 样本 | CNNscore | CNNaffinity | minimizedAffinity (kcal/mol) | minimizedRMSD (Å) |
|---|---:|---:|---:|---:|---:|
| 9WZU_apo | 0 | 0.5564 | 6.6534 | -7.2235 | 2.5045 |
| 9WZU_apo | 1 | 0.2019 | 6.2711 | -6.3935 | 1.8773 |
| 9WZU_apo | 2 | 0.2561 | 6.0074 | -5.2825 | 6.0243 |
| 9WZU_apo | 3 | 0.2450 | 6.4036 | -6.6890 | 3.1496 |
| 9WZU_apo | 4 | 0.2339 | 5.9466 | -2.5996 | 2.4186 |
| apo_8WL6 | 0 | 0.3629 | 6.1169 | -4.9997 | 3.4217 |
| apo_8WL6 | 1 | 0.2169 | 6.1442 | -3.6403 | 2.8570 |
| apo_8WL6 | 2 | 0.3756 | 6.0758 | -5.7186 | 5.6432 |
| apo_8WL6 | 3 | 0.1403 | 6.1806 | -4.8200 | 2.9179 |
| apo_8WL6 | 4 | 0.3799 | 6.0839 | -4.3584 | 3.2646 |
| caspo_T2 | 0 | 0.1451 | 5.4036 | -5.5923 | 3.0048 |
| caspo_T2 | 1 | 0.1311 | 5.5893 | -3.8556 | 3.0439 |
| caspo_T2 | 2 | 0.4038 | 5.5254 | -2.2721 | 4.0791 |
| caspo_T2 | 3 | 0.1370 | 5.7903 | -6.3578 | 2.8299 |
| caspo_T2 | 4 | 0.0391 | 5.6821 | -6.5419 | 4.1529 |

### 7.1 按模板汇总

| 模板 | 平均 CNNscore | 平均 CNNaffinity | 平均 minimizedAffinity | 平均 minimizedRMSD |
|---|---:|---:|---:|---:|
| 9WZU_apo | 0.2987 | 6.2564 | -5.6376 | 3.1949 Å |
| apo_8WL6 | 0.2951 | 6.1203 | -4.7074 | 3.6209 Å |
| caspo_T2 | 0.1712 | 5.5982 | -4.9239 | 3.4221 Å |

caspo_T2 的平均 CNNscore 明显低于两个 apo 模板；同一 Boltz 构象对受体模板状态较敏感，未表现出跨模板一致的高姿势质量。

### 7.2 按 Boltz 样本汇总

| Boltz 样本 | 跨模板平均 CNNscore | 平均 CNNaffinity | 平均 minimizedAffinity | 平均 minimizedRMSD |
|---:|---:|---:|---:|---:|
| 0 | 0.3548 | 6.0580 | -5.9385 | 2.9770 Å |
| 1 | 0.1833 | 6.0016 | -4.6298 | 2.5927 Å |
| 2 | 0.3452 | 5.8695 | -4.4244 | 5.2488 Å |
| 3 | 0.1741 | 6.1248 | -5.9556 | 2.9658 Å |
| 4 | 0.2176 | 5.9042 | -4.4999 | 3.2787 Å |

样本 0 的平均 CNNscore 最高，但从 9WZU 的 0.5564 降至 caspo_T2 的 0.1451，且三个模板的 RMSD 均未小于 2.5 Å。样本 2 的平均 CNNscore 次高，但平均 RMSD 达 5.2488 Å，说明在局部最小化中发生了较大位移。没有一个样本同时满足跨模板高 CNNscore 和低 RMSD。

### 7.3 全局统计

- CNNscore：平均 0.2550，中位数 0.2339，范围 0.0391–0.5564。
- CNNaffinity：平均 5.9916，范围 5.4036–6.6534。
- minimizedAffinity：平均 -5.0897 kcal/mol，范围 -7.2235 至 -2.2721 kcal/mol。
- minimizedRMSD：平均 3.4126 Å，中位数 3.0439 Å，范围 1.8773–6.0243 Å。
- CNNscore ≥ 0.7：0/15。
- minimizedRMSD < 2.5 Å：2/15，均不具备跨模板复现。

CNNaffinity 或 minimizedAffinity 单项较好不能补偿低 CNNscore 和较大构象位移。对本任务最重要的是姿势质量与跨模板一致性，而这两项没有达到预定要求。

## 8. 结果和文件完整性验证

- Boltz：5 个 CIF、5 个 PAE NPZ、1 个 `metrics.json`。
- 对齐输入：15 个 SDF，全部为 64 个原子。
- GNINA：15 个非空 `_min.sdf`、15 个对应 `_min.log`。
- 汇总 CSV：15 行；CNNscore、CNNaffinity、minimizedAffinity、minimizedRMSD 四个关键字段缺失数均为 0。
- 成功作业 `524833` 结束标记：`MODE_A_GNINA_V11_OK`。
- 工作软件 SHA-256 与上传前一致。

## 9. 本地目录和文件说明

完成时模式 A 目录约 203 MB，主要空间来自 Boltz 原始下载归档及 PAE 数组。

| 目录/文件 | 作用 |
|---|---|
| `README_MODE_A.md` | 模式 A 简要说明和运行逻辑 |
| `模式A_工作与结果详细报告.md` | 本报告 |
| `inputs/P38631.fasta` | 官方 FKS1 序列 |
| `config/boltz_payload.yaml` | Boltz 实体、5 样本和软口袋约束 |
| `boltz_runs/` | Boltz 下载器状态、归档及原始输出，约 134 MB |
| `predictions/` | 解包后的 5 CIF、5 PAE NPZ、`metrics.json`，约 69 MB |
| `cofold_poses/` | 15 个模板坐标系中的 ligand-only 输入 SDF |
| `out_gnina/` | 15 个最小化 SDF 和 15 个逐构象日志 |
| `logs/cluster/` | 所有 Slurm 尝试日志和 `sacct_job_history.tsv` |
| `scripts/align_boltz_predictions.ps1` | 5 × 3 PyMOL 批量对齐 |
| `scripts/run_mode_a_gnina.sh` | 独立的本地/WSL 模式 A GNINA 脚本 |
| `scripts/collect_mode_a_scores.py` | 评分、模板汇总、样本汇总和 Boltz 指标汇总 |
| `cluster/*.slurm` | v1.3.2 尝试、v1.1 验证和最终 10 核+A100 作业脚本 |
| `cluster/*.sha256` | 下载过的 GNINA 二进制校验值；大型二进制已删除 |
| `mode_A_scores.csv` | 15 个构象的完整 GNINA 结果 |
| `mode_A_scores_by_template.csv` | 3 个模板的统计汇总 |
| `mode_A_scores_by_sample.csv` | 5 个 Boltz 样本的跨模板统计 |
| `boltz_metrics_summary.csv` | 5 个 Boltz 样本的置信度指标 |
| `SHA256SUMS_CORE.txt` | 排除重复 Boltz 归档、PAE 大文件和缓存后的 78 个核心文件校验值 |

## 10. 建议上传 GitHub 的内容

### 10.1 建议普通 Git 上传

以下文件体积小且是复现、审核和解释结果所必需：

- `README_MODE_A.md`
- `模式A_工作与结果详细报告.md`
- `inputs/P38631.fasta`
- `config/boltz_payload.yaml`
- `scripts/`
- `cluster/*.slurm` 和 `cluster/*.sha256`
- `predictions/metrics.json`
- `predictions/sample_*_predicted_structure.cif`
- `cofold_poses/*.sdf`
- `out_gnina/*_min.sdf`
- `out_gnina/*_min.log`
- 4 个结果 CSV
- `SHA256SUMS_CORE.txt`
- `logs/cluster/sacct_job_history.tsv`
- 成功验证/运行日志 `verify_gnina11_524831.out`、`modeA_gnina11_524833.out`

失败安装日志也很小，可一并上传以保留完整故障处理记录；若只想发布科学结果，可不上传这些失败日志，但本报告中的作业历史应保留。

### 10.2 建议用 Git LFS 或 GitHub Release

`predictions/sample_*_pae.npz` 每个约 13 MB、合计约 65 MB。它们单个没有超过 GitHub 100 MB 硬限制，但会显著增大 Git 历史。建议二选一：

1. 使用 Git LFS 跟踪 `mode_A_boltz_gnina/predictions/*_pae.npz`；或
2. 将 5 个 PAE NPZ 放在 GitHub Release/外部数据归档中，在报告中给出链接。

### 10.3 不建议上传

- `boltz_runs/`：约 134 MB，主要与 `predictions/` 重复；其中 `.boltz-run.json` 是本地下载状态文件，不是科学结果。
- `boltz_runs/**/outputs/archive.tar.gz`：与已解包结果重复。
- Boltz OAuth 配置、Windows 凭据或任何密码/令牌。
- GNINA 二进制、Conda/Micromamba 环境、缓存和临时目录。
- 同一结果的压缩包与解包文件同时上传。

若排除 `boltz_runs/` 和 PAE NPZ，模式 A 的核心可审核结果只有约数 MB，适合普通 Git 上传。

## 11. 局限性与后续建议

1. **没有 AF3 构象。** 本轮只完成 Boltz 分支；如果后续获得 AF3 复合物，应使用同一对齐和 GNINA 流程独立处理，再比较跨方法一致性。
2. **已知口袋软约束。** 这适合检验实验位点内姿势，但不是独立的位点发现证据。
3. **没有显式膜和脂质。** FKS1 是膜蛋白，口袋边缘构象可能受膜环境影响；A3 膜环境 MD 仍是必要的稳定性过滤。
4. **GNINA 版本限制。** 模式 A 使用 v1.1；若模式 B 使用 v1.3.x，不宜把绝对 CNN 数值无条件混在同一阈值体系中。由于本轮最高 CNNscore 仅 0.5564 且跨模板明显下降，负向结论不依赖细小版本差异，但严格比较时仍应统一版本。
5. **A100 未被 GNINA v1.1 使用。** 作业分配了 A100，但程序未检测到 GPU，实际为 CPU CNN。结果完整，性能而非构象数量受到影响；如需严格 GPU/v1.3.x 复算，应使用带新 glibc、CUDA 12 和 cuDNN 9 的 Ubuntu 容器或本地 WSL 环境。
6. **不能从本轮挑一个“最好看”的姿势作为可靠模型。** 样本 0 在单个模板上相对最好，但没有跨模板高分；样本 2 分数尚可但位移过大。选择任一构象进入后续研究都应标记为探索性，而不是模式 A 已验证姿势。

## 12. 可复现命令入口

```powershell
# 重新对齐一个或多个 Boltz CIF
./scripts/align_boltz_predictions.ps1 -PredComplex ./predictions/sample_0_predicted_structure.cif

# 重新生成全部 CSV
python ./scripts/collect_mode_a_scores.py
```

集群脚本：

- `cluster/02_verify_gnina_v11.slurm`
- `cluster/03_mode_a_gnina_v11_gpu.slurm`

成功运行所对应的远程输出仍保留在：

`/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713/out_gnina/`

---

报告生成日期：2026-07-13（Asia/Shanghai）。
