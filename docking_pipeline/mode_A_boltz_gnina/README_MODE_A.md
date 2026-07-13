# 模式 A：Boltz 共折叠构象的 GNINA 最小化与重打分

本目录独立保存模式 A 的输入、Boltz 任务记录、预测结构、模板坐标系下的配体、GNINA 结果与日志，不会与父目录中已完成的模式 B 输出混放。

## 完成状态（2026-07-13）

模式 A 已完成：5 个 Boltz 样本 × 3 个模板，共 15 个 GNINA 最小化/重打分结果，4 个汇总 CSV 均已生成。0/15 的 CNNscore 达到 0.7，且没有构象同时表现出跨模板高 CNNscore 和低 minimizedRMSD，因此本轮没有获得满足预设可靠性标准的共折叠姿势。详细过程、作业记录、完整评分、限制与 GitHub 上传清单见 `模式A_工作与结果详细报告.md`。

## 科学目的

FKS1 的 Papulacandin B 位点由实验确定（S643/R1357 突变、靶点结合及结构数据 B3–B5），不是通过对接发现。这里使用实验位点的软口袋约束生成 Boltz 共折叠构象，再检验这些构象在 9WZU、8WL6 和 caspofungin-bound T2 三个模板中经 GNINA 局部最小化后能否保持稳定，并获得跨模板一致的高 CNN 重打分。

## 目录

- `inputs/P38631.fasta`：UniProt P38631 官方 FKS1 序列（1876 aa）。
- `config/boltz_payload.yaml`：Boltz-2.1 输入；5 个样本，Papulacandin B 为链 B，实验位点为 `force: false` 的软约束。
- `boltz_runs/`：Boltz API 的任务元数据、状态和下载记录。
- `predictions/`：下载并解包的 Boltz 复合物结构。
- `cofold_poses/`：经 `05b_align_pose.pml` 对齐到各实验模板坐标系后的 ligand-only SDF。
- `out_gnina/`：每个构象的 GNINA `--minimize --cnn_scoring rescore` 输出和日志。
- `mode_A_scores.csv`：CNNscore、CNNaffinity、minimizedAffinity 及 GNINA minimizedRMSD 汇总；另有按模板、按样本和 Boltz 指标汇总 CSV。

## 实际执行流程

1. 估算 Boltz 费用并由用户确认后提交后台任务。
2. 下载所有预测样本到 `predictions/`。
3. 对每个样本运行 `scripts/align_boltz_predictions.ps1`，生成 3 个模板坐标系下的构象。
4. 将 15 个构象和 3 个受体上传到 `/home/share/zhangz/Jobs/mode_A_boltz_gnina_20260713`，使用 GNINA v1.1 完成 10 核+A100 Slurm 作业；程序未检测到 GPU，因此实际为 CPU CNN。`scripts/run_mode_a_gnina.sh` 保留为本地/WSL 的模式 A 专用替代入口。
5. 运行 `scripts/collect_mode_a_scores.py` 汇总结果，并撰写最终工作与结果分析文档。

## 解释边界

Papulacandin B 有 64 个重原子，超过 Boltz API ligand–protein binding 指标所支持的少于 50 个原子限制，因此本任务只申请结构预测，不请求该结合指标。高 CNNscore 与低构象位移可支持“该共折叠姿势在已知口袋内自洽且可稳定松弛”，但不能替代结合自由能计算或膜环境 MD。
