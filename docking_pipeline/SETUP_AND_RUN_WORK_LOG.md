# Docking Pipeline Setup and Execution Work Log

Date: 2026-07-13  
Workspace: `C:\Users\dxn-ud\Desktop\Papulacandin_Cowork\28. 对接前准备\docking_pipeline`

## Outcome

The executable parts of `README_DOCKING.md` were completed through `07_collect_scores.py`. Mode B was run with Vina, GNINA de-novo docking, and DiffDock. The combined table contains 186 rows in `docking_scores.csv`.

Mode A was not run because `cofold_poses/` contains no genuine Boltz or AlphaFold 3 ligand poses. These are external scientific inputs and were not fabricated. Once those SDF files are supplied with the README naming convention, `05_dock_gnina.sh` can minimize and rescore them.

## Directories Created

- `templates/`: downloaded source mmCIF files and prepared protein templates.
- `cofold_poses/`: input location for future Boltz/AF3 ligand-only SDF files; currently empty.
- `receptors/`: hydrogenated PDBs, Vina PDBQTs, and DiffDock ESM-segmented PDBs.
- `out_vina/`: Vina poses and logs.
- `out_gnina/`: GNINA poses and logs.
- `out_diffdock/`: downloaded DiffDock poses from the cluster.
- `bin/`: local GNINA executable.
- `DiffDock/`: pinned DiffDock source checkout.
- `cluster/`: Slurm scripts, cluster configuration, source/model archives, wheels, and reproducibility helpers.
- `/home/share/zhangz/Apps`: cluster software, environments, models, caches, and compatibility libraries.
- `/home/share/zhangz/Jobs/docking_pipeline_20260713`: cluster inputs, logs, and DiffDock results.

## Template Preparation

Official structure files were downloaded and retained as:

- `templates/9WZU_source.cif`
- `templates/9PE1_source.cif`
- `templates/8WL6_source.cif`

Chain A protein coordinates were exported with PyMOL as:

- `templates/9WZU_apo.pdb`: 12,445 atoms, residues 141-1861.
- `templates/apo_8WL6.pdb`: 11,922 atoms, residues 146-1845.
- `templates/caspo_T2.pdb`: 13,391 heavy atoms, residues 132-1856; derived from PDB 9PE1.

All 11 configured pocket residues were present in every template. Template-specific boxes were generated and verified:

| Template | Center (x, y, z) | Size (x, y, z) A |
|---|---|---|
| 9WZU_apo | 198.224, 203.560, 210.755 | 29.1, 40.0, 39.9 |
| apo_8WL6 | 126.460, 157.365, 149.119 | 33.1, 40.0, 38.4 |
| caspo_T2 | 154.644, 140.252, 125.101 | 28.9, 40.0, 40.0 |

References: [RCSB 9WZU](https://www.rcsb.org/structure/9WZU), [EMDB/PDB 9PE1](https://www.ebi.ac.uk/emdb/EMD-71550), and the associated [Nature article](https://www.nature.com/articles/s41586-026-10409-7).

## Local Software Environments

Miniforge was installed under WSL at `/root/miniforge3`.

The `dock` environment is `/root/miniforge3/envs/dock` with Python 3.10.20 and the required packages, including RDKit, NumPy, SciPy, pandas, Open Babel, Meeko, Vina, Gemmi, and cuDNN. Vina reports `AutoDock Vina f458505-mod`.

GNINA v1.3.3 was installed as `bin/gnina` from the official release asset. The README's older download URL returned 404, so the current official release was used. SHA-256: `3340c1f49cd3c7c84d8699182a1c6af13c7fa2a22448d1204640446106f72172`. Reference: [GNINA releases](https://github.com/gnina/gnina/releases).

The local DiffDock checkout is pinned to commit `85c49b60d3e0b0182a59ee43a34a6d7036981284`. A local validation environment was created at `/root/miniforge3/envs/diffdock` with Python 3.9.18, PyTorch 1.13.1+cu117, PyG 2.2.0, e3nn 0.5.1, ESM 2.0.0, OpenFold 1.0.0, and RDKit 2022.03.3. Imports, CLI help, and a local RTX 4090 CUDA tensor test passed. Reference: [DiffDock repository](https://github.com/gcorso/DiffDock).

## Receptor and Ligand Preparation

PyMOL hydrogen addition was used because Open Babel pH protonation changed chain/residue identities. Meeko then generated the Vina receptor PDBQTs. All PDBQTs retain the 11 configured pocket residues.

- `receptors/9WZU_apo_H.pdb`: 24,899 atoms; `9WZU_apo.pdbqt` created. A false distant covalent contact at A1446 was excluded only from the Vina PDBQT; it is approximately 34 A from the docking box.
- `receptors/apo_8WL6_H.pdb`: 23,878 atoms; `apo_8WL6.pdbqt` created.
- `receptors/caspo_T2_H.pdb`: 26,776 atoms; `caspo_T2.pdbqt` created.
- `papB.pdbqt`: generated successfully with 32 torsional degrees of freedom reported by the PDBQT.

For DiffDock only, coordinates were not changed, but each long FKS1 chain was split into ESM segments of at most 900 residues because this DiffDock version truncates an individual ESM-2 sequence at 1,022 residues:

- `9WZU_apo_diffdock.pdb`: 900 + 627 residues.
- `apo_8WL6_diffdock.pdb`: 900 + 589 residues.
- `caspo_T2_diffdock.pdb`: 900 + 775 residues.

## Script Corrections and Compatibility Work

- Added reliable WSL Miniforge activation to `01_prep_receptor.sh`, `04_dock_vina.sh`, `05_dock_gnina.sh`, and `06_dock_diffdock.sh`.
- Made receptor preparation resume-safe so validated outputs are not overwritten.
- Normalized CRLF-derived box values before GNINA execution and added resume checks.
- Set all DiffDock Slurm submissions to `--cpus-per-task=10` as requested.
- Added a cluster-specific YAML because this DiffDock version loads YAML after CLI parsing and otherwise overwrites model and sampling paths.
- Changed DiffDock batch size from 6 to 5 so 40 samples divide evenly; the old sampler creates a wrong-sized noise tensor for a remainder batch.
- Replaced RDKit's bulk `Conformer.GetPositions()` bridge with coordinate-equivalent per-atom extraction after a reproducible native segmentation fault. No molecular coordinates or model parameters were changed by this workaround.
- Added explicit GPU availability validation before inference to prevent silent CPU fallback.

The original `00_setup_env.sh` was not used verbatim because several pinned URLs/package combinations are outdated. Its intended environment was implemented and validated with compatible current installers while preserving the README workflow.

## Cluster Public-Software Check and Environment

The module list exposed only CUDA 10.1, but a targeted public-directory check found CUDA 10.2, 11.2, 11.8, 12.2, 12.4, and 12.6 under `/opt/app/cuda`; `/opt/app/cuda/current` points to 11.8. Public CUDA 11.8 supplied `nvcc`, headers, and development libraries, so those were reused instead of downloading another CUDA toolkit.

No public DiffDock, PyTorch, Conda, Miniforge, or compatible ESM cache was found. Project-specific software was therefore placed under `/home/share/zhangz/Apps`:

- `/home/share/zhangz/Apps/envs/diffdock`: cluster-native Python environment.
- `/home/share/zhangz/Apps/DiffDock`: source tree.
- `/home/share/zhangz/Apps/DiffDock_models/v1.1`: score and confidence models.
- `/home/share/zhangz/Apps/torch_cache`: ESM-2 model cache.
- `/home/share/zhangz/Apps/cuda-compat-12-4`: NVIDIA 550-series user-space driver compatibility libraries.

The GPU node exposed A100 device files but not system `libcuda.so`. NVIDIA's official `cuda-compat-12-4` RPM (550.54.15) was extracted without root into Apps to match the node's 550-series kernel driver. A Slurm GPU probe then passed: CUDA available, device `NVIDIA A100-PCIE-40GB`, CUDA tensor sum 28.

The ESM-2 650M model and contact-regression files were downloaded from the official Fair ESM file server, uploaded, and re-hashed on the cluster:

- `esm2_t33_650M_UR50D.pt`: 2,607,531,725 bytes; SHA-256 `c89758ebab0cecd17aebf8fae7a3000f56cac3e544606e1f8780696c00b8cf16`.
- `esm2_t33_650M_UR50D-contact-regression.pt`: 3,687 bytes; SHA-256 `8ffe6edbd4173dc8d45c2cd5cb27d43aad77ec26b4c768200c58ae1f96693575`.

OpenFold was compiled on the cluster against public CUDA 11.8 for A100 `sm_80`. Import validation reported PyTorch 1.13.1+cu117, CUDA ABI 11.7, and PyG 2.2.0.

Final relevant jobs:

- `524669`: cluster environment/OpenFold build; completed.
- `524675`: ESM cache load validation; completed.
- `524684`: A100/PyTorch CUDA probe; completed.
- `524686`: final DiffDock inference, 10 CPUs + 1 A100; completed in 18:36 with exit code 0.

No password or authentication secret was written to the workspace or this log.

## Docking Results

### Vina Mode B

Each template produced 20 modes with exhaustiveness 64 and 16 local CPU cores:

| Template | Best Vina affinity (kcal/mol) |
|---|---:|
| 9WZU_apo | -7.266 |
| apo_8WL6 | -8.153 |
| caspo_T2 | -9.597 |

The large search-volume warning is expected for the README's approximately 29 x 40 x 40 A hypothesis box. These scores are coarse cross-checks, not biological proof.

### GNINA Mode B

Each template produced 20 de-novo poses:

| Template | Top minimized affinity | Top CNNscore | Top CNNaffinity |
|---|---:|---:|---:|
| 9WZU_apo | -7.94953 | 0.50880 | 6.84092 |
| apo_8WL6 | -7.71658 | 0.24219 | 5.70208 |
| caspo_T2 | -11.28892 | 0.60730 | 7.18417 |

### DiffDock Mode B

Job `524686` ran 40 samples per template, 20 configured inference steps (19 actual steps), batch size 5, and generated 41 files per template: one unannotated `rank1.sdf` plus 40 confidence-labelled poses.

| Template | Rank-1 confidence |
|---|---:|
| 9WZU_apo | -2.44 |
| apo_8WL6 | -2.80 |
| caspo_T2 | -2.08 |

Occasional single-sample NaN score warnings were handled by DiffDock's existing small-perturbation fallback; all three complexes completed and wrote all expected files.

## Validation and Deliverables

- `out_vina/`: 3 pose files and 3 logs.
- `out_gnina/`: 3 pose files and 3 logs.
- `out_diffdock/`: 123 SDF files, exactly 41 per template, downloaded from the completed Slurm job.
- `docking_scores.csv`: 186 rows: 3 Vina, 60 GNINA, and 123 DiffDock records.
- All edited shell scripts pass Bash syntax checks; Python scripts compile.

The `rank1.sdf` file has no confidence in its filename, while `rank1_confidence*.sdf` is the scored duplicate. Therefore the collector intentionally contains 123 DiffDock rows rather than 120; use the confidence-labelled rank-1 file for ranking.

## Remaining Scientific Work

1. Supply genuine Boltz/AF3 ligand-only SDF files in `cofold_poses/` to run the README's primary Mode A minimization/rescoring.
2. Inspect whether DiffDock rank-1 poses land in the TM5-TM8 pocket.
3. Superpose and cluster top poses across methods/templates as described in the README, then select no more than three consensus candidates for A3 membrane MD.

Docking outputs should be treated as pose-generation and cross-method consistency evidence, not proof of the binding site or affinity.
