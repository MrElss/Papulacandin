# Open-source docking pipeline (GNINA + Vina + DiffDock) — FKS1 · papulacandin B

The AF3/local docking arms of **Stage A · A2**, entirely open-source (no Schrödinger). Produces poses + scores to combine with the Boltz-2 ensemble into a consensus for **A3 MD**.

## Read this first — two modes

papulacandin B is large and floppy: **64 heavy atoms, 30 rotatable bonds, MW 901**, ~30 Å extended, and the TM5–TM8 pocket box is ~**29 × 40 × 40 Å**. De-novo docking of such a ligand is unreliable. So:

- **Mode A (PRIMARY): GNINA `--minimize` + CNN-rescore of the co-folding poses (Boltz / AF3)** in each template. Reliable for a 900-Da flexible ligand; directly tests whether a predicted pose relaxes to a good CNN score and is consistent across templates.
- **Mode B (secondary): Vina + GNINA de-novo box docking + DiffDock blind docking** — independent searches, treated as coarse cross-checks.

## Files

Scripts: `00_setup_env.sh`, `01_prep_receptor.sh`, `02_make_box.py`, `03_prep_ligand.py`, `04_dock_vina.sh`, `05_dock_gnina.sh`, `06_dock_diffdock.sh`, `07_collect_scores.py`.
Pre-generated & **validated on 9WZU** (ready to use): `box_9WZU_apo.conf`, `papB.sdf`, `papB_ensemble.sdf`, `05b_align_pose.pml`.

## Inputs you provide

- `templates/` — `9WZU_apo.pdb` (you have it), plus `caspo_T2.pdb` and `apo_8WL6.pdb` fetched from RCSB. Run `02_make_box.py` per template. **The box file must be named `box_<template_base>.conf`** (same base as the receptor, e.g. receptor `caspo_T2` → `box_caspo_T2.conf`); the docking scripts now fail loudly if it is missing rather than falling back to another template's box.
- `cofold_poses/` — ligand-only SDFs from your Boltz + AF3 runs, each placed in the target template's frame (see *pose transfer*), named `<template>__<source>.sdf` (e.g. `9WZU_apo__boltz.sdf`).
- GPU recommended (GNINA CNN scoring, DiffDock).

## Run order

```
bash 00_setup_env.sh                                   # once
#  put templates in templates/
bash 01_prep_receptor.sh                               # -> receptors/*_H.pdb (+ .pdbqt)
python 02_make_box.py --pdb templates/caspo_T2.pdb --out box_caspo_T2   # per new template (9WZU done)
python 03_prep_ligand.py                               # papB.sdf (already made; re-run to regen)
#  put co-folding poses in cofold_poses/
bash 05_dock_gnina.sh        # Mode A (primary) + Mode B de-novo
bash 04_dock_vina.sh         # Mode B
bash 06_dock_diffdock.sh     # Mode B (from the DiffDock env)
python 07_collect_scores.py  # -> docking_scores.csv
```

## Pose transfer (Boltz/AF3 complex → template frame) — `05b_align_pose.pml`

Boltz/AF3 give a predicted *complex*. The included `05b_align_pose.pml` superposes the predicted receptor onto an experimental template and writes the ligand out in the template frame (it auto-falls back from `align` to `super` for lower-identity cases):

```
pymol -cq 05b_align_pose.pml -- templates/9WZU_apo.pdb boltz_papB_complex.pdb cofold_poses/9WZU_apo__boltz.sdf
```

Repeat per template × source, naming outputs `<template>__<source>.sdf` (e.g. `9WZU_apo__af3.sdf`) so `05_dock_gnina.sh` and `07_collect_scores.py` keep the Boltz/AF3 origin in the `source` column.

## Reading the scores

- **GNINA** — `CNNscore` (0–1, pose quality), `CNNaffinity` (pKd-like), `minimizedAffinity` (kcal/mol). A co-folding pose that **stays put under `--minimize`** and scores `CNNscore ≳ 0.7` with good `CNNaffinity` **across templates** is strong evidence.
- **Vina** — affinity kcal/mol (coarse here).
- **DiffDock** — rank-1 confidence; check whether rank-1 lands in the TM5–TM8 pocket.

## Consensus → A3

Superpose all top poses (Boltz + AF3 + GNINA-min + Vina + DiffDock) into one frame, cluster by ligand heavy-atom RMSD (< 2.5 Å), keep clusters populated by **≥2 methods on ≥2 templates** → **≤3 candidate models** → **A3 membrane MD** (metastability filter). This upgrades the B1 self-consistency evidence from Boltz-only to cross-method.

## Caveats (keep the reasoning honest)

- Docking treats the receptor as (near-)rigid and **does not model the lipid bilayer**; membrane effects are handled in A3 MD. Sanity-check poses at the lipid-facing edge against the membrane plane.
- 30 rotatable bonds ≈ Vina's 32-torsion ceiling → Vina is a coarse cross-check only.
- The pocket box encodes the **TM5–TM8 hypothesis**; whether that site is correct is settled by mutagenesis / engagement / structure (plan axes B3–B5), not by docking.
