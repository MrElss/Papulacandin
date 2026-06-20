#!/usr/bin/env python3
"""
export_assets.py
================
Bundle every project figure and table into a single, cleanly-organized ZIP so they
can be inserted into slides / documents independently of the PPTX (which may
render embedded objects inconsistently across viewers).

Produces outputs/Papulacandin_deck_assets.zip containing:
  figures/   — all .png figures, phase-prefixed
  tables/    — every .csv  AND a rendered .png image of each table
  reports/   — the markdown / txt write-ups
  INDEX.md   — maps every asset to the phase / slide it supports
"""
import os
import glob
import shutil
import zipfile
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
STAGE = os.path.join(OUT, "deck_assets")
ZIP = os.path.join(OUT, "Papulacandin_deck_assets.zip")


def _fmt(v):
    if isinstance(v, float):
        return f"{v:.3g}"
    s = str(v)
    return s if len(s) <= 42 else s[:39] + "…"


def render_table_png(csv_path, png_path, max_rows=26, max_cols=12):
    df = pd.read_csv(csv_path)
    note = ""
    if df.shape[1] > max_cols:
        df = df.iloc[:, :max_cols]; note += f" (first {max_cols} cols)"
    if df.shape[0] > max_rows:
        df = df.head(max_rows); note += f" (first {max_rows} rows)"
    disp = df.map(_fmt)
    ncol, nrow = disp.shape[1], disp.shape[0]
    fig_w = min(2.0 + 1.55 * ncol, 26)
    fig_h = min(0.9 + 0.32 * (nrow + 1), 22)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=disp.values,
                   colLabels=[textwrap.fill(c, 16) for c in disp.columns],
                   cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.35)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if r == 0:
            cell.set_facecolor("#1F3A5F"); cell.set_text_props(color="white", weight="bold")
        elif r % 2 == 0:
            cell.set_facecolor("#f1f4f8")
    ax.set_title(os.path.basename(csv_path) + note, fontsize=10, weight="bold", pad=10)
    fig.tight_layout()
    fig.savefig(png_path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    if os.path.exists(STAGE):
        shutil.rmtree(STAGE)
    for sub in ("figures", "tables", "reports"):
        os.makedirs(os.path.join(STAGE, sub))

    figs = sorted(glob.glob(os.path.join(OUT, "*.png")))
    for f in figs:
        shutil.copy2(f, os.path.join(STAGE, "figures", os.path.basename(f)))

    csvs = sorted(glob.glob(os.path.join(OUT, "*.csv")))
    rendered = 0
    for c in csvs:
        base = os.path.splitext(os.path.basename(c))[0]
        shutil.copy2(c, os.path.join(STAGE, "tables", os.path.basename(c)))
        try:
            render_table_png(c, os.path.join(STAGE, "tables", base + ".png"))
            rendered += 1
        except Exception as e:
            print(f"[warn] could not render {base}: {e}")

    for r in sorted(glob.glob(os.path.join(OUT, "*.md")) + glob.glob(os.path.join(OUT, "*.txt"))):
        shutil.copy2(r, os.path.join(STAGE, "reports", os.path.basename(r)))

    write_index(figs, csvs)

    with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(STAGE):
            for fn in files:
                fp = os.path.join(root, fn)
                z.write(fp, os.path.relpath(fp, STAGE))
    size_mb = os.path.getsize(ZIP) / 1e6
    print(f"figures: {len(figs)} | tables: {len(csvs)} (rendered {rendered} png) | "
          f"reports: {len(glob.glob(os.path.join(STAGE,'reports','*')))}")
    print(f"Wrote {ZIP}  ({size_mb:.1f} MB)")


INDEX = """# Papulacandin serum-gap — figure & table asset pack

Everything that backs the slide deck, as standalone files. `tables/` holds each
dataset as both the raw `.csv` and a rendered `.png` (drop-in image if a viewer
mangles native tables).

## figures/  (by phase)
- phase1_clogp_vs_serum_mic, phase1_descriptor_boxplots — 2D SAR
- phase2_attrition, phase2_ladder_anchors — where compounds fail / activity ladder
- fig_phase3a_chemspace, fig_phase3a_cv_bars, fig_phase3a_score_dist — external FKS model
- fig_phase3b_coeffs, fig_phase3b_optimism, phase3b_pred_vs_obs — within-series serum model
- fig_phase4_existing_leads, fig_phase4_score_sanity — design leads / score sanity
- phase5_score_distributions, phase5_retrospective — generative library / 2D baseline ρ=0.32
- phase7_retrospective_qm — MMFF-proxy 3D lead (ρ=−0.45)
- phase8_retrospective_crest — real CREST overturns it
- phase9_electronic — electronic descriptors vs serum shift
- phase10_hsa — HSA docking vs serum shift (null)

## tables/  (csv + png)
- serum_gap_pairs — the 24-compound dependent variable
- phase1_descriptor_stats — descriptor split tolerant vs killed
- phase2_activity_ladder — enzyme → serum-free → serum ladder
- phase3a_cv_metrics, phase3a_papulacandin_scores, table_phase3_summary
- phase3b_loo_predictions, phase3b_series_model_coeffs
- table_phase4_existing_leads, table_phase4_candidate_ranking
- phase5_virtual_library — 30 generated analogs + scores
- phase6_qm_descriptors — 12 candidates, real CREST 3D descriptors
- phase7_known_qm_descriptors, phase7_retrospective_qm
- phase8_known_crest_descriptors, phase8_retrospective_crest, **phase8_confound_analysis**
- phase8b_finalist_gfn2_descriptors, phase8b_gfnff_vs_gfn2
- phase9_electronic_descriptors, **phase9_electronic_stats**

## reports/
- SYNTHESIS_phases1-10.md — the master narrative
- phase7–10_findings.md — per-phase write-ups
- design_rules.md, serum_gap_summary.txt, table_phase3_summary.md
"""


def write_index(figs, csvs):
    with open(os.path.join(STAGE, "INDEX.md"), "w") as fh:
        fh.write(INDEX)


if __name__ == "__main__":
    main()
