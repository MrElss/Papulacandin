#!/usr/bin/env python3
"""
phase1_serum_shift_sar.py
=========================
PHASE 1 — Why does serum kill activity? (hypothesis generation)

We have 24 compounds with matched serum-free vs serum MIC (built in
serum_gap_analysis.py). Some keep partial serum activity, most are wiped out.
This script asks: which molecular properties separate the "serum-tolerant"
compounds from the "serum-killed" ones?

The leading mechanistic hypothesis in the literature is that lipophilic
glycolipids bind serum proteins (albumin / lipoproteins), which sequesters the
drug and lowers the free concentration available to reach the fungal cell. If
true, we expect lipophilicity-related descriptors (clogP, lipophilic side-chain
size, low polar surface area) to track with the serum shift.

We test that with simple, transparent statistics appropriate for SMALL data:
  * group medians (tolerant vs killed) for each descriptor
  * Spearman rank correlation of each descriptor with the serum MIC
  * Mann-Whitney U test between the two groups
N is tiny (~24, and only one major chemical series), so we treat every p-value
as HYPOTHESIS-GENERATING, not proof. The point is to find directions, not to
claim significance.

Outputs (analysis/outputs/):
  * phase1_descriptor_stats.csv  — per-descriptor group medians + correlations
  * phase1_clogp_vs_serum_mic.png
  * phase1_descriptor_boxplots.png
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------------------
# Load: serum-gap pairs + molecular descriptors
# --------------------------------------------------------------------------
pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))

DESCRIPTORS = ["clogp", "mw_exact", "tpsa", "hbd", "hba", "rotb", "fsp3",
               "long_acyl_chain_length"]
cm_small = cm[["compound_id"] + DESCRIPTORS].copy()
for c in DESCRIPTORS:
    cm_small[c] = pd.to_numeric(cm_small[c], errors="coerce")

# Drop any descriptor columns already present in `pairs` so the merge does not
# create clogp_x / clogp_y suffixes.
pairs = pairs.drop(columns=[c for c in DESCRIPTORS if c in pairs.columns])
df = pairs.merge(cm_small, on="compound_id", how="left")


# Numeric serum MIC: a censored ">100" is treated as its lower bound (100).
def parse_mic(s):
    s = str(s).strip()
    for r in (">=", "<=", ">", "<", "="):
        if s.startswith(r):
            s = s[len(r):]
            break
    try:
        return float(s)
    except ValueError:
        return np.nan

df["serum_mic_num"] = df["serum_mic_ugml"].apply(parse_mic)
df["log_serum_mic"] = np.log10(df["serum_mic_num"])
df["tolerant"] = (df["serum_active"] == "yes")

# --------------------------------------------------------------------------
# Per-descriptor statistics
# --------------------------------------------------------------------------
records = []
for d in DESCRIPTORS:
    sub = df[["tolerant", "log_serum_mic", d]].dropna()
    if sub[d].nunique() < 2:
        continue
    tol = sub.loc[sub["tolerant"], d]
    kil = sub.loc[~sub["tolerant"], d]
    # Spearman: descriptor vs serum MIC (positive rho => higher descriptor,
    # higher serum MIC => worse serum activity).
    rho, p_rho = stats.spearmanr(sub[d], sub["log_serum_mic"])
    # Mann-Whitney between groups (only if both groups present).
    if len(tol) >= 2 and len(kil) >= 2:
        u, p_u = stats.mannwhitneyu(tol, kil, alternative="two-sided")
    else:
        p_u = np.nan
    records.append({
        "descriptor": d,
        "median_tolerant": round(tol.median(), 3) if len(tol) else np.nan,
        "median_killed": round(kil.median(), 3) if len(kil) else np.nan,
        "n_tolerant": len(tol),
        "n_killed": len(kil),
        "spearman_rho_vs_serumMIC": round(rho, 3),
        "spearman_p": round(p_rho, 4),
        "mannwhitney_p": round(p_u, 4) if not np.isnan(p_u) else np.nan,
    })

stats_df = pd.DataFrame(records).sort_values("spearman_p")
stats_df.to_csv(os.path.join(OUT, "phase1_descriptor_stats.csv"), index=False)

# --------------------------------------------------------------------------
# Plot 1: clogP vs serum MIC (the headline hypothesis)
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4.5))
colors = df["tolerant"].map({True: "#2a9d8f", False: "#e76f51"})
ax.scatter(df["clogp"], df["serum_mic_num"], c=colors, s=60, edgecolor="k", linewidth=0.4)
ax.set_yscale("log")
ax.set_xlabel("clogP (lipophilicity)")
ax.set_ylabel("Serum MIC (ug/mL, log scale; >100 plotted at 100)")
ax.set_title("Phase 1: lipophilicity vs serum activity")
ax.axhline(50, color="grey", ls="--", lw=0.8)
ax.text(ax.get_xlim()[0], 55, "serum-active threshold (50)", fontsize=7, color="grey")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#2a9d8f", label="serum-tolerant"),
                   Patch(color="#e76f51", label="serum-killed")], fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase1_clogp_vs_serum_mic.png"), dpi=130)
plt.close(fig)

# --------------------------------------------------------------------------
# Plot 2: boxplots of the most-separating descriptors by group
# --------------------------------------------------------------------------
top = stats_df.head(4)["descriptor"].tolist()
fig, axes = plt.subplots(1, len(top), figsize=(3.2 * len(top), 4))
if len(top) == 1:
    axes = [axes]
for ax, d in zip(axes, top):
    sub = df[["tolerant", d]].dropna()
    data = [sub.loc[sub["tolerant"], d], sub.loc[~sub["tolerant"], d]]
    ax.boxplot(data, tick_labels=["tolerant", "killed"], widths=0.6)
    ax.set_title(d, fontsize=10)
fig.suptitle("Phase 1: descriptor distributions, serum-tolerant vs serum-killed")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase1_descriptor_boxplots.png"), dpi=130)
plt.close(fig)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
print("PHASE 1 — serum-shift SAR (hypothesis generation)")
print("=" * 60)
print(f"Matched-pair compounds analyzed: {len(df)} "
      f"({df['tolerant'].sum()} tolerant, {(~df['tolerant']).sum()} killed)")
print("\nDescriptors ranked by association with serum MIC (Spearman):\n")
print(stats_df.to_string(index=False))
print("\nWrote: phase1_descriptor_stats.csv, "
      "phase1_clogp_vs_serum_mic.png, phase1_descriptor_boxplots.png")
