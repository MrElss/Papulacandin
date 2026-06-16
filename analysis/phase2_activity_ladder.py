#!/usr/bin/env python3
"""
phase2_activity_ladder.py
=========================
PHASE 2 — The "activity ladder": where in the pipeline do compounds drop off?

A glucan-synthase inhibitor must survive three increasingly stringent stages
to be a useful drug:

  Stage 1  ENZYME    : inhibits beta-1,3-glucan synthase in vitro (IC50)
  Stage 2  CELL      : kills fungi in serum-FREE whole-cell assay (MIC)
  Stage 3  CELL+SERUM: still kills fungi in serum-CONTAINING assay (MIC)

The project's whole thesis is that papulacandins clear Stages 1-2 but fail at
Stage 3. This script makes that attrition explicit and quantitative, both:
  (a) per anchor compound that has data at multiple stages, and
  (b) at the dataset level (how many compounds remain "active" at each stage).

All potencies are converted to micromolar (uM) using the exact MW so enzyme
IC50 and whole-cell MIC are on the same axis. A compound is counted "active"
at a stage if its potency is <= ACTIVE_UM and the value is not censored upward.

Outputs (analysis/outputs/):
  * phase2_activity_ladder.csv  — per-compound potency at each stage (uM)
  * phase2_ladder_anchors.png   — line plot, anchor compounds across stages
  * phase2_attrition.png        — dataset-level attrition bar chart
"""

from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

ACTIVE_UM = 50.0  # potency threshold for "active" (uM)
ALBICANS = {"candida albicans", "c. albicans"}

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------
cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
mw = pd.to_numeric(cm.set_index("compound_id")["mw_exact"], errors="coerce")
name = cm.set_index("compound_id")["canonical_name"]
act = pd.read_csv(os.path.join(CORE, "activity_table.csv"))
enz = pd.read_csv(os.path.join(CORE, "enzyme_assays.csv"))


def ugml_to_um(value, cid):
    """ug/mL -> uM using exact MW (uM = ug/mL * 1000 / MW)."""
    m = mw.get(cid, np.nan)
    if pd.isna(m) or pd.isna(value):
        return np.nan
    return value * 1000.0 / m


def to_num(x):
    return pd.to_numeric(x, errors="coerce")


# --------------------------------------------------------------------------
# Stage 1 — enzyme IC50 (uM). Keep exact "=" values; convert ug/mL to uM.
# --------------------------------------------------------------------------
e = enz[(enz["endpoint_type"] == "IC50") & (enz["relation"].isin(["=", "~"]))].copy()
e["value"] = to_num(e["value"])
stage1 = {}
for cid, g in e.groupby("compound_id"):
    vals = []
    for _, r in g.iterrows():
        u = (r["value"] if str(r["unit"]).strip().lower() == "um"
             else ugml_to_um(r["value"], cid))
        if not pd.isna(u):
            vals.append(u)
    if vals:
        stage1[cid] = float(np.median(vals))

# --------------------------------------------------------------------------
# Stages 2 & 3 — serum-free / serum MIC (uM) against C. albicans.
# --------------------------------------------------------------------------
def is_serum(s):
    s = str(s).strip().lower()
    if s == "none":
        return False
    if any(k in s for k in ("serum", "protein", "fetal", "fbs", "fcs", "albumin")) \
       and "not reported" not in s:
        return True
    return None

mic = act[(act["endpoint_type"] == "MIC")
          & (act["organism"].str.strip().str.lower().isin(ALBICANS))
          & (act["unit"].str.strip().str.lower().isin(["ug/ml", "mcg/ml"]))].copy()
mic["value"] = to_num(mic["endpoint_value"])
mic["serum"] = mic["serum_or_protein"].apply(is_serum)

stage2, stage3 = {}, {}          # exact (uncensored) median potency, uM
stage2_cens, stage3_cens = {}, {}  # whether the representative value is censored
for cid, g in mic.groupby("compound_id"):
    for serum_flag, store, store_c in ((False, stage2, stage2_cens),
                                       (True, stage3, stage3_cens)):
        sub = g[g["serum"] == serum_flag]
        if sub.empty:
            continue
        med = float(sub["value"].median())
        # censored if every measurement uses ">" (i.e. inactive lower bound)
        cens = bool((sub["endpoint_relation"] == ">").all())
        store[cid] = ugml_to_um(med, cid)
        store_c[cid] = cens

# --------------------------------------------------------------------------
# Build per-compound ladder table
# --------------------------------------------------------------------------
all_ids = set(stage1) | set(stage2) | set(stage3)
rows = []
for cid in sorted(all_ids):
    rows.append({
        "compound_id": cid,
        "name": name.get(cid, ""),
        "enzyme_ic50_uM": round(stage1[cid], 3) if cid in stage1 else np.nan,
        "serumfree_mic_uM": round(stage2[cid], 3) if cid in stage2 else np.nan,
        "serumfree_censored": stage2_cens.get(cid, ""),
        "serum_mic_uM": round(stage3[cid], 3) if cid in stage3 else np.nan,
        "serum_censored": stage3_cens.get(cid, ""),
    })
ladder = pd.DataFrame(rows)
ladder.to_csv(os.path.join(OUT, "phase2_activity_ladder.csv"), index=False)

# --------------------------------------------------------------------------
# Plot A — anchor compounds (have >=2 stages) across the ladder
# --------------------------------------------------------------------------
anchors = ladder[(ladder[["enzyme_ic50_uM", "serumfree_mic_uM", "serum_mic_uM"]]
                  .notna().sum(axis=1) >= 2)]
stages = ["enzyme_ic50_uM", "serumfree_mic_uM", "serum_mic_uM"]
labels = ["Enzyme\nIC50", "Serum-free\nMIC", "Serum\nMIC"]
fig, ax = plt.subplots(figsize=(7, 5))
for _, r in anchors.iterrows():
    y = [r[s] for s in stages]
    ax.plot(range(3), y, marker="o", label=f"{r['compound_id']} {str(r['name'])[:18]}")
ax.set_yscale("log")
ax.set_xticks(range(3))
ax.set_xticklabels(labels)
ax.set_ylabel("Potency (uM, log scale; lower = more potent)")
ax.set_title("Phase 2: activity ladder for anchor compounds\n(rising line = losing potency)")
ax.axhline(ACTIVE_UM, color="grey", ls="--", lw=0.8)
ax.legend(fontsize=7, loc="lower right")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase2_ladder_anchors.png"), dpi=130)
plt.close(fig)

# --------------------------------------------------------------------------
# Plot B — dataset-level attrition: # compounds "active" at each stage
# --------------------------------------------------------------------------
def n_active(store, cens=None):
    n = 0
    for cid, v in store.items():
        if pd.isna(v) or v > ACTIVE_UM:
            continue
        if cens is not None and cens.get(cid, False):
            continue
        n += 1
    return n

counts = [
    n_active(stage1),
    n_active(stage2, stage2_cens),
    n_active(stage3, stage3_cens),
]
fig, ax = plt.subplots(figsize=(6, 4.5))
bars = ax.bar(labels, counts, color=["#264653", "#2a9d8f", "#e76f51"])
ax.bar_label(bars)
ax.set_ylabel(f"# compounds active (<= {ACTIVE_UM:g} uM)")
ax.set_title("Phase 2: attrition across assay stringency")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase2_attrition.png"), dpi=130)
plt.close(fig)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
print("PHASE 2 — activity ladder")
print("=" * 60)
print(f"Compounds with enzyme IC50              : {len(stage1)}")
print(f"Compounds with serum-free MIC           : {len(stage2)}")
print(f"Compounds with serum MIC                : {len(stage3)}")
print(f"\nActive (<= {ACTIVE_UM:g} uM) at each stage:")
print(f"  enzyme      : {counts[0]}")
print(f"  serum-free  : {counts[1]}")
print(f"  serum       : {counts[2]}   <-- the translational cliff")
print("\nAnchor compounds (data at >=2 stages):")
print(anchors[["compound_id", "name", "enzyme_ic50_uM",
               "serumfree_mic_uM", "serum_mic_uM"]].to_string(index=False))
print("\nWrote: phase2_activity_ladder.csv, "
      "phase2_ladder_anchors.png, phase2_attrition.png")
