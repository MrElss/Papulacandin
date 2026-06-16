#!/usr/bin/env python3
"""
phase3_figures.py
=================
Presentation-grade figures and tables for Phase 3 (and a cross-phase summary).

Generates:
  * fig_phase3a_cv_bars.png      — random vs scaffold CV vs trivial baselines
  * fig_phase3a_chemspace.png    — t-SNE map: FKS+ / background / papulacandins
  * fig_phase3a_score_dist.png   — FKS-likeness score distributions by group
  * fig_phase3b_coeffs.png       — standardized descriptor coefficients
  * fig_phase3b_optimism.png     — in-sample vs LOO optimism gap
  * table_phase3_summary.csv/.md — one-glance metrics table for slides

Reuses artifacts written by phase3a/phase3b scripts; recomputes ECFP4 for the
chemical-space map. Run phase3a_* and phase3b_* first.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances
import joblib
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
PRE = os.path.join(ROOT, "external", "data", "processed", "pretraining_v0_1")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"figure.dpi": 130, "font.size": 10,
                     "axes.titlesize": 11, "axes.spines.top": False,
                     "axes.spines.right": False})
C_POS, C_BG, C_PAPU = "#1f77b4", "#bbbbbb", "#e63946"
NBITS, RADIUS = 2048, 2

# ==========================================================================
# Figure: Phase 3A cross-validation bar chart
# ==========================================================================
cv = pd.read_csv(os.path.join(OUT, "phase3a_cv_metrics.csv"))
labels = {"random_kfold": "Random CV\n(optimistic)",
          "scaffold_groupkfold": "Scaffold CV\n(honest)",
          "single_descriptor::mol_wt_rdkit": "MW only\n(baseline)",
          "single_descriptor::tpsa_rdkit": "TPSA only\n(baseline)"}
plot = cv[cv["cv"].isin(labels)].copy()
plot["lab"] = plot["cv"].map(labels)
colors = ["#2a9d8f", "#264653", "#adb5bd", "#adb5bd"]
fig, ax = plt.subplots(figsize=(6.5, 4.3))
bars = ax.bar(plot["lab"], plot["roc_auc_mean"], color=colors,
              yerr=plot["roc_auc_std"].fillna(0), capsize=4)
ax.bar_label(bars, fmt="%.3f", padding=3)
ax.set_ylim(0.5, 1.03)
ax.set_ylabel("ROC-AUC")
ax.set_title("Phase 3A: FKS classifier is near-perfect AND near-trivial\n"
             "(simple descriptors already reach ~0.70 -> classes are easy to split)")
ax.axhline(0.5, color="k", lw=0.6, ls=":")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase3a_cv_bars.png"))
plt.close(fig)

# ==========================================================================
# Recompute ECFP4 for papulacandins; load external fingerprints + labels
# ==========================================================================
npz = np.load(os.path.join(PRE, "fks_pretraining_ecfp4_matrix_v0_1.npz"), allow_pickle=True)
Xext = npz["ecfp4_bits"].astype(np.uint8)
yext = npz["y"].astype(int)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)
cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
papu_fps, papu_ids = [], []
for _, r in cm.iterrows():
    smi = r.get("smiles_canonical")
    if isinstance(smi, str) and smi.strip():
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            papu_fps.append(np.array(gen.GetFingerprint(mol), dtype=np.uint8))
            papu_ids.append(r["compound_id"])
Xpapu = np.vstack(papu_fps)

# ==========================================================================
# Figure: chemical-space t-SNE map (Tanimoto/Jaccard distance)
# ==========================================================================
Xall = np.vstack([Xext, Xpapu])
grp = (["pos"] * int((yext == 1).sum()) + ["bg"] * int((yext == 0).sum()))
# reorder external by label so grouping aligns
order = np.argsort(-yext)  # positives first
Xall = np.vstack([Xext[order], Xpapu])
grp = np.array(["FKS-positive"] * (yext == 1).sum()
               + ["background"] * (yext == 0).sum()
               + ["papulacandin (curated)"] * len(Xpapu))

D = pairwise_distances(Xall, metric="jaccard")
emb = TSNE(n_components=2, metric="precomputed", init="random",
           perplexity=30, random_state=42).fit_transform(D)

fig, ax = plt.subplots(figsize=(7, 5.6))
for g, c, z, s in [("background", C_BG, 1, 12),
                   ("FKS-positive", C_POS, 2, 18),
                   ("papulacandin (curated)", C_PAPU, 3, 34)]:
    m = grp == g
    ax.scatter(emb[m, 0], emb[m, 1], c=c, s=s, label=g, zorder=z,
               edgecolor="k" if g.startswith("papu") else "none", linewidth=0.3, alpha=0.85)
ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
ax.set_title("Phase 3A: papulacandins sit OUTSIDE the FKS-positive\n"
             "(echinocandin-dominated) chemical space -> out of applicability domain")
ax.legend(loc="best", fontsize=9)
ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase3a_chemspace.png"))
plt.close(fig)

# ==========================================================================
# Figure: FKS-likeness score distributions
# ==========================================================================
bundle = joblib.load(os.path.join(OUT, "phase3a_fks_model.joblib"))
model = bundle["model"]
score_ext = model.predict_proba(Xext.astype(np.float32))[:, 1]
score_papu = model.predict_proba(Xpapu.astype(np.float32))[:, 1]
fig, ax = plt.subplots(figsize=(6.5, 4.3))
bins = np.linspace(0, 1, 26)
ax.hist(score_ext[yext == 0], bins=bins, color=C_BG, alpha=0.8, label="external background")
ax.hist(score_ext[yext == 1], bins=bins, color=C_POS, alpha=0.7, label="external FKS-positive")
ax.hist(score_papu, bins=bins, color=C_PAPU, alpha=0.8, label="curated papulacandins")
ax.set_xlabel("FKS-likeness score (model P[positive])")
ax.set_ylabel("compound count")
ax.set_title("Phase 3A: papulacandins score LOW despite being real\nglucan-synthase inhibitors (domain mismatch, not inactivity)")
ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase3a_score_dist.png"))
plt.close(fig)

# ==========================================================================
# Figure: Phase 3B standardized coefficients
# ==========================================================================
co = pd.read_csv(os.path.join(OUT, "phase3b_series_model_coeffs.csv"))
co = co.iloc[::-1]
colors = ["#e76f51" if v > 0 else "#2a9d8f" for v in co["std_coefficient"]]
fig, ax = plt.subplots(figsize=(6.5, 4))
ax.barh(co["descriptor"], co["std_coefficient"], color=colors)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("Standardized coefficient on log10(serum MIC)")
ax.set_title("Phase 3B: drivers of serum MIC within the Fusacandin-A series\n"
             "(green = improves serum activity, red = worsens it)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase3b_coeffs.png"))
plt.close(fig)

# ==========================================================================
# Figure: Phase 3B optimism gap
# ==========================================================================
from scipy import stats as _st
loo = pd.read_csv(os.path.join(OUT, "phase3b_loo_predictions.csv"))
rho_loo, p_loo = _st.spearmanr(loo["obs_log_serum_mic"], loo["loo_pred_log_serum_mic"])
# recover in-sample rho from coeffs file note is not stored; recompute quickly
rho_in = 0.654  # reported by phase3b run; shown for the contrast
fig, ax = plt.subplots(figsize=(5.5, 4.3))
b = ax.bar(["In-sample\n(optimistic)", "Leave-one-out\n(honest)"],
           [rho_in, rho_loo], color=["#adb5bd", "#264653"])
ax.bar_label(b, fmt="%.2f")
ax.set_ylim(0, 1)
ax.set_ylabel("Spearman rho (pred vs observed serum MIC)")
ax.set_title("Phase 3B: small-data optimism gap\nhonest validation shows NO robust predictive power yet")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase3b_optimism.png"))
plt.close(fig)

# ==========================================================================
# Cross-phase summary table (for slides)
# ==========================================================================
summary = pd.DataFrame([
    {"phase": "1", "question": "Which properties separate serum-tolerant vs killed?",
     "key_result": "higher MW / lower fsp3 (rigid aromatic) track with serum tolerance; clogP does not",
     "status": "hypothesis (n=24, one series)"},
    {"phase": "2", "question": "Where do compounds fail?",
     "key_result": "potent on enzyme & serum-free cell, lost only in serum (e.g. PapB 1.7->111 uM)",
     "status": "robust per-compound ladder"},
    {"phase": "3A", "question": "Can the external FKS model rank papulacandins?",
     "key_result": "AUC 0.997 but trivial split; papulacandins out-of-domain (scores <=0.39)",
     "status": "honest negative / use as chemical-space ref"},
    {"phase": "3B", "question": "Can we predict serum MIC?",
     "key_result": f"within-series in-sample rho 0.65 -> LOO rho {rho_loo:.2f} (n.s.)",
     "status": "not yet predictive (data-limited)"},
])
summary.to_csv(os.path.join(OUT, "table_phase3_summary.csv"), index=False)
# Hand-roll a GitHub-flavored markdown table (avoids the optional 'tabulate' dep).
cols = list(summary.columns)
md = ["| " + " | ".join(cols) + " |",
      "| " + " | ".join("---" for _ in cols) + " |"]
for _, r in summary.iterrows():
    md.append("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |")
with open(os.path.join(OUT, "table_phase3_summary.md"), "w") as fh:
    fh.write("\n".join(md) + "\n")

print("Wrote 6 figures + summary table to analysis/outputs/:")
for f in ["fig_phase3a_cv_bars.png", "fig_phase3a_chemspace.png",
          "fig_phase3a_score_dist.png", "fig_phase3b_coeffs.png",
          "fig_phase3b_optimism.png", "table_phase3_summary.csv/.md"]:
    print("  -", f)
