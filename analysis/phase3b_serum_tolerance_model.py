#!/usr/bin/env python3
"""
phase3b_serum_tolerance_model.py
================================
PHASE 3B — Modelling serum tolerance (the project's goal-1 deliverable), done
HONESTLY for very small data.

HARD CONSTRAINT WE MUST RESPECT:
Of the 24 compounds with matched serum data, 21 share ONE scaffold (Fusacandin
A); the rest are singletons. So the data CANNOT support a scaffold-general
"serum tolerance" predictor -- holding out the Fusacandin-A scaffold removes
almost every positive example. Pretending otherwise (e.g. random cross-
validation across near-identical analogs) yields a fantasy AUC. The scientific
truth is narrower and still useful: WITHIN the fusacandin-A side-chain series we
can model how the C-6' side chain tunes serum MIC, and we can RANK analogs.

So this script:
  1. Restricts to the one modelable series (Fusacandin A scaffold, n~21).
  2. Fits a small, regularized regression of log10(serum MIC) on the Phase-1
     descriptors (interpretable coefficients, not a black box).
  3. Validates with leave-one-out CV and reports the rank correlation between
     predicted and observed serum MIC -- i.e. "can we rank which analog keeps
     serum activity?" -- plus the in-sample vs LOO gap (the optimism gap).
  4. States the applicability domain explicitly.

Outputs (analysis/outputs/):
  * phase3b_series_model_coeffs.csv   — standardized descriptor coefficients
  * phase3b_loo_predictions.csv       — observed vs LOO-predicted serum MIC
  * phase3b_pred_vs_obs.png
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import LeaveOneOut
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

DESCRIPTORS = ["mw_exact", "fsp3", "tpsa", "hba", "rotb", "clogp"]

# --------------------------------------------------------------------------
# Build the modelling frame: serum MIC + descriptors + scaffold
# --------------------------------------------------------------------------
pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
df = pairs.merge(cm[["compound_id", "parent_scaffold"] + DESCRIPTORS],
                 on="compound_id", suffixes=("", "_cm"))
for c in DESCRIPTORS:
    df[c] = pd.to_numeric(df[c], errors="coerce")


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

# Restrict to the single modelable scaffold series.
SERIES = "Fusacandin A"
series = df[df["parent_scaffold"] == SERIES].dropna(subset=DESCRIPTORS + ["log_serum_mic"]).copy()
print(f"Modelable series '{SERIES}': n={len(series)} "
      f"(of {len(df)} serum-labeled compounds; rest are scaffold singletons)")

X = series[DESCRIPTORS].values
ylog = series["log_serum_mic"].values

# --------------------------------------------------------------------------
# Leave-one-out CV (honest for n~21) vs in-sample (optimistic)
# --------------------------------------------------------------------------
def fit_model():
    return make_pipeline(StandardScaler(), Ridge(alpha=2.0))

# in-sample
m_full = fit_model().fit(X, ylog)
pred_in = m_full.predict(X)

# LOO
loo_pred = np.empty_like(ylog)
for tr, te in LeaveOneOut().split(X):
    loo_pred[te] = fit_model().fit(X[tr], ylog[tr]).predict(X[te])

rho_in, p_in = stats.spearmanr(ylog, pred_in)
rho_loo, p_loo = stats.spearmanr(ylog, loo_pred)

# --------------------------------------------------------------------------
# Standardized coefficients (which descriptor drives serum MIC?)
# --------------------------------------------------------------------------
ridge = m_full.named_steps["ridge"]
coeffs = pd.DataFrame({
    "descriptor": DESCRIPTORS,
    "std_coefficient": np.round(ridge.coef_, 4),
}).assign(direction=lambda d: np.where(
    d["std_coefficient"] > 0, "higher -> WORSE serum (higher MIC)",
    "higher -> BETTER serum (lower MIC)")
).sort_values("std_coefficient", key=np.abs, ascending=False)
coeffs.to_csv(os.path.join(OUT, "phase3b_series_model_coeffs.csv"), index=False)

# --------------------------------------------------------------------------
# Save LOO predictions
# --------------------------------------------------------------------------
out = series[["compound_id", "name", "serum_mic_ugml"]].copy()
out["obs_log_serum_mic"] = np.round(ylog, 3)
out["loo_pred_log_serum_mic"] = np.round(loo_pred, 3)
out = out.sort_values("loo_pred_log_serum_mic")
out.to_csv(os.path.join(OUT, "phase3b_loo_predictions.csv"), index=False)

# --------------------------------------------------------------------------
# Plot: observed vs LOO-predicted
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.scatter(ylog, loo_pred, s=60, edgecolor="k", color="#2a9d8f")
lim = [min(ylog.min(), loo_pred.min()) - 0.1, max(ylog.max(), loo_pred.max()) + 0.1]
ax.plot(lim, lim, "--", color="grey", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("Observed log10 serum MIC")
ax.set_ylabel("LOO-predicted log10 serum MIC")
ax.set_title(f"Phase 3B: within-{SERIES} serum-MIC ranking\n"
             f"LOO Spearman rho={rho_loo:.2f} (p={p_loo:.3f})")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "phase3b_pred_vs_obs.png"), dpi=130)
plt.close(fig)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
print("\nPHASE 3B — within-series serum-tolerance model")
print("=" * 60)
print(f"In-sample (optimistic) Spearman rho : {rho_in:.3f} (p={p_in:.3f})")
print(f"Leave-one-out (honest)  Spearman rho : {rho_loo:.3f} (p={p_loo:.3f})")
print("  -> the gap between these is the small-data optimism you must report.\n")
print("Standardized descriptor coefficients (what tunes serum MIC in this series):")
print(coeffs.to_string(index=False))
print("\nAPPLICABILITY DOMAIN: valid only for Fusacandin-A C-6' side-chain analogs.")
print("Do NOT extrapolate this model to other scaffolds -- no data supports that yet.")
print("\nWrote: phase3b_series_model_coeffs.csv, phase3b_loo_predictions.csv, phase3b_pred_vs_obs.png")
