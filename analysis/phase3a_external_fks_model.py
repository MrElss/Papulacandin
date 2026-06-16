#!/usr/bin/env python3
"""
phase3a_external_fks_model.py
=============================
PHASE 3A — A well-powered "FKS-activity" classifier, transferred to papulacandins.

WHY THIS MODEL IS TRUSTWORTHY (unlike 3B):
The external pretraining set has 818 compounds (321 FKS / beta-1,3-glucan-synthase
POSITIVES vs 497 BACKGROUND decoys/exclusions). That is enough data to train a
real classifier. We learn "what does an FKS-active chemotype look like" from
ECFP4 fingerprints, then SCORE the curated papulacandins with it. High-scoring
papulacandins are the ones whose chemistry most resembles validated FKS actives
-> a prioritization signal for project goal #2 (find FKS1 inhibitors).

THE ONE METHODOLOGICAL TRAP WE AVOID:
Many compounds here are near-duplicate analogs (e.g. echinocandin series). If we
do random cross-validation, train and test folds contain near-identical molecules
and the AUC is falsely high. We therefore cross-validate with GroupKFold grouped
by Murcko scaffold, so an entire scaffold is held out at once. We report BOTH
random-CV and scaffold-CV so the optimism gap is visible.

Outputs (analysis/outputs/):
  * phase3a_cv_metrics.csv          — random vs scaffold CV (AUC, PR-AUC)
  * phase3a_papulacandin_scores.csv — FKS-likeness score for every papulacandin
  * phase3a_fks_model.joblib        — model trained on all 818 (for reuse)
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
import joblib

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
PRE = os.path.join(ROOT, "external", "data", "processed", "pretraining_v0_1")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

RANDOM_STATE = 42
N_FOLDS = 5
NBITS = 2048
RADIUS = 2

# --------------------------------------------------------------------------
# Load pretraining features (ECFP4) + labels, aligned by construction.
# --------------------------------------------------------------------------
npz = np.load(os.path.join(PRE, "fks_pretraining_ecfp4_matrix_v0_1.npz"),
              allow_pickle=True)
X = npz["ecfp4_bits"].astype(np.float32)        # (818, 2048)
y = npz["y"].astype(int)                          # (818,)
sw = npz["sample_weight"].astype(float)           # (818,)
pt = pd.read_csv(os.path.join(PRE, "fks_pretraining_dataset_v0_1.csv"))
# Group label for scaffold-aware CV (fall back to the compound id if missing).
groups = pt["murcko_scaffold_smiles"].fillna(pt["pretraining_compound_id"]).values

print(f"Pretraining matrix: {X.shape}, positives={y.sum()}, negatives={(y==0).sum()}")

# --------------------------------------------------------------------------
# Cross-validation: random vs scaffold-grouped
# --------------------------------------------------------------------------
def make_model():
    return RandomForestClassifier(
        n_estimators=400, max_depth=None, min_samples_leaf=2,
        n_jobs=-1, random_state=RANDOM_STATE, class_weight="balanced")

def run_cv(splits):
    aucs, prs = [], []
    for tr, te in splits:
        m = make_model()
        m.fit(X[tr], y[tr], sample_weight=sw[tr])
        p = m.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        prs.append(average_precision_score(y[te], p))
    return np.array(aucs), np.array(prs)

rand_splits = list(StratifiedKFold(N_FOLDS, shuffle=True, random_state=RANDOM_STATE).split(X, y))
scaf_splits = list(GroupKFold(N_FOLDS).split(X, y, groups))

rand_auc, rand_pr = run_cv(rand_splits)
scaf_auc, scaf_pr = run_cv(scaf_splits)

metrics = pd.DataFrame([
    {"cv": "random_kfold", "roc_auc_mean": rand_auc.mean(), "roc_auc_std": rand_auc.std(),
     "pr_auc_mean": rand_pr.mean(), "note": "optimistic: analogs leak across folds"},
    {"cv": "scaffold_groupkfold", "roc_auc_mean": scaf_auc.mean(), "roc_auc_std": scaf_auc.std(),
     "pr_auc_mean": scaf_pr.mean(), "note": "honest: whole scaffolds held out"},
])

# --- Triviality check: how well do SINGLE crude descriptors separate the
# classes? If a lone descriptor already scores high, the positive-vs-background
# split is "easy" (different chemical universes) rather than a subtle mechanism
# signal. This guards against over-interpreting a near-perfect ECFP4 AUC.
baseline_rows = []
for col in ["mol_wt_rdkit", "tpsa_rdkit", "heavy_atom_count_rdkit", "rotb_rdkit"]:
    v = pd.to_numeric(pt[col], errors="coerce")
    mask = v.notna().values
    auc = roc_auc_score(y[mask], v[mask].values)
    baseline_rows.append({"cv": f"single_descriptor::{col}",
                          "roc_auc_mean": max(auc, 1 - auc), "roc_auc_std": np.nan,
                          "pr_auc_mean": np.nan,
                          "note": "trivial baseline; high value => easy class split"})
metrics = pd.concat([metrics, pd.DataFrame(baseline_rows)], ignore_index=True)
metrics.to_csv(os.path.join(OUT, "phase3a_cv_metrics.csv"), index=False)

# --------------------------------------------------------------------------
# Train final model on ALL pretraining data and persist it.
# --------------------------------------------------------------------------
final = make_model()
final.fit(X, y, sample_weight=sw)
joblib.dump({"model": final, "nbits": NBITS, "radius": RADIUS}, os.path.join(OUT, "phase3a_fks_model.joblib"))

# --------------------------------------------------------------------------
# Featurize curated papulacandins the SAME way (Morgan r=2, 2048 bits) and score.
# --------------------------------------------------------------------------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)

def ecfp4(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return np.array(gen.GetFingerprint(mol), dtype=np.float32)

cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
feats, keep = [], []
for _, r in cm.iterrows():
    fp = ecfp4(r.get("smiles_canonical"))
    if fp is not None:
        feats.append(fp)
        keep.append(r)
Xp = np.vstack(feats)
scores = final.predict_proba(Xp)[:, 1]

papu = pd.DataFrame({
    "compound_id": [r["compound_id"] for r in keep],
    "name": [r["canonical_name"] for r in keep],
    "compound_class": [r["compound_class"] for r in keep],
    "fks_likeness_score": np.round(scores, 4),
})
papu = papu.sort_values("fks_likeness_score", ascending=False)
papu.to_csv(os.path.join(OUT, "phase3a_papulacandin_scores.csv"), index=False)

# --------------------------------------------------------------------------
# Console summary
# --------------------------------------------------------------------------
print("\nPHASE 3A — external FKS-activity classifier")
print("=" * 60)
print(metrics.to_string(index=False))
print(f"\nScored {len(papu)} curated papulacandins. Top FKS-like:")
print(papu.head(10).to_string(index=False))
print("\nWrote: phase3a_cv_metrics.csv, phase3a_papulacandin_scores.csv, phase3a_fks_model.joblib")
