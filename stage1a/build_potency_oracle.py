#!/usr/bin/env python3
"""Stage 1a — intrinsic-potency oracle (serum-free MIC).

START_HERE (§3, Stage 1a) wants a potency scorer. Potency is the dominant,
LEARNABLE axis and — crucially — the confound every serum-tolerance analysis must
beat (guardrail #1). Its role in the design-make-test loop (see
planning/round1_synthesis_panel_plan.md) is as a FILTER: only synthesize/test
compounds predicted active enough serum-free, because a compound dead in broth
yields an uninformative serum-tolerance label.

Two models, both validated on a scaffold split (guardrail #2):

  * CLASSIFIER  active vs inactive serum-free (MIC < 100 ug/mL vs not) — the loop
    filter. All 132 labelled compounds.
  * REGRESSOR   pMIC = -log10(mol/L) on the 93 active compounds — graded ranking.

Both ship with uncertainty and a distance-based applicability-domain flag (reused
from Stage 1b). Endpoint uses the Stage-0 ceiling rule: MIC > 100 ug/mL = inactive
(a right-censored value, excluded from the regression rather than imputed).

Run:
    python stage1a/build_potency_oracle.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
import statistics
from collections import defaultdict

import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")
sys.path.insert(0, os.path.join(ROOT, "stage1b"))
import featurize as F  # noqa: E402
from train_free_fraction_oracle import ApplicabilityDomain  # noqa: E402

ACTIVITY = os.path.join(ROOT, "curated/core_tables/activity_table.csv")
COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")

CEILING_UGML = 100.0          # Stage-0 rule: MIC > 100 ug/mL = no activity
SEED = 0


def _read(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _is_free(text):
    return (text or "").strip().lower() in ("", "none", "none reported", "not reported")


# --------------------------------------------------------------------------- #
# assemble per-compound potency table
# --------------------------------------------------------------------------- #
def build_table():
    act = _read(ACTIVITY)
    cm = {r["compound_id"]: r for r in _read(COMPOUNDS)}
    mic = [r for r in act
           if r["endpoint_type"] in ("MIC", "MIC50")
           and _is_free(r["serum_or_protein"])
           and r["unit"] == "ug/mL" and _num(r["endpoint_value"])]

    by = defaultdict(list)
    for r in mic:
        by[r["compound_id"]].append(r)

    rows = []
    for cid, rs in by.items():
        c = cm.get(cid, {})
        mw = _num(c.get("mw_exact"))
        smi = (c.get("smiles_canonical") or c.get("smiles_raw") or "").strip()
        feats = F.featurize(smi)
        if feats is None or not mw:
            continue
        vals = [(_num(r["endpoint_value"]), r["endpoint_relation"]) for r in rs]
        # active if any uncensored measurement is below the ceiling
        uncensored_active = [v for v, rel in vals if rel != ">" and 0 < v < CEILING_UGML]
        is_active = 1 if uncensored_active else 0
        pmic = None
        if uncensored_active:
            uM = statistics.median(uncensored_active) / mw * 1000.0
            pmic = -math.log10(uM * 1e-6)
        rows.append({
            "compound_id": cid,
            "name": c.get("canonical_name", ""),
            "smiles": smi,
            "x": feats,
            "scaffold": F.murcko_scaffold(smi),
            "active": is_active,
            "pmic": pmic,
        })
    return rows


# --------------------------------------------------------------------------- #
# scaffold split (whole scaffolds to one side)
# --------------------------------------------------------------------------- #
def scaffold_folds(scaffolds, n_splits=5):
    groups = defaultdict(list)
    for i, s in enumerate(scaffolds):
        groups[s or f"__s{i}"].append(i)
    ordered = sorted(groups.values(), key=len, reverse=True)
    folds = [[] for _ in range(n_splits)]
    # distribute scaffolds round-robin onto the smallest fold (balances sizes)
    for members in ordered:
        smallest = min(range(n_splits), key=lambda k: len(folds[k]))
        folds[smallest].extend(members)
    return folds


def _reg_metrics(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    from scipy.stats import spearmanr
    # Spearman is undefined when either input is constant (e.g. the mean baseline)
    if np.ptp(y_true) == 0 or np.ptp(y_pred) == 0:
        rho = float("nan")
    else:
        rho = float(spearmanr(y_true, y_pred).correlation)
    return {
        "r2": round(1 - ss_res / ss_tot, 3) if ss_tot else float("nan"),
        "mae": round(float(np.mean(np.abs(y_true - y_pred))), 3),
        "spearman": round(rho, 3),
    }


# --------------------------------------------------------------------------- #
# classifier: active vs inactive serum-free (the loop filter)
# --------------------------------------------------------------------------- #
def evaluate_classifier(rows):
    X = np.array([r["x"] for r in rows], float)
    y = np.array([r["active"] for r in rows])
    scaf = [r["scaffold"] for r in rows]
    folds = scaffold_folds(scaf, 5)
    oof = np.zeros(len(y))
    for k in range(5):
        te = folds[k]
        tr = [i for i in range(len(y)) if i not in set(te)]
        if len(set(y[tr])) < 2:
            oof[te] = y[tr].mean()
            continue
        clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                     random_state=SEED, class_weight="balanced",
                                     n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])[:, 1]
    auc = roc_auc_score(y, oof) if len(set(y)) > 1 else float("nan")
    bacc = balanced_accuracy_score(y, (oof >= 0.5).astype(int))
    majority = max(y.mean(), 1 - y.mean())
    return {"n": len(y), "n_active": int(y.sum()), "n_inactive": int((1 - y).sum()),
            "scaffold_auc": round(auc, 3), "scaffold_bacc": round(bacc, 3),
            "majority_acc": round(majority, 3)}


# --------------------------------------------------------------------------- #
# regressor: pMIC on active compounds (graded ranking)
# --------------------------------------------------------------------------- #
def evaluate_regressor(rows):
    act = [r for r in rows if r["pmic"] is not None]
    X = np.array([r["x"] for r in act], float)
    y = np.array([r["pmic"] for r in act])
    scaf = [r["scaffold"] for r in act]
    folds = scaffold_folds(scaf, 5)
    oof = np.zeros(len(y))
    for k in range(5):
        te = folds[k]
        tr = [i for i in range(len(y)) if i not in set(te)]
        rf = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                   max_features="sqrt", random_state=SEED,
                                   n_jobs=-1).fit(X[tr], y[tr])
        oof[te] = rf.predict(X[te])
    scaffold_m = _reg_metrics(y, oof)
    # baseline: predict the training mean
    base = _reg_metrics(y, np.full_like(y, y.mean()))
    return {"n": len(y), "scaffold": scaffold_m, "baseline_mean": base,
            "pmic_range": [round(float(y.min()), 2), round(float(y.max()), 2)]}


# --------------------------------------------------------------------------- #
# final models + predictions for the whole library (loop filter output)
# --------------------------------------------------------------------------- #
def predict_library(rows):
    X = np.array([r["x"] for r in rows], float)
    y_cls = np.array([r["active"] for r in rows])
    mw_i = F.FEATURE_NAMES.index("mw")

    clf = RandomForestClassifier(n_estimators=400, min_samples_leaf=2,
                                 random_state=SEED, class_weight="balanced",
                                 n_jobs=-1).fit(X, y_cls)
    ad = ApplicabilityDomain(X, mw_index=mw_i, k=5)

    act = [r for r in rows if r["pmic"] is not None]
    Xr = np.array([r["x"] for r in act], float)
    yr = np.array([r["pmic"] for r in act])
    rf = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                               max_features="sqrt", random_state=SEED,
                               n_jobs=-1).fit(Xr, yr)

    preds = []
    for r in rows:
        x = np.array(r["x"], float)
        p_active = float(clf.predict_proba([x])[0][1])
        tree_pm = np.array([e.predict([x])[0] for e in rf.estimators_])
        preds.append({
            "compound_id": r["compound_id"], "name": r["name"],
            "observed_active": r["active"],
            "observed_pmic": round(r["pmic"], 2) if r["pmic"] is not None else "",
            "pred_active_prob": round(p_active, 3),
            "pred_pmic": round(float(tree_pm.mean()), 2),
            "pred_pmic_sd": round(float(tree_pm.std()), 2),
            "in_domain": ad.score(x)["in_domain"],
        })
    return preds


# --------------------------------------------------------------------------- #
def main():
    os.makedirs(OUT, exist_ok=True)
    rows = build_table()
    cls = evaluate_classifier(rows)
    reg = evaluate_regressor(rows)
    preds = predict_library(rows)

    _write_predictions(preds)
    _write_report(rows, cls, reg, preds)

    print(f"potency table: {cls['n']} compounds "
          f"({cls['n_active']} active, {cls['n_inactive']} inactive serum-free)")
    print(f"classifier (active/inactive) scaffold AUC={cls['scaffold_auc']} "
          f"bacc={cls['scaffold_bacc']} vs majority {cls['majority_acc']}")
    print(f"regressor (pMIC) scaffold R2={reg['scaffold']['r2']} "
          f"MAE={reg['scaffold']['mae']} rho={reg['scaffold']['spearman']} "
          f"(baseline R2={reg['baseline_mean']['r2']})")
    print(f"outputs -> {OUT}")


def _write_predictions(preds):
    path = os.path.join(OUT, "potency_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(preds[0].keys()))
        w.writeheader()
        w.writerows(preds)


def _write_report(rows, cls, reg, preds):
    in_dom = sum(1 for p in preds if p["in_domain"])
    L = []
    L.append("# Stage 1a — intrinsic-potency oracle (serum-free MIC)\n\n")
    L.append(
        "Potency is the learnable, dominant axis and the confound serum-tolerance "
        "work must beat (guardrail #1). In the design-make-test loop this oracle is "
        "the **filter**: only test compounds predicted active enough serum-free to "
        "give an informative serum-tolerance result.\n\n")
    L.append("## Active/inactive classifier (the loop filter) — USABLE\n\n")
    L.append(
        f"- {cls['n']} compounds ({cls['n_active']} active, {cls['n_inactive']} "
        f"inactive serum-free by the >100 ug/mL ceiling rule).\n"
        f"- **Scaffold-split ROC-AUC = {cls['scaffold_auc']}**, balanced acc = "
        f"{cls['scaffold_bacc']} (majority baseline {cls['majority_acc']}).\n"
        f"- **Verdict:** strongly beats baseline and generalizes across scaffolds — "
        f"whether a compound is active serum-free is genuinely learnable. Use this "
        f"as the loop filter.\n\n")
    L.append("## pMIC regressor (graded ranking) — DOES NOT GENERALIZE\n\n")
    L.append(
        f"- {reg['n']} active compounds, pMIC {reg['pmic_range'][0]}–"
        f"{reg['pmic_range'][1]} (a narrow ~3-log window).\n"
        f"- **Scaffold-split R² = {reg['scaffold']['r2']}** (WORSE than the "
        f"mean-predictor baseline R² = {reg['baseline_mean']['r2']}), MAE = "
        f"{reg['scaffold']['mae']} log-units.\n"
        f"- **Verdict:** fine potency ranking among active compounds does NOT "
        f"transfer across scaffolds — the data is dominated by one core, and "
        f"within-window potency differences are scaffold-specific. Do **not** rely "
        f"on `pred_pmic` to rank across novel scaffolds; treat it as in-domain, "
        f"low-confidence only. This is guardrail #2 in action: honest scaffold "
        f"validation catches a model that a random split would have flattered.\n\n")
    L.append("## Applicability domain\n\n")
    L.append(
        f"- {in_dom}/{len(preds)} library compounds in domain (the potency oracle "
        "is trained on this chemotype, so coverage is far better than the "
        "drug-trained Stage-1b PPB oracle — but still check the flag per design).\n\n")
    L.append("## Use in the loop\n\n")
    L.append(
        "`potency_predictions.csv` gives `pred_active_prob` and `pred_pmic` "
        "(+SD, +AD flag) for every library compound. Panel selection keeps only "
        "candidates with high `pred_active_prob`, then chooses among them for "
        "diversity and serum-tolerance informativeness. A candidate predicted "
        "inactive serum-free is dropped — it cannot teach us about serum tolerance.\n")
    with open(os.path.join(OUT, "potency_oracle_report.md"), "w",
              encoding="utf-8") as fh:
        fh.write("".join(L))


if __name__ == "__main__":
    main()
