#!/usr/bin/env python3
"""Stage 1b — free-fraction / plasma-protein-binding (PPB) oracle.

START_HERE (§3, Stage 1b) calls this "the serum-tolerance workhorse ... genuinely
learnable": predict the unbound fraction (fu), because only free, un-degraded drug
acts. The endpoint is **log10(percent unbound)**, the label used by the public
training set; fu = 10**pred / 100.

What this script does, honoring the guardrails:

  1. Train a PPB regressor on a PUBLIC ADME dataset (Fang et al. 2023, Genentech
     "Computational-ADME", human PPB = LOG % unbound, 194 labelled drugs).
  2. Validate it the hard way — a SCAFFOLD-split hold-out, not a random split
     (guardrail #2: cheap validation overstates).
  3. Define an APPLICABILITY DOMAIN (standardized-descriptor kNN distance) and a
     size gate, and flag every prediction in/out of domain (guardrail #4:
     off-the-shelf QSAR does not transfer to this bRo5 chemotype).
  4. Ship predictions with UNCERTAINTY (random-forest tree spread) for the
     papulacandins and echinocandins.
  5. Anchor-check against the in-repo echinocandin PPB (anidulafungin).

The scientifically important output is not a single fu number for each design —
it is the AD verdict: the papulacandin glycolipids (median MW ~930) sit almost
entirely outside the training drugs' MW range (<=666), so the oracle must
EXTRAPOLATE for them. That is precisely why the program plans experimental fu
measurement (guardrails #5, #8); this stage quantifies where computation can and
cannot be trusted.

Run:
    python stage1b/train_free_fraction_oracle.py
Requires: scikit-learn, rdkit, numpy (see requirements.txt). The training CSV is
cached at stage1b/data/ADME_public_set_3521.csv (downloaded once from GitHub).
"""

from __future__ import annotations

import csv
import math
import os
import urllib.request

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

import featurize as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")
DATA = os.path.join(HERE, "data")

ADME_CSV = os.path.join(DATA, "ADME_public_set_3521.csv")
ADME_URL = ("https://raw.githubusercontent.com/molecularinformatics/"
            "Computational-ADME/main/ADME_public_set_3521.csv")
PPB_COL = "LOG PLASMA PROTEIN BINDING (HUMAN) (% unbound)"

CURATED_COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")
EXTERNAL_MATRIX = os.path.join(
    ROOT, "external/data/processed/external_fks_model_ready_v0_1/"
    "external_fks_descriptor_matrix_v0_1.csv")
ECHINO_PPB_SEED = os.path.join(ROOT, "stage0/outputs/echinocandin_free_fraction_seed.csv")

SEED = 0


# --------------------------------------------------------------------------- #
# data loading / featurization
# --------------------------------------------------------------------------- #
def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def ensure_training_csv():
    if os.path.exists(ADME_CSV):
        return
    os.makedirs(DATA, exist_ok=True)
    req = urllib.request.Request(ADME_URL, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        open(ADME_CSV, "wb").write(r.read())


def load_training():
    """Return (X, y, smiles, scaffolds) for the human-PPB labelled drugs."""
    rows = _read(ADME_CSV)
    X, y, smiles, scaffolds = [], [], [], []
    for r in rows:
        val = r.get(PPB_COL, "").strip()
        if not val:
            continue
        feats = F.featurize(r["SMILES"])
        if feats is None:
            continue
        X.append(feats)
        y.append(float(val))                 # log10(% unbound)
        smiles.append(r["SMILES"])
        scaffolds.append(F.murcko_scaffold(r["SMILES"]))
    return np.array(X), np.array(y), smiles, scaffolds


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
def scaffold_split_indices(scaffolds, frac_train=0.8):
    """Group by Murcko scaffold; whole scaffolds go to train or test (no leakage).

    Scaffolds are filled into the train set largest-first until the quota is met;
    the rest form the hold-out. Empty scaffolds are treated as singletons.
    """
    groups = {}
    for i, s in enumerate(scaffolds):
        groups.setdefault(s or f"__singleton_{i}", []).append(i)
    ordered = sorted(groups.values(), key=len, reverse=True)
    n_train_target = int(round(frac_train * len(scaffolds)))
    train_idx, test_idx = [], []
    for members in ordered:
        if len(train_idx) < n_train_target:
            train_idx.extend(members)
        else:
            test_idx.extend(members)
    return sorted(train_idx), sorted(test_idx)


def _metrics(y_true, y_pred):
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    mae = float(np.mean(np.abs(y_true - y_pred)))
    rho = float(spearmanr(y_true, y_pred).correlation)
    return {"n": int(len(y_true)), "r2": round(r2, 3),
            "mae_log": round(mae, 3), "spearman": round(rho, 3)}


def new_model():
    return RandomForestRegressor(
        n_estimators=400, min_samples_leaf=2, max_features="sqrt",
        random_state=SEED, n_jobs=-1)


# --------------------------------------------------------------------------- #
# applicability domain
# --------------------------------------------------------------------------- #
class ApplicabilityDomain:
    """Standardized-descriptor kNN distance + a training MW-range gate.

    A query is IN domain when its mean distance to its k nearest training points
    is within the training distribution (<= 95th percentile of train-to-train kNN
    distances) AND its MW is within the training MW range. Distance-based ADs are
    the honest way to say "this molecule is unlike anything the model was trained
    on" — critical here because random-forest error bars collapse toward the
    training mean under extrapolation and would otherwise look falsely confident.
    """

    def __init__(self, X_train, mw_index, k=5):
        self.scaler = StandardScaler().fit(X_train)
        self.Z = self.scaler.transform(X_train)
        self.k = k
        self.mw_index = mw_index
        self.mw_lo = float(np.min(X_train[:, mw_index]))
        self.mw_hi = float(np.max(X_train[:, mw_index]))
        # train-to-train kNN distance distribution (exclude self)
        d = self._knn_dist(self.Z, exclude_self=True)
        self.threshold = float(np.percentile(d, 95))

    def _knn_dist(self, Z, exclude_self=False):
        out = []
        for z in Z:
            dist = np.sqrt(((self.Z - z) ** 2).sum(axis=1))
            if exclude_self:
                dist = np.sort(dist)[1:self.k + 1]
            else:
                dist = np.sort(dist)[:self.k]
            out.append(float(dist.mean()))
        return np.array(out)

    def score(self, x):
        z = self.scaler.transform([x])[0]
        dist = np.sqrt(((self.Z - z) ** 2).sum(axis=1))
        knn = float(np.sort(dist)[:self.k].mean())
        mw = x[self.mw_index]
        in_mw = self.mw_lo <= mw <= self.mw_hi
        in_domain = (knn <= self.threshold) and in_mw
        return {
            "in_domain": bool(in_domain),
            "knn_distance": round(knn, 3),
            "knn_threshold": round(self.threshold, 3),
            "mw_in_train_range": bool(in_mw),
            "mw_train_range": f"{self.mw_lo:.0f}-{self.mw_hi:.0f}",
        }


# --------------------------------------------------------------------------- #
# prediction with uncertainty
# --------------------------------------------------------------------------- #
def predict_with_uncertainty(model, x):
    """log10(%unbound) prediction + tree-spread SD, converted to fu."""
    preds = np.array([est.predict([x])[0] for est in model.estimators_])
    mean_log = float(preds.mean())
    sd_log = float(preds.std())
    pct_unbound = min(100.0, 10 ** mean_log)     # cap at 100%
    return {
        "pred_log_pct_unbound": round(mean_log, 3),
        "pred_pct_unbound": round(pct_unbound, 4),
        "pred_fu": round(pct_unbound / 100.0, 5),
        "tree_sd_log": round(sd_log, 3),
    }


# --------------------------------------------------------------------------- #
# apply to project compounds
# --------------------------------------------------------------------------- #
def papulacandin_queries():
    out = []
    for r in _read(CURATED_COMPOUNDS):
        smi = (r.get("smiles_canonical") or r.get("smiles_raw") or "").strip()
        feats = F.featurize(smi)
        if feats is None:
            continue
        out.append({"chemotype": "papulacandin", "compound_id": r["compound_id"],
                    "name": r.get("canonical_name", ""), "smiles": smi, "x": feats})
    return out


def echinocandin_queries():
    out = []
    for r in _read(EXTERNAL_MATRIX):
        smi = (r.get("canonical_smiles_rdkit") or "").strip()
        feats = F.featurize(smi)
        if feats is None:
            continue
        out.append({"chemotype": "echinocandin", "compound_id": r["external_compound_id"],
                    "name": r.get("preferred_name", ""), "smiles": smi, "x": feats})
    return out


def anidulafungin_anchor():
    """Observed PPB for anidulafungin from the in-repo seed (external anchor)."""
    rows = _read(ECHINO_PPB_SEED)
    ppb = [float(r["value"]) for r in rows if r["endpoint_type"] == "PPB"]
    return float(np.mean(ppb)) if ppb else None


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    os.makedirs(OUT, exist_ok=True)
    ensure_training_csv()

    X, y, smiles, scaffolds = load_training()
    mw_index = F.FEATURE_NAMES.index("mw")

    # ---- scaffold-split validation (the honest test) ----
    tr, te = scaffold_split_indices(scaffolds, frac_train=0.8)
    m = new_model().fit(X[tr], y[tr])
    scaffold_metrics = _metrics(y[te], m.predict(X[te]))

    # ---- random 5-fold CV for comparison (optimistic baseline) ----
    from sklearn.model_selection import cross_val_predict, KFold
    cv = KFold(n_splits=5, shuffle=True, random_state=SEED)
    cv_pred = cross_val_predict(new_model(), X, y, cv=cv, n_jobs=-1)
    random_metrics = _metrics(y, cv_pred)

    # ---- final model on all training data + applicability domain ----
    final = new_model().fit(X, y)
    ad = ApplicabilityDomain(X, mw_index=mw_index, k=5)

    # ---- predict for project compounds ----
    queries = papulacandin_queries() + echinocandin_queries()
    pred_rows, in_dom = [], 0
    for q in queries:
        pu = predict_with_uncertainty(final, q["x"])
        adv = ad.score(q["x"])
        in_dom += int(adv["in_domain"])
        pred_rows.append({
            "chemotype": q["chemotype"], "compound_id": q["compound_id"],
            "name": q["name"], "mw": round(q["x"][mw_index], 1),
            **pu, **adv,
        })

    # ---- anchor check: anidulafungin ----
    anid = next((q for q in queries if q["name"].upper() == "ANIDULAFUNGIN"), None)
    anchor = {}
    if anid:
        pu = predict_with_uncertainty(final, anid["x"])
        obs_ppb = anidulafungin_anchor()
        anchor = {
            "predicted_pct_unbound": pu["pred_pct_unbound"],
            "predicted_pct_bound": round(100 - pu["pred_pct_unbound"], 2),
            "observed_pct_bound_in_repo": round(obs_ppb, 1) if obs_ppb else None,
            "in_domain": ad.score(anid["x"])["in_domain"],
        }

    _write_predictions(pred_rows)
    _write_ad_report(pred_rows)
    _write_metrics(scaffold_metrics, random_metrics, X, y, ad, pred_rows, anchor)

    print(f"training human-PPB drugs: {len(y)}")
    print(f"scaffold-split hold-out: R2={scaffold_metrics['r2']} "
          f"MAE(log)={scaffold_metrics['mae_log']} rho={scaffold_metrics['spearman']}")
    print(f"random 5-fold (optimistic): R2={random_metrics['r2']}")
    papu = [r for r in pred_rows if r["chemotype"] == "papulacandin"]
    papu_in = sum(1 for r in papu if r["in_domain"])
    print(f"papulacandins in applicability domain: {papu_in}/{len(papu)}")
    if anchor:
        print(f"anidulafungin anchor: pred {anchor['predicted_pct_bound']}% bound "
              f"vs observed {anchor['observed_pct_bound_in_repo']}% "
              f"(in_domain={anchor['in_domain']})")
    print(f"outputs -> {OUT}")


def _write_predictions(rows):
    path = os.path.join(OUT, "free_fraction_predictions.csv")
    fields = ["chemotype", "compound_id", "name", "mw", "pred_log_pct_unbound",
              "pred_pct_unbound", "pred_fu", "tree_sd_log", "in_domain",
              "knn_distance", "knn_threshold", "mw_in_train_range", "mw_train_range"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _write_ad_report(rows):
    path = os.path.join(OUT, "applicability_domain_report.csv")
    agg = {}
    for r in rows:
        c = r["chemotype"]
        a = agg.setdefault(c, {"chemotype": c, "n": 0, "in_domain": 0,
                               "mw_in_range": 0})
        a["n"] += 1
        a["in_domain"] += int(r["in_domain"])
        a["mw_in_range"] += int(r["mw_in_train_range"])
    for a in agg.values():
        a["pct_in_domain"] = round(100 * a["in_domain"] / a["n"], 1)
        a["pct_mw_in_range"] = round(100 * a["mw_in_range"] / a["n"], 1)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["chemotype", "n", "in_domain",
                           "pct_in_domain", "mw_in_range", "pct_mw_in_range"])
        w.writeheader()
        w.writerows(agg.values())


def _write_metrics(scaffold_m, random_m, X, y, ad, pred_rows, anchor):
    papu = [r for r in pred_rows if r["chemotype"] == "papulacandin"]
    echi = [r for r in pred_rows if r["chemotype"] == "echinocandin"]
    papu_in = sum(1 for r in papu if r["in_domain"])
    echi_in = sum(1 for r in echi if r["in_domain"])
    lines = []
    lines.append("# Stage 1b — free-fraction / PPB oracle\n\n")
    lines.append(
        "Endpoint: **log10(% unbound)** (fu = 10^pred / 100). Trained on the "
        "public Fang et al. 2023 Computational-ADME human-PPB set.\n\n")
    lines.append(f"## Training data\n\n- {len(y)} human-PPB labelled drugs "
                 f"(MW {ad.mw_lo:.0f}-{ad.mw_hi:.0f}).\n\n")
    lines.append("## Validation (guardrail #2: validate the hard way)\n\n")
    lines.append(
        f"- **Scaffold-split hold-out (honest):** R² = {scaffold_m['r2']}, "
        f"MAE = {scaffold_m['mae_log']} log-units, Spearman = "
        f"{scaffold_m['spearman']} (n = {scaffold_m['n']}).\n"
        f"- Random 5-fold (optimistic, for contrast): R² = {random_m['r2']}, "
        f"MAE = {random_m['mae_log']}. The gap is the usual scaffold penalty — "
        f"trust the scaffold number.\n\n")
    if anchor:
        lines.append("## Anchor check — anidulafungin\n\n")
        lines.append(
            f"- Predicted {anchor['predicted_pct_bound']}% bound "
            f"({anchor['predicted_pct_unbound']}% unbound); in-repo observed "
            f"~{anchor['observed_pct_bound_in_repo']}% bound. In domain: "
            f"{anchor['in_domain']}. (Literature anidulafungin PPB ≈ 99%; a large "
            f"echinocandin, so an out-of-domain query — read with the AD flag.)\n\n")
    lines.append("## Applicability domain — the decisive result\n\n")
    lines.append(
        f"- **Papulacandins in domain: {papu_in}/{len(papu)} "
        f"({100*papu_in/len(papu):.0f}%).**\n"
        f"- Echinocandins/FKS in domain: {echi_in}/{len(echi)} "
        f"({100*echi_in/len(echi):.0f}%).\n\n")
    lines.append(
        "The papulacandin glycolipids (median MW ~930) sit almost entirely "
        "outside the training drugs' MW range (≤" f"{ad.mw_hi:.0f}), so the oracle "
        "must **extrapolate** for them and its predictions there are not "
        "trustworthy — exactly guardrail #4 (off-the-shelf QSAR does not transfer "
        "to this bRo5 chemotype). This operationalizes guardrails #5 and #8: for "
        "the real design targets, free fraction must be **measured** "
        "(equilibrium dialysis), not predicted. Use the oracle to *rank in-domain* "
        "analogs and to *flag* which designs even fall in a modellable region — "
        "not to assign an fu to a papulacandin.\n\n")
    lines.append("## Files\n\n"
                 "- `free_fraction_predictions.csv` — per-compound fu + uncertainty + AD flag.\n"
                 "- `applicability_domain_report.csv` — in-domain fractions by chemotype.\n")
    with open(os.path.join(OUT, "ppb_oracle_metrics.md"), "w", encoding="utf-8") as fh:
        fh.write("".join(lines))


if __name__ == "__main__":
    main()
