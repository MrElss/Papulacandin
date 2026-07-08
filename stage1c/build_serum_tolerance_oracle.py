#!/usr/bin/env python3
"""Stage 1c — serum-tolerance oracle, and the GATE before generation.

START_HERE (§3, Stage 1c) envisioned a transfer-learning serum-shift oracle:
pretrain on echinocandin serum shifts, fine-tune on the papulacandin pairs, and
"do not proceed to generation until the combined serum-tolerance oracle is
validated on a scaffold-split hold-out."

The data found in Stage 0 forces an honest reframing:

  * There are NO echinocandin serum-shift pairs to pretrain on (serum-free and
    serum MICs never co-occur within a study). The envisioned echinocandin ->
    papulacandin transfer is not supported by the data.
  * The direct serum-tolerance labels are 24 compounds — and 23 of them are a
    single congeneric FUSACANDIN series (C-6' ester analogs from one study), the
    24th being papulacandin B. So out-of-chemotype generalization is essentially
    untestable (n = 1 outside the series).
  * The label is categorical (retained / lost) after the Stage-0 ceiling rule.

This script therefore does the scientifically defensible thing: build the best
classifier the data allows, then VALIDATE IT ADVERSARIALLY and render an explicit
gate decision. The key comparisons (guardrails #1, #2):

  * a POTENCY-ONLY baseline — because retained compounds here are more potent,
    a model must beat "just predict from potency" to add anything;
  * a MAJORITY-class baseline;
  * leave-one-out CV (within-series interpolation) AND scaffold-grouped CV
    (the closest thing to the required scaffold-split hold-out).

The output is a GATE DECISION plus the constructive within-series SAR hypothesis
(which C-6' esters track serum retention) as a testable lead — not a validated
objective to optimize blindly.

Run:
    python stage1c/build_serum_tolerance_oracle.py
"""

from __future__ import annotations

import csv
import os
import sys
import statistics
from collections import defaultdict, Counter

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import LeaveOneOut, LeaveOneGroupOut
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")
sys.path.insert(0, os.path.join(ROOT, "stage1b"))
import featurize as F  # noqa: E402

PAIRS = os.path.join(ROOT, "stage0/outputs/papulacandin_serum_pairs.csv")
COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")
PROXIES = os.path.join(ROOT, "stage0/outputs/invivo_serum_tolerance_proxies.csv")

# Small, interpretable feature set (24 points -> keep it lean). Potency is the
# confound we must beat; the rest capture the C-6' aromatic-ester lever.
STRUCT_FEATURES = ["aromatic_rings", "fsp3", "mw", "clogp", "rotb"]
SEED = 0


def _read(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


# --------------------------------------------------------------------------- #
# assemble the 24-compound modeling table
# --------------------------------------------------------------------------- #
def build_table():
    pairs = _read(PAIRS)
    cm = {r["compound_id"]: r for r in _read(COMPOUNDS)}
    by = defaultdict(list)
    for p in pairs:
        by[p["compound_id"]].append(p)

    rows = []
    for cid, ps in sorted(by.items()):
        inform = [p for p in ps if p["serum_tolerance"] in ("retained", "lost")]
        if not inform:
            continue
        n_ret = sum(1 for p in inform if p["serum_tolerance"] == "retained")
        label = 1 if n_ret > len(inform) - n_ret else 0        # 1 = retained
        pmics = [float(p["intrinsic_potency_pMIC"]) for p in inform
                 if p["intrinsic_potency_pMIC"]]
        pmic = statistics.median(pmics) if pmics else None
        c = cm[cid]
        smi = (c.get("smiles_canonical") or c.get("smiles_raw") or "").strip()
        feats = F.featurize(smi)
        if feats is None or pmic is None:
            continue
        fd = dict(zip(F.FEATURE_NAMES, feats))
        rows.append({
            "compound_id": cid,
            "name": c.get("canonical_name", ""),
            "label": label,                       # 1 retained, 0 lost
            "pmic": pmic,
            "scaffold": F.murcko_scaffold(smi),
            "modification": c.get("modification_summary", ""),
            **{f: fd[f] for f in STRUCT_FEATURES},
        })
    return rows


# --------------------------------------------------------------------------- #
# models + honest validation
# --------------------------------------------------------------------------- #
def _matrix(rows, cols):
    X = np.array([[r[c] for c in cols] for r in rows], float)
    y = np.array([r["label"] for r in rows])
    return X, y


def _loo_auc(X, y):
    """Leave-one-out out-of-fold probabilities -> ROC AUC + balanced accuracy."""
    oof = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(C=1.0, max_iter=1000,
                                               class_weight="balanced"))
        clf.fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])[:, 1]
    return (roc_auc_score(y, oof), balanced_accuracy_score(y, (oof >= 0.5).astype(int)),
            oof)


def _grouped_auc(X, y, groups):
    """Scaffold-grouped out-of-fold AUC (the gate's hard test)."""
    oof = np.zeros(len(y))
    logo = LeaveOneGroupOut()
    for tr, te in logo.split(X, y, groups):
        if len(set(y[tr])) < 2:                 # train fold must have both classes
            oof[te] = y[tr].mean()
            continue
        clf = make_pipeline(StandardScaler(),
                            LogisticRegression(C=1.0, max_iter=1000,
                                               class_weight="balanced"))
        clf.fit(X[tr], y[tr])
        oof[te] = clf.predict_proba(X[te])[:, 1]
    try:
        auc = roc_auc_score(y, oof)
    except ValueError:
        auc = float("nan")
    return auc, oof


def evaluate(rows):
    y = np.array([r["label"] for r in rows])
    groups = np.array([r["scaffold"] for r in rows])
    n_scaffolds = len(set(groups))
    n_fusacandin = sum(1 for r in rows if "fusacandin" in r["name"].lower())

    # baselines and models
    Xp, _ = _matrix(rows, ["pmic"])
    Xfull, _ = _matrix(rows, ["pmic"] + STRUCT_FEATURES)
    Xstruct, _ = _matrix(rows, STRUCT_FEATURES)

    maj = max(np.mean(y), 1 - np.mean(y))       # majority-class accuracy
    pot_auc, pot_bacc, _ = _loo_auc(Xp, y)
    full_auc, full_bacc, full_oof = _loo_auc(Xfull, y)
    struct_auc, struct_bacc, _ = _loo_auc(Xstruct, y)
    grouped_auc, _ = _grouped_auc(Xfull, y, groups)

    return {
        "n": len(y), "n_retained": int(y.sum()), "n_lost": int((1 - y).sum()),
        "n_fusacandin": n_fusacandin, "n_scaffolds": n_scaffolds,
        "majority_acc": round(maj, 3),
        "potency_only_loo_auc": round(pot_auc, 3),
        "potency_only_bacc": round(pot_bacc, 3),
        "structure_only_loo_auc": round(struct_auc, 3),
        "full_loo_auc": round(full_auc, 3),
        "full_bacc": round(full_bacc, 3),
        "full_scaffold_grouped_auc": round(grouped_auc, 3),
        "structure_adds_over_potency": round(full_auc - pot_auc, 3),
        "full_oof": full_oof,
    }


# --------------------------------------------------------------------------- #
# SAR hypothesis + proxy sanity check
# --------------------------------------------------------------------------- #
def sar_summary(rows):
    ret = [r for r in rows if r["label"] == 1]
    lost = [r for r in rows if r["label"] == 0]

    def med(rs, k):
        return round(statistics.median([r[k] for r in rs]), 2)

    return {
        "retained_median_aromatic_rings": med(ret, "aromatic_rings"),
        "lost_median_aromatic_rings": med(lost, "aromatic_rings"),
        "retained_median_fsp3": med(ret, "fsp3"),
        "lost_median_fsp3": med(lost, "fsp3"),
        "retained_median_pmic": med(ret, "pmic"),
        "lost_median_pmic": med(lost, "pmic"),
    }


# --------------------------------------------------------------------------- #
# gate decision
# --------------------------------------------------------------------------- #
def gate_decision(ev):
    """Explicit, rule-based verdict on whether to gate generation on this oracle."""
    beats_majority = ev["full_bacc"] > 0.5 and ev["full_loo_auc"] > 0.5
    beats_potency = ev["structure_adds_over_potency"] > 0.05
    scaffold_diverse = ev["n_scaffolds"] >= 5 and ev["n_fusacandin"] < ev["n"] - 3
    scaffold_validated = ev["full_scaffold_grouped_auc"] >= 0.7

    passed = beats_majority and beats_potency and scaffold_diverse and scaffold_validated
    return {
        "beats_majority_baseline": beats_majority,
        "structure_beats_potency_confound": beats_potency,
        "chemotype_diverse_labels": scaffold_diverse,
        "scaffold_split_validated": scaffold_validated,
        "GATE_PASSED": passed,
    }


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #
def write_report(rows, ev, sar, gate):
    L = []
    L.append("# Stage 1c — serum-tolerance oracle & generation gate\n\n")
    L.append(
        "**Endpoint:** retained vs lost (Stage-0 categorical serum-tolerance call). "
        "**Model:** L2 logistic regression, class-balanced. The point of this "
        "stage is the GATE decision, validated adversarially against a potency "
        "baseline.\n\n")

    L.append("## The label set is one congeneric series\n\n")
    L.append(
        f"- {ev['n']} compounds ({ev['n_retained']} retained, {ev['n_lost']} lost); "
        f"**{ev['n_fusacandin']} are fusacandin C-6' ester analogs** from one "
        f"study. Only {ev['n'] - ev['n_fusacandin']} sits outside that series.\n"
        f"- The envisioned echinocandin->papulacandin transfer is unsupported: no "
        f"echinocandin serum-shift labels exist (Stage 0). So this is a single "
        f"series, not a transfer problem.\n\n")

    L.append("## Validation (adversarial)\n\n")
    L.append("| model | LOO ROC-AUC | balanced acc |\n|---|---|---|\n")
    L.append(f"| majority baseline | – | {ev['majority_acc']} |\n")
    L.append(f"| **potency only** | {ev['potency_only_loo_auc']} | {ev['potency_only_bacc']} |\n")
    L.append(f"| structure only | {ev['structure_only_loo_auc']} | – |\n")
    L.append(f"| potency + structure | {ev['full_loo_auc']} | {ev['full_bacc']} |\n")
    L.append(
        f"\n- Structure adds **{ev['structure_adds_over_potency']:+.3f}** AUC over "
        f"potency alone (guardrail #1 — a serum model must beat the potency "
        f"confound).\n"
        f"- Scaffold-grouped AUC (the gate's hard test): "
        f"**{ev['full_scaffold_grouped_auc']}** — but across only "
        f"{ev['n_scaffolds']} generic scaffolds that all share the fusacandin core, "
        f"so this measures across-ester interpolation, not new-chemotype "
        f"generalization.\n\n")

    L.append("## Within-series SAR hypothesis (constructive)\n\n")
    L.append(
        f"- Retained analogs carry larger extended-aromatic C-6' esters "
        f"(median aromatic rings {sar['retained_median_aromatic_rings']} vs "
        f"{sar['lost_median_aromatic_rings']} for lost; fsp3 "
        f"{sar['retained_median_fsp3']} vs {sar['lost_median_fsp3']}) and are more "
        f"potent (pMIC {sar['retained_median_pmic']} vs {sar['lost_median_pmic']}).\n"
        f"- **Hypothesis (testable):** biphenyl / terphenyl / naphthoyl C-6' esters "
        f"favour serum retention in fusacandins. **Caveat:** confounded with "
        f"potency and unproven beyond this one series — a lead to TEST, not an "
        f"objective to optimize.\n\n")

    L.append("## GATE DECISION\n\n")
    for k, v in gate.items():
        if k == "GATE_PASSED":
            continue
        L.append(f"- {k}: **{v}**\n")
    verdict = "PASSED" if gate["GATE_PASSED"] else "NOT PASSED"
    L.append(f"\n### Gate: **{verdict}** — do "
             f"{'' if gate['GATE_PASSED'] else 'NOT '}proceed to unconstrained "
             f"generation on the basis of this oracle.\n\n")
    if not gate["GATE_PASSED"]:
        L.append(
            "The oracle is, at best, a within-fusacandin-series interpolator whose "
            "signal is entangled with potency; with essentially one scaffold it "
            "cannot be validated for new chemistry (guardrail #8: the binding "
            "constraint is DATA). **Recommendation:** do not gate generation on a "
            "'validated' serum oracle. Either (a) run the Stage-4 design-make-test "
            "loop first to generate a scaffold-diverse serum dataset, or (b) if "
            "generating now, keep it scaffold-constrained to the fusacandin series, "
            "treat the aromatic-ester hypothesis as an experiment, and confirm "
            "serum tolerance in the wet lab rather than trusting the score.\n")
    with open(os.path.join(OUT, "serum_tolerance_gate_report.md"), "w",
              encoding="utf-8") as fh:
        fh.write("".join(L))


def write_predictions(rows, ev):
    path = os.path.join(OUT, "serum_tolerance_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "compound_id", "name", "label_retained", "loo_pred_prob_retained",
            "pmic", "aromatic_rings", "fsp3", "modification"])
        w.writeheader()
        for r, p in zip(rows, ev["full_oof"]):
            w.writerow({
                "compound_id": r["compound_id"], "name": r["name"],
                "label_retained": r["label"],
                "loo_pred_prob_retained": round(float(p), 3),
                "pmic": round(r["pmic"], 2),
                "aromatic_rings": int(r["aromatic_rings"]),
                "fsp3": round(r["fsp3"], 3),
                "modification": r["modification"][:60],
            })


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = build_table()
    ev = evaluate(rows)
    sar = sar_summary(rows)
    gate = gate_decision(ev)
    write_report(rows, ev, sar, gate)
    write_predictions(rows, ev)

    print(f"labelled compounds: {ev['n']} ({ev['n_retained']} retained, "
          f"{ev['n_lost']} lost); {ev['n_fusacandin']} fusacandins")
    print(f"potency-only LOO-AUC: {ev['potency_only_loo_auc']}  |  "
          f"full LOO-AUC: {ev['full_loo_auc']}  "
          f"(structure adds {ev['structure_adds_over_potency']:+.3f})")
    print(f"scaffold-grouped AUC: {ev['full_scaffold_grouped_auc']} "
          f"across {ev['n_scaffolds']} (single-core) scaffolds")
    print(f"GATE PASSED: {gate['GATE_PASSED']}")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
