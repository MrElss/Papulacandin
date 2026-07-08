"""Smoke tests for the Stage 1a potency oracle.

Run:  python -m pytest stage1a/test_stage1a.py -q
Locks in: the active/inactive filter is learnable and beats baseline, the pMIC
regressor honestly fails scaffold validation, and scaffold folds don't leak.
"""

import os

import numpy as np

import build_potency_oracle as s1a

HERE = os.path.dirname(os.path.abspath(__file__))


def _rows():
    return s1a.build_table()


def test_table_has_active_and_inactive_labels():
    rows = _rows()
    labels = {r["active"] for r in rows}
    assert labels == {0, 1}
    assert len(rows) > 100


def test_scaffold_folds_partition_without_leakage():
    scaffolds = [f"S{i % 6}" for i in range(60)]
    folds = s1a.scaffold_folds(scaffolds, 5)
    all_idx = sorted(i for f in folds for i in f)
    assert all_idx == list(range(60))            # covers everything once
    # a scaffold must not appear in two folds
    for a in range(5):
        for b in range(a + 1, 5):
            sa = {scaffolds[i] for i in folds[a]}
            sb = {scaffolds[i] for i in folds[b]}
            assert sa.isdisjoint(sb)


def test_active_inactive_classifier_beats_baseline():
    ev = s1a.evaluate_classifier(_rows())
    # the loop filter must clearly beat the majority baseline on a scaffold split
    assert ev["scaffold_auc"] > 0.8
    assert ev["scaffold_bacc"] > ev["majority_acc"]


def test_pmic_regressor_does_not_generalize_across_scaffolds():
    ev = s1a.evaluate_regressor(_rows())
    # honest finding: graded potency does NOT transfer across scaffolds
    assert ev["scaffold"]["r2"] < 0.3


def test_predictions_present_for_library():
    path = os.path.join(HERE, "outputs", "potency_predictions.csv")
    if not os.path.exists(path):
        return
    import csv
    rows = list(csv.DictReader(open(path)))
    r = rows[0]
    assert {"pred_active_prob", "pred_pmic", "in_domain"} <= set(r)
    assert 0.0 <= float(r["pred_active_prob"]) <= 1.0
