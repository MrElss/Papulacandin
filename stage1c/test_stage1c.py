"""Smoke tests for the Stage 1c serum-tolerance oracle and gate.

Run:  python -m pytest stage1c/test_stage1c.py -q
These lock in the honest conclusions: the label set is one series, potency is a
strong confound, and the gate does NOT pass (so generation must not be blindly
gated on this oracle).
"""

import os

import build_serum_tolerance_oracle as s1c

HERE = os.path.dirname(os.path.abspath(__file__))


def _table():
    return s1c.build_table()


def test_labels_are_one_congeneric_series():
    rows = _table()
    n_fus = sum(1 for r in rows if "fusacandin" in r["name"].lower())
    # essentially a single series -> generalization untestable
    assert n_fus >= len(rows) - 2
    assert len(rows) >= 20


def test_potency_is_a_strong_confound():
    rows = _table()
    ev = s1c.evaluate(rows)
    # potency alone already predicts serum tolerance well (guardrail #1)
    assert ev["potency_only_loo_auc"] > 0.7
    # structure barely adds over potency
    assert ev["structure_adds_over_potency"] < 0.1


def test_gate_does_not_pass():
    rows = _table()
    ev = s1c.evaluate(rows)
    gate = s1c.gate_decision(ev)
    assert gate["GATE_PASSED"] is False
    # it fails for the right reasons: no chemotype diversity / confound not beaten
    assert not gate["chemotype_diverse_labels"] or not gate["structure_beats_potency_confound"]


def test_gate_logic_would_pass_on_ideal_evidence():
    # guard the decision rule itself: strong, diverse, potency-beating evidence passes
    ideal = {
        "full_bacc": 0.85, "full_loo_auc": 0.9, "structure_adds_over_potency": 0.2,
        "n_scaffolds": 8, "n_fusacandin": 5, "n": 40,
        "full_scaffold_grouped_auc": 0.8,
    }
    assert s1c.gate_decision(ideal)["GATE_PASSED"] is True


def test_predictions_written():
    path = os.path.join(HERE, "outputs", "serum_tolerance_predictions.csv")
    if not os.path.exists(path):
        return  # pipeline not yet run
    import csv
    rows = list(csv.DictReader(open(path)))
    assert rows and {"loo_pred_prob_retained", "label_retained"} <= set(rows[0])
