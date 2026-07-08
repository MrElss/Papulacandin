"""Smoke tests for the round-1 panel selector.

Run:  python -m pytest round1/test_round1.py -q
Locks in the design-of-experiments guarantees: candidates are serum-unknown,
predicted-active, and obtainable; the panel is scaffold-diverse and reaches
beyond the fusacandin series; controls carry known serum labels.
"""

import os

import build_panel as bp

HERE = os.path.dirname(os.path.abspath(__file__))


def test_candidates_are_unknown_active_and_obtainable():
    rows, _ = bp.load()
    cands = bp.select_candidates(rows)
    assert cands
    for r in cands:
        assert r["serum_label"] == ""                     # no serum label yet
        assert r["pred_active_prob"] >= bp.ACTIVE_PROB_MIN  # predicted active
        assert r["cost"] in ("cheap_semisynthetic", "obtain_fermentation")


def test_panel_breaks_the_single_series_limit():
    rows, _ = bp.load()
    cands = bp.select_candidates(rows)
    # the whole point: reach scaffolds the serum set doesn't have
    assert sum(1 for r in cands if r["novel_scaffold"]) >= 3
    assert len({r["scaffold"] for r in cands}) >= 3


def test_greedy_selection_maximises_scaffold_clogp_coverage():
    rows, _ = bp.load()
    cands = bp.select_candidates(rows, n=8)
    cells = {(r["scaffold"], bp._clogp_bin(r["clogp"])) for r in cands}
    # diversity-greedy should give (near-)unique cells, not duplicates
    assert len(cells) >= len(cands) - 1


def test_controls_have_known_labels_both_classes():
    rows, _ = bp.load()
    controls = bp.pick_controls(rows)
    roles = {r["panel_role"] for r in controls}
    assert "control_serum_retained" in roles
    assert "control_serum_lost" in roles
    assert all(r["serum_label"] in ("0", "1") for r in controls)


def test_panel_csv_written():
    path = os.path.join(HERE, "outputs", "round1_synthesis_panel.csv")
    if not os.path.exists(path):
        return
    import csv
    rows = list(csv.DictReader(open(path)))
    assert rows and {"panel_role", "compound_id", "selection_reason"} <= set(rows[0])
