"""Smoke tests for the binary serum-activity label set.

Run:  python -m pytest stage0/test_binary_serum_activity.py -q
Locks in the honest design: direct serum-MIC evidence is trustworthy; papulacandin
in-vivo-only proxies are NOT training labels (they'd be likely-false positives).
"""

import os

import build_binary_serum_activity as b

HERE = os.path.dirname(os.path.abspath(__file__))


def test_active_call_uses_the_ceiling_rule():
    assert b._active_call(2.0, "=") == 1
    assert b._active_call(100.0, "=") == 0
    assert b._active_call(50.0, ">") == 0        # non-inhibitory at tested conc
    assert b._active_call(None, "=") is None


def test_observations_pool_multiple_sources():
    obs = (b.papulacandin_serum_mic() + b.echinocandin_serum_mic()
           + b.invivo_proxy_positives() + b.clinical_approved_positives())
    sources = {o["label_source"] for o in obs}
    assert {"papulacandin_serum_mic", "echinocandin_serum_mic",
            "in_vivo_proxy"} <= sources
    # direct serum MICs carry both classes of outcome
    papu = [o for o in obs if o["label_source"] == "papulacandin_serum_mic"]
    assert {o["serum_active"] for o in papu} == {0, 1}


def test_papulacandin_invivo_only_is_not_a_training_label():
    labels = b.consensus(
        b.papulacandin_serum_mic() + b.echinocandin_serum_mic()
        + b.invivo_proxy_positives() + b.clinical_approved_positives())
    weak = [r for r in labels if r["label_tier"] == "proxy_papulacandin_weak"]
    assert weak, "expected papulacandin in-vivo-only compounds"
    # none of them may be marked usable for training (would be false positives)
    assert all(r["usable_for_training"] is False for r in weak)


def test_direct_evidence_is_usable_and_balanced_enough():
    labels = b.consensus(
        b.papulacandin_serum_mic() + b.echinocandin_serum_mic()
        + b.invivo_proxy_positives() + b.clinical_approved_positives())
    train = [r for r in labels if r["usable_for_training"]]
    assert all(r["confidence"] in ("high", "medium") for r in train)
    # direct set has both actives and inactives (not a single-class degenerate set)
    direct = [r for r in labels if r["label_tier"] == "direct_serum"]
    calls = {r["serum_active"] for r in direct}
    assert calls == {0, 1}


def test_chembl_drop_in_is_optional():
    # pipeline must run whether or not the ChEMBL CSV is present
    obs, used = b.chembl_drop_in()
    assert used in (True, False)
    if not used:
        assert obs == []
