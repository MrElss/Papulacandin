"""Smoke tests for the Stage 1b free-fraction / PPB oracle.

Run:  python -m pytest stage1b/test_stage1b.py -q
These lock in the invariants a downstream stage relies on: consistent
featurization across chemotypes, an honest scaffold split (no leakage), a
distance-based applicability domain, and the headline finding that the
papulacandin glycolipids are out-of-domain.
"""

import os

import numpy as np
import pytest

import featurize as F
import train_free_fraction_oracle as s1b

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
# featurization
# --------------------------------------------------------------------------- #
def test_featurize_is_consistent_length_and_handles_large_molecules():
    small = F.featurize("CCO")
    # papulacandin B — a ~900 Da glycolipid must still featurize
    papu_b = ("CCCCCC/C=C/C=C/C(=O)O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)"
              "[C@H]1O")  # truncated stand-in; any large SMILES exercises the path
    big = F.featurize(papu_b)
    assert small is not None and big is not None
    assert len(small) == len(big) == len(F.FEATURE_NAMES)


def test_featurize_returns_none_on_garbage():
    assert F.featurize("not_a_smiles") is None
    assert F.featurize("") is None


# --------------------------------------------------------------------------- #
# scaffold split
# --------------------------------------------------------------------------- #
def test_scaffold_split_has_no_shared_scaffold_between_train_and_test():
    scaffolds = [f"S{i%7}" for i in range(50)]
    tr, te = s1b.scaffold_split_indices(scaffolds, frac_train=0.8)
    assert set(tr).isdisjoint(te)
    assert set(tr) | set(te) == set(range(50))
    train_scaf = {scaffolds[i] for i in tr}
    test_scaf = {scaffolds[i] for i in te}
    assert train_scaf.isdisjoint(test_scaf), "scaffold leaked across the split"


# --------------------------------------------------------------------------- #
# applicability domain
# --------------------------------------------------------------------------- #
def test_applicability_domain_flags_far_points_out():
    rng = np.random.RandomState(0)
    X = rng.normal(0, 1, size=(100, len(F.FEATURE_NAMES)))
    X[:, 0] = rng.uniform(200, 500, size=100)      # MW-like column
    ad = s1b.ApplicabilityDomain(X, mw_index=0, k=5)
    near = X[0].copy()
    far = X[0].copy()
    far[0] = 5000                                   # MW far outside training range
    assert ad.score(near)["in_domain"] is True
    assert ad.score(far)["in_domain"] is False
    assert ad.score(far)["mw_in_train_range"] is False


# --------------------------------------------------------------------------- #
# end-to-end artifacts (require a prior training run)
# --------------------------------------------------------------------------- #
def _load(name):
    import csv
    path = os.path.join(HERE, "outputs", name)
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_papulacandins_are_out_of_domain():
    rows = _load("free_fraction_predictions.csv")
    if rows is None:
        pytest.skip("run train_free_fraction_oracle.py first")
    papu = [r for r in rows if r["chemotype"] == "papulacandin"]
    in_dom = sum(1 for r in papu if r["in_domain"] == "True")
    # the decisive finding: the vast majority of papulacandins cannot be modelled
    assert in_dom / len(papu) < 0.1


def test_predictions_carry_uncertainty_and_ad_flag():
    rows = _load("free_fraction_predictions.csv")
    if rows is None:
        pytest.skip("run train_free_fraction_oracle.py first")
    r = rows[0]
    assert {"pred_fu", "tree_sd_log", "in_domain", "knn_distance"} <= set(r)
    assert 0.0 <= float(r["pred_fu"]) <= 1.0


def test_ad_is_a_valid_reliability_flag_not_a_chemotype_split():
    # the AD's success is that distance predicts error, and it is a size gradient:
    # approved echinocandins (large) are out-of-domain just like papulacandins.
    import analyze_applicability_domain as A
    v = A.ad_validity()
    assert v["in_mae"] < v["out_mae"], "in-domain should be more accurate"
    assert v["spearman"] > 0.15, "distance should track error"
    _, _, _, approved = A.size_breakdown()
    if approved:
        in_dom = sum(1 for r in approved if r["in_domain"] == "True")
        assert in_dom == 0, "approved echinocandins are large -> out-of-domain"
