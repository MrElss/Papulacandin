"""Smoke tests for the Stage 0 serum-shift pipeline.

Run with:  python -m pytest stage0/test_stage0.py -q
These are guardrail checks, not exhaustive unit tests: they lock in the facts a
downstream oracle relies on (24 anchor compounds, correct shift arithmetic,
censoring flags, no fabricated echinocandin pairs).
"""

import csv
import math
import os

import build_serum_shift_table as s0

HERE = os.path.dirname(os.path.abspath(__file__))


# --------------------------------------------------------------------------- #
# builders
# --------------------------------------------------------------------------- #
def test_papulacandin_has_24_anchor_compounds():
    pairs = s0.build_papulacandin_pairs()
    compounds = {p["compound_id"] for p in pairs}
    assert len(compounds) == 24, f"expected 24 anchor compounds, got {len(compounds)}"
    assert len(pairs) > 0


def test_serum_shift_arithmetic_is_correct():
    for p in s0.build_papulacandin_pairs():
        assert p["free_mic"] > 0 and p["serum_mic"] > 0
        expected = round(p["serum_mic"] / p["free_mic"], 4)
        assert math.isclose(p["serum_shift"], expected, rel_tol=1e-6, abs_tol=1e-4)
        assert math.isclose(
            p["log10_serum_shift"], round(math.log10(p["serum_shift"]), 4),
            rel_tol=1e-6, abs_tol=1e-4,
        )


def test_activity_ceiling_rule():
    # point 1: MIC > 100 ug/mL is a categorical INACTIVE call.
    assert s0.activity_call(50.0, "=") == "active"
    assert s0.activity_call(100.0, "=") == "inactive"
    assert s0.activity_call(250.0, ">") == "inactive"
    assert s0.activity_call(40.0, ">") == "ambiguous"   # tested ceiling < 100
    assert s0.activity_call(None, "=") == "unknown"


def test_serum_tolerance_call_is_consistent_with_activity():
    for p in s0.build_papulacandin_pairs():
        if p["serum_tolerance"] == "retained":
            assert p["free_activity"] == "active" and p["serum_activity"] == "active"
        elif p["serum_tolerance"] == "lost":
            assert p["free_activity"] == "active" and p["serum_activity"] == "inactive"
        elif p["serum_tolerance"] == "uninformative":
            assert p["free_activity"] != "active"
    # a meaningful continuous shift only when serum MIC is a real active number
    for p in s0.build_papulacandin_pairs():
        if p["serum_shift_meaningful"]:
            assert p["serum_activity"] == "active"


def test_invivo_proxy_excludes_cellular_only_and_flags_conflicts():
    pairs = s0.build_papulacandin_pairs()
    proxies = s0.build_invivo_proxies(s0._direct_calls_by_compound(pairs))
    assert proxies, "expected in-vivo efficacy proxies"
    # every proxy is in-vivo derived (never cellular-only)
    assert all(p["evidence"] == "in_vivo_efficacy" for p in proxies)
    # canonical echinocandin anchors are present (caspofungin etc.)
    names = {p["compound_name"].upper() for p in proxies}
    assert names & s0.CANONICAL_FKS
    # a proxy conflicting with a direct 'lost'/'mixed' call is down-weighted
    for p in proxies:
        if p["direct_pair_evidence"] in ("lost", "mixed"):
            assert p["confidence"] == "low"


def test_intrinsic_potency_covariate_present():
    # guardrail #1: potency is retained as a covariate for every pair with MW.
    pairs = s0.build_papulacandin_pairs()
    with_potency = [p for p in pairs if p["intrinsic_potency_pMIC"] != ""]
    assert len(with_potency) == len(pairs), "every anchor compound has a known MW"


def test_echinocandin_teacher_is_unpaired_and_has_no_fabricated_shift():
    context, fu = s0.build_echinocandin_teacher()
    assert context, "expected an echinocandin serum-context corpus"
    # the corpus must NOT expose a serum_shift column (no fabricated ratios)
    assert "serum_shift" not in context[0]
    assert {"serum_state", "source_reference"} <= set(context[0])
    # free-fraction seed carries only Fu / PPB endpoints
    assert all(r["endpoint_type"] in ("Fu", "PPB") for r in fu)


# --------------------------------------------------------------------------- #
# emitted artifacts (require a prior `python build_serum_shift_table.py` run)
# --------------------------------------------------------------------------- #
def _load(name):
    path = os.path.join(HERE, "outputs", name)
    if not os.path.exists(path):
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def test_unified_labels_combine_direct_and_proxy():
    rows = _load("serum_tolerance_labels.csv")
    if rows is None:
        return  # pipeline not yet run; builder tests already cover the logic
    sources = {r["label_source"] for r in rows}
    assert sources == {"direct_serum_pair", "in_vivo_proxy"}
    tol = {r["serum_tolerance"] for r in rows}
    assert tol == {"retained", "lost", "retained_presumed"}
    # in-vivo proxies bring the echinocandin anchors into the label set
    assert any(r["chemotype"] == "echinocandin" for r in rows)
