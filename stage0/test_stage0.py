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


def test_censored_serum_mic_is_flagged_as_lower_bound():
    pairs = s0.build_papulacandin_pairs()
    censored = [p for p in pairs if p["serum_censored"]]
    assert censored, "expected some censored (>) serum MICs in the curated data"
    # a serum-only censored pair must be a lower bound on the true shift
    for p in censored:
        if not p["free_censored"]:
            assert p["shift_relation"] == ">="


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


def test_unified_table_holds_only_legitimate_matched_shifts():
    rows = _load("unified_serum_shift_table.csv")
    if rows is None:
        return  # pipeline not yet run; builder tests already cover the logic
    # today only the papulacandin class yields legitimate matched shifts
    assert {r["chemotype"] for r in rows} == {"papulacandin"}
    assert len({r["compound_id"] for r in rows}) == 24
