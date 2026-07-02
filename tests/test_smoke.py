#!/usr/bin/env python3
"""Smoke tests for the Papulacandin / FKS1 serum-gap pipeline.

These are intentionally lightweight: they protect the *negative* scientific
result (documented in analysis/outputs/SYNTHESIS_phases1-11.md) against silent
bit-rot from dependency upgrades or data edits. They do NOT re-validate the
science — only that the curated data is well-formed and the fast, dependency-
light entry points of the pipeline still run end to end and emit the expected
schema on the real curated inputs.

Heavy phases (CREST/xtb/Gaussian QM, Vina docking) are deliberately out of
scope: they need external binaries and hours of compute. Their parsers are
exercised by phase6's own synthetic self-test (`python3 analysis/phase6_qm_layer.py`).

Run:  pytest -q            (from the repo root)
  or:  python3 tests/test_smoke.py
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ANALYSIS = ROOT / "analysis"
CORE = ROOT / "curated" / "core_tables"
OUT = ANALYSIS / "outputs"

# The matched serum-free/serum MIC set is the project's dependent variable and
# is fixed at 24 compounds (13 serum-tolerant, 11 serum-killed). If this count
# moves, a correlation-changing data edit happened and should be reviewed.
N_MATCHED_PAIRS = 24


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _run_phase(script: str) -> None:
    """Execute an analysis script and fail loudly with its output if it errors."""
    proc = subprocess.run(
        [sys.executable, str(ANALYSIS / script)],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    )
    assert proc.returncode == 0, (
        f"{script} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
    )


# --------------------------------------------------------------------------- #
# Curated-data integrity                                                       #
# --------------------------------------------------------------------------- #
def test_core_tables_present_and_nonempty():
    for name in ("compounds_master.csv", "activity_table.csv", "enzyme_assays.csv"):
        rows = _read_csv(CORE / name)
        assert rows, f"{name} has no data rows"


def test_compounds_master_schema():
    rows = _read_csv(CORE / "compounds_master.csv")
    required = {
        "compound_id",
        "canonical_name",
        "smiles_canonical",
        "inchikey",
        "mw_exact",
        "clogp",
        "tpsa",
    }
    missing = required - set(rows[0].keys())
    assert not missing, f"compounds_master.csv missing columns: {missing}"
    # compound_id must be unique and non-empty.
    ids = [r["compound_id"] for r in rows]
    assert all(ids), "blank compound_id present"
    assert len(ids) == len(set(ids)), "duplicate compound_id present"


def test_activity_table_schema():
    rows = _read_csv(CORE / "activity_table.csv")
    required = {"activity_id", "compound_id", "endpoint_type", "endpoint_value", "organism"}
    missing = required - set(rows[0].keys())
    assert not missing, f"activity_table.csv missing columns: {missing}"


def test_activity_compound_ids_resolve():
    compounds = {r["compound_id"] for r in _read_csv(CORE / "compounds_master.csv")}
    for r in _read_csv(CORE / "activity_table.csv"):
        cid = r["compound_id"]
        assert cid in compounds, f"activity row references unknown compound_id {cid}"


def test_canonical_smiles_and_inchikey_roundtrip():
    """A sample of curated SMILES must parse in RDKit and reproduce the stored
    InChIKey — guards against a toolchain change silently altering structures."""
    from rdkit import Chem  # imported lazily so data tests run without RDKit

    rows = [
        r
        for r in _read_csv(CORE / "compounds_master.csv")
        if r.get("smiles_canonical") and r.get("inchikey")
    ]
    assert rows, "no compounds carry both SMILES and InChIKey"
    checked = 0
    for r in rows[:15]:  # a representative sample keeps the test fast
        mol = Chem.MolFromSmiles(r["smiles_canonical"])
        assert mol is not None, f"{r['compound_id']}: SMILES failed to parse"
        got = Chem.MolToInchiKey(mol)
        assert got == r["inchikey"], (
            f"{r['compound_id']}: InChIKey mismatch stored={r['inchikey']} got={got}"
        )
        checked += 1
    assert checked >= 5


# --------------------------------------------------------------------------- #
# Pipeline entry points (fast, dependency-light phases)                        #
# --------------------------------------------------------------------------- #
def test_phase0_serum_gap_analysis():
    _run_phase("serum_gap_analysis.py")
    pairs = _read_csv(OUT / "serum_gap_pairs.csv")
    assert len(pairs) == N_MATCHED_PAIRS, (
        f"expected {N_MATCHED_PAIRS} matched pairs, got {len(pairs)}"
    )
    required = {"compound_id", "serumfree_mic_ugml", "serum_mic_ugml", "serum_shift_fold"}
    assert required <= set(pairs[0].keys())


def test_phase1_serum_shift_sar():
    _run_phase("phase1_serum_shift_sar.py")
    stats = _read_csv(OUT / "phase1_descriptor_stats.csv")
    assert stats, "phase1 produced no descriptor stats"
    required = {"descriptor", "spearman_rho_vs_serumMIC", "spearman_p", "n_tolerant", "n_killed"}
    assert required <= set(stats[0].keys())
    # Sanity: the split must still be 13 tolerant / 11 killed.
    assert stats[0]["n_tolerant"] == "13" and stats[0]["n_killed"] == "11"


def test_phase11_echinocandin_readacross():
    """Cross-chemotype read-across: the echinocandin serum-shift ordering pulled
    from the external FKS corpus must stay caspofungin < anidulafungin < micafungin
    (the robust, literature-consistent signal this phase rests on)."""
    _run_phase("phase11_echinocandin_readacross.py")
    ech = {r["compound"]: r for r in _read_csv(OUT / "phase11_echinocandin_serum_shift.csv")}
    for name in ("CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN"):
        assert name in ech, f"{name} missing from echinocandin read-across"
    folds = {n: float(ech[n]["serum_shift_fold"]) for n in ("CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN")}
    assert folds["CASPOFUNGIN"] < folds["ANIDULAFUNGIN"] <= folds["MICAFUNGIN"], (
        f"echinocandin serum-shift ordering changed: {folds}"
    )
    # The harmonized cross-chemotype table must carry both chemotypes.
    combined = _read_csv(OUT / "phase11_crosschemotype.csv")
    chemotypes = {r["chemotype"] for r in combined}
    assert "papulacandin" in chemotypes and "echinocandin" in chemotypes
    # 24 papulacandin rows should survive into the harmonized table.
    assert sum(r["chemotype"] == "papulacandin" for r in combined) == N_MATCHED_PAIRS


if __name__ == "__main__":
    # Allow running without pytest installed.
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as exc:  # noqa: BLE001 - report and continue
                failures += 1
                print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)
