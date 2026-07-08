"""Tests for the ChEMBL export converter.

Run:  python -m pytest stage0/test_chembl_convert.py -q
Guards the classification rules — especially that serum-FREE assays are not
misread as serum context, and that the ceiling rule sets active/inactive.
"""

import os
import tempfile

import convert_chembl_export as cx


WEB_CSV = (
    "Molecule ChEMBL ID;Molecule Name;Standard Type;Standard Relation;"
    "Standard Value;Standard Units;Assay Description;Document ChEMBL ID\n"
    "CHEMBL1;CASPOFUNGIN;MIC;=;0.5;ug.mL-1;MIC in 50% human serum;DOC1\n"
    "CHEMBL1;CASPOFUNGIN;MIC;>;100;ug.mL-1;MIC in 50% mouse serum;DOC2\n"
    "CHEMBL2;ANIDULAFUNGIN;MIC;=;0.03;ug.mL-1;MIC in RPMI, no serum;DOC3\n"
    "CHEMBL2;ANIDULAFUNGIN;PPB;=;99;%;Plasma protein binding human;DOC4\n"
)


def _write_tmp(text):
    fd, path = tempfile.mkstemp(suffix=".csv")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def test_serum_free_assay_is_not_serum_context():
    assert cx._is_serum("MIC in 50% human serum") is True
    assert cx._is_serum("MIC in RPMI, no serum") is False
    assert cx._is_serum("antifungal activity, serum-free medium") is False
    assert cx._is_serum("plasma protein binding") is True


def test_ceiling_rule_sets_active_inactive():
    path = _write_tmp(WEB_CSV)
    try:
        serum, ppb = cx.convert([path])
    finally:
        os.remove(path)
    calls = {(r["compound_name"], r["serum_active"]) for r in serum}
    # 0.5 ug/mL in serum -> active; >100 -> inactive; the no-serum row excluded
    assert ("CASPOFUNGIN", 1) in calls
    assert ("CASPOFUNGIN", 0) in calls
    assert all("ANIDULAFUNGIN" != r["compound_name"] or r["serum_active"] in (0, 1)
               for r in serum)
    # the serum-free anidulafungin MIC must not appear as a serum label
    assert not any(r["compound_name"] == "ANIDULAFUNGIN" for r in serum)
    # PPB routed to the free-fraction file
    assert any(r["endpoint_type"] == "PPB" for r in ppb)


def test_output_matches_drop_in_schema():
    path = _write_tmp(WEB_CSV)
    try:
        serum, _ = cx.convert([path])
    finally:
        os.remove(path)
    assert serum
    assert {"compound_name", "serum_active", "source_ref"} <= set(serum[0])
