#!/usr/bin/env python3
"""Convert a raw ChEMBL activity export into the Stage-0 drop-in files.

You (the user) download bioactivities for the echinocandins from ChEMBL in your
browser and drop the CSV(s) into stage0/data/ (see CHEMBL_PULL_INSTRUCTIONS.md).
This script reads whatever ChEMBL gave you — it tolerates both the website's
';'-delimited "human header" CSV and the API's ','-delimited lowercase CSV — and
emits:

  stage0/data/chembl_echinocandin_serum.csv   -> binary serum-active labels
        (columns: compound_name, serum_active, source_ref)  picked up automatically
        by build_binary_serum_activity.py
  stage0/data/chembl_free_fraction.csv        -> PPB / Fu records (Stage-1b seed)

Classification rules (the science stays here, not with the download):
  * serum context = the assay description mentions serum / plasma / albumin /
    protein binding.
  * a serum-context MIC below 100 ug/mL (and not ">") = serum-active (1);
    at/above the ceiling or ">" = inactive (0)  [the Stage-0 ceiling rule].
  * PPB / fraction-unbound rows are routed to the free-fraction file.

Run:
    python stage0/convert_chembl_export.py            # reads stage0/data/chembl_raw*.csv
    python stage0/convert_chembl_export.py FILE ...   # or explicit files
"""

from __future__ import annotations

import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

CEILING_UGML = 100.0
MIC_TYPES = {"MIC", "MIC50", "MIC80", "MIC90", "MIC90 ", "MIC100"}
PPB_TYPES = {"PPB", "FU", "FRACTION UNBOUND", "PLASMA PROTEIN BINDING",
             "PROTEIN BINDING", "FB", "PPB (HUMAN)", "PLASMA PROTEIN BINDING (HUMAN)"}
SERUM_WORDS = ("serum", "plasma", "albumin", "protein binding", "protein-binding")
# phrases that mean serum is ABSENT — must not count as serum context
SERUM_NEGATIONS = ("no serum", "without serum", "serum-free", "serum free",
                   "absence of serum", "no plasma", "protein-free", "protein free",
                   "serum-deficient", "minus serum")
# pharmacokinetic phrases: a plasma/serum DRUG-LEVEL reading, not a serum-medium MIC
PK_PHRASES = ("plasma concentration", "serum concentration", "plasma level",
              "serum level", "concentration in", "drug level", "drug concentration",
              "plasma exposure", "pharmacokinet", "plasma auc", "half-life")
UGML_UNITS = {"ug.ml-1", "ug/ml", "µg/ml", "ug ml-1", "mg/l", "mg.l-1"}

# ChEMBL uses different header spellings across web vs API exports; normalise.
HEADER_ALIASES = {
    "molecule_chembl_id": ["molecule chembl id", "molecule_chembl_id", "parent molecule chembl id"],
    "compound_name": ["molecule name", "molecule_pref_name", "compound name", "molecule pref name", "name"],
    "standard_type": ["standard type", "standard_type", "type"],
    "standard_relation": ["standard relation", "standard_relation", "relation"],
    "standard_value": ["standard value", "standard_value", "value"],
    "standard_units": ["standard units", "standard_units", "units"],
    "assay_description": ["assay description", "assay_description", "description"],
    "document_chembl_id": ["document chembl id", "document_chembl_id", "document"],
}


def _sniff_rows(path):
    """Read a ChEMBL CSV that may be ',' or ';' delimited; yield normalised dicts."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.readline()
        delim = ";" if sample.count(";") >= sample.count(",") else ","
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delim)
        lower = {(_norm(k)): k for k in (reader.fieldnames or [])}
        colmap = {}
        for canonical, aliases in HEADER_ALIASES.items():
            for a in aliases:
                if a in lower:
                    colmap[canonical] = lower[a]
                    break
        for row in reader:
            yield {c: (row.get(src, "") or "").strip() for c, src in colmap.items()}


def _norm(s):
    return (s or "").strip().lower().strip('"')


def _num(v):
    try:
        return float(str(v).replace(",", ""))
    except (TypeError, ValueError):
        return None


def _is_serum(desc):
    d = (desc or "").lower()
    if any(neg in d for neg in SERUM_NEGATIONS):
        return False        # e.g. "RPMI, no serum" is a serum-FREE assay
    if any(pk in d for pk in PK_PHRASES):
        return False        # e.g. "Plasma concentration in mouse ..." is a PK reading
    return any(w in d for w in SERUM_WORDS)


def convert(paths):
    serum_rows, ppb_rows = [], []
    seen_serum = set()
    for path in paths:
        for r in _sniff_rows(path):
            stype = (r.get("standard_type") or "").strip().upper()
            name = r.get("compound_name") or r.get("molecule_chembl_id") or ""
            src = r.get("document_chembl_id") or "ChEMBL"
            desc = r.get("assay_description", "")

            if stype in PPB_TYPES:
                ppb_rows.append({
                    "compound_name": name,
                    "molecule_chembl_id": r.get("molecule_chembl_id", ""),
                    "endpoint_type": stype,
                    "value": r.get("standard_value", ""),
                    "unit": r.get("standard_units", ""),
                    "source_ref": src,
                })
                continue

            if stype in MIC_TYPES and _is_serum(desc):
                val = _num(r.get("standard_value"))
                units = _norm(r.get("standard_units"))
                if val is None or val <= 0:
                    continue
                if units and units not in UGML_UNITS:
                    continue  # skip non ug/mL units to keep the ceiling comparable
                rel = r.get("standard_relation", "").strip("'\" ")
                active = 0 if (val >= CEILING_UGML or rel == ">") else 1
                key = (name.upper(), val, active, src)
                if key in seen_serum:
                    continue
                seen_serum.add(key)
                serum_rows.append({
                    "compound_name": name,
                    "compound_id": r.get("molecule_chembl_id", ""),
                    "serum_active": active,
                    "source_ref": f"{src}; {stype} {rel}{val} {r.get('standard_units','')}",
                })
    return serum_rows, ppb_rows


def _write(path, rows):
    if not rows:
        print(f"  (none) {os.path.basename(path)}")
        return
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows)} -> {path}")


def main(argv):
    paths = argv or sorted(glob.glob(os.path.join(DATA, "chembl_raw*.csv")))
    paths = [p for p in paths if os.path.exists(p)]
    if not paths:
        print("No raw ChEMBL export found. Put chembl_raw*.csv in stage0/data/ "
              "(see stage0/data/CHEMBL_PULL_INSTRUCTIONS.md).")
        return 1
    print(f"reading {len(paths)} file(s): {[os.path.basename(p) for p in paths]}")
    serum_rows, ppb_rows = convert(paths)
    os.makedirs(DATA, exist_ok=True)
    _write(os.path.join(DATA, "chembl_echinocandin_serum.csv"), serum_rows)
    _write(os.path.join(DATA, "chembl_free_fraction.csv"), ppb_rows)
    n_pos = sum(1 for r in serum_rows if r["serum_active"] == 1)
    print(f"serum-context MICs: {len(serum_rows)} "
          f"({n_pos} active, {len(serum_rows)-n_pos} inactive); "
          f"PPB/Fu rows: {len(ppb_rows)}")
    print("Now re-run: python stage0/build_binary_serum_activity.py")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
