#!/usr/bin/env python3
"""
serum_gap_analysis.py
=====================
Foundational analysis for the Papulacandin / FKS1 project.

Core scientific question
------------------------
Papulacandin-type compounds are often potent in serum-FREE whole-cell assays
but lose activity in serum-CONTAINING assays (and in vivo). This script builds
the "serum gap" table: for every compound that has whole-cell MIC measurements
both with and without serum, it reports the serum-free MIC, the serum MIC, and
the serum shift (fold loss of activity caused by serum).

This is the dependent variable ("what we ultimately want to predict / design
against") for the whole project, so it is built first, transparently, and with
explicit handling of censored values (e.g. ">100", "<0.03").

Design notes
------------
* Censored values are NOT errors. ">100" means "MIC is at least 100" (inactive
  floor); "<0.03" means "MIC is at most 0.03" (potent ceiling). We keep the
  relation and flag the pair as censored so downstream modeling can treat it
  correctly (e.g. survival/Tobit-style handling) instead of silently dropping
  or mis-reading it.
* We restrict the matched-pair analysis to a single organism (Candida albicans)
  because that is where the serum vs serum-free pairs actually exist in this
  dataset. Mixing organisms would confound the serum effect with species
  differences.
* All MIC values for the serum-paired references are in ug/mL, so the serum
  SHIFT (a within-compound ratio) is unit-independent and robust. For
  cross-compound comparison we also report uM using the parent MW.

Outputs (written to analysis/outputs/)
--------------------------------------
* serum_gap_pairs.csv   : one row per compound with matched serum-free/serum MIC
* serum_gap_summary.txt : human-readable summary of what the data supports

Run:  python3 analysis/serum_gap_analysis.py
"""

from __future__ import annotations
import csv
import os
import statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def load_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def is_serum(serum_field: str) -> bool | None:
    """Classify the serum_or_protein field.

    Returns True (serum present), False (explicitly serum-free), or None
    (unknown / not reported -> excluded from matched pairs).
    """
    s = (serum_field or "").strip().lower()
    if s in ("none", ""):
        # 'none' = serum-free; '' is ambiguous but in this dataset the paired
        # references use 'none' explicitly, so treat blank as unknown.
        return False if s == "none" else None
    if "not reported" in s or "conditions not reported" in s:
        return None
    if any(k in s for k in ("serum", "protein", "fetal", "fbs", "fcs", "albumin")):
        return True
    if s == "none reported":
        return None
    return None


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fold(x):
    return f"{x:.1f}x" if x is not None else "n/a"


# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
compounds = {r["compound_id"]: r for r in load_csv(os.path.join(CORE, "compounds_master.csv"))}
activity = load_csv(os.path.join(CORE, "activity_table.csv"))

# Restrict to whole-cell MIC against Candida albicans (where serum pairs exist).
ALBICANS = {"candida albicans", "c. albicans"}
mic = [
    r for r in activity
    if r["endpoint_type"] == "MIC"
    and (r["organism"] or "").strip().lower() in ALBICANS
    and (r["unit"] or "").strip().lower() in ("ug/ml", "mcg/ml")
]

# --------------------------------------------------------------------------
# Aggregate per compound x serum-condition
# --------------------------------------------------------------------------
# bucket[compound_id][serum_bool] = list of (relation, value_ugmL)
bucket: dict[str, dict[bool, list]] = defaultdict(lambda: defaultdict(list))
for r in mic:
    s = is_serum(r["serum_or_protein"])
    if s is None:
        continue
    val = to_float(r["endpoint_value"])
    if val is None:
        continue
    bucket[r["compound_id"]][s].append((r["endpoint_relation"] or "=", val))


def summarize_condition(entries):
    """Collapse replicate MICs for one condition into a single representative.

    Returns (repr_value, relation, censored_bool, n). We use the median of the
    numeric values; if every measurement is censored in the same direction we
    carry the relation through and flag it.
    """
    if not entries:
        return None
    vals = [v for _, v in entries]
    rels = [rel for rel, _ in entries]
    med = statistics.median(vals)
    # Determine the relation that best represents the median measurement.
    # If all measurements share one censoring relation, keep it.
    uniq_rels = set(rels)
    if uniq_rels == {">"}:
        rel = ">"
    elif uniq_rels == {"<", "<="} or uniq_rels == {"<"} or uniq_rels == {"<="}:
        rel = "<"
    else:
        rel = "="
    censored = rel in (">", "<")
    return (med, rel, censored, len(entries))


# --------------------------------------------------------------------------
# Build the serum-gap table
# --------------------------------------------------------------------------
rows_out = []
for cid, conds in bucket.items():
    free = summarize_condition(conds.get(False, []))
    serum = summarize_condition(conds.get(True, []))
    if free is None or serum is None:
        continue  # only keep compounds with BOTH conditions = a matched pair

    free_v, free_rel, free_cens, free_n = free
    serum_v, serum_rel, serum_cens, serum_n = serum

    # Serum shift = serum MIC / serum-free MIC (>1 means serum hurts potency).
    shift = serum_v / free_v if free_v else None
    # Flag how trustworthy the shift number is given censoring.
    if free_cens or serum_cens:
        shift_note = "lower_bound" if serum_rel == ">" else "censored"
    else:
        shift_note = "exact"

    cm = compounds.get(cid, {})
    rows_out.append({
        "compound_id": cid,
        "name": cm.get("canonical_name", ""),
        "compound_class": cm.get("compound_class", ""),
        "parent_scaffold": cm.get("parent_scaffold", ""),
        "serumfree_mic_ugml": f"{free_rel}{free_v:g}",
        "serumfree_n": free_n,
        "serum_mic_ugml": f"{serum_rel}{serum_v:g}",
        "serum_n": serum_n,
        "serum_shift_fold": round(shift, 2) if shift else "",
        "shift_confidence": shift_note,
        "serum_active": "yes" if (not serum_cens and serum_v <= 50) else "no",
        "clogp": cm.get("clogp", ""),
        "long_acyl_chain_length": cm.get("long_acyl_chain_length", ""),
        "modification_summary": cm.get("modification_summary", ""),
    })

# Sort: serum-active first, then by smallest serum MIC.
def sort_key(r):
    active = 0 if r["serum_active"] == "yes" else 1
    sv = to_float(r["serum_mic_ugml"].lstrip("><=")) or 1e9
    return (active, sv)

rows_out.sort(key=sort_key)

# --------------------------------------------------------------------------
# Write CSV
# --------------------------------------------------------------------------
csv_path = os.path.join(OUT, "serum_gap_pairs.csv")
fieldnames = list(rows_out[0].keys()) if rows_out else []
with open(csv_path, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows_out)

# --------------------------------------------------------------------------
# Human-readable summary
# --------------------------------------------------------------------------
n_pairs = len(rows_out)
n_active = sum(1 for r in rows_out if r["serum_active"] == "yes")
n_total_serum_inactive = sum(1 for r in rows_out if r["serum_mic_ugml"].startswith(">"))

lines = []
lines.append("SERUM GAP ANALYSIS — Candida albicans whole-cell MIC")
lines.append("=" * 60)
lines.append(f"Compounds with matched serum-free AND serum MIC pairs : {n_pairs}")
lines.append(f"  of which retain measurable serum activity (MIC<=50): {n_active}")
lines.append(f"  of which are fully serum-inactivated (serum MIC >X) : {n_total_serum_inactive}")
lines.append("")
lines.append("Interpretation")
lines.append("-" * 60)
lines.append(
    "The matched-pair signal lives almost entirely in two reference series\n"
    "(Yeung 1996 fusacandin analogs; Jackson 1995). For most compounds the\n"
    "serum-free MIC is potent but the serum MIC is censored at '>100' — i.e.\n"
    "serum abolishes activity. A small minority retain partial serum activity\n"
    "(MIC 12.5-50 ug/mL). These partial-retainers are the SAR anchors for the\n"
    "whole project: they are the only positive examples of 'serum-tolerant'\n"
    "chemistry we currently have.")
lines.append("")
lines.append("Top serum-tolerant compounds (the leads to learn from):")
for r in rows_out[:8]:
    lines.append(
        f"  {r['compound_id']}  {r['name'][:38]:38s}  "
        f"free {r['serumfree_mic_ugml']:>7s} -> serum {r['serum_mic_ugml']:>7s}  "
        f"shift {r['serum_shift_fold']} ({r['shift_confidence']})")
lines.append("")
lines.append(f"Full table: {os.path.relpath(csv_path, ROOT)}")

summary = "\n".join(lines)
with open(os.path.join(OUT, "serum_gap_summary.txt"), "w", encoding="utf-8") as fh:
    fh.write(summary + "\n")

print(summary)
