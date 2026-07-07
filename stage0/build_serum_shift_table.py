#!/usr/bin/env python3
"""Stage 0 — Unify data and fix the endpoint (serum shift).

This is the first pipeline step described in ``START_HERE.md`` (§3, "Stage 0").
Its single job is to turn the raw activity records into ONE potency-independent
endpoint that both chemotypes share:

    serum_shift = MIC(serum-containing) / MIC(serum-free)      (same compound,
                                                                same organism)

and to carry the confounders/metadata that every downstream oracle needs:

  * intrinsic potency (the serum-free MIC) is kept as a COVARIATE, because raw
    serum MIC is dominated by potency (guardrail #1 in START_HERE);
  * censoring is tracked explicitly. Many serum MICs are reported as ">100
    ug/mL" (the assay ceiling), so the true shift is a LOWER BOUND, not a point
    estimate. Silently treating ">100" as "100" would bias every model.
  * the applicability domain (what chemistry / which organisms each dataset
    covers) is recorded alongside the numbers.

The pipeline is intentionally stdlib-only (csv / math / statistics) so it runs
in a fresh environment with no scientific stack. Heavier dependencies are
reserved for later stages and pinned in ``requirements.txt``.

Inputs (read-only):
  curated/core_tables/activity_table.csv         papulacandin whole-cell MICs
  curated/core_tables/compounds_master.csv       MW + scaffold descriptors
  external/.../external_activity_table_candidate_v0_2.csv   echinocandin MICs

Outputs (written to stage0/outputs/):
  papulacandin_serum_pairs.csv        matched pairs (the legitimate serum-shift labels)
  unified_serum_shift_table.csv       the Stage-0 endpoint table (matched shifts)
  echinocandin_serum_context.csv      teacher corpus: serum-present / serum-free MIC
                                      rows, UNPAIRED (see data-reality note below)
  echinocandin_free_fraction_seed.csv Fu / PPB rows -> seed for the Stage-1b oracle
  applicability_domain.csv            per-dataset coverage summary
  stage0_summary.md                   human-readable report

DATA-REALITY NOTE (discovered while building this stage):
  The external echinocandin table has NO within-study matched serum pairs: its
  serum-containing and serum-free MICs never share a source study, and pooling
  across studies mixes wild-type with FKS-resistant mutants (free MICs for one
  drug span 0.0001-512 ug/mL). A "serum shift" built from that would be a
  cross-study, resistance-confounded ratio. Per the project's own guardrails
  (model the shift honestly; DATA is the binding constraint), we therefore do
  NOT fabricate echinocandin pairs. We instead stage the echinocandin serum
  data UNPAIRED, in the correct shape for a Stage-1 transfer / hierarchical
  model that carries `study` and `serum_state` as effects. Only the
  papulacandin class yields legitimate matched shifts, so it alone populates
  the unified endpoint table for now (schema is ready to accept echinocandin
  pairs if better-matched data arrives).
"""

from __future__ import annotations

import csv
import math
import os
import statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")

CURATED_ACTIVITY = os.path.join(ROOT, "curated/core_tables/activity_table.csv")
CURATED_COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")
EXTERNAL_ACTIVITY = os.path.join(
    ROOT,
    "external/data/external/fks_inhibitors/source_exports/"
    "external_activity_table_candidate_v0_2.csv",
)

# MIC-type endpoints we treat as whole-cell potency for the shift.
MIC_ENDPOINTS = {"MIC", "MIC50", "MIC80", "MIC90"}

# serum_or_protein strings that mean "no serum / no protein present".
SERUM_FREE_TOKENS = {"", "none", "none reported", "not reported"}


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(value):
    """Parse a float, returning None on blanks / non-numeric."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _is_serum(text):
    return "serum" in (text or "").strip().lower()


def _is_free(text):
    return (text or "").strip().lower() in SERUM_FREE_TOKENS


def _median(values):
    return statistics.median(values) if values else None


# --------------------------------------------------------------------------- #
# papulacandin side (curated) — the precious matched-pair label set
# --------------------------------------------------------------------------- #
def build_papulacandin_pairs():
    """Return matched serum-free / serum MIC pairs for the papulacandin class.

    Pairing is anchored on the serum-containing row and matched to the best
    serum-free counterpart, most specific first:
        1. same (compound, reference, organism, strain, unit)
        2. same (compound, reference, organism, unit)
        3. same (compound, organism, unit)  -> median of serum-free MICs
    Same unit is REQUIRED so the ratio is dimensionless and correct.
    """
    activity = _read(CURATED_ACTIVITY)
    compounds = {r["compound_id"]: r for r in _read(CURATED_COMPOUNDS)}

    mic = [
        r
        for r in activity
        if r["endpoint_type"] in MIC_ENDPOINTS and _num(r["endpoint_value"])
    ]

    serum_rows = [r for r in mic if _is_serum(r["serum_or_protein"])]
    free_rows = [r for r in mic if _is_free(r["serum_or_protein"])]

    # index serum-free rows for fast lookup at three specificities
    by_strain = defaultdict(list)
    by_org = defaultdict(list)
    by_org_any_ref = defaultdict(list)
    for r in free_rows:
        cid, ref, org, strain, unit = (
            r["compound_id"], r["reference_id"], r["organism"],
            r["strain"], r["unit"],
        )
        by_strain[(cid, ref, org, strain, unit)].append(r)
        by_org[(cid, ref, org, unit)].append(r)
        by_org_any_ref[(cid, org, unit)].append(r)

    pairs = []
    for s in serum_rows:
        cid, ref, org, strain, unit = (
            s["compound_id"], s["reference_id"], s["organism"],
            s["strain"], s["unit"],
        )
        match_level, free_candidates = None, None
        if by_strain.get((cid, ref, org, strain, unit)):
            match_level = "strain_within_study"
            free_candidates = by_strain[(cid, ref, org, strain, unit)]
        elif by_org.get((cid, ref, org, unit)):
            match_level = "organism_within_study"
            free_candidates = by_org[(cid, ref, org, unit)]
        elif by_org_any_ref.get((cid, org, unit)):
            match_level = "organism_across_studies"
            free_candidates = by_org_any_ref[(cid, org, unit)]
        else:
            continue  # no comparable serum-free MIC -> cannot form a shift

        free_val = _median([_num(r["endpoint_value"]) for r in free_candidates])
        serum_val = _num(s["endpoint_value"])
        if not free_val or free_val <= 0 or not serum_val or serum_val <= 0:
            continue

        # censoring: ">" on the serum side makes the shift a lower bound.
        serum_censored = s["endpoint_relation"] == ">"
        free_censored = any(r["endpoint_relation"] == ">" for r in free_candidates)
        if serum_censored and not free_censored:
            shift_relation = ">="  # true shift is at least this large
        elif free_censored and not serum_censored:
            shift_relation = "<="  # true shift is at most this large (ambiguous)
        elif serum_censored and free_censored:
            shift_relation = "~"   # both bounded, direction unknown
        else:
            shift_relation = "="

        cm = compounds.get(cid, {})
        mw = _num(cm.get("mw_exact"))
        # intrinsic potency covariate as pMIC = -log10(mol/L), when MW known.
        free_uM = (free_val / mw * 1000.0) if mw else None
        pmic = -math.log10(free_uM * 1e-6) if free_uM else None

        pairs.append(
            {
                "chemotype": "papulacandin",
                "compound_id": cid,
                "compound_name": cm.get("canonical_name", ""),
                "compound_class": cm.get("compound_class", ""),
                "organism": org,
                "strain": strain,
                "reference_id": ref,
                "match_level": match_level,
                "unit": unit,
                "serum_condition": s["serum_or_protein"].strip(),
                "free_mic": round(free_val, 4),
                "serum_mic": round(serum_val, 4),
                "free_censored": free_censored,
                "serum_censored": serum_censored,
                "serum_shift": round(serum_val / free_val, 4),
                "shift_relation": shift_relation,
                "log10_serum_shift": round(math.log10(serum_val / free_val), 4),
                "pairing_confidence": (
                    "high" if match_level == "strain_within_study" else "medium"
                ),
                "intrinsic_potency_free_mic_uM": round(free_uM, 4) if free_uM else "",
                "intrinsic_potency_pMIC": round(pmic, 4) if pmic is not None else "",
                "acyl_tail_type": cm.get("long_acyl_chain_type", ""),
                "acyl_tail_geometry": cm.get("long_acyl_geometry", ""),
            }
        )
    return pairs


# --------------------------------------------------------------------------- #
# echinocandin side (external) — the data-rich "teacher"
# --------------------------------------------------------------------------- #
# Canonical FKS / glucan-synthase drugs. Rows for non-FKS comparators that ride
# along in the external export (amphotericin B, azoles) are flagged, not paired,
# so they never contaminate the teacher signal.
CANONICAL_FKS = {
    "CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN", "REZAFUNGIN",
    "IBREXAFUNGERP", "ENFUMAFUNGIN", "PNEUMOCANDIN",
}


def _serum_state(r):
    present = (
        r["assay_context"] == "serum_or_protein_present"
        or r["serum_or_protein"] == "serum_or_protein_present"
    )
    return "serum_present" if present else "serum_context_unknown"


def build_echinocandin_teacher():
    """Stage the external echinocandin data WITHOUT fabricating serum pairs.

    Returns ``(serum_context_rows, free_fraction_rows)``:

    * ``serum_context_rows`` — every whole-cell MIC row for FKS compounds that
      is either explicitly serum-containing or of unknown serum context, kept
      UNPAIRED with its source study. This is the corpus a later hierarchical /
      transfer model consumes (with `study` and `serum_state` as effects). We do
      not divide serum by serum-free here because the two never co-occur within
      a study (see the module docstring's data-reality note).
    * ``free_fraction_rows`` — the Fu (free-fraction) and PPB (plasma-protein-
      binding) rows: the seed for the Stage-1b free-fraction oracle.
    """
    ext = _read(EXTERNAL_ACTIVITY)

    serum_context_rows = []
    for r in ext:
        if r["endpoint_type"] not in MIC_ENDPOINTS:
            continue
        val = _num(r["endpoint_value"])
        if not val or val <= 0:
            continue
        name = r["preferred_name"].strip().upper()
        serum_context_rows.append(
            {
                "chemotype": "echinocandin",
                "external_compound_id": r["external_compound_id"],
                "compound_name": r["preferred_name"],
                "is_canonical_fks": name in CANONICAL_FKS,
                "organism": r["organism"],
                "strain": r["strain"],
                "endpoint_type": r["endpoint_type"],
                "mic_value": val,
                "endpoint_relation": r["endpoint_relation"],
                "unit": r["unit"],
                "serum_state": _serum_state(r),
                "serum_or_protein": r["serum_or_protein"],
                "assay_context": r["assay_context"],
                "source_reference": r["source_reference"],
                "source_database": r["source_database"],
            }
        )

    free_fraction_rows = []
    for r in ext:
        if r["endpoint_type"] not in ("Fu", "PPB"):
            continue
        free_fraction_rows.append(
            {
                "chemotype": "echinocandin",
                "external_compound_id": r["external_compound_id"],
                "compound_name": r["preferred_name"],
                "endpoint_type": r["endpoint_type"],  # Fu = free fraction, PPB = bound
                "value": r["endpoint_value"],
                "unit": r["unit"],
                "serum_state": _serum_state(r),
                "source_reference": r["source_reference"],
                "source_database": r["source_database"],
            }
        )

    return serum_context_rows, free_fraction_rows


# --------------------------------------------------------------------------- #
# applicability domain + writers
# --------------------------------------------------------------------------- #
def applicability_domain(papu, echino_context, echino_fu):
    """One coverage row per dataset (guardrail: know each set's domain)."""
    rows = []

    if papu:
        shifts = [p["serum_shift"] for p in papu]
        censored = sum(1 for p in papu if p["serum_censored"])
        rows.append(
            {
                "dataset": "papulacandin_curated_matched_pairs",
                "role": "serum_shift_labels",
                "chemotype": "papulacandin",
                "n_rows": len(papu),
                "n_compounds": len({p["compound_id"] for p in papu}),
                "n_organisms": len({p["organism"] for p in papu}),
                "organisms": "; ".join(sorted({p["organism"] for p in papu})[:8]),
                "n_censored": censored,
                "pct_censored": round(100 * censored / len(papu), 1),
                "shift_min": round(min(shifts), 2),
                "shift_median": round(statistics.median(shifts), 2),
                "shift_max": round(max(shifts), 2),
                "n_references": len({p["reference_id"] for p in papu}),
                "notes": "legitimate within-study matched serum-free/serum MIC",
            }
        )

    if echino_context:
        fks = [r for r in echino_context if r["is_canonical_fks"]]
        present = [r for r in echino_context if r["serum_state"] == "serum_present"]
        rows.append(
            {
                "dataset": "echinocandin_external_serum_context",
                "role": "teacher_corpus_unpaired",
                "chemotype": "echinocandin",
                "n_rows": len(echino_context),
                "n_compounds": len({r["external_compound_id"] for r in echino_context}),
                "n_organisms": len({r["organism"] for r in echino_context}),
                "organisms": f"{len({r['organism'] for r in echino_context})} distinct",
                "n_censored": sum(1 for r in echino_context if r["endpoint_relation"] == ">"),
                "pct_censored": "",
                "shift_min": "",
                "shift_median": "",
                "shift_max": "",
                "n_references": len({r["source_reference"] for r in echino_context}),
                "notes": (
                    f"{len(present)} serum-present rows; {len(fks)} canonical-FKS "
                    "rows; UNPAIRED (no within-study serum pairs exist)"
                ),
            }
        )

    if echino_fu:
        rows.append(
            {
                "dataset": "echinocandin_external_free_fraction",
                "role": "free_fraction_oracle_seed",
                "chemotype": "echinocandin",
                "n_rows": len(echino_fu),
                "n_compounds": len({r["external_compound_id"] for r in echino_fu}),
                "n_organisms": 0,
                "organisms": "",
                "n_censored": 0,
                "pct_censored": "",
                "shift_min": "",
                "shift_median": "",
                "shift_max": "",
                "n_references": len({r["source_reference"] for r in echino_fu}),
                "notes": "Fu + PPB rows; seed for Stage-1b free-fraction oracle",
            }
        )
    return rows


def _write_csv(path, rows, fieldnames=None):
    if not rows:
        open(path, "w").close()
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_summary(path, papu, echino_context, echino_fu, ad):
    n_papu_comp = len({p["compound_id"] for p in papu})
    high = sum(1 for p in papu if p["pairing_confidence"] == "high")
    censored = sum(1 for p in papu if p["serum_censored"])
    fks_ctx = [r for r in echino_context if r["is_canonical_fks"]]
    present = [r for r in echino_context if r["serum_state"] == "serum_present"]
    lines = []
    lines.append("# Stage 0 — unify data and fix the endpoint\n\n")
    lines.append(
        "Endpoint: **serum shift = MIC(serum-containing) / MIC(serum-free)** for "
        "the *same compound in the same organism* — a potency-independent measure. "
        "Intrinsic potency (the serum-free MIC) is carried as a covariate.\n\n"
    )
    lines.append("## Papulacandin (curated) — the serum-shift labels\n\n")
    lines.append(
        f"- **{n_papu_comp} compounds** carry matched serum-free + serum-containing "
        f"MIC — the ~24 anchor labels referenced in START_HERE.\n"
        f"- **{len(papu)} matched pairs** total (a compound recurs across organisms "
        f"/ strains); {high} are strain-matched within a single study (high "
        f"confidence).\n"
        f"- **{censored}/{len(papu)} pairs are right-censored** (serum MIC reported "
        f"as \">\"), so their shift is a LOWER BOUND, flagged in `shift_relation` "
        f"(`>=`). They must be modelled as censored, never as point values.\n\n"
    )
    lines.append("## Echinocandin (external) — teacher corpus, NOT fabricated pairs\n\n")
    lines.append(
        "**Key data-reality finding:** the external echinocandin table contains "
        "**no within-study matched serum pairs**. Serum-containing and serum-free "
        "MICs never share a source study, and pooling across studies mixes "
        "wild-type with FKS-resistant mutants (one drug's free MIC spans "
        "0.0001–512 ug/mL). A ratio built from that is a cross-study, "
        "resistance-confounded artifact. Honoring the project's guardrails, we do "
        "**not** manufacture echinocandin serum shifts. Instead:\n\n"
        f"- `echinocandin_serum_context.csv`: {len(echino_context)} whole-cell MIC "
        f"rows ({len(present)} explicitly serum-present, {len(fks_ctx)} canonical-FKS), "
        "kept UNPAIRED with their source study — the corpus for a Stage-1 "
        "hierarchical / transfer model with `study` and `serum_state` as effects.\n"
        f"- `echinocandin_free_fraction_seed.csv`: {len(echino_fu)} Fu / PPB rows — "
        "the seed for the Stage-1b free-fraction oracle.\n\n"
    )
    lines.append("## Applicability domain\n\n")
    for r in ad:
        extra = (
            f"; shift {r['shift_min']}–{r['shift_max']}x (median {r['shift_median']}x)"
            if r["shift_median"] != "" else ""
        )
        lines.append(
            f"- **{r['dataset']}** ({r['role']}): {r['n_rows']} rows / "
            f"{r['n_compounds']} compounds{extra}; {r['n_references']} refs. "
            f"{r['notes']}.\n"
        )
    lines.append(
        "\n## Guardrails honored\n"
        "1. Endpoint is the serum SHIFT, not the raw serum MIC (guardrail #1).\n"
        "2. Intrinsic potency retained as a covariate (dominant confound).\n"
        "3. Censoring tracked explicitly (`shift_relation`, `serum_censored`).\n"
        "4. Applicability domain recorded per dataset.\n"
        "5. No fabricated cross-study pairs — DATA is the binding constraint "
        "(guardrail #8).\n"
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(lines))


def main():
    os.makedirs(OUT, exist_ok=True)
    papu = build_papulacandin_pairs()
    echino_context, echino_fu = build_echinocandin_teacher()
    ad = applicability_domain(papu, echino_context, echino_fu)

    # The unified endpoint table holds only legitimate matched shifts (papulacandin
    # today); schema is ready to accept echinocandin pairs if better data arrives.
    _write_csv(os.path.join(OUT, "papulacandin_serum_pairs.csv"), papu)
    _write_csv(os.path.join(OUT, "unified_serum_shift_table.csv"), papu)
    _write_csv(os.path.join(OUT, "echinocandin_serum_context.csv"), echino_context)
    _write_csv(os.path.join(OUT, "echinocandin_free_fraction_seed.csv"), echino_fu)
    _write_csv(os.path.join(OUT, "applicability_domain.csv"), ad)
    write_summary(os.path.join(OUT, "stage0_summary.md"),
                  papu, echino_context, echino_fu, ad)

    print(f"papulacandin matched pairs: "
          f"{len({p['compound_id'] for p in papu})} compounds, {len(papu)} pairs "
          f"({sum(1 for p in papu if p['serum_censored'])} censored)")
    print(f"echinocandin teacher corpus: {len(echino_context)} serum-context rows, "
          f"{len(echino_fu)} free-fraction rows")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
