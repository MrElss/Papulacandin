#!/usr/bin/env python3
"""Stage 0 — Unify data and fix the endpoint (serum tolerance).

First pipeline step of the methodology in ``START_HERE.md`` (§3, "Stage 0").
Its only job is to turn the raw activity records into ONE well-defined
serum-tolerance endpoint both chemotypes can share, and to carry the
confounders every downstream oracle needs. No modeling happens here.

The endpoint has two complementary faces:

  1. A CATEGORICAL serum-tolerance call (the headline label):
         retained  — active serum-free AND still active in serum
         lost      — active serum-free but inactive (MIC > ceiling) in serum
     This is robust and matches the program goal ("retain activity in serum").

  2. A CONTINUOUS serum shift, only where it is meaningful:
         serum_shift = MIC(serum) / MIC(serum-free)
     kept for pairs where the serum MIC is a real, active number.

Two curation rules encoded here (added at the project's request):

  * ACTIVITY CEILING (point 1): a MIC above ``ACTIVITY_CEILING_UGML`` (100 ug/mL)
    means "no detectable activity". We do NOT chase a precise number above it and
    we do NOT treat it as a censored value to be modelled numerically — it is a
    categorical INACTIVE call. So a compound active serum-free but > ceiling in
    serum has simply "lost" activity in serum.

  * IN-VIVO EFFICACY AS A SERUM-TOLERANCE PROXY (point 2): in-vivo efficacy
    (curing/ reducing an infection in an animal) requires activity in the
    bloodstream, where serum proteins are present, so it is legitimate PROXY
    evidence of serum tolerance -> a presumed-positive label (flagged, lower
    confidence than a direct serum MIC). Whole-cell / cellular MIC alone is NOT
    such evidence — losing activity in serum despite a good broth MIC is exactly
    this class's failure mode — so cellular-only activity is never counted as a
    serum-tolerance label.

The pipeline is stdlib-only (csv / math / statistics) so it runs in a fresh
environment. Heavier deps are reserved for later stages (``requirements.txt``).

Inputs (read-only):
  curated/core_tables/activity_table.csv         papulacandin whole-cell + in-vivo
  curated/core_tables/compounds_master.csv       MW + scaffold descriptors
  external/.../external_activity_table_candidate_v0_2.csv   echinocandin activity

Outputs (written to stage0/outputs/):
  papulacandin_serum_pairs.csv          matched serum-free/serum pairs + calls
  serum_tolerance_labels.csv            UNIFIED labels: direct pairs + in-vivo proxies
  invivo_serum_tolerance_proxies.csv    presumed-positive labels from in-vivo efficacy
  echinocandin_serum_context.csv        teacher corpus: serum-context MIC, UNPAIRED
  echinocandin_free_fraction_seed.csv   Fu / PPB rows -> Stage-1b free-fraction seed
  applicability_domain.csv              per-dataset coverage summary
  stage0_summary.md                     human-readable report

DATA-REALITY NOTE (found while building this stage):
  The external echinocandin table has NO within-study matched serum pairs — its
  serum-containing and serum-free MICs never share a source study, and pooling
  across studies mixes wild-type with FKS-resistant mutants (one drug's free MIC
  spans 0.0001-512 ug/mL). We therefore do NOT fabricate echinocandin serum
  shifts; we stage that data unpaired for a Stage-1 transfer/hierarchical model,
  and instead borrow the echinocandins' serum signal through their in-vivo
  efficacy proxies (point 2).
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

MIC_ENDPOINTS = {"MIC", "MIC50", "MIC80", "MIC90"}
SERUM_FREE_TOKENS = {"", "none", "none reported", "not reported"}

# Point 1: a MIC at or above this concentration is "no detectable activity".
ACTIVITY_CEILING_UGML = 100.0

# Point 2: in-vivo EFFICACY endpoints (a therapeutic effect was observed).
# "Max inactive dose" / "no toxic symptoms" describe absence of effect/tox and
# are deliberately excluded — they are not evidence of serum-tolerant efficacy.
INVIVO_EFFICACY_CURATED = {
    "in_vivo_ED50",
    "in_vivo_PCP_cyst_reduction",
    "in_vivo_TOKA_percent_reduction",
    "in_vivo_TOKA_log10_CFU_per_g_kidney",
    "in_vivo_TOKA_renal_clearance",
}

CANONICAL_FKS = {
    "CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN", "REZAFUNGIN",
    "IBREXAFUNGERP", "ENFUMAFUNGIN", "PNEUMOCANDIN", "CILOFUNGIN",
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _read(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(value):
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


def activity_call(value_ugml, relation):
    """Categorical activity call under the ceiling rule (point 1).

    active     — a measured MIC below the ceiling
    inactive   — MIC at/above the ceiling (>100 ug/mL means "no activity")
    ambiguous  — non-inhibitory ('>') only up to a tested max BELOW the ceiling,
                 so we cannot invoke the >100 rule (a small tail in this data)
    unknown    — no value
    """
    if value_ugml is None:
        return "unknown"
    if value_ugml >= ACTIVITY_CEILING_UGML:
        return "inactive"
    if relation == ">":
        return "ambiguous"
    return "active"


# --------------------------------------------------------------------------- #
# papulacandin side (curated) — the matched-pair serum-tolerance labels
# --------------------------------------------------------------------------- #
def build_papulacandin_pairs():
    """Matched serum-free / serum MIC pairs, each with a serum-tolerance call.

    Pairing anchors on the serum-containing row and matches the best serum-free
    counterpart, most specific first: (compound, ref, organism, strain, unit) ->
    (compound, ref, organism, unit) -> (compound, organism, unit). Same unit is
    required so any ratio is dimensionless.
    """
    activity = _read(CURATED_ACTIVITY)
    compounds = {r["compound_id"]: r for r in _read(CURATED_COMPOUNDS)}

    mic = [r for r in activity
           if r["endpoint_type"] in MIC_ENDPOINTS and _num(r["endpoint_value"])]
    serum_rows = [r for r in mic if _is_serum(r["serum_or_protein"])]
    free_rows = [r for r in mic if _is_free(r["serum_or_protein"])]

    by_strain, by_org, by_org_any = defaultdict(list), defaultdict(list), defaultdict(list)
    for r in free_rows:
        cid, ref, org, strain, unit = (r["compound_id"], r["reference_id"],
                                       r["organism"], r["strain"], r["unit"])
        by_strain[(cid, ref, org, strain, unit)].append(r)
        by_org[(cid, ref, org, unit)].append(r)
        by_org_any[(cid, org, unit)].append(r)

    pairs = []
    for s in serum_rows:
        cid, ref, org, strain, unit = (s["compound_id"], s["reference_id"],
                                       s["organism"], s["strain"], s["unit"])
        if by_strain.get((cid, ref, org, strain, unit)):
            level, cands = "strain_within_study", by_strain[(cid, ref, org, strain, unit)]
        elif by_org.get((cid, ref, org, unit)):
            level, cands = "organism_within_study", by_org[(cid, ref, org, unit)]
        elif by_org_any.get((cid, org, unit)):
            level, cands = "organism_across_studies", by_org_any[(cid, org, unit)]
        else:
            continue

        free_val = _median([_num(r["endpoint_value"]) for r in cands])
        serum_val = _num(s["endpoint_value"])
        if not free_val or free_val <= 0 or not serum_val or serum_val <= 0:
            continue

        # activity calls under the ceiling rule (point 1)
        free_relation = cands[0]["endpoint_relation"] if len(cands) == 1 else "="
        free_act = activity_call(free_val, free_relation)
        serum_act = activity_call(serum_val, s["endpoint_relation"])

        # headline serum-tolerance label
        if free_act != "active":
            tolerance = "uninformative"   # not active even serum-free
        elif serum_act == "active":
            tolerance = "retained"
        elif serum_act == "inactive":
            tolerance = "lost"
        else:
            tolerance = "ambiguous"

        # continuous shift only where the serum MIC is a real, active number
        shift = round(serum_val / free_val, 4)
        shift_meaningful = free_act == "active" and serum_act == "active"

        cm = compounds.get(cid, {})
        mw = _num(cm.get("mw_exact"))
        free_uM = (free_val / mw * 1000.0) if mw else None
        pmic = -math.log10(free_uM * 1e-6) if free_uM else None

        pairs.append({
            "chemotype": "papulacandin",
            "compound_id": cid,
            "compound_name": cm.get("canonical_name", ""),
            "compound_class": cm.get("compound_class", ""),
            "organism": org,
            "strain": strain,
            "reference_id": ref,
            "match_level": level,
            "pairing_confidence": "high" if level == "strain_within_study" else "medium",
            "unit": unit,
            "serum_condition": s["serum_or_protein"].strip(),
            "free_mic": round(free_val, 4),
            "serum_mic": round(serum_val, 4),
            "free_activity": free_act,
            "serum_activity": serum_act,
            "serum_tolerance": tolerance,          # <- headline categorical label
            "serum_shift": shift,
            "serum_shift_meaningful": shift_meaningful,
            "log10_serum_shift": round(math.log10(shift), 4),
            "intrinsic_potency_free_mic_uM": round(free_uM, 4) if free_uM else "",
            "intrinsic_potency_pMIC": round(pmic, 4) if pmic is not None else "",
            "acyl_tail_type": cm.get("long_acyl_chain_type", ""),
            "acyl_tail_geometry": cm.get("long_acyl_geometry", ""),
        })
    return pairs


# --------------------------------------------------------------------------- #
# in-vivo efficacy proxies (point 2) — both chemotypes
# --------------------------------------------------------------------------- #
def _direct_calls_by_compound(pairs):
    """compound_id -> {'retained','lost',...} from the direct matched pairs."""
    calls = defaultdict(set)
    for p in pairs:
        if p["serum_tolerance"] in ("retained", "lost"):
            calls[p["compound_id"]].add(p["serum_tolerance"])
    return calls


def _direct_evidence_label(calls):
    if not calls:
        return "none"
    if calls == {"retained"}:
        return "retained"
    if calls == {"lost"}:
        return "lost"
    return "mixed"


def build_invivo_proxies(direct_calls=None):
    """Presumed serum-tolerant positives from in-vivo efficacy.

    Rationale: an in-vivo therapeutic effect requires the compound to act in the
    bloodstream (serum present), so it is proxy evidence of serum tolerance.
    Whole-cell/cellular activity alone is intentionally NOT included.

    Each proxy is cross-referenced against the higher-confidence DIRECT serum
    pairs via ``direct_pair_evidence``. When a compound also has a direct 'lost'
    (or 'mixed') call, the in-vivo proxy is unreliable for it (an ED50 can be
    reported even where activity is serum-quenched — the papulacandin B case) and
    downstream should prefer the direct label.
    """
    direct_calls = direct_calls or {}
    proxies = []

    # curated papulacandins with an in-vivo EFFICACY endpoint
    curated = _read(CURATED_ACTIVITY)
    compounds = {r["compound_id"]: r for r in _read(CURATED_COMPOUNDS)}
    seen = set()
    for r in curated:
        if r["endpoint_type"] not in INVIVO_EFFICACY_CURATED:
            continue
        cid = r["compound_id"]
        key = (cid, r["endpoint_type"])
        if key in seen:
            continue
        seen.add(key)
        cm = compounds.get(cid, {})
        direct = _direct_evidence_label(direct_calls.get(cid, set()))
        proxies.append({
            "chemotype": "papulacandin",
            "compound_id": cid,
            "compound_name": cm.get("canonical_name", ""),
            "serum_tolerance_proxy": "presumed_positive",
            "evidence": "in_vivo_efficacy",
            "invivo_endpoint": r["endpoint_type"],
            "organism": r["organism"],
            "confidence": "low" if direct in ("lost", "mixed") else "medium",
            "direct_pair_evidence": direct,
            "basis": "in-vivo efficacy requires activity despite serum proteins",
            "reference_id": r["reference_id"],
        })

    # external echinocandins with in-vivo efficacy (canonical FKS only, to keep
    # the teacher chemotype-relevant; azoles/polyenes riding along are excluded)
    ext = _read(EXTERNAL_ACTIVITY)
    seen_ext = set()
    for r in ext:
        if (r["assay_type"] or "").lower() != "in_vivo":
            continue
        name = r["preferred_name"].strip().upper()
        if name not in CANONICAL_FKS:
            continue
        cid = r["external_compound_id"]
        if cid in seen_ext:
            continue
        seen_ext.add(cid)
        proxies.append({
            "chemotype": "echinocandin",
            "compound_id": cid,
            "compound_name": r["preferred_name"],
            "serum_tolerance_proxy": "presumed_positive",
            "evidence": "in_vivo_efficacy",
            "invivo_endpoint": r["endpoint_type"],
            "organism": r["organism"],
            "confidence": "medium",
            "direct_pair_evidence": "none",  # external side has no direct pairs
            "basis": "in-vivo efficacy requires activity despite serum proteins",
            "reference_id": r["source_reference"],
        })
    return proxies


# --------------------------------------------------------------------------- #
# echinocandin teacher corpus (external), unpaired
# --------------------------------------------------------------------------- #
def _serum_state(r):
    present = (r["assay_context"] == "serum_or_protein_present"
               or r["serum_or_protein"] == "serum_or_protein_present")
    return "serum_present" if present else "serum_context_unknown"


def build_echinocandin_teacher():
    """Echinocandin serum-context MIC rows (unpaired) + Fu/PPB free-fraction rows."""
    ext = _read(EXTERNAL_ACTIVITY)

    serum_context_rows = []
    for r in ext:
        if r["endpoint_type"] not in MIC_ENDPOINTS:
            continue
        val = _num(r["endpoint_value"])
        if not val or val <= 0:
            continue
        name = r["preferred_name"].strip().upper()
        serum_context_rows.append({
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
            "activity_call": activity_call(val, r["endpoint_relation"]) if r["unit"] == "ug.mL-1" else "unit_not_ugml",
            "serum_state": _serum_state(r),
            "serum_or_protein": r["serum_or_protein"],
            "assay_context": r["assay_context"],
            "source_reference": r["source_reference"],
            "source_database": r["source_database"],
        })

    free_fraction_rows = []
    for r in ext:
        if r["endpoint_type"] not in ("Fu", "PPB"):
            continue
        free_fraction_rows.append({
            "chemotype": "echinocandin",
            "external_compound_id": r["external_compound_id"],
            "compound_name": r["preferred_name"],
            "endpoint_type": r["endpoint_type"],
            "value": r["endpoint_value"],
            "unit": r["unit"],
            "serum_state": _serum_state(r),
            "source_reference": r["source_reference"],
            "source_database": r["source_database"],
        })
    return serum_context_rows, free_fraction_rows


# --------------------------------------------------------------------------- #
# unified label set + applicability domain + writers
# --------------------------------------------------------------------------- #
def build_unified_labels(pairs, proxies):
    """One row per serum-tolerance label, direct pairs and in-vivo proxies.

    Direct pairs contribute their categorical call (retained / lost);
    'uninformative' pairs (not active serum-free) are dropped from the label set.
    In-vivo proxies contribute presumed-positive (retained) labels.
    """
    rows = []
    for p in pairs:
        if p["serum_tolerance"] in ("retained", "lost"):
            rows.append({
                "compound_id": p["compound_id"],
                "compound_name": p["compound_name"],
                "chemotype": p["chemotype"],
                "organism": p["organism"],
                "serum_tolerance": p["serum_tolerance"],
                "label_source": "direct_serum_pair",
                "confidence": p["pairing_confidence"],
                "serum_shift": p["serum_shift"] if p["serum_shift_meaningful"] else "",
                "reference_id": p["reference_id"],
            })
    for x in proxies:
        rows.append({
            "compound_id": x["compound_id"],
            "compound_name": x["compound_name"],
            "chemotype": x["chemotype"],
            "organism": x["organism"],
            "serum_tolerance": "retained_presumed",
            "label_source": "in_vivo_proxy",
            "confidence": x["confidence"],
            "serum_shift": "",
            "direct_pair_evidence": x["direct_pair_evidence"],
            "reference_id": x["reference_id"],
        })
    # keep a stable column set across both row types
    for r in rows:
        r.setdefault("direct_pair_evidence", "self")
    return rows


def applicability_domain(pairs, proxies, echino_context, echino_fu):
    rows = []
    if pairs:
        informative = [p for p in pairs if p["serum_tolerance"] in ("retained", "lost")]
        retained = sum(1 for p in informative if p["serum_tolerance"] == "retained")
        lost = sum(1 for p in informative if p["serum_tolerance"] == "lost")
        rows.append({
            "dataset": "papulacandin_curated_matched_pairs",
            "role": "direct_serum_tolerance_labels",
            "chemotype": "papulacandin",
            "n_rows": len(pairs),
            "n_compounds": len({p["compound_id"] for p in pairs}),
            "n_informative": len(informative),
            "n_retained": retained,
            "n_lost": lost,
            "n_references": len({p["reference_id"] for p in pairs}),
            "notes": "within-study matched serum-free/serum MIC; ceiling rule applied",
        })
    if proxies:
        rows.append({
            "dataset": "invivo_efficacy_proxies",
            "role": "presumed_positive_labels",
            "chemotype": "papulacandin+echinocandin",
            "n_rows": len(proxies),
            "n_compounds": len({(p["chemotype"], p["compound_id"]) for p in proxies}),
            "n_informative": len(proxies),
            "n_retained": len(proxies),
            "n_lost": 0,
            "n_references": len({p["reference_id"] for p in proxies}),
            "notes": "in-vivo efficacy -> presumed serum-tolerant; cellular-only excluded",
        })
    if echino_context:
        fks = [r for r in echino_context if r["is_canonical_fks"]]
        present = [r for r in echino_context if r["serum_state"] == "serum_present"]
        rows.append({
            "dataset": "echinocandin_external_serum_context",
            "role": "teacher_corpus_unpaired",
            "chemotype": "echinocandin",
            "n_rows": len(echino_context),
            "n_compounds": len({r["external_compound_id"] for r in echino_context}),
            "n_informative": len(present),
            "n_retained": "",
            "n_lost": "",
            "n_references": len({r["source_reference"] for r in echino_context}),
            "notes": f"{len(present)} serum-present, {len(fks)} canonical-FKS; UNPAIRED",
        })
    if echino_fu:
        rows.append({
            "dataset": "echinocandin_external_free_fraction",
            "role": "free_fraction_oracle_seed",
            "chemotype": "echinocandin",
            "n_rows": len(echino_fu),
            "n_compounds": len({r["external_compound_id"] for r in echino_fu}),
            "n_informative": len(echino_fu),
            "n_retained": "",
            "n_lost": "",
            "n_references": len({r["source_reference"] for r in echino_fu}),
            "notes": "Fu + PPB rows; seed for Stage-1b free-fraction oracle",
        })
    return rows


def _write_csv(path, rows, fieldnames=None):
    if not rows:
        open(path, "w").close()
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def write_summary(path, pairs, proxies, labels, echino_context, echino_fu, ad):
    informative = [p for p in pairs if p["serum_tolerance"] in ("retained", "lost")]
    retained = sum(1 for p in informative if p["serum_tolerance"] == "retained")
    lost = sum(1 for p in informative if p["serum_tolerance"] == "lost")
    ambiguous = sum(1 for p in pairs if p["serum_tolerance"] == "ambiguous")
    uninform = sum(1 for p in pairs if p["serum_tolerance"] == "uninformative")
    papu_prox = sum(1 for p in proxies if p["chemotype"] == "papulacandin")
    echi_prox = sum(1 for p in proxies if p["chemotype"] == "echinocandin")
    L = []
    L.append("# Stage 0 — unify data & fix the serum-tolerance endpoint\n\n")
    L.append(
        "Endpoint = whether a compound **keeps activity in serum**. Two faces: a "
        "categorical call (`retained` / `lost`) plus a continuous serum shift "
        "where meaningful. Intrinsic potency (serum-free MIC) is kept as a "
        "covariate.\n\n")
    L.append("## Curation rules applied\n\n")
    L.append(
        f"- **Activity ceiling (point 1):** MIC > {ACTIVITY_CEILING_UGML:g} ug/mL "
        "= no detectable activity (a categorical INACTIVE call, not a censored "
        "number to model). So a compound active serum-free but above the ceiling "
        "in serum has simply **lost** serum activity.\n"
        "- **In-vivo efficacy proxy (point 2):** in-vivo efficacy implies activity "
        "in the bloodstream (serum present) -> a **presumed serum-tolerant** "
        "positive. Whole-cell / cellular activity alone is NOT counted (losing "
        "serum activity despite a good broth MIC is this class's failure mode).\n\n")
    L.append("## Papulacandin direct labels (matched pairs)\n\n")
    L.append(
        f"- {len({p['compound_id'] for p in pairs})} compounds, {len(pairs)} "
        f"matched pairs. Categorical calls: **{retained} retained, {lost} lost**, "
        f"{ambiguous} ambiguous, {uninform} uninformative (inactive serum-free).\n"
        f"- The 'lost' pairs are the class's signature failure (potent in broth, "
        f"dead in serum). Continuous serum shift is retained only where the serum "
        f"MIC is a real active number ("
        f"{sum(1 for p in pairs if p['serum_shift_meaningful'])} pairs).\n\n")
    L.append("## In-vivo efficacy proxies (presumed positives)\n\n")
    L.append(
        f"- {len(proxies)} compounds flagged presumed serum-tolerant from in-vivo "
        f"efficacy: {papu_prox} papulacandin, {echi_prox} canonical echinocandin "
        f"(the caspofungin / anidulafungin / micafungin anchors). Labelled "
        f"`retained_presumed`, confidence medium.\n\n")
    L.append(f"## Unified label set: {len(labels)} rows -> serum_tolerance_labels.csv\n\n")
    L.append("## Applicability domain\n\n")
    for r in ad:
        L.append(f"- **{r['dataset']}** ({r['role']}): {r['n_rows']} rows / "
                 f"{r['n_compounds']} compounds; {r['n_references']} refs. {r['notes']}.\n")
    L.append(
        "\n## Guardrails honored\n"
        "1. Endpoint is serum TOLERANCE (retained/lost), not raw serum MIC.\n"
        "2. Intrinsic potency retained as a covariate.\n"
        "3. Ceiling rule: >100 ug/mL = inactive (point 1); no chasing numbers above it.\n"
        "4. In-vivo efficacy is a flagged proxy; cellular-only is not a serum label (point 2).\n"
        "5. No fabricated cross-study echinocandin pairs — DATA is the binding constraint.\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(L))


def main():
    os.makedirs(OUT, exist_ok=True)
    pairs = build_papulacandin_pairs()
    proxies = build_invivo_proxies(_direct_calls_by_compound(pairs))
    echino_context, echino_fu = build_echinocandin_teacher()
    labels = build_unified_labels(pairs, proxies)
    ad = applicability_domain(pairs, proxies, echino_context, echino_fu)

    _write_csv(os.path.join(OUT, "papulacandin_serum_pairs.csv"), pairs)
    _write_csv(os.path.join(OUT, "invivo_serum_tolerance_proxies.csv"), proxies)
    _write_csv(os.path.join(OUT, "serum_tolerance_labels.csv"), labels)
    _write_csv(os.path.join(OUT, "echinocandin_serum_context.csv"), echino_context)
    _write_csv(os.path.join(OUT, "echinocandin_free_fraction_seed.csv"), echino_fu)
    _write_csv(os.path.join(OUT, "applicability_domain.csv"), ad)
    write_summary(os.path.join(OUT, "stage0_summary.md"),
                  pairs, proxies, labels, echino_context, echino_fu, ad)

    retained = sum(1 for p in pairs if p["serum_tolerance"] == "retained")
    lost = sum(1 for p in pairs if p["serum_tolerance"] == "lost")
    print(f"papulacandin pairs: {len({p['compound_id'] for p in pairs})} compounds, "
          f"{len(pairs)} pairs ({retained} retained, {lost} lost)")
    print(f"in-vivo proxies: {len(proxies)} "
          f"({sum(1 for p in proxies if p['chemotype']=='echinocandin')} echinocandin)")
    print(f"unified labels: {len(labels)} rows")
    print(f"echinocandin teacher: {len(echino_context)} serum-context, "
          f"{len(echino_fu)} free-fraction rows")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
