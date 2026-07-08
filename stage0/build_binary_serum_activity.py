#!/usr/bin/env python3
"""Stage 0 (binary reframe) — pooled serum-ACTIVE / INACTIVE label set.

Motivated by the project decision (and confirmed by the literature — see
planning/serum_mechanism_evidence.md) that serum inactivation is not reducible to
free fraction + stability: a direct serum->target effect exists that no descriptor
predicts. The robust move is therefore to *measure/collect the serum outcome
directly* and, for discovery, to represent it as a coarse BINARY label —
"retains activity in serum?" — which lets us pool far more, and more diverse,
data than the quantitative serum-shift endpoint could.

Endpoint: ``serum_active`` in {1, 0}. A MIC below the ceiling in a serum-containing
assay = active; at/above the ceiling (Stage-0 rule: >100 ug/mL = no activity) =
inactive. In-vivo efficacy and clinical use are additional (presumed-positive)
evidence, since serum in-vitro activity tracks animal-model outcomes (Nasar 2013).

Pooled sources
  1. papulacandin serum-containing MIC        (curated)          direct, high
  2. echinocandin serum-present MIC           (external)         direct, high
  3. in-vivo efficacy proxies                 (both chemotypes)  presumed, medium
  4. approved-echinocandin clinical annotation                    presumed, high
  5. ChEMBL echinocandin serum records        (optional drop-in) direct, high
     -> place a CSV at stage0/data/chembl_echinocandin_serum.csv with columns
        compound_name, serum_active (1/0), source_ref; ingested automatically
        when present. (Pending the ChEMBL connector; the pipeline runs without it.)

Outputs
  binary_serum_activity_observations.csv   one row per observation
  binary_serum_activity_labels.csv         one row per compound (consensus)
  binary_serum_activity_summary.md         coverage + the data-expansion win
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
DATA = os.path.join(HERE, "data")

CURATED_ACTIVITY = os.path.join(ROOT, "curated/core_tables/activity_table.csv")
CURATED_COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")
EXTERNAL_ACTIVITY = os.path.join(
    ROOT, "external/data/external/fks_inhibitors/source_exports/"
    "external_activity_table_candidate_v0_2.csv")
INVIVO_PROXIES = os.path.join(OUT, "invivo_serum_tolerance_proxies.csv")
CHEMBL_DROP_IN = os.path.join(DATA, "chembl_echinocandin_serum.csv")

CEILING_UGML = 100.0
APPROVED_ECHINOCANDINS = {
    "CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN", "REZAFUNGIN",
}


def _read(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _active_call(value, relation):
    """Binary serum-activity from a MIC (ceiling rule). None if unusable."""
    if value is None or value <= 0:
        return None
    if value >= CEILING_UGML or relation == ">":
        return 0
    return 1


# --------------------------------------------------------------------------- #
# observation collectors (one dict per observation)
# --------------------------------------------------------------------------- #
def _obs(chemotype, cid, name, active, source, confidence, detail):
    return {"chemotype": chemotype, "compound_id": cid, "compound_name": name,
            "serum_active": active, "label_source": source,
            "confidence": confidence, "detail": detail}


def papulacandin_serum_mic():
    act = _read(CURATED_ACTIVITY)
    cm = {r["compound_id"]: r for r in _read(CURATED_COMPOUNDS)}
    out = []
    for r in act:
        if r["endpoint_type"] not in ("MIC", "MIC50"):
            continue
        if "serum" not in (r["serum_or_protein"] or "").lower():
            continue
        val = _num(r["endpoint_value"])
        call = _active_call(val, r["endpoint_relation"])
        if call is None:
            continue
        name = cm.get(r["compound_id"], {}).get("canonical_name", "")
        out.append(_obs("papulacandin", r["compound_id"], name, call,
                        "papulacandin_serum_mic", "high",
                        f"{r['serum_or_protein'].strip()}; {r['organism']}; "
                        f"MIC {r['endpoint_relation']}{val} ug/mL"))
    return out


def echinocandin_serum_mic():
    ext = _read(EXTERNAL_ACTIVITY)
    out = []
    for r in ext:
        if r["endpoint_type"] not in ("MIC", "MIC50", "MIC80", "MIC90"):
            continue
        present = (r["assay_context"] == "serum_or_protein_present"
                   or r["serum_or_protein"] == "serum_or_protein_present")
        if not present or r["unit"] != "ug.mL-1":
            continue
        val = _num(r["endpoint_value"])
        call = _active_call(val, r["endpoint_relation"])
        if call is None:
            continue
        out.append(_obs("echinocandin", r["external_compound_id"],
                        r["preferred_name"], call, "echinocandin_serum_mic",
                        "high", f"serum-present; {r['organism']}; "
                        f"MIC {r['endpoint_relation']}{val} ug/mL"))
    return out


def invivo_proxy_positives():
    if not os.path.exists(INVIVO_PROXIES):
        return []
    out = []
    for r in _read(INVIVO_PROXIES):
        out.append(_obs(r["chemotype"], r["compound_id"], r["compound_name"],
                        1, "in_vivo_proxy", "medium",
                        f"{r['invivo_endpoint']}; {r.get('organism', '')}"))
    return out


def clinical_approved_positives():
    """Approved echinocandins are serum-active by definition (used in bloodstream)."""
    ext = _read(EXTERNAL_ACTIVITY)
    seen, out = set(), []
    for r in ext:
        name = (r["preferred_name"] or "").strip().upper()
        if name in APPROVED_ECHINOCANDINS and name not in seen:
            seen.add(name)
            out.append(_obs("echinocandin", r["external_compound_id"],
                            r["preferred_name"], 1, "clinical_approved", "high",
                            "approved echinocandin; active in the bloodstream"))
    return out


def chembl_drop_in():
    """Ingest ChEMBL echinocandin serum records if the drop-in CSV is present."""
    if not os.path.exists(CHEMBL_DROP_IN):
        return [], False
    out = []
    for r in _read(CHEMBL_DROP_IN):
        call = _num(r.get("serum_active"))
        if call not in (0.0, 1.0):
            continue
        out.append(_obs("echinocandin", r.get("compound_id", ""),
                        r.get("compound_name", ""), int(call), "chembl_serum",
                        "high", r.get("source_ref", "ChEMBL")))
    return out, True


# --------------------------------------------------------------------------- #
# consensus per compound
# --------------------------------------------------------------------------- #
def consensus(observations):
    # Key on normalized NAME so the same drug from different sources merges
    # (ChEMBL uses CHEMBL... ids, the in-repo tables use EXT-FKS.../PAPU... ids —
    # keying on id would split e.g. caspofungin into two phantom compounds).
    by = defaultdict(list)
    for o in observations:
        name_key = (o["compound_name"] or o["compound_id"]).strip().upper()
        by[(o["chemotype"], name_key)].append(o)

    # high-confidence = direct serum MIC / ChEMBL serum / clinical-approved
    DIRECT_SOURCES = {"papulacandin_serum_mic", "echinocandin_serum_mic",
                      "chembl_serum", "clinical_approved"}
    rows = []
    for (chemotype, key), obs in by.items():
        name = obs[0]["compound_name"]
        cid = obs[0]["compound_id"]
        direct = [o for o in obs if o["label_source"] in DIRECT_SOURCES]
        deciding = direct if direct else obs
        n_pos = sum(1 for o in deciding if o["serum_active"] == 1)
        n_neg = sum(1 for o in deciding if o["serum_active"] == 0)
        call = 1 if n_pos > n_neg else (0 if n_neg > n_pos
                                        else (1 if n_pos else 0))

        # Trustworthiness. Direct serum evidence is usable. An in-vivo/clinical
        # proxy is a reliable POSITIVE for echinocandins (serum-active drugs whose
        # in-vitro serum activity tracks animal outcomes, Nasar 2013) but NOT for
        # papulacandins, whose class phenotype is serum LOSS despite reported
        # in-vivo ED50s (papulacandin B is the archetype). So a papulacandin known
        # only from an in-vivo proxy has UNKNOWN serum activity — kept as a test
        # candidate, not a training label.
        if direct:
            tier, usable, caveat = "direct_serum", True, ""
        elif chemotype == "echinocandin":
            tier, usable, caveat = "proxy_echinocandin", True, "in-vivo/clinical positive"
        else:
            tier, usable, caveat = ("proxy_papulacandin_weak", False,
                                    "in-vivo only; serum activity UNKNOWN for this "
                                    "class — a test candidate, not a training label")
        rows.append({
            "chemotype": chemotype, "compound_id": cid, "compound_name": name,
            "serum_active": call,
            "label_tier": tier,
            "usable_for_training": usable,
            "confidence": "high" if direct else "medium",
            "n_active_obs": sum(1 for o in obs if o["serum_active"] == 1),
            "n_inactive_obs": sum(1 for o in obs if o["serum_active"] == 0),
            "mixed": bool(n_pos and n_neg),
            "sources": ";".join(sorted({o["label_source"] for o in obs})),
            "caveat": caveat,
        })
    return rows


# --------------------------------------------------------------------------- #
def _write(path, rows, fields=None):
    if not rows:
        open(path, "w").close()
        return
    fields = fields or list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _summary(path, obs, labels, chembl_used):
    def n_cpd(rows, chemo=None):
        return len({(r["chemotype"], r["compound_id"] or r["compound_name"])
                    for r in rows if chemo is None or r["chemotype"] == chemo})
    from collections import Counter
    src = Counter(o["label_source"] for o in obs)
    train = [r for r in labels if r["usable_for_training"]]
    tpos = sum(1 for r in train if r["serum_active"] == 1)
    tneg = sum(1 for r in train if r["serum_active"] == 0)
    weak = [r for r in labels if not r["usable_for_training"]]
    L = []
    L.append("# Stage 0 (binary reframe) — serum-active / inactive label set\n\n")
    L.append(
        "Endpoint reframed to BINARY `serum_active` (retains activity in serum: "
        "yes/no) to pool more, and more diverse, data than the quantitative "
        "serum-shift endpoint could — the right call for a discovery goal, and "
        "robust to the fact that serum inactivation is not fully computable "
        "(planning/serum_mechanism_evidence.md).\n\n")
    L.append("## Trustworthy training set (use this)\n\n")
    L.append(
        f"- **{len(train)} compounds with usable labels** ({tpos} serum-active, "
        f"{tneg} inactive) — direct serum-MIC evidence plus echinocandin "
        f"in-vivo/clinical positives.\n"
        f"- Direct serum-MIC compounds: "
        f"{sum(1 for r in labels if r['label_tier']=='direct_serum')} "
        f"({n_cpd([r for r in labels if r['label_tier']=='direct_serum'], 'papulacandin')} "
        f"papulacandin, "
        f"{n_cpd([r for r in labels if r['label_tier']=='direct_serum'], 'echinocandin')} "
        f"echinocandin).\n"
        f"- vs the **24** matched-pair compounds the quantitative endpoint could "
        f"use — a modest but real gain, and now the endpoint is the "
        f"discovery-relevant binary form.\n\n")
    L.append("## Caveated (NOT training labels)\n\n")
    L.append(
        f"- **{len(weak)} papulacandin compounds known only from an in-vivo "
        f"proxy**: serum activity is UNKNOWN for this class (in-vivo ED50s coexist "
        f"with serum loss — papulacandin B is the archetype). These are **test "
        f"candidates** for the wet-lab panel, not positive training labels — "
        f"pooling them would inject ~{len(weak)} likely-false positives.\n"
        f"- {sum(1 for r in labels if r['mixed'])} compounds are mixed (active in "
        f"some serum assays, lost in others).\n\n")
    chembl_obs = src.get("chembl_serum", 0)
    if chembl_used:
        # richest-evidence echinocandins after the ChEMBL merge
        rich = sorted(
            (r for r in labels if r["chemotype"] == "echinocandin"
             and "chembl_serum" in r["sources"]),
            key=lambda r: r["n_active_obs"] + r["n_inactive_obs"], reverse=True)[:3]
        rich_txt = "; ".join(
            f"{r['compound_name']} {r['n_active_obs']}+{r['n_inactive_obs']}"
            for r in rich)
        L.append("## What the ChEMBL pull added (ingested)\n\n")
        L.append(
            f"- **{chembl_obs} serum-context observations** merged. Its value is "
            "DEPTH, not width: the echinocandin labels are now data-rich, robust "
            f"calls rather than single points ({rich_txt} active+inactive obs).\n"
            "- It did **not** add many new distinct compounds — serum-context "
            "antifungal data clusters on the approved echinocandins, which we "
            "already had. The path to more *distinct* serum-labelled compounds "
            "remains the wet-lab panel (papulacandins).\n"
            "- Bonus: `stage0/data/chembl_free_fraction.csv` PPB/Fu rows enlarge the "
            "Stage-1b free-fraction seed beyond the original anidulafungin-only set.\n"
            f"- {len(obs)} total observations. By source: "
            + ", ".join(f"{k} {v}" for k, v in src.most_common()) + ".\n\n")
    else:
        L.append("## Where the real expansion comes from: ChEMBL\n\n")
        L.append(
            f"- Only {n_cpd([r for r in labels if r['label_tier']=='direct_serum'], 'echinocandin')} "
            "echinocandins have direct serum-MIC evidence in-repo. The binary "
            "endpoint can absorb many more echinocandin/FKS serum-shift records — "
            "the high-value ChEMBL pull (see drop-in below).\n"
            f"- {len(obs)} total observations. By source: "
            + ", ".join(f"{k} {v}" for k, v in src.most_common()) + ".\n\n")
    L.append("## ChEMBL drop-in\n\n")
    L.append(
        f"- ChEMBL echinocandin serum records: "
        f"{'INGESTED' if chembl_used else 'not present yet'}. "
        "Place a CSV at `stage0/data/chembl_echinocandin_serum.csv` "
        "(columns: compound_name, serum_active, source_ref) and re-run to merge "
        "them automatically. (Pending the ChEMBL connector.)\n\n")
    L.append("## Use\n\n")
    L.append(
        "`binary_serum_activity_labels.csv` is the pooled training/target set for a "
        "serum-active classifier and for panel selection. Keep the potency oracle "
        "(Stage 1a) as the filter and control potency as a covariate — the binary "
        "label is still potency-correlated within a series. `evidence_level` and "
        "`confidence` let you weight direct serum MICs above in-vivo/clinical "
        "proxies.\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("".join(L))


def main():
    os.makedirs(OUT, exist_ok=True)
    chembl_obs, chembl_used = chembl_drop_in()
    observations = (papulacandin_serum_mic() + echinocandin_serum_mic()
                    + invivo_proxy_positives() + clinical_approved_positives()
                    + chembl_obs)
    labels = consensus(observations)

    _write(os.path.join(OUT, "binary_serum_activity_observations.csv"), observations)
    _write(os.path.join(OUT, "binary_serum_activity_labels.csv"), labels)
    _summary(os.path.join(OUT, "binary_serum_activity_summary.md"),
             observations, labels, chembl_used)

    train = [r for r in labels if r["usable_for_training"]]
    tpos = sum(1 for r in train if r["serum_active"] == 1)
    weak = sum(1 for r in labels if not r["usable_for_training"])
    print(f"observations: {len(observations)} (chembl drop-in: "
          f"{'yes' if chembl_used else 'pending'})")
    print(f"trustworthy training labels: {len(train)} "
          f"({tpos} serum-active, {len(train)-tpos} inactive) vs 24 in the "
          f"quantitative endpoint")
    print(f"caveated papulacandin in-vivo-only test candidates (not labels): {weak}")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
