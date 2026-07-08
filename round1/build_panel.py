#!/usr/bin/env python3
"""Round 1 — informative synthesis/test panel selector (design of experiments).

The Stage-1c gate did not pass because the serum labels are one fusacandin
series (guardrail #8: the binding constraint is DATA). This selector turns the
oracles into a concrete "what to test first" list that maximally *generates* the
missing data, per planning/round1_synthesis_panel_plan.md.

Objective — maximise information about serum tolerance per cost, NOT predicted
potency (that would re-select the class's known serum failures). Concretely:

  * FILTER by the Stage-1a potency oracle: only compounds predicted active
    serum-free (a broth-dead compound can't demonstrate serum tolerance).
  * NOVELTY: only compounds with NO serum label yet (testing knowns wastes a slot).
  * DIVERSITY: span scaffolds and clogP (the mechanism-key lipophilicity /
    amphiphile-sequestration axis, planning/serum_mechanism_evidence.md), with a
    bonus for scaffolds ABSENT from the current serum set — those directly test
    whether serum tolerance generalises beyond fusacandins.
  * CHEAP FIRST: prefer obtainable-by-fermentation naturals and cheap
    semisynthetic analogs (synthesis_feasibility.csv) — round 1 needs no new
    chemistry, just the serum assay on a diverse obtainable set.
  * CONTROLS: include known serum-RETAINED and serum-LOST papulacandins plus an
    external serum-active echinocandin, so the round-1 assay is calibrated.

Output: a ranked panel with the design rationale, ready to hand to the wet lab
(protein-adjusted MIC ± albumin; equilibrium-dialysis fu). The measurements then
retrain the oracles (Stage 4 loop).

Run:
    python round1/build_panel.py
"""

from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "outputs")
sys.path.insert(0, os.path.join(ROOT, "stage1b"))
import featurize as F  # noqa: E402

COMPOUNDS = os.path.join(ROOT, "curated/core_tables/compounds_master.csv")
SYNTH = os.path.join(ROOT, "curated/core_tables/synthesis_feasibility.csv")
POTENCY = os.path.join(ROOT, "stage1a/outputs/potency_predictions.csv")
FU = os.path.join(ROOT, "stage1b/outputs/free_fraction_predictions.csv")
SERUM = os.path.join(ROOT, "stage0/outputs/binary_serum_activity_labels.csv")

N_CANDIDATES = 12
ACTIVE_PROB_MIN = 0.5


def _read(p):
    with open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _clogp_bin(x):
    if x is None:
        return "NA"
    return f"clogP<{int(x)+1}"          # 1-unit bins


# --------------------------------------------------------------------------- #
# assemble
# --------------------------------------------------------------------------- #
def load():
    cm = {r["compound_id"]: r for r in _read(COMPOUNDS)}
    synth = defaultdict(set)
    for r in _read(SYNTH):
        synth[r["compound_id"]].add(r["route_type"])
    pot = {r["compound_id"]: r for r in _read(POTENCY)}
    fu = {r["compound_id"]: r for r in _read(FU)}

    # serum-labelled papulacandins (direct evidence) — name-keyed, with call
    serum_by_name = {}
    for r in _read(SERUM):
        if r["chemotype"] == "papulacandin" and r["label_tier"] == "direct_serum":
            serum_by_name[r["compound_name"].upper()] = r["serum_active"]

    # scaffolds already represented in the serum set (mostly the fusacandin core)
    labelled_scaffolds = set()
    for cid, c in cm.items():
        if c["canonical_name"].upper() in serum_by_name:
            labelled_scaffolds.add(_scaffold(c))

    rows = []
    for cid, c in cm.items():
        if cid not in synth or cid not in pot:
            continue
        name = c["canonical_name"]
        pa = _num(pot[cid]["pred_active_prob"])
        routes = synth[cid]
        cost = ("cheap_semisynthetic" if "semisynthetic" in routes else
                "obtain_fermentation" if "fermentation_isolation" in routes else
                "synthesis")
        scaf = _scaffold(c)
        rows.append({
            "compound_id": cid, "name": name,
            "pred_active_prob": pa,
            "pred_pmic": _num(pot[cid]["pred_pmic"]),
            "clogp": _num(c.get("clogp")), "mw": _num(c.get("mw_exact")),
            "aromatic_rings": _aro(c),
            "scaffold": scaf,
            "cost": cost,
            "serum_label": serum_by_name.get(name.upper(), ""),  # '', '1', '0'
            "novel_scaffold": scaf not in labelled_scaffolds and scaf != "",
            "fu_pred": fu.get(cid, {}).get("pred_fu", ""),
            "fu_in_domain": fu.get(cid, {}).get("in_domain", ""),
        })
    return rows, labelled_scaffolds


def _scaffold(c):
    smi = (c.get("smiles_canonical") or c.get("smiles_raw") or "").strip()
    return F.murcko_scaffold(smi) if smi else ""


def _aro(c):
    smi = (c.get("smiles_canonical") or c.get("smiles_raw") or "").strip()
    feats = F.featurize(smi) if smi else None
    return feats[F.FEATURE_NAMES.index("aromatic_rings")] if feats else None


# --------------------------------------------------------------------------- #
# greedy diversity selection
# --------------------------------------------------------------------------- #
def select_candidates(rows, n=N_CANDIDATES):
    pool = [r for r in rows
            if r["serum_label"] == ""                       # serum-unknown
            and r["pred_active_prob"] is not None
            and r["pred_active_prob"] >= ACTIVE_PROB_MIN      # predicted active
            and r["cost"] in ("cheap_semisynthetic", "obtain_fermentation")]
    chosen, covered = [], set()
    while pool and len(chosen) < n:
        def gain(r):
            cell = (r["scaffold"], _clogp_bin(r["clogp"]))
            new_cell = cell not in covered
            return (
                new_cell,                                    # covers a new scaffold×clogP cell
                r["novel_scaffold"],                         # scaffold absent from serum set
                r["cost"] == "obtain_fermentation",          # cheapest (no synthesis)
                r["pred_active_prob"],                       # cleaner serum readout
            )
        best = max(pool, key=gain)
        chosen.append(best)
        covered.add((best["scaffold"], _clogp_bin(best["clogp"])))
        pool.remove(best)
    for r in chosen:
        r["panel_role"] = ("candidate_novel_scaffold" if r["novel_scaffold"]
                           else "candidate_diversity")
        r["selection_reason"] = _reason(r)
    return chosen


def _reason(r):
    bits = []
    if r["novel_scaffold"]:
        bits.append("scaffold absent from serum set (tests generalisation)")
    bits.append(f"clogP {r['clogp']:.1f}" if r["clogp"] is not None else "clogP NA")
    bits.append(r["cost"].replace("_", " "))
    bits.append(f"pAct {r['pred_active_prob']:.2f}")
    return "; ".join(bits)


def pick_controls(rows):
    """Known serum-retained + serum-lost papulacandins to calibrate the assay."""
    retained = [r for r in rows if r["serum_label"] == "1"]
    lost = [r for r in rows if r["serum_label"] == "0"]
    controls = []
    for r in sorted(retained, key=lambda z: -(z["pred_active_prob"] or 0))[:2]:
        r = dict(r); r["panel_role"] = "control_serum_retained"
        r["selection_reason"] = "known serum-RETAINED — positive assay control"
        controls.append(r)
    for r in sorted(lost, key=lambda z: -(z["pred_active_prob"] or 0))[:2]:
        r = dict(r); r["panel_role"] = "control_serum_lost"
        r["selection_reason"] = "known serum-LOST — negative assay control"
        controls.append(r)
    return controls


# --------------------------------------------------------------------------- #
def write_panel(panel):
    os.makedirs(OUT, exist_ok=True)
    fields = ["rank", "panel_role", "compound_id", "name", "scaffold_short",
              "novel_scaffold", "cost", "pred_active_prob", "pred_pmic", "clogp",
              "mw", "aromatic_rings", "fu_pred", "fu_in_domain", "serum_label",
              "selection_reason"]
    with open(os.path.join(OUT, "round1_synthesis_panel.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i, r in enumerate(panel, 1):
            w.writerow({
                "rank": i, "panel_role": r["panel_role"],
                "compound_id": r["compound_id"], "name": r["name"],
                "scaffold_short": (r["scaffold"] or "")[:24],
                "novel_scaffold": r.get("novel_scaffold", ""),
                "cost": r["cost"],
                "pred_active_prob": round(r["pred_active_prob"], 2) if r["pred_active_prob"] is not None else "",
                "pred_pmic": r["pred_pmic"], "clogp": round(r["clogp"], 2) if r["clogp"] is not None else "",
                "mw": round(r["mw"]) if r["mw"] is not None else "",
                "aromatic_rings": int(r["aromatic_rings"]) if r["aromatic_rings"] is not None else "",
                "fu_pred": r["fu_pred"], "fu_in_domain": r["fu_in_domain"],
                "serum_label": r["serum_label"], "selection_reason": r["selection_reason"],
            })


def write_report(candidates, controls, labelled_scaffolds):
    cand_scaf = {r["scaffold"] for r in candidates}
    novel = [r for r in candidates if r["novel_scaffold"]]
    clogps = [r["clogp"] for r in candidates if r["clogp"] is not None]
    L = []
    L.append("# Round 1 — informative serum-test panel\n\n")
    L.append(
        "Turns the oracles into a concrete first experiment that *generates* the "
        "scaffold-diverse serum data the Stage-1c gate lacked. Objective: maximum "
        "information about serum tolerance per cost — not predicted potency.\n\n")
    L.append("## The panel\n\n")
    L.append(
        f"- **{len(candidates)} test candidates** + **{len(controls)} controls** "
        f"(known serum-retained / serum-lost, to calibrate the assay).\n"
        f"- **{len(novel)} candidates sit on scaffolds ABSENT from the current "
        f"serum set** — these directly test whether serum tolerance generalises "
        f"beyond the fusacandin series (the gate's blind spot).\n"
        f"- Candidates span {len(cand_scaf)} scaffolds and clogP "
        f"{min(clogps):.1f}–{max(clogps):.1f} (the mechanism-key "
        f"lipophilicity/sequestration axis).\n"
        f"- All are obtainable without new synthesis (fermentation naturals + "
        f"cheap semisynthetic), so round 1 is just the serum assay on a diverse "
        f"obtainable set — the cheapest way to break the single-series limit.\n\n")
    L.append("## Why this design (and not 'make the predicted-best')\n\n")
    L.append(
        "Ranking by predicted potency would re-pick potent compounds that may die "
        "in serum — the class's 50-year failure. Instead the potency oracle is a "
        "FILTER (all candidates are predicted active serum-free, so a serum result "
        "is interpretable), and selection maximises scaffold + clogP diversity so "
        "the resulting labels can finally *separate serum tolerance from potency* "
        "and *generalise past fusacandins*.\n\n")
    L.append("## Assay & loop\n\n")
    L.append(
        "Measure serum-relevant endpoints directly: protein-adjusted MIC ± albumin "
        "(binary retained/lost + shift) and equilibrium-dialysis fu. Feed results "
        "back to Stage 0 (binary labels) and re-validate the Stage-1c gate. Round 2 "
        "then uses generative decoration to probe the clogP / aromatic-ester lever "
        "with potency-matched pairs.\n\n")
    L.append("## Note on the free-fraction column\n\n")
    L.append(
        "`fu_pred` is shown for reference but most candidates are OUT of the "
        "Stage-1b oracle's domain (large glycolipids) — treat fu as measure-not-"
        "predict (guardrails #4/#5). The equilibrium-dialysis fu in round 1 is what "
        "calibrates it.\n")
    with open(os.path.join(OUT, "round1_panel_report.md"), "w", encoding="utf-8") as fh:
        fh.write("".join(L))


def main():
    rows, labelled_scaffolds = load()
    candidates = select_candidates(rows)
    controls = pick_controls(rows)
    panel = candidates + controls
    write_panel(panel)
    write_report(candidates, controls, labelled_scaffolds)

    novel = sum(1 for r in candidates if r["novel_scaffold"])
    print(f"candidates: {len(candidates)} ({novel} on scaffolds absent from the "
          f"serum set); controls: {len(controls)}")
    print(f"candidate scaffolds: {len({r['scaffold'] for r in candidates})}; "
          f"clogP span "
          f"{min(r['clogp'] for r in candidates if r['clogp'] is not None):.1f}–"
          f"{max(r['clogp'] for r in candidates if r['clogp'] is not None):.1f}")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
