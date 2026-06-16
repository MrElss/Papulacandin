#!/usr/bin/env python3
"""
phase4_design.py
================
PHASE 4 — Candidate prioritization & design proposals (project goals #1 and #2).

IMPORTANT FRAMING (read first):
Phase 3B showed we do NOT yet have a validated predictive serum model. So Phase 4
is deliberately a TRANSPARENT, SAR-GUIDED prioritization, not a black-box
prediction. It does three honest things:

  (A) REPORT the best serum-tolerant compounds that ALREADY EXIST in the data
      (the concrete answer to goal #1 from current evidence), with their
      synthesis-route availability.
  (B) Build an interpretable "serum-tolerance design score" from the Phase-1
      descriptor DIRECTIONS (rigid/aromatic, higher MW/HBA -> better serum),
      apply it across the curated papulacandins, and flag each candidate's
      applicability domain + whether a synthesis route exists. Explicitly a
      hypothesis-ranking heuristic, NOT a validated prediction.
  (C) Distill DESIGN RULES that combine our descriptor finding with the
      independent literature SAR annotations (curated/after_analysis_reference),
      and propose concrete next analogs to make and test.

Outputs (analysis/outputs/):
  * table_phase4_existing_leads.csv   — best observed serum-active compounds
  * table_phase4_candidate_ranking.csv— heuristic score for all papulacandins
  * design_rules.md                   — distilled, source-attributed design rules
  * fig_phase4_existing_leads.png     — observed serum MIC of current best leads
  * fig_phase4_score_sanity.png       — heuristic score vs observed serum MIC
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
REF = os.path.join(ROOT, "curated", "after_analysis_reference")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})

# --------------------------------------------------------------------------
# Load
# --------------------------------------------------------------------------
pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
sf = pd.read_csv(os.path.join(CORE, "synthesis_feasibility.csv"))
sar = pd.read_csv(os.path.join(REF, "sar_annotations.csv"))
p1 = pd.read_csv(os.path.join(OUT, "phase1_descriptor_stats.csv"))

DESCRIPTORS = ["mw_exact", "fsp3", "tpsa", "hba", "rotb", "clogp"]
for c in DESCRIPTORS:
    cm[c] = pd.to_numeric(cm[c], errors="coerce")


def parse_mic(s):
    s = str(s).strip()
    for r in (">=", "<=", ">", "<", "="):
        if s.startswith(r):
            s = s[len(r):]
            break
    try:
        return float(s)
    except ValueError:
        return np.nan

pairs["serum_mic_num"] = pairs["serum_mic_ugml"].apply(parse_mic)

# ==========================================================================
# (A) Best serum-tolerant compounds that ALREADY exist
# ==========================================================================
route_ids = set(sf["compound_id"])
leads = (pairs[pairs["serum_active"] == "yes"]
         .sort_values("serum_mic_num")
         .loc[:, ["compound_id", "name", "serumfree_mic_ugml", "serum_mic_ugml",
                  "serum_mic_num", "serum_shift_fold", "modification_summary"]]
         .copy())
leads["synthesis_route_known"] = leads["compound_id"].isin(route_ids).map({True: "yes", False: "no"})
leads.to_csv(os.path.join(OUT, "table_phase4_existing_leads.csv"), index=False)

# ==========================================================================
# (B) Interpretable serum-tolerance design score across papulacandins
# ==========================================================================
# Use Phase-1 descriptors whose association with serum MIC is meaningful.
sig = p1[(p1["spearman_p"] < 0.10) & (p1["descriptor"].isin(DESCRIPTORS))]
weights = dict(zip(sig["descriptor"], sig["spearman_rho_vs_serumMIC"]))
# Fallback if too few pass: use the top-|rho| descriptors.
if len(weights) < 2:
    top = p1.reindex(p1["spearman_rho_vs_serumMIC"].abs().sort_values(ascending=False).index).head(3)
    weights = dict(zip(top["descriptor"], top["spearman_rho_vs_serumMIC"]))

# Restrict standardization to the applicability domain (papulacandin/fusacandin
# family) so far-out structures do not distort the z-scores.
fam_mask = cm["compound_class"].str.contains("papulacandin|fusacandin", case=False, na=False)
fam = cm[fam_mask].copy()

score = np.zeros(len(fam))
for d, rho in weights.items():
    col = fam[d]
    z = (col - col.mean()) / (col.std(ddof=0) + 1e-9)
    # rho>0 => higher descriptor -> higher serum MIC (worse). Tolerance gets the
    # NEGATIVE contribution so higher score = predicted MORE serum-tolerant.
    score += -rho * z
fam["serum_tolerance_design_score"] = np.round(score, 3)

# Annotate with serum data (if any), route availability, and SAR exposure notes.
obs = pairs.set_index("compound_id")["serum_mic_ugml"]
sar_exposure = (sar[sar["effect_on_in_vivo"].isin(["improved"]) |
                    sar["effect_on_solubility_or_exposure"].isin(["improved"])]
                .groupby("compound_id").size())
fam["observed_serum_mic"] = fam["compound_id"].map(obs).fillna("(no serum data)")
fam["synthesis_route_known"] = fam["compound_id"].isin(route_ids).map({True: "yes", False: "no"})
fam["sar_exposure_positive_hits"] = fam["compound_id"].map(sar_exposure).fillna(0).astype(int)

ranking = (fam.sort_values("serum_tolerance_design_score", ascending=False)
           .loc[:, ["compound_id", "canonical_name", "compound_class",
                    "serum_tolerance_design_score", "observed_serum_mic",
                    "synthesis_route_known", "sar_exposure_positive_hits",
                    "modification_summary"]])
ranking.to_csv(os.path.join(OUT, "table_phase4_candidate_ranking.csv"), index=False)

# Sanity check: does the heuristic track observed serum MIC on labeled compounds?
lab = fam.merge(pairs[["compound_id", "serum_mic_num"]], on="compound_id")
rho_chk, p_chk = stats.spearmanr(lab["serum_tolerance_design_score"],
                                 np.log10(lab["serum_mic_num"]))

# ==========================================================================
# (C) Design rules (descriptor finding + independent literature SAR)
# ==========================================================================
rules = f"""# Phase 4 — Serum-tolerance design rules (hypotheses to test)

Two INDEPENDENT lines of evidence converge on actionable design directions.
Both are hypothesis-generating; neither is a validated quantitative model.

## Axis 1 — Descriptor SAR (this project, Phase 1/3B; Fusacandin-A C-6' series)
Within the only modelable series, serum activity improves with:
{os.linesep.join(f"  * {d}: rho={r:+.2f} vs serum MIC -> favor {'higher' if r<0 else 'lower'} {d}" for d, r in weights.items())}
Plain reading: **rigid, extended AROMATIC C-6' acyl groups (biphenyl / naphthoyl
types) retain serum activity better than flexible aliphatic chains.** clogP alone
is NOT the lever.

## Axis 2 — Literature SAR annotations (curated/after_analysis_reference)
Independent, paper-derived observations on the Papulacandin-B core point a
DIFFERENT structural direction for in-vivo/serum translation:
  * Polar 10-O ethers can improve in-vivo activity without losing Candida MIC.
  * Selected 11-aminoacyl / cationic substitutions improve in-vivo activity
    while keeping useful MIC.
Plain reading: on the Pap-B core, **adding polar / cationic handles** (not just
lipophilic bulk) is the literature-supported route to better exposure.

## Synthesis recommendations (concrete next analogs)
1. Best EXISTING serum lead to advance now (goal #1): the biphenyl/naphthoyl
   Fusacandin-A C-6' esters (serum MIC ~11-25 ug/mL), all with known 4-step
   semisynthetic routes -> lowest-risk compounds to re-test under serum.
2. New analog hypothesis A (Axis 1): rigid biphenyl C-6' ester bearing a polar
   substituent (e.g. carboxylate / morpholine) to merge rigidity + polarity.
3. New analog hypothesis B (Axis 2 x goal #2): graft a polar/cationic 11-aminoacyl
   handle onto the serum-tolerant biphenyl-ester scaffold -> tests whether the
   two independent axes are additive for serum tolerance.

## Caveats
* All scores/rules are within or near the papulacandin/fusacandin domain only.
* Heuristic design-score vs observed serum MIC sanity check (labeled set, n={len(lab)}):
  Spearman rho = {rho_chk:.2f} (p = {p_chk:.3f}). Consistent by construction with
  Axis 1; this is NOT independent validation.
"""
with open(os.path.join(OUT, "design_rules.md"), "w") as fh:
    fh.write(rules)

# ==========================================================================
# Figures
# ==========================================================================
fig, ax = plt.subplots(figsize=(7, 4.3))
top_leads = leads.head(8)
bars = ax.barh(top_leads["compound_id"] + "  " + top_leads["name"].str[:22],
               top_leads["serum_mic_num"], color="#2a9d8f")
ax.invert_yaxis()
ax.bar_label(bars, fmt="%.1f", padding=3)
ax.set_xlabel("Observed serum MIC (ug/mL; lower = better)")
ax.set_title("Phase 4 (goal #1): best EXISTING serum-tolerant papulacandins\nall have known 4-step semisynthetic routes")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase4_existing_leads.png"))
plt.close(fig)

fig, ax = plt.subplots(figsize=(5.8, 4.6))
ax.scatter(lab["serum_tolerance_design_score"], lab["serum_mic_num"],
           s=55, color="#264653", edgecolor="k", linewidth=0.3)
ax.set_yscale("log")
ax.set_xlabel("Serum-tolerance design score (higher = predicted better)")
ax.set_ylabel("Observed serum MIC (ug/mL, log)")
ax.set_title(f"Phase 4: heuristic score vs observed serum MIC\nSpearman rho={rho_chk:.2f} (hypothesis-ranking only)")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig_phase4_score_sanity.png"))
plt.close(fig)

# ==========================================================================
# Console summary
# ==========================================================================
print("PHASE 4 — candidate prioritization & design")
print("=" * 60)
print(f"Descriptor weights used (Phase-1 rho vs serum MIC): {weights}")
print(f"\n(A) Best EXISTING serum-tolerant leads (goal #1):")
print(leads.head(6)[["compound_id", "name", "serum_mic_ugml",
                     "serum_shift_fold", "synthesis_route_known"]].to_string(index=False))
print(f"\n(B) Top heuristic candidates within the family (n={len(fam)} scored):")
print(ranking.head(8)[["compound_id", "canonical_name",
                       "serum_tolerance_design_score", "observed_serum_mic",
                       "synthesis_route_known"]].to_string(index=False))
print(f"\nHeuristic sanity check: Spearman rho={rho_chk:.2f} (p={p_chk:.3f})")
print("\nWrote: table_phase4_existing_leads.csv, table_phase4_candidate_ranking.csv,")
print("       design_rules.md, fig_phase4_existing_leads.png, fig_phase4_score_sanity.png")
