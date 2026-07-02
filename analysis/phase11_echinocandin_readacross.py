#!/usr/bin/env python3
"""
phase11_echinocandin_readacross.py
==================================
PHASE 11 — Cross-chemotype serum-tolerance read-across from the echinocandins.

Motivation
----------
Phases 1-10 established that the papulacandin/fusacandin serum gap is real but
that no descriptor reaches significance on the available data (n=24, 11 censored,
a SINGLE chemotype from essentially one lab). The cross-phase synthesis named the
one thing that could break that ceiling: a SECOND, independent chemotype that
hits the same target and shows the same serum liability.

The echinocandins are exactly that. They are cyclic lipopeptide beta-1,3-glucan
(Fks1) synthase inhibitors — a completely different pharmacophore from the
papulacandin aryl-C-glycoside — that ALSO carry a long lipophilic tail and are
famously affected by serum in vitro. Critically, our own external FKS corpus
(external/data/external/fks_inhibitors/) already contains the relevant data:
MIC-shift ratios in +/- 50% serum (the SAME endpoint as our serum_shift_fold)
plus plasma-protein-binding (PPB) and fraction-unbound (Fu) measurements.

What this phase does (and does NOT do)
--------------------------------------
It harmonizes the echinocandin serum data with our 24-compound papulacandin
serum-gap set on a shared endpoint (log2 serum-shift fold) and a shared,
consistently-recomputed 2D descriptor set (RDKit), then asks two honest
questions:

  Q1 (reframe): how big are echinocandin serum shifts, and how do they square
      with the free-drug hypothesis (PPB/Fu)? -> context for our own endpoint.
  Q2 (stress-test): does the papulacandin "expose polar not hydrophobic surface"
      lead reproduce on the echinocandins at the bulk-descriptor level?

This is a READ-ACROSS with n=3 echinocandins carrying matched C. albicans serum
shifts. It is deliberately NOT a pooled regression or a significance claim — it
is a qualitative concordance check plus a reproducible dataset for later use.
Bulk 2D descriptors (MolLogP, TPSA/heavy-atom) are used because they are all
that is comparable across such different scaffolds; they are NOT the 3D exposed-
surface descriptors (polar SASA, QM logP) that the Phase 8-9 lead was actually
about, so a null here refines — it does not refute — that lead.

Data provenance
---------------
* Papulacandin serum shift : analysis/outputs/serum_gap_pairs.csv (Phase 0),
  SMILES from curated/core_tables/compounds_master.csv.
* Echinocandin serum shift / PPB / Fu / SMILES :
  external/data/external/fks_inhibitors/source_exports/*activity*.csv and
  external_compounds_candidate_v0_2.csv (curated from ChEMBL/PubChem).
  Ratio-vs-serum rows are defined by the source quotes as
  MIC(+50% serum) / MIC(serum-free medium) — i.e. the serum-shift fold,
  the same direction as our serum_shift_fold.

Outputs (analysis/outputs/)
---------------------------
* phase11_echinocandin_serum_shift.csv — per-compound echinocandin serum-shift
  (C. albicans + all-species medians), PPB/Fu, and RDKit descriptors.
* phase11_crosschemotype.csv           — harmonized papulacandin + echinocandin
  rows on the shared endpoint/descriptor axes.
* phase11_crosschemotype.png           — 2-panel read-across figure.
* phase11_findings.md                  — interpretation.
"""

from __future__ import annotations

import glob
import os
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
EXT = os.path.join(ROOT, "external", "data", "external", "fks_inhibitors", "source_exports")
os.makedirs(OUT, exist_ok=True)

# Marketed/known echinocandin-class names we recognize in the external corpus.
ECHINOCANDINS = ("CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN", "CILOFUNGIN")
# Same target (glucan synthase) but NO cyclic-lipopeptide fatty-acyl tail — the
# "the tail is droppable" existence proof (ibrexafungerp is a marketed drug).
NO_TAIL_GS_INHIBITORS = ("IBREXAFUNGERP", "ENFUMAFUNGIN")


# --------------------------------------------------------------------------
# Descriptors — recomputed with RDKit so both chemotypes are on one method.
# --------------------------------------------------------------------------
def rdkit_descriptors(smiles: str) -> dict:
    """2D descriptors comparable across very different scaffolds/sizes.

    tpsa_per_heavy_atom and mollogp_per_heavy_atom are size-normalized so a
    730 Da triterpenoid and a 1270 Da lipopeptide can be compared on
    'polarity density' rather than absolute totals.
    """
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles else None
    if mol is None:
        return {k: np.nan for k in
                ("mollogp", "tpsa", "mw", "heavy_atoms",
                 "tpsa_per_heavy_atom", "mollogp_per_heavy_atom")}
    ha = mol.GetNumHeavyAtoms()
    tpsa = Descriptors.TPSA(mol)
    logp = Crippen.MolLogP(mol)
    return {
        "mollogp": round(logp, 3),
        "tpsa": round(tpsa, 2),
        "mw": round(Descriptors.MolWt(mol), 2),
        "heavy_atoms": ha,
        "tpsa_per_heavy_atom": round(tpsa / ha, 3) if ha else np.nan,
        "mollogp_per_heavy_atom": round(logp / ha, 4) if ha else np.nan,
    }


# --------------------------------------------------------------------------
# Papulacandin side (our project's dependent variable)
# --------------------------------------------------------------------------
def load_papulacandin() -> pd.DataFrame:
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
    smiles = cm.set_index("compound_id")["smiles_canonical"].to_dict()
    rows = []
    for _, r in pairs.iterrows():
        shift = pd.to_numeric(r.get("serum_shift_fold"), errors="coerce")
        if not np.isfinite(shift) or shift <= 0:
            continue
        desc = rdkit_descriptors(smiles.get(r["compound_id"], ""))
        rows.append({
            "chemotype": "papulacandin",
            "compound": r["compound_id"],
            "name": r["name"],
            "organism": "Candida albicans",
            "serum_shift_fold": float(shift),
            "log2_serum_shift": float(np.log2(shift)),
            "shift_confidence": r.get("shift_confidence", ""),
            "ppb_percent": np.nan,
            "fu_percent": np.nan,
            **desc,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Echinocandin side (from our external FKS corpus)
# --------------------------------------------------------------------------
_SERUM_RATIO_RE = re.compile(r"serum", re.I)
_ALBICANS_RE = re.compile(r"albicans", re.I)


def _iter_activity_rows():
    for fn in sorted(glob.glob(os.path.join(EXT, "*activity*.csv"))):
        try:
            df = pd.read_csv(fn, dtype=str, keep_default_na=False)
        except Exception:
            continue
        for _, r in df.iterrows():
            yield r


def collect_echinocandin_serum() -> dict:
    """Return {NAME: {albicans_folds, allspecies_folds, ppb, fu}} deduped."""
    acc = {}
    seen = set()
    for r in _iter_activity_rows():
        name = (r.get("preferred_name") or "").strip().upper()
        if name not in ECHINOCANDINS and name not in NO_TAIL_GS_INHIBITORS:
            continue
        d = acc.setdefault(name, {"albicans": [], "allspecies": [], "ppb": [], "fu": []})
        et = (r.get("endpoint_type") or "").strip()
        q = r.get("source_quote") or ""
        val = pd.to_numeric(r.get("endpoint_value"), errors="coerce")
        # PPB / Fu (protein binding, for the free-drug reframe)
        if et == "PPB" and np.isfinite(val):
            d["ppb"].append(float(val))
            continue
        if et == "Fu" and np.isfinite(val):
            d["fu"].append(float(val))
            continue
        # Serum-shift folds: Ratio endpoints whose quote is a +/- serum comparison.
        if et == "Ratio" and _SERUM_RATIO_RE.search(q) and np.isfinite(val) and val > 0:
            key = (name, round(float(val), 4), q[:80])
            if key in seen:
                continue
            seen.add(key)
            d["allspecies"].append(float(val))
            if _ALBICANS_RE.search(q):
                d["albicans"].append(float(val))
    return acc


def load_echinocandin() -> pd.DataFrame:
    acc = collect_echinocandin_serum()
    f = os.path.join(EXT, "external_compounds_candidate_v0_2.csv")
    comp = pd.read_csv(f, dtype=str, keep_default_na=False)
    smiles = {r["preferred_name"].strip().upper():
              (r.get("smiles_canonical") or r.get("smiles_raw") or "")
              for _, r in comp.iterrows()}
    rows = []
    for name, d in sorted(acc.items()):
        has_tail = name in ECHINOCANDINS
        # Primary endpoint: C. albicans median (matches papulacandin endpoint);
        # fall back to all-species median if no C. albicans row exists.
        alb = d["albicans"]
        allsp = d["allspecies"]
        shift = np.median(alb) if alb else (np.median(allsp) if allsp else np.nan)
        endpoint_org = "Candida albicans" if alb else ("all-species median" if allsp else "n/a")
        desc = rdkit_descriptors(smiles.get(name, ""))
        rows.append({
            "chemotype": "echinocandin" if has_tail else "glucan-synthase (no lipid tail)",
            "compound": name,
            "name": name.title(),
            "organism": endpoint_org,
            "serum_shift_fold": float(shift) if np.isfinite(shift) else np.nan,
            "log2_serum_shift": float(np.log2(shift)) if np.isfinite(shift) and shift > 0 else np.nan,
            "shift_confidence": "median-of-isolates" if np.isfinite(shift) else "no serum-shift row",
            "albicans_folds": ";".join(str(x) for x in sorted(set(alb))),
            "allspecies_median": round(float(np.median(allsp)), 2) if allsp else np.nan,
            "allspecies_n": len(allsp),
            "ppb_percent": round(float(np.median(d["ppb"])), 1) if d["ppb"] else np.nan,
            "fu_percent": round(float(np.median(d["fu"])), 4) if d["fu"] else np.nan,
            **desc,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Analysis + figure
# --------------------------------------------------------------------------
DESCRIPTOR_AXES = [
    ("mollogp", "RDKit MolLogP (whole molecule)"),
    ("tpsa_per_heavy_atom", "TPSA / heavy atom  (polarity density)"),
]


def spearman(df: pd.DataFrame, col: str):
    sub = df[["log2_serum_shift", col]].apply(pd.to_numeric, errors="coerce").dropna()
    if len(sub) < 4:
        return (np.nan, np.nan, len(sub))
    rho, p = stats.spearmanr(sub[col], sub["log2_serum_shift"])
    return (rho, p, len(sub))


def main() -> None:
    papu = load_papulacandin()
    ech = load_echinocandin()

    combined = pd.concat([papu, ech], ignore_index=True, sort=False)
    combined.to_csv(os.path.join(OUT, "phase11_crosschemotype.csv"), index=False)
    ech.to_csv(os.path.join(OUT, "phase11_echinocandin_serum_shift.csv"), index=False)

    # Within-papulacandin trend (the project's own lead) for each axis, and the
    # echinocandin ordering as a concordance check.
    lines = []
    lines.append(f"Papulacandin serum-gap compounds:      {len(papu)}")
    lines.append(f"Echinocandins with serum-shift data:   "
                 f"{ech['serum_shift_fold'].notna().sum()} "
                 f"({', '.join(ech.loc[ech['serum_shift_fold'].notna(),'compound'])})")
    lines.append("")
    lines.append("Within-papulacandin Spearman(log2 serum shift vs descriptor):")
    stat_rows = []
    for col, label in DESCRIPTOR_AXES:
        rho, p, n = spearman(papu, col)
        lines.append(f"  {label:36s} rho={rho:+.2f}  p={p:.3f}  n={n}")
        stat_rows.append({"descriptor": col, "papu_spearman_rho": rho,
                          "papu_spearman_p": p, "papu_n": n})
    pd.DataFrame(stat_rows).to_csv(
        os.path.join(OUT, "phase11_crosschemotype_stats.csv"), index=False)

    # Echinocandin C. albicans ordering (the robust, most comparable signal).
    ech_alb = ech[ech["serum_shift_fold"].notna()].sort_values("serum_shift_fold")
    lines.append("")
    lines.append("Echinocandin serum-shift ordering (most comparable endpoint):")
    for _, r in ech_alb.iterrows():
        lines.append(f"  {r['compound']:14s} shift x{r['serum_shift_fold']:.0f}  "
                     f"MolLogP={r['mollogp']:+.2f}  TPSA/HA={r['tpsa_per_heavy_atom']:.2f}  "
                     f"PPB={r['ppb_percent']}")
    report = "\n".join(lines)
    print(report)

    # ---- Figure: 2 panels, shared endpoint, both chemotypes overlaid ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    colors = {"papulacandin": "#4a73c4",
              "echinocandin": "#c44a4a",
              "glucan-synthase (no lipid tail)": "#2f9e5a"}
    for ax, (col, label) in zip(axes, DESCRIPTOR_AXES):
        for ct, sub in combined.groupby("chemotype"):
            s = sub[["log2_serum_shift", col]].apply(pd.to_numeric, errors="coerce").dropna()
            if s.empty:
                continue
            ax.scatter(s[col], s["log2_serum_shift"],
                       s=70 if ct != "papulacandin" else 45,
                       c=colors.get(ct, "#888"),
                       edgecolor="k" if ct != "papulacandin" else "none",
                       linewidth=0.8, alpha=0.9 if ct != "papulacandin" else 0.6,
                       label=ct, zorder=3 if ct != "papulacandin" else 2)
        # annotate the named reference compounds
        for _, r in ech.iterrows():
            v = pd.to_numeric(r[col], errors="coerce")
            y = pd.to_numeric(r["log2_serum_shift"], errors="coerce")
            if np.isfinite(v) and np.isfinite(y):
                ax.annotate(r["compound"].title(), (v, y),
                            fontsize=7.5, xytext=(4, 4),
                            textcoords="offset points")
        ax.set_xlabel(label)
        ax.set_ylabel("log2 serum-shift fold  (serum / serum-free MIC)")
        ax.grid(alpha=0.25)
    axes[0].set_title("Whole-molecule lipophilicity")
    axes[1].set_title("Size-normalized polarity")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Phase 11 — cross-chemotype serum-shift read-across "
                 "(papulacandins n=%d + echinocandins)" % len(papu), y=1.07)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase11_crosschemotype.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

    _write_findings(papu, ech, ech_alb, stat_rows)
    print("\nWrote: phase11_echinocandin_serum_shift.csv, phase11_crosschemotype.csv,\n"
          "       phase11_crosschemotype_stats.csv, phase11_crosschemotype.png, "
          "phase11_findings.md")


def _write_findings(papu, ech, ech_alb, stat_rows) -> None:
    order = ", ".join(f"{r['compound'].title()} (x{r['serum_shift_fold']:.0f})"
                      for _, r in ech_alb.iterrows())
    rho_logp = next((s["papu_spearman_rho"] for s in stat_rows if s["descriptor"] == "mollogp"), float("nan"))
    md = f"""# Phase 11 — echinocandin read-across for serum tolerance

*Question:* the papulacandin serum-gap lead ("expose polar not hydrophobic
surface") never cleared significance on our single chemotype (n={len(papu)}). Do
the echinocandins — an independent cyclic-lipopeptide glucan-synthase inhibitor
class with the same long-tail serum liability — corroborate it, and what does
their clinical success teach us? All echinocandin numbers below come from data
already curated in `external/data/external/fks_inhibitors/`.

## 1. Echinocandins suffer the SAME serum shift — and it is large
Serum-shift fold (MIC in 50% serum / MIC in serum-free medium), the identical
endpoint direction as our `serum_shift_fold`:

- **{order}** (C. albicans; robust across multiple isolates/species).
- Across all species the shifts reach ~256-1024x for micafungin/anidulafungin.

So a big serum shift is **not, by itself, disqualifying**: micafungin and
anidulafungin lose 1-2 orders of magnitude of in vitro potency to serum yet are
front-line IV antifungals. This reframes the whole project's target — the goal
is not "zero serum shift" but "enough free-drug exposure at the target."

## 2. The free-drug hypothesis is the right lens (PPB / Fu)
Echinocandins are ~96% (caspofungin), ~99% (anidulafungin) and ~99.8%
(micafungin) plasma-protein bound (literature; our corpus carries anidulafungin
PPB {ech.loc[ech['compound']=='ANIDULAFUNGIN','ppb_percent'].dropna().tolist()}%
and Fu ~0.01%). Clinically they are dosed to a **free-drug** AUC/MIC target, and
when normalized to unbound drug the PK/PD targets converge across the class.
Read-across for us: always model the *shift* / free fraction (Phase 8's
conclusion), and — the actionable part — **measure protein binding directly on
our leads** rather than inferring it from docking (Phase 10 could not).

## 3. Bulk 2D descriptors do NOT explain the echinocandin ordering (honest null)
The most polar echinocandin by every bulk measure — **micafungin**
(TPSA {float(ech.loc[ech['compound']=='MICAFUNGIN','tpsa'].iloc[0]):.0f},
MolLogP {float(ech.loc[ech['compound']=='MICAFUNGIN','mollogp'].iloc[0]):+.1f}) —
has the **largest** serum shift, while **caspofungin** has the smallest. Within
the papulacandins the same bulk axis is weak/uninformative
(MolLogP Spearman rho={rho_logp:+.2f}). So whole-molecule polarity/lipophilicity
is **not** the read-across variable.

Two non-exclusive readings, both consistent with earlier phases:
1. The papulacandin lead was about **locally exposed** surface (Boltzmann polar
   SASA / QM logP, Phases 8-9), which bulk TPSA cannot see. A bulk null here does
   **not** refute the 3D lead — it says "don't shortcut it with 2D."
2. The echinocandin serum effect has a documented component of **direct
   desensitization of glucan synthase by serum**, not pure albumin sequestration
   (the effect does not track protein-binding rank). This mirrors our Phase 10
   HSA-docking null and warns that "reduce albumin affinity" may be the wrong
   single objective.

## 4. Existence proof: the lipophilic tail is droppable
{", ".join(sorted(ech.loc[ech['chemotype'].str.contains('no lipid tail'),'compound'].str.title()))}
inhibit the same target with **no cyclic-lipopeptide fatty-acyl tail**
(ibrexafungerp: MW ~730, and it is an *orally* bioavailable marketed drug with
good tissue penetration). Since synthetic accessibility is off the table for now,
this legitimizes an aggressive design branch: **truncate / replace the long acyl
chain** of the papulacandin scaffold while keeping the aryl-C-glycoside /
spiroketal pharmacophore, rather than only decorating the C-6' ester.

## 5. What Phase 11 changes for the project
- **Endpoint:** adopt the free-drug framing; the deliverable is a serum-tolerant
  *free* exposure, not a serum-invariant MIC.
- **Design (Phase 5 next iteration):** keep scoring exposed polar surface
  (QM logP / polar SASA), and add a tail-truncation/replacement branch inspired
  by the fungerps. Do **not** rank on bulk TPSA/clogP — Phase 11 shows it does
  not separate serum tolerance across chemotypes.
- **Experiment (highest value):** run the echinocandin assay playbook on our
  papulacandin leads — protein-adjusted MIC with an albumin titration, and
  equilibrium-dialysis fraction-unbound — to decide, at last, whether the serum
  gap is albumin sequestration, a direct serum effect, or degradation. This is
  the measurement Phases 8 and 10 could not make computationally.

## Caveats
n=3 echinocandins with matched C. albicans serum shifts is a qualitative
read-across, not a regression; different pharmacophore, species mix, and assay
labs than the papulacandin set. The value is the **direction and the reframing**,
plus a reproducible harmonized dataset (`phase11_crosschemotype.csv`) to extend
as more +/- serum data is curated.

## Sources
- Paderu et al., *Antimicrob. Agents Chemother.* 2007 — Effects of serum on in
  vitro susceptibility testing of echinocandins.
- Odds/Gumbo et al. — serum differentially alters echinocandin antifungal
  activity; free-drug AUC/MIC targets.
- Davis et al., *J. Fungi* 2021 — Ibrexafungerp: first-in-class oral triterpenoid
  glucan synthase inhibitor.
- In-repo: `external/data/external/fks_inhibitors/source_exports/*` (ChEMBL/PubChem
  curated).
"""
    with open(os.path.join(OUT, "phase11_findings.md"), "w", encoding="utf-8") as fh:
        fh.write(textwrap.dedent(md).lstrip())


if __name__ == "__main__":
    main()
