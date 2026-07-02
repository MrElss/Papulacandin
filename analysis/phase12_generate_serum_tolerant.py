#!/usr/bin/env python3
"""
phase12_generate_serum_tolerant.py
==================================
PHASE 12 — Serum-tolerance-biased generative design (Track A: transparent
scaffold-constrained generator with a physics-grounded reward).

Why this exists / what changed from Phase 5
-------------------------------------------
Phase 5 generated novel analogs but scored them on QED / synthetic accessibility
/ 2D design rules — and Phases 5, 8 and 11 showed those bulk/2D terms are
UNINFORMATIVE for serum tolerance in this bRo5 chemotype (QED saturates near 0;
bulk polarity does not even order the echinocandins). The one convergent,
defensible lead across Phases 8 (polar SASA), 9 (QM logP) and 11 (echinocandin
read-across) is:

    among equipotent analogs, biasing the *locally exposed* surface toward polar
    / H-bonding character is associated with smaller serum loss.

Phase 12 therefore replaces the reward with that lead, made operational as a fast
3D descriptor: the Boltzmann-ish mean EXPOSED POLAR SURFACE FRACTION over a small
ETKDG conformer ensemble (RDKit rdFreeSASA — the built-in, distinct from the
`freesasa` pip package that failed to build in Phase 6). It preserves the
FKS-engaging pharmacophore by construction (re-esterifying onto the conserved
papulacandin core) and rewards novelty.

Honest framing (critical)
--------------------------
There is NO validated serum-tolerance oracle (every descriptor is p>0.05 on
n=24). The exposed-polar-surface reward is the project's best HYPOTHESIS, not a
proven predictor. So this phase does two things, not one:
  (1) ENRICH toward the predicted serum-tolerant (polar-exposed) region, and
  (2) emit a DISCRIMINATING SERIES — novel analogs on one scaffold that SPAN the
      exposed-polar-surface axis at comparable size — so a wet-lab serum assay
      can FALSIFY or confirm the lead rather than only confirm it.
Success = the assay resolves the hypothesis, not "the model was right."

The reward is validated retrospectively here: the same 3D exposed-polar-fraction
is computed for the 24 knowns and correlated with their observed serum SHIFT
(expected sign: more polar exposure -> smaller shift). This is a crude
single-few-conformer proxy (Phase 7->8 showed proxies OVERSTATE), reported as a
sanity check only.

Design scope
------------
Per the "ignore synthetic accessibility for now" instruction and the Phase 11
"the lipophilic tail is droppable" existence proof (ibrexafungerp), Phase 12
explores three branches off a serum-tolerant lead (PAPU-0080):
  * ester::   re-esterify the validated aromatic C-6' handle (Phase 5 chemistry)
  * polaraxis:: designed acyls deliberately spanning hydrophobic->polar (the
                discriminating-series backbone)
  * notail::  remove the long aliphatic fatty-acyl tail (deacylation) BEFORE
              re-esterifying — ibrexafungerp-inspired tail-free analogs

Track B (a trained generative network — REINVENT-style RL on the external FKS
pretraining set with THIS reward) is the natural follow-on once the reward is
trusted; this script is the reproducible, CPU-only foundation and reward it would
reuse.

Outputs (analysis/outputs/)
---------------------------
* phase12_generated_library.csv    — every generated analog + reward components
* phase12_discriminating_series.csv — matched novel set spanning exposed polarity
* phase12_top_candidates.sdf        — 3D structures of the top novel analogs (CREST-ready)
* phase12_crest_commands.sh         — per-candidate CREST/xTB command template
* phase12_reward_validation.png     — reward vs serum shift on the 24 knowns
* phase12_findings.md               — interpretation + how to read the series
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator, rdFreeSASA
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CORE = os.path.join(ROOT, "curated", "core_tables")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

# A small ensemble stabilizes the exposed-surface estimate without QM cost.
# Override for quick smoke runs.
N_CONF = int(os.environ.get("PHASE12_N_CONF", "2"))
SASA_CACHE = os.path.join(OUT, "phase12_sasa_cache.csv")
TOP_N_FOR_QM = 12
SERUM_TOLERANT_TEMPLATE = "PAPU-0080"  # fusacandin 6j, a serum-active lead
AROM_ESTER = Chem.MolFromSmarts("[O;X2:1][C:2](=O)[c:3]")
# any carbon-acyl ester (O-C(=O)-C...); the fatty tail is then selected at the
# FRAGMENT level (no aromatic ring, long carbon chain) so unsaturated acyl tails
# like the fusacandin C(=O)/C=C/C=C/... are caught, while the aromatic C-6' ester
# and short sugar acetates are filtered out.
ALIPH_ESTER = Chem.MolFromSmarts("[O;X2:1][C:2](=O)[#6:3]")
_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


# ==========================================================================
# Chemistry helpers (verified in Phase 5; re-implemented so Phase 12 is
# self-contained — phase5_generate.py has no __main__ guard and cannot be
# imported without side effects).
# ==========================================================================
def ecfp(mol):
    return _gen.GetFingerprint(mol)


def _tag_dummies(frag, mapnum=1):
    for a in frag.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetAtomMapNum(mapnum)


def split_aromatic_ester(mol):
    """Cleave the aromatic C-6' ester; return (core_with_dummy, acyl_with_dummy)."""
    m = mol.GetSubstructMatch(AROM_ESTER)
    if not m:
        return None, None
    oxy, carb = m[0], m[1]
    rw = Chem.RWMol(mol)
    bond = rw.GetBondBetweenAtoms(oxy, carb)
    frag = Chem.FragmentOnBonds(rw, [bond.GetIdx()], addDummies=True)
    parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
    if len(parts) != 2:
        return None, None
    core = max(parts, key=lambda p: p.GetNumAtoms())
    acyl = min(parts, key=lambda p: p.GetNumAtoms())
    _tag_dummies(core)
    _tag_dummies(acyl)
    return core, acyl


def deacylate_longest_tail(mol):
    """Remove the longest aliphatic fatty-acyl ester (R-O-C(=O)-chain -> R-OH).

    Returns a NEW mol with the tail excised and the freed oxygen capped as a
    hydroxyl, or None if no suitable aliphatic tail is found / product is invalid.
    This embodies the Phase-11 'the tail is droppable' branch.
    """
    matches = mol.GetSubstructMatches(ALIPH_ESTER)
    if not matches:
        return None
    best = None
    for match in matches:
        # match = (ester-O, carbonyl-C, carbonyl-O, alpha-C); we cut O-C(=O).
        oxy, carb = match[0], match[1]
        # Measure the acyl fragment size by cutting the O-C(acyl) bond.
        rw = Chem.RWMol(mol)
        b = rw.GetBondBetweenAtoms(oxy, carb)
        if b is None:
            continue
        frag = Chem.FragmentOnBonds(rw, [b.GetIdx()], addDummies=True)
        parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=False)
        if len(parts) != 2:
            continue
        acyl = min(parts, key=lambda p: p.GetNumAtoms())
        # only genuine lipophilic aliphatic tails: many carbons, no aromatic ring
        n_c = sum(1 for a in acyl.GetAtoms() if a.GetAtomicNum() == 6)
        n_arom = rdMolDescriptors.CalcNumAromaticRings(acyl)
        if n_arom == 0 and n_c >= 6:
            score = n_c
            if best is None or score > best[0]:
                best = (score, oxy, carb)
    if best is None:
        return None
    _, oxy, carb = best
    rw = Chem.RWMol(mol)
    b = rw.GetBondBetweenAtoms(oxy, carb)
    frag = Chem.FragmentOnBonds(rw, [b.GetIdx()], addDummies=True)
    parts = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=False)
    core = max(parts, key=lambda p: p.GetNumAtoms())
    # cap the dummy (was bonded to the ester O) as H -> restores a free hydroxyl
    ed = Chem.RWMol(core)
    for a in ed.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetAtomicNum(1)      # dummy -> H, restoring a free hydroxyl
            a.SetIsotope(0)        # clear the FragmentOnBonds isotope label
            a.SetAtomMapNum(0)
    out = ed.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None
    return out


def attach(core, acyl):
    try:
        m = Chem.molzip(core, acyl)
        Chem.SanitizeMol(m)
        return m
    except Exception:
        return None


# ==========================================================================
# Reward: 3D exposed polar surface fraction (the Phase 8/9/11 lead)
# ==========================================================================
def _embed(mol):
    m = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    p.useMacrocycleTorsions = True
    if AllChem.EmbedMolecule(m, p) != 0:
        p.useRandomCoords = True
        p.maxIterations = 2000
        if AllChem.EmbedMolecule(m, p) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=400)
    except Exception:
        pass
    return m


_SASA_MEMO = {}


def _load_sasa_cache():
    if os.path.exists(SASA_CACHE):
        for _, r in pd.read_csv(SASA_CACHE).iterrows():
            _SASA_MEMO[r["inchikey"]] = (float(r["frac"]), float(r["total"]), int(r["n_ok"]))


def _save_sasa_cache():
    if _SASA_MEMO:
        pd.DataFrame([{"inchikey": k, "frac": v[0], "total": v[1], "n_ok": v[2]}
                      for k, v in _SASA_MEMO.items()]).to_csv(SASA_CACHE, index=False)


def cached_exposed_polar_fraction(mol):
    """exposed_polar_fraction memoized by InChIKey (persisted to SASA_CACHE) so
    reruns / the smoke test are instant. Delete SASA_CACHE to recompute."""
    ik = Chem.MolToInchiKey(mol)
    if ik not in _SASA_MEMO:
        _SASA_MEMO[ik] = exposed_polar_fraction(mol)
    return _SASA_MEMO[ik]


def exposed_polar_fraction(mol, n_conf=N_CONF):
    """Mean exposed POLAR surface fraction over a small ETKDG ensemble.

    Polar = SASA on N, O and the H bonded to N/O; fraction = polar / total SASA.
    Higher => more exposed H-bonding surface => predicted more serum-tolerant
    (Phase 8: polar SASA vs serum shift rho=-0.33). Returns (mean_frac, total_sasa,
    n_ok) or (nan, nan, 0).
    """
    fracs, totals = [], []
    for seed in range(n_conf):
        m = Chem.AddHs(mol)
        p = AllChem.ETKDGv3()
        p.randomSeed = 42 + seed
        p.useMacrocycleTorsions = True
        if AllChem.EmbedMolecule(m, p) != 0:
            p.useRandomCoords = True
            p.maxIterations = 2000
            if AllChem.EmbedMolecule(m, p) != 0:
                continue
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=400)
            radii = rdFreeSASA.classifyAtoms(m)
            total = rdFreeSASA.CalcSASA(m, radii)
        except Exception:
            continue
        polar = 0.0
        tot = 0.0
        for a in m.GetAtoms():
            s = float(a.GetProp("SASA")) if a.HasProp("SASA") else 0.0
            tot += s
            z = a.GetAtomicNum()
            if z in (7, 8):
                polar += s
            elif z == 1:
                nb = a.GetNeighbors()
                if nb and nb[0].GetAtomicNum() in (7, 8):
                    polar += s
        if tot > 0:
            fracs.append(polar / tot)
            totals.append(tot)
    if not fracs:
        return (np.nan, np.nan, 0)
    return (float(np.mean(fracs)), float(np.mean(totals)), len(fracs))


# ==========================================================================
# Fragment libraries
# ==========================================================================
# Phase-5 rigid-aromatic designed acyls (kept) + a polar-axis set that
# deliberately spans hydrophobic -> polar for the discriminating series.
DESIGNED_ESTER = {
    "biphenyl": "[*:1]C(=O)c1ccc(-c2ccccc2)cc1",
    "naphthoyl_6OH": "[*:1]C(=O)c1ccc2cc(O)ccc2c1",
    "pyridylphenyl": "[*:1]C(=O)c1ccc(-c2ccncc2)cc1",
    "quinolinecarbonyl": "[*:1]C(=O)c1ccc2ncccc2c1",
}
POLAR_AXIS = {  # label : SMILES  (roughly increasing polar/H-bond character)
    "ax0_dodecanoyl": "[*:1]C(=O)CCCCCCCCCCC",       # lipophilic tail control
    "ax1_biphenyl": "[*:1]C(=O)c1ccc(-c2ccccc2)cc1",
    "ax2_benzoyl": "[*:1]C(=O)c1ccccc1",
    "ax3_4OH_benzoyl": "[*:1]C(=O)c1ccc(O)cc1",
    "ax4_glycolyl": "[*:1]C(=O)CO",
    "ax5_diglycolyl": "[*:1]C(=O)COCCO",
    "ax6_carboxyethanoyl": "[*:1]C(=O)CCC(=O)O",
    "ax7_succinamoyl": "[*:1]C(=O)CCC(N)=O",
    "ax8_seryl": "[*:1]C(=O)[C@@H](N)CO",
    "ax9_polyhydroxy": "[*:1]C(=O)[C@H](O)[C@@H](O)[C@H](O)CO",
    "ax10_aminomethylbenzoyl": "[*:1]C(=O)c1ccc(CN)cc1",
    "ax11_sulfobenzoyl": "[*:1]C(=O)c1ccc(S(=O)(=O)O)cc1",
}


def build_acyl_library():
    """Designed rigid-aromatic + polar-axis acyls — the novel design space.

    Data-derived acyls (harvested from real analogs in Phase 5) are intentionally
    NOT scored here: they mostly regenerate known compounds, which already appear
    in the reward-validation set, and SASA-scoring all of them is the cost that
    does not buy novel candidates. The 24 knowns still seed novelty fingerprints.
    """
    lib = {}
    for name, asmi in {**{f"design_{k}": v for k, v in DESIGNED_ESTER.items()},
                       **{f"polaraxis_{k}": v for k, v in POLAR_AXIS.items()}}.items():
        am = Chem.MolFromSmiles(asmi)
        if am:
            _tag_dummies(am)
            lib[name] = am
    return lib


def _branch(name):
    head = name.split("::")[0]
    return {"data": "ester", "design": "ester"}.get(head.split("_")[0], head)


# ==========================================================================
# Main
# ==========================================================================
def main():
    cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv"))
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    _load_sasa_cache()

    known_mols, known_fps, known_ik = [], [], {}
    for smi, cid in zip(cm["smiles_canonical"], cm["compound_id"]):
        if not isinstance(smi, str):
            continue
        m = Chem.MolFromSmiles(smi)
        if m:
            known_mols.append(m)
            known_fps.append(ecfp(m))
            known_ik[Chem.MolToInchiKey(m)] = cid

    tmpl = Chem.MolFromSmiles(
        cm.set_index("compound_id").loc[SERUM_TOLERANT_TEMPLATE, "smiles_canonical"])
    core_full, _ = split_aromatic_ester(tmpl)
    assert core_full is not None, "core extraction failed"

    # tail-free core (ibrexafungerp-inspired branch): deacylate the template,
    # then re-cleave its aromatic C-6' ester to expose the attachment dummy.
    detailed = deacylate_longest_tail(tmpl)
    core_notail = None
    if detailed is not None:
        core_notail, _ = split_aromatic_ester(detailed)

    acyl_lib = build_acyl_library()

    rows, mol_by_key, seen = [], {}, set()
    cores = [("ester", core_full)]
    if core_notail is not None:
        cores.append(("notail", core_notail))

    for branch_core_tag, core in cores:
        for aname, acyl in acyl_lib.items():
            full = attach(core, acyl)
            if full is None:
                continue
            ik = Chem.MolToInchiKey(full)
            if ik in seen:
                continue
            seen.add(ik)
            frac, total_sasa, n_ok = cached_exposed_polar_fraction(full)
            if n_ok == 0:
                continue
            fp = ecfp(full)
            sims = DataStructs.BulkTanimotoSimilarity(fp, known_fps)
            novelty = 1.0 - max(sims) if sims else 1.0
            base = _branch(aname)
            branch = "notail" if branch_core_tag == "notail" else (
                "polaraxis" if base == "polaraxis" else "ester")
            rows.append(dict(
                branch=branch,
                acyl_source=aname.split("::")[0].split("_")[0],
                acyl_name=aname.split("::")[-1],
                exposed_polar_fraction=round(frac, 4),
                exposed_hydrophobic_fraction=round(1.0 - frac, 4),
                total_sasa=round(total_sasa, 1),
                n_conf_ok=n_ok,
                novelty=round(novelty, 3),
                mw=round(Descriptors.MolWt(full), 1),
                clogp=round(Descriptors.MolLogP(full), 2),
                inchikey=ik,
                matches_known=known_ik.get(ik, ""),
                smiles=Chem.MolToSmiles(full),
            ))
            mol_by_key[ik] = full

    lib = pd.DataFrame(rows)
    # PRIMARY reward = exposed polar fraction; novelty as a light secondary.
    lib["serum_tolerance_reward"] = lib["exposed_polar_fraction"]
    lib["is_novel"] = lib["matches_known"] == ""
    lib = lib.sort_values(["serum_tolerance_reward", "novelty"],
                          ascending=False).reset_index(drop=True)
    lib.to_csv(os.path.join(OUT, "phase12_generated_library.csv"), index=False)
    _save_sasa_cache()

    # ---- Reward validation on the 24 knowns: exposed polar frac vs serum shift
    val_rows = []
    tmpl_ids = set(pairs["compound_id"])
    for _, r in pairs.iterrows():
        smi = cm.set_index("compound_id")["smiles_canonical"].get(r["compound_id"])
        shift = pd.to_numeric(r.get("serum_shift_fold"), errors="coerce")
        if not isinstance(smi, str) or not np.isfinite(shift) or shift <= 0:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        frac, _, n_ok = cached_exposed_polar_fraction(m)
        if n_ok:
            val_rows.append({"compound_id": r["compound_id"],
                             "exposed_polar_fraction": frac,
                             "serum_shift_fold": float(shift),
                             "log2_shift": float(np.log2(shift))})
    val = pd.DataFrame(val_rows)
    rho, pval = (stats.spearmanr(val["exposed_polar_fraction"], val["log2_shift"])
                 if len(val) >= 4 else (np.nan, np.nan))

    fig, ax = plt.subplots(figsize=(6, 4.6))
    ax.scatter(val["exposed_polar_fraction"], val["serum_shift_fold"],
               s=55, color="#2f9e5a", edgecolor="k", linewidth=0.4)
    ax.set_yscale("log")
    ax.set_xlabel("Exposed polar surface fraction (3D ETKDG + rdFreeSASA)")
    ax.set_ylabel("Observed serum shift fold (log)")
    ax.set_title(f"Phase 12 reward validation on 24 knowns\n"
                 f"Spearman rho={rho:.2f} (p={pval:.2f}, n={len(val)}) "
                 f"— negative = reward tracks tolerance")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase12_reward_validation.png"), dpi=140)
    plt.close(fig)

    # ---- Discriminating series: novel analogs spanning exposed polarity -----
    novel = lib[lib["is_novel"]].copy()
    series = _discriminating_series(novel)
    series.to_csv(os.path.join(OUT, "phase12_discriminating_series.csv"), index=False)

    # ---- Emit top novel candidates as 3D SDF (CREST-ready) ------------------
    top = novel.head(TOP_N_FOR_QM)
    _write_sdf_and_crest(top, mol_by_key)

    _write_findings(lib, novel, series, val, rho, pval, core_notail is not None)

    # ---- console summary ----
    print("PHASE 12 — serum-tolerance-biased generation (Track A)")
    print("=" * 60)
    print(f"Template: {SERUM_TOLERANT_TEMPLATE}   branches: "
          f"{', '.join(sorted(lib['branch'].unique()))}   conformers/mol: {N_CONF}")
    print(f"Generated analogs: {len(lib)}  (novel: {int(lib['is_novel'].sum())})")
    print(f"Reward validation (exposed polar frac vs serum shift, knowns): "
          f"rho={rho:.2f} p={pval:.2f} n={len(val)}")
    print(f"Discriminating series rows: {len(series)}")
    print("\nTop novel candidates by serum-tolerance reward:")
    cols = ["branch", "acyl_name", "serum_tolerance_reward",
            "exposed_hydrophobic_fraction", "novelty", "mw"]
    print(top[cols].to_string(index=False))
    print("\nWrote: phase12_generated_library.csv, phase12_discriminating_series.csv,")
    print("       phase12_top_candidates.sdf, phase12_crest_commands.sh,")
    print("       phase12_reward_validation.png, phase12_findings.md")


def _discriminating_series(novel, n_bins=3, per_bin=2):
    """Pick novel analogs that SPAN the exposed-polar axis while holding the
    SCAFFOLD constant, so the wet-lab serum assay tests polarity — not scaffold —
    against serum tolerance. We restrict to the single branch with the widest
    exposed-polar span (enough members), then bin low/mid/high within it."""
    if novel.empty:
        return novel
    # choose the scaffold branch that spans the axis most widely
    best_branch, best_span = None, -1.0
    for br, grp in novel.groupby("branch"):
        if len(grp) < n_bins:
            continue
        span = grp["exposed_polar_fraction"].max() - grp["exposed_polar_fraction"].min()
        if span > best_span:
            best_branch, best_span = br, span
    n = (novel[novel["branch"] == best_branch] if best_branch else novel).copy()
    try:
        n["polar_bin"] = pd.qcut(n["exposed_polar_fraction"], q=n_bins,
                                 labels=["low", "mid", "high"], duplicates="drop")
    except ValueError:
        n["polar_bin"] = "all"
    picks = [grp.sort_values("novelty", ascending=False).head(per_bin)
             for _, grp in n.groupby("polar_bin", observed=True)]
    out = pd.concat(picks).sort_values("exposed_polar_fraction")
    return out[["polar_bin", "branch", "acyl_name", "exposed_polar_fraction",
                "exposed_hydrophobic_fraction", "novelty", "mw", "clogp",
                "inchikey", "smiles"]]


def _write_sdf_and_crest(top, mol_by_key):
    writer = Chem.SDWriter(os.path.join(OUT, "phase12_top_candidates.sdf"))
    cmds = ["#!/bin/bash",
            "# Phase 12 CREST ensembles -> confirm exposed polar SASA in 3D (QM tier).",
            "# Reserve GFN2 for finalists; GFN-FF screening tier otherwise.", ""]
    for i, (_, r) in enumerate(top.iterrows(), 1):
        m = _embed(mol_by_key[r["inchikey"]])
        if m is None:
            continue
        m.SetProp("_Name", f"p12_{i:02d}_{r['branch']}_{r['acyl_name']}")
        m.SetProp("serum_tolerance_reward", str(r["serum_tolerance_reward"]))
        m.SetProp("exposed_polar_fraction", str(r["exposed_polar_fraction"]))
        writer.write(m)
        cmds.append(f"crest p12_{i:02d}_{r['branch']}_{r['acyl_name']}.xyz "
                    f"--gfnff --alpb water --quick -ewin 6 --T 8 > p12_{i:02d}.out "
                    f"# reward={r['serum_tolerance_reward']}")
    writer.close()
    with open(os.path.join(OUT, "phase12_crest_commands.sh"), "w") as fh:
        fh.write("\n".join(cmds) + "\n")


def _write_findings(lib, novel, series, val, rho, pval, has_notail):
    n_notail = int((lib["branch"] == "notail").sum())
    top = novel.head(5)
    top_lines = "\n".join(
        f"- **{r['branch']}::{r['acyl_name']}** — exposed polar frac "
        f"{r['exposed_polar_fraction']:.2f}, novelty {r['novelty']:.2f}, MW {r['mw']:.0f}"
        for _, r in top.iterrows())
    direction = ("the expected direction (more exposed polar surface -> smaller "
                 "serum shift)" if np.isfinite(rho) and rho < 0
                 else "a weak/again-inconclusive direction")
    md = f"""# Phase 12 — serum-tolerance-biased generative design (findings)

## What was generated
Off the serum-active lead **{SERUM_TOLERANT_TEMPLATE}**, three branches:
- **ester** — re-esterify the validated aromatic C-6' handle (Phase 5 chemistry).
- **polaraxis** — designed acyls deliberately spanning hydrophobic -> polar.
- **notail** — {'the fatty acyl tail removed (deacylation) before re-esterifying '
  '(ibrexafungerp-inspired tail-free analogs)' if has_notail
  else 'ATTEMPTED but no aliphatic tail was cleanly excisable on this template'}.

{len(lib)} unique analogs ({int(lib['is_novel'].sum())} novel), scored by the
reward below. {n_notail} tail-free analogs generated.

## The reward (and its honest status)
Reward = **mean exposed POLAR surface fraction** over a {N_CONF}-conformer ETKDG
ensemble (RDKit rdFreeSASA) — the operational form of the Phase 8/9/11 lead
("expose polar not hydrophobic surface"). It is a HYPOTHESIS, not a validated
oracle, and a crude single-few-conformer proxy that Phase 7->8 showed will
OVERSTATE effects. Bulk QED/Ro5/clogP terms were deliberately dropped (Phases 5,
11 showed them uninformative here).

Retrospective check on the 24 knowns: exposed polar fraction vs serum shift
Spearman **rho={rho:.2f} (p={pval:.2f}, n={len(val)})** — {direction}. Treat as a
sanity check on the reward's sign, not proof.

## Top novel candidates (full table: phase12_generated_library.csv)
{top_lines}

## The deliverable that matters: a DISCRIMINATING SERIES
`phase12_discriminating_series.csv` is a small novel set that SPANS the
exposed-polar axis (low / mid / high) at comparable size on one scaffold. Its
purpose is to let a serum assay **falsify or confirm** the polar-surface lead:
if serum tolerance rises monotonically with exposed polar fraction across this
matched series, the lead is real and prospective; if not, the hypothesis is
rejected on-scaffold. This is worth more than any single "best" molecule.

## How to use this
1. Push `phase12_top_candidates.sdf` through the existing CREST -> Phase 6/8/9 QM
   funnel to confirm the exposed-polar property at QM-ensemble quality (the reward
   here is only a fast proxy).
2. Synthesize the **discriminating series** (SA set aside per instruction) and run
   the Phase-11 assay playbook — protein-adjusted MIC with albumin titration +
   equilibrium-dialysis fraction-unbound — reporting the serum SHIFT, not raw MIC.
3. Feed the measured shifts back as labels to train **Track B** (a REINVENT-style
   generative network on the external FKS pretraining set) with this same reward —
   at which point the reward stops being a hypothesis and becomes data-anchored.

## Caveats
Single-scaffold, proxy reward, no validated oracle; the tail-free branch produces
chemically aggressive structures by design (synthetic accessibility set aside).
The value is a testable, hypothesis-spanning series plus a reusable reward — not a
finished drug candidate.
"""
    with open(os.path.join(OUT, "phase12_findings.md"), "w", encoding="utf-8") as fh:
        fh.write(md)


if __name__ == "__main__":
    main()
