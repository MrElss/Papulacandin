#!/usr/bin/env python3
"""
phase7_retrospective_qm.py
==========================
PHASE 7 — Retrospective test of the Phase-6 3D hypothesis on REAL serum data.

THE QUESTION
------------
Phase 5 showed the flat 2D design-rule score only weakly tracks observed serum
MIC across the known fusacandin C-6' ester series (Spearman rho = 0.32, n.s.).
Phase 6 proposed the discriminating signal is 3D/conformational: the
Boltzmann-weighted exposed HYDROPHOBIC surface (SASA) and its ensemble SPREAD
(a true rigidity metric), not 2D topology.

This module tests that proposal head-to-head, on the SAME compounds and the
SAME endpoint the 0.32 baseline used: the Yeung-1996 fusacandin analog series
(6a-6u) plus the reference natural products, all of which carry a matched
serum-free / serum whole-cell MIC against C. albicans.

METHOD (and its honest limits)
------------------------------
The 12 Phase-5 *novel candidates* were profiled from cluster CREST/GFN-FF
ensembles. The KNOWN compounds here do NOT have CREST ensembles yet (that is
the rigorous follow-up this script is meant to justify). To get a fast,
self-consistent 3D readout NOW, each known compound's ensemble is built with
RDKit ETKDGv3 + MMFF94 (the same machinery phase6's self-test uses), written
in CREST's exact xyz format, and pushed through the IDENTICAL phase6 SASA /
shape descriptor pipeline. Because every compound is treated by one method on
one shared scaffold, RANK correlations are internally valid even though the
absolute SASA values are MMFF- not GFN-FF-derived. A positive result here is
the green light to spend cluster time CREST-ing the known set for confirmation.

OUTPUTS
-------
  outputs/phase7_known_qm_descriptors.csv  — 3D descriptors per known compound
  outputs/phase7_retrospective_qm.csv      — descriptors + serum MIC + stats
  outputs/phase7_retrospective_qm.png      — 3D descriptor vs serum MIC, with
                                             the 2D baseline (rho=0.32) for scale
"""
from __future__ import annotations
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem

import phase6_qm_layer as p6

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
SDF_DIR = os.path.join(HERE, "..", "curated", "structures", "sdf")
N_CONFS = 40
RANDOM_SEED = 0xC0FFEE
BASELINE_RHO = 0.32   # Phase-5 2D design-rule score vs serum MIC (n.s.)


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


def find_sdf(compound_id):
    hits = glob.glob(os.path.join(SDF_DIR, f"{compound_id}_*.sdf"))
    return hits[0] if hits else None


def build_ensemble_xyz(sdf_path, name, out_xyz, n_confs=N_CONFS):
    """RDKit ETKDGv3 + MMFF94 ensemble, written in CREST crest_conformers.xyz
    format (energy on the comment line, Hartree, lowest first) so the phase6
    parser/descriptor code consumes it unchanged. Returns (mol_with_Hs, n)."""
    mol = Chem.SDMolSupplier(sdf_path, removeHs=False)[0]
    if mol is None:
        return None, 0
    mol = Chem.AddHs(mol, addCoords=False)
    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    params.useRandomCoords = True
    params.pruneRmsThresh = 0.5
    params.numThreads = 0
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not cids:                                   # large flexible bRo5: retry harder
        params.maxIterations = 2000
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not cids:
        return None, 0
    energies = []
    keep = []
    props = AllChem.MMFFGetMoleculeProperties(mol)
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid) if props else None
        if ff is None:
            continue
        ff.Minimize(maxIts=1000)
        energies.append(ff.CalcEnergy())
        keep.append(cid)
    if not keep:
        return None, 0
    order = np.argsort(energies)
    with open(out_xyz, "w") as fh:
        for rank in order:
            cid = keep[int(rank)]
            conf = mol.GetConformer(cid)
            fake_hartree = -2000.0 + energies[int(rank)] / p6.KCAL_PER_HARTREE
            fh.write(f"{mol.GetNumAtoms()}\n{fake_hartree:.8f}\n")
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                fh.write(f"{atom.GetSymbol():2s} {pos.x:14.8f} {pos.y:14.8f} {pos.z:14.8f}\n")
    return mol, len(keep)


def compute_known_descriptors():
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    pairs["serum_mic_num"] = pairs["serum_mic_ugml"].apply(parse_mic)
    ens_dir = os.path.join(OUT, "qm_runs_known_rdkit")
    os.makedirs(ens_dir, exist_ok=True)
    rows = []
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        sdf = find_sdf(cid)
        if sdf is None:
            print(f"[skip] no SDF for {cid}")
            continue
        out_xyz = os.path.join(ens_dir, f"{cid}_ensemble.xyz")
        mol, n = build_ensemble_xyz(sdf, cid, out_xyz)
        if mol is None or n < 4:
            print(f"[skip] ensemble failed for {cid} ({r['name']})")
            continue
        polar_mask, elements_ref = p6.classify_polarity(mol)
        frames = p6.parse_crest_ensemble(out_xyz)
        desc = p6.ensemble_descriptors(frames, polar_mask, elements_ref, cid)
        desc.pop("weights", None)
        desc.pop("rel_kcal", None)
        desc["compound_id"] = cid
        desc["name"] = r["name"]
        rows.append(desc)
        print(f"[ok] {cid:10s} {r['name']:24s} n={n:3d} "
              f"hSASA={desc['hydrophobic_sasa_mean']:6.1f} "
              f"hSASA_sd={desc['hydrophobic_sasa_std']:4.1f} "
              f"hfrac={desc['hydrophobic_fraction_mean']:.3f}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase7_known_qm_descriptors.csv"), index=False)
    return df.merge(pairs, on=["compound_id", "name"], how="left")


DESCRIPTORS = [
    ("hydrophobic_sasa_mean", "Boltzmann-wtd hydrophobic SASA (A^2)"),
    ("hydrophobic_fraction_mean", "Hydrophobic fraction of total SASA"),
    ("hydrophobic_sasa_std", "Hydrophobic SASA ensemble spread (rigidity)"),
    ("polar_sasa_mean", "Boltzmann-wtd polar SASA (A^2)"),
    ("rg_mean", "Radius of gyration (A)"),
]


def correlate_and_plot(m):
    y = np.log10(m["serum_mic_num"].values)
    print("\n" + "=" * 68)
    print(f"RETROSPECTIVE: 3D descriptor vs log10(serum MIC), n={len(m)}")
    print(f"2D baseline to beat: Spearman rho = {BASELINE_RHO:+.2f} (n.s.)")
    print("(positive rho = descriptor rises as the compound gets more"
          " serum-INACTIVATED, i.e. the descriptor flags serum liability)")
    print("=" * 68)
    results = []
    for col, label in DESCRIPTORS:
        rho, p = stats.spearmanr(m[col].values, y)
        results.append((col, label, rho, p))
        flag = "  <-- beats 2D baseline" if abs(rho) > abs(BASELINE_RHO) else ""
        print(f"  {col:28s} rho={rho:+.2f}  p={p:.3f}{flag}")

    # serum-active vs serum-killed group separation (Mann-Whitney)
    print("\nGroup separation (serum_active yes vs no), Mann-Whitney U:")
    grp = m.groupby("serum_active")
    for col, label in DESCRIPTORS:
        if {"yes", "no"} <= set(grp.groups):
            a = grp.get_group("yes")[col].values
            b = grp.get_group("no")[col].values
            u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
            print(f"  {col:28s} tolerant={a.mean():6.1f}  killed={b.mean():6.1f}"
                  f"  p={pu:.3f}")

    # best descriptor figure, with the 2D baseline panel for scale
    best = max(results, key=lambda t: abs(t[2]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    colors = m["serum_active"].map({"yes": "#2a9d8f", "no": "#e76f51"})
    ax.scatter(m[best[0]], m["serum_mic_num"], c=colors, s=60,
               edgecolor="k", linewidth=0.4)
    ax.set_yscale("log")
    ax.set_xlabel(best[1])
    ax.set_ylabel("Observed serum MIC (ug/mL, log)")
    ax.set_title(f"Phase 7 (3D, RDKit ensemble): best descriptor\n"
                 f"{best[0]}  Spearman rho={best[2]:+.2f} (p={best[3]:.3f}, n={len(m)})")
    for _, r in m.iterrows():
        ax.annotate(r["name"].replace("Fusacandin analog ", ""),
                    (r[best[0]], r["serum_mic_num"]), fontsize=6,
                    xytext=(2, 2), textcoords="offset points")

    ax = axes[1]
    ax.bar([d[0].replace("_mean", "").replace("hydrophobic", "hphob")
            for d in results], [abs(d[2]) for d in results],
           color="#264653")
    ax.axhline(abs(BASELINE_RHO), color="#e76f51", ls="--",
               label=f"2D baseline |rho|={abs(BASELINE_RHO):.2f}")
    ax.set_ylabel("|Spearman rho| vs serum MIC")
    ax.set_title("3D descriptors vs the 2D baseline")
    ax.tick_params(axis="x", rotation=35, labelsize=7)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase7_retrospective_qm.png"), dpi=140)
    plt.close(fig)

    out = m.copy()
    out.to_csv(os.path.join(OUT, "phase7_retrospective_qm.csv"), index=False)
    return results, best


if __name__ == "__main__":
    print("Building RDKit ensembles + 3D descriptors for known serum-gap set")
    print("(MMFF94 proxy for CREST; rank-correlation is the readout) ...\n")
    m = compute_known_descriptors()
    results, best = correlate_and_plot(m)
    print("\n" + "=" * 68)
    print(f"VERDICT: best 3D descriptor = {best[0]} (|rho|={abs(best[2]):.2f}) "
          f"vs 2D baseline |rho|={abs(BASELINE_RHO):.2f}")
    print("Wrote outputs/phase7_known_qm_descriptors.csv, "
          "phase7_retrospective_qm.csv, phase7_retrospective_qm.png")
    print("=" * 68)
