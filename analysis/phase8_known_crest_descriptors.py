#!/usr/bin/env python3
"""
phase8_known_crest_descriptors.py
=================================
PHASE 8 — confirm (or revise) the Phase-7 retrospective with REAL CREST/GFN-FF
ensembles instead of the MMFF proxy.

Phase 7 used fast RDKit-MMFF ensembles for the 24 known serum-gap compounds and
found one 3D descriptor that beat the 2D baseline: Boltzmann-weighted hydrophobic
SASA, Spearman rho = -0.45 (p=0.029) vs log10(serum MIC), where the 2D
design-rule score managed only rho=0.32 (n.s.). The sign was the surprise:
*more* exposed hydrophobic surface tracked *better* serum activity, i.e. aromatic
extent/size (confounded with intrinsic potency), not albumin avoidance; the
rigidity/SASA-spread idea did not validate.

Phase 8 re-runs the IDENTICAL phase6 descriptor engine on the real CREST GFN-FF
conformer ensembles (analysis/outputs/qm_runs_known/PAPU-*/crest_conformers.xyz)
and re-tests the same correlations, so the MMFF-proxy conclusion can be checked
at proper QM-ensemble quality on the same compounds and endpoint.

Note on PAPU-0078 (6h): its CREST run (the --noreftopo retry) yielded a single
conformer (search did not sample; structure verified intact). It carries no
ensemble/flexibility information, so the PRIMARY analysis excludes it (n=23);
a secondary pass reports n=24 with 6h as a single-point for transparency.

OUTPUTS
  outputs/phase8_known_crest_descriptors.csv  — real-CREST 3D descriptors
  outputs/phase8_retrospective_crest.csv      — descriptors + serum MIC + flags
  outputs/phase8_retrospective_crest.png      — real CREST vs MMFF proxy vs 2D
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

import phase6_qm_layer as p6

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
SDF_DIR = os.path.join(HERE, "..", "curated", "structures", "sdf")
KNOWN_ROOT = os.path.join(OUT, "qm_runs_known")
BASELINE_RHO = 0.32     # Phase-5 2D design-rule score vs serum MIC (n.s.)
PROXY_RHO = -0.45       # Phase-7 MMFF-proxy hydrophobic_sasa_mean
MIN_CONF = 4            # below this an ensemble is treated as degenerate

DESCRIPTORS = [
    ("hydrophobic_sasa_mean", "Boltzmann-wtd hydrophobic SASA (A^2)"),
    ("hydrophobic_fraction_mean", "Hydrophobic fraction of total SASA"),
    ("hydrophobic_sasa_std", "Hydrophobic SASA ensemble spread (rigidity)"),
    ("polar_sasa_mean", "Boltzmann-wtd polar SASA (A^2)"),
    ("rg_mean", "Radius of gyration (A)"),
]


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


def find_sdf(cid):
    hits = glob.glob(os.path.join(SDF_DIR, f"{cid}_*.sdf"))
    return hits[0] if hits else None


def folder_for(cid):
    hits = glob.glob(os.path.join(KNOWN_ROOT, f"{cid}_*"))
    return hits[0] if hits else None


DESC_CACHE = "phase8_known_crest_descriptors.csv"


def compute(force=False):
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    pairs["serum_mic_num"] = pairs["serum_mic_ugml"].apply(parse_mic)
    pairs["serumfree_mic_num"] = pairs["serumfree_mic_ugml"].apply(parse_mic)
    cache = os.path.join(OUT, DESC_CACHE)
    if os.path.exists(cache) and not force:
        print(f"[cache] loading descriptors from {DESC_CACHE} "
              "(delete it or pass force=True to recompute SASA)")
        df = pd.read_csv(cache)
        return df.merge(pairs, on=["compound_id", "name"], how="left")
    rows = []
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        sdf, folder = find_sdf(cid), folder_for(cid)
        if not sdf or not folder:
            print(f"[skip] missing SDF/folder for {cid}")
            continue
        xyz = os.path.join(folder, "crest_conformers.xyz")
        if not os.path.exists(xyz):
            print(f"[skip] no CREST ensemble for {cid}")
            continue
        mol = Chem.AddHs(Chem.SDMolSupplier(sdf, removeHs=False)[0], addCoords=False)
        polar_mask, elements_ref = p6.classify_polarity(mol)
        frames = p6.parse_crest_ensemble(xyz)
        if frames[0]["elements"] != elements_ref:
            print(f"[warn] atom-order mismatch for {cid} "
                  f"({len(frames[0]['elements'])} vs {len(elements_ref)}); skipping")
            continue
        desc = p6.ensemble_descriptors(frames, polar_mask, elements_ref, cid)
        desc.pop("weights", None)
        desc.pop("rel_kcal", None)
        desc["compound_id"] = cid
        desc["name"] = r["name"]
        desc["degenerate"] = desc["n_conformers"] < MIN_CONF
        rows.append(desc)
        tag = "  *DEGENERATE (single conf)" if desc["degenerate"] else ""
        print(f"[ok] {cid:10s} {r['name']:24s} nconf={desc['n_conformers']:4d} "
              f"hSASA={desc['hydrophobic_sasa_mean']:6.1f} "
              f"hSASA_sd={desc['hydrophobic_sasa_std']:4.1f} "
              f"hfrac={desc['hydrophobic_fraction_mean']:.3f}{tag}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase8_known_crest_descriptors.csv"), index=False)
    return df.merge(pairs, on=["compound_id", "name"], how="left")


def spearman_block(m, label):
    y = np.log10(m["serum_mic_num"].values)
    print(f"\n--- {label} (n={len(m)}) ---")
    res = []
    for col, lab in DESCRIPTORS:
        rho, p = stats.spearmanr(m[col].values, y)
        res.append((col, lab, rho, p))
        flag = "  <-- beats 2D |0.32|" if abs(rho) > abs(BASELINE_RHO) else ""
        print(f"  {col:28s} rho={rho:+.2f}  p={p:.3f}{flag}")
    return res


def partial_spearman(x, y, z):
    """Spearman of x vs y controlling for z (rank-residual partial correlation)."""
    from scipy.stats import rankdata
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    rx = xr - np.polyval(np.polyfit(zr, xr, 1), zr)
    ry = yr - np.polyval(np.polyfit(zr, yr, 1), zr)
    return stats.pearsonr(rx, ry)


def confound_analysis(clean):
    """The decisive test: is the descriptor->serum-MIC link just intrinsic
    potency? (1) how strongly serum-free MIC alone tracks serum MIC; (2) does
    any descriptor survive controlling for serum-free MIC; (3) does any
    descriptor track the pure serum SHIFT (potency-independent)."""
    sf = clean["serumfree_mic_num"].values
    sm = np.log10(clean["serum_mic_num"].values)
    shift = np.log10(clean["serum_shift_fold"].values)
    rho_sf, p_sf = stats.spearmanr(sf, clean["serum_mic_num"].values)
    print("\n" + "=" * 70)
    print("CONFOUND CHECK — is the 3D signal just intrinsic potency?")
    print("=" * 70)
    print(f"serum-free MIC vs serum MIC: Spearman rho={rho_sf:+.2f} (p={p_sf:.3f})"
          "  <- potency dominates serum MIC")
    print("\nPartial Spearman (descriptor vs serum MIC | controlling serum-free MIC):")
    print("Direct Spearman (descriptor vs serum SHIFT fold | potency-independent):")
    out = []
    for col, _ in DESCRIPTORS:
        pr, pp = partial_spearman(clean[col].values, sm, np.log10(sf))
        sr, sp = stats.spearmanr(clean[col].values, shift)
        out.append(dict(descriptor=col, partial_rho_vs_serumMIC=round(pr, 2),
                        partial_p=round(pp, 3), rho_vs_shift=round(sr, 2),
                        shift_p=round(sp, 3)))
        print(f"  {col:28s} partial={pr:+.2f} (p={pp:.3f}) | shift={sr:+.2f} (p={sp:.3f})")
    pd.DataFrame(out).to_csv(os.path.join(OUT, "phase8_confound_analysis.csv"), index=False)
    return rho_sf, out


def main():
    print("PHASE 8 — real CREST/GFN-FF descriptors for the known serum-gap set\n")
    m = compute()
    m.to_csv(os.path.join(OUT, "phase8_retrospective_crest.csv"), index=False)

    clean = m[~m["degenerate"]].copy()
    n_deg = int(m["degenerate"].sum())
    primary_label = ("PRIMARY: all ensembles (none degenerate)" if n_deg == 0
                     else f"PRIMARY: healthy ensembles only (excl. {n_deg} degenerate)")
    print("\n" + "=" * 70)
    print("RETROSPECTIVE vs log10(serum MIC) — real CREST ensembles")
    print(f"benchmarks: 2D baseline rho={BASELINE_RHO:+.2f} (n.s.); "
          f"Phase-7 MMFF proxy hydrophobic_sasa_mean rho={PROXY_RHO:+.2f}")
    print("=" * 70)
    res_primary = spearman_block(clean, primary_label)
    res_all = spearman_block(m, "SECONDARY: all compounds (6h as single-point)")
    confound_analysis(clean)

    # group separation on primary set
    print("\nGroup separation (serum_active yes vs no), Mann-Whitney U [primary]:")
    grp = clean.groupby("serum_active")
    if {"yes", "no"} <= set(grp.groups):
        for col, _ in DESCRIPTORS:
            a = grp.get_group("yes")[col].values
            b = grp.get_group("no")[col].values
            u, pu = stats.mannwhitneyu(a, b, alternative="two-sided")
            print(f"  {col:28s} tolerant={a.mean():6.1f}  killed={b.mean():6.1f}  p={pu:.3f}")

    # figure: best CREST descriptor vs MIC, + bar of |rho| (CREST vs proxy vs 2D)
    best = max(res_primary, key=lambda t: abs(t[2]))
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7))
    ax = axes[0]
    colors = clean["serum_active"].map({"yes": "#2a9d8f", "no": "#e76f51"})
    ax.scatter(clean[best[0]], clean["serum_mic_num"], c=colors, s=60,
               edgecolor="k", linewidth=0.4)
    ax.set_yscale("log")
    ax.set_xlabel(best[1])
    ax.set_ylabel("Observed serum MIC (ug/mL, log)")
    ax.set_title(f"Phase 8 (real CREST/GFN-FF): best descriptor\n"
                 f"{best[0]}  rho={best[2]:+.2f} (p={best[3]:.3f}, n={len(clean)})")
    for _, r in clean.iterrows():
        ax.annotate(str(r["name"]).replace("Fusacandin analog ", ""),
                    (r[best[0]], r["serum_mic_num"]), fontsize=6,
                    xytext=(2, 2), textcoords="offset points")
    ax = axes[1]
    labels = [d[0].replace("_mean", "").replace("hydrophobic", "hphob") for d in res_primary]
    ax.bar(labels, [abs(d[2]) for d in res_primary], color="#264653")
    ax.axhline(abs(BASELINE_RHO), color="#e76f51", ls="--", label=f"2D baseline |rho|={abs(BASELINE_RHO):.2f}")
    ax.axhline(abs(PROXY_RHO), color="#e9c46a", ls=":", label=f"MMFF proxy |rho|={abs(PROXY_RHO):.2f}")
    ax.set_ylabel("|Spearman rho| vs serum MIC (real CREST)")
    ax.set_title("Real CREST descriptors vs MMFF proxy & 2D baseline")
    ax.tick_params(axis="x", rotation=35, labelsize=7)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase8_retrospective_crest.png"), dpi=140)
    plt.close(fig)

    print("\n" + "=" * 70)
    print(f"VERDICT: best real-CREST descriptor = {best[0]} (rho={best[2]:+.2f}); "
          f"MMFF proxy was {PROXY_RHO:+.2f}, 2D baseline {BASELINE_RHO:+.2f}")
    print("Wrote phase8_known_crest_descriptors.csv, phase8_retrospective_crest.csv,"
          " phase8_retrospective_crest.png")
    print("=" * 70)


if __name__ == "__main__":
    main()
