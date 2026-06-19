#!/usr/bin/env python3
"""
phase9_electronic.py  (Phase 9 — electronic descriptors, analysis stage)
========================================================================
Parse the GFN2-xTB single points produced by run_xtb_electronic.sbatch and test
whether ELECTRONIC / solvation descriptors predict serum tolerance where shape
did not (Phase 8). Mirrors the Phase-8 statistics exactly so results are
directly comparable.

Per conformer it reads <conf>_water.out and <conf>_octanol.out and extracts:
  * dipole moment (Debye)                  -- overall polarity
  * HOMO-LUMO gap (eV)                      -- electronic softness/hardness
  * molecular polarizability alpha(0) (au)  -- dispersion / hydrophobic binding capacity
  * G_solv(water), G_solv(octanol) (Eh)     -- solvation free energies
  * logP_xtb = (Gsolv_water - Gsolv_oct)/(2.303 RT)  -- QM hydrophobicity
Then Boltzmann-averages each over the compound's populated conformers (weights
from weights.csv) and correlates vs serum data:
  - Spearman vs log10(serum MIC)
  - partial Spearman controlling for serum-free MIC (the Phase-8 confound)
  - Spearman vs the pure serum SHIFT fold

OUTPUTS
  outputs/phase9_electronic_descriptors.csv
  outputs/phase9_electronic_stats.csv
  outputs/phase9_electronic.png
"""
from __future__ import annotations
import os
import re
import glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
XTB_ROOT = os.path.join(OUT, "qm_runs_known_xtb")
HARTREE_KCAL = 627.5094740631
RT_KCAL = 1.987204e-3 * 298.15           # 2.303*RT denominator below
LOGP_DENOM = 2.302585 * RT_KCAL          # kcal/mol


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


def _search_float(text, pattern, group=1):
    m = re.search(pattern, text)
    return float(m.group(group)) if m else np.nan


def parse_xtb_out(path):
    """Extract dipole (Debye), HOMO-LUMO gap (eV), polarizability alpha(0) (au),
    and G_solv (Eh) from a GFN2-xTB --alpb output. Robust to format drift via
    targeted regexes; returns NaN for anything missing."""
    if not os.path.exists(path):
        return {}
    txt = open(path, errors="ignore").read()
    d = {}
    # dipole: the "full:" line under "molecular dipole:", last column = tot (Debye).
    # non-greedy .*? skips the header + "q only:" lines robustly.
    m = re.search(r"molecular dipole:.*?full:\s+"
                  r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt, re.S)
    d["dipole_D"] = float(m.group(4)) if m else np.nan
    d["homo_lumo_gap_eV"] = _search_float(txt, r"HOMO-LUMO\s+GAP\s+([-\d.]+)\s+eV")
    if np.isnan(d["homo_lumo_gap_eV"]):
        d["homo_lumo_gap_eV"] = _search_float(txt, r"HL-Gap\s+([-\d.]+)\s*Eh")  # fallback (Eh)
    d["alpha_au"] = _search_float(txt, r"Mol\.\s*α\(0\)\s*/au\s*:?\s*([-\d.]+)")
    if np.isnan(d["alpha_au"]):
        d["alpha_au"] = _search_float(txt, r"Mol\.\s*alpha\(0\).*?([-\d.]+)")
    d["gsolv_Eh"] = _search_float(txt, r"->\s*Gsolv\s+([-\d.]+)\s+Eh")
    return d


def boltzmann_avg(values, weights):
    v = np.array(values, float)
    w = np.array(weights, float)
    ok = ~np.isnan(v)
    if not ok.any():
        return np.nan
    w = w[ok] / w[ok].sum()
    return float((v[ok] * w).sum())


def compute():
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    pairs["serum_mic_num"] = pairs["serum_mic_ugml"].apply(parse_mic)
    pairs["serumfree_mic_num"] = pairs["serumfree_mic_ugml"].apply(parse_mic)
    rows = []
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        cdir = os.path.join(XTB_ROOT, cid)
        wfile = os.path.join(cdir, "weights.csv")
        if not os.path.exists(wfile):
            print(f"[skip] no xtb inputs for {cid}")
            continue
        wdf = pd.read_csv(wfile)
        dip, gap, alpha, gsolv_w, gsolv_o, logp, wts = [], [], [], [], [], [], []
        for _, cr in wdf.iterrows():
            conf = cr["conformer"]
            w = parse_xtb_out(os.path.join(cdir, f"{conf}_water.out"))
            o = parse_xtb_out(os.path.join(cdir, f"{conf}_octanol.out"))
            if not w:
                continue
            dip.append(w.get("dipole_D", np.nan))
            gap.append(w.get("homo_lumo_gap_eV", np.nan))
            alpha.append(w.get("alpha_au", np.nan))
            gw, go = w.get("gsolv_Eh", np.nan), o.get("gsolv_Eh", np.nan)
            gsolv_w.append(gw)
            gsolv_o.append(go)
            logp.append((gw - go) * HARTREE_KCAL / LOGP_DENOM
                        if not (np.isnan(gw) or np.isnan(go)) else np.nan)
            wts.append(cr["pop"])
        if not wts:
            print(f"[skip] no parsed xtb outputs for {cid}")
            continue
        rows.append(dict(
            compound_id=cid, name=r["name"], n_conf_used=len(wts),
            dipole_D=boltzmann_avg(dip, wts),
            homo_lumo_gap_eV=boltzmann_avg(gap, wts),
            alpha_au=boltzmann_avg(alpha, wts),
            gsolv_water_kcal=boltzmann_avg(gsolv_w, wts) * HARTREE_KCAL,
            logP_xtb=boltzmann_avg(logp, wts),
        ))
        print(f"[ok] {cid:10s} dipole={rows[-1]['dipole_D']:.2f}D "
              f"gap={rows[-1]['homo_lumo_gap_eV']:.2f}eV "
              f"logP={rows[-1]['logP_xtb']:.2f} "
              f"Gsolv(w)={rows[-1]['gsolv_water_kcal']:.1f}kcal")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase9_electronic_descriptors.csv"), index=False)
    return df.merge(pairs, on=["compound_id", "name"], how="left")


DESCRIPTORS = [
    ("dipole_D", "Molecular dipole (Debye)"),
    ("logP_xtb", "QM logP (octanol/water solvation)"),
    ("gsolv_water_kcal", "Aqueous solvation free energy (kcal/mol)"),
    ("alpha_au", "Molecular polarizability alpha(0) (au)"),
    ("homo_lumo_gap_eV", "HOMO-LUMO gap (eV)"),
]


def partial_spearman(x, y, z):
    from scipy.stats import rankdata
    xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
    rx = xr - np.polyval(np.polyfit(zr, xr, 1), zr)
    ry = yr - np.polyval(np.polyfit(zr, yr, 1), zr)
    return stats.pearsonr(rx, ry)


def main():
    print("PHASE 9 — GFN2-xTB electronic/solvation descriptors vs serum data\n")
    m = compute()
    if not len(m):
        print("No parsed data yet. Run run_xtb_electronic.sbatch first.")
        return
    sm = np.log10(m["serum_mic_num"].values)
    sf = np.log10(m["serumfree_mic_num"].values)
    shift = np.log10(m["serum_shift_fold"].values)
    print(f"\n{'descriptor':22s} {'rho(serumMIC)':>14s} {'partial|pot':>12s} "
          f"{'rho(shift)':>11s}")
    out = []
    for col, _ in DESCRIPTORS:
        x = m[col].values
        r_mic, p_mic = stats.spearmanr(x, sm)
        r_par, p_par = partial_spearman(x, sm, sf)
        r_sh, p_sh = stats.spearmanr(x, shift)
        out.append(dict(descriptor=col, rho_serumMIC=round(r_mic, 2), p_serumMIC=round(p_mic, 3),
                        partial_rho=round(r_par, 2), partial_p=round(p_par, 3),
                        rho_shift=round(r_sh, 2), p_shift=round(p_sh, 3)))
        print(f"{col:22s} {r_mic:+.2f} (p={p_mic:.2f})  {r_par:+.2f} (p={p_par:.2f})  "
              f"{r_sh:+.2f} (p={p_sh:.2f})")
    res = pd.DataFrame(out)
    res.to_csv(os.path.join(OUT, "phase9_electronic_stats.csv"), index=False)

    # figure: |rho vs shift| per descriptor (the potency-free endpoint), with the
    # Phase-8 best shape lead (polar SASA, |rho|=0.33) for scale
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    labels = [d[0] for d in DESCRIPTORS]
    ax.bar(labels, [abs(o["rho_shift"]) for o in out], color="#2a6f97")
    ax.axhline(0.33, color="#e76f51", ls="--", label="Phase-8 best shape lead (polar SASA |rho|=0.33)")
    ax.set_ylabel("|Spearman rho| vs serum SHIFT")
    ax.set_title(f"Phase 9 electronic descriptors vs serum tolerance (n={len(m)})")
    ax.tick_params(axis="x", rotation=25, labelsize=8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase9_electronic.png"), dpi=140)
    plt.close(fig)
    best = max(out, key=lambda o: abs(o["rho_shift"]))
    print(f"\nBest electronic descriptor vs serum shift: {best['descriptor']} "
          f"(rho={best['rho_shift']:+.2f}, p={best['p_shift']:.3f})")
    print("Wrote phase9_electronic_descriptors.csv, phase9_electronic_stats.csv, phase9_electronic.png")


if __name__ == "__main__":
    main()
