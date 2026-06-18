#!/usr/bin/env python3
"""
phase8b_finalist_gfn2.py
========================
Recompute the 3 funnel finalists' 3D descriptors from their GFN2-RERANKED
ensembles (crest --screen output, analysis/outputs/qm_runs/cand0*_gfn2/
crest_ensemble.xyz), and regenerate their Gaussian DFT inputs from the GFN2
Boltzmann populations -- so the expensive B3LYP single points run on the
conformers GFN2 (not GFN-FF) says are populated.

Compares each finalist's GFN-FF-tier descriptors (phase6_qm_descriptors.csv)
against the GFN2-reranked values, to show how much the better energies moved the
populations / SASA.

OUTPUTS
  outputs/phase8b_finalist_gfn2_descriptors.csv   — GFN2 descriptors (3 finalists)
  outputs/phase8b_gfnff_vs_gfn2.csv               — side-by-side delta
  outputs/qm_runs_gaussian_gfn2/<cand>/*.gjf      — DFT inputs from GFN2 pops
"""
from __future__ import annotations
import os
import glob
import pandas as pd
from rdkit import Chem
import phase6_qm_layer as p6

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
SDF = os.path.join(OUT, "phase5_top_candidates.sdf")
GFN2_GAUSS = os.path.join(OUT, "qm_runs_gaussian_gfn2")

FINALISTS = ["cand01_quinolinecarbonyl", "cand02_naphthoyl_6OH", "cand03_pyridylphenyl"]


def main():
    refs = {m.GetProp("_Name"): m for m in Chem.SDMolSupplier(SDF, removeHs=False) if m}
    gfnff = pd.read_csv(os.path.join(OUT, "phase6_qm_descriptors.csv")).set_index("compound")
    rows = []
    for name in FINALISTS:
        folder = os.path.join(OUT, "qm_runs", f"{name}_gfn2")
        ens = os.path.join(folder, "crest_ensemble.xyz")
        if not os.path.exists(ens):
            print(f"[skip] no crest_ensemble.xyz for {name}")
            continue
        mol = refs[name]
        polar_mask, elements_ref = p6.classify_polarity(mol)
        frames = p6.parse_crest_ensemble(ens)
        if frames[0]["elements"] != elements_ref:
            print(f"[warn] atom-order mismatch for {name}; skipping")
            continue
        desc = p6.ensemble_descriptors(frames, polar_mask, elements_ref, name)
        weights = desc.pop("weights")
        desc.pop("rel_kcal", None)
        rows.append(desc)
        # regenerate DFT inputs from the GFN2 populations
        gdir = os.path.join(GFN2_GAUSS, name)
        gpaths, gweights = p6.write_gaussian_inputs(frames, weights, name, gdir)
        print(f"[ok] {name}: {desc['n_conformers']} GFN2 conformers -> "
              f"{len(gpaths)} Gaussian inputs (cum pop {sum(gweights):.2f})")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase8b_finalist_gfn2_descriptors.csv"), index=False)

    # side-by-side GFN-FF vs GFN2 for the populated descriptors
    cmp_cols = ["n_conformers", "hydrophobic_sasa_mean", "hydrophobic_sasa_std",
                "polar_sasa_mean", "hydrophobic_fraction_mean", "rg_mean"]
    print("\nGFN-FF tier  ->  GFN2-reranked (Boltzmann-weighted):")
    cmp_rows = []
    for _, r in df.iterrows():
        name = r["compound"]
        rec = {"compound": name}
        for c in cmp_cols:
            a = gfnff.loc[name, c] if name in gfnff.index else float("nan")
            b = r[c]
            rec[f"{c}__gfnff"] = a
            rec[f"{c}__gfn2"] = b
            if c != "n_conformers":
                print(f"  {name:26s} {c:26s} {a:8.2f} -> {b:8.2f}  (d={b-a:+.2f})")
        cmp_rows.append(rec)
    pd.DataFrame(cmp_rows).to_csv(os.path.join(OUT, "phase8b_gfnff_vs_gfn2.csv"), index=False)
    print("\nWrote phase8b_finalist_gfn2_descriptors.csv, phase8b_gfnff_vs_gfn2.csv, "
          "and qm_runs_gaussian_gfn2/<cand>/*.gjf")


if __name__ == "__main__":
    main()
