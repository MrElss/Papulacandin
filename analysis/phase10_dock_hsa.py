#!/usr/bin/env python3
"""
phase10_dock_hsa.py
===================
Phase 10 — does human serum albumin (HSA) binding explain the serum gap?
Direct mechanistic test of the recurring Phase 8/9 lead (exposed polar vs
hydrophobic surface), by docking each known compound to HSA and correlating
binding strength with the serum SHIFT.

METHOD (and why this specific method)
-------------------------------------
These ligands are 1000-1200 Da glycolipids with ~38 rotatable bonds -- far
outside the regime where flexible docking is reliable (Vina ~15 torsions) and
too large to fit rigidly into HSA's buried Sudlow drug pockets (rigid pocket
docking clashes, +40 kcal/mol). So we use RIGID ENSEMBLE SURFACE DOCKING:
  * receptor: HSA chain A (PDB 1AO6, defatted), rigid.
  * ligands : the top Boltzmann CREST/xTB conformers (real QM geometries) kept
              RIGID -- conformational flexibility is sampled by the ensemble,
              not by Vina's (unreliable here) torsion search.
  * boxes   : large (40 A) boxes over the two Sudlow regions (subdomains IIA &
              IIIA); a big box lets the rigid amphiphile settle into the nearest
              surface groove rather than being forced into a buried pocket.
  * score   : best (most negative) Vina affinity over {conformers x sites} =
              the compound's strongest accessible HSA association.

HONEST CONTROL
--------------
A rigid surface-contact score partly tracks ligand SIZE. So the load-bearing
test is NOT the raw correlation but the PARTIAL correlation of HSA score vs the
serum shift controlling for size/potency (polarizability alpha(0), and serum-free
MIC). If HSA binding predicts the shift only because bigger binds more, it adds
nothing over Phase 8/9; if it survives the control, it is genuine mechanism.

OUTPUTS
  outputs/phase10_docking/phase10_hsa_scores.csv      (raw per conf/site)
  outputs/phase10_docking/phase10_hsa_per_compound.csv
  outputs/phase10_findings via stdout + phase10_hsa.png
"""
from __future__ import annotations
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
DOCK = os.path.join(OUT, "phase10_docking")
SDF_DIR = os.path.join(HERE, "..", "curated", "structures", "sdf")
XTB_ROOT = os.path.join(OUT, "qm_runs_known_xtb")
RECEPTOR = os.path.join(DOCK, "hsa_receptor.pdbqt")
SITES = {"siteI_IIA": [30.84, 35.88, 33.75], "siteII_IIIA": [14.12, 27.99, 18.63]}
BOX = [40.0, 40.0, 40.0]
TOP_CONFS = 3
EXHAUST = 8


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


def rigid_pdbqt(cid, conf_xyz, out_path):
    """Map a CREST/xTB conformer onto the bonded SDF (AddHs order preserved),
    Meeko-prepare, then flatten to a rigid (TORSDOF 0) PDBQT."""
    sdf = glob.glob(os.path.join(SDF_DIR, f"{cid}_*.sdf"))[0]
    mol = Chem.AddHs(Chem.SDMolSupplier(sdf, removeHs=False)[0], addCoords=True)
    lines = open(conf_xyz).read().splitlines()[2:]
    xyz = np.array([[float(x) for x in l.split()[1:4]] for l in lines])
    if mol.GetNumAtoms() != len(xyz):
        return False
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, p.tolist())
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    s, ok, _ = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(mol)[0])
    if not ok:
        return False
    atoms = [l for l in s.splitlines() if l.startswith(("ATOM", "HETATM"))]
    with open(out_path, "w") as fh:
        fh.write("ROOT\n" + "\n".join(atoms) + "\nENDROOT\nTORSDOF 0\n")
    return True


def prepare_ligands():
    pairs = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv"))
    lig = {}
    for _, r in pairs.iterrows():
        cid = r["compound_id"]
        wfile = os.path.join(XTB_ROOT, cid, "weights.csv")
        if not os.path.exists(wfile):
            continue
        w = pd.read_csv(wfile).sort_values("pop", ascending=False).head(TOP_CONFS)
        ld = os.path.join(DOCK, "ligands", cid)
        os.makedirs(ld, exist_ok=True)
        paths = []
        for _, cr in w.iterrows():
            cx = os.path.join(XTB_ROOT, cid, f"{cr['conformer']}.xyz")
            pq = os.path.join(ld, f"{cr['conformer']}.pdbqt")
            if os.path.exists(cx) and rigid_pdbqt(cid, cx, pq):
                paths.append(pq)
        if paths:
            lig[cid] = paths
    return lig, pairs


def dock_all(lig):
    rows = []
    for site, center in SITES.items():
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(RECEPTOR)
        v.compute_vina_maps(center=center, box_size=BOX)
        for cid, paths in lig.items():
            for pq in paths:
                v.set_ligand_from_file(pq)
                try:
                    v.dock(exhaustiveness=EXHAUST, n_poses=3)
                    aff = v.energies(n_poses=1)[0][0]
                except Exception as e:
                    aff = np.nan
                rows.append(dict(compound_id=cid, site=site,
                                 conformer=os.path.basename(pq).replace(".pdbqt", ""),
                                 affinity=float(aff)))
                print(f"  {cid} {site} {os.path.basename(pq):14s} {aff:7.2f}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DOCK, "phase10_hsa_scores.csv"), index=False)
    return df


def analyze(df, pairs):
    pairs["serum_mic_num"] = pairs["serum_mic_ugml"].apply(parse_mic)
    pairs["serumfree_mic_num"] = pairs["serumfree_mic_ugml"].apply(parse_mic)
    # best (strongest) HSA association per compound, over conformers x sites
    best = df.groupby("compound_id")["affinity"].min().rename("hsa_best_affinity")
    mean = df.groupby("compound_id")["affinity"].mean().rename("hsa_mean_affinity")
    m = pairs.merge(best, on="compound_id").merge(mean, on="compound_id")
    # bring in size/electronic covariate (polarizability) if available
    ph9 = os.path.join(OUT, "phase9_electronic_descriptors.csv")
    if os.path.exists(ph9):
        m = m.merge(pd.read_csv(ph9)[["compound_id", "alpha_au"]], on="compound_id", how="left")
    m.to_csv(os.path.join(DOCK, "phase10_hsa_per_compound.csv"), index=False)

    shift = np.log10(m["serum_shift_fold"].values)
    smic = np.log10(m["serum_mic_num"].values)
    sfree = np.log10(m["serumfree_mic_num"].values)

    def sp(a, b):
        r, p = stats.spearmanr(a, b)
        return r, p

    def partial(x, y, z):
        from scipy.stats import rankdata
        xr, yr, zr = rankdata(x), rankdata(y), rankdata(z)
        rx = xr - np.polyval(np.polyfit(zr, xr, 1), zr)
        ry = yr - np.polyval(np.polyfit(zr, yr, 1), zr)
        return stats.pearsonr(rx, ry)

    print("\n" + "=" * 68)
    print(f"PHASE 10 — HSA rigid-ensemble docking vs serum data (n={len(m)})")
    print("(more negative affinity = stronger HSA binding; sequestration model:")
    print(" stronger binding -> larger serum shift -> negative rho expected)")
    print("=" * 68)
    aff = m["hsa_best_affinity"].values
    for lab, y in [("serum SHIFT", shift), ("serum MIC", smic), ("serum-free MIC", sfree)]:
        r, p = sp(aff, y)
        print(f"  HSA best-affinity vs {lab:16s}: rho={r:+.2f} (p={p:.3f})")
    print("\nKEY CONTROLS (does HSA binding predict the shift beyond size?):")
    for cov_name, cov in [("serum-free MIC", sfree)]:
        r, p = partial(aff, shift, cov)
        print(f"  partial(HSA vs shift | {cov_name:14s}): rho={r:+.2f} (p={p:.3f})")
    if "alpha_au" in m and m["alpha_au"].notna().all():
        r, p = partial(aff, shift, m["alpha_au"].values)
        print(f"  partial(HSA vs shift | polarizability ): rho={r:+.2f} (p={p:.3f})")
        rr, pp = sp(aff, m["alpha_au"].values)
        print(f"  [check] HSA affinity vs polarizability(size): rho={rr:+.2f} (p={pp:.3f})")

    # figure
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    colors = m["serum_active"].map({"yes": "#2a9d8f", "no": "#e76f51"})
    ax.scatter(aff, m["serum_shift_fold"], c=colors, s=60, edgecolor="k", linewidth=0.4)
    ax.set_yscale("log")
    ax.set_xlabel("Best HSA Vina affinity (kcal/mol, more neg = stronger)")
    ax.set_ylabel("Serum shift fold (log)")
    r, p = sp(aff, shift)
    ax.set_title(f"Phase 10: HSA binding vs serum shift\nSpearman rho={r:+.2f} (p={p:.3f}, n={len(m)})")
    for _, rr in m.iterrows():
        ax.annotate(str(rr["name"]).replace("Fusacandin analog ", ""),
                    (rr["hsa_best_affinity"], rr["serum_shift_fold"]),
                    fontsize=6, xytext=(2, 2), textcoords="offset points")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "phase10_hsa.png"), dpi=140)
    plt.close(fig)
    print("\nWrote phase10_docking/phase10_hsa_scores.csv, phase10_hsa_per_compound.csv, "
          "outputs/phase10_hsa.png")


if __name__ == "__main__":
    print("Preparing rigid ligand ensembles ...")
    lig, pairs = prepare_ligands()
    print(f"  {len(lig)} compounds, {sum(len(v) for v in lig.values())} conformers")
    print("Docking into HSA (2 Sudlow-region boxes) ...")
    df = dock_all(lig)
    analyze(df, pairs)
