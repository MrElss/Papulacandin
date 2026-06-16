#!/usr/bin/env python3
"""
phase6_qm_layer.py
===================
PHASE 6 (Tier 4) — QM integration layer: turn CREST conformer ensembles into
the 3D/bRo5-appropriate descriptors that Phase 5 showed 2D rules cannot
capture, and prepare/parse the Gaussian DFT follow-up on the populated
conformers.

WHY THIS EXISTS
---------------
Phase 5 found: (a) QED/Ro5 are saturated near zero for this whole bRo5
chemotype (uninformative), and (b) flat 2D design-rule scores only weakly
track observed serum MIC retrospectively (Spearman rho=0.32, n.s.). Both
point at the same gap: the discriminating signal is conformational/3D, not
2D-topological. This module computes the natural 3D analogue of the Phase-1
flexibility hypothesis directly from a CREST ensemble:
    - Boltzmann-weighted exposed HYDROPHOBIC vs POLAR solvent-accessible
      surface area (SASA) -> the structural correlate of serum-protein
      (e.g. albumin) binding.
    - Ensemble SPREAD of hydrophobic SASA (weighted std) -> a genuine 3D
      flexibility metric: does the molecule sample conformers that expose
      large hydrophobic patches, even if the global minimum does not?
    - Radius of gyration / asphericity -> overall 3D shape/compactness.

No external SASA dependency: freesasa failed to build in this sandbox, so
SASA is computed in-house with a vectorized Shrake-Rupley algorithm (Bondi
radii, 1.4 A water probe, golden-spiral point sets) — pure numpy + RDKit,
portable to whatever machine parses the real CREST output later.

PIPELINE STAGES
----------------
1. parse_crest_ensemble()  — read CREST's multi-frame XYZ (energy in the
   comment line, Hartree, lowest first — standard CREST output format),
   recover Boltzmann weights at 298.15 K from the GFN2/GFN-FF energies.
2. ensemble_descriptors()  — per-conformer SASA breakdown + shape, then
   Boltzmann-weighted aggregation -> one row per compound.
3. write_gaussian_inputs() — pick the top Boltzmann-populated conformers
   (cumulative population >= POP_CUTOFF, capped at MAX_DFT_CONFS) and emit
   Gaussian .gjf single-point inputs (DFT charges/dipole) for the funnel's
   final, most expensive step.
4. parse_gaussian_log()    — extract SCF energy, dipole, ESP (Merz-Kollman)
   charges from a returned .gjf log; Boltzmann-average back up to the
   compound level using the CREST-derived weights.

A self-test at the bottom builds a SYNTHETIC ensemble (RDKit multi-conformer
embedding + MMFF energies) for one real Phase-5 candidate, in CREST's exact
file format, and runs the full parse -> descriptor -> Gaussian-input path on
it end-to-end. This validates the code now, before any real CREST/Gaussian
output exists. Swap in real `crest_conformers.xyz` paths once your fast
screening-tier runs finish (see phase5_crest_commands.sh for candidate
naming) and re-run with REAL_RUN=True.
"""
from __future__ import annotations
import os
import re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)

KCAL_PER_HARTREE = 627.5094740631
RT_298 = 1.987204e-3 * 298.15          # kcal/mol
PROBE_RADIUS = 1.4                      # water probe, Angstrom
N_SASA_POINTS = 194                     # golden-spiral points per atom sphere
POP_CUTOFF = 0.50                       # cumulative Boltzmann population for DFT subset
MAX_DFT_CONFS = 5

# Bondi-type van der Waals radii (Angstrom)
VDW_RADII = {"H": 1.10, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
             "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}
POLAR_ELEMENTS = {"O", "N", "S", "P"}

# ==========================================================================
# Stage 0: atom polarity classification from a reference RDKit mol
# ==========================================================================
def classify_polarity(mol):
    """Return bool array: True = polar (heteroatom O/N/S/P, or H bonded to one)."""
    mol = Chem.AddHs(mol) if mol.GetNumAtoms() == sum(
        1 for a in mol.GetAtoms() if a.GetAtomicNum() != 1) else mol
    polar = np.zeros(mol.GetNumAtoms(), dtype=bool)
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        if sym in POLAR_ELEMENTS:
            polar[atom.GetIdx()] = True
        elif sym == "H":
            nbrs = atom.GetNeighbors()
            if nbrs and nbrs[0].GetSymbol() in POLAR_ELEMENTS:
                polar[atom.GetIdx()] = True
    elements = [a.GetSymbol() for a in mol.GetAtoms()]
    return polar, elements

# ==========================================================================
# Stage 1: vectorized Shrake-Rupley SASA
# ==========================================================================
def _golden_sphere_points(n):
    i = np.arange(0, n)
    phi = np.arccos(1 - 2 * (i + 0.5) / n)
    golden = np.pi * (1 + 5 ** 0.5)
    theta = golden * i
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return np.stack([x, y, z], axis=1)

_SPHERE = _golden_sphere_points(N_SASA_POINTS)

def per_atom_sasa(coords, elements, probe=PROBE_RADIUS):
    """Shrake-Rupley SASA per atom (Angstrom^2). coords: (N,3), elements: list[str]."""
    n = len(elements)
    radii = np.array([VDW_RADII.get(e, 1.70) + probe for e in elements])
    sasa = np.zeros(n)
    for i in range(n):
        d = coords - coords[i]
        dist = np.linalg.norm(d, axis=1)
        cand = np.where((dist > 0) & (dist < (radii + radii[i])))[0]
        pts = coords[i] + _SPHERE * radii[i]
        buried = np.zeros(len(pts), dtype=bool)
        for j in cand:
            dd = np.linalg.norm(pts - coords[j], axis=1)
            buried |= dd < radii[j]
            if buried.all():
                break
        exposed_frac = 1.0 - buried.mean()
        sasa[i] = exposed_frac * 4 * np.pi * radii[i] ** 2
    return sasa

def shape_descriptors(coords, masses):
    """Radius of gyration and asphericity from the mass-weighted gyration tensor."""
    com = np.average(coords, axis=0, weights=masses)
    x = coords - com
    gyr = np.einsum("i,ia,ib->ab", masses, x, x) / masses.sum()
    eigval = np.sort(np.linalg.eigvalsh(gyr))[::-1]  # descending
    rg = np.sqrt(eigval.sum())
    asphericity = eigval[0] - 0.5 * (eigval[1] + eigval[2])
    return float(rg), float(asphericity)

ATOMIC_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
               "F": 18.998, "P": 30.974, "S": 32.06, "Cl": 35.45,
               "Br": 79.90, "I": 126.90}

# ==========================================================================
# Stage 2: CREST ensemble parsing (standard crest_conformers.xyz format)
# ==========================================================================
def parse_crest_ensemble(xyz_path):
    """Parse a CREST multi-frame XYZ. Comment line of each frame is the
    absolute electronic energy in Hartree (standard crest_conformers.xyz
    convention; CREST sorts frames lowest-energy first and preserves the
    atom order of the input structure across all frames)."""
    frames = []
    with open(xyz_path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        natoms = int(lines[i].strip())
        energy = float(re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", lines[i + 1])[0])
        elements, coords = [], []
        for j in range(i + 2, i + 2 + natoms):
            parts = lines[j].split()
            elements.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        frames.append(dict(elements=elements, coords=np.array(coords), energy_hartree=energy))
        i += 2 + natoms
    return frames

def boltzmann_weights(energies_hartree, T_kelvin=298.15):
    e_kcal = np.array(energies_hartree) * KCAL_PER_HARTREE
    rel = e_kcal - e_kcal.min()
    rt = 1.987204e-3 * T_kelvin
    w = np.exp(-rel / rt)
    return w / w.sum(), rel

# ==========================================================================
# Stage 3: per-compound Boltzmann-weighted descriptor aggregation
# ==========================================================================
def ensemble_descriptors(frames, polar_mask, elements_ref, compound_name):
    elements = frames[0]["elements"]
    if elements != elements_ref:
        raise ValueError(
            f"{compound_name}: CREST atom order/elements do not match the "
            f"reference RDKit mol — re-export the input xyz from the same "
            f"SDF used to build phase5_top_candidates.sdf.")
    energies = [f["energy_hartree"] for f in frames]
    weights, rel_kcal = boltzmann_weights(energies)
    masses = np.array([ATOMIC_MASS.get(e, 12.0) for e in elements])

    hphob, polar_sasa, total_sasa, rg_list, asph_list = [], [], [], [], []
    for f in frames:
        sasa = per_atom_sasa(f["coords"], elements)
        hphob.append(sasa[~polar_mask].sum())
        polar_sasa.append(sasa[polar_mask].sum())
        total_sasa.append(sasa.sum())
        rg, asph = shape_descriptors(f["coords"], masses)
        rg_list.append(rg)
        asph_list.append(asph)
    hphob, polar_sasa, total_sasa = map(np.array, (hphob, polar_sasa, total_sasa))
    rg_arr, asph_arr = np.array(rg_list), np.array(asph_list)

    def wmean(x):
        return float(np.sum(weights * x))

    def wstd(x):
        m = wmean(x)
        return float(np.sqrt(np.sum(weights * (x - m) ** 2)))

    return dict(
        compound=compound_name,
        n_conformers=len(frames),
        energy_window_kcal=round(float(rel_kcal.max()), 2),
        hydrophobic_sasa_mean=round(wmean(hphob), 1),
        hydrophobic_sasa_std=round(wstd(hphob), 1),     # ensemble flexibility metric
        polar_sasa_mean=round(wmean(polar_sasa), 1),
        total_sasa_mean=round(wmean(total_sasa), 1),
        hydrophobic_fraction_mean=round(wmean(hphob / total_sasa), 3),
        rg_mean=round(wmean(rg_arr), 2),
        rg_std=round(wstd(rg_arr), 2),
        asphericity_mean=round(wmean(asph_arr), 2),
        weights=weights, rel_kcal=rel_kcal,  # kept for write_gaussian_inputs
    )

# ==========================================================================
# Stage 4: Gaussian input generation for the top Boltzmann-populated conformers
# ==========================================================================
GAUSSIAN_ROUTE = "#P B3LYP/6-31G(d) SCRF=(PCM,Solvent=Water) Pop=MK NoSymm"

def write_gaussian_inputs(frames, weights, compound_name, out_dir,
                           charge=0, mult=1, route=GAUSSIAN_ROUTE):
    """Emit .gjf inputs for the conformers covering >= POP_CUTOFF cumulative
    Boltzmann population (capped at MAX_DFT_CONFS) — the funnel's most
    expensive step should only ever touch a handful of structures."""
    os.makedirs(out_dir, exist_ok=True)
    order = np.argsort(weights)[::-1]
    cum, chosen = 0.0, []
    for idx in order:
        chosen.append(idx)
        cum += weights[idx]
        if cum >= POP_CUTOFF or len(chosen) >= MAX_DFT_CONFS:
            break
    paths = []
    for rank, idx in enumerate(chosen, 1):
        f = frames[idx]
        fname = f"{compound_name}_conf{rank:02d}_pop{weights[idx]:.2f}.gjf"
        fpath = os.path.join(out_dir, fname)
        with open(fpath, "w") as fh:
            fh.write(f"%chk={compound_name}_conf{rank:02d}.chk\n%nprocshared=8\n%mem=16GB\n")
            fh.write(f"{route}\n\n{compound_name} conformer {rank} (Boltzmann pop={weights[idx]:.3f})\n\n")
            fh.write(f"{charge} {mult}\n")
            for el, xyz in zip(f["elements"], f["coords"]):
                fh.write(f"{el:2s} {xyz[0]:14.8f} {xyz[1]:14.8f} {xyz[2]:14.8f}\n")
            fh.write("\n")
        paths.append(fpath)
    return paths, [weights[i] for i in chosen]

# ==========================================================================
# Stage 5: Gaussian log parser (run after DFT jobs return)
# ==========================================================================
def parse_gaussian_log(log_path):
    """Extract final SCF energy (Hartree), dipole magnitude (Debye), and
    Merz-Kollman ESP charges from a Gaussian .log/.out file."""
    with open(log_path) as fh:
        text = fh.read()
    scf = re.findall(r"SCF Done:\s*E\([^)]*\)\s*=\s*([-\d.]+)", text)
    dip = re.findall(r"Dipole moment \(field-independent basis, Debye\):\s*\n\s*X=\s*([-\d.]+)\s*Y=\s*([-\d.]+)\s*Z=\s*([-\d.]+)\s*Tot=\s*([-\d.]+)", text)
    charges = []
    m = re.search(r"Charges from ESP fit.*?\n(.*?)\n\s*Sum of ESP charges", text, re.S)
    if m:
        for line in m.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) >= 3 and parts[0].isdigit():
                charges.append((parts[1], float(parts[2])))
    return dict(
        scf_energy_hartree=float(scf[-1]) if scf else np.nan,
        dipole_debye=float(dip[-1][3]) if dip else np.nan,
        esp_charges=charges,
    )

# ==========================================================================
# Orchestration: run over a directory of {compound_name}/crest_conformers.xyz
# ==========================================================================
def run_qm_layer(sdf_path, crest_root, out_csv, gaussian_out_root, real_run=True):
    """sdf_path: phase5_top_candidates.sdf (reference topology/atom order).
       crest_root: directory containing one subfolder per compound name with
       a crest_conformers.xyz inside (the candidate's fast-screening-tier
       output)."""
    refs = {m.GetProp("_Name"): m for m in Chem.SDMolSupplier(sdf_path, removeHs=False) if m}
    rows = []
    for name, mol in refs.items():
        xyz_path = os.path.join(crest_root, name, "crest_conformers.xyz")
        if not os.path.exists(xyz_path):
            if real_run:
                print(f"[skip] no CREST ensemble yet for {name} ({xyz_path})")
            continue
        polar_mask, elements_ref = classify_polarity(mol)
        frames = parse_crest_ensemble(xyz_path)
        desc = ensemble_descriptors(frames, polar_mask, elements_ref, name)
        weights = desc.pop("weights")
        desc.pop("rel_kcal")
        rows.append(desc)
        gdir = os.path.join(gaussian_out_root, name)
        gpaths, gweights = write_gaussian_inputs(frames, weights, name, gdir)
        print(f"[ok] {name}: {desc['n_conformers']} conformers -> "
              f"{len(gpaths)} Gaussian inputs (cum pop {sum(gweights):.2f})")
    df = pd.DataFrame(rows)
    if len(df):
        df.to_csv(out_csv, index=False)
    return df

# ==========================================================================
# Self-test: synthetic CREST-format ensemble for one real Phase-5 candidate
# ==========================================================================
def _self_test():
    print("=" * 70)
    print("SELF-TEST: synthetic ensemble (RDKit ETKDG + MMFF), CREST xyz format")
    print("Purpose: validate parse -> descriptor -> Gaussian-input code paths")
    print("BEFORE any real CREST output exists. Not a QM result.")
    print("=" * 70)
    sdf_path = os.path.join(OUT, "phase5_top_candidates.sdf")
    mol = next(m for m in Chem.SDMolSupplier(sdf_path, removeHs=False) if m)
    name = mol.GetProp("_Name")
    polar_mask, elements_ref = classify_polarity(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 1
    params.useRandomCoords = True
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=12, params=params))
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid)
        ff.Minimize(maxIts=500)
        energies.append(ff.CalcEnergy())
    order = np.argsort(energies)

    mock_xyz = os.path.join(OUT, f"_selftest_{name}_crest_conformers.xyz")
    with open(mock_xyz, "w") as fh:
        for rank in order:
            conf = mol.GetConformer(int(cids[rank]))
            # fabricate a plausible Hartree-scale energy from MMFF kcal/mol so the
            # parser's unit handling is exercised exactly as it will be on real data
            fake_hartree = -2000.0 + energies[rank] / KCAL_PER_HARTREE
            fh.write(f"{mol.GetNumAtoms()}\n{fake_hartree:.8f}\n")
            for atom in mol.GetAtoms():
                p = conf.GetAtomPosition(atom.GetIdx())
                fh.write(f"{atom.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}\n")

    frames = parse_crest_ensemble(mock_xyz)
    desc = ensemble_descriptors(frames, polar_mask, elements_ref, name)
    weights = desc.pop("weights")
    desc.pop("rel_kcal")
    print(f"\nParsed {desc['n_conformers']} synthetic conformers for {name}")
    for k, v in desc.items():
        print(f"  {k:28s} {v}")

    gdir = os.path.join(OUT, "_selftest_gaussian_inputs", name)
    paths, gweights = write_gaussian_inputs(frames, weights, name, gdir)
    print(f"\nWrote {len(paths)} Gaussian .gjf inputs (cumulative pop {sum(gweights):.2f}):")
    for p, w in zip(paths, gweights):
        print(f"  pop={w:.2f}  {os.path.relpath(p, OUT)}")
    print("\nSelf-test passed: parse -> SASA/shape -> Boltzmann aggregate -> "
          "Gaussian input all executed without error.")
    return desc

if __name__ == "__main__":
    desc = _self_test()
    print("\n" + "=" * 70)
    print("NEXT STEP (real data): once your fast screening-tier CREST runs")
    print("finish, lay out crest_conformers.xyz under e.g.")
    print("  analysis/outputs/qm_runs/<candidate_name>/crest_conformers.xyz")
    print("and call:")
    print("  run_qm_layer('outputs/phase5_top_candidates.sdf',")
    print("               'outputs/qm_runs', 'outputs/phase6_qm_descriptors.csv',")
    print("               'outputs/qm_runs_gaussian')")
    print("=" * 70)
