"""Shared RDKit featurization for the Stage 1b free-fraction / PPB oracle.

One descriptor function is used for the training drugs, the papulacandins and the
echinocandins so the model sees every molecule on the same basis. The descriptor
set is deliberately physicochemical (not a fingerprint): these are the axes that
govern plasma-protein binding (size, lipophilicity, polar surface, H-bonding,
shape) and they remain computable for the large bRo5 glycolipids in this project.
"""

from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def _heteroatoms(m):
    return sum(1 for a in m.GetAtoms() if a.GetAtomicNum() not in (1, 6))


def _frac_aromatic_atoms(m):
    n = m.GetNumHeavyAtoms()
    return (sum(1 for a in m.GetAtoms() if a.GetIsAromatic()) / n) if n else 0.0


# (name, callable) — order is the feature-vector order.
DESCRIPTORS = [
    ("mw", Descriptors.MolWt),
    ("clogp", Crippen.MolLogP),
    ("tpsa", rdMolDescriptors.CalcTPSA),
    ("hbd", rdMolDescriptors.CalcNumHBD),
    ("hba", rdMolDescriptors.CalcNumHBA),
    ("rotb", rdMolDescriptors.CalcNumRotatableBonds),
    ("fsp3", rdMolDescriptors.CalcFractionCSP3),
    ("aromatic_rings", rdMolDescriptors.CalcNumAromaticRings),
    ("aliphatic_rings", rdMolDescriptors.CalcNumAliphaticRings),
    ("ring_count", rdMolDescriptors.CalcNumRings),
    ("heavy_atoms", lambda m: m.GetNumHeavyAtoms()),
    ("heteroatoms", _heteroatoms),
    ("labute_asa", rdMolDescriptors.CalcLabuteASA),
    ("frac_aromatic_atoms", _frac_aromatic_atoms),
]

FEATURE_NAMES = [name for name, _ in DESCRIPTORS]


def mol_from_smiles(smiles):
    if not smiles or not str(smiles).strip():
        return None
    return Chem.MolFromSmiles(str(smiles).strip())


def featurize(smiles):
    """Return the descriptor vector for a SMILES, or None if unparseable."""
    mol = mol_from_smiles(smiles)
    if mol is None:
        return None
    out = []
    for _, fn in DESCRIPTORS:
        try:
            out.append(float(fn(mol)))
        except Exception:
            return None
    return out


def murcko_scaffold(smiles):
    """Generic (atom-agnostic) Murcko scaffold SMILES, for scaffold splitting."""
    mol = mol_from_smiles(smiles)
    if mol is None:
        return ""
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        generic = MurckoScaffold.MakeScaffoldGeneric(scaf)
        return Chem.MolToSmiles(generic)
    except Exception:
        return ""
