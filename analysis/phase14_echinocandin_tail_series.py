#!/usr/bin/env python3
"""
phase14_echinocandin_tail_series.py
===================================
PHASE 14 — Round-1 synthesis shortlist, redirected by echinocandin SAR.

Why this differs from Phase 13
------------------------------
Phase 13 chased exposed POLARITY on the tail; GFN2 (Step 4) killed those leads
(the flexible polar heads fold and bury). Meanwhile the echinocandins — the same
FKS/glucan-synthase target, the same long-tail liability, but with real clinical
serum data (our Phase 11 corpus) — point at a DIFFERENT axis:

  * the lipid tail is REQUIRED for potency (it anchors the drug to the membrane;
    removing/over-polarizing it kills activity), so the tail cannot be made polar
    or short without losing the antifungal;
  * among tails that keep potency, the marketed drug with the SMALLEST serum shift
    (caspofungin, ~2-11x) carries a BRANCHED SATURATED fatty tail, while the
    largest shifts (anidulafungin ~16-512x, micafungin ~64-1024x) carry RIGID
    EXTENDED AROMATIC tails (alkoxy-terphenyl, diphenyl-isoxazole).

The native papulacandin tail is a linear CONJUGATED POLYENE — rigid, planar,
extended — i.e. structurally in the anidulafungin/micafungin "high serum binding"
regime, NOT the caspofungin "flexible saturated" regime. So the echinocandin-
grounded hypothesis is:

  >> at matched chain length (hold lipophilicity/potency ~constant), replace the
  >> rigid conjugated polyene with a FLEXIBLE SATURATED / BRANCHED chain to move
  >> from the anidulafungin-like toward the caspofungin-like regime.

Crucially this is NOT a bulk-clogP play: the saturated analogs are MORE lipophilic
than the native polyene (conjugation lowers clogP), yet the flexible saturated
tail is predicted to bind serum protein/membrane LESS tightly because it lacks the
rigid planar shape — consistent with Phase 11's finding that bulk clogP does not
explain the serum shift. This rigidity axis was never tested in Phase 13 because
the polarity-biased proxy filtered saturated tails out.

The series (all C16, so length/lipophilicity is held; only rigidity varies)
--------------------------------------------------------------------------
  native C16 polyene (4 C=C, rigid/conjugated) -- reference baseline
  C16:1 palmitoleoyl (1 cis kink, fluid)       -- intermediate rigidity control
  C16:0 palmitoyl    (saturated, flexible)     -- the clean single test
  branched saturated (caspofungin-like)        -- closest analog to the low-shift drug

Synthetic accessibility (the practical answer to "which first")
---------------------------------------------------------------
These are ALSO the easiest tails to install: palmitic / palmitoleic / branched
fatty acids are cheap commodities put on in ONE esterification, vs the native
conjugated-polyene tail (a hard multi-step synthesis) or Phase-13's charged/PEG
tails. So the best scientific bet and the easiest chemistry coincide.

Outputs (analysis/outputs/)
---------------------------
* phase14_tail_series.csv         — the 4 targets + tail rigidity descriptors + rationale
* phase14_top_candidates.sdf      — 3D structures (CREST-ready; names match QM dirs)
* phase14_qm_runs/<name>/<name>.xyz + run_crest.sbatch  — optional QM confirmation
* phase14_findings.md             — rationale, priority order, and the free-drug caveat
"""
from __future__ import annotations

import os
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phase13_fatty_tail_optimization as p13  # noqa: E402 (guarded __main__)

p12 = p13.p12
OUT = p13.OUT
CORE = p13.CORE
QM_DIR = os.path.join(OUT, "phase14_qm_runs")
TEMPLATE = p13.TEMPLATE  # PAPU-0080

# label : (SMILES or None for native, priority, rationale, synthesis note)
LADDER = {
    "native_C16_polyene": (
        None, "reference",
        "rigid conjugated tetraene (4 C=C); the anidulafungin/micafungin-like "
        "high-serum-binding regime — this is what we are trying to beat.",
        "you likely already have the parent; hard (conjugated-polyene synthesis)."),
    "C16_0_palmitoyl": (
        "[*:1]C(=O)CCCCCCCCCCCCCCC", "MAKE FIRST",
        "saturated, fully flexible, same C16 length -> holds lipophilicity/potency, "
        "removes the rigidity/conjugation echinocandin SAR blames for serum binding. "
        "The single cleanest test of the hypothesis.",
        "trivial: palmitic acid (commodity), one esterification."),
    "branched_sat_caspofungin_like": (
        "[*:1]C(=O)CCCCCCCCCCCC(C)CC", "MAKE FIRST",
        "branched saturated ~C16 — the closest structural analog to caspofungin's "
        "tail, the marketed echinocandin with the SMALLEST serum shift.",
        "easy: branched (anteiso/iso) fatty acid, one esterification."),
    "C16_1_palmitoleoyl": (
        "[*:1]C(=O)CCCCCCC/C=C\\CCCCCC", "second",
        "one cis double bond (kinked, fluid) — intermediate rigidity control that "
        "separates 'conjugation/rigidity' from 'saturation': if it behaves like "
        "palmitoyl, the culprit is the CONJUGATION, not merely the double bonds.",
        "easy: palmitoleic/oleic acid, one esterification."),
}


def tail_descriptors(acyl):
    nc = sum(1 for a in acyl.GetAtoms() if a.GetAtomicNum() == 6)
    ndb = sum(1 for b in acyl.GetBonds()
              if b.GetBondType() == Chem.BondType.DOUBLE
              and b.GetBeginAtom().GetAtomicNum() == 6
              and b.GetEndAtom().GetAtomicNum() == 6)
    return dict(tail_n_carbon=nc, tail_CC_double_bonds=ndb,
                tail_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(acyl),
                tail_clogp=round(Crippen.MolLogP(acyl), 2))


def main():
    cm = pd.read_csv(os.path.join(CORE, "compounds_master.csv")).set_index("compound_id")
    tmpl = Chem.MolFromSmiles(cm.loc[TEMPLATE, "smiles_canonical"])
    core, native = p13.cleave_longest_fatty_tail(tmpl)
    assert core is not None, "core/tail cleavage failed"

    rows, mol_by_name = [], {}
    for name, (smi, priority, rationale, synth) in LADDER.items():
        acyl = native if smi is None else Chem.MolFromSmiles(smi)
        for a in acyl.GetAtoms():
            if a.GetAtomicNum() == 0:
                a.SetAtomMapNum(1)
        full = p12.attach(core, acyl)
        if full is None:
            print(f"[skip] attach failed for {name}")
            continue
        td = tail_descriptors(acyl)
        rows.append(dict(name=name, priority=priority,
                         mw=round(Descriptors.MolWt(full), 1),
                         whole_clogp=round(Crippen.MolLogP(full), 2),
                         **td, rationale=rationale, synthesis=synth,
                         smiles=Chem.MolToSmiles(full)))
        mol_by_name[name] = full

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "phase14_tail_series.csv"), index=False)
    _write_qm_inputs(df, mol_by_name)
    _write_findings(df)

    print("PHASE 14 — echinocandin-motivated tail series (round-1 shortlist)")
    print("=" * 64)
    print("All C16 (length/lipophilicity matched); only rigidity/saturation varies.\n")
    show = ["name", "priority", "tail_CC_double_bonds", "tail_rotatable_bonds",
            "tail_clogp", "whole_clogp", "mw"]
    print(df[show].to_string(index=False))
    print("\nNote: native polyene has the LOWEST tail clogP yet the WORST serum "
          "behavior -> it is rigidity/shape, not bulk lipophilicity (Phase 11).")
    print(f"\nWrote: phase14_tail_series.csv, phase14_top_candidates.sdf, "
          f"phase14_qm_runs/, phase14_findings.md")


def _write_qm_inputs(df, mol_by_name):
    os.makedirs(QM_DIR, exist_ok=True)
    writer = Chem.SDWriter(os.path.join(OUT, "phase14_top_candidates.sdf"))
    for i, r in df.iterrows():
        name = r["name"]
        m = p12._embed(mol_by_name[name])
        if m is None:
            continue
        m.SetProp("_Name", name)
        writer.write(m)
        d = os.path.join(QM_DIR, name)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f"{name}.xyz"), "w") as fh:
            fh.write(Chem.MolToXYZBlock(m))
    writer.close()
    with open(os.path.join(QM_DIR, "run_crest.sbatch"), "w") as fh:
        fh.write(p13.RUN_CREST_SBATCH.replace("Phase 13", "Phase 14").replace("p13_", "p14_"))


def _write_findings(df):
    def row(name):
        r = df[df["name"] == name].iloc[0]
        return (f"- **{name}** [{r['priority']}] — C={r['tail_n_carbon']}, "
                f"C=C={r['tail_CC_double_bonds']}, rot={r['tail_rotatable_bonds']}, "
                f"tail clogP {r['tail_clogp']}, whole clogP {r['whole_clogp']}\n"
                f"    - *why:* {r['rationale']}\n    - *synthesis:* {r['synthesis']}")
    order = ["C16_0_palmitoyl", "branched_sat_caspofungin_like",
             "C16_1_palmitoleoyl", "native_C16_polyene"]
    body = "\n".join(row(n) for n in order if (df["name"] == n).any())
    md = f"""# Phase 14 — which tail to synthesize first (echinocandin-guided)

## The redirect
Phase 13 optimized the tail for exposed POLARITY; GFN2 (Step 4) showed those
leads were force-field artifacts. The echinocandins point at a different, better-
evidenced axis: the tail is **required for potency** (membrane anchor), and among
potency-keeping tails the low-serum-shift drug (caspofungin) is **branched
saturated**, while the high-shift drugs (anidulafungin, micafungin) are **rigid
extended aromatics**. The native papulacandin tail is a rigid conjugated polyene —
in the high-shift regime. So: at matched length, **de-rigidify the tail**
(saturate / branch) rather than polarize it.

Key point: this is NOT bulk lipophilicity — the saturated analogs have HIGHER
clogP than the native polyene, yet are predicted to bind serum LESS because they
lack the rigid planar shape (consistent with Phase 11: bulk clogP does not track
serum shift).

## Priority order (best scientific bet AND easiest chemistry coincide)
{body}

## Why this order
1. **C16:0 palmitoyl first** — the single cleanest test (only rigidity changes vs
   native), and the easiest to make (palmitic acid, one esterification). If serum
   tolerance improves here, the rigidity hypothesis is validated on-scaffold.
2. **Branched saturated (caspofungin-like)** — the most direct read-across from the
   lowest-serum-shift marketed drug; tells you if branching adds over plain
   saturation.
3. **C16:1 palmitoleoyl** — separates conjugation from unsaturation (one isolated
   double bond, not conjugated). Cheap; run if palmitoyl is encouraging.
4. **native C16 polyene** — the reference; you must measure its serum shift in YOUR
   assay to interpret the rest (you likely already have the parent).

## Read-out and the free-drug caveat (from echinocandins)
Measure the serum SHIFT (protein-adjusted MIC ± albumin) AND fraction-unbound
(equilibrium dialysis) on this matched set. Interpret via the FREE drug: a tail
that is somewhat more protein-bound can still win if it stays potent (echinocandins
are ~96-99.8% bound yet effective). Do not reject an analog on total-drug MIC
alone. And keep serum-free potency in the readout — a de-rigidified tail must not
have quietly lost antifungal activity (the tail is essential for it).

## Optional QM pre-check (same funnel; Step-4 taught us to)
3D inputs are under `phase14_qm_runs/<name>/`. If you want a compute check before
synthesis, CREST them (GFN-FF screen) and — because Step 4 showed GFN-FF can
mislead — GFN2 re-rank, then compare exposed hydrophobic surface across the
rigidity ladder. But given the descriptor's poor track record for tails, the
decisive evidence is the wet-lab serum shift of these four.
"""
    with open(os.path.join(OUT, "phase14_findings.md"), "w", encoding="utf-8") as fh:
        fh.write(md)


if __name__ == "__main__":
    main()
