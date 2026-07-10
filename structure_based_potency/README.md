# Structure-based potency + property-prediction (session artifacts)

AI-orchestrated work toward a **serum-robust papulacandin**, covering (1) the verified
FKS1 target, (2) a rational-design analog set, (3) a physicochemical property panel, and
(4) a pocket-constrained FKS1 binding-affinity screen. All predictions are **orientation-grade**
(model-internal scales, out-of-domain scaffold) unless stated otherwise.

## 1. Target — FKS1 (verified)

- *S. cerevisiae* β-1,3-glucan synthase Fks1, **UniProt P38631**, **1876 aa**.
- Sequence integrity **checksum-verified**: CRC64 `AD4B4CB8CB28B5D8` (matches UniProt). File: `../target/fks1_clean.seq`.

### Experimental structure landscape (verified on RCSB)
| PDB | State | Study |
|---|---|---|
| 8JZN | apo, catalytic region (2.47 Å) | You/Zhao, *Sci Adv* 2023 |
| 7YUY | **S643P** echinocandin-resistant mutant | Hu, *Nature* 2023 |
| 8WL6 / 8WLA | FKS1–**Rho1** complex (activated) | Li, *Nat Commun* 2025 |
| **9PE1** | **caspofungin–Rho1–long-glucan ternary** (active-state inhibitor template) | Ren, *Nature* 2026 |
| **9WZU** | **enfumafungin-bound** (membrane-anchored amphiphile template) | You, *Nat Commun* 2026 |

### Binding pocket (papulacandin working hypothesis)
Shared TM5–TM8 amphiphile/glucan-gate site, from the caspofungin (9PE1) and enfumafungin (9WZU) contacts:
**E635, Y638, F639, S643, W695 (TM5/6); R1357, L1360, S1361, I1364, V1365, I1368 (HS2/TM8).**
0-indexed for Boltz: `[634, 637, 638, 642, 694, 1356, 1359, 1360, 1363, 1364, 1367]`.

## 2. Rational designs (`design_smiles.csv`, built by `build_designs.py`)

Papulacandin B has two acyl positions (pos1 = long-tail ester on the glucoside core; pos2 = short-tail
ester on the galactose). Native rebuild is byte-identical to verified papulacandin B (CID 6436198, C47H64O17).
The set varies three axes:
- **Tail shape:** native rigid polyene → one-kink → palmitoyl (flexible) → caspofungin branched → anidulafungin/micafungin (rigid aromatic, grafted from the real echinocandin side chains).
- **Linkage:** ester → **amide** (esterase-hardened bioisostere).
- **Position:** pos1 only vs **dual** (both acyl positions).

## 3. Property panel (`rational_design_panel.csv` / `.png`; `build_panel2.py`)

Boltz **ADME-v1** (permeability, lipophilicity) + **Inductive Bio** LogD@7.4. Key reads:
- Native papulacandin B = **permeability floor**; every de-rigidified tail raises it (delivery-device mechanism).
- **Ester→amide swap raises predicted permeability at ~constant LogD** — esterase hardening is "free" on these axes.
- **Rigid aromatic grafts rank below flexible aliphatic tails** (shape effect).
- **Caveat:** ADME-v1 / Inductive Bio do **not** provide fraction-unbound, albumin binding, serum-shift, or potency. Those require (fu) an in-domain model built from the echinocandin PPB/Fu data in `external/` + equilibrium dialysis; (albumin) MD/free-energy; (serum-shift) the wet-lab RPMI±50% serum assay; (potency) the screen in §4.

## 4. FKS1 affinity screen (potency term)

Boltz `small_molecule_screen`, **pocket-constrained** to the residues above, **built-in PAINS/druglikeness filter disabled** (it rejects this bRo5 chemotype).
- Job: `sm_scr_ZM9jvNzI0fbBbcJk5ePR` (workspace `ws_8knSS6VFwWhir9lE0F1d`).
- Share link: https://lab.boltz.bio/share/shr_1M31aOmbXK846vrwt-899MP2U03VbY35k696JEEreR4
- **Controls:** ibrexafungerp + caspofungin (positive, bind this site) vs fluconazole (negative). Validity check = positives high, negative low.
- **Caveat:** the screen folds FKS1 from sequence (no experimental template). Top hits should be re-scored with a **9WZU-template-anchored** `structure_and_binding` run and validated against resistance genetics + MD.

## Reproduce
`build_designs.py` → design SMILES; `build_panel.py` / `build_panel2.py` → property panels + figures.
Predictions were run via the Boltz and Inductive Bio MCP tools (see job IDs above).
