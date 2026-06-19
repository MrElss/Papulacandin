# Papulacandin / Fusacandin serum gap — cross-phase synthesis (Phases 1–10)

*Question:* why do a few papulacandin-class FKS1 inhibitors keep antifungal
activity in serum while structurally close analogs lose it — and can we compute a
rule to design serum-tolerant analogs?

*Answer, in one line:* serum **potency** is governed by intrinsic potency ×
molecular size/lipophilicity; serum **tolerance** (the free→serum MIC shift) has, at
most, a weak association with exposing **polar rather than hydrophobic** surface —
a lead seen consistently across independent methods but never statistically
significant on the available data. No computable descriptor or docking model
predicts serum tolerance on this dataset; the binding constraint is the **data**,
not the methods.

---

## 1. The dataset (the dependent variable)
- 24 compounds with matched serum-free **and** serum whole-cell MIC vs *C. albicans*
  (Yeung-1996 fusacandin analogs 6a–6u + Fusacandin A/B + Papulacandin B).
- 13/24 retain measurable serum activity; **11/24 serum MICs are censored at 100**,
  and values are tied on {12.5, 25, 50, 100}. This censoring/tying is the ceiling on
  every correlation in the project.
- Endpoints used: raw **serum MIC**, and the **serum SHIFT fold** (serum/serum-free)
  — the potency-independent measure of true serum tolerance.

## 2. Phase-by-phase

| Phase | What | Headline result |
|---|---|---|
| **1–4** | Curate serum-gap table; 2D SAR; interpretable design-rule score | Best 2D predictor (rigid-aromatic + polar design score) vs serum MIC: **Spearman ρ = 0.32, n.s.** Motivated a 3D approach. |
| **5** | Scaffold-constrained generative design (cleave/re-esterify the C-6′ aromatic ester) → 12 novel candidates | Whole chemotype is **bRo5** (QED ≈ 0.01–0.03; 0/30 pass Lipinski) → 2D drug-likeness is uninformative; 3D needed. |
| **6** | QM funnel infrastructure: CREST ensemble parser, in-house Shrake–Rupley SASA, Boltzmann weighting, Gaussian DFT I/O | Pipeline + per-conformer 3D descriptors; DFT inputs for populated conformers. |
| **7** | Retrospective on 24 knowns with **MMFF** ensembles (fast proxy) | Hydrophobic SASA vs serum MIC **ρ = −0.45 (p = 0.029)** — beat 2D. Promising enough to justify real CREST. |
| **8** | Same retrospective with **real CREST/GFN-FF** ensembles + confound analysis | Signal **shrinks and fails to confirm**: ρ = −0.31 (p = 0.14). Decisive: serum-free MIC alone tracks serum MIC at **ρ = 0.79**; partialling potency out collapses hydrophobic SASA to **ρ = 0.02**. Rigidity hypothesis dies (ρ = +0.17). Only coherent lead: **polar SASA vs shift ρ = −0.33** (p = 0.12). |
| **9** | GFN2-xTB **electronic / solvation** descriptors (dipole, gap, polarizability, QM logP) | Polarizability is the biggest raw serum-MIC correlate (ρ = −0.54, p = 0.01) but is a **size proxy** (shift ρ = −0.02). **QM logP** tracks hydrophobic SASA (ρ = +0.76) and gives **shift ρ = +0.30** — independently corroborating the polar-surface lead. Two descriptor families converge; neither significant. |
| **10** | Explicit **HSA docking** (rigid ensemble surface docking; whole-molecule flexible docking is intractable at ~38 torsions) | **Null:** HSA affinity vs serum shift **ρ = +0.22 (p = 0.30)**, sign opposite to sequestration; flat under potency/size controls. Not a size proxy (ρ = −0.15 vs polarizability), so the model captured real surface association that simply doesn't track tolerance. |

## 3. The one convergent, defensible finding
Three independent computational families — 3D shape (polar SASA), QM solvation
(logP), and weakly the dipole — agree on a single qualitative rule and disagree
with nothing:

> **Among equipotent analogs, biasing the exposed surface toward polar / H-bonding
> character (higher polar SASA, lower QM logP) is associated with smaller serum
> loss.** Effect size |ρ| ≈ 0.30–0.33, p ≈ 0.12–0.15 — directional, not proven.

Equally important negatives, each now rigorously established rather than assumed:
- **Raw serum MIC is the wrong endpoint** — ~ potency-dominated (ρ = 0.79 with
  serum-free MIC); always model the *shift* or use serum-free MIC as a covariate.
- **The Phase-1 rigidity hypothesis is not supported** (SASA-spread shift ρ ≈ −0.1).
- **Bulk lipophilicity/size (clogP, Rg, polarizability) tracks potency, not
  tolerance** — the Phase-7 hydrophobic-SASA "hit" was a potency artifact.
- **Explicit HSA drug-site docking carries no signal** for these oversized ligands.

## 4. Design output (carried forward with caveats)
The 12 novel Phase-5 candidates were taken through the full funnel; the 3 finalists
by multi-objective score — **cand01 quinolinecarbonyl, cand02 naphthoyl-6-OH,
cand03 pyridylphenyl** (extended rigid aromatics) — have GFN2-reranked ensembles and
ready Gaussian DFT inputs. Given Phase 8–10, finalist *ranking* is unchanged but the
*rationale* shifts: prefer designs that add **polar/H-bonding surface** on the
variable C-6′ group, not just rigid aromatic bulk. This is a hypothesis to test, not
a validated predictor.

## 5. Limitations (honest)
- **n = 24, 11 censored, single chemotype/lab** — three convergent descriptors at
  |ρ| ≈ 0.30, p ≈ 0.13 is exactly what a real-but-modest effect looks like under this
  much censoring; it cannot be pushed to significance with more modeling.
- MMFF→CREST showed proxy ensembles **overstate** effects (−0.45 → −0.31): trust the
  QM-ensemble numbers, treat any fast-proxy result as a screen only.
- Vina is not parameterized for 1000–1200 Da amphiphiles; the Phase-10 null cannot
  *exclude* albumin sequestration, only say this proxy is uninformative.

## 6. Recommended next steps (now experimental, not computational)
1. **Acquire uncensored serum MICs** on more analogs — ideally a second chemotype —
   so the polar-surface lead can be tested at adequate power. This is the single
   biggest lever; everything else is capped by the data.
2. **Direct HSA-binding assay** (equilibrium dialysis or fluorescence-probe
   displacement at Sudlow I/II) to ground-truth the sequestration mechanism docking
   could not resolve; include the fatty-acid sites given the lipid tails.
3. If pursued in silico: model the **fatty-acid binding sites / other serum
   carriers** (AAG, lipoproteins) rather than more drug-site docking, and run the
   prepared **finalist DFT** single points to add ESP descriptors — but treat these
   as exploratory, not decision-grade, until (1) lands.

## 7. Asset map (all under `analysis/`)
- Data & 2D: `outputs/serum_gap_pairs.csv`, `phase1–4` outputs, `outputs/design_rules.md`
- Generative: `phase5_generate.py`, `outputs/phase5_*`
- QM funnel: `phase6_qm_layer.py`; CREST inputs/results in `outputs/qm_runs/`,
  `outputs/qm_runs_known/`; GFN2 finalist reranks in `outputs/qm_runs/cand0*_gfn2/`
- Retrospectives: `phase7_*`, `phase8_*` (+ `phase8_findings.md`, `phase8_confound_analysis.csv`)
- Electronics: `gen_known_xtb_inputs.py`, `phase9_electronic.py`, `outputs/phase9_*`
- Docking: `phase10_dock_hsa.py`, `outputs/phase10_*`
- Per-phase write-ups: `outputs/phase7_findings.md` … `phase10_findings.md`; this file:
  `outputs/SYNTHESIS_phases1-10.md`
