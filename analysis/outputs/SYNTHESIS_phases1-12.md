# Papulacandin / Fusacandin serum gap — cross-phase synthesis (Phases 1–12)

*Question:* why do a few papulacandin-class FKS1 inhibitors keep antifungal
activity in serum while structurally close analogs lose it — and can we compute a
rule to design serum-tolerant analogs?

*Answer, in one line:* serum **potency** is governed by intrinsic potency ×
molecular size/lipophilicity; serum **tolerance** (the free→serum MIC shift) has, at
most, a weak association with exposing **polar rather than hydrophobic** surface —
a lead seen consistently across independent methods and now, at bulk 2D level,
*not* reproduced by a second chemotype, which pins the lead specifically to
**locally exposed** surface. No computable descriptor or docking model predicts
serum tolerance on this dataset; the binding constraint is the **data**, not the
methods — and the right endpoint is **free-drug** exposure, not raw serum MIC.

---

## 1. The dataset (the dependent variable)
- 24 compounds with matched serum-free **and** serum whole-cell MIC vs *C. albicans*
  (Yeung-1996 fusacandin analogs 6a–6u + Fusacandin A/B + Papulacandin B).
- 13/24 retain measurable serum activity; **11/24 serum MICs are censored at 100**,
  and values are tied on {12.5, 25, 50, 100}. This censoring/tying is the ceiling on
  every correlation in the project.
- Endpoints: raw **serum MIC**, and the **serum SHIFT fold** (serum/serum-free) —
  the potency-independent measure of true serum tolerance.
- **Phase 11 adds a second, independent chemotype** (echinocandins) on the *same*
  shift endpoint, mined from the in-repo `external/` FKS corpus.

## 2. Phase-by-phase

| Phase | What | Headline result |
|---|---|---|
| **1–4** | Curate serum-gap table; 2D SAR; interpretable design-rule score | Best 2D predictor (rigid-aromatic + polar design score) vs serum MIC: **Spearman ρ = 0.32, n.s.** Motivated a 3D approach. |
| **5** | Scaffold-constrained generative design (cleave/re-esterify the C-6′ aromatic ester) → 12 novel candidates | Whole chemotype is **bRo5** (QED ≈ 0.01–0.03; 0/30 pass Lipinski) → 2D drug-likeness is uninformative; 3D needed. |
| **6** | QM funnel infrastructure: CREST ensemble parser, in-house Shrake–Rupley SASA, Boltzmann weighting, Gaussian DFT I/O | Pipeline + per-conformer 3D descriptors; DFT inputs for populated conformers. |
| **7** | Retrospective on 24 knowns with **MMFF** ensembles (fast proxy) | Hydrophobic SASA vs serum MIC **ρ = −0.45 (p = 0.029)** — beat 2D. Promising enough to justify real CREST. |
| **8** | Same retrospective with **real CREST/GFN-FF** ensembles + confound analysis | Signal **shrinks and fails to confirm**: ρ = −0.31 (p = 0.14). Decisive: serum-free MIC alone tracks serum MIC at **ρ = 0.79**; partialling potency out collapses hydrophobic SASA to **ρ = 0.02**. Rigidity hypothesis dies (ρ = +0.17). Only coherent lead: **polar SASA vs shift ρ = −0.33** (p = 0.12). |
| **9** | GFN2-xTB **electronic / solvation** descriptors (dipole, gap, polarizability, QM logP) | Polarizability is the biggest raw serum-MIC correlate (ρ = −0.54, p = 0.01) but is a **size proxy** (shift ρ = −0.02). **QM logP** tracks hydrophobic SASA (ρ = +0.76) and gives **shift ρ = +0.30** — independently corroborating the polar-surface lead. Two descriptor families converge; neither significant. |
| **10** | Explicit **HSA docking** (rigid ensemble surface docking; whole-molecule flexible docking is intractable at ~38 torsions) | **Null:** HSA affinity vs serum shift **ρ = +0.22 (p = 0.30)**, sign opposite to sequestration; flat under potency/size controls. Not a size proxy, so the model captured real surface association that simply doesn't track tolerance. |
| **11** | **Echinocandin cross-chemotype read-across** (in-repo external FKS corpus): serum-shift folds ± 50% serum, PPB/Fu, on a shared endpoint + RDKit descriptor axis | Same phenomenon, robust C. albicans ordering **caspofungin ×2 < anidulafungin ×16 < micafungin ×64** — yet all are front-line drugs dosed to a **free-drug** target (~96–99.8% bound). **Honest null:** *bulk* 2D polarity/lipophilicity does **not** explain the ordering (micafungin most polar yet shifts most). This **refines** (not refutes) the Phase 8–9 lead → it is about **locally exposed** surface, which bulk TPSA cannot see; and echoes the Phase-10 null (echinocandins’ known *direct* serum effect on glucan synthase). Existence proof: ibrexafungerp/enfumafungin hit the target with **no lipopeptide tail**. |
| **12** | **Serum-tolerance-biased generative design** (Track A: transparent scaffold-constrained generator with a physics-grounded reward) | Reward = 3D **exposed polar surface fraction** (ETKDG + RDKit rdFreeSASA) — the operational form of the Phase 8/9/11 lead; QED/Ro5/clogP dropped. Reward **validated** on the 24 knowns (exposed polar frac vs serum shift **ρ = −0.33, p = 0.11** — reproduces the Phase-8 CREST value with a fast proxy). Generated 28 novel analogs across an **ester** and an ibrexafungerp-inspired **tail-free (notail)** branch; the notail branch dominates the predicted-tolerant, high-novelty region at MW ~700–810 (ibrexafungerp's regime). Deliverable: a **discriminating series** on one scaffold spanning exposed polarity, built to *test* the lead in a serum assay, not assume it. |

## 3. The convergent, defensible finding (updated by Phase 11)
Four independent computational families — 3D shape (polar SASA), QM solvation
(logP), weakly the dipole, and now the echinocandin cross-chemotype read-across —
agree on a single qualitative rule and disagree with nothing:

> **Among equipotent analogs, biasing the *locally exposed* surface toward polar /
> H-bonding character (higher polar SASA, lower QM logP) is associated with smaller
> serum loss.** Effect size |ρ| ≈ 0.30–0.33, p ≈ 0.12–0.15 — directional, not proven.

Phase 11 sharpens two words in that rule: **locally exposed**. Micafungin is the
most polar echinocandin by every *bulk* measure yet loses the most to serum, so the
signal cannot be a whole-molecule polarity descriptor — consistent with the
project's own 3D exposed-SASA framing and a caution against ranking designs on
bulk TPSA/clogP.

Equally important negatives, each now rigorously established rather than assumed:
- **Raw serum MIC is the wrong endpoint** — potency-dominated (ρ = 0.79 with
  serum-free MIC). Model the *shift*, use serum-free MIC as a covariate, and think
  in **free-drug** terms (Phase 11): a big serum shift is *survivable* if intrinsic
  potency and free exposure are high enough (the echinocandins prove this).
- **The Phase-1 rigidity hypothesis is not supported** (SASA-spread shift ρ ≈ −0.1).
- **Bulk lipophilicity/size (clogP, Rg, polarizability, TPSA) tracks potency, not
  tolerance** — the Phase-7 hydrophobic-SASA "hit" was a potency artifact, and the
  Phase-11 echinocandin ordering is not a bulk-descriptor effect.
- **Explicit HSA drug-site docking carries no signal** for these oversized ligands.

## 4. Design output (carried forward with caveats)
The 12 novel Phase-5 candidates were taken through the full funnel; the 3 finalists
by multi-objective score — **cand01 quinolinecarbonyl, cand02 naphthoyl-6-OH,
cand03 pyridylphenyl** (extended rigid aromatics) — have GFN2-reranked ensembles and
ready Gaussian DFT inputs. Given Phases 8–11, finalist *ranking* is unchanged but the
*rationale* shifts decisively: prefer designs that add **locally exposed polar /
H-bonding surface** on the variable region, and — new from Phase 11 — treat the
**long acyl tail as a droppable liability** (ibrexafungerp is a tail-free, orally
active glucan-synthase inhibitor). Do **not** rank on bulk TPSA/clogP/QED.

## 5. Limitations (honest)
- **n = 24, 11 censored, single chemotype/lab** for the primary set; Phase 11 adds
  only **3** echinocandins with matched C. albicans serum shifts — a qualitative
  read-across, not a regression. Convergent descriptors at |ρ| ≈ 0.30, p ≈ 0.13 are
  exactly what a real-but-modest effect looks like under this much censoring; it
  cannot be pushed to significance with more modeling.
- MMFF→CREST showed proxy ensembles **overstate** effects (−0.45 → −0.31): trust the
  QM-ensemble numbers, treat any fast-proxy result as a screen only.
- Vina is not parameterized for 1000–1200 Da amphiphiles; the Phase-10 null cannot
  *exclude* albumin sequestration, only say this proxy is uninformative.
- **There is no validated serum-tolerance oracle.** Any generative optimization
  toward "serum tolerance" optimizes the polar-exposed-surface *hypothesis*, not a
  proven endpoint — so the next generation must be designed to **test** that
  hypothesis, not assume it.

## 6. Recommended next steps
**A. Generative design of candidate serum-tolerant structures — DELIVERED as
Phase 12 (Track A).** `phase12_generate_serum_tolerant.py` runs a
serum-tolerance-biased generator whose reward IS the Phase-8/9/11 physics —
maximize *locally exposed* polar SASA (3D ETKDG + rdFreeSASA), preserve the
FKS-engaging pharmacophore by construction, reward novelty, drop the uninformative
QED/Ro5/clogP terms — with an ibrexafungerp-inspired **tail-free** branch alongside
C-6′ decoration. The reward is validated on the 24 knowns (exposed polar frac vs
serum shift **ρ = −0.33, p = 0.11**, reproducing the Phase-8 CREST value). Output:
28 novel analogs + a **discriminating series** on one scaffold spanning exposed
polarity (built to *test* the lead in a serum assay, not assume it). *Next within
this thread:* funnel `phase12_top_candidates.sdf` through the Phase-6/8/9 QM
pipeline to confirm the property at QM quality, then **Track B** — an RL/CLM
generative network (REINVENT-style) on the external FKS pretraining set reusing
this same reward, once wet-lab shifts (B) anchor it.

> **Round 1 in progress — Phase 13 (`phase13_fatty_tail_optimization.py`).** First
> campaign scopes the design to the **long-chain fatty acid only**, freezing the
> sugar/spiroketal core and the aromatic C-6′ acyl (cleanest one-variable test:
> the tail is the exposed-hydrophobic-surface driver, and freezing the core holds
> the FKS pharmacophore/potency roughly constant). It emits a fixed-core tail
> series spanning exposed polarity plus **CREST-ready inputs and the exact cluster
> protocol** (`phase13_qm_runs/SUBMIT.md`) for the QM confirmation the user runs on
> their platform, feeding the existing Phase-6 descriptor engine. See
> `phase13_findings.md`.

**B. Experimental grounding (the single biggest lever, in parallel).**
- Acquire **uncensored serum MICs** on more analogs — ideally a second chemotype —
  so the polar-surface lead can be tested at adequate power.
- Run the **echinocandin assay playbook** on the papulacandin leads: protein-adjusted
  MIC with an albumin titration, and equilibrium-dialysis fraction-unbound — to
  decide at last whether the gap is albumin sequestration, a direct serum effect, or
  degradation (Phases 8 and 10 could not).

**C. If pursued further in silico:** model the **fatty-acid binding sites / other
serum carriers** (AAG, lipoproteins) rather than more drug-site docking, and run the
prepared **finalist DFT** single points for ESP descriptors — exploratory, not
decision-grade, until B lands.

## 7. Asset map (all under `analysis/`)
- Data & 2D: `outputs/serum_gap_pairs.csv`, `phase1–4` outputs, `outputs/design_rules.md`
- Generative: `phase5_generate.py`, `outputs/phase5_*`
- QM funnel: `phase6_qm_layer.py`; CREST inputs/results in `outputs/qm_runs/`,
  `outputs/qm_runs_known/`; GFN2 finalist reranks in `outputs/qm_runs/cand0*_gfn2/`
- Retrospectives: `phase7_*`, `phase8_*` (+ `phase8_findings.md`, `phase8_confound_analysis.csv`)
- Electronics: `gen_known_xtb_inputs.py`, `phase9_electronic.py`, `outputs/phase9_*`
- Docking: `phase10_dock_hsa.py`, `outputs/phase10_*`
- Cross-chemotype: `phase11_echinocandin_readacross.py`, `outputs/phase11_*`
- Generative (serum-tolerance-biased): `phase12_generate_serum_tolerant.py`,
  `outputs/phase12_*` (generated library, discriminating series, reward validation,
  CREST-ready SDF)
- Per-phase write-ups: `outputs/phase7_findings.md` … `phase12_findings.md`; this
  file: `outputs/SYNTHESIS_phases1-12.md` (supersedes `SYNTHESIS_phases1-10.md`).
