# Papulacandin / FKS1 serum-gap project

A computational medicinal-chemistry study of the papulacandin/fusacandin family
of antifungal natural products (β-1,3-glucan / FKS1 synthase inhibitors),
built around one question:

> **Why do a few papulacandin-class analogs keep their antifungal activity in
> blood serum while structurally near-identical analogs lose it — and can we
> compute a rule to design serum-tolerant analogs?**

**Answer, in one line** (full version in
[`analysis/outputs/SYNTHESIS_phases1-10.md`](analysis/outputs/SYNTHESIS_phases1-10.md)):
serum *potency* is governed by intrinsic potency × molecular size/lipophilicity;
serum *tolerance* (the free→serum MIC shift) has at most a weak, non-significant
association with exposing **polar rather than hydrophobic** surface — a lead that
recurs across independent methods but never clears significance on the available
data. **No computable descriptor or docking model predicts serum tolerance on
this dataset; the binding constraint is the data (n = 24, 11 censored, single
chemotype), not the methods.**

This repository is one plank of a larger program (see
[`research_goal.html`](research_goal.html), in Chinese) to create original
antifungal drug candidates for invasive fungal infections.

## Repository map

| Path | What it holds |
|---|---|
| [`curated/`](curated/) | Clean, human-QC'd source data: `compounds_master.csv` (138 compounds), `activity_table.csv` (1042 MIC/assay records), enzyme assays, synthesis-feasibility notes, and structure files (SDF/MOL/CDX). Packaged "data-only" so a downstream tool can re-derive SAR independently. |
| [`external/`](external/) | External FKS / glucan-synthase inhibitor datasets (source exports, curation notes, processed model-ready + pretraining matrices). |
| [`analysis/`](analysis/) | The pipeline: self-contained, ordered Python scripts `phase1…phase10` writing to `analysis/outputs/`. See [`analysis/README.md`](analysis/README.md) for the phase-by-phase narrative. |
| [`analysis/outputs/`](analysis/outputs/) | Results: CSVs, figures, QM run artifacts, per-phase `*_findings.md`, the cross-phase [`SYNTHESIS_phases1-10.md`](analysis/outputs/SYNTHESIS_phases1-10.md), and slide decks. |
| [`tests/`](tests/) | Smoke tests (curated-data integrity + fast pipeline entry points). |
| [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) | Data sources, curation dates, toolchain versions, endpoint/censoring conventions. |

## The pipeline at a glance

`serum_gap_analysis.py` (dependent variable) → **1** 2D SAR → **2** activity
ladder → **3** external FKS model + within-series serum-tolerance model →
**4** interpretable design rules → **5** scaffold-constrained generative design
(12 novel candidates) → **6** QM funnel (CREST ensembles, in-house Shrake–Rupley
SASA, Boltzmann weighting, Gaussian DFT I/O) → **7–8** retrospective validation
(MMFF proxy → real CREST + confound analysis) → **9** GFN2-xTB electronics →
**10** explicit HSA docking → **11** echinocandin cross-chemotype read-across
(reuses the `external/` FKS corpus to reframe the endpoint via the free-drug
hypothesis and stress-test the polar-surface lead against a second chemotype).

Phases 0–9 run on a normal workstation. Phases that generate raw QM/docking data
need external binaries (CREST, xtb, Gaussian, AutoDock Vina); the Python scripts
otherwise only *parse* those tools' outputs.

## Quickstart

```bash
# 1. Install the core analysis stack (validated versions, see requirements.txt)
python3 -m pip install -r requirements.txt

# 2. Run the fast entry points (read curated/, write analysis/outputs/)
python3 analysis/serum_gap_analysis.py     # builds the matched serum-gap table
python3 analysis/phase1_serum_shift_sar.py # 2D serum-shift SAR

# 3. Run the smoke tests
python3 -m pytest tests/ -q
#   or without pytest:  python3 tests/test_smoke.py
```

Optional docking / slide-deck extras and the external-binary requirements are
listed in [`requirements-optional.txt`](requirements-optional.txt).
The full run order and every output is documented in
[`analysis/README.md`](analysis/README.md).

## Reading order for the science

1. [`analysis/outputs/SYNTHESIS_phases1-10.md`](analysis/outputs/SYNTHESIS_phases1-10.md) — the payoff: the whole story and the honest limitations.
2. `analysis/outputs/serum_gap_summary.txt` — the dependent variable and the serum-tolerant leads.
3. `analysis/outputs/phase8_findings.md` … `phase10_findings.md` — the decisive confound analysis and null results.
4. [`analysis/README.md`](analysis/README.md) — how each phase was run (incl. the real CREST/xtb settings used).

## Limitations (in one breath)

n = 24 matched pairs, 11 serum MICs censored at `>100`, values tied on a few
levels, mostly one chemotype from one lab. Three convergent descriptors at
|ρ| ≈ 0.30 (p ≈ 0.13) is what a real-but-modest effect looks like under this much
censoring — it cannot be pushed to significance with more modeling. The
recommended next steps are **experimental** (uncensored serum MICs on more
analogs; a direct HSA-binding assay), not more computation.
