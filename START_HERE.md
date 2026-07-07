# START HERE — Papulacandin serum-tolerant antifungal project (fresh restart)

*Read this first. It is the complete handoff for a new working session. The repo
was deliberately reset to a clean slate: only the raw **databases** are kept. The
previous exploratory analysis pipeline was removed from the working tree (it
remains in git history — see the last section — if you ever need to consult it).*

---

## 1. The goal (one paragraph)

Design **novel small-molecule antifungal candidates that stay active in blood
serum**, in the **papulacandin / fusacandin** structural class, which inhibit the
fungal enzyme **β-1,3-glucan synthase (FKS1)** — the machine that builds the
fungal cell wall (humans have no cell wall, so it is a low-toxicity target). The
class is potent in a test tube but **loses activity in serum**, which is why it
never became a drug. The objective is to use **AI + a generative model** to
propose candidates that close that "serum gap," and to do it as a disciplined
**design → make → test → learn** loop, not a one-shot prediction.

The overarching program context is in `research_goal.html` (in Chinese).

---

## 2. What is in this repo right now (your only assets)

Two curated databases and their own documentation. **Read the READMEs inside each
folder first — they describe every table.**

### `curated/` — the papulacandin (target) chemotype
- `curated/README_MAIN_CURATED_DATA_ONLY_FOR_CLAUDE_CODE.md` — start here.
- `curated/core_tables/compounds_master.csv` — **138 compounds** (SMILES, InChIKey, descriptors).
- `curated/core_tables/activity_table.csv` — **1,042 activity records** (whole-cell MIC with/without serum, etc.).
- `curated/core_tables/enzyme_assays.csv` — enzyme-level (target) assays.
- `curated/core_tables/synthesis_feasibility.csv` — **which analogs are makeable** (use to prioritize designs).
- `curated/core_tables/structure_quality_notes.tsv`, `structure_qc_manual_review.csv` — structure QC.
- `curated/after_analysis_reference/sar_annotations.csv` — literature SAR (reference only).
- `curated/structures/{sdf,mol,cdx}/` — 3-D/2-D structures.
- **The scarce but precious label:** only ~**24 compounds** have *matched* serum-free **and** serum MIC — the direct "serum tolerance" data.

### `external/` — the data-rich related chemotypes (your "teacher")
- `external/README_DATA_ONLY_CLEAN_FOR_CLAUDE_CODE.md` — start here.
- `external/data/external/fks_inhibitors/source_exports/` — **echinocandins and other FKS/glucan-synthase inhibitors** (caspofungin, anidulafungin, micafungin, ibrexafungerp, …) from ChEMBL/PubChem: thousands of MICs, **~279 serum-context rows, plasma-protein-binding (PPB) and free-fraction (Fu) values**, across many species/strains. Look at the `serum_or_protein` / `assay_context` columns and the `PPB` / `Fu` endpoint types.
- `external/data/processed/pretraining_v0_1/` — a **pretraining corpus** (FKS-relevant chemistry) for a generative/predictive model.
- `external/data/processed/external_fks_model_ready_v0_1/` — model-ready descriptor matrices.

---

## 3. The methodology to follow (already decided — this is the plan)

**Do not try to train a "serum-MIC predictor" from 24 compounds — it will fail.**
Instead, **factor serum tolerance into learnable pieces** and use the echinocandin
data as a teacher:

> serum tolerance ≈ **intrinsic potency** × **free fraction (fᵤ)** × **stability**
> (the "free-drug" principle: only unbound, un-degraded drug acts).

### Stage 0 — Unify data & fix the endpoint (before any model)
- Merge `curated` + `external` on a common endpoint: the **serum shift = serum MIC ÷ serum-free MIC** (a *potency-independent* measure). Pool papulacandins **and** echinocandins.
- Always keep **intrinsic potency as a covariate** (it is the dominant confound — see guardrails).
- Record the **applicability domain** of every dataset (what chemistry it covers).

### Stage 1 — Build & VALIDATE the oracles (the crux; do this before generating)
- **Potency oracle:** structure-based scoring against the **glucan-synthase target structure** *(a cryo-EM structure of fungal β-1,3-glucan synthase with caspofungin was reported in Nature, 2026 — locate the PDB and verify)*, plus a QSAR on the pooled MIC data.
- **Free-fraction / protein-binding oracle:** train on the echinocandin **PPB/Fu** rows + public ADMET protein-binding datasets. *This is the serum-tolerance workhorse and it is genuinely learnable.*
- **Serum-shift oracle by transfer learning:** **pretrain on the ~279 echinocandin serum-shift rows, then fine-tune on the 24 papulacandin pairs.**
- **Physics features as a backstop** where data is thin: 3-D exposed polar vs greasy surface (conformer-ensemble SASA), computed logP — but only after retrospective validation (see guardrails).
- Every oracle ships with **uncertainty** and an **applicability-domain flag**.
- **Gate:** do not proceed to generation until the combined serum-tolerance oracle is validated on a scaffold-split hold-out.

### Stage 2 — Generative model
- **Pretrain** on `pretraining_v0_1` + broad chemistry.
- Because these are large "beyond-Rule-of-5" glycolipids with a fixed target-engaging core and a few modifiable positions (the fatty tail, the C-6′ aromatic ester, sugar hydroxyls), use **scaffold-constrained / fragment-decoration generation on the fixed core** — *not* de-novo whole-molecule generation. More tractable, more synthesizable, keeps the pharmacophore intact.
- **Optimize** against the Stage-1 multi-objective score (potency × fᵤ × stability × novelty) via RL (e.g. REINVENT-style) or GFlowNets, with **synthesizability** (cross-checked against `synthesis_feasibility.csv`) and **applicability-domain** penalties in the reward.

### Stage 3 — In-silico triage funnel
Validity/ADMET/applicability filters → multi-objective ranking → diversity clustering → **expensive physics/QM re-scoring on the top few dozen** (to catch oracle artifacts) → a **small, diverse, synthesizable, informative** shortlist.

### Stage 4 — Design–make–test–learn loop (the actual engine)
Serum labels are scarce, so round 1 also *generates* data: synthesize a **diverse discriminating series** spanning the physicochemical axes → measure serum-relevant assays **directly** (protein-adjusted MIC ± albumin; equilibrium-dialysis fᵤ) → **retrain** oracles → **regenerate**. Iterate. The generator is one gear in this loop.

---

## 4. Concrete first tasks for this new session

1. **Read** both database READMEs and profile the tables (columns, row counts, missingness). Confirm the ~24 matched serum pairs and the echinocandin serum/PPB/Fu rows.
2. **Stage 0:** build the unified, potency-covariate-aware serum-shift table (papulacandin + echinocandin pooled).
3. **Stage 1a:** locate/verify the glucan-synthase target structure; stand up a potency scorer.
4. **Stage 1b:** train the **free-fraction / PPB** oracle (the most data-rich, highest-value piece).
5. **Stage 1c:** the transfer-learning serum-shift oracle (echinocandin → papulacandin); validate on a scaffold-split hold-out. **Gate before Stage 2.**
6. Only then build the scaffold-constrained generator and run the funnel.

Set up the environment first (see §6) and re-establish reproducibility scaffolding (pinned deps, a couple of smoke tests) as you build — the previous project's was removed with the old pipeline.

---

## 5. Hard-won guardrails — decisions already learned; **do not re-derive these**

The previous exploratory work (now removed) established these the hard way. Honor
them from day one to save months:

1. **Model the serum SHIFT, never the raw serum MIC.** Raw serum MIC is dominated by intrinsic potency (potency alone tracked serum MIC at Spearman ρ≈0.79). Any descriptor correlated with raw serum MIC is probably just tracking potency. **Always control for potency (partial correlation / covariate).**
2. **Cheap methods overstate.** A fast 3-D descriptor looked strong (ρ≈−0.45) but shrank to non-significant at accurate quality and collapsed to ρ≈0.02 once potency was controlled. **Validate any descriptor/oracle at accurate quality, on a scaffold-split hold-out, with an applicability-domain limit — before trusting or optimizing it.**
3. **Do not design before the oracle is validated.** The prior project generated molecules first and validated the scoring later; that caused wasted rounds. **Validate → then generate.**
4. **An off-the-shelf external QSAR does NOT transfer to this chemotype.** A classifier trained on external FKS inhibitors scored all papulacandins as "inactive" (out-of-domain). Use **structure-based** potency + **transfer learning**, not a naive external QSAR, and always check applicability domain.
5. **Docking to serum albumin (HSA) was uninformative** for these large amphiphiles (Vina is not calibrated for ~1,000–1,200 Da). **Measure free fraction experimentally**; don't rely on docking to predict protein binding here.
6. **2-D drug-likeness (QED/Rule-of-5) is uninformative** for this bRo5 class (QED ≈ 0.01–0.03; nothing passes Ro5). Don't use it to rank designs.
7. **The tail is the key lever, and shape matters more than greasiness.** Approved-drug (echinocandin) evidence: the least serum-affected drug (caspofungin) has a **flexible/branched** tail; the worst have **rigid/aromatic** tails. The papulacandin native tail is a **rigid polyene** (unfavorable). The tail is *required for potency* (do not simply remove it) but is *modifiable* (ibrexafungerp works with no lipopeptide tail). A promising, untested design hypothesis: **de-rigidify the tail at constant length** (e.g. saturated/branched instead of the rigid polyene).
8. **The binding constraint is DATA, not methods.** With ~24 labels, no amount of computation yields a validated predictor. **Plan the wet-lab loop from the start**; computation's job is to *aim* the experiments.

*(These are guidance, not results to reproduce. Treat them as priors to test, not dogma — but do not waste effort rediscovering points 1–6, which are well established.)*

---

## 6. Environment & tooling (suggested)

Python 3.11+, with (minimum): `rdkit`, `pandas`, `numpy`, `scipy`, `scikit-learn`,
`matplotlib`. For generation/QM later: a generative framework (e.g. REINVENT or a
graph/diffusion model), and optionally `xtb`/`crest` for physics descriptors and a
docking engine for structure-based potency. Pin versions in a `requirements.txt`
you create, and keep a couple of smoke tests as you build.

---

## 7. Repo & git logistics

- **Kept:** `curated/`, `external/` (the databases), `research_goal.html`, this file, `README.md`.
- **Removed from the working tree** (recoverable from git history): the old `analysis/` pipeline (Phases 0–14: scripts, QM runs, figures, decks, tutorials, summaries), `tests/`, `.github/`, and the old scaffolding/provenance/requirements files.
- **To recover anything old:** `git log --oneline` and `git show <commit>:<path>`, or `git checkout <commit> -- <path>`. Nothing is permanently lost.
- **Branch:** this clean state is on the working branch. To make it the state a fresh clone sees, **merge it into the default branch (`main`) before starting the new session**, or open the new session on this branch.
