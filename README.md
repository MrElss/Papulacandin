# Papulacandin — serum-tolerant antifungal design (AI + generative, fresh start)

Goal: use **AI + a generative model** to design novel **papulacandin/fusacandin**
antifungal candidates (β-1,3-glucan / **FKS1** synthase inhibitors) that **retain
activity in blood serum** — the flaw that stopped this class from becoming drugs.

This repository has been reset to a clean slate for a fresh start. It currently
contains only the **raw databases** needed to begin:

- **`curated/`** — the papulacandin chemotype: 138 compounds, 1,042 activity
  records (incl. the scarce ~24 matched serum-free/serum MIC pairs), structures,
  synthesis-feasibility. See its internal README.
- **`external/`** — the data-rich related chemotypes (echinocandins & other FKS
  inhibitors) with serum-shift, plasma-protein-binding (PPB) and free-fraction
  data, plus a pretraining corpus. See its internal README.
- `research_goal.html` — the overarching program objectives (in Chinese).

## ▶ New here? Read [`START_HERE.md`](START_HERE.md) first.

It is the full handoff: the goal, what each database contains, the decided
methodology (factor serum tolerance into learnable oracles → validate → generate
→ design–make–test loop), the concrete first tasks, and the **hard-won guardrails**
to avoid repeating known dead ends.

*The previous exploratory analysis (Phases 0–14) was removed from the working tree
but remains in git history if ever needed.*
