# Plan — AI-guided synthesis prioritization for round-1 serum data

*Feasibility analysis of the proposal "use algorithms to build a descriptor SAR
and tell me which compounds to synthesize," and the recommended refinement.
Grounded in the data as of Stages 0–1c.*

## The problem this solves

Stage 1c's gate did **not** pass: the serum-tolerance labels are one congeneric
fusacandin series, the signal is ~90% potency, and no oracle can be validated for
new chemistry. The binding constraint is **data** (guardrail #8) — specifically,
**scaffold-diverse serum data**. Papulacandins are natural products and expensive
to make, so the question is which compounds to synthesize first.

## What the data says (grounding)

- **Potency labels:** 51 compounds have serum-free MIC — enough for a potency SAR.
- **Serum labels:** 24 compounds, but 23 are one fusacandin series → the gate gap.
- **Diversity exists and is cheap:** the full library spans 34 generic scaffolds;
  the **affordable makeable pool** (short aryl-glycoside route + semisynthetic
  single-position analogs) is **44 compounds / 16 scaffolds**, clogP −1.6→13.5,
  rotatable bonds 4→35, MW 284→1180. **23 of them have no serum data yet.**

## Feasibility verdict

| element of the proposal | verdict |
|---|---|
| descriptor SAR for **potency** | ✅ feasible & useful (Stage 1a) |
| descriptor SAR for **serum tolerance** (predictive) | ❌ not from current data (repeats the Stage 1c failure) |
| "prioritize compounds to synthesize" | ✅ right question — but change the objective |
| involve AI + a generative model | ✅ achievable honestly (below) |

**The trap:** prioritizing the *predicted-most-active* optimizes potency — and this
class's failure is that potent compounds die in serum. That re-selects the known
dead ends. Flip the objective: prioritize for **information about serum tolerance
per synthesis dollar**, using potency as a *filter*, not the selector.

## Recommended plan (reframed, still AI + generative)

1. **Potency oracle (Stage 1a).** Structure-based + descriptor QSAR on the 51
   serum-free MIC labels, scaffold-split validated. Role in the loop: keep only
   candidates predicted **active enough serum-free** — a compound dead in broth
   yields an uninformative serum label.
2. **Generative scaffold-constrained decoration.** Not de-novo whole molecules;
   decorate the FIXED core at the modifiable positions (C-6′ aromatic ester, acyl
   tail, sugar hydroxyls) to enumerate a virtual **makeable** library. This is the
   generative-model component, and it respects the pharmacophore.
3. **Synthesizability + cost filter.** Cross-check `synthesis_feasibility.csv`;
   prefer the cheap aryl-glycoside / single-position semisynthetic routes over
   total synthesis or fermentation. Maximize information *per dollar*.
4. **Active-learning panel selection (the core AI step).** Choose a small panel
   that is (a) makeable & cheap, (b) scaffold- and physicochemically **diverse**
   (fix the gate gap), (c) **potency-matched but differing in the hypothesized
   serum drivers** — tail rigidity, aromatic-ester size, logP — so results
   *decouple* serum tolerance from potency, and (d) where the oracles are most
   uncertain.
5. **Measure → retrain → regenerate.** Serum-shift MIC ± albumin and
   equilibrium-dialysis fu on the panel → retrain oracles → regenerate. The
   design–make–test loop (START_HERE Stage 4). Round 1 both tests candidates and
   *manufactures the scaffold-diverse dataset the gate needs.*

## Testable hypotheses to build the panel around

- **Aromatic-ester lever (from Stage 1c):** larger biphenyl/terphenyl/naphthoyl
  C-6′ esters associate with serum retention — but confounded with potency.
  Design potency-matched pairs that vary ester size to test it cleanly.
- **Tail-shape lever (guardrail #7):** de-rigidify the acyl tail at constant
  length (saturated/branched vs the native rigid polyene) — untested.

## Why this satisfies "AI + generative" honestly

Supervised potency oracle (AI) + scaffold-constrained generative decoration
(generative) + active-learning experimental design (AI) — none of which pretends
to predict an endpoint we cannot yet predict. The novelty is using AI to design
the *most informative makeable experiment*, not to guess the winner.

## First concrete build options

- **A. Stage 1a potency oracle** — the filter everything else needs.
- **B. Round-1 panel selector** — rank the 44 makeable compounds (23 unlabelled)
  by information × diversity × makeability, output a prioritized synthesis list.
- **C. Generative decoration prototype** — enumerate makeable core-decorations and
  score them, producing the candidate pool for B.
