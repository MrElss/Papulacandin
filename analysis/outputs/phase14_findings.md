# Phase 14 — which tail to synthesize first (echinocandin-guided)

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
- **C16_0_palmitoyl** [MAKE FIRST] — C=16, C=C=0, rot=14, tail clogP 5.54, whole clogP 3.32
    - *why:* saturated, fully flexible, same C16 length -> holds lipophilicity/potency, removes the rigidity/conjugation echinocandin SAR blames for serum binding. The single cleanest test of the hypothesis.
    - *synthesis:* trivial: palmitic acid (commodity), one esterification.
- **branched_sat_caspofungin_like** [MAKE FIRST] — C=16, C=C=0, rot=13, tail clogP 5.4, whole clogP 3.17
    - *why:* branched saturated ~C16 — the closest structural analog to caspofungin's tail, the marketed echinocandin with the SMALLEST serum shift.
    - *synthesis:* easy: branched (anteiso/iso) fatty acid, one esterification.
- **C16_1_palmitoleoyl** [second] — C=16, C=C=1, rot=13, tail clogP 5.32, whole clogP 3.09
    - *why:* one cis double bond (kinked, fluid) — intermediate rigidity control that separates 'conjugation/rigidity' from 'saturation': if it behaves like palmitoyl, the culprit is the CONJUGATION, not merely the double bonds.
    - *synthesis:* easy: palmitoleic/oleic acid, one esterification.
- **native_C16_polyene** [reference] — C=16, C=C=4, rot=10, tail clogP 3.62, whole clogP 1.39
    - *why:* rigid conjugated tetraene (4 C=C); the anidulafungin/micafungin-like high-serum-binding regime — this is what we are trying to beat.
    - *synthesis:* you likely already have the parent; hard (conjugated-polyene synthesis).

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
