# Phase 4 — Serum-tolerance design rules (hypotheses to test)

Two INDEPENDENT lines of evidence converge on actionable design directions.
Both are hypothesis-generating; neither is a validated quantitative model.

## Axis 1 — Descriptor SAR (this project, Phase 1/3B; Fusacandin-A C-6' series)
Within the only modelable series, serum activity improves with:
  * mw_exact: rho=-0.59 vs serum MIC -> favor higher mw_exact
  * fsp3: rho=+0.58 vs serum MIC -> favor lower fsp3
  * hba: rho=-0.39 vs serum MIC -> favor higher hba
  * tpsa: rho=-0.35 vs serum MIC -> favor higher tpsa
Plain reading: **rigid, extended AROMATIC C-6' acyl groups (biphenyl / naphthoyl
types) retain serum activity better than flexible aliphatic chains.** clogP alone
is NOT the lever.

## Axis 2 — Literature SAR annotations (curated/after_analysis_reference)
Independent, paper-derived observations on the Papulacandin-B core point a
DIFFERENT structural direction for in-vivo/serum translation:
  * Polar 10-O ethers can improve in-vivo activity without losing Candida MIC.
  * Selected 11-aminoacyl / cationic substitutions improve in-vivo activity
    while keeping useful MIC.
Plain reading: on the Pap-B core, **adding polar / cationic handles** (not just
lipophilic bulk) is the literature-supported route to better exposure.

## Synthesis recommendations (concrete next analogs)
1. Best EXISTING serum lead to advance now (goal #1): the biphenyl/naphthoyl
   Fusacandin-A C-6' esters (serum MIC ~11-25 ug/mL), all with known 4-step
   semisynthetic routes -> lowest-risk compounds to re-test under serum.
2. New analog hypothesis A (Axis 1): rigid biphenyl C-6' ester bearing a polar
   substituent (e.g. carboxylate / morpholine) to merge rigidity + polarity.
3. New analog hypothesis B (Axis 2 x goal #2): graft a polar/cationic 11-aminoacyl
   handle onto the serum-tolerant biphenyl-ester scaffold -> tests whether the
   two independent axes are additive for serum tolerance.

## Caveats
* All scores/rules are within or near the papulacandin/fusacandin domain only.
* Heuristic design-score vs observed serum MIC sanity check (labeled set, n=24):
  Spearman rho = -0.57 (p = 0.004). Consistent by construction with
  Axis 1; this is NOT independent validation.
