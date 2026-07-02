# Phase 11 — echinocandin read-across for serum tolerance

*Question:* the papulacandin serum-gap lead ("expose polar not hydrophobic
surface") never cleared significance on our single chemotype (n=24). Do
the echinocandins — an independent cyclic-lipopeptide glucan-synthase inhibitor
class with the same long-tail serum liability — corroborate it, and what does
their clinical success teach us? All echinocandin numbers below come from data
already curated in `external/data/external/fks_inhibitors/`.

## 1. Echinocandins suffer the SAME serum shift — and it is large
Serum-shift fold (MIC in 50% serum / MIC in serum-free medium), the identical
endpoint direction as our `serum_shift_fold`:

- **Caspofungin (x2), Anidulafungin (x16), Micafungin (x64)** (C. albicans; robust across multiple isolates/species).
- Across all species the shifts reach ~256-1024x for micafungin/anidulafungin.

So a big serum shift is **not, by itself, disqualifying**: micafungin and
anidulafungin lose 1-2 orders of magnitude of in vitro potency to serum yet are
front-line IV antifungals. This reframes the whole project's target — the goal
is not "zero serum shift" but "enough free-drug exposure at the target."

## 2. The free-drug hypothesis is the right lens (PPB / Fu)
Echinocandins are ~96% (caspofungin), ~99% (anidulafungin) and ~99.8%
(micafungin) plasma-protein bound (literature; our corpus carries anidulafungin
PPB [84.0]%
and Fu ~0.01%). Clinically they are dosed to a **free-drug** AUC/MIC target, and
when normalized to unbound drug the PK/PD targets converge across the class.
Read-across for us: always model the *shift* / free fraction (Phase 8's
conclusion), and — the actionable part — **measure protein binding directly on
our leads** rather than inferring it from docking (Phase 10 could not).

## 3. Bulk 2D descriptors do NOT explain the echinocandin ordering (honest null)
The most polar echinocandin by every bulk measure — **micafungin**
(TPSA 510,
MolLogP -3.9) —
has the **largest** serum shift, while **caspofungin** has the smallest. Within
the papulacandins the same bulk axis is weak/uninformative
(MolLogP Spearman rho=+0.15). So whole-molecule polarity/lipophilicity
is **not** the read-across variable.

Two non-exclusive readings, both consistent with earlier phases:
1. The papulacandin lead was about **locally exposed** surface (Boltzmann polar
   SASA / QM logP, Phases 8-9), which bulk TPSA cannot see. A bulk null here does
   **not** refute the 3D lead — it says "don't shortcut it with 2D."
2. The echinocandin serum effect has a documented component of **direct
   desensitization of glucan synthase by serum**, not pure albumin sequestration
   (the effect does not track protein-binding rank). This mirrors our Phase 10
   HSA-docking null and warns that "reduce albumin affinity" may be the wrong
   single objective.

## 4. Existence proof: the lipophilic tail is droppable
Enfumafungin, Ibrexafungerp
inhibit the same target with **no cyclic-lipopeptide fatty-acyl tail**
(ibrexafungerp: MW ~730, and it is an *orally* bioavailable marketed drug with
good tissue penetration). Since synthetic accessibility is off the table for now,
this legitimizes an aggressive design branch: **truncate / replace the long acyl
chain** of the papulacandin scaffold while keeping the aryl-C-glycoside /
spiroketal pharmacophore, rather than only decorating the C-6' ester.

## 5. What Phase 11 changes for the project
- **Endpoint:** adopt the free-drug framing; the deliverable is a serum-tolerant
  *free* exposure, not a serum-invariant MIC.
- **Design (Phase 5 next iteration):** keep scoring exposed polar surface
  (QM logP / polar SASA), and add a tail-truncation/replacement branch inspired
  by the fungerps. Do **not** rank on bulk TPSA/clogP — Phase 11 shows it does
  not separate serum tolerance across chemotypes.
- **Experiment (highest value):** run the echinocandin assay playbook on our
  papulacandin leads — protein-adjusted MIC with an albumin titration, and
  equilibrium-dialysis fraction-unbound — to decide, at last, whether the serum
  gap is albumin sequestration, a direct serum effect, or degradation. This is
  the measurement Phases 8 and 10 could not make computationally.

## Caveats
n=3 echinocandins with matched C. albicans serum shifts is a qualitative
read-across, not a regression; different pharmacophore, species mix, and assay
labs than the papulacandin set. The value is the **direction and the reframing**,
plus a reproducible harmonized dataset (`phase11_crosschemotype.csv`) to extend
as more +/- serum data is curated.

## Sources
- Paderu et al., *Antimicrob. Agents Chemother.* 2007 — Effects of serum on in
  vitro susceptibility testing of echinocandins.
- Odds/Gumbo et al. — serum differentially alters echinocandin antifungal
  activity; free-drug AUC/MIC targets.
- Davis et al., *J. Fungi* 2021 — Ibrexafungerp: first-in-class oral triterpenoid
  glucan synthase inhibitor.
- In-repo: `external/data/external/fks_inhibitors/source_exports/*` (ChEMBL/PubChem
  curated).
