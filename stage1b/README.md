# Stage 1b — free-fraction / plasma-protein-binding (PPB) oracle

START_HERE (§3) calls this "the serum-tolerance workhorse … genuinely learnable":
only unbound drug acts, so the free fraction **fu** is a core factor of serum
tolerance (`serum tolerance ≈ potency × fu × stability`). This stage predicts
**log10(% unbound)** (fu = 10^pred / 100).

## Run

```bash
python stage1b/train_free_fraction_oracle.py      # trains + predicts + validates
python -m pytest stage1b/test_stage1b.py -q        # 6 smoke tests
```

Requires `scikit-learn`, `rdkit`, `scipy`, `numpy`. The training CSV is cached at
`stage1b/data/ADME_public_set_3521.csv` (fetched once from the Genentech
Computational-ADME GitHub repo; the script re-downloads if absent).

## Method

- **Training data:** public Fang et al. 2023 *Computational-ADME* set — **194
  human-PPB drugs** with SMILES and `LOG % unbound` (MW 151–666).
- **Features:** 14 physicochemical RDKit descriptors (size, logP, TPSA, H-bonding,
  shape) computed by one shared `featurize.py` for training drugs, papulacandins
  **and** echinocandins — so every molecule is on the same basis.
- **Model:** random forest; uncertainty = spread across trees.
- **Validation (guardrail #2):** a **Murcko scaffold split**, not a random split.
- **Applicability domain (guardrail #4):** standardized-descriptor kNN distance
  (≤ 95th-percentile of train-to-train distance) **and** a training-MW-range gate.

## Results

| validation | R² | MAE (log) | Spearman |
|------------|----|-----------|----------|
| scaffold-split hold-out (honest) | ~0.39 | ~0.60 | ~0.66 |
| random 5-fold (optimistic) | ~0.45 | — | — |

PPB is learnable in-domain, at modest accuracy given 194 compounds. The
scaffold–random gap is the usual optimism of random splits — trust the scaffold
number.

**Anchor check — anidulafungin:** predicted ~96% bound vs ~88% observed in-repo
(literature ≈ 99%); right regime, but flagged **out-of-domain** (a large
echinocandin) — read with the AD flag.

## The decisive finding (why this matters)

| chemotype | in applicability domain |
|-----------|-------------------------|
| **papulacandin** | **2 / 137 (1.5%)** |
| echinocandin / FKS | 91 / 321 (28%) |

The papulacandin glycolipids (median MW ~930) sit almost entirely **outside** the
training drugs' range (≤666 Da; the typical query's kNN distance is ~13 vs a
threshold of ~3). A PPB model built on normal drugs must **extrapolate** for them,
and random-forest error bars collapse toward the training mean under
extrapolation — which is why the distance-based AD flag, not the tree spread, is
the trustworthy guide here.

### Is the AD "successful" because more echinocandins pass it?

No — that reading is a trap; see `outputs/applicability_domain_analysis.md`
(`analyze_applicability_domain.py`). The echinocandin/papulacandin gap is a
**molecular-size gradient**, not chemotype discrimination: the AD admits small
molecules of either class (MW<500: 73% of echinocandins, 22% of papulacandins)
and rejects large ones of either class (MW>700: 0% of both). Every **approved**
echinocandin — caspofungin, anidulafungin, micafungin, cilofungin, ibrexafungerp,
the serum-tolerant drugs we'd use as the *teacher* — is out-of-domain, exactly
like the papulacandins. The in-domain echinocandins are small synthetic analogs,
not the serum-tolerance exemplars.

The AD *does* work, but the evidence is the error test, not the head-count: on
held-out labelled data, in-domain MAE = 0.47 vs out-of-domain 0.80 log-units, MAE
rises monotonically across distance quartiles (0.32→0.43→0.49→0.69), and
Spearman(distance, |error|) = +0.33. So the AD is a valid **reliability flag** —
and it says target compounds *and* the large-molecule teachers are both
out-of-domain for a drug-trained PPB oracle. (This does not threaten the strategy:
Stage 1c's serum transfer-learning uses the echinocandins' own measured serum /
in-vivo data, not this oracle.)

This is exactly **guardrail #4** (off-the-shelf QSAR does not transfer to this
bRo5 chemotype) made quantitative, and it operationalizes **guardrails #5 and #8**:
for the real design targets, free fraction must be **measured** (equilibrium
dialysis), not predicted. Use this oracle to (a) **rank in-domain** analogs and
(b) **flag** which designs even fall in a modellable region — never to assign an
fu to a full-size papulacandin.

## Outputs (`stage1b/outputs/`)

| file | what it is |
|------|------------|
| `free_fraction_predictions.csv` | per-compound fu + tree-SD uncertainty + AD flag, all papulacandins & echinocandins |
| `applicability_domain_report.csv` | in-domain fractions by chemotype |
| `ppb_oracle_metrics.md` | validation, anchor check, and the AD finding |
| `applicability_domain_analysis.md` | why the echinocandin in-domain rate is a size effect, and the error-vs-distance validation of the AD |

## Next

- **Stage 1a** potency oracle (structure-based + pooled-MIC QSAR).
- **Stage 1c** transfer serum-tolerance oracle, gated on scaffold-split validation.
- Feed this fu prediction (in-domain only) as a covariate into 1c, and use the
  AD flag to decide which designs need **experimental** fu before trusting them.
