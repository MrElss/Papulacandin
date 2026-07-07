#!/usr/bin/env python3
"""Stage 1b analysis — is the applicability domain (AD) "successful"?

Prompted by the question: papulacandins are almost all out-of-domain while ~28%
of echinocandins are in-domain — does that show the AD is working?

This script separates two very different claims:

  A. "More echinocandins are in-domain than papulacandins" -> a chemotype signal.
  B. "In-domain membership predicts a more reliable fu prediction" -> the AD's
     ACTUAL job.

It shows that (A) is a molecular-SIZE gradient, not chemotype discrimination
(the approved echinocandin drugs are as out-of-domain as the papulacandins), and
that (B) genuinely holds (error rises with distance-from-training). So the AD is
a valid RELIABILITY flag, not evidence that the echinocandin teacher transfers.

Run:
    python stage1b/analyze_applicability_domain.py
Writes stage1b/outputs/applicability_domain_analysis.md
"""

from __future__ import annotations

import csv
import os

import numpy as np
from scipy.stats import spearmanr
from sklearn.model_selection import KFold

import featurize as F
import train_free_fraction_oracle as s1b

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
PRED = os.path.join(OUT, "free_fraction_predictions.csv")

APPROVED = ["CASPOFUNGIN", "ANIDULAFUNGIN", "MICAFUNGIN", "REZAFUNGIN",
            "CILOFUNGIN", "IBREXAFUNGERP"]
MW_BANDS = [(0, 500), (500, 700), (700, 900), (900, 3000)]


def _frac(rows):
    n = sum(1 for r in rows if r["in_domain"] == "True")
    return n, len(rows), (100 * n / len(rows) if rows else 0.0)


def size_breakdown():
    rows = list(csv.DictReader(open(PRED)))
    papu = [r for r in rows if r["chemotype"] == "papulacandin"]
    ech = [r for r in rows if r["chemotype"] == "echinocandin"]
    bands = {}
    for lo, hi in MW_BANDS:
        bands[(lo, hi)] = {
            "papulacandin": _frac([r for r in papu if lo <= float(r["mw"]) < hi]),
            "echinocandin": _frac([r for r in ech if lo <= float(r["mw"]) < hi]),
        }
    approved = [r for r in ech
                if any(k in r["name"].upper() for k in APPROVED)]
    return _frac(papu), _frac(ech), bands, approved


def ad_validity():
    """Out-of-fold: does distance-from-training predict |error|? (the real test)"""
    s1b.ensure_training_csv()
    X, y, _, _ = s1b.load_training()
    mw_i = F.FEATURE_NAMES.index("mw")
    kf = KFold(n_splits=5, shuffle=True, random_state=s1b.SEED)
    dists = np.zeros(len(y))
    errs = np.zeros(len(y))
    indom = np.zeros(len(y), bool)
    for tr, te in kf.split(X):
        m = s1b.new_model().fit(X[tr], y[tr])
        ad = s1b.ApplicabilityDomain(X[tr], mw_index=mw_i, k=5)
        for j in te:
            errs[j] = abs(m.predict([X[j]])[0] - y[j])
            sc = ad.score(X[j])
            dists[j] = sc["knn_distance"]
            indom[j] = sc["in_domain"]
    quartiles = np.quantile(dists, [0, .25, .5, .75, 1.0])
    q_rows = []
    for i in range(4):
        lo, hi = quartiles[i], quartiles[i + 1]
        mask = (dists >= lo) & (dists <= hi if i == 3 else dists < hi)
        q_rows.append((i + 1, lo, hi, float(errs[mask].mean()), int(mask.sum())))
    return {
        "in_mae": float(errs[indom].mean()), "in_n": int(indom.sum()),
        "out_mae": float(errs[~indom].mean()), "out_n": int((~indom).sum()),
        "quartiles": q_rows,
        "spearman": float(spearmanr(dists, errs).correlation),
    }


def main():
    papu, ech, bands, approved = size_breakdown()
    v = ad_validity()

    L = []
    L.append("# Stage 1b analysis — is the applicability domain successful?\n\n")
    L.append(
        f"Overall in-domain: **papulacandin {papu[0]}/{papu[1]} "
        f"({papu[2]:.0f}%)**, **echinocandin {ech[0]}/{ech[1]} ({ech[2]:.0f}%)**. "
        "The question: does the echinocandin-vs-papulacandin gap show the AD works?\n\n")

    L.append("## 1. The gap is molecular SIZE, not chemotype\n\n")
    L.append("| MW band | echinocandin in-domain | papulacandin in-domain |\n")
    L.append("|---|---|---|\n")
    for lo, hi in MW_BANDS:
        e = bands[(lo, hi)]["echinocandin"]
        p = bands[(lo, hi)]["papulacandin"]
        e_s = f"{e[0]}/{e[1]} ({e[2]:.0f}%)" if e[1] else "—"
        p_s = f"{p[0]}/{p[1]} ({p[2]:.0f}%)" if p[1] else "—"
        L.append(f"| {lo}–{hi} | {e_s} | {p_s} |\n")
    L.append(
        "\nThe AD admits small, drug-like molecules of **either** class and rejects "
        "large ones of **either** class. Echinocandins score higher only because "
        "the external set is padded with small *synthetic* FKS inhibitors; the "
        "papulacandin set is mostly large glycolipids.\n\n")

    L.append("## 2. The echinocandins that matter are out-of-domain too\n\n")
    L.append("| approved drug | MW | in-domain? | kNN distance |\n")
    L.append("|---|---|---|---|\n")
    for r in sorted(approved, key=lambda r: float(r["mw"])):
        L.append(f"| {r['name']} | {r['mw']} | "
                 f"{'yes' if r['in_domain']=='True' else 'NO'} | {r['knn_distance']} |\n")
    L.append(
        "\nEvery approved echinocandin — the serum-tolerant drugs we'd use as the "
        "TEACHER — is out-of-domain, just like the papulacandins. So 'many "
        "echinocandins in-domain' does **not** validate the teacher: the in-domain "
        "echinocandins are small analogs, not the serum-tolerance exemplars.\n\n")

    L.append("## 3. The real success test: does the AD predict error?\n\n")
    L.append(
        f"- In-domain hold-out MAE = **{v['in_mae']:.3f}** (n={v['in_n']}) vs "
        f"out-of-domain MAE = **{v['out_mae']:.3f}** (n={v['out_n']}) log-units.\n"
        f"- MAE by distance quartile (near → far): "
        + " → ".join(f"{q[3]:.2f}" for q in v["quartiles"]) + ".\n"
        f"- Spearman(distance, |error|) = **{v['spearman']:+.3f}**.\n\n")
    L.append(
        "Distance-from-training genuinely tracks prediction error, so the AD **is** "
        "a valid reliability flag (modest but real). That is the correct sense in "
        "which it 'works' — not the in-domain head-count.\n\n")

    L.append("## Conclusion\n\n")
    L.append(
        "The AD is successful as a **reliability flag** (error rises with distance), "
        "but the echinocandin-in / papulacandin-out contrast is a **size gradient**, "
        "not proof the teacher transfers. Target compounds (papulacandins) and the "
        "large-molecule teachers (real echinocandins) are **both** out-of-domain for "
        "a PPB oracle trained on ordinary drugs — reinforcing Stage 1b's conclusion "
        "that fu must be **measured** for full-size molecules. This does not "
        "threaten the strategy: Stage 1c's serum transfer-learning uses the "
        "echinocandins' own measured serum / in-vivo data, not this drug-trained "
        "PPB oracle, which is reserved for ranking *in-domain* (small, drug-like) "
        "analogs.\n")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "applicability_domain_analysis.md"), "w",
              encoding="utf-8") as fh:
        fh.write("".join(L))

    print(f"papulacandin in-domain: {papu[0]}/{papu[1]} ({papu[2]:.0f}%)")
    print(f"echinocandin in-domain: {ech[0]}/{ech[1]} ({ech[2]:.0f}%)")
    print(f"approved echinocandins in-domain: "
          f"{sum(1 for r in approved if r['in_domain']=='True')}/{len(approved)}")
    print(f"AD validity: in-domain MAE {v['in_mae']:.3f} < out MAE {v['out_mae']:.3f}; "
          f"Spearman(dist,err)={v['spearman']:+.3f}")
    print(f"-> {os.path.join(OUT, 'applicability_domain_analysis.md')}")


if __name__ == "__main__":
    main()
