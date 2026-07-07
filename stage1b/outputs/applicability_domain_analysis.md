# Stage 1b analysis — is the applicability domain successful?

Overall in-domain: **papulacandin 2/137 (1%)**, **echinocandin 91/321 (28%)**. The question: does the echinocandin-vs-papulacandin gap show the AD works?

## 1. The gap is molecular SIZE, not chemotype

| MW band | echinocandin in-domain | papulacandin in-domain |
|---|---|---|
| 0–500 | 30/41 (73%) | 2/9 (22%) |
| 500–700 | 61/121 (50%) | 0/17 (0%) |
| 700–900 | 0/13 (0%) | 0/33 (0%) |
| 900–3000 | 0/146 (0%) | 0/78 (0%) |

The AD admits small, drug-like molecules of **either** class and rejects large ones of **either** class. Echinocandins score higher only because the external set is padded with small *synthetic* FKS inhibitors; the papulacandin set is mostly large glycolipids.

## 2. The echinocandins that matter are out-of-domain too

| approved drug | MW | in-domain? | kNN distance |
|---|---|---|---|
| IBREXAFUNGERP | 730.1 | NO | 6.956 |
| CILOFUNGIN | 1030.1 | NO | 19.824 |
| CASPOFUNGIN | 1093.3 | NO | 23.556 |
| ANIDULAFUNGIN | 1140.3 | NO | 20.237 |
| MICAFUNGIN | 1270.3 | NO | 27.368 |

Every approved echinocandin — the serum-tolerant drugs we'd use as the TEACHER — is out-of-domain, just like the papulacandins. So 'many echinocandins in-domain' does **not** validate the teacher: the in-domain echinocandins are small analogs, not the serum-tolerance exemplars.

## 3. The real success test: does the AD predict error?

- In-domain hold-out MAE = **0.467** (n=184) vs out-of-domain MAE = **0.796** (n=10) log-units.
- MAE by distance quartile (near → far): 0.32 → 0.43 → 0.49 → 0.69.
- Spearman(distance, |error|) = **+0.332**.

Distance-from-training genuinely tracks prediction error, so the AD **is** a valid reliability flag (modest but real). That is the correct sense in which it 'works' — not the in-domain head-count.

## Conclusion

The AD is successful as a **reliability flag** (error rises with distance), but the echinocandin-in / papulacandin-out contrast is a **size gradient**, not proof the teacher transfers. Target compounds (papulacandins) and the large-molecule teachers (real echinocandins) are **both** out-of-domain for a PPB oracle trained on ordinary drugs — reinforcing Stage 1b's conclusion that fu must be **measured** for full-size molecules. This does not threaten the strategy: Stage 1c's serum transfer-learning uses the echinocandins' own measured serum / in-vivo data, not this drug-trained PPB oracle, which is reserved for ranking *in-domain* (small, drug-like) analogs.
