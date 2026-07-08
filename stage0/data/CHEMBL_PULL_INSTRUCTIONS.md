# How to pull the ChEMBL echinocandin data (for you to download)

The ChEMBL connector isn't reachable from my side, but you can reach ChEMBL in a
browser. **You download; I parse.** You do NOT need to classify anything — just
grab the raw bioactivity tables and drop them in. The parser
(`stage0/convert_chembl_export.py`) does the serum-detection and active/inactive
calls.

## What we need

Bioactivity records for the echinocandins / FKS inhibitors — I'll keep only the
two useful kinds automatically:
- **serum-context MICs** (assays whose description mentions serum / plasma /
  albumin / protein binding), and
- **PPB / fraction-unbound (Fu)** records.

**Compounds to grab (each as its own file is fine):**
`caspofungin`, `anidulafungin`, `micafungin`, `rezafungin`, `ibrexafungerp`
(and any other FKS/glucan-synthase inhibitor you like).

## Easiest way — ChEMBL website (GUI)

1. Go to <https://www.ebi.ac.uk/chembl/>.
2. Search a drug name (e.g. *caspofungin*) and open its **compound report card**.
3. Scroll to the **Activities** / **Bioactivities** table.
4. Click the table's **download** button → choose **CSV** (or TSV). Save the whole
   table — no need to filter; I'll filter for serum/PPB.
5. Repeat for each compound.

## Precise way — API URL in your browser (optional)

If you'd rather use the API, first get each compound's ChEMBL ID (shown on its
report card, e.g. `CHEMBL…`), then open this URL (it downloads a CSV):

```
https://www.ebi.ac.uk/chembl/api/data/activity.csv?molecule_chembl_id=<CHEMBL_ID>&limit=1000
```

Repeat per compound. (One combined file or several separate files — both work.)

## What to name them and where to put them

- Save each file with a name that **starts with `chembl_raw`**, e.g.
  `chembl_raw_caspofungin.csv`, `chembl_raw_anidulafungin.csv`.
- Put them in **`stage0/data/`** in the repo and push/upload them to the
  `claude/repository-review-ta7yof` branch.

The files can be the website's `;`-delimited export **or** the API's
`,`-delimited CSV — the parser auto-detects both, and tolerates ChEMBL's
different header spellings.

## Then tell me — I'll run

```
python stage0/convert_chembl_export.py          # raw -> chembl_echinocandin_serum.csv (+ chembl_free_fraction.csv)
python stage0/build_binary_serum_activity.py    # merges into the binary label set
```

The serum records drop straight into the binary serum-active/inactive label set
(expanding the echinocandin side, which currently has only 6 direct compounds),
and the PPB/Fu rows extend the Stage-1b free-fraction seed. No code changes
needed on your end.

## Notes

- Don't worry about units or "which assays count" — the parser keeps ug/mL MICs,
  applies the >100 µg/mL = inactive rule, and skips serum-*free* assays
  automatically.
- If a download is huge, that's fine; the parser keeps only the relevant rows.
