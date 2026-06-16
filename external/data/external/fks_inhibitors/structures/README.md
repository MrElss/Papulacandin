# External Structures

Store machine-readable external compound structures here.

Preferred order:

1. Import SDF or canonical SMILES from ChEMBL, PubChem, BindingDB, supplier, or paper SI.
2. Standardize with RDKit.
3. Only hand-fix structures that are missing, ambiguous, or clearly wrong.

Do not manually draw every external compound. Manual drawing is reserved for
important compounds whose database structures are unavailable or suspicious.

Suggested subfolders when the collection grows:

```text
sdf/
mol/
smiles/
qc_reports/
```

