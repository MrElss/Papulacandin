"""
Construct the rational-design panel by editing the validated papulacandin B scaffold.
Two acyl positions: pos1 = long-tail ester on the glucoside core; pos2 = short-tail
ester on the galactose C6. We vary (a) the pos1 acyl (echinocandin side chains +
saturation series), (b) ester->amide linkage, (c) single vs dual-position edits.
No RDKit here, so we self-check paren/bracket balance and ring-digit parity, and
verify the 'native' rebuild equals the known-good papulacandin B SMILES.
"""

# --- validated reference (from earlier ADME/Inductive Bio runs) ---
PAPB_KNOWN = "CCC(C)CC/C=C/C=C(C)/C(O)C/C=C/C=C/C(=O)O[C@@H]1[C@@H](O)[C@@]2(OCc3cc(O)cc(O)c23)O[C@H](CO)[C@H]1O[C@@H]4O[C@H](COC(=O)/C=C/C=C/C=C/C(O)CC)[C@H](O)[C@H](O)[C@H]4O"

# glycoside core body; {p2link} = 'O'(ester) or 'N'(amide); {p2acyl} carbonyl-first
def glyco(p2link, p2acyl):
    return (f"[C@@H]1[C@@H](O)[C@@]2(OCc3cc(O)cc(O)c23)O[C@H](CO)[C@H]1O"
            f"[C@@H]4O[C@H](C{p2link}{p2acyl})[C@H](O)[C@H](O)[C@H]4O")

# pos1 aliphatic acyls (end in C(=O)); pos1 aromatic acyls use {X}=link+glyco, rings 5-7
NATIVE_LONG = "CCC(C)CC/C=C/C=C(C)/C(O)C/C=C/C=C/C(=O)"
CASPO       = "CC[C@H](C)C[C@H](C)CCCCCCCCC(=O)"          # branched saturated (flexible)
PALMITOYL   = "CCCCCCCCCCCCCCCC(=O)"                        # C16 straight (flexible)
ONEKINK     = "CCCCCC/C=C\\CCCCCCCC(=O)"                    # C16:1 one cis kink
ANIDULA_T   = "CCCCCOc5ccc(-c6ccc(-c7ccc(C(=O){X})cc7)cc6)cc5"   # pentyloxy-terphenyl
MICA_T      = "CCCCCOc5ccc(-c6cc(-c7ccc(C(=O){X})cc7)no6)cc5"    # pentyloxy-phenyl-isoxazole-phenyl

# pos2 acyls (carbonyl-first)
NATIVE_SHORT = "C(=O)/C=C/C=C/C=C/C(O)CC"
CASPO_P2     = "C(=O)CCCCCCCCCC(C)CC(C)CC"                  # branched saturated (achiral form)

def ali(p1acyl, p1link, p2link, p2acyl):    # aliphatic pos1
    return p1acyl + p1link + glyco(p2link, p2acyl)
def aro(p1tmpl, p1link, p2link, p2acyl):    # aromatic pos1
    return p1tmpl.replace("{X}", p1link + glyco(p2link, p2acyl))

IBREX = "CC(C)[C@@H](C)[C@@]1(C)CC[C@]2(C)[C@H]3CC[C@@H]4[C@@]5(COC[C@@]4(C)[C@@H](OC[C@](C)(N)C(C)(C)C)[C@H](n4ncnc4-c4ccncc4)C5)C3=CC[C@@]2(C)[C@@H]1C(=O)O"
FLUC  = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"

designs = [
    ("papB_native",         ali(NATIVE_LONG, "O", "O", NATIVE_SHORT), "reference: native rigid polyene tail, di-ester"),
    ("A1_caspo_ester",      ali(CASPO,      "O", "O", NATIVE_SHORT),  "pos1 caspofungin branched-sat tail, ester"),
    ("A2_palmitoyl_ester",  ali(PALMITOYL,  "O", "O", NATIVE_SHORT),  "pos1 palmitoyl straight-sat tail, ester"),
    ("A3_onekink_ester",    ali(ONEKINK,    "O", "O", NATIVE_SHORT),  "pos1 C16:1 one-kink control, ester"),
    ("G_anidula_ester",     aro(ANIDULA_T,  "O", "O", NATIVE_SHORT),  "pos1 anidulafungin terphenyl (rigid aromatic), ester"),
    ("G_mica_ester",        aro(MICA_T,     "O", "O", NATIVE_SHORT),  "pos1 micafungin isoxazole-diaryl (rigid aromatic), ester"),
    ("M_palmitoyl_amide",   ali(PALMITOYL,  "N", "O", NATIVE_SHORT),  "pos1 palmitoyl, ESTER->AMIDE (esterase-hardened)"),
    ("M_caspo_amide",       ali(CASPO,      "N", "O", NATIVE_SHORT),  "pos1 caspo branched, ESTER->AMIDE"),
    ("M_anidula_amide",     aro(ANIDULA_T,  "N", "O", NATIVE_SHORT),  "pos1 anidulafungin terphenyl, ESTER->AMIDE"),
    ("D_caspo_both_ester",  ali(CASPO,      "O", "O", CASPO_P2),      "caspo branched at BOTH positions, di-ester"),
    ("D_caspo_both_amide",  ali(CASPO,      "N", "N", CASPO_P2),      "caspo branched at BOTH positions, di-AMIDE (fully hardened)"),
    ("ref_ibrexafungerp",   IBREX, "reference oral drug"),
    ("ref_fluconazole",     FLUC,  "reference small oral drug"),
]

def sane(smi):
    issues = []
    if smi.count("(") != smi.count(")"): issues.append("paren")
    if smi.count("[") != smi.count("]"): issues.append("bracket")
    for d in "123456789":
        # ring-closure digits: count only those not inside brackets-as-charges (none here)
        if smi.count(d) % 2 != 0: issues.append(f"ring{d}")
    return issues or "OK"

print(f"native rebuild matches known papulacandin B: {designs[0][1] == PAPB_KNOWN}\n")
import csv
out = "/home/claude/papulacandin_work/design_smiles.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["id","smiles","rationale","sanity","nC","nO","nN"])
    for name, smi, note in designs:
        nC = smi.count("C") + smi.count("c")
        nO = smi.count("O") + smi.count("o")
        nN = smi.count("N") + smi.count("n")
        print(f"{name:22} sane={str(sane(smi)):22} C{nC} O{nO} N{nN}")
        w.writerow([name, smi, note, sane(smi), nC, nO, nN])
print("\nwrote", out)
