"""
Assemble the first calibrated property panel for the papulacandin tail-analog set.
Data are from Boltz ADME-v1 (permeability, lipophilicity, solubility-confidence)
and Inductive Bio (LogD @ pH7.4, most-acidic pKa, most-basic pKa).
Scales for ADME-v1 permeability/lipophilicity are model-internal (unitless);
LogD is in real log units. Values flagged OOD = out-of-domain (extrapolation).
"""
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# id, role, ADME perm, ADME lipo, ADME sol-conf, LogD, LogD_flag, apKa, bpKa, note
rows = [
    ("papB_native",      "reference (native)", 0.0134, 0.367, "high",  1.39, "OOD/low-conf", 8.21, 2.0, "C18 tetraene hydroxy tail; 2 esters"),
    ("papB_A1_caspo",    "analog (branched-sat)", 0.304, 0.854, "med", 2.81, "in-domain",    8.19, 2.0, "caspofungin's own branched-saturated acyl grafted on"),
    ("papB_A2_palmitoyl","analog (flexible-sat)", 0.287, 1.113, "med", 2.98, "in-domain",    8.14, 2.0, "straight C16 palmitoyl; fully flexible"),
    ("papB_A3_one_kink", "analog (mono-unsat)",   0.245, 0.723, "med", 2.76, "in-domain",    8.11, 2.0, "C16:1 single cis kink (rigidity control)"),
    ("ibrexafungerp",    "reference (oral drug)", 0.093, 2.025, "med", 2.96, "OOD/low-conf", 3.90, 7.84, "approved oral GS inhibitor; low-affinity/reversible PPB"),
    ("caspofungin",      "reference (IV drug)",   None,  None,  None,  0.47, "in-domain",    7.60, 8.45, "ADME failed: too ionizable (252 microstates > 128 cap)"),
    ("fluconazole",      "reference (small oral)",1.605, 0.650, "high",0.73, "in-domain",    11.9, 3.09, "small well-absorbed oral anchor"),
]

# --- CSV ---
csv_path = "/home/claude/papulacandin_work/papulacandin_property_panel.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id","role","ADME_permeability(model units)","ADME_lipophilicity(model units)",
                "ADME_solubility_confidence","LogD_pH7.4","LogD_flag","acidic_pKa","basic_pKa","note"])
    for r in rows:
        w.writerow(r)
print("wrote", csv_path)

# --- Figure: LogD vs passive permeability ---
plot = [r for r in rows if r[2] is not None]  # exclude caspofungin (no ADME)
colors = {
    "papB_native": "#C0392B",
    "papB_A1_caspo": "#1F6FB2",
    "papB_A2_palmitoyl": "#2E86C1",
    "papB_A3_one_kink": "#5DADE2",
    "ibrexafungerp": "#1E8449",
    "fluconazole": "#7F8C8D",
}
labels = {
    "papB_native": "papulacandin B (native)",
    "papB_A1_caspo": "A1 branched-sat",
    "papB_A2_palmitoyl": "A2 palmitoyl",
    "papB_A3_one_kink": "A3 one-kink",
    "ibrexafungerp": "ibrexafungerp",
    "fluconazole": "fluconazole",
}

fig, ax = plt.subplots(figsize=(9.5, 6.6))
for r in plot:
    cid, logd, perm = r[0], r[5], r[2]
    ax.scatter(logd, perm, s=170, c=colors[cid], edgecolor="black", linewidth=0.8,
               zorder=3, marker=("*" if cid=="ibrexafungerp" else "o"))
    dy = 1.14 if cid != "papB_native" else 0.80
    ax.annotate(labels[cid], (logd, perm), xytext=(6, 8), textcoords="offset points",
                fontsize=10, fontweight=("bold" if cid.startswith("papB") else "normal"))

ax.set_yscale("log")
ax.set_xlabel("Predicted LogD @ pH 7.4  (Inductive Bio; real log units)", fontsize=11)
ax.set_ylabel("Predicted passive permeability  (Boltz ADME-v1; model units, log scale)", fontsize=11)
ax.set_title("First property panel: papulacandin tail analogs vs reference drugs",
             fontsize=13, fontweight="bold")
# Shade the ibrexafungerp-like lipophilicity band
ax.axvspan(2.6, 3.1, color="#2ECC71", alpha=0.08, zorder=0)
ax.text(2.85, ax.get_ylim()[0]*1.4, "ibrexafungerp\nLogD band", ha="center",
        fontsize=9, color="#1E8449")
ax.grid(True, which="both", axis="both", alpha=0.25)
ax.margins(x=0.15)

note = ("Native papulacandin: lowest permeability (relies on tail for membrane access).\n"
        "Tail analogs shift up ~18-23x in permeability and into ibrexafungerp's LogD band.\n"
        "Scales are model-internal & the papulacandin scaffold is out-of-domain -> directional only.\n"
        "Neither model predicts fraction-unbound, metabolic stability, or albumin/serum binding.")
ax.text(0.015, 0.02, note, transform=ax.transAxes, fontsize=8.3, va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", fc="#FBFCFC", ec="#BDC3C7"))

plt.tight_layout()
png = "/home/claude/papulacandin_work/papulacandin_property_panel.png"
plt.savefig(png, dpi=160)
print("wrote", png)
