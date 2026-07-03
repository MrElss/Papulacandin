#!/usr/bin/env python3
"""
make_deck_figures.py
====================
Generate a cohesive set of professional, presentation-grade figures for the
promotion deck. All share one academic style and palette. Real project data is
used where it exists (serum-gap table, Phase-13/14 outputs); schematic/framework
figures (pipeline, QM funnel) are drawn cleanly with matplotlib patches.

Outputs -> analysis/outputs/deck_figures/*.png  (white bg, 200 dpi)
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
FIG = os.path.join(OUT, "deck_figures")
os.makedirs(FIG, exist_ok=True)

INK = "#14233A"; BLUE = "#2F56A6"; TEAL = "#1F8A8A"; CORAL = "#D84A3C"
GREEN = "#2E8B57"; AMBER = "#D88A2B"; GREY = "#5B6472"; LIGHT = "#EEF2F8"

mpl.rcParams.update({
    "figure.dpi": 200, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 12, "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GREY, "axes.labelcolor": INK, "text.color": INK,
    "axes.titlecolor": INK, "xtick.color": GREY, "ytick.color": GREY,
    "axes.titlesize": 15, "axes.titleweight": "bold", "axes.labelsize": 13,
})


def _parse_mic(s):
    s = str(s).strip(); rel = ""
    for r in (">=", "<=", ">", "<", "="):
        if s.startswith(r):
            rel = r; s = s[len(r):]; break
    try:
        return float(s), rel
    except ValueError:
        return np.nan, rel


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, facecolor="white"); plt.close(fig)
    print("wrote", os.path.relpath(p, HERE))


# ---------------------------------------------------------------- 1. serum gap
def fig_serum_gap():
    df = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv")).drop_duplicates("compound_id")
    rows = []
    for _, r in df.iterrows():
        f, _ = _parse_mic(r["serumfree_mic_ugml"])
        sv, rel = _parse_mic(r["serum_mic_ugml"])
        if np.isfinite(f) and np.isfinite(sv):
            rows.append((f, sv, rel, r["serum_active"]))
    rows.sort(key=lambda x: (x[3] != "yes", x[0]))
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    for i, (f, sv, rel, act) in enumerate(rows):
        col = GREEN if act == "yes" else CORAL
        ax.plot([f, sv], [i, i], color=col, lw=2.2, alpha=0.55, zorder=1)
        ax.scatter([f], [i], s=44, color=BLUE, zorder=3, edgecolor="white", lw=0.6)
        marker = ">" if rel == ">" else "o"
        ax.scatter([sv], [i], s=60, color=col, marker=marker, zorder=3,
                   edgecolor="white", lw=0.6)
    ax.set_xscale("log")
    ax.axvline(100, ls="--", color=GREY, lw=1)
    ax.text(102, len(rows) - 1.5, "assay ceiling\n(activity lost)", color=GREY,
            fontsize=10, va="top")
    ax.set_xlabel("MIC — dose needed to stop the fungus  (µg/mL, log scale)")
    ax.set_yticks([])
    ax.set_ylabel("24 matched compounds")
    ax.set_title("The serum gap: blood serum abolishes activity for most analogs")
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=BLUE, markersize=9, label="serum-free MIC (potent)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=GREEN, markersize=9, label="serum MIC — tolerant (13)"),
           Line2D([0], [0], marker=">", color="w", markerfacecolor=CORAL, markersize=10, label="serum MIC — activity lost (11)")]
    ax.legend(handles=leg, loc="upper left", frameon=True, fontsize=10.5)
    ax.annotate("", xy=(11.5, 3), xytext=(0.85, 3),
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.6))
    ax.text(3.2, 3.6, "“serum shift” = the gap we must close", color=INK, fontsize=11,
            style="italic", ha="center")
    save(fig, "fig_serum_gap.png")


# ------------------------------------------------------------- 2. pipeline
def _rbox(ax, x, y, w, h, text, fc, tc="white", fs=11, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                fc=fc, ec="none", zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="bold" if bold else "normal", zorder=3)


def fig_pipeline():
    fig, ax = plt.subplots(figsize=(11.5, 4.6)); ax.axis("off")
    ax.set_xlim(0, 12); ax.set_ylim(0, 5)
    acts = [("ACT I", "Data &\nfirst clues", "Phases 1–4", BLUE),
            ("ACT II", "Go 3-D\n(quantum)", "Phases 5–9", TEAL),
            ("ACT III", "Learn from\napproved drugs", "Phases 10–11", AMBER),
            ("ACT IV", "AI design +\nstress-test", "Phases 12–14", GREEN)]
    w, gap = 2.5, 0.55; x = 0.3
    for i, (act, body, ph, col) in enumerate(acts):
        _rbox(ax, x, 2.6, w, 1.7, "", col)
        ax.text(x + w / 2, 3.95, act, ha="center", color="white", fontsize=13, fontweight="bold")
        ax.text(x + w / 2, 3.25, body, ha="center", color="white", fontsize=11.5)
        ax.text(x + w / 2, 2.35, ph, ha="center", color=col, fontsize=10.5, fontweight="bold")
        if i < 3:
            ax.add_patch(FancyArrowPatch((x + w + 0.04, 3.45), (x + w + gap - 0.04, 3.45),
                                         arrowstyle="-|>", mutation_scale=18, color=GREY, lw=2))
        x += w + gap
    # narrowing funnel underlay
    ax.add_patch(Polygon([(0.3, 1.9), (11.85, 1.9), (9.6, 0.7), (2.55, 0.7)],
                         closed=True, fc=LIGHT, ec="none", zorder=1))
    ax.text(6.07, 1.28, "broad & cheap  →  narrow & high-confidence", ha="center",
            color=INK, fontsize=12, fontstyle="italic")
    ax.text(6.07, 4.75, "A 14-stage AI + quantum-chemistry discovery pipeline",
            ha="center", color=INK, fontsize=15, fontweight="bold")
    save(fig, "fig_pipeline.png")


# ------------------------------------------------------------- 3. QM funnel
def fig_qm_funnel():
    fig, ax = plt.subplots(figsize=(7.6, 6.2)); ax.axis("off")
    ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    tiers = [("AI molecule generation", "novel structures", 9.0, 9.6, BLUE),
             ("2-D rules & drug-likeness", "cheap filter", 7.2, 7.6, TEAL),
             ("Cheap 3-D  (GFN-FF)", "conformer screen", 5.4, 5.9, AMBER),
             ("Accurate 3-D  (GFN2)", "re-rank; catch artifacts", 3.6, 4.4, CORAL),
             ("Quantum DFT", "finalists only", 1.8, 2.8, INK)]
    top = 9.6
    for i, (name, sub, yc, halfw, col) in enumerate(tiers):
        w_top = 9.2 * (1 - i * 0.17); w_bot = 9.2 * (1 - (i + 1) * 0.17)
        yt, yb = 9.7 - i * 1.75, 9.7 - (i + 1) * 1.75 + 0.25
        cx = 5.0
        ax.add_patch(Polygon([(cx - w_top / 2, yt), (cx + w_top / 2, yt),
                              (cx + w_bot / 2, yb), (cx - w_bot / 2, yb)],
                             closed=True, fc=col, ec="white", lw=1.5, zorder=2, alpha=0.92))
        ax.text(cx, (yt + yb) / 2 + 0.12, name, ha="center", va="center",
                color="white", fontsize=12.5, fontweight="bold", zorder=3)
        ax.text(cx, (yt + yb) / 2 - 0.32, sub, ha="center", va="center",
                color="white", fontsize=10, zorder=3, alpha=0.95)
    ax.annotate("", xy=(9.4, 1.3), xytext=(9.4, 9.5),
                arrowprops=dict(arrowstyle="-|>", color=GREY, lw=2))
    ax.text(9.75, 5.4, "increasing accuracy & cost", rotation=90, va="center",
            color=GREY, fontsize=11)
    ax.text(5, 0.35, "Each tier is more accurate and more expensive — so we\nspend"
            " the costly methods only on the few survivors.", ha="center",
            color=INK, fontsize=10.5, fontstyle="italic")
    ax.set_title("The staged quantum-chemistry funnel", y=0.99)
    save(fig, "fig_qm_funnel.png")


# ------------------------------------------------------------- 4. honesty ladder
def fig_honesty_ladder():
    labels = ["2-D\nchemistry", "Cheap 3-D\n(MMFF)", "Accurate 3-D\n(real CREST)",
              "…controlling\nfor potency"]
    vals = [0.32, 0.45, 0.31, 0.02]
    sig = [False, True, False, False]
    cols = [GREY, AMBER, BLUE, CORAL]
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    bars = ax.bar(labels, vals, color=cols, width=0.62, zorder=3, edgecolor="white")
    ax.axhline(0.40, ls="--", color=GREY, lw=1)
    ax.text(3.4, 0.415, "≈ significance\nfloor (n=24)", color=GREY, fontsize=9.5, va="bottom", ha="right")
    notes = ["weak\n(promising)", "STRONG\n(tempting!)", "shrinks\n(n.s.)", "≈ 0\n(artifact)"]
    for b, v, nt, c in zip(bars, vals, notes, cols):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"|ρ| = {v:.2f}",
                ha="center", fontsize=11.5, fontweight="bold", color=INK)
        ax.text(b.get_x() + b.get_width() / 2, v / 2, nt, ha="center", va="center",
                color="white", fontsize=10.5, fontweight="bold")
    ax.add_patch(FancyArrowPatch((1.15, 0.47), (3.0, 0.06), connectionstyle="arc3,rad=-0.25",
                                 arrowstyle="-|>", mutation_scale=20, color=CORAL, lw=2.2, zorder=5))
    ax.text(2.5, 0.36, "more rigor →\nsignal collapses", color=CORAL, fontsize=11,
            fontweight="bold", ha="center")
    ax.set_ylabel("apparent |correlation| with serum tolerance")
    ax.set_ylim(0, 0.52)
    ax.set_title("The “honesty ladder”: each rise in rigor shrank the easy signal")
    save(fig, "fig_honesty_ladder.png")


# ------------------------------------------------------------- 5. echinocandins
def fig_echinocandin():
    fig, ax = plt.subplots(figsize=(9.4, 5.8))
    ax.add_patch(Rectangle((0, 1), 0.5, 29, fc=GREEN, alpha=0.07, zorder=0))
    ax.add_patch(Rectangle((0.5, 30), 0.5, 1470, fc=CORAL, alpha=0.08, zorder=0))
    # drug: (x rigidity, low, high, color, label, text-xy in data coords, ha)
    drugs = [(0.18, 2, 11, GREEN, "Caspofungin\n(branched, flexible)", (0.31, 14), "left"),
             (0.55, 16, 32, INK, "Native papulacandin\ntail (rigid polyene)", (0.23, 78), "left"),
             (0.78, 16, 512, CORAL, "Anidulafungin\n(rigid, aromatic)", (0.49, 200), "left"),
             (0.90, 64, 1024, CORAL, "Micafungin\n(rigid, aromatic)", (0.60, 600), "left")]
    for x, lo, hi, col, lab, (tx, ty), ha in drugs:
        gm = (lo * hi) ** 0.5
        ax.errorbar([x], [gm], yerr=[[gm - lo], [hi - gm]], fmt="o", ms=13, color=col,
                    ecolor=col, elinewidth=2, capsize=5, zorder=4, mec="white", mew=1)
        ax.annotate(lab, xy=(x, gm), xytext=(tx, ty), ha=ha, va="center",
                    fontsize=10.5, fontweight="bold", color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=1, alpha=0.55))
    ax.set_yscale("log")
    ax.set_xlim(0, 1); ax.set_ylim(1, 1500)
    ax.set_xticks([0.15, 0.85]); ax.set_xticklabels(["flexible / saturated", "rigid / aromatic"], fontsize=12)
    ax.set_xlabel("tail shape  →")
    ax.set_ylabel("serum shift  (activity lost, log)")
    ax.set_title("Approved-drug evidence: rigid tails lose serum activity")
    ax.text(0.06, 1.35, "GOOD zone", color=GREEN, fontsize=11.5, fontweight="bold", ha="left")
    ax.text(0.56, 1050, "BAD zone", color=CORAL, fontsize=11.5, fontweight="bold", ha="left")
    save(fig, "fig_echinocandin.png")


# ------------------------------------------------------------- 6. round-1 result
def fig_round1():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.2, 5.0), gridspec_kw={"width_ratios": [1, 1]})
    # left: attrition — decreasing centered blocks, big count inside, label below
    axL.axis("off"); axL.set_xlim(0, 10); axL.set_ylim(0, 10)
    blocks = [(12, "12", "AI-designed tails", BLUE, 8.1),
              (2, "2", "pass cheap 3-D screen", AMBER, 5.0),
              (0, "0", "survive accurate GFN2", CORAL, 1.9)]
    maxw = 8.6
    for cnt, num, lab, col, yc in blocks:
        w = maxw * (cnt / 12.0) if cnt > 0 else 1.3
        axL.add_patch(FancyBboxPatch((5 - w / 2, yc - 0.6), w, 1.2,
                                     boxstyle="round,pad=0.02,rounding_size=0.12",
                                     fc=col, ec="none"))
        axL.text(5, yc, num, ha="center", va="center", color="white",
                 fontsize=21, fontweight="bold")
        axL.text(5, yc - 1.15, lab, ha="center", va="center", color=INK,
                 fontsize=11.5, fontweight="bold")
    axL.annotate("", xy=(5, 6.0), xytext=(5, 6.9), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=2))
    axL.annotate("", xy=(5, 2.9), xytext=(5, 3.8), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=2))
    axL.set_title("12 → 2 → 0", fontsize=14)
    # right: the artifact — sulfonate GFN-FF vs GFN2 hydrophobic fraction
    g = pd.read_csv(os.path.join(OUT, "phase13_gfn2_ranking.csv")).set_index("finalist")
    r = g.loc["t01_C8_omega_sulfonate"]
    axR.bar(["cheap 3-D\n(GFN-FF)", "accurate\n(GFN2)"],
            [r["gfnff_hydrophobic_fraction"], r["gfn2_hydrophobic_fraction"]],
            color=[AMBER, CORAL], width=0.55, zorder=3, edgecolor="white")
    axR.axhline(0.58, ls="--", color=INK, lw=1.4)
    axR.text(1.45, 0.585, "native baseline", color=INK, fontsize=10, ha="right", va="bottom")
    for i, v in enumerate([r["gfnff_hydrophobic_fraction"], r["gfn2_hydrophobic_fraction"]]):
        axR.text(i, v + 0.006, f"{v:.2f}", ha="center", fontweight="bold", fontsize=12)
    axR.annotate("looks better\nthan native", (0, 0.534), xytext=(-0.15, 0.50),
                 fontsize=10, color=GREEN, ha="center")
    axR.annotate("actually WORSE\n(folds & buries)", (1, 0.636), xytext=(1.0, 0.66),
                 fontsize=10, color=CORAL, ha="center", fontweight="bold")
    axR.set_ylim(0.45, 0.70); axR.set_ylabel("exposed greasy surface (fraction)")
    axR.set_title("The 'winner' was a force-field artifact")
    fig.suptitle("Round 1: the quantum funnel caught a false positive before synthesis",
                 fontsize=15, fontweight="bold", y=1.02)
    save(fig, "fig_round1.png")


# ------------------------------------------------------------- 7. tail ladder
def fig_tail_ladder():
    df = pd.read_csv(os.path.join(OUT, "phase14_tail_series.csv"))
    order = ["native_C16_polyene", "C16_1_palmitoleoyl", "C16_0_palmitoyl",
             "branched_sat_caspofungin_like"]
    disp = {"native_C16_polyene": ("Native polyene", "the problem", CORAL, ""),
            "C16_1_palmitoleoyl": ("Palmitoleoyl (1 kink)", "control", TEAL, ""),
            "C16_0_palmitoyl": ("Palmitoyl (saturated)", "★ make first", GREEN, "★"),
            "branched_sat_caspofungin_like": ("Branched saturated", "★ make first", GREEN, "★")}
    df = df.set_index("name").loc[order]
    fig, ax = plt.subplots(figsize=(9.8, 5.2))
    y = np.arange(len(order))[::-1]
    cvals = df["tail_CC_double_bonds"].values
    cols = [disp[n][2] for n in order]
    ax.barh(y, cvals, color=cols, height=0.55, zorder=3, edgecolor="white")
    for yi, n in zip(y, order):
        name, tag, col, star = disp[n]
        ax.text(-0.12, yi, name, ha="right", va="center", fontsize=12, fontweight="bold", color=INK)
        ax.text(df.loc[n, "tail_CC_double_bonds"] + 0.12, yi, tag, va="center",
                fontsize=11, color=col, fontweight="bold")
    ax.set_yticks([])
    ax.set_xlabel("tail rigidity  =  number of C=C double bonds  (fewer = more flexible)")
    ax.set_xlim(0, 5)
    ax.set_title("Round-2 design: de-rigidify the tail (same length, keep potency)")
    ax.add_patch(FancyArrowPatch((3.9, 3), (3.9, 0.2), arrowstyle="-|>",
                                 mutation_scale=20, color=INK, lw=2))
    ax.text(4.15, 1.6, "de-rigidify", rotation=90, va="center", color=INK,
            fontsize=11, fontweight="bold")
    ax.text(2.35, 0.55, "All four tails are C16 (same length &\ngreasiness) — only the SHAPE changes.",
            fontsize=10.5, fontstyle="italic", color=GREY, ha="center", va="center")
    save(fig, "fig_tail_ladder.png")


# ------------------------------------------------------------- 8. potency confound
def fig_confound():
    df = pd.read_csv(os.path.join(OUT, "serum_gap_pairs.csv")).drop_duplicates("compound_id")
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.4, 5.2), gridspec_kw={"width_ratios": [1, 1]})
    # A: serum-free vs serum MIC
    for _, r in df.iterrows():
        f, _ = _parse_mic(r["serumfree_mic_ugml"]); sv, rel = _parse_mic(r["serum_mic_ugml"])
        if np.isfinite(f) and np.isfinite(sv):
            col = GREEN if r["serum_active"] == "yes" else CORAL
            axA.scatter([f], [sv], s=70, color=col, alpha=0.8, edgecolor="white", lw=0.6,
                        marker=("^" if rel == ">" else "o"), zorder=3)
    lims = [0.3, 200]
    axA.plot(lims, lims, ls="--", color=GREY, lw=1, zorder=1)
    axA.set_xscale("log"); axA.set_yscale("log"); axA.set_xlim(*lims); axA.set_ylim(*lims)
    axA.set_xlabel("serum-FREE MIC  (potency, log)")
    axA.set_ylabel("serum MIC  (log)")
    axA.set_title("Serum outcome ≈ intrinsic potency", fontsize=13.5)
    axA.text(0.5, 130, "ρ = +0.79\n(potency dominates)", color=INK, fontsize=12,
             fontweight="bold", va="top")
    # B: partial correlation after controlling potency
    cf = pd.read_csv(os.path.join(OUT, "phase8_confound_analysis.csv"))
    names = {"hydrophobic_sasa_mean": "greasy surface", "hydrophobic_fraction_mean": "greasy fraction",
             "hydrophobic_sasa_std": "rigidity", "polar_sasa_mean": "POLAR surface", "rg_mean": "size (Rg)"}
    cf = cf[cf["descriptor"].isin(names)].copy()
    cf["label"] = cf["descriptor"].map(names)
    cf = cf.sort_values("partial_rho_vs_serumMIC")
    cols = [GREEN if l == "POLAR surface" else GREY for l in cf["label"]]
    axB.barh(cf["label"], cf["partial_rho_vs_serumMIC"], color=cols, zorder=3, edgecolor="white")
    axB.axvline(0, color=INK, lw=1)
    axB.set_xlim(-0.5, 0.5)
    axB.set_xlabel("correlation after removing potency  (partial ρ)")
    axB.set_title("After removing potency: signal ≈ 0", fontsize=13.5)
    labs = list(cf["label"])
    axB.annotate("greasy surface — the\nPhase-7 lead → now ≈ 0",
                 xy=(0.02, labs.index("greasy surface")), xytext=(0.14, labs.index("greasy surface") - 0.9),
                 color=GREY, fontsize=9.5, arrowprops=dict(arrowstyle="->", color=GREY, lw=1))
    axB.annotate("← polar surface:\nthe one weak survivor", xy=(-0.31, labs.index("POLAR surface")),
                 xytext=(0.03, labs.index("POLAR surface") + 0.15), color=GREEN, fontsize=10,
                 fontweight="bold", va="center")
    fig.suptitle("The decisive check: the exciting 3-D signal was mostly just potency",
                 fontsize=15, fontweight="bold", y=1.02)
    save(fig, "fig_confound.png")


# ------------------------------------------------------------- 9. all-12 tail ranking
def fig_tail_ranking():
    df = pd.read_csv(os.path.join(OUT, "phase13_qm_ranking.csv")).sort_values(
        "hydrophobic_fraction_mean", ascending=False)
    native = 0.58
    fig, ax = plt.subplots(figsize=(10.2, 5.8))
    y = np.arange(len(df))
    cols = [GREEN if v < native else GREY for v in df["hydrophobic_fraction_mean"]]
    ax.barh(y, df["hydrophobic_fraction_mean"], color=cols, zorder=3, edgecolor="white")
    ax.set_yticks(y); ax.set_yticklabels(df["tail_name"], fontsize=10.5)
    ax.axvline(native, ls="--", color=CORAL, lw=1.8)
    ax.text(native + 0.004, len(df) - 0.5, "native\nbaseline", color=CORAL, fontsize=10.5, va="top")
    ax.set_xlim(0.45, 0.72)
    ax.set_xlabel("exposed GREASY surface (fraction) — lower is better")
    ax.set_title("Round-1 screen: only 2 of 12 tails beat native (cheap 3-D tier)")
    for yi, v in zip(y, df["hydrophobic_fraction_mean"]):
        if v < native:
            ax.text(v - 0.004, yi, "✓", color="white", fontsize=13, va="center", ha="right", fontweight="bold")
    save(fig, "fig_tail_ranking.png")


# ------------------------------------------------------------- 10. electronics (Phase 9)
def fig_electronics():
    df = pd.read_csv(os.path.join(OUT, "phase9_electronic_stats.csv"))
    names = {"alpha_au": "polarizability α", "logP_xtb": "QM logP", "dipole_D": "dipole",
             "homo_lumo_gap_eV": "HOMO–LUMO gap", "gsolv_water_kcal": "solvation"}
    df = df[df["descriptor"].isin(names)].copy(); df["label"] = df["descriptor"].map(names)
    df = df.sort_values("rho_serumMIC")
    fig, ax = plt.subplots(figsize=(10.0, 5.4))
    yy = np.arange(len(df)); h = 0.38
    ax.barh(yy + h / 2, df["rho_serumMIC"], height=h, color=BLUE, label="vs serum MIC (raw)", zorder=3, edgecolor="white")
    ax.barh(yy - h / 2, df["rho_shift"], height=h, color=TEAL, label="vs serum SHIFT (potency-free)", zorder=3, edgecolor="white")
    ax.set_yticks(yy); ax.set_yticklabels(df["label"], fontsize=11)
    ax.axvline(0, color=INK, lw=1); ax.set_xlim(-0.7, 0.5)
    ax.set_xlabel("Spearman correlation")
    ax.set_title("Electronics: α is a size/potency proxy; QM logP corroborates polar-surface")
    ax.legend(loc="lower left", fontsize=10, frameon=True)
    # annotate the two key points
    a = df[df["label"] == "polarizability α"].iloc[0]
    ax.annotate("big raw signal (−0.54)\nbut a SIZE proxy (shift ≈ 0)", (a["rho_serumMIC"], list(df["label"]).index("polarizability α") + h / 2),
                xytext=(-0.66, list(df["label"]).index("polarizability α") - 1.2), fontsize=9.5, color=GREY,
                arrowprops=dict(arrowstyle="->", color=GREY, lw=1))
    save(fig, "fig_electronics.png")


# ------------------------------------------------------------- 11. evidence base
def fig_dataset():
    fig, ax = plt.subplots(figsize=(11.2, 4.6)); ax.axis("off")
    ax.set_xlim(0, 12); ax.set_ylim(0, 5)
    def stat(x, w, num, lab, col):
        ax.add_patch(FancyBboxPatch((x, 1.5), w, 2.4, boxstyle="round,pad=0.02,rounding_size=0.08",
                                    fc=col, ec="none"))
        ax.text(x + w / 2, 3.1, num, ha="center", va="center", color="white", fontsize=30, fontweight="bold")
        ax.text(x + w / 2, 2.05, lab, ha="center", va="center", color="white", fontsize=12)
    stat(0.3, 2.5, "138", "curated\ncompounds", BLUE)
    stat(3.05, 2.5, "1,042", "activity\nrecords", TEAL)
    stat(5.8, 2.5, "24", "matched serum\npairs (the target)", INK)
    # split of the 24
    ax.add_patch(FancyBboxPatch((8.6, 2.75), 3.1, 1.15, boxstyle="round,pad=0.02,rounding_size=0.08", fc=GREEN, ec="none"))
    ax.text(10.15, 3.32, "13 serum-tolerant", ha="center", va="center", color="white", fontsize=13, fontweight="bold")
    ax.add_patch(FancyBboxPatch((8.6, 1.5), 3.1, 1.15, boxstyle="round,pad=0.02,rounding_size=0.08", fc=CORAL, ec="none"))
    ax.text(10.15, 2.07, "11 serum-killed", ha="center", va="center", color="white", fontsize=13, fontweight="bold")
    ax.add_patch(FancyArrowPatch((8.35, 2.7), (8.58, 2.7), arrowstyle="-|>", mutation_scale=15, color=GREY, lw=2))
    ax.text(6.07, 0.75, "+ an independent echinocandin drug corpus (3 approved drugs, 279 serum measurements) "
            "imported for cross-validation", ha="center", fontsize=11, color=INK, fontstyle="italic")
    ax.text(6.07, 4.55, "The evidence base — curated from decades of scattered literature",
            ha="center", fontsize=15, fontweight="bold", color=INK)
    save(fig, "fig_dataset.png")


if __name__ == "__main__":
    fig_serum_gap()
    fig_pipeline()
    fig_qm_funnel()
    fig_honesty_ladder()
    fig_echinocandin()
    fig_round1()
    fig_tail_ladder()
    fig_confound()
    fig_tail_ranking()
    fig_electronics()
    fig_dataset()
    print("\nAll deck figures written to", os.path.relpath(FIG, HERE))
