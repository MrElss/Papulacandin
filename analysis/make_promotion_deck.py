#!/usr/bin/env python3
"""
make_promotion_deck.py
======================
Build a professional, figure-rich promotion-report deck (16:9) for the
Papulacandin serum-gap program (Phases 1-14), for a general-faculty audience.

Embeds the cohesive figure set from make_deck_figures.py (run that first) plus a
real result plot; framework accents (chips, cards, roadmap loop, contributions
grid) are drawn with shapes.

Output: analysis/outputs/Papulacandin_promotion_report.pptx
"""
import os
import struct

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")
FIG = os.path.join(OUT, "deck_figures")

INK = RGBColor(0x14, 0x23, 0x3A); BLUE = RGBColor(0x2F, 0x56, 0xA6)
TEAL = RGBColor(0x1F, 0x8A, 0x8A); CORAL = RGBColor(0xD8, 0x4A, 0x3C)
GREEN = RGBColor(0x2E, 0x8B, 0x57); AMBER = RGBColor(0xD8, 0x8A, 0x2B)
GREY = RGBColor(0x5B, 0x64, 0x72); LIGHT = RGBColor(0xF2, 0xF5, 0xF9)
CARD2 = RGBColor(0xE8, 0xEE, 0xF6); LINE = RGBColor(0xD3, 0xDA, 0xE4)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)
SW, SH = prs.slide_width, prs.slide_height
BLANK = prs.slide_layouts[6]
ML = Inches(0.75); CW = SW - 2 * ML
_page = 0


def _solid(sh, c): sh.fill.solid(); sh.fill.fore_color.rgb = c; sh.line.fill.background()
def _bg(s, c=WHITE): s.background.fill.solid(); s.background.fill.fore_color.rgb = c


def rect(s, l, t, w, h, c, shape=MSO_SHAPE.RECTANGLE, line=None, lw=1.0):
    sh = s.shapes.add_shape(shape, int(l), int(t), int(w), int(h)); _solid(sh, c)
    if line is not None:
        sh.line.color.rgb = line; sh.line.width = Pt(lw)
        sh.line.fill.solid(); sh.line.fill.fore_color.rgb = line
    sh.shadow.inherit = False
    return sh


def txt(s, l, t, w, h, runs, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP):
    tb = s.shapes.add_textbox(int(l), int(t), int(w), int(h)); tf = tb.text_frame
    tf.word_wrap = True; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Pt(2); tf.margin_top = tf.margin_bottom = Pt(1)
    for i, para in enumerate(runs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align; p.space_after = Pt(4)
        for (text, size, bold, color, *rest) in para:
            r = p.add_run(); r.text = text; r.font.size = Pt(size)
            r.font.bold = bold; r.font.color.rgb = color; r.font.name = "Calibri"
            if rest and rest[0]:
                r.font.italic = True
    return tb


def chevron(s, l, t, w=Inches(0.34), h=Inches(0.34), c=GREY):
    rect(s, l, t, w, h, c, shape=MSO_SHAPE.CHEVRON)


def card(s, l, t, w, h, fill=LIGHT, line=LINE):
    return rect(s, l, t, w, h, fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE, line=line)


def _png_size(path):
    with open(path, "rb") as f:
        head = f.read(26)
    return struct.unpack(">II", head[16:24])  # PNG width,height


def fig(s, name, l, t, boxw, boxh):
    """Fit-and-center an image within (l,t,boxw,boxh)."""
    path = os.path.join(FIG, name) if not os.path.isabs(name) else name
    if not os.path.exists(path):
        path2 = os.path.join(OUT, name)
        path = path2 if os.path.exists(path2) else path
    iw, ih = _png_size(path); ar = iw / ih; boxar = boxw / boxh
    if ar > boxar:
        w = boxw; h = boxw / ar
    else:
        h = boxh; w = boxh * ar
    x = l + (boxw - w) / 2; y = t + (boxh - h) / 2
    s.shapes.add_picture(path, int(x), int(y), width=int(w), height=int(h))


def base(s):
    _bg(s); rect(s, 0, 0, SW, Inches(0.16), INK)
    rect(s, 0, SH - Inches(0.16), SW, Inches(0.16), BLUE)


def header(kicker, title):
    global _page
    s = prs.slides.add_slide(BLANK); base(s)
    txt(s, ML, Inches(0.40), CW, Inches(0.30), [[(kicker.upper(), 13, True, BLUE)]])
    txt(s, ML, Inches(0.68), CW, Inches(0.70), [[(title, 26, True, INK)]])
    rect(s, ML, Inches(1.44), Inches(1.6), Pt(3), TEAL)
    _page += 1
    txt(s, SW - Inches(1.1), SH - Inches(0.5), Inches(0.8), Inches(0.3),
        [[(str(_page), 10, False, GREY)]], align=PP_ALIGN.RIGHT)
    return s


def strip(s, text_runs, top=Inches(6.35), fill=INK):
    card(s, ML, top, CW, Inches(0.9), fill, line=None)
    txt(s, ML + Inches(0.4), top + Inches(0.06), CW - Inches(0.8), Inches(0.78),
        text_runs, anchor=MSO_ANCHOR.MIDDLE)


# =========================================================== 1 TITLE
s = prs.slides.add_slide(BLANK); _bg(s, INK)
rect(s, 0, 0, SW, Inches(0.22), BLUE); rect(s, 0, SH - Inches(0.22), SW, Inches(0.22), TEAL)
txt(s, Inches(0.9), Inches(1.6), Inches(11.5), Inches(0.4),
    [[("PROMOTION REVIEW  ·  RESEARCH SUMMARY", 14, True, RGBColor(0x9C, 0xB4, 0xD8))]])
txt(s, Inches(0.9), Inches(2.15), Inches(11.6), Inches(2.4),
    [[("Designing Serum-Tolerant Antifungal Candidates", 40, True, WHITE)],
     [("An AI- and Quantum-Chemistry–Guided Drug Discovery Program", 21, False, TEAL)],
     [("on the Papulacandin / Fusacandin class (β-1,3-glucan / FKS1 inhibitors)", 16, False, RGBColor(0xC8, 0xD2, 0xDC))]])
txt(s, Inches(0.9), Inches(5.7), Inches(11.6), Inches(1.0),
    [[("Presenter:  ______________________          Date:  ____________", 14, False, RGBColor(0xC8, 0xD2, 0xDC))],
     [("A complete, reproducible pipeline — from scattered literature to a ready-to-test design.", 13, False, RGBColor(0x9C, 0xB4, 0xD8), True)]])
_page = 1

# =========================================================== 2 PROBLEM
s = header("Background · the problem", "Invasive fungal infections: an urgent, under-served threat")
chips = [("Rising incidence", "hundreds of thousands of deaths / year", CORAL),
         ("Very few drug classes", "limited options, serious side effects", AMBER),
         ("Resistance spreading", "existing drugs losing ground", BLUE)]
cw = Inches(3.87); x = ML; y = Inches(1.85)
for head_, sub, col in chips:
    card(s, x, y, cw, Inches(1.65))
    rect(s, x, y, Inches(0.14), Inches(1.65), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x + Inches(0.35), y + Inches(0.26), cw - Inches(0.5), Inches(1.2),
        [[(head_, 17, True, INK)], [(sub, 13, False, GREY)]])
    x += cw + Inches(0.26)
card(s, ML, Inches(3.95), CW, Inches(2.35), CARD2, line=None)
txt(s, ML + Inches(0.4), Inches(4.2), CW - Inches(0.8), Inches(1.95),
    [[("The opportunity — and the catch", 17, True, BLUE)],
     [("The ", 15, False, INK), ("papulacandin", 15, True, INK),
      (" natural products block the enzyme that builds the fungal cell wall (β-1,3-glucan "
       "synthase / FKS1). Humans have no cell wall → an attractive, low-toxicity target.", 15, False, INK)],
     [("❗  The catch: potent in a test tube, but ", 15, False, CORAL),
      ("activity is lost in blood serum", 15, True, CORAL),
      (" — so they never became drugs.", 15, False, CORAL)],
     [("If we understand and fix that, we revive an entire drug class.", 14, True, INK)]])

# =========================================================== EVIDENCE BASE (fig)
s = header("Background · the data", "The evidence base: decades of scattered literature, curated")
fig(s, "fig_dataset.png", ML, Inches(1.75), CW, Inches(3.95))
strip(s, [[("A clean, provenance-tracked benchmark — itself a lasting resource for the field.",
            15, True, WHITE)]], top=Inches(6.05))

# =========================================================== 3 SERUM GAP (fig)
s = header("The scientific question", "The serum gap, quantified: blood switches most analogs off")
fig(s, "fig_serum_gap.png", ML, Inches(1.6), CW, Inches(4.35))
strip(s, [[("“serum shift”  =  dose needed in serum  ÷  dose needed in broth", 16, True, WHITE),
           ("      — the one number we set out to lower.", 14, False, RGBColor(0xC8, 0xD2, 0xDC))]])

# =========================================================== 4 STRATEGY (fig)
s = header("Strategy", "A 14-stage funnel: start broad and cheap, narrow to confident designs")
fig(s, "fig_pipeline.png", ML, Inches(1.6), CW, Inches(4.5))
strip(s, [[("Each stage checks our reasoning before the next — and the whole platform is reusable.",
            15, True, WHITE)]])

# =========================================================== 5 METHODS (fig+text)
s = header("Methods · the algorithm", "The staged quantum-chemistry funnel")
fig(s, "fig_qm_funnel.png", ML, Inches(1.55), Inches(6.7), Inches(5.4))
tx = ML + Inches(7.0)
txt(s, tx, Inches(2.0), Inches(4.9), Inches(4.6),
    [[("How it works", 17, True, BLUE)],
     [("①  AI invents novel molecules.", 15, True, INK)],
     [("②  Cheap filters remove the obvious losers.", 15, False, INK)],
     [("③  Accurate quantum chemistry computes the real 3-D shapes — and re-checks the "
       "cheap results.", 15, False, INK)],
     [("④  The costliest method (DFT) touches only the few survivors.", 15, False, INK)],
     [("", 8, False, INK)],
     [("Same measurement engine at every tier → a fair, honest comparison.", 14, True, TEAL)]])

# =========================================================== 6 HONESTY LADDER (fig) money
s = header("Key methodological result", "The “honesty ladder”: more rigor shrank the easy signal")
fig(s, "fig_honesty_ladder.png", ML, Inches(1.55), CW, Inches(4.35))
strip(s, [[("A less careful analysis would have made a confident — and WRONG — claim at every rung. ",
            14, True, WHITE), ("Knowing what is NOT true is a real, durable result.", 14, False, RGBColor(0xC8, 0xD2, 0xDC))]])

# =========================================================== CONFOUND (fig)
s = header("Result · the decisive check", "Why the exciting signal collapsed: it was mostly potency")
fig(s, "fig_confound.png", ML, Inches(1.6), CW, Inches(4.45))
strip(s, [[("Serum outcome is governed by intrinsic potency (ρ=0.79); remove it and the 3-D signal vanishes — ",
            14, True, WHITE), ("so we model the potency-free SHIFT.", 14, False, RGBColor(0xC8, 0xD2, 0xDC))]])

# =========================================================== ELECTRONICS (fig+text)
s = header("Result · independent check", "A different property (electronics) — the same honest story")
fig(s, "fig_electronics.png", ML, Inches(1.6), Inches(7.5), Inches(5.0))
tx = ML + Inches(7.75)
txt(s, tx, Inches(2.05), Inches(4.3), Inches(4.4),
    [[("What it shows", 16, True, BLUE)],
     [("•  Polarizability looks like the strongest signal — but it is a ", 13.5, False, INK),
      ("size / potency proxy", 13.5, True, INK), (" (vanishes on the shift).", 13.5, False, INK)],
     [("•  ", 13.5, False, INK), ("QM logP", 13.5, True, TEAL),
      (" independently corroborates the polar-surface lead (+0.30 on the shift).", 13.5, False, INK)],
     [("", 8, False, INK)],
     [("Two unrelated descriptor families converge on the same weak, honest conclusion.",
       13.5, True, INK)]])

# =========================================================== 7 ECHINOCANDIN (fig+text)
s = header("Learning from approved drugs", "Echinocandins: same target, real clinical data")
fig(s, "fig_echinocandin.png", ML, Inches(1.55), Inches(6.9), Inches(5.5))
tx = ML + Inches(7.15)
txt(s, tx, Inches(1.95), Inches(4.75), Inches(1.0),
    [[("Two lessons that redirected our design:", 16, True, BLUE)]])
c1 = card(s, tx, Inches(2.75), Inches(4.75), Inches(1.5), LIGHT, line=LINE)
rect(s, tx, Inches(2.75), Inches(0.12), Inches(1.5), CORAL, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, tx + Inches(0.3), Inches(2.9), Inches(4.3), Inches(1.25),
    [[("1.  The tail is REQUIRED for potency", 14.5, True, INK)],
     [("can’t just make it water-friendly or cut it off, or the drug stops working.", 13, False, GREY)]])
c2 = card(s, tx, Inches(4.45), Inches(4.75), Inches(1.5), LIGHT, line=LINE)
rect(s, tx, Inches(4.45), Inches(0.12), Inches(1.5), GREEN, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
txt(s, tx + Inches(0.3), Inches(4.6), Inches(4.3), Inches(1.25),
    [[("2.  SHAPE matters more than greasiness", 14.5, True, INK)],
     [("flexible tail = less serum loss. Our native tail is rigid & flat → the ‘bad’ zone.", 13, False, GREY)]])

# =========================================================== 8 AI DESIGN + REWARD (fig+text)
s = header("AI-guided design", "We built an AI generator — and validated its scoring on real data")
fig(s, os.path.join(OUT, "phase12_reward_validation.png"), ML, Inches(1.65), Inches(6.4), Inches(4.9))
tx = ML + Inches(6.8)
txt(s, tx, Inches(2.0), Inches(5.0), Inches(4.6),
    [[("The AI, done honestly", 17, True, BLUE)],
     [("•  An AI generator invents novel, on-scaffold molecules.", 14.5, False, INK)],
     [("•  Its ‘good design’ reward was ", 14.5, False, INK),
      ("checked against real measured data", 14.5, True, INK),
      (" (left) before we trusted it.", 14.5, False, INK)],
     [("•  It is built to produce a ", 14.5, False, INK),
      ("discriminating series", 14.5, True, TEAL),
      (" — molecules that TEST the idea, not just assume it.", 14.5, False, INK)],
     [("", 8, False, INK)],
     [("Left: the reward tracks true serum tolerance on the 24 known compounds "
       "(a weak-but-real trend, honestly reported).", 12.5, False, GREY, True)]])

# =========================================================== TAIL RANKING (fig)
s = header("Result · round 1 (detail)", "AI-designed tails, ranked by the cheap 3-D screen")
fig(s, "fig_tail_ranking.png", ML, Inches(1.55), CW, Inches(4.5))
strip(s, [[("Only 2 of 12 beat the native tail here — and the accurate re-check (next slide) removed even those.",
            14.5, True, WHITE)]])

# =========================================================== 9 ROUND-1 RESULT (fig)
s = header("Result · round 1", "The quantum funnel caught a false positive before synthesis")
fig(s, "fig_round1.png", ML, Inches(1.6), CW, Inches(4.55))
strip(s, [[("The accurate method exposed the cheap method’s ‘winner’ as an artifact — ", 14, True, WHITE),
           ("caught before any lab time was spent.", 14, False, RGBColor(0xC8, 0xD2, 0xDC))]])

# =========================================================== 10 SHORTLIST (fig+text)
s = header("The deliverable", "A ranked, easy-to-make shortlist — redirected by the drug lesson")
fig(s, "fig_tail_ladder.png", ML, Inches(1.55), Inches(7.0), Inches(4.6))
tx = ML + Inches(7.2)
txt(s, tx, Inches(1.85), Inches(4.7), Inches(0.6), [[("Make first (best bet):", 16, True, BLUE)]])
for i, (nm, sub, col) in enumerate([("Palmitoyl (saturated)", "cleanest single test", GREEN),
                                    ("Branched saturated", "caspofungin-like", GREEN),
                                    ("One-kink control", "make second", TEAL)]):
    yy = Inches(2.5) + i * Inches(0.82)
    card(s, tx, yy, Inches(4.7), Inches(0.72), RGBColor(0xE7, 0xF4, 0xEA) if col == GREEN else LIGHT, line=col)
    txt(s, tx + Inches(0.28), yy + Inches(0.05), Inches(4.3), Inches(0.6),
        [[(nm, 14, True, INK), (f"   ·  {sub}", 12.5, False, GREY)]], anchor=MSO_ANCHOR.MIDDLE)
strip(s, [[("Happy coincidence:  ", 15, True, TEAL),
           ("these are also the EASIEST to synthesize", 15, True, WHITE),
           (" — cheap fatty acids, one step. Best science = cheapest chemistry.", 14, False, RGBColor(0xC8, 0xD2, 0xDC))]])

# =========================================================== 11 CONTRIBUTIONS (grid)
s = header("Contributions", "What this program delivered — a reusable platform, not a one-off")
cards = [("📊  Benchmark dataset", "138 compounds, 1,042 records,\nfull provenance — a field resource"),
         ("🔬  Validated QM method", "multi-tier quantum pipeline +\nin-house 3-D surface engine"),
         ("🤖  AI molecule generator", "novel-structure design with a\nreward checked against real data"),
         ("💊  Cross-drug read-across", "clinical echinocandin evidence\nimported into the design logic"),
         ("⚙️  Reproducible engineering", "14 documented stages, automated\ntests, cluster run protocols"),
         ("🗺️  Roadmap + shortlist", "the exact experiment that\nunlocks predictive AI")]
cols = 3; gap = Inches(0.25); cw = (CW - gap * (cols - 1)) / cols; ch = Inches(2.05)
for i, (head_, body) in enumerate(cards):
    r, c = divmod(i, cols)
    x = ML + c * (cw + gap); y = Inches(1.9) + r * (ch + Inches(0.22))
    card(s, x, y, cw, ch, LIGHT, line=LINE)
    rect(s, x, y, cw, Pt(4), BLUE)
    txt(s, x + Inches(0.25), y + Inches(0.22), cw - Inches(0.5), ch - Inches(0.4),
        [[(head_, 15, True, INK)], [(body, 12.5, False, GREY)]])

# =========================================================== 12 LESSON + ROADMAP
s = header("The key lesson & the road ahead", "Where AI helps — and the loop that gets us there")
card(s, ML, Inches(1.85), CW, Inches(1.5), CARD2, line=None)
txt(s, ML + Inches(0.4), Inches(2.0), CW - Inches(0.8), Inches(1.25),
    [[("This is NOT “AI failed.”  ", 15, True, BLUE), ("It is a precise statement of the critical path:", 15, False, INK)],
     [("AI drug design needs a measured ‘answer key’. Physics-based guesses alone were shown — rigorously — to mislead.", 14, False, INK)],
     [("We now know exactly which few experiments unlock the AI — and the AI is already built to use them.", 14, True, GREEN)]])
loop = [("1.  Synthesize the\n4-molecule shortlist", BLUE), ("2.  Measure serum shift\n+ protein binding", AMBER),
        ("3.  Train the AI\non the results", GREEN), ("4.  AI designs\nround 2 (wider)", TEAL)]
n = len(loop); gap = Inches(0.55); cw = (CW - gap * (n - 1)) / n; y = Inches(4.1); x = ML
for i, (t, col) in enumerate(loop):
    card(s, x, y, cw, Inches(1.5), LIGHT, line=col)
    rect(s, x, y, Inches(0.12), Inches(1.5), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x + Inches(0.28), y, cw - Inches(0.4), Inches(1.5), [[(t, 13.5, True, INK)]], anchor=MSO_ANCHOR.MIDDLE)
    x += cw
    if i < n - 1:
        chevron(s, x + Inches(0.12), y + Inches(0.55), w=Inches(0.3), h=Inches(0.4), c=col)
    x += gap
txt(s, ML, Inches(5.85), CW, Inches(0.5),
    [[("↻  A virtuous loop: each round of data makes the AI a better designer.  ", 13.5, True, TEAL),
      ("The immediate next step needs no more computing.", 13.5, False, GREY)]], align=PP_ALIGN.CENTER)

# =========================================================== 13 SIGNIFICANCE
s = header("Significance", "Why this matters — and the one-line takeaway")
pts = [("Revived a stalled, high-value drug class as a tractable design problem", "urgent unmet clinical need"),
       ("Built an end-to-end AI + quantum-chemistry discovery platform", "reusable across programs"),
       ("Demonstrated rigorous self-correction (the ‘honesty ladder’)", "trustworthy, independent science"),
       ("Integrated chemistry + machine learning + quantum physics + pharmacology", "genuinely interdisciplinary"),
       ("Defined the exact experiment that unlocks predictive AI", "de-risks future investment")]
y = Inches(1.85)
for head_, sub in pts:
    rect(s, ML, y + Inches(0.06), Inches(0.16), Inches(0.16), TEAL, shape=MSO_SHAPE.OVAL)
    txt(s, ML + Inches(0.35), y, CW - Inches(0.4), Inches(0.6),
        [[(head_ + "  ", 14.5, True, INK), ("— " + sub, 13, False, GREY, True)]])
    y += Inches(0.62)
card(s, ML, Inches(5.2), CW, Inches(1.5), INK, line=None)
txt(s, ML + Inches(0.45), Inches(5.35), CW - Inches(0.9), Inches(1.25),
    [[("“We set out to teach a computer to design a serum-tolerant antifungal. We built the platform, "
       "curated the data, imported lessons from approved drugs, and — most importantly — learned to tell "
       "a real signal from a flattering illusion.", 14, False, WHITE, True)],
     [("The computer is ready; it awaits one clean experiment, designed down to the exact four molecules.”",
       14, True, TEAL, True)]])

# =========================================================== 14 THANK YOU
s = prs.slides.add_slide(BLANK); _bg(s, INK)
rect(s, 0, 0, SW, Inches(0.22), BLUE); rect(s, 0, SH - Inches(0.22), SW, Inches(0.22), TEAL)
txt(s, Inches(0.9), Inches(2.7), Inches(11.6), Inches(1.5),
    [[("Thank you", 44, True, WHITE)], [("Questions & discussion welcome", 20, False, TEAL)]])
txt(s, Inches(0.9), Inches(6.2), Inches(11.6), Inches(0.6),
    [[("Full technical detail, code, and reproducible reports: project repository (Phases 1–14).",
       13, False, RGBColor(0x9C, 0xB4, 0xD8), True)]])

path = os.path.join(OUT, "Papulacandin_promotion_report.pptx")
prs.save(path)
print(f"Wrote {path}  ({len(prs.slides._sldIdLst)} slides)")
