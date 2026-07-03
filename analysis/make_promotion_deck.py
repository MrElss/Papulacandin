#!/usr/bin/env python3
"""
make_promotion_deck.py
======================
Build a professional, polished promotion-report deck (16:9) summarizing the
Papulacandin serum-gap program (Phases 1-14) for a general-faculty audience.

Self-contained: all framework diagrams (funnel, "honesty ladder", echinocandin
quadrant, shortlist ladder, roadmap loop, deliverables grid) are drawn with
shapes; one real result figure (phase13_qm_ranking.png) is embedded.

Output: analysis/outputs/Papulacandin_promotion_report.pptx
"""
import os

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "outputs")

# ---- palette -------------------------------------------------------------
INK   = RGBColor(0x14, 0x23, 0x3A)   # near-black navy (titles/text)
BLUE  = RGBColor(0x2F, 0x56, 0xA6)   # primary accent
TEAL  = RGBColor(0x1F, 0x8A, 0x8A)   # secondary
CORAL = RGBColor(0xD8, 0x4A, 0x3C)   # negative / problem
GREEN = RGBColor(0x2E, 0x8B, 0x57)   # positive
AMBER = RGBColor(0xD8, 0x8A, 0x2B)   # caution / tempting
GREY  = RGBColor(0x5B, 0x64, 0x72)   # secondary text
LIGHT = RGBColor(0xF2, 0xF5, 0xF9)   # card background
CARD2 = RGBColor(0xE8, 0xEE, 0xF6)   # alt card
LINE  = RGBColor(0xD3, 0xDA, 0xE4)   # hairlines
WHITE = RGBColor(0xFF, 0xFF, 0xFF)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
SW, SH = prs.slide_width, prs.slide_height
BLANK = prs.slide_layouts[6]
ML = Inches(0.75)                    # left margin
CW = SW - 2 * ML                     # content width

_page = 0


def _solid(shape, color):
    shape.fill.solid(); shape.fill.fore_color.rgb = color
    shape.line.fill.background()


def _bg(slide, color=WHITE):
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = color


def rect(slide, l, t, w, h, color, shape=MSO_SHAPE.RECTANGLE, line=None, lw=1.0):
    s = slide.shapes.add_shape(shape, l, t, w, h)
    _solid(s, color)
    if line is not None:
        s.line.color.rgb = line; s.line.width = Pt(lw); s.line.fill.solid(); s.line.fill.fore_color.rgb = line
    s.shadow.inherit = False
    return s


def txt(slide, l, t, w, h, runs, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, wrap=True):
    """runs: list of paragraphs; each paragraph is list of (text, size, bold, color, italic)."""
    tb = slide.shapes.add_textbox(l, t, w, h)
    tf = tb.text_frame; tf.word_wrap = wrap; tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Pt(2); tf.margin_top = tf.margin_bottom = Pt(1)
    for i, para in enumerate(runs):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(4)
        for (text, size, bold, color, *rest) in para:
            r = p.add_run(); r.text = text
            r.font.size = Pt(size); r.font.bold = bold; r.font.color.rgb = color
            r.font.name = "Calibri"
            if rest and rest[0]:
                r.font.italic = True
    return tb


def base(slide):
    _bg(slide)
    rect(slide, 0, 0, SW, Inches(0.16), INK)                 # top accent
    rect(slide, 0, SH - Inches(0.16), SW, Inches(0.16), BLUE)  # bottom accent


def header(slide, kicker, title):
    base(slide)
    txt(slide, ML, Inches(0.42), CW, Inches(0.32),
        [[(kicker.upper(), 13, True, BLUE)]])
    txt(slide, ML, Inches(0.70), CW, Inches(0.75),
        [[(title, 27, True, INK)]])
    rect(slide, ML, Inches(1.48), Inches(1.6), Pt(3), TEAL)
    global _page
    _page += 1
    txt(slide, SW - Inches(1.1), SH - Inches(0.52), Inches(0.8), Inches(0.3),
        [[(str(_page), 10, False, GREY)]], align=PP_ALIGN.RIGHT)


def new(kicker, title):
    s = prs.slides.add_slide(BLANK); header(s, kicker, title); return s


def chevron(slide, l, t, w=Inches(0.34), h=Inches(0.34), color=GREY):
    rect(slide, l, t, w, h, color, shape=MSO_SHAPE.CHEVRON)


def card(slide, l, t, w, h, fill=LIGHT, line=LINE):
    return rect(slide, l, t, w, h, fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE, line=line, lw=1.0)


# =========================================================================
# 1 — TITLE
# =========================================================================
s = prs.slides.add_slide(BLANK); _bg(s, INK)
rect(s, 0, 0, SW, Inches(0.22), BLUE)
rect(s, 0, SH - Inches(0.22), SW, Inches(0.22), TEAL)
txt(s, Inches(0.9), Inches(1.7), Inches(11.5), Inches(0.4),
    [[("PROMOTION REVIEW  ·  RESEARCH SUMMARY", 14, True, RGBColor(0x9C, 0xB4, 0xD8))]])
txt(s, Inches(0.9), Inches(2.2), Inches(11.6), Inches(2.4),
    [[("Designing Serum-Tolerant Antifungal Candidates", 40, True, WHITE)],
     [("An AI- and Quantum-Chemistry–Guided Drug Discovery Program", 21, False, TEAL)],
     [("on the Papulacandin / Fusacandin class (β-1,3-glucan / FKS1 inhibitors)", 16, False, RGBColor(0xC8, 0xD2, 0xDC))]])
txt(s, Inches(0.9), Inches(5.7), Inches(11.6), Inches(1.0),
    [[("Presenter:  ______________________          Date:  ____________", 14, False, RGBColor(0xC8, 0xD2, 0xDC))],
     [("A complete, reproducible pipeline — from scattered literature to a ready-to-test design.", 13, False, RGBColor(0x9C, 0xB4, 0xD8), True)]])
_page = 1

# =========================================================================
# 2 — THE PROBLEM
# =========================================================================
s = new("The problem", "Invasive fungal infections: an urgent, under-served threat")
chips = [("Rising incidence", "hundreds of thousands of deaths / year", CORAL),
         ("Very few drug classes", "limited options, serious side effects", AMBER),
         ("Resistance spreading", "existing drugs losing ground", BLUE)]
cw = Inches(3.9); gap = Inches(0.25); x = ML; y = Inches(1.9)
for head, sub, col in chips:
    c = card(s, x, y, cw, Inches(1.7))
    rect(s, x, y, Inches(0.14), Inches(1.7), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x + Inches(0.35), y + Inches(0.28), cw - Inches(0.5), Inches(1.2),
        [[(head, 17, True, INK)], [(sub, 13, False, GREY)]])
    x += cw + gap
card(s, ML, Inches(4.0), CW, Inches(2.4), CARD2, line=None)
txt(s, ML + Inches(0.4), Inches(4.25), CW - Inches(0.8), Inches(2.0),
    [[("The opportunity — and the catch", 17, True, BLUE)],
     [("The ", 15, False, INK), ("papulacandin", 15, True, INK),
      (" natural products block the enzyme that builds the fungal cell wall "
       "(β-1,3-glucan synthase / FKS1). Humans have no cell wall → an attractive, "
       "low-toxicity target.", 15, False, INK)],
     [("❗  The catch: these molecules are potent in a test tube but ", 15, False, CORAL),
      ("lose activity in blood serum", 15, True, CORAL),
      (" — so they never became drugs.", 15, False, CORAL)],
     [("If we understand and fix that, we revive an entire drug class.", 14, True, INK)]])

# =========================================================================
# 3 — THE PUZZLE + METRIC
# =========================================================================
s = new("The scientific question", "Why does blood switch the drug off — and can we design around it?")
# two panels
p1 = card(s, ML, Inches(1.9), Inches(5.7), Inches(2.5), RGBColor(0xE7, 0xF4, 0xEA), line=GREEN)
txt(s, ML + Inches(0.35), Inches(2.15), Inches(5.1), Inches(2.0),
    [[("🧪  In clean broth", 18, True, GREEN)],
     [("Drug + fungus  →  fungus dies", 15, False, INK)],
     [("Potent. Low dose works.", 14, False, GREY, True)]])
p2 = card(s, ML + Inches(6.0), Inches(1.9), Inches(5.85), Inches(2.5), RGBColor(0xFB, 0xE8, 0xE6), line=CORAL)
txt(s, ML + Inches(6.35), Inches(2.15), Inches(5.2), Inches(2.0),
    [[("🩸  In blood serum", 18, True, CORAL)],
     [("Same drug + fungus  →  fungus survives", 15, False, INK)],
     [("Activity lost. The clinical failure mode.", 14, False, GREY, True)]])
# metric band
card(s, ML, Inches(4.7), CW, Inches(1.9), INK, line=None)
txt(s, ML + Inches(0.45), Inches(4.95), CW - Inches(0.9), Inches(1.5),
    [[("We defined ONE number to explain and to design against — the ", 15, False, WHITE),
      ("serum shift", 16, True, TEAL)],
     [("serum shift  =  (dose needed in serum)  ÷  (dose needed in broth)", 17, True, WHITE)],
     [("shift ≈ 1  →  serum has no effect (ideal)      large shift  →  serum wrecks the drug", 14, False, RGBColor(0xC8, 0xD2, 0xDC))],
     [("Goal: lower the serum shift without losing potency.", 14, True, TEAL)]])

# =========================================================================
# 4 — STRATEGY FUNNEL
# =========================================================================
s = new("Strategy", "A staged funnel: start broad and cheap, narrow to high-confidence designs")
acts = [("ACT I", "Data & first clues", "Curate literature → 2-D patterns\n& interpretable design rules", "Phases 1–4", BLUE),
        ("ACT II", "Go 3-D (quantum)", "AI generates molecules;\nfull 3-D shape & electronics", "Phases 5–9", TEAL),
        ("ACT III", "Learn from drugs", "Test 'albumin sponge';\nborrow echinocandin lessons", "Phases 10–11", AMBER),
        ("ACT IV", "AI design + test", "AI designs candidates;\nquantum funnel filters them", "Phases 12–14", GREEN)]
n = len(acts); gap = Inches(0.2)
cw = (CW - gap * (n - 1)) / n
x = ML; y = Inches(2.2)
for i, (act, head, body, ph, col) in enumerate(acts):
    c = card(s, x, y, cw, Inches(2.9))
    rect(s, x, y, cw, Inches(0.5), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x, y + Inches(0.04), cw, Inches(0.42), [[(act, 14, True, WHITE)]], align=PP_ALIGN.CENTER)
    txt(s, x + Inches(0.18), y + Inches(0.62), cw - Inches(0.36), Inches(2.1),
        [[(head, 15, True, INK)],
         [(body, 12.5, False, GREY)],
         [(ph, 11, True, col)]])
    x += cw
    if i < n - 1:
        chevron(s, x - Inches(0.02), y + Inches(1.25), w=Inches(0.24), h=Inches(0.34), color=GREY)
    x += gap
card(s, ML, Inches(5.5), CW, Inches(1.0), CARD2, line=None)
txt(s, ML + Inches(0.4), Inches(5.62), CW - Inches(0.8), Inches(0.8),
    [[("Outcome:  ", 15, True, BLUE),
      ("a ranked, synthesizable shortlist + a precise experimental roadmap — and a reusable platform for future questions.", 15, False, INK)]],
    anchor=MSO_ANCHOR.MIDDLE)

# =========================================================================
# 5 — THE HONESTY LADDER (money slide)
# =========================================================================
s = new("The intellectual backbone", "The “honesty ladder”: more rigor → the easy signal shrank")
txt(s, ML, Inches(1.65), CW, Inches(0.5),
    [[("At every step we used a more rigorous method — and reported what was actually there, "
       "not the flattering early number.", 14, False, GREY)]])
rungs = [("2-D chemistry", 0.30, "weak — promising", GREY),
         ("Cheap 3-D (approximate)", 0.62, "looks strong — TEMPTING", AMBER),
         ("Accurate quantum (real 3-D shapes)", 0.30, "shrinks — corrected", BLUE),
         ("…controlling for plain potency", 0.06, "≈ zero — artifact removed", CORAL),
         ("Best quantum re-check of AI designs", 0.05, "no winner — decision → experiment", CORAL)]
y = Inches(2.35); rowh = Inches(0.62); barx = ML + Inches(3.95); barmax = Inches(4.35)
for label, frac, verdict, col in rungs:
    txt(s, ML + Inches(0.35), y, Inches(3.5), rowh, [[(label, 13, True, INK)]], anchor=MSO_ANCHOR.MIDDLE)
    rect(s, barx, y + Inches(0.10), barmax, Inches(0.34), RGBColor(0xEC, 0xEF, 0xF3))  # track
    rect(s, barx, y + Inches(0.10), Emu(int(barmax * frac)), Inches(0.34), col,
         shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, barx + barmax + Inches(0.12), y, Inches(2.75), rowh,
        [[(verdict, 11.5, False, col)]], anchor=MSO_ANCHOR.MIDDLE)
    y += rowh
# arrow rail (increasing rigor, downward)
rect(s, ML + Inches(0.02), Inches(2.42), Pt(3), Inches(2.85), LINE)
txt(s, ML - Inches(0.05), Inches(2.35), Inches(0.35), Inches(3.0),
    [[("▼", 12, True, GREY)]], align=PP_ALIGN.CENTER)
card(s, ML, Inches(5.75), CW, Inches(0.85), INK, line=None)
txt(s, ML + Inches(0.4), Inches(5.85), CW - Inches(0.8), Inches(0.7),
    [[("A less careful analysis would have made a confident — and wrong — claim at every rung. "
       "Knowing precisely what is NOT true is a real, durable result.", 14, True, WHITE)]],
    anchor=MSO_ANCHOR.MIDDLE)

# =========================================================================
# 6 — ECHINOCANDIN INSIGHT (quadrant)
# =========================================================================
s = new("Learning from approved drugs", "Echinocandins: same target, real clinical data — shape beats greasiness")
# quadrant box
qx, qy, qw, qh = ML, Inches(1.95), Inches(5.6), Inches(4.4)
rect(s, qx, qy, qw, qh, RGBColor(0xF7, 0xF9, 0xFC), line=LINE, shape=MSO_SHAPE.RECTANGLE)
rect(s, qx, qy + qh/2, qw, Pt(1.2), LINE)
rect(s, qx + qw/2, qy, Pt(1.2), qh, LINE)
txt(s, qx, qy + qh + Inches(0.02), qw, Inches(0.3),
    [[("flexible / saturated tail        →        rigid / aromatic tail", 11, True, GREY)]], align=PP_ALIGN.CENTER)
txt(s, qx - Inches(0.15), qy, Inches(0.3), qh,
    [[("loses ← serum activity → keeps", 11, True, GREY)]], anchor=MSO_ANCHOR.MIDDLE)
def dot(cx, cy, label, col):
    d = Inches(0.26)
    rect(s, qx + Emu(int(qw*cx)) - d/2, qy + Emu(int(qh*(1-cy))) - d/2, d, d, col, shape=MSO_SHAPE.OVAL)
    txt(s, qx + Emu(int(qw*cx)) - Inches(0.9), qy + Emu(int(qh*(1-cy))) + Inches(0.12), Inches(1.8), Inches(0.3),
        [[(label, 10.5, True, col)]], align=PP_ALIGN.CENTER)
dot(0.20, 0.85, "Caspofungin ✓", GREEN)
dot(0.80, 0.35, "Anidulafungin", CORAL)
dot(0.78, 0.12, "Micafungin", CORAL)
dot(0.72, 0.20, "Our native tail", INK)
# right column: table + lessons
tx = ML + Inches(6.0); tw = Inches(5.85)
txt(s, tx, Inches(1.95), tw, Inches(0.4), [[("Approved-drug pattern (from our own data)", 15, True, BLUE)]])
rows = [("Caspofungin", "branched, flexible", "small  (~2–11×)", GREEN),
        ("Anidulafungin", "rigid, aromatic", "large  (~16–512×)", CORAL),
        ("Micafungin", "rigid, aromatic", "large  (~64–1024×)", CORAL)]
yy = Inches(2.45)
for name, tail, shift, col in rows:
    card(s, tx, yy, tw, Inches(0.62), LIGHT, line=LINE)
    rect(s, tx, yy, Inches(0.1), Inches(0.62), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, tx + Inches(0.3), yy + Inches(0.06), tw - Inches(0.5), Inches(0.5),
        [[(f"{name}", 13, True, INK), (f"   ·   {tail}   ·   serum loss ", 12, False, GREY),
          (shift, 12.5, True, col)]], anchor=MSO_ANCHOR.MIDDLE)
    yy += Inches(0.72)
card(s, tx, Inches(4.75), tw, Inches(1.6), CARD2, line=None)
txt(s, tx + Inches(0.3), Inches(4.9), tw - Inches(0.6), Inches(1.35),
    [[("Two lessons that redirected our design:", 13.5, True, BLUE)],
     [("1.  The tail is REQUIRED for potency — can’t just make it water-friendly.", 13, False, INK)],
     [("2.  SHAPE matters more than greasiness — flexible tail = less serum loss.", 13, False, INK)],
     [("Our native tail is rigid & flat → in the ‘bad’ zone.", 12.5, True, CORAL)]])

# =========================================================================
# 7 — ROUND 1 RESULT (embed real figure) + funnel-caught-artifact
# =========================================================================
s = new("Round 1 result", "AI designed the tails; the quantum funnel caught a false positive")
steps = [("12 designed tails", GREY), ("Cheap screen: 2 look best", AMBER),
         ("Accurate quantum: 0 of 3 hold up", CORAL)]
x = ML; bw = Inches(3.55); y = Inches(1.85)
for i, (t, col) in enumerate(steps):
    card(s, x, y, bw, Inches(0.8), LIGHT, line=col)
    txt(s, x + Inches(0.2), y, bw - Inches(0.4), Inches(0.8), [[(t, 13.5, True, col)]], anchor=MSO_ANCHOR.MIDDLE)
    x += bw
    if i < len(steps) - 1:
        chevron(s, x - Inches(0.05), y + Inches(0.24), w=Inches(0.28), h=Inches(0.32), color=GREY)
    x += Inches(0.05)
img = os.path.join(OUT, "phase13_qm_ranking.png")
if os.path.exists(img):
    s.shapes.add_picture(img, ML, Inches(2.95), height=Inches(3.3))
txt(s, ML + Inches(6.7), Inches(3.0), Inches(5.1), Inches(3.4),
    [[("What happened", 15, True, BLUE)],
     [("The cheap method’s ‘winner’ (a water-friendly group) ", 13, False, INK),
      ("folds up and hides in accurate 3-D", 13, True, CORAL),
      (" — an artifact.", 13, False, INK)],
     [("We caught it BEFORE any lab synthesis.", 13.5, True, GREEN)],
     [("", 6, False, INK)],
     [("Conclusion: the computed surface score does not, by itself, "
       "reliably rank these tails. The decision moves to experiment.", 13, False, INK)]])

# =========================================================================
# 8 — THE SHORTLIST (what to make first)
# =========================================================================
s = new("The deliverable", "A ranked, easy-to-make shortlist — redirected by the drug lesson")
txt(s, ML, Inches(1.65), CW, Inches(0.5),
    [[("New idea (from echinocandins): don’t make the tail water-friendly — ", 14, False, INK),
      ("make it less rigid", 14, True, GREEN),
      (".  Same length (keeps potency), only the shape changes.", 14, False, INK)]])
native = card(s, ML, Inches(2.6), Inches(2.9), Inches(2.2), RGBColor(0xFB, 0xE8, 0xE6), line=CORAL)
txt(s, ML + Inches(0.25), Inches(2.8), Inches(2.5), Inches(1.9),
    [[("Native tail", 15, True, CORAL)], [("rigid, flat (a polyene)", 12.5, False, GREY)],
     [("= the problem", 12.5, True, CORAL)], [("reference baseline", 11.5, False, GREY, True)]])
chevron(s, ML + Inches(2.95), Inches(3.5), w=Inches(0.3), h=Inches(0.4), color=GREY)
targets = [("Palmitoyl", "straight, flexible", "★ make first", GREEN),
           ("Branched saturated", "caspofungin-like", "★ make first", GREEN),
           ("One-kink version", "rigidity control", "make second", TEAL)]
x = ML + Inches(3.5); tw = Inches(2.85)
for name, sub, tag, col in targets:
    card(s, x, Inches(2.6), tw, Inches(2.2), RGBColor(0xE7, 0xF4, 0xEA) if col == GREEN else LIGHT, line=col)
    txt(s, x + Inches(0.22), Inches(2.8), tw - Inches(0.44), Inches(1.9),
        [[(name, 15, True, INK)], [(sub, 12.5, False, GREY)],
         [("", 4, False, INK)], [(tag, 12.5, True, col)]])
    x += tw + Inches(0.15)
card(s, ML, Inches(5.15), CW, Inches(1.3), INK, line=None)
txt(s, ML + Inches(0.4), Inches(5.3), CW - Inches(0.8), Inches(1.05),
    [[("A happy coincidence:  ", 15, True, TEAL),
      ("these are also the EASIEST to synthesize", 15, True, WHITE),
      (" — cheap, off-the-shelf fatty acids in one step,", 14, False, RGBColor(0xC8, 0xD2, 0xDC))],
     [("unlike the original rigid tail or exotic water-friendly tails. ", 14, False, RGBColor(0xC8, 0xD2, 0xDC)),
      ("Best scientific bet = cheapest chemistry.", 14, True, WHITE)]])

# =========================================================================
# 9 — WHAT I BUILT
# =========================================================================
s = new("Contributions", "What this program delivered — a reusable platform, not a one-off")
cards = [("📊  Benchmark dataset", "138 compounds, 1,042 records,\nfull provenance — a field resource"),
         ("🔬  Validated QM method", "multi-tier quantum pipeline +\nin-house 3-D surface engine"),
         ("🤖  AI molecule generator", "novel-structure design with a\nreward checked against real data"),
         ("💊  Cross-drug read-across", "clinical echinocandin evidence\nimported into the design logic"),
         ("⚙️  Reproducible engineering", "14 documented stages, automated\ntests, cluster run protocols"),
         ("🗺️  Roadmap + shortlist", "the exact experiment that\nunlocks predictive AI")]
cols = 3; gap = Inches(0.25)
cw = (CW - gap * (cols - 1)) / cols; ch = Inches(2.1)
for i, (head, body) in enumerate(cards):
    r, c = divmod(i, cols)
    x = ML + c * (cw + gap); y = Inches(1.95) + r * (ch + Inches(0.25))
    card(s, x, y, cw, ch, LIGHT, line=LINE)
    rect(s, x, y, cw, Pt(4), BLUE, shape=MSO_SHAPE.RECTANGLE)
    txt(s, x + Inches(0.25), y + Inches(0.22), cw - Inches(0.5), ch - Inches(0.4),
        [[(head, 15, True, INK)], [(body, 12.5, False, GREY)]])

# =========================================================================
# 10 — WHERE AI HELPS + ROADMAP
# =========================================================================
s = new("The key lesson & the road ahead", "Where AI helps — and the loop that gets us there")
# lesson strip
card(s, ML, Inches(1.9), CW, Inches(1.55), CARD2, line=None)
txt(s, ML + Inches(0.4), Inches(2.05), CW - Inches(0.8), Inches(1.3),
    [[("This is NOT “AI failed.”  ", 15, True, BLUE),
      ("It is a precise statement of the critical path:", 15, False, INK)],
     [("AI drug design needs a measured ‘answer key’. Physics-based guesses alone were shown "
       "— rigorously — to mislead.", 14, False, INK)],
     [("We now know exactly which few experiments unlock the AI — and the AI is already built to use them.", 14, True, GREEN)]])
# roadmap loop
loop = [("1.  Synthesize the\n4-molecule shortlist", BLUE),
        ("2.  Measure serum shift\n+ protein binding", AMBER),
        ("3.  Train the AI\non the results", GREEN),
        ("4.  AI designs\nround 2 (wider)", TEAL)]
n = len(loop); gap = Inches(0.55); cw = (CW - gap * (n - 1)) / n; y = Inches(4.2)
x = ML
for i, (t, col) in enumerate(loop):
    card(s, x, y, cw, Inches(1.5), LIGHT, line=col)
    rect(s, x, y, Inches(0.12), Inches(1.5), col, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    txt(s, x + Inches(0.28), y, cw - Inches(0.4), Inches(1.5), [[(t, 13.5, True, INK)]], anchor=MSO_ANCHOR.MIDDLE)
    x += cw
    if i < n - 1:
        chevron(s, x + Inches(0.12), y + Inches(0.55), w=Inches(0.3), h=Inches(0.4), color=col)
    x += gap
txt(s, ML, Inches(5.95), CW, Inches(0.5),
    [[("↻  A virtuous loop: each round of data makes the AI a better designer.  ", 13.5, True, TEAL),
      ("Immediate next step needs no more computing.", 13.5, False, GREY)]], align=PP_ALIGN.CENTER)

# =========================================================================
# 11 — SIGNIFICANCE / CLOSE
# =========================================================================
s = new("Significance", "Why this matters for the field — and the one-line takeaway")
pts = [("Revived a stalled, high-value drug class as a tractable design problem",
        "addresses an urgent unmet clinical need"),
       ("Built an end-to-end AI + quantum-chemistry discovery platform",
        "reusable across natural-product programs"),
       ("Demonstrated rigorous self-correction (the ‘honesty ladder’)",
        "the hallmark of trustworthy, independent science"),
       ("Integrated chemistry + machine learning + quantum physics + pharmacology",
        "genuinely interdisciplinary"),
       ("Defined the exact experiment that unlocks predictive AI",
        "de-risks and directs future investment")]
y = Inches(1.9)
for head, sub in pts:
    rect(s, ML, y + Inches(0.06), Inches(0.16), Inches(0.16), TEAL, shape=MSO_SHAPE.OVAL)
    txt(s, ML + Inches(0.35), y, CW - Inches(0.4), Inches(0.62),
        [[(head + "  ", 14.5, True, INK), ("— " + sub, 13, False, GREY, True)]])
    y += Inches(0.66)
card(s, ML, Inches(5.4), CW, Inches(1.4), INK, line=None)
txt(s, ML + Inches(0.45), Inches(5.55), CW - Inches(0.9), Inches(1.15),
    [[("“We set out to teach a computer to design a serum-tolerant antifungal. We built the platform, "
       "curated the data, imported lessons from approved drugs, and — most importantly — learned to tell "
       "a real signal from a flattering illusion.", 14, False, WHITE, True)],
     [("The computer is ready; it awaits one clean experiment, which we have designed down to the exact four molecules.”",
       14, True, TEAL, True)]])

# =========================================================================
# 12 — THANK YOU
# =========================================================================
s = prs.slides.add_slide(BLANK); _bg(s, INK)
rect(s, 0, 0, SW, Inches(0.22), BLUE)
rect(s, 0, SH - Inches(0.22), SW, Inches(0.22), TEAL)
txt(s, Inches(0.9), Inches(2.7), Inches(11.6), Inches(1.5),
    [[("Thank you", 44, True, WHITE)],
     [("Questions & discussion welcome", 20, False, TEAL)]])
txt(s, Inches(0.9), Inches(6.2), Inches(11.6), Inches(0.6),
    [[("Full technical detail, code, and reproducible reports: project repository (Phases 1–14).",
       13, False, RGBColor(0x9C, 0xB4, 0xD8), True)]])

path = os.path.join(OUT, "Papulacandin_promotion_report.pptx")
prs.save(path)
print(f"Wrote {path}  ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)")
