#!/usr/bin/env python3
"""Publication figures for the pgx-suite validation manuscript.

Data are hard-coded from docs/AxiomValidation.md (the numbers are script-generated
upstream by build_axiom_concordance.py; reproduced here for figure rendering).

Run:  python3 paper/figures.py      # writes PNGs (+SVG) to paper/figures/
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(OUT, exist_ok=True)

# Okabe-Ito colorblind-safe palette
CB = {"blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
      "vermillion": "#D55E00", "sky": "#56B4E9", "grey": "#999999"}
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150})


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)


def fig2_waterfall():
    """Adjudication waterfall: 73.0% raw -> 85.8% adjudicated."""
    fig, ax = plt.subplots(figsize=(7, 4))
    labels = ["Concordant\n(raw)", "+ array_FN", "+ array_auth.", "Residual"]
    vals = [241, 30, 12, 47]
    colors = [CB["green"], CB["sky"], CB["blue"], CB["grey"]]
    bottom = 0
    for lab, v, c in zip(labels, vals, colors):
        ax.bar(0.5, v, bottom=bottom, width=0.6, color=c, label=f"{lab} ({v})")
        ax.text(0.5, bottom + v / 2, str(v), ha="center", va="center",
                color="white", fontweight="bold")
        bottom += v
    ax.axhline(241, color=CB["green"], ls="--", lw=1)
    ax.axhline(283, color=CB["blue"], ls="--", lw=1)
    ax.text(0.95, 241, "  73.0% (241/330)", va="center", color=CB["green"])
    ax.text(0.95, 283, "  85.8% (283/330)", va="center", color=CB["blue"])
    ax.set_xlim(0, 1.6)
    ax.set_ylim(0, 340)
    ax.set_xticks([])
    ax.set_ylabel("Gene-calls (of 330)")
    ax.set_title("Verdict concordance with Axiom array, before and after adjudication")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    save(fig, "figure2_waterfall")


def fig3_per_gene():
    """Per-gene concordance, sorted descending."""
    data = [("ABCG2", 100), ("CYP2B6", 100), ("CYP3A4", 100), ("CYP3A5", 100),
            ("NAT2", 100), ("TPMT", 100), ("CYP2C19", 92), ("CYP2C9", 92),
            ("NUDT15", 92), ("SLCO1B1", 92), ("CYP4F2", 58), ("VKORC1", 54),
            ("UGT1A1", 42), ("CYP2D6", 17), ("DPYD", 8), ("RYR1", 0)]
    data.sort(key=lambda x: x[1], reverse=True)
    genes = [d[0] for d in data]
    vals = [d[1] for d in data]
    colors = [CB["green"] if v >= 90 else CB["orange"] if v >= 50 else CB["vermillion"]
              for v in vals]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(genes, vals, color=colors)
    for i, v in enumerate(vals):
        ax.text(v + 1, i, f"{v}%", va="center", fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 108)
    ax.set_xlabel("Raw verdict concordance (%)")
    ax.set_title("Per-gene concordance (16 array-overlapping genes)")
    ax.axvline(90, color=CB["grey"], ls=":", lw=1)
    save(fig, "figure3_per_gene")


def fig4_platform():
    """Platform equivalence: MGI vs Illumina."""
    fig, ax = plt.subplots(figsize=(5, 4))
    plats = ["MGI", "Illumina"]
    pct = [73.9, 72.1]
    n = ["122/165", "119/165"]
    bars = ax.bar(plats, pct, width=0.5, color=[CB["blue"], CB["orange"]])
    for b, p, nn in zip(bars, pct, n):
        ax.text(b.get_x() + b.get_width() / 2, p + 0.5, f"{p}%\n({nn})",
                ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 85)
    ax.set_ylabel("Raw verdict concordance (%)")
    ax.set_title("Platform equivalence (1.8 pp apart)")
    save(fig, "figure4_platform")


def fig1_schematic():
    """Four-tier reconciliation engine schematic."""
    fig, ax = plt.subplots(figsize=(6, 6.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis("off")
    tiers = [
        ("All caller diplotypes + phenotypes", CB["grey"], 11),
        ("Tier 1 — Coverage gate\n(NO_CALL if depth < 10x)", CB["sky"], 9),
        ("Tier 2 — Synonym collapse\n(allele_synonyms.json)", CB["green"], 7),
        ("Tier 3 — Authoritative caller\n(Cyrius / PharmCAT / VCF-Check)", CB["blue"], 5),
        ("Tier 4 — Phenotype consensus\n(>=2 callers agree on CPIC phenotype)", CB["orange"], 3),
        ("One verdict per gene\n(+ authority label)", CB["vermillion"], 1),
    ]
    for text, color, y in tiers:
        box = FancyBboxPatch((1.5, y - 0.55), 7, 1.1, boxstyle="round,pad=0.1",
                             fc=color, ec="black", alpha=0.85)
        ax.add_patch(box)
        ax.text(5, y, text, ha="center", va="center", color="white",
                fontsize=9.5, fontweight="bold")
    for y0, y1 in [(11, 9), (9, 7), (7, 5), (5, 3), (3, 1)]:
        ax.add_patch(FancyArrowPatch((5, y0 - 0.55), (5, y1 + 0.55),
                     arrowstyle="-|>", mutation_scale=15, color="black"))
    ax.set_title("Reconciliation engine", fontsize=12, fontweight="bold")
    save(fig, "figure1_schematic")


if __name__ == "__main__":
    fig1_schematic()
    fig2_waterfall()
    fig3_per_gene()
    fig4_platform()
    # sanity: adjudication categories must sum to 330
    assert 241 + 30 + 12 + 6 + 34 + 7 == 330, "Table 4 categories must sum to 330"
    assert 241 + 30 + 12 == 283, "Adjudicated concordant must be 283"
    print(f"Wrote figures 1-4 (png+svg) to {OUT}/  [sanity checks passed]")
