#!/usr/bin/env python3
"""Create a publication-ready verification summary figure.

Outputs: verification_figure.png (300 dpi) and verification_figure.pdf (vector).
The error ratios are computed from the reported benchmark errors and frozen
acceptance tolerances in the manuscript.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatterMathtext, LogLocator

OUT = Path(__file__).resolve().parent
PNG = OUT / "verification_figure.png"
PDF = OUT / "verification_figure.pdf"

# Okabe-Ito-inspired, colorblind-friendly colours.
INK = "#263238"
MUTED = "#56636B"
GRID = "#D9E1E5"
BLUE = "#0072B2"
ORANGE = "#D55E00"
PALE = "#EDF3F5"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": INK,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

# Reported relative errors and their predeclared acceptance tolerances.
labels = [
    "Transport: relative L₂ error",
    "Transport: maximum error",
    "Elastic bar: displacement L₂ error",
    "Elastic bar: equilibrium residual",
]
errors = np.array([1.8605e-5, 2.0245e-5, 1.3165e-15, 3.3284e-14])
tolerances = np.array([0.02, 0.03, 0.005, 1.0e-10])
ratios = errors / tolerances
colors = [BLUE, BLUE, ORANGE, ORANGE]

fig = plt.figure(figsize=(7.4, 3.45), constrained_layout=False)
grid = fig.add_gridspec(
    nrows=1, ncols=2,
    width_ratios=(1.45, 1.0),
    left=0.19, right=0.98, bottom=0.20, top=0.82, wspace=0.32,
)
ax = fig.add_subplot(grid[0, 0])
summary = fig.add_subplot(grid[0, 1])

# Panel (a): error-to-tolerance ratios. Values below one satisfy the limit.
y = np.arange(len(labels))[::-1]
ax.set_xscale("log")
ax.set_xlim(1e-14, 2.0)
ax.set_ylim(-0.55, 3.55)
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8.4)
ax.xaxis.set_major_locator(LogLocator(base=10, numticks=7))
ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
ax.grid(axis="x", which="major", color=GRID, linewidth=0.7, zorder=0)
ax.grid(axis="x", which="minor", color=GRID, linewidth=0.35, alpha=0.55, zorder=0)
ax.axvline(1.0, color=INK, linestyle=(0, (4, 3)), linewidth=1.2, zorder=1)
ax.text(0.985, 1.015, "acceptance limit", transform=ax.get_xaxis_transform(),
        ha="right", va="bottom", fontsize=8, color=MUTED)

for yi, value, color in zip(y, ratios, colors):
    ax.hlines(yi, 1e-14, value, color=color, linewidth=2.0, alpha=0.65, zorder=2)
    ax.scatter(value, yi, s=38, color=color, edgecolor="white", linewidth=0.8, zorder=3)
    ax.annotate(f"{value:.2e}", (value, yi), xytext=(6, 0),
                textcoords="offset points", va="center", ha="left",
                fontsize=8, color=INK)

ax.set_xlabel("Relative error / declared tolerance", fontsize=8.7, labelpad=7)
ax.set_title("(a) Analytical limiting cases", loc="left", fontsize=10.2,
             fontweight="bold", pad=10, color=INK)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_color("#7D8B91")
    ax.spines[spine].set_linewidth(0.8)
ax.tick_params(axis="y", length=0, pad=7, labelsize=8.4)
ax.tick_params(axis="x", labelsize=7.8, width=0.7)

# Panel (b): verification facts, grouped with clear typographic hierarchy.
summary.set_xlim(0, 1)
summary.set_ylim(0, 1)
summary.axis("off")
summary.set_title("(b) Reaction and conservation checks", loc="left",
                  fontsize=10.2, fontweight="bold", pad=10, color=INK)

summary.text(0.02, 0.78, "6 / 6", fontsize=14, fontweight="bold", color=BLUE,
             ha="left", va="center")
summary.text(0.45, 0.80, "reaction vectors balanced", fontsize=8.5,
             fontweight="bold", color=INK, ha="left", va="center")
summary.text(0.45, 0.66, "for elements and charge", fontsize=8.3,
             color=MUTED, ha="left", va="center")
summary.plot([0.0, 1.0], [0.55, 0.55], color=GRID, linewidth=0.8, clip_on=False)

summary.text(0.02, 0.41, "12 / 12", fontsize=14, fontweight="bold", color=BLUE,
             ha="left", va="center")
summary.text(0.49, 0.44, "baseline cases passed", fontsize=8.5,
             fontweight="bold", color=INK, ha="left", va="center")
summary.text(0.49, 0.30, "stepwise element, charge,", fontsize=8.0,
             color=MUTED, ha="left", va="center")
summary.text(0.49, 0.19, "and pore-volume checks", fontsize=8.0,
             color=MUTED, ha="left", va="center")

summary.text(0.00, 0.13, "Defect corrected", fontsize=8.4,
             fontweight="bold", color=ORANGE, ha="left", va="center")
summary.text(0.00, 0.025,
             "Negative-inventory clipping was detected,\ncorrected, and the baseline cases rerun.",
             fontsize=7.3, color=MUTED, ha="left", va="center", linespacing=1.15)

fig.savefig(PNG, dpi=300, bbox_inches="tight", pad_inches=0.08)
fig.savefig(PDF, format="pdf", bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"Saved {PNG}")
print(f"Saved {PDF}")
