#!/usr/bin/env python3
"""Publication-ready rank-count comparison for paired LHS samples."""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
PNG = OUT / "rank_counts_figure.png"
PDF = OUT / "rank_counts_figure.pdf"

INK = "#263238"
MUTED = "#56636B"
GRID = "#D9E1E5"
BLUE = "#0072B2"
ORANGE = "#D55E00"

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

n_samples = 32
rank_labels = ["Numerical rank = 3", "Practical rank = 1"]
y = np.array([1, 0])

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.25), sharex=True, sharey=True)
fig.subplots_adjust(left=0.18, right=0.98, bottom=0.25, top=0.78, wspace=0.20)
fig.suptitle("Rank classification across the paired coverage design",
             fontsize=10.8, fontweight="bold", color=INK, y=0.96)

for ax, title in zip(axes, ["Kozeny–Carman", "Power law"]):
    ax.barh(y, [n_samples, n_samples], height=0.52,
            color=[BLUE, ORANGE], alpha=0.82, zorder=2)
    for yi in y:
        ax.text(n_samples - 0.35, yi, "32/32  ·  100%", ha="right", va="center",
                fontsize=8.2, fontweight="bold", color="white", zorder=3)
    ax.set_title(title, fontsize=9.5, fontweight="bold", color=INK, pad=9)
    ax.set_xlim(0, 34)
    ax.set_ylim(-0.65, 1.65)
    ax.set_yticks(y, rank_labels, fontsize=8.2)
    ax.set_xticks([0, 8, 16, 24, 32])
    ax.grid(axis="x", color=GRID, linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#7D8B91")
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="y", length=0, pad=7)
    ax.tick_params(axis="x", labelsize=7.8, width=0.7)

axes[0].set_ylabel("Rank classification", fontsize=8.5, labelpad=9)
fig.supxlabel("Paired LHS design samples (of 32)", fontsize=8.6, y=0.075)

fig.savefig(PNG, dpi=300, bbox_inches="tight", pad_inches=0.08)
fig.savefig(PDF, format="pdf", bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"Saved {PNG}")
print(f"Saved {PDF}")
