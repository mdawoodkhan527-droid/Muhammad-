#!/usr/bin/env python3
"""Publication-ready summary of the benchmark-Ba localisation screen."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT = Path(__file__).resolve().parent
PNG = OUT / "barite_localisation_figure.png"
PDF = OUT / "barite_localisation_figure.pdf"

INK = "#263238"
MUTED = "#56636B"
GRID = "#D9E1E5"
BLUE = "#0072B2"
ORANGE = "#D55E00"
PALE_BLUE = "#DDEBF2"
PALE_ORANGE = "#F7E6D4"

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

bone_range = (0.870, 0.921)
matrix_range = (0.890, 0.931)

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.55), sharex=True)
fig.subplots_adjust(left=0.09, right=0.98, bottom=0.28, top=0.78, wspace=0.24)
fig.suptitle("Benchmark-Ba barite-localisation screen", fontsize=10.8,
             fontweight="bold", color=INK, y=0.96)

panels = [
    (axes[0], "(a) Bone-proxy detection", bone_range, 0.80, BLUE, PALE_BLUE,
     "Required minimum: 0.80", "≥"),
    (axes[1], "(b) Matrix-proxy detection", matrix_range, 0.20, ORANGE, PALE_ORANGE,
     "Required maximum: 0.20", "≤"),
]
for ax, title, interval, threshold, color, pale, threshold_label, relation in panels:
    ax.set_title(title, loc="left", fontsize=9.3, fontweight="bold", pad=9, color=INK)
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.5, 0.8)
    ax.set_yticks([])
    ax.set_xticks([0, .2, .4, .6, .8, 1.0], ["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
    ax.grid(axis="x", color=GRID, linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)
    if relation == "≥":
        ax.axvspan(threshold, 1.0, color=pale, alpha=0.85, zorder=0)
    else:
        ax.axvspan(0, threshold, color=pale, alpha=0.85, zorder=0)
    ax.axvline(threshold, color=INK, linestyle=(0, (4, 2.5)), linewidth=1.15, zorder=2)
    lo, hi = interval
    ax.hlines(0.15, lo, hi, color=color, linewidth=4.0, zorder=3)
    ax.vlines([lo, hi], 0.02, 0.28, color=color, linewidth=1.4, zorder=3)
    ax.scatter([lo, hi], [0.15, 0.15], s=28, color=color,
               edgecolor="white", linewidth=0.7, zorder=4)
    ax.text((lo + hi) / 2, 0.43, f"{lo:.3f}–{hi:.3f}", ha="center", va="center",
            fontsize=8.8, fontweight="bold", color=INK)
    ax.text(threshold + (0.015 if relation == "≥" else 0.02), 0.68,
            threshold_label, ha="left", va="center", fontsize=7.4, color=MUTED,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor=GRID, linewidth=0.7))
    ax.set_xlabel("Prior-averaged detection probability", fontsize=8.2, labelpad=6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#7D8B91")
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="x", labelsize=7.7, width=0.7)
    ax.tick_params(axis="y", length=0)

fig.text(0.5, 0.105,
         "Across three benchmark-Ba model forms, no variant meets both criteria. "
         "The low-Ba stress case succeeds only outside the defensible screening range "
         "and within a narrow source-rate window.",
         ha="center", va="center", fontsize=7.7, color=MUTED, wrap=True)

fig.savefig(PNG, dpi=300, bbox_inches="tight", pad_inches=0.08)
fig.savefig(PDF, format="pdf", bbox_inches="tight", pad_inches=0.08)
plt.close(fig)
print(f"Saved {PNG}")
print(f"Saved {PDF}")
