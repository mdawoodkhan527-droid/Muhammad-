#!/usr/bin/env python3
"""Create a publication-ready workflow diagram for the observation-operator study.

The script writes workflow_observation_operator.png (300 dpi) and
workflow_observation_operator.pdf (vector PDF) beside this script.
"""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


# Output files are written beside this script so it remains portable.
OUTPUT_DIR = Path(__file__).resolve().parent
PNG_PATH = OUTPUT_DIR / "workflow_observation_operator.png"
PDF_PATH = OUTPUT_DIR / "workflow_observation_operator.pdf"

# Muted, colorblind-friendly academic palette.
BOX_FACE = "#E8F0F2"
BOX_EDGE = "#2F4F4F"
ARROW_COLOR = "#333333"
TITLE_COLOR = "#1F2933"
OUTPUT_COLOR = "#52616B"

# Use common sans-serif fonts with a reliable Matplotlib fallback.
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def build_flowchart() -> tuple[plt.Figure, plt.Axes]:
    """Draw the left-to-right workflow and return the figure and axes."""
    fig, ax = plt.subplots(figsize=(7.0, 3.0))
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.98)
    ax.set_xlim(0, 7)
    ax.set_ylim(0, 3)
    ax.axis("off")

    # The title is wrapped intentionally to preserve a comfortable font size.
    ax.text(
        3.5,
        2.68,
        "No fitting: admissibility and information content are assessed\nbefore calibration",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        color=TITLE_COLOR,
        linespacing=1.15,
    )

    labels = [
        "Published\nevidence\n(21)\nobservation\ngroups",
        "Synthetic\n1D model\ntransport +\nreactions",
        "Verification\n+ one-way\nmechanics",
        "Observation\noperator",
        "Structural\ngate\nthen SVD",
    ]

    # Uniform box geometry keeps the flow visually balanced.
    box_width = 1.18
    box_height = 1.08
    gap = 0.22
    first_x = 0.11
    box_y = 1.06
    center_y = box_y + box_height / 2
    xs = [first_x + i * (box_width + gap) for i in range(len(labels))]

    for x, label in zip(xs, labels):
        box = FancyBboxPatch(
            (x, box_y),
            box_width,
            box_height,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.2,
            edgecolor=BOX_EDGE,
            facecolor=BOX_FACE,
            antialiased=True,
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(
            x + box_width / 2,
            center_y,
            label,
            ha="center",
            va="center",
            fontsize=10.0,
            color=TITLE_COLOR,
            linespacing=1.18,
            zorder=3,
        )

    # Draw arrows in the open gaps, from box edge to box edge.
    for left_x, right_x in zip(xs[:-1], xs[1:]):
        arrow = FancyArrowPatch(
            (left_x + box_width + 0.015, center_y),
            (right_x - 0.015, center_y),
            arrowstyle="-|>",
            mutation_scale=9,
            linewidth=1.5,
            color=ARROW_COLOR,
            shrinkA=0,
            shrinkB=0,
            zorder=4,
        )
        ax.add_patch(arrow)

    ax.text(
        3.5,
        0.43,
        "Output: verified benchmark, rank result, and measurements needed to add information",
        ha="center",
        va="center",
        fontsize=10,
        fontstyle="italic",
        color=OUTPUT_COLOR,
    )

    return fig, ax


def main() -> None:
    """Generate PNG and vector PDF outputs."""
    fig, _ = build_flowchart()
    fig.savefig(PNG_PATH, dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(PDF_PATH, format="pdf", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"Saved PNG: {PNG_PATH}")
    print(f"Saved PDF: {PDF_PATH}")


if __name__ == "__main__":
    main()
