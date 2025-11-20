# tails_heads_violin.py
# Template: plots ONLY the violin plots of phage TAIL and HEAD size distributions.
# Output: figure_violin_tails_heads.png
# ────────────────────────────────────────────────────────────────────────────────
# HOW TO USE
# 1) Replace the lists under DATA with your measured sizes in nanometers (nm).
# 2) Run:  python tails_heads_violin.py
# 3) The figure will be saved next to the script.

from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from typing import Dict, List, Tuple

# ================================ CONFIG =======================================
TITLE              = "(A) Phage Tail & Head Distributions"
YLABEL             = "Size (nm)"
OUTPUT_PNG         = "figure_violin_tails_heads.png"
LOG_SCALE_Y        = True
POINTS_SHOW        = True
POINTS_SIZE        = 18
POINTS_ALPHA       = 0.9
POINTS_JITTER      = 0.08          # horizontal jitter for scatter points
VIOLIN_ALPHA       = 0.35
VIOLIN_LINEWIDTH   = 1.2
MEDIAN_LINEWIDTH   = 2.0
XTICK_ROTATE_DEG   = 0

# Colors by family (feel free to change)
COL_PODO = "#4F81BD"   # blue
COL_MYO  = "#9BBB59"   # green
COL_SIPHO= "#C0504D"   # red

# ================================ DATA =========================================
# Enter sizes in nanometers (nm). Replace the example arrays below with your data.
# If a family lacks data for a category, use an empty list [].

DATA: Dict[str, Dict[str, List[float]]] = {
    "Podoviridae": {
        "Tail (nm)": [20.0, 15.0, 12.0, 10.0, 18.0, 13.0, 17.0, 46.0],             # ← replace
        "Head (nm)": [60.0, 70.0, 63.0, 70.0, 61.0, 60.0, 67.0, 60.0],             # ← replace
    },
    "Myoviridae": {
        "Tail (nm)": [113.0, 170.0, 135.0, 120.0, 150.0, 170.0, 160.0, 53.0],      # ← replace
        "Head (nm)": [85.0, 60.0, 60.0, 94.0, 64.0, 144.5, 85.1],                  # ← replace
    },
    "Siphoviridae": {
        "Tail (nm)": [150.0, 151.0, 180.0, 135.0, 98.0, 210.0, 100.0, 140.0, 355], # ← replace
        "Head (nm)": [60.0, 60.0, 80.0, 60.0, 55.0, 53.0, 60.0, 81.0],             # ← replace
    },
}

FAMILY_ORDER: List[Tuple[str, str]] = [
    ("Podoviridae", COL_PODO),
    ("Myoviridae",  COL_MYO),
    ("Siphoviridae",COL_SIPHO),
]
CATEGORY_ORDER = ["Tail (nm)", "Head (nm)"]  # left→right within each family

# ============================== PLOT LOGIC =====================================
def _violin(ax, data: List[np.ndarray], positions: List[float], color: str):
    v = ax.violinplot(
        data, positions=positions, widths=0.7,
        showmeans=False, showmedians=False, showextrema=False
    )
    for b in v["bodies"]:
        b.set_facecolor(color)
        b.set_alpha(VIOLIN_ALPHA)
        b.set_edgecolor("k")
        b.set_linewidth(VIOLIN_LINEWIDTH)
    # draw medians as thick lines
    for xs, ys in zip(positions, data):
        if len(ys) == 0:
            continue
        q2 = np.median(ys)
        ax.plot([xs - 0.3, xs + 0.3], [q2, q2], lw=MEDIAN_LINEWIDTH, color="k")

def _scatter(ax, y: np.ndarray, x: float, color: str):
    if len(y) == 0:
        return
    jitter = (np.random.rand(len(y)) - 0.5) * 2 * POINTS_JITTER
    ax.scatter(
        np.full_like(y, x, dtype=float) + jitter, y,
        s=POINTS_SIZE, alpha=POINTS_ALPHA, c=color, edgecolors="k", linewidths=0.4, zorder=3
    )

def main():
    plt.close("all")
    fig, ax = plt.subplots(figsize=(8.5, 4.5))

    x_positions = []
    x_labels    = []
    x_current   = 1.0
    family_patches = []

    for fam, color in FAMILY_ORDER:
        fam_block_positions = []
        for cat in CATEGORY_ORDER:
            y = np.asarray(DATA.get(fam, {}).get(cat, []), dtype=float)
            _violin(ax, [y], [x_current], color=color)
            if POINTS_SHOW:
                _scatter(ax, y, x_current, color)
            x_positions.append(x_current)
            x_labels.append(cat)
            fam_block_positions.append(x_current)
            x_current += 1.2  # spacing within the family

        # legend handle per family
        family_patches.append(Patch(facecolor=color, alpha=VIOLIN_ALPHA, label=fam))
        x_current += 0.8  # extra gap between families

    ax.set_xticks(x_positions)
    ax.set_xticklabels(x_labels, rotation=XTICK_ROTATE_DEG)
    ax.set_ylabel(YLABEL)
    ax.set_title(TITLE)

    if LOG_SCALE_Y:
        ax.set_yscale("log")

    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    ax.legend(handles=family_patches, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300)
    plt.show()

if __name__ == "__main__":
    np.random.seed(7)  # deterministic jitter
    main()
