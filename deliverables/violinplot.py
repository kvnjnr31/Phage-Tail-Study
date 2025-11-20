"""
This program makes special violin plots — a graph that shows how data are spread out.
Each “violin” appears as a wide, vertical shape that gets thicker where more values occur.
It can handle very small or very large numbers safely by using a log scale, to facilitate
data that span many orders of magnitude. The function 'plot_violins' draws several violins
side by side, labels and colors them for comparing groups of measurements on the same chart.
"""

from __future__ import annotations
from typing import Iterable, Optional, Sequence, Dict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext, LogLocator
from violin import draw_violin

def _to_1d_float_array(vals) -> np.ndarray:
    """Coerce to 1D float array; return empty array if not iterable."""
    try:
        arr = np.asarray(vals, dtype=float).ravel()
        if arr.ndim != 1:
            arr = arr.flatten()
        return arr
    except Exception:
        return np.array([], dtype=float)

def plot_violins(
    ax: plt.Axes,
    series: Sequence[Optional[Iterable[float]]],
    labels: Sequence[Optional[str]],
    *,
    colors: Optional[Sequence[Optional[str]]] = None,
    positions: Optional[Sequence[float]] = None,
    logy: bool = True,
    width: float = 0.9,
    title: str = "",
    xlabel: str = "",
    ylabel: str = "",
    y_floor: float = 1e-9,
    show_points: bool = True,
    
) -> Dict[str, float]:
    """Draw a row of violins with log-safe rendering."""
    assert len(series) == len(labels)
    n = len(series)
    if positions is None:
        positions = list(np.arange(1, n + 1, dtype=float))
    if colors is None:
        colors = [None] * n

    label_to_x: Dict[str, float] = {}
    xpos_iter = iter(positions)

    for idx, (vals, lab, col) in enumerate(zip(series, labels, colors)):
        x = next(xpos_iter)

        # Treat Ellipsis as spacer
        if lab is Ellipsis:
            lab = None
        if vals is Ellipsis:
            vals = None

        if lab is not None:
            label_to_x[lab] = x
        if vals is None:
            continue

        v = _to_1d_float_array(vals)
        v = v[np.isfinite(v) & (v > 0)]
        if v.size == 0:
            continue

        draw_violin(ax, v, x, facecolor=col, edgecolor=col,
                    logy=logy, y_floor=y_floor)

    ax.set_xticks(positions)
    ax.set_xticklabels([l if l is not None else "" for l in labels], rotation=20, ha="right")

    if logy:
        ax.set_yscale("log")
        ax.set_ylim(y_floor, ax.get_ylim()[1])
        ax.yaxis.set_major_locator(LogLocator(base=10.0))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.tick_params(axis="both", labelsize=10)
    return label_to_x
