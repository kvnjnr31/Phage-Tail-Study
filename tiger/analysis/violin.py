"""Log-safe violin plots for adsorption rates and phage sizes.

Merges the 2025 helpers ``violin.py`` (one violin) and ``violinplot.py`` (a row of
violins) from the Phage-Tail-Study deliverables (originals in
``archive/violin_study_2025/``). Behaviour is unchanged: the kernel density is
estimated in log10 space when ``logy`` is set, groups with fewer than two distinct
values are drawn as a small box, and each violin carries jittered raw points and a
median / interquartile bar.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatterMathtext, LogLocator
from scipy.stats import gaussian_kde

FAMILY_COLORS = {"Podoviridae": "#4F81BD", "Myoviridae": "#9BBB59", "Siphoviridae": "#C0504D"}
HOST_ABBREV = {"Enterobacteriaceae": "Entero.", "Synechococcaceae": "Syn.",
               "Vibrionaceae": "Vib.", "Pseudomonadaceae": "Psd."}


@dataclass
class ViolinStyle:
    facealpha: float = 0.35
    edgealpha: float = 0.9
    linewidth: float = 1.2
    point_size: float = 22.0
    cap_width: float = 0.12
    box_linewidth: float = 1.25


def _positive(vals) -> np.ndarray:
    try:
        v = np.asarray(vals, dtype=float).ravel()
    except Exception:
        return np.array([], dtype=float)
    return v[np.isfinite(v) & (v > 0)]


def _draw_points(ax, y, x, width, s: ViolinStyle):
    rng = np.random.default_rng(7)
    jitter = (rng.random(y.size) - 0.5) * (width * 0.6)
    ax.scatter(np.full_like(y, x) + jitter, y, s=s.point_size,
               c="k", edgecolors="k", linewidths=0.3, zorder=3)


def _draw_box(ax, y, x, s: ViolinStyle):
    q1, med, q3 = np.percentile(y, [25, 50, 75]) if y.size >= 3 else (y.min(), np.median(y), y.max())
    ax.plot([x - s.cap_width, x + s.cap_width], [med, med], color="k", lw=s.box_linewidth, zorder=4)
    ax.plot([x, x], [q1, q3], color="k", lw=s.box_linewidth, zorder=4)


def draw_violin(ax: plt.Axes, y: Iterable[float], x: float, *, width: float = 0.9,
                facecolor: Optional[str] = None, edgecolor: Optional[str] = None,
                alpha: float = 0.35, show_points: bool = True, show_box: bool = True,
                style: ViolinStyle = ViolinStyle(), logy: bool = False, y_floor: float = 1e-9):
    """Draw one violin at position ``x``; the density is computed in log10 space if ``logy``."""
    y = _positive(y)
    if y.size == 0:
        return
    fc, ec = facecolor or "C0", edgecolor or "C0"
    if np.unique(y).size >= 2 and np.nanstd(y) > 0:
        if logy:
            ly = np.log10(np.clip(y, y_floor, None))
            grid = np.linspace(ly.min(), ly.max(), 200)
            dens = gaussian_kde(ly)(grid)
            dens = dens / dens.max() * (width / 2)
            for sign in (-1, 1):
                ax.fill_betweenx(10 ** grid, x + sign * dens, x, color=fc, alpha=alpha,
                                 edgecolor=ec, linewidth=style.linewidth)
        else:
            body = ax.violinplot(y, positions=[x], widths=width, showmeans=False,
                                 showmedians=False, showextrema=False)["bodies"][0]
            body.set_facecolor(fc); body.set_edgecolor(ec)
            body.set_alpha(alpha); body.set_linewidth(style.linewidth)
    else:
        y0 = float(np.median(y))
        hh = max(1e-12, 0.03 * max(y0, 1.0))
        ax.add_patch(plt.Polygon([[x - width / 8, y0 - hh], [x + width / 8, y0 - hh],
                                  [x + width / 8, y0 + hh], [x - width / 8, y0 + hh]],
                                 facecolor=fc, edgecolor=ec, alpha=alpha, linewidth=style.linewidth))
    if show_points:
        _draw_points(ax, y, x, width * 0.45, style)
    if show_box:
        _draw_box(ax, y, x, style)


def plot_violins(ax: plt.Axes, series: Sequence[Optional[Iterable[float]]],
                 labels: Sequence[Optional[str]], *, colors: Optional[Sequence[Optional[str]]] = None,
                 positions: Optional[Sequence[float]] = None, logy: bool = True, width: float = 0.9,
                 title: str = "", xlabel: str = "", ylabel: str = "", y_floor: float = 1e-9,
                 show_points: bool = True) -> Dict[str, float]:
    """Draw a row of violins. ``None`` in ``series`` leaves a gap. Returns label -> x position."""
    assert len(series) == len(labels)
    n = len(series)
    positions = list(np.arange(1, n + 1, dtype=float)) if positions is None else list(positions)
    colors = [None] * n if colors is None else colors
    label_to_x: Dict[str, float] = {}
    for vals, lab, col, x in zip(series, labels, colors, positions):
        if lab:
            label_to_x[lab] = x
        if vals is None:
            continue
        draw_violin(ax, vals, x, facecolor=col, edgecolor=col, logy=logy, y_floor=y_floor,
                    show_points=show_points, width=width)
    ax.set_xticks(positions)
    ax.set_xticklabels([l or "" for l in labels], rotation=20, ha="right")
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
