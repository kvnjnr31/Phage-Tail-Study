"""
This program provides lightweight helper functions to be used by data_violinplots.py 
"""
 

from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, Optional
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


@dataclass
class ViolinStyle:
    facealpha: float = 0.35
    edgealpha: float = 0.9
    linewidth: float = 1.2
    point_size: float = 22.0
    cap_width: float = 0.12
    box_linewidth: float = 1.25


def _draw_points(ax: plt.Axes, y: np.ndarray, x: float, width: float, s: ViolinStyle):
    if y.size == 0:
        return
    rng = np.random.default_rng(7)
    jitter = (rng.random(y.size) - 0.5) * (width * 0.6)
    ax.scatter(np.full_like(y, x) + jitter, y, s=s.point_size,
               c="k", edgecolors="k", linewidths=0.3, zorder=3)


def _draw_box(ax: plt.Axes, y: np.ndarray, x: float, s: ViolinStyle):
    if y.size == 0:
        return
    q1, med, q3 = np.percentile(y, [25, 50, 75]) if y.size >= 3 else (y.min(), np.median(y), y.max())
    ax.plot([x - s.cap_width, x + s.cap_width], [med, med],
            color="k", lw=s.box_linewidth, zorder=4)
    ax.plot([x, x], [q1, q3], color="k", lw=s.box_linewidth, zorder=4)


def draw_violin(
    ax: plt.Axes,
    y: Iterable[float],
    x: float,
    *,
    width: float = 0.9,
    facecolor: Optional[str] = None,
    edgecolor: Optional[str] = None,
    alpha: float = 0.35,
    show_points: bool = True,
    show_box: bool = True,
    style: ViolinStyle = ViolinStyle(),
    logy: bool = False,
    y_floor: float = 1e-9,
):
    """Draw a single violin at position x with proper log support."""
    y = np.asarray([v for v in y if np.isfinite(v) and v > 0], dtype=float)
    if y.size == 0:
        return

    unique_vals = np.unique(y)
    can_kde = (unique_vals.size >= 2) and (np.nanstd(y) > 0)

    if can_kde:
        if logy:
            ly = np.log10(np.clip(y, y_floor, None))
            grid = np.linspace(ly.min(), ly.max(), 200)
            kde = gaussian_kde(ly)
            dens = kde(grid)
            dens = dens / dens.max() * (width / 2)
            for sign in (-1, 1):
                ax.fill_betweenx(10 ** grid,
                                 x + sign * dens,
                                 x,
                                 color=facecolor or "C0",
                                 alpha=alpha,
                                 edgecolor=edgecolor or "C0",
                                 linewidth=style.linewidth)
        else:
            parts = ax.violinplot(y, positions=[x], widths=width,
                                  showmeans=False, showmedians=False, showextrema=False)
            body = parts["bodies"][0]
            if facecolor is not None:
                body.set_facecolor(facecolor)
            if edgecolor is not None:
                body.set_edgecolor(edgecolor)
            body.set_alpha(alpha)
            body.set_linewidth(style.linewidth)
    else:
        y0 = float(np.median(y))
        hh = max(1e-12, 0.03 * (np.nanmax(y) - np.nanmin(y) if can_kde else max(y0, 1.0)))
        poly = plt.Polygon([[x - width/8, y0 - hh], [x + width/8, y0 - hh],
                            [x + width/8, y0 + hh], [x - width/8, y0 + hh]],
                           facecolor=facecolor or "C0", edgecolor=edgecolor or "C0",
                           alpha=alpha, linewidth=style.linewidth)
        ax.add_patch(poly)

    if show_points:
        _draw_points(ax, y, x, width * 0.45, style)
    if show_box:
        _draw_box(ax, y, x, style)
