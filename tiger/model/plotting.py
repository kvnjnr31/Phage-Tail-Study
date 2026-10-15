"""Plots for the receptor simulation.

Ports simulation_plotter.m, sim_plotter.m (radial concentration profile) and the
trajectory figure at the end of encounter_rates_main.m.
"""
from __future__ import annotations

import numpy as np

from .geometry import plot_capsule


def radial_concentration(pos: np.ndarray, rmax: float = 2.5, bins: int = 40):
    """Shell-normalised radial concentration of phage around the origin (sim_plotter.m)."""
    r = np.linalg.norm(pos, axis=1)
    r = r[r <= rmax]
    n, edges = np.histogram(r, bins=bins)
    centers = 0.5 * (edges[1:] + edges[:-1])
    half = (edges[1] - edges[0]) / 2
    vol = (4 / 3) * np.pi * ((centers + half) ** 3 - (centers - half) ** 3)
    return centers, n / (max(n.sum(), 1) * vol)


def plot_simulation(res, axis_lim: float = 5.0, every: int | None = None):
    """Four-panel summary like simulation_plotter.m. Returns the matplotlib figure."""
    import matplotlib.pyplot as plt

    p = res.params
    fig = plt.figure(figsize=(14, 7))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    plot_capsule(ax, p["Rb"], p["bac_length"])
    ax.scatter(*res.receptors.T, s=4, c="k")
    ax.scatter(*res.phage.pos.T, s=20, c="r", marker="D")
    ax.set(xlim=(-axis_lim, axis_lim), ylim=(-axis_lim, axis_lim), zlim=(-axis_lim, axis_lim),
           xlabel="x", ylabel="y", zlabel="z", title="Simulation")
    ax.view_init(24, 43)

    n = len(res.kernel)
    every = every or max(n // 20, 1)
    t = np.arange(n)[::every]
    a1 = fig.add_subplot(2, 4, 3)
    a1.plot(t, res.kernel[::every], "k")
    a1.set(title="Encounter rate", xlabel="step", ylabel="kernel")
    a2 = fig.add_subplot(2, 4, 4)
    a2.plot(t, res.n_bound[::every], "k")
    a2.set(title="Bound phage", xlabel="step")
    a3 = fig.add_subplot(2, 4, 7)
    a3.plot(t, res.n_attached[::every], "ko")
    a3.set(title="Attached phage", xlabel="step")
    a4 = fig.add_subplot(2, 4, 8)
    a4.plot(t, res.tail_median[::every], "ko")
    a4.set(title="Tail length (median)", xlabel="step", ylabel="um")
    fig.tight_layout()
    return fig


def plot_walk(res, start: int = 10, axis_lim: float = 5.0):
    """Receptors plus the surface trajectory of phage 0 (end of encounter_rates_main.m)."""
    import matplotlib.pyplot as plt

    p = res.params
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    plot_capsule(ax, p["Rb"], p["bac_length"])
    ax.scatter(*res.receptors.T, s=4, c="k")
    w = res.walk[start:]
    ax.plot(*w.T, "m", lw=2)
    ax.scatter(*w[0], c="g", s=40)
    ax.scatter(*w[-1], c="r", s=40)
    ax.set(xlim=(-axis_lim, axis_lim), ylim=(-axis_lim, axis_lim), zlim=(-axis_lim, axis_lim))
    ax.view_init(24, 43)
    return fig
