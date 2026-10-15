"""Capsule (spherocylinder) geometry for a rod-shaped bacterium.

The bacterium is a cylinder of radius ``Rb`` along the y-axis for |y| <= L
("body"), capped by hemispheres centred at y = +L ("head") and y = -L ("tail"),
where ``L`` is ``bac_length`` in the MATLAB code (half the cylinder length).

Ports the surface pieces that were repeated inline across
Langevin_fun.m, RATTLE.m, noise_gen.m and overlap.m in Phage_Bacteria_2/Simulation:

* region classification (body / head / tail),
* the level-set function phi(r) = |r - c(r)|^2 - R^2 and its gradient 2 (r - c),
* projection of a point back onto an (offset) capsule surface,
* inside/overlap tests and a surface mesh for plotting.
"""
from __future__ import annotations

import numpy as np

BODY, HEAD, TAIL = 0, 1, 2


def region(points: np.ndarray, L: float) -> np.ndarray:
    """Return 0 (body, |y|<=L), 1 (head, y>L) or 2 (tail, y<-L) for each point."""
    y = np.asarray(points)[..., 1]
    out = np.full(y.shape, BODY, dtype=int)
    out[y > L] = HEAD
    out[y < -L] = TAIL
    return out


def axis_point(points: np.ndarray, L: float) -> np.ndarray:
    """Closest point on the capsule's central segment (0, clip(y,-L,L), 0)."""
    p = np.asarray(points, dtype=float)
    c = np.zeros_like(p)
    c[..., 1] = np.clip(p[..., 1], -L, L)
    return c


def level_set(points: np.ndarray, R, L: float):
    """phi = |r - c|^2 - R^2 and grad phi = 2 (r - c), as in the MATLAB code.

    ``R`` may be a scalar or a per-point array (the MATLAB walk uses an offset
    radius Rb + (tail + Rp)/2 that differs per phage).
    """
    p = np.asarray(points, dtype=float)
    d = p - axis_point(p, L)
    phi = np.sum(d * d, axis=-1) - np.asarray(R) ** 2
    return phi, 2.0 * d


def outward_normal(points: np.ndarray, L: float) -> np.ndarray:
    """Unit outward normal of the capsule through each point."""
    d = np.asarray(points, dtype=float) - axis_point(points, L)
    n = np.linalg.norm(d, axis=-1, keepdims=True)
    n[n == 0] = 1.0
    return d / n


def distance_to_axis(points: np.ndarray, L: float) -> np.ndarray:
    """Distance from each point to the central segment; < R means inside capsule R."""
    p = np.asarray(points, dtype=float)
    return np.linalg.norm(p - axis_point(p, L), axis=-1)


def project_to_surface(points: np.ndarray, R, L: float, exact: bool = True) -> np.ndarray:
    """Project points onto the capsule of radius R (scalar or per point).

    ``exact=False`` reproduces the single Newton step used in Langevin_fun.m:
    r2 = r1 - phi * grad / |grad|^2, which lands slightly outside the surface.
    ``exact=True`` (default) places the point exactly at radius R from the axis,
    which is what that Newton step converges to.
    """
    p = np.asarray(points, dtype=float)
    if exact:
        c = axis_point(p, L)
        d = p - c
        n = np.linalg.norm(d, axis=-1, keepdims=True)
        n[n == 0] = 1.0
        return c + np.asarray(R)[..., None] * d / n if np.ndim(R) else c + R * d / n
    phi, g = level_set(p, R, L)
    gg = np.sum(g * g, axis=-1)
    gg[gg == 0] = 1.0
    return p - (phi / gg)[..., None] * g


def inside(points: np.ndarray, R, L: float) -> np.ndarray:
    """True where the point lies strictly inside the capsule of radius R (overlap.m)."""
    return distance_to_axis(points, L) < R


def surface_area(R: float, L: float) -> float:
    """Capsule area: cylinder (2 pi R * 2L) + sphere (4 pi R^2)."""
    return 2 * np.pi * R * (2 * L) + 4 * np.pi * R ** 2


def capsule_mesh(R: float, L: float, n: int = 60):
    """(X, Y, Z) arrays for plotting the capsule with matplotlib plot_surface."""
    u = np.linspace(0, 2 * np.pi, n)
    # body
    yb = np.linspace(-L, L, 2)
    U, Yb = np.meshgrid(u, yb)
    body = (R * np.cos(U), Yb, R * np.sin(U))
    # caps
    v = np.linspace(0, np.pi / 2, n // 2)
    U2, V = np.meshgrid(u, v)
    head = (R * np.cos(V) * np.cos(U2), L + R * np.sin(V), R * np.cos(V) * np.sin(U2))
    tail = (R * np.cos(V) * np.cos(U2), -L - R * np.sin(V), R * np.cos(V) * np.sin(U2))
    return body, head, tail


def plot_capsule(ax, R: float, L: float, color="tab:blue", alpha=0.3, n=60):
    """Draw the capsule on a 3D matplotlib axis (replaces the cylinder/sphere surf calls)."""
    for X, Y, Z in capsule_mesh(R, L, n):
        ax.plot_surface(X, Y, Z, color=color, alpha=alpha, linewidth=0)
    return ax
