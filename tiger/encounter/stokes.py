"""Stokes flow around a translating sphere (the swimming bacterium).

Ports the velocity formulas repeated in encounter_rates_2D.m,
encounter_rates_hydro/Langevin_fun.m, surface_walks/Langevin_fun.m and stokes_vel.m.

For a sphere of radius Rb moving with velocity U, the disturbance velocity at
r (relative to the sphere centre) is

  u_s(r) = (3 Rb / 4 r) (1 + Rb^2 / 3 r^2) U
         + (3 Rb / 4 r) (1 / r^2 - Rb^2 / r^4) r (U . r)

and in the sphere's frame a passive tracer moves with U - u_s (MATLAB:
``x - (v_comp + x_comp) dt + U dt``).
"""
from __future__ import annotations

import numpy as np


def stokes_disturbance(r: np.ndarray, U: np.ndarray, Rb: float, mask=None) -> np.ndarray:
    """u_s at points r (N x d) for sphere velocity U (d,). ``mask`` zeros disabled axes."""
    r = np.asarray(r, float)
    U = np.asarray(U, float)
    if mask is not None:
        r = r * mask
    rr = np.linalg.norm(r, axis=-1, keepdims=True)
    rr = np.where(rr == 0, np.inf, rr)
    dot = r @ U
    a = (3 * Rb / (4 * rr)) * (1 + Rb ** 2 / (3 * rr ** 2))
    b = (3 * Rb / (4 * rr)) * (1 / rr ** 2 - Rb ** 2 / rr ** 4)
    return a * U + b * r * dot[..., None]


def comoving_velocity(r: np.ndarray, U: np.ndarray, Rb: float) -> np.ndarray:
    """Tracer velocity in the sphere's frame: U - u_s (the flow past a fixed sphere)."""
    return np.asarray(U, float) - stokes_disturbance(r, U, Rb)


def stokes_field_2d(v_x: float, v_y: float, Rb: float, LX: float, LY: float, h: float = 0.01):
    """Speed field around the sphere on a grid (stokes_vel.m).

    Returns (S, X, Y, UX, UY): speed magnitude and components of
    (v - u_s) with r clipped to >= Rb inside the body, as in the MATLAB.
    """
    x = np.arange(-LX / 2, LX / 2 + h / 2, h)
    y = np.arange(-LY / 2, LY / 2 + h / 2, h)
    X, Y = np.meshgrid(x, y)
    r = np.sqrt(X ** 2 + Y ** 2)
    r[r < Rb] = Rb
    dot = v_x * X + v_y * Y
    a = (3 * Rb / (4 * r)) * (1 + Rb ** 2 / (3 * r ** 2))
    b = (3 * Rb / (4 * r)) * (1 / r ** 2 - Rb ** 2 / r ** 4)
    UX = v_x - (a * v_x + b * X * dot)
    UY = v_y - (a * v_y + b * Y * dot)
    return np.hypot(UX, UY), X, Y, UX, UY


def plot_stokes_field(v_x=30.0, v_y=0.0, Rb=1.0, LX=10.0, LY=10.0, h=0.05):
    """Contour + streamlines + centre-line speed profiles (stokes_vel.m figure)."""
    import matplotlib.pyplot as plt

    S, X, Y, UX, UY = stokes_field_2d(v_x, v_y, Rb, LX, LY, h)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(1, 2, 1)
    ax.contourf(X, Y, S, cmap="gray")
    ax.streamplot(X, Y, -UX, -UY, density=1.2, color="k", linewidth=0.7)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(Rb * np.cos(th), Rb * np.sin(th), "k", lw=3)
    ax.set_aspect("equal")
    ax2 = fig.add_subplot(2, 2, 2)
    j0 = np.argmin(np.abs(X[0]))
    ax2.plot(Y[:, j0][Y[:, j0] > 0] - Rb, S[:, j0][Y[:, j0] > 0], "k")
    ax2.set(xlabel="distance (um)", ylabel="|v_F|")
    ax3 = fig.add_subplot(2, 2, 4)
    i0 = np.argmin(np.abs(Y[:, 0]))
    ax3.plot(X[i0][X[i0] > 0] - Rb, S[i0][X[i0] > 0], "k")
    ax3.set(xlabel="distance (um)", ylabel="|v_F|")
    fig.tight_layout()
    return fig
