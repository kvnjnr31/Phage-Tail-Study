"""RATTLE integrator for a phage constrained to the bacterium surface.

Port of RATTLE.m (Phage_Bacteria_2/Simulation): velocity-Verlet with a
holonomic constraint g(x) = |x - c(x)|^2 - R^2 = 0 (stay on the capsule),
driven by Lennard-Jones forces from nearby receptors.

The MATLAB file was a draft that never ran correctly, and the encounter
driver never called it. This version implements the same scheme properly.
Fixes relative to the MATLAB:
* both Newton loops were ``while r_mag < eta`` with r_mag starting at 2, so
  they never executed; here they iterate until the residual < tol;
* the receptor distance used ``(xm(3)*MatB(:,3))`` instead of a difference;
* the force was +dU/dx; the force is -dU/dx;
* ``Rb`` was not passed in, and ``x_lambda``/``pm_plus_half`` were undefined
  when the loop did not run;
* the velocity stage used the old normal ``nm`` instead of ``nm_plus``.

Units: um, s, kg (kBT = 300 * 1.38e-23 * 1e12 kg um^2 / s^2), as in the MATLAB.
"""
from __future__ import annotations

import numpy as np

from .geometry import level_set

KBT_UM = 300 * 1.38e-23 * 1e12


def lj_force(x: np.ndarray, receptors: np.ndarray, sigma: float = 1e-4,
             epsilon: float = 0.5 * KBT_UM, rc: float = 1.0) -> np.ndarray:
    """Total Lennard-Jones force on a point x from receptors within cutoff rc.

    U(r) = 4 eps [(sigma/r)^12 - (sigma/r)^6]; returns F = -sum dU/dx.
    """
    d = x - receptors
    r = np.linalg.norm(d, axis=1)
    m = (r < rc) & (r > 0)
    if not np.any(m):
        return np.zeros(3)
    r, d = r[m], d[m]
    sr6 = (sigma / r) ** 6
    dUdr = 4 * epsilon * (-12 * sr6 * sr6 + 6 * sr6) / r
    return -np.sum((dUdr / r)[:, None] * d, axis=0)


def rattle_step(x: np.ndarray, p: np.ndarray, receptors: np.ndarray, Rb: float, L: float,
                dt: float, mass: float = 0.1e-15, tol: float = 1e-12, max_iter: int = 50,
                force=lj_force):
    """Advance one particle (position x, momentum p) by one RATTLE step.

    Returns (x_new, p_new). The new position satisfies g(x_new) = 0 to ``tol``
    and the new momentum is tangent to the surface.
    """
    x = np.asarray(x, float)
    p = np.asarray(p, float)
    f = force(x, receptors)
    _, n0 = level_set(x, Rb, L)

    # position stage: solve for lambda so that g(x_new) = 0
    base = x + (dt / mass) * (p + 0.5 * dt * f)
    step = -(dt * dt / (2 * mass)) * n0
    lam = 0.0
    for _ in range(max_iter):
        xn = base + lam * step
        g, gn = level_set(xn, Rb, L)
        if abs(g) < tol:
            break
        dg = float(gn @ step)
        if dg == 0:
            break
        lam -= g / dg
    xn = base + lam * step
    p_half = p + 0.5 * dt * (f - lam * n0)

    # momentum stage: choose mu so that grad g(x_new) . p_new = 0
    fn = force(xn, receptors)
    _, n1 = level_set(xn, Rb, L)
    q = p_half + 0.5 * dt * fn
    mu = float(n1 @ q) / (0.5 * dt * float(n1 @ n1))
    p_new = q - 0.5 * dt * mu * n1
    return xn, p_new


def rattle(pos: np.ndarray, mom: np.ndarray, receptors: np.ndarray, Rb: float, L: float,
           dt: float, **kw):
    """Apply :func:`rattle_step` to every row of pos/mom (N x 3). Returns new arrays."""
    X = np.empty_like(pos, dtype=float)
    Pm = np.empty_like(mom, dtype=float)
    for i in range(pos.shape[0]):
        X[i], Pm[i] = rattle_step(pos[i], mom[i], receptors, Rb, L, dt, **kw)
    return X, Pm
