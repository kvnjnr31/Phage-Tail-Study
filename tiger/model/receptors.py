"""Receptor placement on the capsule surface.

Port of init_receptors.m and receptor_overlap.m (Phage_Bacteria_2/Simulation).

Receptors are scattered inside a helical strip that winds along the bacterium:
the strip has pitch angle ``m`` and half-thickness ``delta``. A receptor at
helical angle ``theta`` sits at height

    y = m * theta + (a + delta) / cos(m) - L - Rb,   a ~ U(-delta, delta)

and is wrapped onto the cylinder (|y| <= L) or onto a hemispherical cap.
Candidates closer than 2*Rr to an existing receptor are rejected and redrawn.
"""
from __future__ import annotations

import numpy as np


def _helix_point(rng, delta, m, Rb, L):
    theta_max = 2 * (L + Rb - delta / np.cos(m)) / m
    a = 2 * delta * rng.random() - delta
    theta = theta_max * rng.random()
    y = m * theta + (a + delta) / np.cos(m) - L - Rb
    if abs(y) < L:  # body receptor
        return np.array([Rb * np.cos(theta), y, Rb * np.sin(theta)])
    yc = y - L if y > L else y + L  # head (y > L) or tail cap
    rho = np.sqrt(max(Rb * Rb - yc * yc, 0.0))
    return np.array([rho * np.cos(theta), y, rho * np.sin(theta)])


def receptor_overlap(MatB: np.ndarray, k: int, Rr: float) -> bool:
    """True if receptor k lies within 2*Rr of any earlier receptor (receptor_overlap.m)."""
    if k == 0:
        return False
    d = np.linalg.norm(MatB[:k] - MatB[k], axis=1)
    return bool(np.any(d < 2 * Rr))


def init_receptors(B: int, delta: float, m: float, Rb: float, L: float, Rr: float,
                   rng: np.random.Generator | None = None, max_tries: int = 100000) -> np.ndarray:
    """Place B non-overlapping receptors in a helical strip; returns (B, 3) positions.

    Parameters mirror init_receptors.m: B receptors, strip half-thickness delta,
    pitch angle m (radians), bacterium radius Rb, half cylinder length L,
    receptor radius Rr.
    """
    rng = np.random.default_rng() if rng is None else rng
    MatB = np.zeros((B, 3))
    for k in range(B):
        for _ in range(max_tries):
            MatB[k] = _helix_point(rng, delta, m, Rb, L)
            if not receptor_overlap(MatB, k, Rr):
                break
        else:
            raise RuntimeError(f"could not place receptor {k} without overlap; "
                               "reduce B or Rr, or widen the strip")
    return MatB


def receptor_density(MatB: np.ndarray, Rb: float, L: float) -> float:
    """Receptors per unit surface area of the capsule."""
    from .geometry import surface_area
    return len(MatB) / surface_area(Rb, L)
