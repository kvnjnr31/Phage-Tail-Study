"""Tail-receptor binding, box boundaries and overlap checks.

Ports binding_fun.m, boundarys_receptor.m and overlap.m (Phage_Bacteria_2/Simulation).
"""
from __future__ import annotations

import warnings

import numpy as np

from .geometry import distance_to_axis
from .phage import PhageState


def binding_fun(ph: PhageState, MatB: np.ndarray, Rr: float, rng: np.random.Generator,
                rule: str = "matlab"):
    """Attempt receptor binding for every unbound phage (binding_fun.m).

    For phage k, receptors within the tail-cone radius Rcone = tail + Rp are found
    and their covered fraction is P = n_in_cone * (Rr / Rcone)^2.

    ``rule="matlab"`` (default) reproduces the MATLAB test ``rand > P and P > 0``:
    binding happens with probability 1 - P whenever at least one receptor is in
    the cone.
    ``rule="slides"`` follows the update slides ("probability of NOT interacting
    is the receptor-free area fraction; bind if rand > P_free"), i.e. bind with
    probability P (the covered fraction).

    The bound receptor is chosen uniformly among those in the cone.

    Returns (Nenc, Nact, ph) where Nenc is the number of new bindings and Nact
    is the number of phage with the attached flag set (MATLAB counts
    ``MatP(:,9) > 0``; its comment calls these "active (unbound)").
    """
    Nact = int(np.count_nonzero(ph.attached))
    Nenc = 0
    for k in np.flatnonzero(~ph.bound):
        Rcone = ph.tail_length[k] + ph.head_radius[k]
        d = np.linalg.norm(MatB - ph.pos[k], axis=1)
        inside = np.flatnonzero(d < Rcone)
        if inside.size == 0:
            continue
        P = inside.size * (Rr / Rcone) ** 2
        u = rng.random()
        hit = (u > P) if rule == "matlab" else (u > 1.0 - P)
        if hit:
            ph.receptor[k] = inside[rng.integers(inside.size)]
            Nenc += 1
    return Nenc, Nact, ph


def boundarys_receptor(ph: PhageState, x_right: float, y_right: float, z_right: float,
                       Rp: float, rng: np.random.Generator, tail_mutation: float = 0.01):
    """Box boundary handling with tail-length mutation (boundarys_receptor.m).

    A phage leaving through a face is re-inserted just inside the opposite face
    (coordinate -> -(edge - Rp) * sign). Unbound phage that cross get a tail
    length mutation e ~ U(-tail_mutation, tail_mutation) clipped to [0, 1].
    Head radius mutation is multiplied by 0 in the MATLAB code, so heads are
    only clipped ([0, 1] on x faces, [0.05, 0.5] on y/z faces), as there.
    With ``tail_mutation=0`` traits are left untouched (no clipping), which is
    what tail-length sweeps need.

    Fixes relative to the MATLAB:
    * the loop tested ``MatP(kk,8)`` (row kk) instead of ``MatP(jj(kk),8)``;
    * the "reactivate" line set ``MatP(kk,8) = 1``, which would bind phage kk to
      receptor 1. Here a phage re-entering the box is simply detached and unbound.
    """
    for axis, edge, (hlo, hhi) in ((0, x_right, (0.0, 1.0)), (1, y_right, (0.05, 0.5)),
                                   (2, z_right, (0.05, 0.5))):
        jj = np.flatnonzero(np.abs(ph.pos[:, axis]) >= edge)
        if jj.size == 0:
            continue
        ph.pos[jj, axis] = -(edge - Rp) * np.sign(ph.pos[jj, axis])
        ph.prev[jj] = 0.0
        ph.prev[jj, axis] = 1.0
        ub = jj[~ph.bound[jj]]
        if ub.size and tail_mutation == 0.0:
            ph.attached[ub] = False
        elif ub.size:
            e = 2 * tail_mutation * rng.random(ub.size) - tail_mutation
            ph.tail_length[ub] = np.clip(ph.tail_length[ub] + e, 0.0, 1.0)
            ph.head_radius[ub] = np.clip(ph.head_radius[ub], hlo, hhi)
            ph.attached[ub] = False
    return ph


def overlap(ph: PhageState, Rb: float, L: float, warn: bool = True) -> np.ndarray:
    """Indices of phage whose centre is inside the bacterium (overlap.m).

    The MATLAB version popped up a 3D figure for each offender; this returns the
    indices and (optionally) warns.
    """
    bad = np.flatnonzero(distance_to_axis(ph.pos, L) < Rb)
    if warn and bad.size:
        warnings.warn(f"{bad.size} phage overlap the bacterium", RuntimeWarning, stacklevel=2)
    return bad
