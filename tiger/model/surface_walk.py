"""Phage motion: bulk CTRW jumps and random walks on the bacterium surface.

Ports Langevin_fun.m (surface walk with tether forces) and noise_gen.m
(bulk jumps + attachment) from Phage_Bacteria_2/Simulation.

Noise convention (all MATLAB files): each Cartesian kick is
sqrt(24 D0 dt) * (U - 1/2), U ~ U(0,1), whose variance is 2 D0 dt.
CTRW waiting times: tau = tmin * U^(-1/nu)  (Pareto; nu = 1 gives normal diffusion).
"""
from __future__ import annotations

import numpy as np

from .geometry import axis_point, distance_to_axis, level_set
from .phage import PhageState

KB = 1.38e-23


def stokes_einstein_D(radius_um: float, T: float = 300.0, eta: float = 1e-3) -> float:
    """Diffusion constant in um^2/s for a sphere of radius ``radius_um`` (um).

    Same as encounter_rates_main.m: D0 = 1e12 * kBT / (6 pi eta R).
    """
    return 1e12 * KB * T / (6 * np.pi * eta * radius_um * 1e-6)


def _uniform_kick(rng, n, D0, dt):
    return np.sqrt(24 * D0 * dt) * (rng.random((n, 3)) - 0.5)


def surface_step(ph: PhageState, MatB: np.ndarray, Rb: float, D0: float, k_bind: float,
                 L: float, dt: float, rng: np.random.Generator,
                 force_scale: float = 0.0, exact_projection: bool = True) -> PhageState:
    """One step of the surface random walk for attached phage (Langevin_fun.m).

    1. Build a random tangent vector at each attached phage: v = n x w, w ~ U(0,1)^3,
       scaled to length sqrt(24 D0 dt) (U - 1/2) and rotated about the normal by a
       random angle (Rodrigues rotation).
    2. Bound phage feel a linear tether force from their receptor,
       F = k_bind * x * e,  x = |r_rec - r| - tail - Rp,  only when x > 0,
       projected onto the tangent plane.
    3. The kicked point is projected back onto the offset capsule of radius
       Rb + (tail + Rp)/2.

    ``force_scale`` multiplies the tether force before it is added as a
    displacement. The MATLAB line is ``r1 = r0 + v + 0*F_proj`` (tether force
    switched off), so the default is 0.0. Set force_scale > 0 to turn it on.

    Fixes relative to the MATLAB:
    * the Heaviside used ``x/abs(x)`` (a matrix division for vectors); here it is
      elementwise ``x > 0``;
    * the Rodrigues rotation used the un-normalised gradient as the axis; here the
      unit normal is used so the kick length is not rescaled by |grad phi|.
    """
    idx = np.flatnonzero(ph.attached)
    if idx.size == 0:
        return ph
    r0 = ph.pos[idx]
    N = idx.size
    _, grad = level_set(r0, Rb, L)
    n = grad / np.linalg.norm(grad, axis=1, keepdims=True)

    # random tangent kick
    w = rng.random((N, 3))
    v = np.cross(n, w)
    vn = np.linalg.norm(v, axis=1, keepdims=True)
    vn[vn == 0] = 1.0
    length = np.sqrt(24 * D0 * dt) * (rng.random((N, 1)) - 0.5)
    v = length * v / vn
    th = 2 * np.pi * rng.random((N, 1))
    v = v * np.cos(th) + np.cross(n, v) * np.sin(th) + n * np.sum(n * v, 1, keepdims=True) * (1 - np.cos(th))

    # tether force from bound receptor
    F_proj = np.zeros((N, 3))
    rec = ph.receptor[idx]
    b = rec >= 0
    if np.any(b) and force_scale != 0.0:
        vb = MatB[rec[b]] - r0[b]
        vmag = np.linalg.norm(vb, axis=1, keepdims=True)
        vmag[vmag == 0] = 1.0
        x = vmag[:, 0] - ph.tail_length[idx][b] - ph.head_radius[idx][b]
        F = (k_bind * np.clip(x, 0, None))[:, None] * vb / vmag
        nb = n[b]
        F_proj[b] = F - np.sum(F * nb, 1, keepdims=True) * nb

    r1 = r0 + v + force_scale * F_proj

    # project back to offset surface radius Rb + (tail + Rp)/2
    R_off = Rb + (ph.tail_length[idx] + ph.head_radius[idx]) / 2
    if exact_projection:
        c = axis_point(r1, L)
        d = r1 - c
        dn = np.linalg.norm(d, axis=1, keepdims=True)
        dn[dn == 0] = 1.0
        r2 = c + R_off[:, None] * d / dn
    else:
        phi, g = level_set(r1, R_off, L)
        r2 = r1 - (phi / np.sum(g * g, 1))[:, None] * g

    ph.prev[idx] = 0.0  # MATLAB sets the "drift velocity" columns to zero here
    ph.pos[idx] = r2
    return ph


def bulk_jumps(ph: PhageState, D0: float, dt: float, tmin: float, nu: float, time: float,
               Rb: float, L: float, rng: np.random.Generator, max_tries: int = 100000) -> PhageState:
    """CTRW jumps for detached phage, with attachment on contact (noise_gen.m).

    Detached phage whose waiting time has elapsed draw a new waiting time
    tau = tmin * U^(-1/nu) and make one uniform jump. Jumps that would put the
    phage within Rb + Rp of the bacterium axis are redrawn; if any redraw was
    needed (the phage "hit" the surface) the phage becomes attached.
    """
    det = np.flatnonzero(~ph.attached)
    ready = det[time - ph.t_wait[det] >= 0]
    if ready.size == 0:
        return ph
    U1 = rng.random(ready.size)
    ph.t_wait[ready] = time + tmin * (1.0 / U1) ** (1.0 / nu)

    pending = ready.copy()
    new = np.empty((ready.size, 3))
    tries = np.zeros(ready.size, dtype=int)
    slot = np.arange(ready.size)
    for _ in range(max_tries):
        cand = ph.pos[pending] + _uniform_kick(rng, pending.size, D0, dt)
        Rmin = Rb + ph.head_radius[pending]
        ok = distance_to_axis(cand, L) > Rmin
        tries[slot] += 1
        new[slot[ok]] = cand[ok]
        pending, slot = pending[~ok], slot[~ok]
        if pending.size == 0:
            break
    else:
        raise RuntimeError("bulk_jumps: could not find a non-overlapping jump")

    ph.prev[ready] = ph.pos[ready]
    ph.pos[ready] = new
    ph.attached[ready[tries > 1]] = True
    return ph
