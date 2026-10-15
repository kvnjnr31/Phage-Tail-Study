"""Phage-bacterium encounter-rate simulators (point bacterium, spherical capture).

Ports from Phage_Bacteria_1:
* encounter_rates_3D.m (+ Langevin_fun_3D, noise_gen_3D, encounter_fun_3D,
  boundarys_3D) - "Code for paper/encounter_rates_algorithm" and the copy in
  "Encounter Rates/encounter_rates_3D" (identical apart from a rescaling of v0/omega)
* encounter_rates_ballistic.m
* encounter_rates_2D.m (+ boundarys, encounter_fun) - 2D with Stokes flow
* encounter_rates_stratified (encounter_rates.m, noise_gen.m, Langevin_fun.m,
  encounter_fun.m, boundarys.m, diffusion_vs_depth.m) - mucus layer with
  depth-dependent T4 diffusion and a run-and-tumble bacterium

All return a dict with the per-step encounter kernel
``kernel[n] = (Nenc / Nact) * V / dt`` (um^3/s), where Nenc counts active phage
that came within Rcrit this step and Nact the active phage before the step.
"""
from __future__ import annotations

import numpy as np

from .noise import pareto_wait, uniform_kick
from .stokes import stokes_disturbance


def _init_outside(rng, P, box, Rp, center, Rmin, dim=3):
    lo, hi = -box / 2 + Rp, box / 2 - Rp
    out = np.empty((P, dim))
    k = 0
    while k < P:
        c = rng.uniform(lo, hi, size=(max(2 * (P - k), 1000), dim))
        c = c[np.linalg.norm(c - center[:dim], axis=1) >= Rmin][: P - k]
        out[k:k + len(c)] = c
        k += len(c)
    return out


# ----------------------------------------------------------------------------- 3D
def boundarys_3d(pos, edge, Rad, rng):
    """boundarys_3D.m: periodic wrap on the crossing axis, randomise the other two.

    Randomised coordinates are drawn as 2 (edge - Rad) U - edge, as in the MATLAB.
    """
    for ax in range(3):
        jj = np.flatnonzero(np.abs(pos[:, ax]) >= edge)
        if jj.size == 0:
            continue
        pos[jj, ax] -= 2 * edge * np.sign(pos[jj, ax])
        for o in range(3):
            if o != ax:
                pos[jj, o] = 2 * (edge - Rad) * rng.random(jj.size) - edge
    return pos


def encounter_rates_3d(P=10_000, nsteps=1000, D0=2.54, nu=1.0, tmin=0.001, dt=0.01, v0=30.0,
                       omega=100.0, Rb=1.0, Rp=0.1, LX=10.0, diffusion_scale=0.0,
                       react_dist=5.0, seed=None):
    """Ballistic bacterium through diffusing phage in a periodic box (encounter_rates_3D.m).

    The bacterium moves at (-v0, 0, 0) and never tumbles (omega is only used by
    the theory curve). Phage that come within Rcrit = Rb + Rp are counted and
    deactivated, and reactivated once farther than ``react_dist``.

    ``diffusion_scale`` multiplies the phage noise. In noise_gen_3D.m the noise is
    ``0*sqrt(24*D0*dt)*U`` (phage diffusion switched off), so the default is 0.
    Set 1.0 to turn phage diffusion on.
    """
    rng = np.random.default_rng(seed)
    half = LX / 2
    V = LX ** 3
    Rcrit = Rb + Rp
    bac = np.zeros(3)
    vel = np.array([-v0, 0.0, 0.0])
    pos = _init_outside(rng, P, LX, Rp, bac, Rb + Rp)
    t_wait = np.zeros(P)
    active = np.ones(P, bool)
    kernel = np.zeros(nsteps - 1)
    time = 0.0
    for n in range(nsteps - 1):
        kk = np.flatnonzero(time - t_wait >= 0)
        w = np.zeros((P, 3))
        w[kk] = diffusion_scale * uniform_kick(rng, (kk.size, 3), D0, dt)
        t_wait[kk] = time + pareto_wait(rng, kk.size, tmin, nu)
        pos += w
        bac = bac + vel * dt
        Nact = active.sum()
        dist = np.linalg.norm(pos - bac, axis=1)
        hit = dist <= Rcrit
        Nenc = np.count_nonzero(active & hit)
        active[hit] = False
        active[dist > react_dist] = True
        kernel[n] = Nenc / Nact * V / dt if Nact else np.nan
        boundarys_3d(pos, half, Rp, rng)
        bac = boundarys_3d(bac[None, :], half, Rb, rng)[0]
        time += dt
    return dict(kernel=kernel, dt=dt, Rcrit=Rcrit, v0=v0, D0=D0, nu=nu, omega=omega)


def encounter_rates_ballistic(P=50_000, nsteps=100, v=30.0, Rb=1.0, Rp=0.1, dt=0.01, LX=50.0,
                              seed=None):
    """Pure advection of phage past a stationary bacterium (encounter_rates_ballistic.m).

    Expected kernel: pi (Rb + Rp)^2 v. Note the MATLAB called boundarys_3D with the
    wrong argument list; here only phage boundaries are applied.
    """
    rng = np.random.default_rng(seed)
    half = LX / 2
    V = LX ** 3
    Rcrit = Rb + Rp
    pos = _init_outside(rng, P, LX, Rp, np.zeros(3), Rcrit)
    active = np.ones(P, bool)
    vel = np.array([-v, 0.0, 0.0])
    kernel = np.zeros(nsteps - 1)
    for n in range(nsteps - 1):
        pos += vel * dt
        Nact = active.sum()
        dist = np.linalg.norm(pos, axis=1)
        hit = dist <= Rcrit
        Nenc = np.count_nonzero(active & hit)
        active[hit] = False
        active[dist > 5] = True
        kernel[n] = Nenc / Nact * V / dt if Nact else np.nan
        boundarys_3d(pos, half, Rp, rng)
    return dict(kernel=kernel, dt=dt, Rcrit=Rcrit, theory=np.pi * Rcrit ** 2 * v)


# ----------------------------------------------------------------------------- 2D
def encounter_rates_2d(P=1000, nsteps=1000, v=0.0, Rb=1.0, Rp=0.1, dt=0.01, LX=10.0,
                       D0=2.54, nu=0.82, tmin=0.001, seed=None):
    """2D phage diffusion in Stokes flow past a disc (encounter_rates_2D.m).

    Every phage jumps every step (the MATLAB waiting-time test is multiplied by 0);
    jumps landing inside Rb + Rp are redrawn. Phage are then advected by the
    Stokes flow for a sphere moving at (-v, 0). Contact deactivates a phage; it is
    reactivated only when it leaves the box and is re-inserted at x = +LX/2.
    V = LX^3 as in the MATLAB.
    """
    rng = np.random.default_rng(seed)
    half = LX / 2
    V = LX ** 3
    Rcrit = Rb + Rp
    U = np.array([-v, 0.0])
    pos = _init_outside(rng, P, LX, Rp, np.zeros(3), Rcrit, dim=2)
    active = np.ones(P, bool)
    kernel = np.zeros(nsteps - 1)
    for n in range(nsteps - 1):
        pending = np.arange(P)
        new = pos.copy()
        while pending.size:
            cand = pos[pending] + uniform_kick(rng, (pending.size, 2), D0, dt)
            ok = np.linalg.norm(cand, axis=1) >= Rcrit
            new[pending[ok]] = cand[ok]
            pending = pending[~ok]
        pos = new
        pos = pos - stokes_disturbance(pos, U, Rb) * dt + U * dt
        Nact = active.sum()
        hit = np.linalg.norm(pos, axis=1) <= Rcrit
        Nenc = np.count_nonzero(active & hit)
        active[hit] = False
        kernel[n] = Nenc / Nact * V / dt if Nact else np.nan
        out = np.flatnonzero((np.abs(pos[:, 0]) >= half) | (np.abs(pos[:, 1]) >= half))
        pos[out, 0] = half
        pos[out, 1] = 2 * (half - 0.1) * rng.random(out.size) - half
        active[out] = True
    return dict(kernel=kernel, dt=dt, Rcrit=Rcrit, theory_ballistic=np.pi * Rcrit ** 2 * v)


# ---------------------------------------------------------------------- stratified
_D_MUCIN = np.array([3.64, 3.14, 2.63, 2.54, 1.37, 0.74])
_ALPHA_WT = np.array([1.0237, 0.92924, 0.82187, 0.81797, 0.91178, 0.85635])
_TMIN_WT = np.array([-0.0003, 0.0006, 0.0015, 0.0015, 0.0008, 0.0012])
_ALPHA_HOC = np.array([0.99, 0.99, 1.0, 1.0, 0.99, 0.89])
_MUCIN_EDGES = np.array([0, 0.2, 0.6, 1, 2, 4, 5])


def diffusion_vs_depth(mucin: np.ndarray, ptype: int = 0):
    """(D0, nu, tmin) for T4 wild type (ptype=1) or T4 hoc- (ptype=0) vs mucin % (diffusion_vs_depth.m).

    Mucin outside [0, 5) gives zeros, as in the MATLAB.
    """
    mucin = np.asarray(mucin, float)
    k = np.digitize(mucin, _MUCIN_EDGES) - 1
    ok = (k >= 0) & (k < 6)
    D = np.zeros_like(mucin); a = np.zeros_like(mucin); tm = np.zeros_like(mucin)
    alpha = _ALPHA_WT if ptype == 1 else _ALPHA_HOC
    tmin0 = _TMIN_WT if ptype == 1 else np.zeros(6)
    D[ok], a[ok], tm[ok] = _D_MUCIN[k[ok]], alpha[k[ok]], tmin0[k[ok]]
    return D, a, tm


def encounter_rates_stratified(P=100_000, nsteps=10_000, ptype=0, omega=1.0, dt=0.01, v0=30.0,
                               Rb=1.0, Rp=0.1, LX=10.0, thickness=None, seed=None):
    """Run-and-tumble bacterium in a mucus layer with depth-dependent phage diffusion.

    Assembles the encounter_rates_stratified/ functions into a working driver
    (the MATLAB encounter_rates.m called mismatched helper signatures and could
    not run). The mucin concentration at height z is 3 (z_top - z) / thickness;
    each phage's (D0, nu, tmin) comes from :func:`diffusion_vs_depth`. The
    bacterium starts at the top wall (z = LX/2), picks a new random 3D direction
    at exponential times (rate omega), and reflects off the top wall.
    ``thickness`` defaults to LX (the MATLAB variable was never defined).
    """
    rng = np.random.default_rng(seed)
    half = LX / 2
    V = LX ** 3
    thickness = LX if thickness is None else thickness
    Rcrit = Rb + Rp
    bac = np.array([0.0, 0.0, half])
    bac_prev = bac.copy()
    vel = np.array([-v0, 0.0, 0.0])
    t_tumble = 0.0
    pos = _init_outside(rng, P, LX, Rp, bac, Rcrit)
    prev = pos.copy()
    t_wait = np.zeros(P)
    active = np.ones(P, bool)
    kernel = np.zeros(nsteps - 1)
    time = 0.0
    for n in range(nsteps - 1):
        kk = np.flatnonzero(time - t_wait >= 0)
        mucin = (3 / thickness) * (half - pos[kk, 2])
        D0, nu, tmin = diffusion_vs_depth(mucin, ptype)
        w = np.zeros((P, 3))
        w[kk] = np.sqrt(24 * D0 * dt)[:, None] * (rng.random((kk.size, 3)) - 0.5)
        with np.errstate(divide="ignore", invalid="ignore"):
            tau = np.where(nu > 0, tmin * (1 / rng.random(kk.size)) ** (1 / np.where(nu > 0, nu, 1)), dt)
        t_wait[kk] = time + tau
        if time >= t_tumble:
            th, al = 2 * np.pi * rng.random(2)
            t_tumble = time - np.log(rng.random()) / omega
            vel = v0 * np.array([np.cos(al) * np.cos(th), np.cos(al) * np.sin(th), np.sin(al)])
        prev[:] = pos
        pos += w
        bac_prev = bac.copy()
        bac = bac + vel * dt

        Nact = active.sum()
        dist = np.linalg.norm(pos - bac, axis=1)
        hit = dist <= Rcrit
        Nenc = np.count_nonzero(active & hit)
        active[hit] = False
        active[dist >= half] = True
        kernel[n] = Nenc / Nact * V / dt if Nact else np.nan

        for ax in (0, 1):  # periodic in x, y
            jj = np.abs(pos[:, ax]) >= half
            pos[jj, ax] -= 2 * half * np.sign(pos[jj, ax]); active[jj] = True
            if abs(bac[ax]) >= half:
                bac[ax] -= 2 * half * np.sign(bac[ax])
        jj = np.abs(pos[:, 2]) >= half  # walls in z: undo the step
        pos[jj, 2] = prev[jj, 2]; active[jj] = True
        if abs(bac[2]) > half:
            bac = bac_prev.copy(); vel[2] = -vel[2]
        time += dt
    return dict(kernel=kernel, dt=dt, Rcrit=Rcrit)
