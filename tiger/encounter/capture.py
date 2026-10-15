"""Escape / capture times for diffusing phage.

Ports Phage_Bacteria_1/Encounter Rates/Diffusion Times:
diffusion_times.m + capture_fun.m (3D release from the bacterium centre) and
mean_time_to_capture.m (1D CTRW escape from |x| < LX).
"""
from __future__ import annotations

import numpy as np

from .noise import pareto_wait, tmin_for_nu, uniform_kick
from .theory import mean_capture_time, mfpt_1d_ctrw


def diffusion_times(P=50_000, nsteps=1000, Rb=1.0, Rp=0.1, dt=0.01, D0=2.54, nu=1.0, seed=None):
    """Release P phage at the origin; record when each first reaches Rcrit = Rb + Rp.

    Every phage jumps every step (the MATLAB waiting-time test is multiplied by 0).
    The recorded time is ``Time - dt`` at detection, as in capture_fun.m.
    Returns dict with per-phage capture times (NaN if never captured), the running
    mean and the theory value (Rcrit^2 / D)^(1/nu), D = 6 D0 / dt^(nu-1).
    """
    rng = np.random.default_rng(seed)
    Rcrit = Rb + Rp
    pos = np.zeros((P, 3))
    t_cap = np.full(P, np.nan)
    active = np.ones(P, bool)
    mean_t = []
    time = 0.0
    for _ in range(nsteps - 1):
        if not active.any():
            break
        pos += uniform_kick(rng, (P, 3), D0, dt)
        out = active & (np.linalg.norm(pos, axis=1) >= Rcrit)
        t_cap[out] = time - dt
        active[out] = False
        mean_t.append(np.nanmean(t_cap) if np.any(~active) else np.nan)
        time += dt
    return dict(capture_times=t_cap, running_mean=np.array(mean_t), dt=dt,
                theory=mean_capture_time(Rcrit, D0, nu, dt))


def mean_time_to_capture(P=1000, LX=10.0, dt=0.01, nu=0.82, delta=1.0, max_steps=10_000_000, seed=None):
    """1D CTRW with Pareto waiting times; time to leave |x| < LX (mean_time_to_capture.m).

    D0 = delta^2 / (2 dt), tmin = dt * sinc(nu)^(1/nu). Returns capture times and the
    theoretical MFPT  LX^2 / (2 nu D0 (1 - (tmin/(dt+tmin))^nu)).
    """
    rng = np.random.default_rng(seed)
    D0 = delta ** 2 / (2 * dt)
    tmin = tmin_for_nu(dt, nu)
    x = np.zeros(P)
    t_wait = pareto_wait(rng, P, tmin, nu)
    active = np.ones(P, bool)
    t_cap = np.full(P, np.nan)
    time = 0.0
    for _ in range(max_steps):
        if not active.any():
            break
        kk = np.flatnonzero(time - t_wait >= 0)
        x[kk] += uniform_kick(rng, kk.size, D0, dt)
        t_wait[kk] = time + pareto_wait(rng, kk.size, tmin, nu)
        out = active & (np.abs(x) >= LX)
        t_cap[out] = time
        active[out] = False
        time += dt
    return dict(capture_times=t_cap, mean=np.nanmean(t_cap), theory=mfpt_1d_ctrw(LX, nu, D0, tmin, dt))
