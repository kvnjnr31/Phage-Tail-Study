"""Mean-squared-displacement simulators for phage (CTRW) and bacteria (run-and-tumble).

Ports Phage_Bacteria_1/bacteria_phage_simulators:
track_gen.m, phage_msd.m, phage_msd_ctrw.m, diffuser_msd.m (CTRW phage),
bacteria_msd.m, propeller_msd.m (2D run-and-tumble), FIG2A_LEFT_PLOT.m (mucin MSD
figure) and reflecting.m. The MATLAB versions stored every track in a J x 2T
matrix and looped in parfor; these are vectorised over trials.
"""
from __future__ import annotations

from math import gamma

import numpy as np

KB = 1.38e-23


def phage_D0(Rp=60e-9, T=300.0, eta=1e-3):
    """Stokes-Einstein D in m^2/s (phage_msd.m defaults: Rp = 60 nm)."""
    return KB * T / (6 * np.pi * eta * Rp)


def ctrw_tracks(nu, D0, dt, steps, J, x0=0.0, y0=0.0, offset_dt=True, seed=None):
    """2D CTRW tracks with waiting times dt (Gamma(1-nu) U)^(-1/nu) (track_gen.m / phage_msd_ctrw.m).

    Returns array (J, steps, 2). ``offset_dt=True`` subtracts dt from each waiting
    time (phage_msd.m, phage_msd_ctrw.m, diffuser_msd.m); track_gen.m did not.
    For nu -> 1, Gamma(1-nu) diverges and every step jumps.
    """
    rng = np.random.default_rng(seed)
    xy = np.zeros((J, 2)); xy[:, 0] = x0; xy[:, 1] = y0
    wt = np.zeros(J)
    out = np.empty((J, steps, 2))
    g = gamma(1 - nu) if nu < 1 else np.inf
    s = np.sqrt(24 * D0 * dt)
    for n in range(steps):
        t = n * dt
        jump = t >= wt
        k = np.count_nonzero(jump)
        if k:
            xy[jump] += s * (rng.random((k, 2)) - 0.5)
            if np.isfinite(g):
                tau = dt * (g * rng.random(k)) ** (-1.0 / nu)
            else:
                tau = np.zeros(k)
            wt[jump] = t + tau - (dt if offset_dt else 0.0)
        out[:, n] = xy
    return out


def ensemble_msd(tracks, dt, origin=None):
    """Ensemble MSD vs lag, measured from the origin (phage_*) or from the first point.

    phage_msd.m measures |r(h) - 0|^2; bacteria_msd.m measures |r(h) - r(1)|^2.
    Returns (tau, msd) for lags 1..T-1.
    """
    ref = tracks[:, :1, :] if origin is None else np.asarray(origin)
    d = tracks[:, 1:, :] - ref
    msd = np.mean(np.sum(d * d, axis=2), axis=0)
    tau = dt * np.arange(1, tracks.shape[1])
    return tau, msd


def fit_power_law(tau, msd, idx=None):
    """Fit log MSD = log K + alpha log tau. Returns (alpha, K)."""
    idx = slice(None) if idx is None else idx
    p = np.polyfit(np.log(tau[idx]), np.log(msd[idx]), 1)
    return p[0], np.exp(p[1])


def phage_msd_ctrw(nu, dt=0.0435, steps=10_000, J=1000, Rp=60e-9, T=300.0, eta=1e-3, seed=None):
    """phage_msd_ctrw.m: returns (tau, msd) in m^2 for exponent nu."""
    D0 = phage_D0(Rp, T, eta)
    tr = ctrw_tracks(nu, D0, dt, steps, J, seed=seed)
    return ensemble_msd(tr, dt, origin=np.zeros(2))


def phage_msd_sweep(nus=None, dt=0.0435, steps=10_000, J=1000, D0=None, sample=slice(4, None, 10), seed=None):
    """phage_msd.m / diffuser_msd.m: fitted diffusion exponent alpha for each CTRW nu.

    Returns dict(nu, alpha, D_fit, theory_prefactor) where the theory MSD is
    4 D0 / (Gamma(nu+1) dt^(nu-1)) * t^nu.
    """
    nus = np.linspace(1e-6, 2, 20) if nus is None else np.asarray(nus)
    D0 = phage_D0() if D0 is None else D0
    alphas, Ds = [], []
    for i, nu in enumerate(nus):
        if nu >= 1:  # Gamma(1-nu) undefined/negative: MATLAB produced complex times
            alphas.append(np.nan); Ds.append(np.nan); continue
        tau, msd = ensemble_msd(ctrw_tracks(nu, D0, dt, steps, J, seed=None if seed is None else seed + i),
                                dt, origin=np.zeros(2))
        a, K = fit_power_law(tau, msd, sample)
        alphas.append(a); Ds.append(K / 4)
    pref = 4 * D0 / (np.array([gamma(n + 1) for n in nus]) * dt ** (nus - 1))
    return dict(nu=nus, alpha=np.array(alphas), D_fit=np.array(Ds), theory_prefactor=pref)


def run_tumble_tracks(rate, v, dt, steps, N, mode="propeller", seed=None):
    """2D run-and-tumble tracks (propeller_msd.m / bacteria_msd.m).

    mode="propeller": constant time step dt, new uniform direction at exponential
    times with rate ``rate`` (propeller_msd.m).
    mode="jump": each step is a full run of exponential duration (bacteria_msd.m,
    which redrew the direction every step and moved v * run_time).
    Returns (N, steps, 2).
    """
    rng = np.random.default_rng(seed)
    xy = np.zeros((N, 2))
    out = np.empty((N, steps, 2))
    th = 2 * np.pi * rng.random(N)
    T = -np.log(rng.random(N)) / rate
    t = 0.0
    for j in range(steps):
        if mode == "propeller":
            tumble = t >= T
            k = np.count_nonzero(tumble)
            th[tumble] = 2 * np.pi * rng.random(k)
            T[tumble] = t - np.log(rng.random(k)) / rate
            step = v * dt
        else:
            th = 2 * np.pi * rng.random(N)
            step = v * (-np.log(rng.random(N)) / rate)
        xy[:, 0] += step * np.cos(th)
        xy[:, 1] += step * np.sin(th)
        out[:, j] = xy
        t += dt
    return out


def run_tumble_D(v, rate, dims=2):
    """Long-time diffusivity of a run-and-tumble walker with uniform reorientation.

    MSD -> 2 dims D t with D = v^2 / (dims * rate). propeller_msd.m compared
    against 2 v^2 / rate (= MSD slope in 2D, i.e. 4D).
    """
    return v * v / (dims * rate)


# mucin MSD parameters (FIG2A_LEFT_PLOT.m): MSD(t) = 4 D_alpha t^alpha
MUCIN_PERCENT = np.array([1e-6, 0.002, 0.006, 0.01, 0.02, 0.04])
MUCIN_D_ALPHA = np.array([15.34, 9.5003, 4.7175, 4.1145, 3.5975, 1.6845]) * 1e-12 / 4
MUCIN_ALPHA = np.array([1.0237, 0.92924, 0.82187, 0.81797, 0.91178, 0.85635])


def fig2a_mucin_msd(idx_use=(0, 1, 3, 5), steps=10_000, trials=1000, dt=0.0435, t_ref=1.0, seed=1):
    """FIG2A_LEFT_PLOT.m: experimental power laws vs CTRW simulations scaled to match at t_ref.

    Returns a matplotlib figure (MSD in um^2 vs t in s, log-log).
    """
    import matplotlib.pyplot as plt

    t_line = np.logspace(-1, 1, 100)
    t_marks = np.logspace(-1, 1, 12)
    fig, ax = plt.subplots(figsize=(6, 6))
    markers = "osd^"
    for k, i in enumerate(idx_use):
        a, Da = MUCIN_ALPHA[i], MUCIN_D_ALPHA[i]
        ax.loglog(t_line, 4 * Da * t_line ** a * 1e12, "k--", lw=2)
        nu = min(a, 0.999)
        tau, msd = phage_msd_ctrw(nu, dt, steps, trials, seed=seed + k)
        f = lambda t: np.exp(np.interp(np.log(t), np.log(tau), np.log(msd)))
        scale = 4 * Da * t_ref ** a / f(t_ref)
        ax.loglog(t_marks, scale * f(t_marks) * 1e12, markers[k % 4], mfc="w", mec="k", ms=8,
                  label=f"{MUCIN_PERCENT[i]:g} % mucin")
    ax.set(xlabel="t (s)", ylabel="MSD (um^2)")
    ax.legend(frameon=False)
    return fig


def reflect(z_new, z_old, z_upper):
    """reflecting.m: reflect a step that crosses |z| = z_upper. Returns (z_new, z_old)."""
    if abs(z_new) >= z_upper:
        z_new = z_old - (z_new - z_old)
        z_old = z_upper * np.sign(z_new)
    else:
        z_old = z_new
    return z_new, z_old
