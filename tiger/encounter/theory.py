"""Analytic encounter kernels and small analytic models.

Ports plot_fun.m (encounter_rates_algorithm / encounter_rates_hydro),
Diffusion Times/plot_fun.m and distro_plotter.m, sherwood_numbers.m, and
adsorption_model.m + my_model.m (Phage_Bacteria_1).
"""
from __future__ import annotations

import numpy as np


def effective_D(t, v0, omega, D0, nu, dt):
    """Run-and-tumble + (sub)diffusive phage diffusivity used in plot_fun.m.

    Dt(t) = 2 (v0^2 / 3 omega) t + 2 D0 t^nu / dt^(nu-1)        (MSD-like term)
    D(t)  = 2 v0^2 / 3 omega + 2 nu D0 t^(nu-1) / dt^(nu-1)     (its derivative)
    """
    t = np.asarray(t, float)
    Dt = 2 * (v0 ** 2 / (3 * omega)) * t + 2 * D0 * t ** nu / dt ** (nu - 1)
    D = 2 * v0 ** 2 / (3 * omega) + 2 * nu * D0 * t ** (nu - 1) / dt ** (nu - 1)
    return Dt, D


def diffusive_kernel(t, Rcrit, v0, omega, D0, nu, dt):
    """Time-dependent Smoluchowski kernel J(t) = 4 pi Rcrit D (1 + Rcrit^2 / sqrt(pi Dt t)).

    Note: the MATLAB uses Rcrit^2 / sqrt(pi Dt t); the classical transient term is
    Rcrit / sqrt(pi D t). Kept as in the MATLAB.
    """
    t = np.asarray(t, float)
    Dt, D = effective_D(t, v0, omega, D0, nu, dt)
    C = 1 + Rcrit ** 2 / np.sqrt(np.pi * Dt * t)
    return 4 * np.pi * Rcrit * D * C


def ballistic_kernel(Rcrit, v):
    """J_ball = pi Rcrit^2 v."""
    return np.pi * Rcrit ** 2 * v


def to_mL_per_min(kernel_um3_per_s):
    """MATLAB plot scaling 6 * kernel / 1e2, labelled mL/min in the figures."""
    return 6 * np.asarray(kernel_um3_per_s) / 1e2


def plot_kernel(kernel, Rcrit, v0, D0, dt, nu, omega, ax=None, scale=1.0):
    """Simulated kernel vs diffusive and ballistic theory (plot_fun.m)."""
    import matplotlib.pyplot as plt

    ax = ax or plt.gca()
    n = len(kernel)
    t = dt * np.arange(1, n + 1)
    ax.plot(t, scale * np.asarray(kernel), "k", lw=2, label="simulation")
    ax.plot(t, scale * diffusive_kernel(t, Rcrit, v0, omega, D0, nu, dt), "b--", lw=2, label="diffusive")
    ax.axhline(scale * ballistic_kernel(Rcrit, v0), color="r", ls="--", lw=2, label="ballistic")
    ax.set(xlabel="time (s)", ylabel="kernel")
    ax.legend()
    return ax


def mean_capture_time(Rcrit, D0, nu, dt):
    """Diffusion Times/plot_fun.m: t = (Rcrit^2 / D)^(1/nu), D = 6 D0 / dt^(nu-1)."""
    D = 6 * D0 / dt ** (nu - 1)
    return (Rcrit ** 2 / D) ** (1 / nu)


def mfpt_1d_ctrw(LX, nu, D0, tmin, dt):
    """Mean first passage out of |x| < LX for the 1D CTRW (mean_time_to_capture.m)."""
    return LX ** 2 / (2 * nu * D0 * (1 - (tmin / (dt + tmin)) ** nu))


def capture_ratio_curve(D0=2.54, R=1.1, t=None):
    """distro_plotter.m: cumulative ratio of two Gaussian-kernel integrals."""
    t = np.arange(0.01, 1000.0, 0.01) if t is None else np.asarray(t)
    g = np.exp(-(R ** 2) / (2 * D0 * t)) / np.sqrt(4 * np.pi * D0 * t)
    return t, np.cumsum(t * g) / np.cumsum(g)


def sherwood(U, L, D):
    """Sherwood number Sh = (1 + (1 + Pe)^(1/3)) / 2 with Pe = U L / D (sherwood_numbers.m)."""
    Pe = np.asarray(U) * L / D
    return 0.5 * (1 + (1 + Pe) ** (1 / 3))


def sherwood_kernel(U, L, D):
    """J = 4 pi D Sh (per unit radius, as plotted in sherwood_numbers.m)."""
    return 4 * np.pi * D * sherwood(U, L, D)


# T4 wild-type and T4 hoc- diffusion coefficients at increasing mucin (sherwood_numbers.m)
D_T4_WT = np.array([3.84, 2.38, 1.18, 1.03, 0.90, 0.42])
D_T4_HOC = np.array([3.64, 3.14, 2.63, 2.54, 1.37, 0.74])


def two_state_adsorption(t_span=(0.01, 600.0), P0=(0.0, 1e4), D0=2.56, tauD=1e-2, alpha=1.0,
                         vB=30.0, Rb=1.1, Ra=1.0, c=5e-5, t_eval=None):
    """Reversible/irreversible adsorption ODE (adsorption_model.m + my_model.m).

    P1 = reversibly bound, P2 = free phage:
      dP1/dt = c [pi vB Rb^2 P2 - 4 pi D(t) Rb P1 - 4 pi D0 Ra P1]
      dP2/dt = c [-pi vB Rb^2 P2 + 4 pi D(t) Rb P1]
    with D(t) = 6 alpha K t^(alpha-1), K = D0 / tauD^(alpha-1).
    Returns scipy's OdeResult (t in seconds).
    """
    from scipy.integrate import solve_ivp

    K = D0 / tauD ** (alpha - 1)

    def rhs(t, P):
        D = 6 * alpha * K * t ** (alpha - 1)
        return c * np.array([np.pi * vB * Rb ** 2 * P[1] - 4 * np.pi * D * Rb * P[0] - 4 * np.pi * D0 * Ra * P[0],
                             -np.pi * vB * Rb ** 2 * P[1] + 4 * np.pi * D * Rb * P[0]])

    return solve_ivp(rhs, t_span, list(P0), method="RK45", t_eval=t_eval, rtol=1e-6)
