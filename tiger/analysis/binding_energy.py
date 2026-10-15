"""Tail-receptor binding potential, rupture and escape (first-passage) times.

Ports Phage_Bacteria_2/Data/first_passage.m and data_analyzer.m.

The gp37-LPS force-distance curve (Brzozowska et al. 2019; data/binding_force_curve_gp37_LPS.csv)
is linear up to the first rupture, F = k_bind x with k_bind ~ -65 pN/nm. The
potential U(x) = -int F dx is truncated where F drops below the rupture force,
mirrored to x < 0, and the mean escape time over the barrier is the
Kramers/first-passage integral

  tau = (1/D) int_{-xA}^{xA} exp(-U(y)/kT) [ int_{-xA}^{y} exp(U(z)/kT) dz ] dy

with D the Stokes-Einstein constant of a 0.1 um head.
"""
from __future__ import annotations

import numpy as np

KBT = 300 * 1.38e-23


def truncated_potential(kbind=-65e-3, f_rup_pN=-79.0, x_max_nm=50.0, n=10_000):
    """Mirrored potential from a linear force truncated at the rupture force.

    kbind in N/m (pN/nm * 1e-3), rupture force in pN. Returns (d, U) in m and J.
    """
    x = np.linspace(0, x_max_nm, n) * 1e-9
    dx = x[1] - x[0]
    F = kbind * x
    F[F < f_rup_pN * 1e-12] = 0.0
    U = np.cumsum(-F[:-1] * dx)
    d = np.concatenate([-x[1:][::-1], x[:-1]])
    UU = np.concatenate([U[::-1], U])
    return d, UU


def first_passage(kbind=-65e-3, f_rup_pN=-79.0, head_radius_m=0.1e-6, eta=1e-3, n=10_000):
    """Mean escape time (s) over the truncated binding potential (first_passage.m)."""
    d, UU = truncated_potential(kbind, f_rup_pN, n=n)
    dx = d[1] - d[0]
    Uesc = UU.max() / KBT
    xA = np.min(np.abs(d[UU / KBT >= Uesc]))
    D = KBT / (6 * np.pi * eta * head_radius_m)
    Ut = UU[np.abs(d) <= xA] / KBT
    I = np.cumsum(np.exp(Ut) * dx)
    return float(np.sum(np.exp(-Ut) * I * dx) / D)


def escape_times(rupture_forces_pN=None, kbind=-65e-3, **kw):
    """first_passage over a range of rupture forces (data_analyzer.m bottom panel)."""
    f = np.linspace(-50, -120, 20) if rupture_forces_pN is None else np.asarray(rupture_forces_pN)
    return f, np.array([first_passage(kbind, fr, **kw) for fr in f])


def potential_from_force_curve(force_data: np.ndarray):
    """Integrate a measured (distance nm, force pN) curve into a mirrored potential (J)."""
    dist = force_data[:, 0] * 1e-9
    force = force_data[:, 1] * 1e-12
    dU = -force[:-1] * np.diff(dist)
    U = np.cumsum(dU)
    d = np.concatenate([-dist[:-1][::-1], dist[:-1]])
    UU = np.concatenate([U[::-1], U])
    return d, UU


def fit_sech_potential(force_data: np.ndarray):
    """Fit U/kT = a (1 - 2 e^{b x} / (e^{2 b x} + 1)) = a (1 - sech(b x)) to the mirrored potential.

    x in nm. Returns (a, b), the data (d_nm, U_kT) and the energy/distance at
    the first rupture (max |F|).
    """
    from scipy.optimize import curve_fit

    d, UU = potential_from_force_curve(force_data)
    x, y = d / 1e-9, UU / KBT
    f = lambda x, a, b: a * (1 - 1 / np.cosh(b * x))
    (a, b), cov = curve_fit(f, x, y, p0=[y.max(), 0.1], maxfev=20000)
    j = np.argmin(force_data[:, 1])
    return dict(a=a, b=b, cov=cov, d_nm=x, U_kT=y,
                U_rupture_kT=UU[len(UU) // 2 + j] / KBT if j < len(UU) // 2 else np.nan,
                d_rupture_nm=force_data[j, 0])


def plot_binding_analysis(force_data: np.ndarray, kbind=-65e-3, f_rup=-79.0):
    """Three-panel figure like data_analyzer.m: force curve, truncated potential, escape times."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 7))
    ax1 = fig.add_subplot(1, 2, 1)
    xs = np.arange(-2, 18)
    ax1.plot(xs, -65 * xs - 55, "k--", lw=2)
    ax1.plot(force_data[:, 0], force_data[:, 1], "o", mfc="0.5", mec="k")
    ax1.set(xlabel="distance (nm)", ylabel="force (pN)")
    ax2 = fig.add_subplot(2, 2, 2)
    d, UU = truncated_potential(kbind, f_rup, n=1000)
    x = np.linspace(0, 50, 1000) * 1e-9
    U0 = np.cumsum(-(kbind * x)[:-1] * (x[1] - x[0]))
    UU0 = np.concatenate([U0[::-1], U0])
    m = np.abs(d) < 2e-9
    ax2.plot(d[m] / 1e-9, UU0[m] / KBT, "k--", label="(1/2) k x^2")
    ax2.plot(d[m] / 1e-9, UU[m] / KBT, "k", label="truncated")
    ax2.set(xlabel="x (nm)", ylabel="U / kT"); ax2.legend()
    ax3 = fig.add_subplot(2, 2, 4)
    f, tau = escape_times(kbind=kbind)
    ax3.plot(f, tau / 60, "ks")
    ax3.set(xlabel="rupture force (pN)", ylabel="mean escape time (min)")
    fig.tight_layout()
    return fig


def schwartz_saturation(schwartz_data: np.ndarray, c=1.25e-2):
    """data_converter.m: adsorption vs receptor number with k = k_max (1 - exp(-c n)).

    ``schwartz_data`` columns: receptors per cell, rate (1e-10 mL/s units as in the
    MATLAB, converted with 3600e-10 to mL/hr). Returns (n, rates, model).
    """
    rates = schwartz_data[:, 1] * 3600e-10
    kmax = rates.max()
    n = np.arange(1, int(schwartz_data[:, 0].max()) + 1)
    return schwartz_data[:, 0], rates, kmax * (1 - np.exp(-n * c)), n
