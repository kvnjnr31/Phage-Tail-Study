"""Adsorption-rate theory: diffusive and ballistic limits, temperature scaling, steric factors.

Ports (Phage_Bacteria_2): TheoreticalAdsorbtion.m, PhageTailExperiment.m,
DiffusionConstants.m, LinearModeling.m, LinearModelingMoldovan.m,
Data Modeling/phage_characteristics.m.

Most functions work in CGS-hour units like the MATLAB: lengths in cm, D in
cm^2/hr, rates in mL/hr (= cm^3/hr). kBT is written as
1.38e-23 * T * 3600^2 * 1e4 (J -> g cm^2 / hr^2) and viscosity as
eta * 3600 / 100 (Pa s -> g / (cm hr)).
"""
from __future__ import annotations

import numpy as np

KB = 1.38e-23

# water viscosity (Pa s) at the temperatures used throughout the MATLAB
TEMPS_C = np.array([1, 3, 4, 10, 15, 20, 24, 30, 35, 37, 40], float)
ETA_WATER = np.array([1.73, 1.62, 1.57, 1.31, 1.14, 1.00, 0.91, 0.80, 0.72, 0.69, 0.65]) * 1e-3
ETA25 = 0.89e-3

# Moldovan et al. 2007 lambda adsorption constants (mL/s in the source, x3600 -> mL/hr)
MOLDOVAN_K_MALTOSE = np.array([7.4e-13, 6.0e-12, 3.6e-12, 1.3e-11, 3.56e-11, 1.4e-11, 2.2e-11,
                               1.6e-11, 2.7e-11, 2.6e-11, 2.9e-11]) * 3600
MOLDOVAN_K_GLUCOSE = np.array([9.6e-13, 3.0e-12, 2.3e-12, 1.6e-11, 2.0e-11, 1.2e-11, 3.3e-11,
                               1.5e-11, 2.0e-11, 1.6e-11, 1.6e-11]) * 3600

# E. coli used in the MATLAB: 3.8 x 0.8 um rod, effective radius sqrt(a b / 4)
ECOLI_A_CM = 0.5 * 3.8e-4
ECOLI_B_CM = 0.5 * 0.8e-4
ECOLI_RB_CM = np.sqrt(ECOLI_A_CM * ECOLI_B_CM / 4)


def kbT_cgs_hr(T_K):
    return KB * np.asarray(T_K) * 3600 ** 2 * 1e4


def eta_cgs_hr(eta_pa_s):
    return np.asarray(eta_pa_s) * 3600 / 1e2


def eta_water(T_C):
    """Interpolated water viscosity (Pa s) from the MATLAB table."""
    return np.interp(T_C, TEMPS_C, ETA_WATER)


def stokes_einstein(R, T_K, eta_pa_s, cgs_hr=True):
    """D = kBT / (6 pi eta R). R in cm (cgs_hr) -> cm^2/hr, or SI R in m -> m^2/s."""
    if cgs_hr:
        return kbT_cgs_hr(T_K) / (6 * np.pi * eta_cgs_hr(eta_pa_s) * R)
    return KB * T_K / (6 * np.pi * eta_pa_s * R)


def hydrodynamic_radius(D, T_K, eta_pa_s):
    """Invert Stokes-Einstein (cgs-hr): Rp = kBT / (6 pi eta D)."""
    return kbT_cgs_hr(T_K) / (6 * np.pi * eta_cgs_hr(eta_pa_s) * D)


def scale_D(D_ref, T_K, eta, T_ref=298.15, eta_ref=ETA25):
    """D(T) = D_ref (T / T_ref) (eta_ref / eta)  (DiffusionConstants.m, LinearModeling.m)."""
    return D_ref * (np.asarray(T_K) / T_ref) * (eta_ref / np.asarray(eta))


def k_diffusive(D, R):
    """Smoluchowski: k = 4 pi D R."""
    return 4 * np.pi * np.asarray(D) * np.asarray(R)


def k_ballistic(R, v):
    """k = pi R^2 v."""
    return np.pi * np.asarray(R) ** 2 * np.asarray(v)


def k_run_tumble(R, v, tau):
    """Diffusive limit with bacterial run-and-tumble D_B = v^2 tau / 3 (DiffusionConstants.m).

    The MATLAB plotted 4 * 3 * pi * z * D_B; this returns 4 pi R D_B.
    """
    return 4 * np.pi * np.asarray(R) * np.asarray(v) ** 2 * np.asarray(tau) / 3


def k_prolate(a, b, D):
    """Diffusive rate to a prolate spheroid (TheoreticalAdsorbtion.m): 4 pi a D / ln(2a/b)."""
    return 4 * np.pi * a * D / np.log(2 * a / b)


def capture_radius(Rb, Rp, tail):
    """z = Rb + Rp + tail (all MATLAB scripts)."""
    return Rb + Rp + tail


def sherwood_swim(Rb, v, D):
    """Pe = Rb v / D, Sh = (1 + (1 + 2 Pe)^(1/3)) / 2 (LinearModelingMoldovan.m)."""
    Pe = Rb * np.asarray(v) / D
    return 0.5 * (1 + (1 + 2 * Pe) ** (1 / 3))


# ---------------------------------------------------------------- scripts as functions
def theoretical_adsorption(Dp=6.2e-12, vb=30e-6, L=150e-9, H=60e-9):
    """TheoreticalAdsorbtion.m, lambda on E. coli (SI in, mL/hr out).

    Returns dict with k_diff, k_ball, steric factors psi = (k_exp / k_theory) (z / b)
    for the Moldovan maltose and glucose data, and temperatures.
    """
    a, b = 0.5 * 3.8e-6, 0.5 * 0.8e-6
    b = np.sqrt(a * b / 4)
    z = b + L + H
    kdiff = 4 * np.pi * z * Dp * 1e6 * 3600
    kball = np.pi * z ** 2 * vb * 1e6 * 3600
    out = dict(T_C=TEMPS_C, k_diff=kdiff, k_ball=kball, k_exp_maltose=MOLDOVAN_K_MALTOSE,
               k_exp_glucose=MOLDOVAN_K_GLUCOSE)
    for name, k in (("maltose", MOLDOVAN_K_MALTOSE), ("glucose", MOLDOVAN_K_GLUCOSE)):
        out[f"psi_diff_{name}"] = (k / kdiff) * (z / b)
        out[f"psi_ball_{name}"] = (k / kball) * (z / b)
    return out


def kdiff_size_map(T_C=24.0, Dp_range_cm=(30e-7, 200e-7), tail_range_cm=(5e-7, 1000e-7), n=400):
    """PhageTailExperiment.m figure 2: k_diff over capsid diameter x tail length.

    Returns (capsid_diam_nm, tail_nm, kdiff_mL_hr) grids.
    """
    T = T_C + 273.15
    kbT = kbT_cgs_hr(T)
    eta = eta_cgs_hr(eta_water(T_C))
    Dp = np.linspace(*Dp_range_cm, n)
    rho = np.linspace(*tail_range_cm, n)
    DP, RHO = np.meshgrid(Dp, rho)
    RP = 0.5 * DP
    D = kbT / (6 * np.pi * eta * RP)
    k = 4 * np.pi * D * (ECOLI_RB_CM + RP + RHO)
    return DP * 1e7, RHO * 1e7, k


def kdiff_head_vs_tail(T_C=24.0, Rp=0.39e-5, rho=150e-7, eps=None):
    """PhageTailExperiment.m figure 1: grow the tail (Rp fixed) vs grow the head by eps (cm)."""
    eps = np.linspace(-1e-6, 1e-6, 1000) if eps is None else eps
    T = T_C + 273.15
    kbT, eta = kbT_cgs_hr(T), eta_cgs_hr(eta_water(T_C))
    z = ECOLI_RB_CM + Rp + rho + eps
    k_tail = 4 * np.pi * kbT / (6 * np.pi * eta * Rp) * z
    k_head = 4 * np.pi * kbT / (6 * np.pi * eta * (Rp + eps)) * z
    return eps * 1e7, k_tail, k_head


# phage geometry used in PhageTailExperiment.m (capsid radius, tail length; cm)
PHAGE_GEOMETRY_CM = {
    "Lambda": (0.5 * 63e-7, 150e-7), "T1": (0.39e-5, 150e-7), "T5": (0.5 * 65e-5, 150e-7),
    "T2": (27.5e-7, 290e-7), "T4": (27.5e-7, 290e-7), "T6": (27.5e-7, 290e-7),
    "T7": (0.5 * 60e-7, 290e-7), "T3": (0.5 * 60e-7, 290e-7), "J1": (27.5e-7, 290e-7),
}


def diffusion_constants(maeda_speed: np.ndarray, maeda_tumble: np.ndarray, temps_c=None):
    """DiffusionConstants.m: lambda D(T), E. coli run-and-tumble D(T) and rate limits.

    ``maeda_speed`` and ``maeda_tumble`` are the (T_C, value) arrays from
    data/bacteria_characteristics.xlsx (tiger.analysis.maeda_arrays): speed in um/s, tumble rate in 1/s.
    Returns a dict of arrays (cm^2/hr and mL/hr) on the splined temperatures.
    """
    from scipy.interpolate import PchipInterpolator

    temps_c = np.array([1, 3, 4, 10, 15, 20, 24, 30, 35, 37, 39], float) if temps_c is None else temps_c
    T = TEMPS_C + 273.15
    DLam = 6.2e-8 * 3600
    DLams = scale_D(DLam, T, ETA_WATER)
    Rp = hydrodynamic_radius(DLams, T, ETA_WATER)[0]
    z = ECOLI_RB_CM + Rp + 150e-7
    Tm, vel = maeda_speed[:, 0], maeda_speed[:, 1] * 1e-6
    tumb = maeda_tumble[:, 1]
    vels = PchipInterpolator(Tm, vel)(temps_c)
    tumbs = PchipInterpolator(maeda_tumble[:, 0], tumb)(temps_c)
    vbs = vels * 1e2
    taus = 1 / tumbs
    D_ecoli = (1 / 3) * taus * vbs ** 2 * 3600
    return dict(T_C=temps_c, speed_um_s=vels * 1e6, tumble_rate=tumbs, D_ecoli=D_ecoli,
                T_lambda_C=TEMPS_C, D_lambda=DLams, z=z,
                k_runtumble=4 * 3 * np.pi * z * D_ecoli,  # as plotted in the MATLAB
                k_lambda_diff=4 * np.pi * z * DLams,
                k_ball=np.pi * z ** 2 * vbs * 3600)


def linear_model_moldovan(maeda_speed: np.ndarray, which="maltose", drop=None):
    """LinearModelingMoldovan.m: regress Moldovan k_exp on k_diff(T) (x1e7), dropping the outlier.

    Outliers: index 4 (15 C) for maltose, 6 (24 C) for glucose, as in the MATLAB.
    Returns dict with predictor, response, OLS result (statsmodels if available,
    else numpy) and the Sherwood numbers.
    """
    from scipy.interpolate import PchipInterpolator

    T = TEMPS_C + 273.15
    D = scale_D(6.2e-8 * 3600, T, ETA_WATER)
    Rp = hydrodynamic_radius(D, T, ETA_WATER)
    z = ECOLI_RB_CM + Rp
    kdiff = 4 * np.pi * D * z
    vb = PchipInterpolator(maeda_speed[:, 0], maeda_speed[:, 1] * 1e-6)(TEMPS_C) * 1e2 * 3600
    Sh = sherwood_swim(ECOLI_RB_CM, vb, D)
    k = MOLDOVAN_K_MALTOSE if which == "maltose" else MOLDOVAN_K_GLUCOSE
    drop = (4 if which == "maltose" else 6) if drop is None else drop
    keep = np.ones(len(k), bool); keep[drop] = False
    X, y = kdiff[keep] * 1e7, k[keep] * 1e7
    from .stats import ols
    fit = ols(y, X[:, None], names=["X"])
    return dict(X=X, y=y, outlier=(kdiff[drop] * 1e7, k[drop] * 1e7), fit=fit, k_ball=np.pi * z ** 2 * vb,
                sherwood=Sh, k_diff=kdiff)


def linear_modeling_literature():
    """LinearModeling.m: k_max * Rb/(Rb+z) vs measured k for T2, lambda, T4, T3 at 37 C."""
    corr = (310 / 298) * (ETA25 / 0.69e-3)
    D = np.array([2.8e-8, 6.2e-8, 3.6e-8, 1.19e-7]) * 3600 * corr
    k = np.array([1.26e-7, 9.36e-8, 1.38e-7, 1.8e-7])
    Rb = 0.5e-4
    Rp = kbT_cgs_hr(298) / (6 * np.pi * eta_cgs_hr(ETA25) * D)
    rho = np.array([113, 150, 113, 150]) * 1e-7
    z = Rb + Rp + rho
    kmax = 4 * np.pi * D * z
    X = kmax * Rb / (Rb + z)
    return dict(phage=["T2", "Lambda", "T4", "T3"], X=X, k=k, D=D, Rp=Rp)


def family_diffusion(T=300.0, eta=1e-3):
    """phage_characteristics.m: mean tail/head (nm) and Stokes-Einstein D (m^2/s) per tail family."""
    fam = {
        "Myoviridae": ([113, 170, 135, 120, 150, 170, 160, 53.4], [111, 85, 60, 60, 94, 64, 144.5, 85]),
        "Siphoviridae": ([151, 180, 135, 980, 210, 100, 140, 355], [60, 80, 60, 56, 55, 53, 60, 81]),
        "Podoviridae": ([20, 15, 12, 10, 18, 13, 17, 46], [60, 70, 63, 70, 59, 145, 60, 67, 60, 54]),
    }
    out = {}
    for k, (tails, heads) in fam.items():
        h = np.mean(heads)
        out[k] = dict(tail_nm=np.mean(tails), head_nm=h,
                      D=KB * T / (6 * np.pi * eta * h * 1e-9 / 2))
    return out
