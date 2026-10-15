"""Random-walk noise shared by all encounter simulators.

* ``uniform_kick``: sqrt(24 D0 dt) * (U - 1/2) per axis (variance 2 D0 dt), the
  increment used in every MATLAB Langevin/noise_gen file.
* ``pareto_wait``: tau = tmin * U^(-1/nu)  (noise_gen*.m, diffusion_times.m).
* ``ml_wait``: tau = dt * (Gamma(1 - nu) U)^(-1/nu)  (track_gen.m, phage_msd*.m),
  the heavy-tailed CTRW waiting time used for the subdiffusion MSD runs.
* ``tmin_for_nu``: tmin = dt * sinc(nu)^(1/nu)  (diffusion_times.m), using
  MATLAB's normalised sinc(x) = sin(pi x)/(pi x).
"""
from __future__ import annotations

import numpy as np
from math import gamma


def uniform_kick(rng: np.random.Generator, shape, D0, dt):
    return np.sqrt(24 * np.asarray(D0) * dt) * (rng.random(shape) - 0.5)


def pareto_wait(rng: np.random.Generator, n: int, tmin, nu):
    U = rng.random(n)
    return np.asarray(tmin) * (1.0 / U) ** (1.0 / np.asarray(nu))


def ml_wait(rng: np.random.Generator, n: int, dt: float, nu: float):
    U = rng.random(n)
    return dt * (gamma(1 - nu) * U) ** (-1.0 / nu)


def tmin_for_nu(dt: float, nu: float) -> float:
    return dt * np.sinc(nu) ** (1.0 / nu)
