"""Background flows that advect free (detached) phage around the bacterium.

A flow is any object with ``velocity(pos, t) -> (N, 3)`` in um/s, expressed in
the bacterium's frame of the receptor model (capsule along y, centred at the
origin). Two are provided:

* :class:`StokesSwimmerFlow` - analytic Stokes flow past a sphere of radius R
  translating with velocity U (the hydro encounter model), in the body frame:
  u = -u_s(r) + ... i.e. fluid streams past the body at -U far away.
* :class:`TigerLBMFlow` - the velocity field computed by the TIGER lattice-
  Boltzmann swimmer (models/tiger_lbm, hybrid_simulation.py), loaded from
  ``velocity_field_small.npy`` / ``velocity_field.npy`` and interpolated
  trilinearly. Lattice units are converted with the calibration in
  hybrid_simulation.py: dx = Rb_um / Rb_lu, dt_unit = f_bundle_lb / f_bundle_Hz.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .encounter.stokes import stokes_disturbance


class NoFlow:
    def velocity(self, pos, t):
        return np.zeros_like(pos)


class StokesSwimmerFlow:
    """Flow past a sphere swimming with velocity ``U`` (um/s), seen from the body frame.

    Far from the body the fluid moves at -U; on the sphere surface it is zero
    (no slip). ``R`` is the effective hydrodynamic radius (um); for the capsule
    of the receptor model an equal-volume sphere is a reasonable choice.
    """

    def __init__(self, U=(0.0, -30.0, 0.0), R=1.2):
        self.U = np.asarray(U, float)
        self.R = float(R)

    def velocity(self, pos, t):
        # lab-frame disturbance of a sphere moving at U, minus U -> body frame.
        # Inside r < R (between the capsule and the effective sphere) the
        # exterior solution is meaningless; the fluid is taken to be at rest there.
        u = stokes_disturbance(pos, self.U, self.R) - self.U
        u[np.linalg.norm(pos, axis=1) < self.R] = 0.0
        return u


def equal_volume_radius(Rb, L):
    """Radius of the sphere with the capsule's volume (cylinder of length 2L + two caps)."""
    V = np.pi * Rb ** 2 * 2 * L + 4 / 3 * np.pi * Rb ** 3
    return (3 * V / (4 * np.pi)) ** (1 / 3)


class TigerLBMFlow:
    """Interpolated TIGER LBM velocity field in the bacterium frame.

    Parameters
    ----------
    field_path : path to velocity_field(_small).npy with shape (T, nx, ny, nz, 3)
                 in lattice units (lu / tu).
    path_path  : optional bacterium_path.npy (T_full, 3); the body centre used to
                 shift the grid into the body frame (default: grid centre of
                 the first frame, from the config: [48, 32, 32]).
    Rb_lu, Rb_um : body radius in lattice and physical units (config: 10, 1.0).
    f_bundle_lb, f_bundle_hz : measured bundle frequency (1/tu) and target (Hz)
                 (config: 1.6978e-2, 100) giving seconds per lattice time unit.
    frame      : which saved frame to use ("mean" averages all frames, an int
                 picks one). The LBM run spans ~ms, far shorter than the phage
                 time step, so a time-averaged (quasi-steady) field is the default.
    axis_map   : how LBM axes map to the receptor-model axes. The LBM swimmer
                 swims along x; the receptor capsule's long axis is y, so the
                 default maps (lbm x, y, z) -> (model y, z, x).
    speed_scale: extra multiplier on the converted velocity (1 = calibrated).
    target_speed: if given (um/s), rescale the field so that the median fluid
                 speed in the shell 1.5-3 um from the body centre equals it.
                 Needed in practice: with the saved run (Oct 2026 check) the
                 calibrated field has a median speed of ~470 um/s and a 99th
                 percentile of ~26,000 um/s inside the grid, far above E. coli's
                 ~30 um/s swim speed, so only the *shape* of the field is usable.

    Outside the LBM grid (about +-4.8 x +-3.2 x +-3.2 um) the velocity is zero.
    """

    def __init__(self, field_path, path_path=None, Rb_lu=10.0, Rb_um=1.0, f_bundle_lb=1.6978e-2,
                 f_bundle_hz=100.0, frame="mean", axis_map=(1, 2, 0), center_lu=(48, 32, 32),
                 speed_scale=1.0, target_speed=None, mmap=True):
        from scipy.interpolate import RegularGridInterpolator

        F = np.load(field_path, mmap_mode="r" if mmap else None)
        if frame == "mean":
            u = np.asarray(F.mean(axis=0))
        else:
            u = np.asarray(F[int(frame)])
        nx, ny, nz, _ = u.shape
        self.dx_um = Rb_um / Rb_lu
        self.dt_unit_s = abs(f_bundle_lb) / f_bundle_hz
        self.vel_scale = self.dx_um / self.dt_unit_s * speed_scale  # (lu/tu) -> um/s
        c = np.asarray(center_lu, float)
        if path_path is not None and Path(path_path).exists():
            c = np.load(path_path)[0]
        self.axis_map = tuple(axis_map)
        # grid coordinates in um relative to body centre (LBM axes)
        self.axes = [(np.arange(n) - c[i]) * self.dx_um for i, n in enumerate((nx, ny, nz))]
        self._interp = RegularGridInterpolator(self.axes, u * self.vel_scale, bounds_error=False,
                                               fill_value=0.0)
        self.extent_um = [(a[0], a[-1]) for a in self.axes]
        self.rescale = 1.0
        if target_speed is not None:
            rng = np.random.default_rng(0)
            g = rng.uniform(-3, 3, size=(20000, 3))
            r = np.linalg.norm(g, axis=1)
            sh = g[(r > 1.5) & (r < 3)]
            med = np.median(np.linalg.norm(self.velocity(sh, 0.0), axis=1))
            if med > 0:
                self.rescale = target_speed / med

    def velocity(self, pos, t):
        pos = np.asarray(pos, float)
        # model -> LBM coordinates: model axis axis_map[i] holds LBM axis i
        q = np.empty_like(pos)
        for i_lbm, i_model in enumerate(self.axis_map):
            q[:, i_lbm] = pos[:, i_model]
        u_lbm = self._interp(q) * getattr(self, "rescale", 1.0)
        out = np.empty_like(u_lbm)
        for i_lbm, i_model in enumerate(self.axis_map):
            out[:, i_model] = u_lbm[:, i_lbm]
        return out

    def describe(self) -> str:
        return (f"TIGER LBM flow: dx = {self.dx_um:.3f} um, 1 lu/tu = {self.vel_scale:.1f} um/s, "
                f"rescale = {self.rescale:.3g}, "
                f"extent (um) = {[tuple(round(float(v), 2) for v in e) for e in self.extent_um]}")
