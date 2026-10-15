#!/usr/bin/env python3
"""
Berg-style recovery of rotation rate and effective drag coefficients
from a rear-view velocity slice (back_slice.npy).

Input file:
    data/numerical/back_slice.npy

Expected shape:
    (T, ny, nz, 3)
        T  : number of time steps
        ny : grid points in y
        nz : grid points in z
        3  : velocity components (ux, uy, uz)

What this script estimates:
    - Mean tangential speed |v_t|
    - Rotation frequency f (1/time)
    - Angular speed omega (rad/time)
    - Effective perpendicular drag coefficient ξ_perp
    - Effective parallel drag coefficient ξ_parallel
      (assuming ξ_perp ≈ 2 ξ_parallel as in Berg-style scaling)

Usage example (from bacteria_ABM root):

    python recovery_script.py --dt 0.2 --r 8.0 --dr 2.0

You can adjust r and dr based on where the swirl is strongest.
"""

import os
import argparse
import numpy as np


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def load_back_slice(path: str) -> np.ndarray:
    arr = np.load(path)
    if arr.ndim != 4 or arr.shape[-1] != 3:
        raise ValueError(f"Expected shape (T, ny, nz, 3), got {arr.shape}")
    return arr


def estimate_center(back_slice: np.ndarray) -> tuple:
    """
    Estimate the swirl center in the y–z plane by locating the
    time-averaged kinetic-energy maximum.

    back_slice: (T, ny, nz, 3)
    returns (cy, cz) as floats
    """
    speed2 = np.sum(back_slice**2, axis=-1)   # (T, ny, nz)
    mean_ke = speed2.mean(axis=0)            # (ny, nz)
    iy, iz = np.unravel_index(np.argmax(mean_ke), mean_ke.shape)
    return float(iy), float(iz)


def sample_ring_indices(ny, nz, center, r, tol):
    """
    Get indices in a y–z ring of radius r ± tol around 'center'.

    center: (cy, cz) in index space
    r     : radius
    tol   : half-thickness of ring
    """
    cy, cz = center
    ys = np.arange(ny)
    zs = np.arange(nz)
    Y, Z = np.meshgrid(ys, zs, indexing="ij")

    dy = Y - cy
    dz = Z - cz
    rad = np.sqrt(dy**2 + dz**2)

    mask = (rad >= r - tol) & (rad <= r + tol)
    iy_list, iz_list = np.where(mask)

    indices = [(int(iy), int(iz)) for iy, iz in zip(iy_list, iz_list)]
    rad_vals = rad[mask]
    return indices, rad_vals


def extract_tangential_velocity(back_slice, indices, center):
    """
    Compute mean tangential velocity around x̂ at the ring indices,
    for each time step.

    back_slice: (T, ny, nz, 3)
    indices   : list of (iy, iz)
    center    : (cy, cz)
    returns v_t_mean: (T,)
    """
    T, ny, nz, _ = back_slice.shape
    cy, cz = center
    v_t_mean = np.zeros(T, dtype=float)

    for t in range(T):
        vt_list = []

        for iy, iz in indices:
            dy = iy - cy
            dz = iz - cz
            r = np.sqrt(dy**2 + dz**2) + 1e-12

            # radial and tangential unit vectors in y–z
            e_r = np.array([0.0, dy/r, dz/r])
            e_phi = np.array([0.0, -dz/r, dy/r])  # tangential around x̂

            u = back_slice[t, iy, iz, :]   # (ux, uy, uz)
            vt = np.dot(u, e_phi)
            vt_list.append(vt)

        v_t_mean[t] = float(np.mean(vt_list)) if vt_list else 0.0

    return v_t_mean


def estimate_rotation_from_vt(v_t, r):
    """
    Berg-style circular motion:
        C = 2π r,  v ≈ mean |v_t|
        f = v / (2π r),  omega = v / r
    """
    v_mag = np.abs(v_t)
    v_mean = float(np.mean(v_mag))
    if r <= 0:
        return 0.0, 0.0, v_mean
    f = v_mean / (2.0 * np.pi * r)
    omega = v_mean / r
    return f, omega, v_mean


def estimate_drag_coefficients(v_t, r, dt, rho=1.0, dr=1.0):
    """
    Crude effective perpendicular drag coefficient via:

        volume of annulus per unit length in x:
            V ≈ 2π r dr
        m ≈ ρ V
        F ≈ m dv/dt
        ξ_perp ≈ F / v

    This is not exact Stokes theory, but provides a consistent
    "Berg-style" effective coefficient from the simulation.
    """
    volume = 2.0 * np.pi * r * dr
    m = rho * volume

    dv = np.diff(v_t)
    dv_dt = dv / dt
    v_mid = 0.5 * (v_t[1:] + v_t[:-1])

    mask = np.abs(v_mid) > 1e-8
    if not np.any(mask):
        return 0.0, 0.0

    F_est = m * dv_dt[mask]
    xi_perp_vals = F_est / v_mid[mask]

    xi_perp = float(np.median(xi_perp_vals))
    xi_parallel = 0.5 * xi_perp   # Berg-style scaling
    return xi_perp, xi_parallel


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Berg-style recovery from back_slice.npy (rear-view velocity slice)."
    )
    parser.add_argument(
        "--path",
        default=os.path.join("data", "numerical", "back_slice.npy"),
        help="Path to back_slice.npy (default: data/numerical/back_slice.npy)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        required=True,
        help="Time step used in the simulation (simulation units).",
    )
    parser.add_argument(
        "--r",
        type=float,
        required=True,
        help="Ring radius in grid units (choose where swirl is strongest).",
    )
    parser.add_argument(
        "--dr",
        type=float,
        default=1.0,
        help="Ring half-thickness (default 1.0).",
    )
    parser.add_argument(
        "--rho",
        type=float,
        default=1.0,
        help="Fluid density (default 1.0 in LBM units).",
    )
    args = parser.parse_args()

    back_slice = load_back_slice(args.path)
    T, ny, nz, _ = back_slice.shape
    print(f"Loaded back_slice: T={T}, ny={ny}, nz={nz}")

    # 1) Estimate swirl center from mean kinetic energy
    center = estimate_center(back_slice)
    print(f"Estimated center (cy, cz) ≈ ({center[0]:.2f}, {center[1]:.2f})")

    # 2) Sample ring at radius r ± dr/2
    indices, _ = sample_ring_indices(ny, nz, center, args.r, tol=args.dr / 2.0)
    if not indices:
        raise RuntimeError("No grid points found in the requested ring; adjust r/dr.")

    print(f"Ring points found: {len(indices)} at radius ≈ {args.r}")

    # 3) Extract tangential velocity over time
    v_t = extract_tangential_velocity(back_slice, indices, center)

    # 4) Berg-style rotation estimates
    f, omega, v_mean = estimate_rotation_from_vt(v_t, args.r)

    # 5) Effective drag coefficients
    xi_perp, xi_parallel = estimate_drag_coefficients(
        v_t, args.r, args.dt, rho=args.rho, dr=args.dr
    )

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------
    print("\n=== Berg-style recovery results ===")
    print(f"Mean |v_t|          : {v_mean:.4e} [grid units / time]")
    print(f"Rotation frequency f: {f:.4e} [1/time]")
    print(f"Angular speed ω     : {omega:.4e} [rad/time]")
    print(f"ξ_perp (effective)  : {xi_perp:.4e} [force / (velocity)]")
    print(f"ξ_parallel (eff)    : {xi_parallel:.4e} [force / (velocity)]")


if __name__ == "__main__":
    main()
