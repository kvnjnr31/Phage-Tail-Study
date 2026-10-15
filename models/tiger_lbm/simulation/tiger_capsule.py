"""TIGER swimmer with a solid capsule body and a thrusting flagellar bundle (D3Q19 LBM).

The cell is the capsule used by the phage model (``tiger.model``): radius R and
cylinder half-length L (default 0.8 and 1.9 um, total length 5.4 um). It is a
no-slip solid that rotates about its long axis (moving-wall half-way bounce-back).
Behind it a helical bundle pushes on the fluid with

  * thrust  F_t (along -x, pushing fluid backward; the cell swims toward +x), and
  * torque  Gamma about +x (the bundle spin).

The simulation runs in the frame of the cell in a periodic box. In steady state the
wall drag on the fluid balances the thrust, so the cell plus bundle is force-free,
and the far fluid streams past the cell at -U: this is the body-frame flow of a
free swimmer with speed U. The body rotation rate is set from torque balance,
omega_body = -Gamma / C_rot with C_rot = mu (8 pi R^2 L + 8 pi R^3) (cylinder + two
hemispherical caps), so the swimmer is also (approximately) torque-free.

At Re << 1 the flow is linear in the forcing, so fields are reported in lattice units
and in units of the measured swim speed U (``u/U``); multiply by the target swim
speed (e.g. 30 um/s) to get physical velocities. Lengths: dx = config ``dx_um``.

Run from models/tiger_lbm:

    python simulation/tiger_capsule.py                          # config: simulation/tiger_capsule_config.yaml
    python simulation/tiger_capsule.py --steps 400              # quick check
    python simulation/plot_flow_slice.py --run data/numerical/capsule_thrust

Outputs (``--out``, default data/numerical/capsule_thrust): u_final.npy, u_mean.npy (time
average over the last ``avg`` steps), solid.npy (capsule mask), path.npy, meta.npy,
geometry.npy [R, L, x_centre, y_centre, z_centre] (cells), history.csv (U vs step).
"""
import argparse
import os
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

C = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
              [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 0, 1], [-1, 0, -1],
              [1, 0, -1], [-1, 0, 1], [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1]])
W = np.array([1 / 3] + [1 / 18] * 6 + [1 / 36] * 12)
OPP = np.array([np.flatnonzero((C == -ci).all(1))[0] for ci in C])
CS2 = 1 / 3


def capsule_mask(shape, centre, R, L):
    """Solid nodes of a capsule with axis along x: distance to the segment |x - xc| <= L is <= R."""
    X, Y, Z = np.indices(shape, dtype=float)
    ax = np.clip(X - centre[0], -L, L)
    d2 = (X - centre[0] - ax) ** 2 + (Y - centre[1]) ** 2 + (Z - centre[2]) ** 2
    return d2 <= R * R


def deposit(F, pts, vecs, sigma):
    """Spread point forces onto the periodic grid with a normalised Gaussian kernel (+-2 sigma)."""
    shape = F.shape[:3]
    r = int(np.ceil(2 * sigma))
    offs = np.arange(-r, r + 1)
    for p, v in zip(pts, vecs):
        i0 = np.round(p).astype(int)
        gx, gy, gz = (i0[k] + offs for k in range(3))
        wx, wy, wz = (np.exp(-0.5 * (g - p[k]) ** 2 / sigma ** 2) for k, g in enumerate((gx, gy, gz)))
        wgt = wx[:, None, None] * wy[None, :, None] * wz[None, None, :]
        wgt /= wgt.sum()
        ix = np.ix_(gx % shape[0], gy % shape[1], gz % shape[2])
        for k in range(3):
            F[..., k][ix] += wgt * v[k]


def bundle_points(cfg, centre, t):
    """Helix points behind the rear pole and their thrust + torque forces."""
    b = cfg["bundle"]
    R, L = cfg["body"]["R"], cfg["body"]["L"]
    x0 = centre[0] - L - R - b["gap"]
    xs = x0 - np.linspace(0, b["length"], b["n_points"])
    k = 2 * np.pi / b["pitch"]
    th = k * (xs - x0) + b["omega"] * t
    r = b["radius"]
    pts = np.column_stack([xs, centre[1] + r * np.cos(th), centre[2] + r * np.sin(th)])
    n = len(xs)
    thrust = np.array([-b["thrust"] / n, 0.0, 0.0])
    f_tan = b["torque"] / (n * r)                       # tangential force giving torque Gamma about +x
    vecs = np.column_stack([np.full(n, thrust[0]), -f_tan * np.sin(th), f_tan * np.cos(th)])
    return pts, vecs


def run(cfg, steps, avg, out, log_every=200):
    g, body, fl = cfg["grid"], cfg["body"], cfg["fluid"]
    shape = (g["nx"], g["ny"], g["nz"])
    tau = fl["tau"]
    nu = CS2 * (tau - 0.5)
    mu = nu  # rho = 1
    centre = np.array(body.get("centre") or [shape[0] * body.get("centre_x_frac", 0.6), shape[1] / 2, shape[2] / 2])
    R, L = body["R"], body["L"]
    solid = capsule_mask(shape, centre, R, L)
    fluid = ~solid
    C_rot = mu * (8 * np.pi * R ** 2 * L + 8 * np.pi * R ** 3)
    omega_body = -cfg["bundle"]["torque"] / C_rot          # rotation about x (rad / step)

    # links fluid(x) -> solid(x + c_i) and the moving-wall correction at the link midpoint
    X = np.indices(shape, dtype=float)
    links = []
    for i in range(1, 19):
        hit = fluid & np.roll(solid, -C[i], axis=(0, 1, 2))
        if not hit.any():
            continue
        idx = np.nonzero(hit)
        ym = X[1][idx] + 0.5 * C[i][1] - centre[1]
        zm = X[2][idx] + 0.5 * C[i][2] - centre[2]
        uw = np.stack([np.zeros_like(ym), -omega_body * zm, omega_body * ym], 1)   # omega x r, omega along x
        corr = 2 * W[i] * (uw @ C[i]) / CS2
        links.append((i, idx, corr))

    # start from rest, or from a uniform stream -U0 past the cell. The uniform (mean-flow) mode is
    # the slowest to relax (time ~ fluid mass x U / thrust, thousands of steps), so starting near
    # the expected swim speed cuts the run several-fold without changing the steady state.
    U0 = float(cfg["run"].get("init_U", 0.0))
    u0 = np.zeros(shape + (3,))
    u0[..., 0] = -U0
    u0[solid] = 0.0
    cu0 = u0 @ C.T
    f = W * (1 + cu0 / CS2 + 0.5 * cu0 ** 2 / CS2 ** 2 - 0.5 * (u0 * u0).sum(-1)[..., None] / CS2)
    pref = 1 - 1 / (2 * tau)
    mean = np.zeros(shape + (3,))
    n_avg, hist, t0 = 0, [], time.time()
    for step in range(steps):
        F = np.zeros(shape + (3,))
        pts, vecs = bundle_points(cfg, centre, step)
        deposit(F, pts, vecs, cfg["bundle"]["sigma"])
        F[solid] = 0.0
        rho = f.sum(-1)
        u = (f @ C + 0.5 * F) / rho[..., None]
        u[solid] = 0.0
        cu = u @ C.T
        uu = (u * u).sum(-1)[..., None]
        feq = W * rho[..., None] * (1 + cu / CS2 + 0.5 * cu ** 2 / CS2 ** 2 - 0.5 * uu / CS2)
        cF = F @ C.T
        S = pref * W * ((cF - (u * F).sum(-1)[..., None]) / CS2 + cu * cF / CS2 ** 2)
        fpost = f - (f - feq) / tau + S
        for i in range(19):
            f[..., i] = np.roll(fpost[..., i], C[i], axis=(0, 1, 2))
        for i, idx, corr in links:
            f[..., OPP[i]][idx] = fpost[..., i][idx] - rho[idx] * corr
        f[solid] = W                                     # keep solid nodes at rest (not used)
        if step >= steps - avg:
            mean += u
            n_avg += 1
        if step % log_every == 0 or step == steps - 1:
            U = -u[..., 0][fluid].mean()
            hist.append((step, U, np.abs(u).max()))
            msg = ""
            if len(hist) >= 3:   # Aitken extrapolation of the mean-flow relaxation
                (_, a1, _), (_, a2, _), (_, a3, _) = hist[-3:]
                den = a1 + a3 - 2 * a2
                if abs(den) > 1e-15:
                    msg = f"  (extrapolated U = {(a1 * a3 - a2 * a2) / den:.4e})"
            print(f"[{step:05d}] U = {U:.4e}  max|u| = {np.abs(u).max():.3e}  {time.time() - t0:.0f} s{msg}", flush=True)

    os.makedirs(out, exist_ok=True)
    um = mean / max(n_avg, 1)
    np.save(os.path.join(out, "u_final.npy"), u.astype(np.float32))
    np.save(os.path.join(out, "u_mean.npy"), um.astype(np.float32))
    np.save(os.path.join(out, "solid.npy"), solid)
    np.save(os.path.join(out, "path.npy"), np.tile(centre, (steps, 1)))
    np.save(os.path.join(out, "geometry.npy"), np.array([R, L, *centre]))
    dx_m = cfg["units"]["dx_um"] * 1e-6
    np.save(os.path.join(out, "meta.npy"), np.array([dx_m, np.nan, 1.0, R, n_avg]))
    np.savetxt(os.path.join(out, "history.csv"), np.array(hist), delimiter=",", header="step,U,max_abs_u", comments="")
    U = -um[..., 0][fluid].mean()
    Re = U * 2 * R / nu
    print(f"saved {out}\nswim speed U = {U:.4e} lu/step (Re = {Re:.3f}); body omega = {omega_body:.3e} rad/step; "
          f"surface speed / U = {abs(omega_body) * R / U:.2f}")
    return um, U


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "tiger_capsule_config.yaml"))
    ap.add_argument("--steps", type=int, default=None)
    ap.add_argument("--avg", type=int, default=None)
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "numerical", "capsule_thrust"))
    a = ap.parse_args()
    with open(a.config) as fh:
        cfg = yaml.safe_load(fh)
    steps = a.steps or cfg["run"]["steps"]
    avg = min(a.avg if a.avg is not None else cfg["run"]["avg"], steps)
    run(cfg, steps, avg, a.out)
