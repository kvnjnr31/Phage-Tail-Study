"""Plot 2D slices of a TIGER flow field saved by ``run_tiger.py``.

(A) x-y plane through the body centre (swimming axis horizontal): colour is |u|
    including the out-of-plane component, white lines are in-plane streamlines.
(B) y-z cross-section through the body centre: the swirl around the cell.
Also prints |u| on spherical shells around the body and the mean axial flow,
which shows whether the cell is actually pushing fluid along its axis.

Run from models/tiger_lbm:

    python simulation/plot_flow_slice.py                         # data/numerical/run, time-averaged field
    python simulation/plot_flow_slice.py --field final
    python simulation/plot_flow_slice.py --run data/numerical/capsule_thrust --frame lab

Runs from ``tiger_capsule.py`` (they contain geometry.npy) are drawn with the capsule outline,
the solid masked out, and velocities in units of the swim speed U. ``--frame body`` shows the
flow seen by the cell (far field = -U); ``--frame lab`` adds U and shows the disturbance the
swimming cell makes in still water.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--run", default=os.path.join(ROOT, "data", "numerical", "run"), help="folder written by run_tiger.py")
ap.add_argument("--field", choices=["mean", "final"], default="mean")
ap.add_argument("--out", default=None, help="output PNG (default: figures/tiger_flow_slice_<field>.png)")
ap.add_argument("--title", default="TIGER lattice Boltzmann flow")
ap.add_argument("--frame", choices=["body", "lab"], default="body", help="capsule runs only")
a = ap.parse_args()
out = a.out or os.path.join(ROOT, "figures", f"tiger_flow_slice_{a.field}.png")
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

U = np.load(os.path.join(a.run, f"u_{a.field}.npy")).astype(float)
geo_file = os.path.join(a.run, "geometry.npy")
capsule = os.path.exists(geo_file)
unit_label = "lattice units"
if capsule:
    Rc, Lc = np.load(geo_file)[:2]
    solid = np.load(os.path.join(a.run, "solid.npy"))
    Uswim = -U[..., 0][~solid].mean()
    if a.frame == "lab":
        U[..., 0] += Uswim
    U[solid] = np.nan
    U /= Uswim
    unit_label = "u / U$_{swim}$"
path = np.load(os.path.join(a.run, "path.npy"))
dx_m, dt_s, dt_lu, Rb, navg = np.load(os.path.join(a.run, "meta.npy"))
nx, ny, nz, _ = U.shape
c = path[-int(navg):].mean(0) if a.field == "mean" and navg > 0 else path[-1]
ix, iy, iz = (int(round(v)) for v in c)
dx_um = dx_m * 1e6
x = (np.arange(nx) - c[0]) * dx_um
y = (np.arange(ny) - c[1]) * dx_um
z = (np.arange(nz) - c[2]) * dx_um
th = np.linspace(0, 2 * np.pi, 200)
ring = (Rb * dx_um * np.cos(th), Rb * dx_um * np.sin(th))
if capsule:   # capsule outline in the x-y plane, and its circular cross-section
    ph = np.linspace(-np.pi / 2, np.pi / 2, 60)
    xs = np.r_[Lc + Rc * np.cos(ph), -Lc - Rc * np.cos(ph[::-1]), Lc + Rc * np.cos(ph[0])] * dx_um
    ys = np.r_[Rc * np.sin(ph), Rc * np.sin(ph[::-1]), Rc * np.sin(ph[0])] * dx_um
    outline_xy = (xs, ys)
    ring = (Rc * dx_um * np.cos(th), Rc * dx_um * np.sin(th))
else:
    outline_xy = ring

fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.5, 1]))
S = U[:, :, iz, :]
ax = axs[0]
im = ax.pcolormesh(x, y, np.linalg.norm(S, axis=-1).T, shading="auto", cmap="viridis")
ax.streamplot(x, y, np.nan_to_num(S[:, :, 0].T), np.nan_to_num(S[:, :, 1].T), color="w", density=1.4,
              linewidth=0.6, arrowsize=0.8)
body_lab = (f"capsule (R = {Rc * dx_um:.2f} µm, length {2 * (Lc + Rc) * dx_um:.1f} µm)" if capsule
            else f"body (R$_B$ = {Rb * dx_um:g} µm)")
ax.plot(*outline_xy, "r-", lw=1.5, label=body_lab)
ax.set(xlabel="x (µm), swimming axis", ylabel="y (µm)", aspect="equal", title="(A) x–y plane through the body centre")
ax.legend(loc="lower left", fontsize=8, framealpha=0.7)
fig.colorbar(im, ax=ax, label=f"|u| incl. out-of-plane ({unit_label})")

S = U[ix, :, :, :]
ax = axs[1]
im = ax.pcolormesh(y, z, np.linalg.norm(S, axis=-1).T, shading="auto", cmap="viridis")
ax.streamplot(y, z, np.nan_to_num(S[:, :, 1].T), np.nan_to_num(S[:, :, 2].T), color="w", density=1.2,
              linewidth=0.6, arrowsize=0.8)
ax.plot(*ring, "r-", lw=1.5)
ax.set(xlabel="y (µm)", ylabel="z (µm)", aspect="equal", title="(B) y–z cross-section through the body")
fig.colorbar(im, ax=ax, label=f"|u| ({unit_label})")
lab = f"time average over the last {int(navg)} steps" if a.field == "mean" else "final step"
if capsule:
    lab += f", {a.frame} frame"
fig.suptitle(f"{a.title} — {lab}", fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=200)

# diagnostics: decay with distance and net axial flow
mag = np.linalg.norm(np.nan_to_num(U), axis=-1)
G = np.indices(mag.shape).transpose(1, 2, 3, 0) - c
if capsule:   # distance from the capsule surface
    G[..., 0] -= np.clip(G[..., 0], -Lc, Lc)
    r = np.linalg.norm(G, axis=-1) - Rc
    shells, what = (1.5, 3, 6, 12), "from the surface"
else:
    r = np.linalg.norm(G, axis=-1)
    shells, what = (Rb, 1.5 * Rb, 2 * Rb, 3 * Rb), "from the centre"
print(f"saved {out}")
if capsule:
    print(f"swim speed U = {Uswim:.4e} lu/step; values below in units of U ({a.frame} frame)")
print(f"max |u| = {mag.max():.3e};  mean u_x = {np.nanmean(U[..., 0]):.3e};  "
      f"body drift = {np.linalg.norm(path[-1] - path[0]):.2e} cells")
for R in shells:
    shell = np.abs(r - R) < 0.75
    print(f"  {R * dx_um:.2f} µm {what}: median |u| = {np.median(mag[shell]):.3e}")
