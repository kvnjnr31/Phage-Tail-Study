"""Plot 2D slices of a TIGER flow field saved by ``run_tiger.py``.

(A) x-y plane through the body centre (swimming axis horizontal): colour is |u|
    including the out-of-plane component, white lines are in-plane streamlines.
(B) y-z cross-section through the body centre: the swirl around the cell.
Also prints |u| on spherical shells around the body and the mean axial flow,
which shows whether the cell is actually pushing fluid along its axis.

Run from models/tiger_lbm:

    python simulation/plot_flow_slice.py                         # data/numerical/run, time-averaged field
    python simulation/plot_flow_slice.py --field final
    python simulation/plot_flow_slice.py --run data/numerical/thrust --out figures/thrust_slice.png
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
a = ap.parse_args()
out = a.out or os.path.join(ROOT, "figures", f"tiger_flow_slice_{a.field}.png")
os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

U = np.load(os.path.join(a.run, f"u_{a.field}.npy")).astype(float)
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

fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw=dict(width_ratios=[1.5, 1]))
S = U[:, :, iz, :]
ax = axs[0]
im = ax.pcolormesh(x, y, np.linalg.norm(S, axis=-1).T, shading="auto", cmap="viridis")
ax.streamplot(x, y, S[:, :, 0].T, S[:, :, 1].T, color="w", density=1.4, linewidth=0.6, arrowsize=0.8)
ax.plot(*ring, "r-", lw=1.5, label=f"body (R$_B$ = {Rb * dx_um:g} µm)")
ax.set(xlabel="x (µm), swimming axis", ylabel="y (µm)", aspect="equal", title="(A) x–y plane through the body centre")
ax.legend(loc="lower left", fontsize=8, framealpha=0.7)
fig.colorbar(im, ax=ax, label="|u| incl. out-of-plane (lattice units)")

S = U[ix, :, :, :]
ax = axs[1]
im = ax.pcolormesh(y, z, np.linalg.norm(S, axis=-1).T, shading="auto", cmap="viridis")
ax.streamplot(y, z, S[:, :, 1].T, S[:, :, 2].T, color="w", density=1.2, linewidth=0.6, arrowsize=0.8)
ax.plot(*ring, "r-", lw=1.5)
ax.set(xlabel="y (µm)", ylabel="z (µm)", aspect="equal", title="(B) y–z cross-section through the body")
fig.colorbar(im, ax=ax, label="|u| (lattice units)")
lab = f"time average over the last {int(navg)} steps" if a.field == "mean" else "final step"
fig.suptitle(f"{a.title} — {lab}", fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=200)

# diagnostics: decay with distance and net axial flow
mag = np.linalg.norm(U, axis=-1)
r = np.linalg.norm(np.indices(mag.shape).transpose(1, 2, 3, 0) - c, axis=-1)
print(f"saved {out}")
print(f"max |u| = {mag.max():.3e} lu;  domain-mean u_x = {U[..., 0].mean():.3e} lu;  "
      f"body drift = {np.linalg.norm(path[-1] - path[0]):.2e} cells")
for R in (Rb, 1.5 * Rb, 2 * Rb, 3 * Rb):
    shell = np.abs(r - R) < 1
    print(f"  r = {R * dx_um:.1f} µm: median |u| = {np.median(mag[shell]):.3e} lu")
