# phage_bulk_surface_exchange.py
# Simulates phage that diffuse in the fluid around a spherical bacterium,
# ADSORB (stick) when they contact the surface, then SURFACE-WALK for a
# random residence time before DESORBING (knock-off) back into the fluid.
#
# Output: CSV of trajectories + state, and a quick figure.

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ============================================================================
#                               C O N F I G
# Move/adjust these only. CLI flags will override these defaults at runtime.
# ============================================================================
CONFIG = {
    # ---- simulation scale/time ----
    "n_phage": 1000,
    "steps": 10000,
    "dt": 0.01,

    # ---- geometry (units of length) ----
    "rb": 1.0,        # bacterium radius
    "rp": 0.05,       # phage center offset from the cell surface
    "r_init": 1.8,    # initial radius as multiple of (Rb+Rp)

    # ---- motion/kinetics ----
    "sigma_bulk": 0.06,   # 3D bulk step std (length units)
    "sigma_surf": 0.06,   # surface step std (radians on unit sphere)
    "k_off": 0.2,         # surface desorption rate (1/s)

    # ---- small epsilons ----
    "eps_contact_push": 1e-6,  # push to surface shell after deep penetration
    "eps_desorb": 3e-3,        # outward nudge on desorption

    # ---- RNG ----
    "seed": 7,

    # ---- IO ----
    "png": "data/image/exchange_demo.png",
    "csv": "data/numerical/phage_surface_walk_trajectories.csv",

    # ---- quick viz styling ----
    "trail_frames": 200,   # last N frames to draw per trajectory
    "figsize": (7.5, 7.0),
    "lw_bulk": 0.8,
    "lw_surf": 1.4,
    "alpha_bulk": 0.5,
    "alpha_surf": 0.9,
    "scatter_s_bulk": 10,
    "scatter_s_surf": 16,
    "view_elev": 20,
    "view_azim": 30,
}
# ============================================================================

# ---------- helpers ----------
BULK, SURF = 0, 1

def to_unit(v):
    return v / np.clip(np.linalg.norm(v, axis=-1, keepdims=True), 1e-12, None)

def random_unit_vectors(n):
    v = np.random.normal(size=(n, 3))
    return to_unit(v)

def random_tangent_step(u, sigma):
    """
    u: (n,3) unit vectors on sphere (current surface positions' directions)
    returns new unit vectors after a small tangent random step, reprojected to sphere.
    """
    n = u.shape[0]
    t_helper = np.tile(np.array([1.0, 0.0, 0.0]), (n,1))
    collinear = np.abs(np.sum(t_helper*u, axis=1)) > 0.9
    t_helper[collinear] = np.array([0.0, 1.0, 0.0])
    t1 = np.cross(u, t_helper); t1 = to_unit(t1)
    t2 = np.cross(u, t1)
    a = np.random.normal(scale=sigma, size=(n,1))
    b = np.random.normal(scale=sigma, size=(n,1))
    v = u + a*t1 + b*t2
    return to_unit(v)

# ---------- simulation ----------
def simulate_exchange(n_phage, steps, dt, Rb, Rp, sigma_bulk, sigma_surf, k_off,
                      r_init, eps_contact_push=1e-6, eps_desorb=3e-3, seed=None):
    if seed is not None and seed >= 0:
        np.random.seed(seed)

    R = Rb + Rp               # contact radius for phage centers
    prob_off = 1.0 - np.exp(-k_off*dt)

    # initialize in BULK at radius r_init*R with random directions
    u0 = random_unit_vectors(n_phage)
    pos = (r_init * R) * u0
    state = np.full(n_phage, BULK, dtype=int)

    # storage
    traj = np.zeros((steps, n_phage, 3), float)
    states = np.zeros((steps, n_phage), int)
    traj[0] = pos; states[0] = state

    for t in range(1, steps):
        # --- BULK walkers: 3D Gaussian step ---
        bulk_idx = np.where(state == BULK)[0]
        if bulk_idx.size:
            step3d = np.random.normal(scale=sigma_bulk, size=(bulk_idx.size, 3))
            pos[bulk_idx] += step3d

            # If they cross inside the bacterium (|x| < R), reflect to the shell at radius R
            r_now = np.linalg.norm(pos[bulk_idx], axis=1)
            deep = np.where(r_now < R)[0]
            if deep.size:
                u = to_unit(pos[bulk_idx][deep])
                pos[bulk_idx][deep] = (R + eps_contact_push) * u

            # Adsorb upon contact with the surface (first time r <= R)
            r_now = np.linalg.norm(pos[bulk_idx], axis=1)
            hit = np.where(r_now <= R)[0]
            if hit.size:
                u_hit = to_unit(pos[bulk_idx][hit])
                pos[bulk_idx][hit] = R * u_hit
                state[bulk_idx[hit]] = SURF

        # --- SURF walkers: tangent step + possible desorption ---
        surf_idx = np.where(state == SURF)[0]
        if surf_idx.size:
            u = to_unit(pos[surf_idx])                 # directions on the sphere
            u_new = random_tangent_step(u, sigma_surf) # surface move
            pos[surf_idx] = R * u_new

            # Desorb with prob_off
            off_mask = np.random.rand(surf_idx.size) < prob_off
            if np.any(off_mask):
                u_off = to_unit(pos[surf_idx][off_mask])
                pos[surf_idx][off_mask] = (R + eps_desorb) * u_off
                state[surf_idx[off_mask]] = BULK

        traj[t] = pos
        states[t] = state

    time = np.arange(steps) * dt
    return traj, states, time, R

# ---------- IO ----------
def save_csv(traj, states, time, out_csv):
    """
    Write only PHAGE-IDS that EVER walked on the surface (state==1 at least once).
    We still write ALL their timesteps (including their bulk intervals) so downstream
    scripts can choose to filter by state without losing time alignment if needed.
    """
    steps, n, _ = traj.shape
    ever_surf = (states == 1).any(axis=0)  # shape (n,)
    keep_ids = np.where(ever_surf)[0]

    rows = []
    for t in range(steps):
        for i in keep_ids:
            x,y,z = traj[t, i]
            rows.append({
                "t": float(time[t]), "id": int(i),
                "x": float(x), "y": float(y), "z": float(z),
                "state": int(states[t, i])  # 0 bulk, 1 surf
            })
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)

# ---------- quick viz ----------
def quick_figure(traj, states, R, png_path,
                 trail_frames, figsize, lw_bulk, lw_surf,
                 alpha_bulk, alpha_surf, scatter_s_bulk, scatter_s_surf,
                 view_elev, view_azim):
    plt.rcParams.update({
        "figure.figsize": figsize,
        "axes.titlesize": 16,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11
    })

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # draw bacterium surface
    u = np.linspace(0, 2*np.pi, 120)
    v = np.linspace(0, np.pi, 60)
    xs = R * np.outer(np.cos(u), np.sin(v))
    ys = R * np.outer(np.sin(u), np.sin(v))
    zs = R * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, alpha=0.15, linewidth=0, shade=True)

    # plot last N frames to reduce clutter
    T = traj.shape[0]
    t0 = max(0, T - trail_frames)
    end_state = states[-1]
    bulk_ids = np.where(end_state == BULK)[0]
    surf_ids = np.where(end_state == SURF)[0]

    for i in bulk_ids:
        ax.plot(traj[t0:, i, 0], traj[t0:, i, 1], traj[t0:, i, 2],
                lw=lw_bulk, alpha=alpha_bulk)
    for i in surf_ids:
        ax.plot(traj[t0:, i, 0], traj[t0:, i, 1], traj[t0:, i, 2],
                lw=lw_surf, alpha=alpha_surf)

    ax.scatter(traj[-1, bulk_ids, 0], traj[-1, bulk_ids, 1], traj[-1, bulk_ids, 2],
               s=scatter_s_bulk, label="bulk")
    ax.scatter(traj[-1, surf_ids, 0], traj[-1, surf_ids, 1], traj[-1, surf_ids, 2],
               s=scatter_s_surf, label="surface")

    lim = 2.2*R
    ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("Bulk ⇄ Surface Exchange")
    ax.legend(loc="upper left")
    ax.view_init(elev=view_elev, azim=view_azim)
    plt.tight_layout()

    if png_path:
        out_png = Path(png_path)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_png, dpi=220)
    else:
        plt.show()
    plt.close(fig)

# ---------- CLI ----------
def parse_args():
    p = argparse.ArgumentParser(description="Phage bulk–surface exchange on a spherical bacterium.")
    p.add_argument("--n_phage", type=int, default=CONFIG["n_phage"])
    p.add_argument("--steps", type=int, default=CONFIG["steps"])
    p.add_argument("--dt", type=float, default=CONFIG["dt"])

    p.add_argument("--rb", type=float, default=CONFIG["rb"], help="Bacterium radius.")
    p.add_argument("--rp", type=float, default=CONFIG["rp"], help="Phage center offset from the cell surface.")
    p.add_argument("--r_init", type=float, default=CONFIG["r_init"], help="Initial radius as multiple of (Rb+Rp).")

    p.add_argument("--sigma_bulk", type=float, default=CONFIG["sigma_bulk"], help="3D bulk step std (units of length).")
    p.add_argument("--sigma_surf", type=float, default=CONFIG["sigma_surf"], help="Surface step std (radians on unit sphere).")
    p.add_argument("--k_off", type=float, default=CONFIG["k_off"], help="Surface desorption rate (1/s).")

    p.add_argument("--eps_contact_push", type=float, default=CONFIG["eps_contact_push"],
                   help="Outward push to the shell when a bulk walker penetrates inside the cell.")
    p.add_argument("--eps_desorb", type=float, default=CONFIG["eps_desorb"],
                   help="Outward nudge size on desorption to avoid instant re-adsorption.")

    p.add_argument("--png", type=str, default=CONFIG["png"])
    p.add_argument("--csv", type=str, default=CONFIG["csv"])
    p.add_argument("--seed", type=int, default=CONFIG["seed"])

    # quick viz styling
    p.add_argument("--trail_frames", type=int, default=CONFIG["trail_frames"])
    p.add_argument("--fig_w", type=float, default=CONFIG["figsize"][0])
    p.add_argument("--fig_h", type=float, default=CONFIG["figsize"][1])
    p.add_argument("--lw_bulk", type=float, default=CONFIG["lw_bulk"])
    p.add_argument("--lw_surf", type=float, default=CONFIG["lw_surf"])
    p.add_argument("--alpha_bulk", type=float, default=CONFIG["alpha_bulk"])
    p.add_argument("--alpha_surf", type=float, default=CONFIG["alpha_surf"])
    p.add_argument("--scatter_s_bulk", type=float, default=CONFIG["scatter_s_bulk"])
    p.add_argument("--scatter_s_surf", type=float, default=CONFIG["scatter_s_surf"])
    p.add_argument("--view_elev", type=float, default=CONFIG["view_elev"])
    p.add_argument("--view_azim", type=float, default=CONFIG["view_azim"])
    return p.parse_args()

def main():
    args = parse_args()
    traj, states, time, R = simulate_exchange(
        n_phage=args.n_phage, steps=args.steps, dt=args.dt,
        Rb=args.rb, Rp=args.rp,
        sigma_bulk=args.sigma_bulk, sigma_surf=args.sigma_surf,
        k_off=args.k_off, r_init=args.r_init,
        eps_contact_push=args.eps_contact_push,
        eps_desorb=args.eps_desorb,
        seed=args.seed
    )

    if args.csv:
        save_csv(traj, states, time, args.csv)
    if args.png:
        quick_figure(
            traj, states, R, args.png,
            trail_frames=args.trail_frames,
            figsize=(args.fig_w, args.fig_h),
            lw_bulk=args.lw_bulk, lw_surf=args.lw_surf,
            alpha_bulk=args.alpha_bulk, alpha_surf=args.alpha_surf,
            scatter_s_bulk=args.scatter_s_bulk, scatter_s_surf=args.scatter_s_surf,
            view_elev=args.view_elev, view_azim=args.view_azim
        )

    # quick stats
    surf_frac = states.mean(axis=0).mean()
    print(f"Saved: {args.csv} and {args.png}")
    print(f"Avg. fraction of time on surface (across walkers): {surf_frac:.3f}")

if __name__ == "__main__":
    main()
