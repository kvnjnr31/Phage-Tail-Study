# phage_surface_walk.py
# Tail-anchored surface walking on a single spherical bacterium.
# All tunable parameters live in the CONFIG block below. CLI flags (optional)
# will override these defaults at runtime.

# ============================== CONFIG (EDIT HERE) ==============================

CONFIG = {
    # ---- Simulation / Time ----
    "n_phage": 100,            # number of walkers
    "steps": 100,             # time steps
    "dt": 0.01,               # time step size (s)
    "seed": 7,                # RNG seed; set <0 for nondeterministic

    # ---- Geometry (Bacterium + Phage) ----
    "rb": 1.0,                # bacterium radius
    "head_radius": 0.08,      # phage head radius
    "tail_length": 12,        # beads in tail
    "tail_spacing": 0.06,     # arc spacing between tail beads along surface

    # ---- Dynamics (Surface Walk) ----
    "sigma": 0.06,            # tangent-plane Gaussian step scale at contact point

    # ---- Outputs ----
    "png": "data/image/phage_surface_walk.png",
    "csv": "data/numerical/phage_surface_walk_trajectories.csv",
    "mp4": "",                # set e.g. "data/phage_surface_walk.mp4" to render (requires ffmpeg)

    # ---- Plot / Render ----
    "figure_size": (7, 7),
    "dpi": 220,
    "view_elev": 18,
    "view_azim": 30,
}

# ================================================================================
# Internals below (no need to edit)
# ================================================================================

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from PhageSimpleTail import PhageSimpleTail  # beads-and-stick geometry

# ------------------------- utils ------------------------- #
def ensure_parent_dir(path: str):
    if path:
        d = os.path.dirname(path)
        if d and not os.path.exists(d):
            os.makedirs(d, exist_ok=True)

def random_unit_vectors(n: int) -> np.ndarray:
    v = np.random.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v

def orthonormal_basis(u: np.ndarray):
    """Return two orthonormal tangent vectors (t1, t2) at unit normal u on S^2."""
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, u)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    t1 = np.cross(u, tmp); t1 /= np.linalg.norm(t1)
    t2 = np.cross(u, t1);  t2 /= np.linalg.norm(t2)
    return t1, t2

def random_tangent_step(u: np.ndarray, sigma: float) -> np.ndarray:
    """
    Small Gaussian step in the tangent plane at u, then reproject to S^2.
    Returns the next unit normal u'.
    """
    t1, t2 = orthonormal_basis(u)
    a, b = np.random.normal(scale=sigma), np.random.normal(scale=sigma)
    v = u + a * t1 + b * t2
    return v / np.linalg.norm(v)

def rotate_about_axis(vec: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Rodrigues’ rotation formula: rotate 'vec' around unit 'axis' by angle 'theta' (radians).
    """
    axis = axis / max(np.linalg.norm(axis), 1e-12)
    c, s = np.cos(theta), np.sin(theta)
    return (vec * c
            + np.cross(axis, vec) * s
            + axis * (np.dot(axis, vec) * (1 - c)))

def tail_beads_on_surface(contact_u: np.ndarray, trail_dir: np.ndarray,
                          rb: float, tail_length: int, tail_spacing: float) -> np.ndarray:
    """
    Lay tail beads *on the bacterium surface* (radius rb) along the great-circle
    behind the contact point 'contact_u'. 'trail_dir' is a unit tangent direction
    pointing “backwards”; bead k sits at arc s_k = k * tail_spacing => angle theta_k = s_k / rb.
    """
    n = np.cross(contact_u, trail_dir)
    n_norm = np.linalg.norm(n)
    if n_norm < 1e-12:
        t1, _ = orthonormal_basis(contact_u)
        n = np.cross(contact_u, t1)
        n /= max(np.linalg.norm(n), 1e-12)
    else:
        n /= n_norm

    beads = []
    for k in range(1, tail_length + 1):
        theta = (k * tail_spacing) / rb  # radians of arc along the surface
        dir_k = rotate_about_axis(contact_u, n, -theta)  # rotate backwards along trail
        beads.append(rb * dir_k)
    return np.array(beads)

# ----------------------------- Core simulation ----------------------------- #
def simulate_surface_walk(
    n_phage: int,
    steps: int,
    dt: float,
    rb: float,
    head_radius: float,
    tail_length: int,
    tail_spacing: float,
    sigma: float,
    seed: int | None
):
    """
    Tail-anchored model:
      - The *contact point* of the tail is parameterized by its outward normal u on S^2.
      - u performs a tangent-plane random walk (Brownian on the sphere).
      - Head sits at (rb + head_radius) * u.
      - Tail beads lie ON the surface (radius rb) trailing behind the contact point
        along the local great circle defined by a unit tangent 'trail_dir'.
    """
    if seed is not None and seed >= 0:
        np.random.seed(seed)

    # Initialize contact normals and initial trail directions
    u = random_unit_vectors(n_phage)           # surface normals at contact
    trail_dirs = np.zeros_like(u)              # persistent tangent trail directions
    for i in range(n_phage):
        trail_dirs[i], _ = orthonormal_basis(u[i])

    R_head = rb + head_radius

    # Create phage containers (we use class for geometry fields; dynamics here)
    phage = [
        PhageSimpleTail(
            head_pos=R_head * u[i],
            head_radius=head_radius,
            tail_length=tail_length,
            tail_spacing=tail_spacing
        ) for i in range(n_phage)
    ]

    # Storage (we log head centers for MSD/coverage; add tail dump if needed)
    traj = np.zeros((steps, n_phage, 3), dtype=float)
    traj[0] = R_head * u
    time = np.arange(steps) * dt

    for t in range(1, steps):
        for i in range(n_phage):
            u_prev = u[i]
            u[i] = random_tangent_step(u_prev, sigma)

            # Compute actual tangent step direction to set the trail direction
            step_vec = u[i] - u_prev
            step_vec -= np.dot(step_vec, u[i]) * u[i]  # project to tangent
            norm = np.linalg.norm(step_vec)
            trail_dir = (-step_vec / norm) if norm > 1e-12 else trail_dirs[i]
            trail_dirs[i] = trail_dir  # persist

            # Update head position (rides above contact)
            head_pos = R_head * u[i]
            phage[i].head_pos = head_pos

            # Update tail beads ON the surface trailing the contact
            phage[i].tail_coords = tail_beads_on_surface(
                contact_u=u[i],
                trail_dir=trail_dir,
                rb=rb,
                tail_length=tail_length,
                tail_spacing=tail_spacing
            )

            traj[t, i] = head_pos

    return traj, time

# ----------------------------- I/O helpers ----------------------------- #
def save_csv(traj: np.ndarray, time: np.ndarray, path: str, dpi: int = 220):
    """
    Write trajectories in the same structure as phage_bulk_surface_exchange.py:
      columns: t, id, x, y, z, state
    Here, every row is a surface sample, so state=1 for all rows.
    """
    ensure_parent_dir(path)

    steps, n_phage, _ = traj.shape
    rows = []
    for t in range(steps):
        for i in range(n_phage):
            rows.append({
                "t": float(time[t]),
                "id": int(i),
                "x": float(traj[t, i, 0]),
                "y": float(traj[t, i, 1]),
                "z": float(traj[t, i, 2]),
                "state": 1  # surface
            })

    # keep column order identical to the exchange writer
    pd.DataFrame(rows, columns=["t", "id", "x", "y", "z", "state"]).to_csv(path, index=False)

def plot_png(traj: np.ndarray, rb: float, head_radius: float, png_path: str,
             figure_size=(7,7), dpi=220, view_elev=18, view_azim=30):
    if png_path:
        ensure_parent_dir(png_path)
    plt.rcParams.update({
        "figure.figsize": figure_size,
        "axes.titlesize": 16,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
    })

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # bacterium sphere
    u_s = np.linspace(0, 2*np.pi, 120)
    v_s = np.linspace(0, np.pi, 60)
    xs = rb * np.outer(np.cos(u_s), np.sin(v_s))
    ys = rb * np.outer(np.sin(u_s), np.sin(v_s))
    zs = rb * np.outer(np.ones_like(u_s), np.cos(v_s))
    ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)

    # trajectories (head centers)
    steps, n_phage, _ = traj.shape
    for i in range(n_phage):
        ax.plot(traj[:, i, 0], traj[:, i, 1], traj[:, i, 2], linewidth=1.0)
        ax.scatter(traj[-1, i, 0], traj[-1, i, 1], traj[-1, i, 2], s=12)

    ax.set_title("Tail-Anchored Phage Surface Walking (head follows tail)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    lim = (rb + head_radius) * 1.2
    ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
    ax.view_init(elev=view_elev, azim=view_azim)
    plt.tight_layout()
    if png_path:
        plt.savefig(png_path, dpi=dpi)
    else:
        plt.show()
    plt.close(fig)

def render_mp4(traj: np.ndarray, rb: float, head_radius: float, dt: float, mp4_path: str,
               figure_size=(7,7), fps: int = 30):
    if not mp4_path:
        return
    ensure_parent_dir(mp4_path)

    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111, projection='3d')

    # Precompute sphere
    u_s = np.linspace(0, 2*np.pi, 120)
    v_s = np.linspace(0, np.pi, 60)
    xs = rb * np.outer(np.cos(u_s), np.sin(v_s))
    ys = rb * np.outer(np.sin(u_s), np.sin(v_s))
    zs = rb * np.outer(np.ones_like(u_s), np.cos(v_s))

    steps, n_phage, _ = traj.shape

    def init():
        ax.clear()
        ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)
        lim = (rb + head_radius) * 1.2
        ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.view_init(elev=CONFIG["view_elev"], azim=CONFIG["view_azim"])
        return []

    def update(frame):
        ax.clear()
        ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)
        for i in range(n_phage):
            ax.plot(traj[:frame+1, i, 0], traj[:frame+1, i, 1], traj[:frame+1, i, 2], linewidth=1.0)
            ax.scatter(traj[frame, i, 0], traj[frame, i, 1], traj[frame, i, 2], s=12)
        lim = (rb + head_radius) * 1.2
        ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.set_title(f"Tail-Anchored Surface Walk (t = {frame*dt:.2f}s)")
        ax.view_init(elev=CONFIG["view_elev"], azim=CONFIG["view_azim"])
        return []

    anim = FuncAnimation(fig, update, init_func=init, frames=steps, interval=1000/fps, blit=False)
    writer = FFMpegWriter(fps=fps, bitrate=2400)
    anim.save(mp4_path, writer=writer)
    plt.close(fig)

# ----------------------------- CLI (optional) ----------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        description="Tail-anchored phage surface-walking; CLI overrides CONFIG defaults."
    )
    # simulation / time
    p.add_argument("--n_phage", type=int, default=None)
    p.add_argument("--steps", type=int, default=None)
    p.add_argument("--dt", type=float, default=None)
    p.add_argument("--seed", type=int, default=None)

    # geometry
    p.add_argument("--rb", type=float, default=None)
    p.add_argument("--head_radius", type=float, default=None)
    p.add_argument("--tail_length", type=int, default=None)
    p.add_argument("--tail_spacing", type=float, default=None)

    # dynamics
    p.add_argument("--sigma", type=float, default=None)

    # outputs
    p.add_argument("--png", type=str, default=None)
    p.add_argument("--csv", type=str, default=None)
    p.add_argument("--mp4", type=str, default=None)

    # plot/render
    p.add_argument("--figure_size", type=str, default=None, help="e.g. '7,7'")
    p.add_argument("--dpi", type=int, default=None)
    p.add_argument("--view_elev", type=float, default=None)
    p.add_argument("--view_azim", type=float, default=None)
    return p

def apply_overrides(cfg: dict, args: argparse.Namespace) -> dict:
    new_cfg = dict(cfg)
    for k, v in vars(args).items():
        if v is None:
            continue
        if k == "figure_size" and isinstance(v, str):
            try:
                parts = [float(x) for x in v.split(",")]
                if len(parts) == 2:
                    new_cfg[k] = (parts[0], parts[1])
                    continue
            except Exception:
                pass
        new_cfg[k] = v
    return new_cfg

# -------------------------------- MAIN -------------------------------- #
def main():
    # Parse CLI (optional) and apply on top of CONFIG
    parser = build_parser()
    args = parser.parse_args()
    cfg = apply_overrides(CONFIG, args)

    # Run sim
    seed = cfg["seed"] if (cfg["seed"] is None or cfg["seed"] >= 0) else None
    traj, time = simulate_surface_walk(
        n_phage=cfg["n_phage"],
        steps=cfg["steps"],
        dt=cfg["dt"],
        rb=cfg["rb"],
        head_radius=cfg["head_radius"],
        tail_length=cfg["tail_length"],
        tail_spacing=cfg["tail_spacing"],
        sigma=cfg["sigma"],
        seed=seed
    )

    # Save outputs
    if cfg["csv"]:
        save_csv(traj, time, cfg["csv"])
    plot_png(
        traj, cfg["rb"], cfg["head_radius"], cfg["png"],
        figure_size=tuple(cfg["figure_size"]), dpi=cfg["dpi"],
        view_elev=cfg["view_elev"], view_azim=cfg["view_azim"]
    )
    if cfg["mp4"]:
        try:
            render_mp4(
                traj, cfg["rb"], cfg["head_radius"], cfg["dt"], cfg["mp4"],
                figure_size=tuple(cfg["figure_size"]), fps=30
            )
        except Exception as e:
            print("MP4 render skipped (ffmpeg missing or other error):", e)

if __name__ == "__main__":
    main()
