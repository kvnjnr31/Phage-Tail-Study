# phage_surface_walk.py
# Tail-anchored surface walking with receptor adsorption:
# - Tail contact point performs Brownian motion on the bacterium surface.
# - Receptors are absorbing patches on the sphere (geodesic capture).
# - Binding probability per contact is modulated by Mg2+ (schwartz-style effect).
#
# Supports tail-length sweeps (podoviridae/myoviridae/siphoviridae analogs) and
# receptor-density sweeps to reproduce sublinear adsorption-vs-density and the
# reduced Mg2+ sensitivity at high receptor density.
#
# Example:
#   python phage_surface_walk.py --n_phage 200 --steps 6000 --dt 0.002 \
#     --rb 1.0 --head_radius 0.08 --tail_length 4 --tail_spacing 0.06 \
#     --sigma 0.04 --receptor_count 30 --r_cap 0.015 \
#     --p_bind 1.0 --mg_mM 0.2 --mg_Km_mM 0.1 --mg_alpha 1.0 \
#     --csv data/traj.csv --events data/events.csv --png data/plot.png
#
# Base surface-walk and beads-and-stick geometry follow our previous version and the
# minimal tail class in PhageSimpleTail.py (head + linear tail).  (See citations below.)

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from PhageSimpleTail import PhageSimpleTail  # beads-and-stick geometry

# ------------------------- Sphere & tangent helpers ------------------------- #
def random_unit_vectors(n: int) -> np.ndarray:
    v = np.random.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v

def orthonormal_basis(u: np.ndarray):
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(tmp, u)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    t1 = np.cross(u, tmp); t1 /= np.linalg.norm(t1)
    t2 = np.cross(u, t1);  t2 /= np.linalg.norm(t2)
    return t1, t2

def random_tangent_step(u: np.ndarray, sigma: float) -> np.ndarray:
    t1, t2 = orthonormal_basis(u)
    a, b = np.random.normal(scale=sigma), np.random.normal(scale=sigma)
    v = u + a * t1 + b * t2
    return v / np.linalg.norm(v)

def rotate_about_axis(vec: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    axis = axis / max(np.linalg.norm(axis), 1e-12)
    c, s = np.cos(theta), np.sin(theta)
    return (vec * c
            + np.cross(axis, vec) * s
            + axis * (np.dot(axis, vec) * (1 - c)))

# ----------------------------- Tail geometry ----------------------------- #
def tail_beads_on_surface(contact_u: np.ndarray, trail_dir: np.ndarray,
                          rb: float, tail_length: int, tail_spacing: float) -> np.ndarray:
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
        theta = (k * tail_spacing) / rb  # radians of arc
        dir_k = rotate_about_axis(contact_u, n, -theta)
        beads.append(rb * dir_k)
    return np.array(beads)

# ----------------------------- Receptors & binding ----------------------------- #
def place_receptors_on_sphere(rb: float, n_receptors: int, rng=np.random) -> np.ndarray:
    """
    Return unit directions r̂_j to receptors (positions are rb * r̂_j).
    Uniform on sphere by sampling cos(theta) ~ U[-1,1], phi ~ U[0,2π).
    """
    z = rng.uniform(-1.0, 1.0, size=n_receptors)
    phi = rng.uniform(0.0, 2*np.pi, size=n_receptors)
    x = np.sqrt(1 - z*z) * np.cos(phi)
    y = np.sqrt(1 - z*z) * np.sin(phi)
    return np.stack([x, y, z], axis=1)

def mg_binding_factor(mg_mM: float, Km_mM: float, alpha: float) -> float:
    """
    Hill-like Mg2+ modulation for per-contact binding probability (0..1 scale).
    """
    if Km_mM <= 0:
        return 1.0
    frac = mg_mM / (mg_mM + Km_mM)
    return float(np.clip(frac**alpha, 0.0, 1.0))

def geodesic_angle(u: np.ndarray, v: np.ndarray) -> float:
    return float(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0)))

# ----------------------------- Core simulation ----------------------------- #
def simulate_surface_walk_adsorption(
    n_phage: int,
    steps: int,
    dt: float,
    rb: float,
    head_radius: float,
    tail_length: int,
    tail_spacing: float,
    sigma: float,
    receptor_count: int,
    r_cap: float,              # capture radius (linear, same units as rb)
    p_bind: float,             # baseline per-contact binding probability
    mg_mM: float,
    mg_Km_mM: float,
    mg_alpha: float,
    seed: int | None
):
    """
    Tail-anchored surface diffusion with absorbing receptors.
      - Contact point unit vector u_i(t) diffuses on S^2.
      - Head sits at (rb + head_radius)*u_i(t).
      - Receptors are at directions rhat_j; a contact occurs when geodesic angle <= r_cap / rb.
      - Conditional on contact, binding occurs with probability p_bind_eff = p_bind * f_Mg.
      - Upon binding, walker i stops and its binding time is recorded.
    Returns:
      traj[T,N,3] (head centers for unbound walkers; bound ones freeze at binding frame),
      time[T],
      events: DataFrame with columns [id, t_bind, step_bind, receptor_id] (NaN if unbound).
    """
    rng = np.random.default_rng(seed)
    # Initial contact directions
    u = random_unit_vectors(n_phage)
    R_head = rb + head_radius

    # Phage geometry containers
    phage = [
        PhageSimpleTail(
            head_pos=R_head * u[i],
            head_radius=head_radius,
            tail_length=tail_length,
            tail_spacing=tail_spacing
        ) for i in range(n_phage)
    ]

    # Receptors
    rhat = place_receptors_on_sphere(rb, receptor_count, rng=rng)  # unit directions
    theta_cap = r_cap / max(rb, 1e-12)

    # Mg2+ modulation
    p_eff = float(np.clip(p_bind, 0.0, 1.0)) * mg_binding_factor(mg_mM, mg_Km_mM, mg_alpha)
    p_eff = float(np.clip(p_eff, 0.0, 1.0))

    # Storage
    traj = np.zeros((steps, n_phage, 3), dtype=float)
    traj[0] = R_head * u
    time = np.arange(steps) * dt

    bound = np.zeros(n_phage, dtype=bool)
    t_bind = np.full(n_phage, np.nan)
    step_bind = np.full(n_phage, -1, dtype=int)
    rec_hit = np.full(n_phage, -1, dtype=int)

    # persistent trail directions
    trail = np.zeros_like(u)
    for i in range(n_phage):
        trail[i], _ = orthonormal_basis(u[i])

    for t in range(1, steps):
        for i in range(n_phage):
            if bound[i]:
                traj[t, i] = traj[t-1, i]
                continue

            u_prev = u[i]
            u[i] = random_tangent_step(u_prev, sigma)

            # trail direction = opposite of the actual step
            step_vec = u[i] - u_prev
            step_vec -= np.dot(step_vec, u[i]) * u[i]
            nrm = np.linalg.norm(step_vec)
            if nrm > 1e-12:
                trail[i] = -step_vec / nrm

            # Update geometry (head follows tail)
            head_pos = R_head * u[i]
            phage[i].head_pos = head_pos
            phage[i].tail_coords = tail_beads_on_surface(u[i], trail[i], rb, tail_length, tail_spacing)
            traj[t, i] = head_pos

            # --- adsorption test (tail tip contact ≈ u[i]) ---
            # Check geodesic distance from contact to all receptors
            dots = rhat @ u[i]
            dots = np.clip(dots, -1.0, 1.0)
            ang = np.arccos(dots)
            idx = np.where(ang <= theta_cap)[0]
            if idx.size > 0:
                # Encountered ≥1 receptors this step; bind with Bernoulli(p_eff)
                if rng.random() < p_eff:
                    bound[i] = True
                    t_bind[i] = time[t]
                    step_bind[i] = t
                    rec_hit[i] = int(idx[0])

    events = pd.DataFrame({
        "id": np.arange(n_phage, dtype=int),
        "t_bind": t_bind,
        "step_bind": step_bind,
        "receptor_id": rec_hit
    })
    return traj, time, events, rhat

# ----------------------------- I/O helpers ----------------------------- #
def save_traj_csv(traj: np.ndarray, time: np.ndarray, path: str):
    rows = []
    T, N, _ = traj.shape
    for t in range(T):
        for i in range(N):
            rows.append({"t": float(time[t]), "id": int(i),
                         "x": float(traj[t, i, 0]),
                         "y": float(traj[t, i, 1]),
                         "z": float(traj[t, i, 2])})
    pd.DataFrame(rows).to_csv(path, index=False)

def save_events_csv(events: pd.DataFrame, path: str):
    events.to_csv(path, index=False)

def plot_snapshot_and_survival(traj: np.ndarray, rb: float, head_radius: float,
                               events: pd.DataFrame, png_path: str):
    import matplotlib.gridspec as gridspec
    plt.rcParams.update({
        "figure.figsize": (10, 5.6),
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })
    fig = plt.figure()
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1.0], wspace=0.25)

    # 3D snapshot of final positions
    ax = fig.add_subplot(gs[0], projection='3d')
    u_s = np.linspace(0, 2*np.pi, 120)
    v_s = np.linspace(0, np.pi, 60)
    xs = rb * np.outer(np.cos(u_s), np.sin(v_s))
    ys = rb * np.outer(np.sin(u_s), np.sin(v_s))
    zs = rb * np.outer(np.ones_like(u_s), np.cos(v_s))
    ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)

    T, N, _ = traj.shape
    for i in range(N):
        ax.plot(traj[:, i, 0], traj[:, i, 1], traj[:, i, 2], linewidth=0.8)
        ax.scatter(traj[-1, i, 0], traj[-1, i, 1], traj[-1, i, 2], s=10)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("Surface Walk Trajectories")
    lim = (rb + head_radius) * 1.2
    ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])

    # Simple survival curve S(t) = P(unbound by t)
    ax2 = fig.add_subplot(gs[1])
    tb = events["t_bind"].to_numpy()
    tb = tb[~np.isnan(tb)]
    if tb.size > 0:
        tb_sorted = np.sort(tb)
        n0 = events.shape[0]
        surv = 1.0 - np.arange(1, tb_sorted.size + 1) / n0
        ax2.step(tb_sorted, surv, where="post")
        ax2.set_ylim(0, 1)
    ax2.set_title("Survival (Unbound Fraction)")
    ax2.set_xlabel("time")
    ax2.set_ylabel("S(t)")
    plt.tight_layout()
    if png_path:
        plt.savefig(png_path, dpi=220)
    else:
        plt.show()
    plt.close(fig)

# -------------------------------- CLI -------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(description="Tail-anchored phage surface-walking with receptor adsorption.")
    # sim/time
    p.add_argument("--n_phage", type=int, default=200)
    p.add_argument("--steps", type=int, default=6000)
    p.add_argument("--dt", type=float, default=0.002)
    p.add_argument("--seed", type=int, default=7)

    # geometry
    p.add_argument("--rb", type=float, default=1.0, help="Bacterium radius")
    p.add_argument("--head_radius", type=float, default=0.08, help="Phage head radius")
    p.add_argument("--tail_length", type=int, default=12, help="Tail beads (family proxy)")
    p.add_argument("--tail_spacing", type=float, default=0.06, help="Arc spacing along surface")

    # diffusion
    p.add_argument("--sigma", type=float, default=0.04, help="Tangent-plane step scale (controls Ds)")

    # receptors
    p.add_argument("--receptor_count", type=int, default=30, help="Number of receptors on the cell")
    p.add_argument("--r_cap", type=float, default=0.015, help="Capture radius (linear units)")

    # chemistry (Mg2+)
    p.add_argument("--p_bind", type=float, default=1.0, help="Baseline per-contact bind probability")
    p.add_argument("--mg_mM", type=float, default=0.2, help="Mg2+ concentration (mM)")
    p.add_argument("--mg_Km_mM", type=float, default=0.1, help="Half-effect conc. (mM)")
    p.add_argument("--mg_alpha", type=float, default=1.0, help="Hill exponent for Mg2+ effect")

    # outputs
    p.add_argument("--png", type=str, default="data/phage_surface_walk.png")
    p.add_argument("--csv", type=str, default="data/phage_surface_walk_trajectories.csv")
    p.add_argument("--events", type=str, default="data/phage_adsorption_events.csv")
    p.add_argument("--mp4", type=str, default="")
    return p.parse_args()

def main():
    args = parse_args()
    seed = None if (args.seed is None or args.seed < 0) else args.seed

    traj, time, events, rhat = simulate_surface_walk_adsorption(
        n_phage=args.n_phage,
        steps=args.steps,
        dt=args.dt,
        rb=args.rb,
        head_radius=args.head_radius,
        tail_length=args.tail_length,
        tail_spacing=args.tail_spacing,
        sigma=args.sigma,
        receptor_count=args.receptor_count,
        r_cap=args.r_cap,
        p_bind=args.p_bind,
        mg_mM=args.mg_mM,
        mg_Km_mM=args.mg_Km_mM,
        mg_alpha=args.mg_alpha,
        seed=seed
    )

    if args.csv:
        save_traj_csv(traj, time, args.csv)
    if args.events:
        save_events_csv(events, args.events)
    if args.png:
        plot_snapshot_and_survival(traj, args.rb, args.head_radius, events, args.png)
    if args.mp4:
        try:
            # simple orbiting camera animation (optional)
            fig = plt.figure(figsize=(7,7))
            ax = fig.add_subplot(111, projection='3d')
            u_s = np.linspace(0, 2*np.pi, 120)
            v_s = np.linspace(0, np.pi, 60)
            xs = args.rb * np.outer(np.cos(u_s), np.sin(v_s))
            ys = args.rb * np.outer(np.sin(u_s), np.sin(v_s))
            zs = args.rb * np.outer(np.ones_like(u_s), np.cos(v_s))
            T, N, _ = traj.shape
            def init():
                ax.clear()
                ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)
                lim = (args.rb + args.head_radius) * 1.2
                ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
                ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
                return []
            def update(frame):
                ax.clear()
                ax.plot_surface(xs, ys, zs, alpha=0.18, linewidth=0, shade=True)
                for i in range(N):
                    ax.plot(traj[:frame+1, i, 0], traj[:frame+1, i, 1], traj[:frame+1, i, 2], linewidth=0.8)
                lim = (args.rb + args.head_radius) * 1.2
                ax.set_xlim([-lim, lim]); ax.set_ylim([-lim, lim]); ax.set_zlim([-lim, lim])
                ax.view_init(elev=18, azim=30 + 0.2*frame)
                return []
            anim = FuncAnimation(fig, update, init_func=init, frames=T, interval=20, blit=False)
            anim.save(args.mp4, writer=FFMpegWriter(fps=40, bitrate=3000))
            plt.close(fig)
        except Exception as e:
            print("MP4 render skipped:", e)

if __name__ == "__main__":
    main()
