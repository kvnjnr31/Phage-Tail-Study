# validate_phage_surface_walk.py
# Verifies SURFACE-WALK diagnostics from phage_surface_walk_trajectories.csv
# MSD now uses only sustained surface segments with duration >= --min_surface_dur (sec).

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ---------- style ----------
RC = {
    "figure.figsize": (7.2, 4.4),
    "axes.titlesize": 18,
    "axes.titlepad": 12,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 13,
    "font.family": "DejaVu Sans",
}
plt.rcParams.update(RC)

def _apply_tick_style(ax, nbins=5):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nbins))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=nbins))
    ax.tick_params(axis="both", which="major", length=3.0, width=0.9)

# ---------- helpers ----------
def to_unit(v):
    return v / np.clip(np.linalg.norm(v, axis=-1, keepdims=True), 1e-12, None)

def angular_distance(u, v):
    c = np.sum(u * v, axis=-1)
    return np.arccos(np.clip(c, -1.0, 1.0))

def spherical_coords(xyz):
    r = np.linalg.norm(xyz, axis=-1)
    theta = np.arccos(np.clip(xyz[..., 2] / np.clip(r, 1e-12, None), -1, 1))
    phi = (np.arctan2(xyz[..., 1], xyz[..., 0]) + 2*np.pi) % (2*np.pi)
    return r, theta, phi

def equal_area_bins(n_theta=12, n_phi=24):
    return np.linspace(-1, 1, n_theta+1), np.linspace(0, 2*np.pi, n_phi+1)

def segment_by_time(t, dt, tol=1.5):
    if len(t) == 0:
        return []
    segs = []
    start = 0
    for k in range(1, len(t)):
        if (t[k] - t[k-1]) > tol*dt + 1e-12:
            segs.append((start, k))
            start = k
    segs.append((start, len(t)))
    return segs

# ---------- single-panel plotters ----------
def plot_step_angle_hist(step_angles, ax):
    ax.hist(step_angles, bins=60, density=True)
    ax.set_title("Step Angle Distribution (surface only)", pad=12)
    ax.set_xlabel("angle (rad)")
    ax.set_ylabel("pdf")
    _apply_tick_style(ax)

def plot_angular_msd_from_segments(segments_uv, dt, ax, kfit=50, marker_size=6, extend_fit=True, max_lag=400):
    # accumulate sums and counts per lag
    lag_max = 0
    sums = {}
    counts = {}
    for U in segments_uv:
        m = len(U)
        if m < 2:
            continue
        Lmax = min(max_lag, m-1)
        lag_max = max(lag_max, Lmax)
        for L in range(1, Lmax+1):
            a = angular_distance(U[L:], U[:-L])**2
            sums[L] = sums.get(L, 0.0) + float(np.sum(a))
            counts[L] = counts.get(L, 0) + a.size

    if lag_max == 0:
        ax.text(0.5, 0.5, "Insufficient sustained surface data", ha="center", va="center", transform=ax.transAxes)
        return np.nan, np.nan, np.array([]), np.array([])

    lags = np.arange(1, lag_max+1)
    msd = np.array([sums[L] / counts[L] for L in lags])
    t = lags * dt

    kfit = max(2, min(kfit, len(lags)))
    slope, intercept = np.polyfit(t[:kfit], msd[:kfit], 1)
    D_hat = slope / 4.0

    ax.plot(t, msd, marker='o', ms=marker_size, lw=1.3, label="MSD (surface, sustained)")
    if extend_fit:
        t_fit = np.linspace(t.min(), t.max(), 200)
        ax.plot(t_fit, slope * t_fit + intercept, lw=2.25, label=f"fit (k={kfit})")
    else:
        ax.plot(t[:kfit], slope * t[:kfit] + intercept, lw=2.25, label=f"fit (k={kfit})")

    ax.set_title("Angular MSD on Sphere (sustained surface-only)", pad=12)
    ax.set_xlabel("lag (s)")
    ax.set_ylabel(r"$\langle \alpha^2 \rangle$")
    ax.legend()
    _apply_tick_style(ax)
    return slope, D_hat, lags, msd

def plot_equal_area_heat(mu, phi, ax, n_theta=12, n_phi=24):
    mu_edges, phi_edges = equal_area_bins(n_theta, n_phi)
    H, _, _ = np.histogram2d(mu, phi, bins=[mu_edges, phi_edges])
    im = ax.imshow(
        H.T, origin="lower", aspect="auto",
        extent=[mu_edges[0], mu_edges[-1], phi_edges[0], phi_edges[-1]]
    )
    ax.set_title(r"Equal-Area Occupancy ($\mu=\cos\theta$ vs $\phi$) — surface", pad=12)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\phi$")
    _apply_tick_style(ax)
    return H, im

# ---------- three-panel composer ----------
def summary_three_panel(step_angles, dt, segments_uv, mu_ss, phi_ss, outpng, n_walkers, sigma_for_caption, kfit=50):
    fig, axes = plt.subplots(3, 1, figsize=(8.4, 11.4))
    plt.subplots_adjust(top=0.92, hspace=0.42)

    # Panel A: MSD from sustained segments
    plot_angular_msd_from_segments(segments_uv, dt, axes[0], kfit=kfit, marker_size=6, extend_fit=True)

    # Panel B: Equal-area heatmap (all surface rows)
    _, im = plot_equal_area_heat(mu_ss, phi_ss, axes[1], n_theta=12, n_phi=24)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="4%", pad=0.10)
    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label("count")

    # Panel C: Step-angle PDF (all surface steps)
    plot_step_angle_hist(step_angles, axes[2])

    fig.suptitle(
        rf"Surface Random Walk Diagnostics  (walkers={n_walkers},  $\sigma$={sigma_for_caption})",
        fontsize=18, y=0.98
    )
    fig.savefig(outpng, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(description="Validate phage SURFACE-walk CSV with sustained-segment MSD.")
    parser.add_argument("--csv", type=str, default="data/numerical/phage_surface_walk_trajectories.csv")
    parser.add_argument("--rb", type=float, default=1.0)
    parser.add_argument("--rp", type=float, default=0.05)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--outdir", type=str, default="data/image")
    parser.add_argument("--max_ids", type=int, default=0)
    parser.add_argument("--kfit", type=int, default=50, help="# initial MSD lags used for linear fit.")
    parser.add_argument("--sigma", type=float, default=None, help="Sim σ for caption; if omitted, estimated from steps.")
    parser.add_argument("--min_surface_dur", type=float, default=4.0,
                        help="Minimum continuous on-surface duration (sec) required for MSD segments.")
    args = parser.parse_args()

    csv_path = Path.cwd() / args.csv
    if not csv_path.exists():
        print(f"[ERROR] CSV file not found in current directory: {csv_path.name}")
        return

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Load & FILTER to surface walkers and surface rows
    df = pd.read_csv(csv_path)
    surface_ids = df.groupby("id")["state"].max()
    surface_ids = surface_ids[surface_ids == 1].index.values
    df = df[df["id"].isin(surface_ids)].copy()
    df = df[df["state"] == 1].copy()

    ids_sorted = np.sort(df["id"].unique())
    if args.max_ids and args.max_ids < len(ids_sorted):
        ids_sorted = ids_sorted[:args.max_ids]
        df = df[df["id"].isin(ids_sorted)]

    # geometry (surface rows only)
    R_target = args.rb + args.rp
    r = np.linalg.norm(df[["x","y","z"]].values, axis=1)
    geom_max = float(np.max(np.abs(r - R_target))) if len(r) else np.nan
    geom_mean = float(np.mean(np.abs(r - R_target))) if len(r) else np.nan

    # Build segments per id (surface-only), then keep only sustained segments
    dt = args.dt
    sustained_uv = []
    all_step_angles = []  # still from ALL surface steps for PDF
    kept, dropped = 0, 0

    for pid in ids_sorted:
        sub = df[df["id"] == pid].sort_values("t")
        t = sub["t"].values
        xyz = sub[["x","y","z"]].values
        u = to_unit(xyz)

        # all steps (for PDF)
        if len(u) >= 2:
            all_step_angles.append(angular_distance(u[1:], u[:-1]))

        # segment by time; keep only those meeting duration threshold
        for a, b in segment_by_time(t, dt, tol=1.5):
            if b - a >= 2:
                dur = t[b-1] - t[a]
                if dur >= args.min_surface_dur - 1e-12:
                    sustained_uv.append(u[a:b]); kept += 1
                else:
                    dropped += 1

    step_angles = np.concatenate(all_step_angles) if len(all_step_angles) else np.array([])

    # Equal-area occupancy map — use all surface samples
    xyz_surface = df[["x","y","z"]].values
    _, theta_ss, phi_ss = spherical_coords(xyz_surface)
    mu_ss = np.cos(theta_ss)

    # Sigma caption
    sigma_caption = f"{np.std(step_angles):.3f} (est.)" if (args.sigma is None and step_angles.size) else (
        f"{args.sigma:.3f}" if args.sigma is not None else "n/a"
    )

    # Figures
    msd_png = outdir / "angular_msd_surface_sustained.png"
    fig1, ax1 = plt.subplots(); plt.subplots_adjust(top=0.92)
    plot_angular_msd_from_segments(sustained_uv, dt, ax1, kfit=args.kfit, marker_size=6, extend_fit=True)
    fig1.savefig(msd_png, dpi=300, bbox_inches="tight"); plt.close(fig1)

    heat_png = outdir / "equal_area_heatmap_surface.png"
    fig2, ax2 = plt.subplots()
    _, im = plot_equal_area_heat(mu_ss, phi_ss, ax2, n_theta=12, n_phi=24)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="4%", pad=0.10)
    cbar = fig2.colorbar(im, cax=cax, orientation="vertical"); cbar.set_label("count")
    fig2.savefig(heat_png, dpi=300, bbox_inches="tight"); plt.close(fig2)

    step_png = outdir / "step_angle_hist_surface.png"
    fig3, ax3 = plt.subplots(); plot_step_angle_hist(step_angles, ax3)
    fig3.savefig(step_png, dpi=300, bbox_inches="tight"); plt.close(fig3)

    summary_png = outdir / "summary_3panel_surface_sustained.png"
    summary_three_panel(step_angles, dt, sustained_uv, mu_ss, phi_ss,
                        outpng=summary_png, n_walkers=len(ids_sorted),
                        sigma_for_caption=sigma_caption, kfit=args.kfit)

    # Report
    print("\n")
    report = f"""
    VALIDATION REPORT (Surface-only; sustained MSD)
    ===============================================
    Walkers with surface data        : {len(ids_sorted)}
    Target Radius (Rb+Rp)            : {R_target:.6f}

    Radius Geometry (surface rows)
        max  |r - R|                : {geom_max:.3e}
        mean |r - R|                : {geom_mean:.3e}

    Sustained segment threshold      : {args.min_surface_dur:.3f} s
    Segments kept / dropped          : {kept} / {dropped}

    Figures
        Angular MSD (sustained)     : {msd_png}
        Equal-Area Heatmap          : {heat_png}
        Step PDF                    : {step_png}
        Summary 3×1 (sustained MSD) : {summary_png}
    """
    (outdir/"validation_report_surface_sustained.txt").write_text(report.strip()+"\n", encoding="utf-8")
    print(report)

if __name__ == "__main__":
    main()
