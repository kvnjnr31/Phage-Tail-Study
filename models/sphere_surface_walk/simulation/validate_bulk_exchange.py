# validate_bulk_exchange.py
# Kinetic diagnostics for phage_bulk_surface_exchange.py outputs.
# Creates:
#   - contact_fraction_over_time.png
#   - residence_time_hist.png
#   - bulk_exchange_summary.csv (per-phage table)
#   - bulk_exchange_residence_times.csv (all segments)
#   - summary_4panel_bulk.png  [2x2 summary panel]
#
# Input CSV schema (required): t,id,x,y,z,state
# Notes:
#   * Exchange CSV contains only phage that ever touched the surface.
#   * dt is inferred from timestamps unless overridden.

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ---------------------- plotting style ----------------------
RC = {
    "figure.figsize": (7.6, 4.4),
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
}
plt.rcParams.update(RC)

# ---------------------- utils ----------------------
def ensure_dir(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)

def infer_dt(df: pd.DataFrame) -> float:
    dts = []
    for _, sub in df.groupby("id"):
        t = np.asarray(sub["t"])
        if len(t) > 1:
            d = np.diff(np.sort(t))
            d = d[d > 0]
            if d.size:
                dts.append(np.median(d))
    return float(np.median(dts)) if dts else np.nan

def segment_residence_times(t: np.ndarray, state: np.ndarray) -> list:
    """Return list of contiguous SURF (state==1) run durations for one id."""
    if t.size == 0:
        return []
    idx = np.argsort(t)
    t = t[idx]
    s = state[idx].astype(int)
    durs = []
    run = False
    a = 0
    for k in range(len(s)):
        if s[k] == 1 and not run:
            run = True
            a = k
        end_series = (k == len(s) - 1)
        if run and (s[k] == 0 or end_series):
            b = k if s[k] == 0 else k
            if b > a:
                dur = float(np.diff(t[a:b+1]).sum())
                if dur > 0:
                    durs.append(dur)
            run = False
    return durs

def adsorption_desorption_counts(state: np.ndarray, t: np.ndarray) -> tuple[int,int]:
    idx = np.argsort(t)
    s = state[idx].astype(int)
    ads = int(np.sum((s[1:] == 1) & (s[:-1] == 0)))
    des = int(np.sum((s[1:] == 0) & (s[:-1] == 1)))
    return ads, des

def survival_from_samples(x: np.ndarray):
    """Empirical survival S(t)=1-CDF(t) for positive samples x."""
    if x.size == 0:
        return np.array([]), np.array([])
    x = np.sort(x)
    n = x.size
    y = 1.0 - np.arange(1, n + 1) / n
    return x, y

# ---------------------- plots (single) ----------------------
def plot_contact_timeseries(df: pd.DataFrame, out_png: Path):
    ts = (df.groupby("t")["state"].mean()).reset_index()
    fig, ax = plt.subplots()
    ax.plot(ts["t"].values, ts["state"].values, lw=1.8)
    ax.set_title("Population Contact Fraction (tracked subset)")
    ax.set_xlabel("time")
    ax.set_ylabel("fraction on surface")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.grid(True, alpha=0.3)
    ensure_dir(out_png)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return float(ts["state"].mean()), float(ts["state"].median()), ts

def plot_residence_hist(res_times: np.ndarray, out_png: Path):
    fig, ax = plt.subplots()
    if res_times.size:
        ax.hist(res_times, bins=50, density=True)
        ax.set_xscale("log")
    ax.set_title("Residence-Time Distribution (surface)")
    ax.set_xlabel("residence time")
    ax.set_ylabel("pdf")
    ax.grid(True, which="both", alpha=0.3)
    ensure_dir(out_png)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------------------- 2x2 summary panel ----------------------
def make_panel(ts_df: pd.DataFrame,
               per_df: pd.DataFrame,
               res_df: pd.DataFrame,
               out_png: Path,
               suptitle: str = "Bulk Exchange Diagnostics — Summary"):
    ensure_dir(out_png)
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    plt.subplots_adjust(wspace=0.28, hspace=0.36, top=0.90)

    # Panel A: contact fraction over time
    ax = axs[0, 0]
    if len(ts_df):
        ax.plot(ts_df["t"].values, ts_df["state"].values, lw=1.8)
        ax.set_title("A) Contact Fraction Over Time")
        ax.set_xlabel("time")
        ax.set_ylabel("fraction on surface")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.grid(True, alpha=0.3)

    # Panel B: per-phage contact % histogram
    ax = axs[0, 1]
    if len(per_df):
        vals = per_df["pct_surface"].values
        ax.hist(vals, bins=30, range=(0, 1), density=True)
        ax.axvline(np.mean(vals), ls="--", lw=1.6, label=f"mean={np.mean(vals):.2f}")
        ax.axvline(np.median(vals), ls=":", lw=1.6, label=f"median={np.median(vals):.2f}")
        ax.set_title("B) Per-Phage Contact Percentage")
        ax.set_xlabel("fraction of time on surface")
        ax.set_ylabel("pdf")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3)

    # Panel C: residence time histogram (log-x)
    ax = axs[1, 0]
    if len(res_df):
        x = res_df["residence_time"].values
        if x.size:
            ax.hist(x, bins=50, density=True)
            ax.set_xscale("log")
        ax.set_title("C) Residence-Time Distribution")
        ax.set_xlabel("residence time")
        ax.set_ylabel("pdf")
        ax.grid(True, which="both", alpha=0.3)

    # Panel D: survival with exponential fit
    ax = axs[1, 1]
    if len(res_df):
        x = res_df["residence_time"].values
        x = x[x > 0]
        xs, S = survival_from_samples(x)
        if xs.size:
            ax.plot(xs, S, lw=1.8, label="empirical S(t)")
            lam = 1.0 / np.mean(x)  # MLE for exponential
            tfit = np.linspace(xs.min(), xs.max(), 300)
            ax.plot(tfit, np.exp(-lam * tfit), lw=1.8, ls="--", label=f"exp fit λ={lam:.3g}")
            ax.set_xscale("log")
            ax.set_yscale("log")
        ax.set_title("D) Residence-Time Survival")
        ax.set_xlabel("t")
        ax.set_ylabel("S(t) = 1 − CDF(t)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(frameon=False)

    fig.suptitle(suptitle, fontsize=17)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------------------- main ----------------------
def main():
    ap = argparse.ArgumentParser(description="Bulk-exchange kinetics validation and summary panel.")
    ap.add_argument("--csv", type=str, default="data/numerical/phage_surface_walk_trajectories.csv",
                    help="CSV produced by phage_bulk_surface_exchange.py")
    ap.add_argument("--outdir_img", type=str, default="data/image",
                    help="Directory for figures")
    ap.add_argument("--summary_csv", type=str, default="data/numerical/bulk_exchange_summary.csv",
                    help="Per-phage summary table output")
    ap.add_argument("--residence_csv", type=str, default="data/numerical/bulk_exchange_residence_times.csv",
                    help="All residence segments table")
    ap.add_argument("--panel_png", type=str, default="data/image/summary_4panel_bulk.png",
                    help="2x2 summary panel image")
    ap.add_argument("--dt", type=float, default=None, help="Override dt (otherwise inferred from t).")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    if not csv_path.exists():
        print(f"[ERROR] CSV not found: {csv_path}")
        return

    df = pd.read_csv(csv_path)
    required = {"t","id","x","y","z","state"}
    if not required.issubset(df.columns):
        print("[ERROR] CSV must have columns: t,id,x,y,z,state")
        return

    # infer dt if not given (reported for reference)
    dt = args.dt if args.dt is not None else infer_dt(df)

    # 1) contact fraction series
    out_ts_png = Path(args.outdir_img) / "contact_fraction_over_time.png"
    mean_frac, median_frac, ts_df = plot_contact_timeseries(df, out_ts_png)

    # 2) per-phage stats + residence segments
    per_rows = []
    all_residence_rows = []

    for pid, sub in df.groupby("id"):
        t = np.asarray(sub["t"])
        s = np.asarray(sub["state"]).astype(int)

        pct_surface = float(np.mean(s == 1)) if s.size else np.nan
        n_ads, n_des = adsorption_desorption_counts(s, t)
        res_list = segment_residence_times(t, s)
        res_arr = np.array(res_list, float) if res_list else np.array([], float)

        per_rows.append({
            "id": int(pid),
            "pct_surface": pct_surface,
            "n_adsorptions": int(n_ads),
            "n_desorptions": int(n_des),
            "n_residence_segments": int(res_arr.size),
            "total_surface_time": float(res_arr.sum()) if res_arr.size else 0.0,
            "mean_res_time": float(res_arr.mean()) if res_arr.size else np.nan,
            "median_res_time": float(np.median(res_arr)) if res_arr.size else np.nan,
            "max_res_time": float(res_arr.max()) if res_arr.size else np.nan,
        })

        for r in res_list:
            all_residence_rows.append({"id": int(pid), "residence_time": float(r)})

    per_df = pd.DataFrame(per_rows).sort_values("id")
    res_df = pd.DataFrame(all_residence_rows)

    # write tables
    summary_csv = Path(args.summary_csv)
    residence_csv = Path(args.residence_csv)
    ensure_dir(summary_csv)
    ensure_dir(residence_csv)
    per_df.to_csv(summary_csv, index=False)
    res_df.to_csv(residence_csv, index=False)

    # single plots
    out_res_png = Path(args.outdir_img) / "residence_time_hist.png"
    plot_residence_hist(res_df["residence_time"].values if len(res_df) else np.array([]), out_res_png)

    # panel: summary_4panel_bulk.png
    make_panel(ts_df, per_df, res_df, Path(args.panel_png),
               suptitle=f"Bulk Exchange Diagnostics — Summary  (tracked={per_df['id'].nunique()}, dt={dt:.4g})")

    # console report
    n_ids = per_df["id"].nunique() if len(per_df) else 0
    n_segments = len(res_df)
    mean_res = float(res_df["residence_time"].mean()) if n_segments else np.nan
    med_res  = float(res_df["residence_time"].median()) if n_segments else np.nan

    print("\nBULK-EXCHANGE KINETICS REPORT")
    print("================================")
    print(f"CSV: {csv_path}")
    print(f"Inferred dt: {dt if args.dt is None else f'{args.dt} (override)'}")
    print(f"Tracked phage (ever-surface): {n_ids}")
    print(f"Population contact fraction over time (mean / median): {mean_frac:.3f} / {median_frac:.3f}")
    print(f"Residence segments: {n_segments}  |  mean={mean_res:.4g}, median={med_res:.4g}")
    print("Saved figures:")
    print(f"  - {out_ts_png}")
    print(f"  - {out_res_png}")
    print(f"  - {args.panel_png}")
    print("Saved tables:")
    print(f"  - {summary_csv}")
    print(f"  - {residence_csv}\n")

if __name__ == "__main__":
    main()
