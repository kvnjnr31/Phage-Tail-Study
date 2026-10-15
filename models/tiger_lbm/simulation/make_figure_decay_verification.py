# make_figure_decay_r_over_b_1_to_100.py
# -------------------------------------------------------------------
# Compute ⟨|u|⟩(r) on the YZ plane (⊥ x̂) for r/b ∈ [1, 100],
# save CSV + (optional) log–log plot with slope guides.
#
# Run:
#   python simulation/make_figure_decay_r_over_b_1_to_100.py
# -------------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
import csv

# =========================== USER SETTINGS ===========================

# Averaging convention:
MAGNITUDE_AFTER_AVG = False   # False → ⟨|u|⟩ (paper-like), True → |⟨u⟩|

# Frames to include (None = all):
AVG_LAST_K = 300

# Distance normalization b:
NORM_BY = "bundle"            # "bundle" (b=0.85*Rb), "Rb", or a float (cells)
BUNDLE_FACTOR = 0.85          # if NORM_BY == "bundle": b = BUNDLE_FACTOR * Rb

# Binning across r/b ∈ [1, 100]:
NBINS = 40                    # log-spaced bin count

# What to do beyond domain limits:
CLIP_TO_DOMAIN     = True     # True → bins outside the domain become NaN
EXTRAPOLATE_TAIL   = False    # if CLIP_TO_DOMAIN is False, extrapolate with power laws
TAIL_SLOPE_FAR     = -2.0     # slope for far field extrapolation
TAIL_SLOPE_NEAR    = -3.0     # optional inner extrapolation if needed

# Plot?
MAKE_PLOT   = True
FIGSIZE     = (6.4, 5.0)
DPI         = 300
TITLE       = r"Decay of time-averaged flow speed on the YZ plane"
XLAB        = r"distance: $r/b$"
YLAB        = r"$\langle |\,\mathbf{u}\,| \rangle_t$  (lattice units)"

# ====================================================================

def detect_project_root(start=None):
    if start is None:
        start = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cur = os.path.abspath(start)
    for _ in range(6):
        if os.path.isdir(os.path.join(cur, "data", "numerical")):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return os.path.abspath(start)

PROJ     = detect_project_root()
NUM_DIR  = os.path.join(PROJ, "data", "numerical")
IMG_DIR  = os.path.join(PROJ, "data", "image")
CFG_PATH = os.path.join(PROJ, "simulation", "tiger_simulation_config.yaml")
os.makedirs(IMG_DIR, exist_ok=True)

def load_Rb_from_config(default=10.0):
    try:
        import yaml
        with open(CFG_PATH, "r") as f:
            cfg = yaml.safe_load(f)
        return float(cfg.get("bacterium", {}).get("Rb", default))
    except Exception:
        return float(default)

def find_saved_outputs():
    def p(n): return os.path.join(NUM_DIR, n)
    files = {
        "velocity_field":       p("velocity_field.npy"),
        "velocity_field_small": p("velocity_field_small.npy"),
        "bacterium_path":       p("bacterium_path.npy"),
    }
    for k,v in list(files.items()):
        if not os.path.exists(v): files[k] = None
    return files

def load_field_and_center():
    files = find_saved_outputs()
    if files["velocity_field"]:
        U = np.load(files["velocity_field"], allow_pickle=True)
    elif files["velocity_field_small"]:
        U = np.load(files["velocity_field_small"], allow_pickle=True)
    else:
        raise FileNotFoundError("No velocity_field(.npy/_small.npy) found.")
    if files["bacterium_path"] is not None:
        center = np.asarray(np.load(files["bacterium_path"]), float)[-1]
    else:
        nx, ny, nz = U.shape[1:4]
        center = np.array([nx/2, ny/2, nz/2], float)
    return U, center

def time_average(U):
    if AVG_LAST_K is None or AVG_LAST_K >= U.shape[0]:
        Usub = U
    else:
        Usub = U[-AVG_LAST_K:]
    if MAGNITUDE_AFTER_AVG:
        u_mean = Usub.mean(axis=0)                 # (nx,ny,nz,3)
        v = np.sqrt((u_mean**2).sum(axis=-1))      # |⟨u⟩|
    else:
        spd = np.sqrt((Usub**2).sum(axis=-1))      # |u|
        v = spd.mean(axis=0)                        # ⟨|u|⟩
    return v

def get_b_norm(Rb):
    if isinstance(NORM_BY, (int, float)):
        return float(NORM_BY)
    if str(NORM_BY).lower() == "rb":
        return float(Rb)
    return float(BUNDLE_FACTOR) * float(Rb)

def radial_profile_fixed_range(v3d, center, b_norm, rmin_b=1.0, rmax_b=100.0, nbins=40):
    """
    Median ⟨|u|⟩ vs r/b in log-spaced bins from rmin_b to rmax_b on the YZ plane.
    Returns (r_centers, v_vals, mask_valid) with NaNs for bins outside domain if CLIP_TO_DOMAIN.
    """
    nx, ny, nz = v3d.shape
    ix = int(np.clip(np.round(center[0]), 0, nx-1))
    yc, zc = float(center[1]), float(center[2])

    V = v3d[ix, :, :]
    yy = np.arange(ny)[:, None]
    zz = np.arange(nz)[None, :]
    R_cells = np.sqrt((yy - yc)**2 + (zz - zc)**2)
    r_over_b = R_cells / b_norm

    edges = np.logspace(np.log10(rmin_b), np.log10(rmax_b), nbins+1)
    centers = np.sqrt(edges[:-1] * edges[1:])
    vals = np.full(nbins, np.nan, dtype=float)
    valid = np.zeros(nbins, dtype=bool)

    # Domain max in r/b (on this plane)
    rmax_domain = np.nanmax(r_over_b) if np.isfinite(r_over_b).any() else 0.0

    for k in range(nbins):
        lo, hi = edges[k], edges[k+1]
        if lo > rmax_domain:
            continue  # beyond domain
        mask = (r_over_b >= lo) & (r_over_b < min(hi, rmax_domain))
        if np.any(mask):
            vals[k] = np.median(V[mask])
            valid[k] = True

    # Optional extrapolation
    if not CLIP_TO_DOMAIN and not valid.all():
        # Far-tail extrapolation with slope = TAIL_SLOPE_FAR anchored at last valid
        if valid.any():
            i_last = np.where(valid)[0][-1]
            r0, v0 = centers[i_last], vals[i_last]
            for i in range(i_last+1, nbins):
                vals[i] = v0 * (centers[i]/r0)**(TAIL_SLOPE_FAR)
                valid[i] = True
        # Inner extrapolation with slope = TAIL_SLOPE_NEAR anchored at first valid
        if valid.any():
            i_first = np.where(valid)[0][0]
            r0, v0 = centers[i_first], vals[i_first]
            for i in range(0, i_first):
                vals[i] = v0 * (centers[i]/r0)**(TAIL_SLOPE_NEAR)
                valid[i] = True

    return centers, vals, valid

def main():
    Rb = load_Rb_from_config(default=10.0)
    b_norm = get_b_norm(Rb)

    U, center = load_field_and_center()
    v3d = time_average(U)

    r_c, v_c, valid = radial_profile_fixed_range(
        v3d, center, b_norm, rmin_b=1.0, rmax_b=100.0, nbins=NBINS
    )

    # Save CSV (r/b, velocity, valid_flag)
    out_csv = os.path.join(IMG_DIR, "decay_profile_r_over_b_1_to_100.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r_over_b", "velocity", "in_domain"])
        for r, v, ok in zip(r_c, v_c, valid):
            w.writerow([f"{r:.8g}", f"{v:.8g}" if np.isfinite(v) else "", int(ok)])
    print("Saved:", out_csv)

    if MAKE_PLOT:
        fig, ax = plt.subplots(figsize=FIGSIZE)
        # Plot valid as points; invalid (NaN) won’t show
        ax.loglog(r_c[valid], v_c[valid], "o", color="k", markersize=4.0, label="simulation")
        # Indicate out-of-domain bins (if any)
        if (~valid).any():
            ax.loglog(r_c[~valid], np.full((~valid).sum(), np.nan), "o", color="gray", alpha=0.4)

        # Optional slope guides through a mid anchor (only if we have valid data across).
        if valid.any():
            r_valid, v_valid = r_c[valid], v_c[valid]
            r0 = np.exp((np.log(r_valid.min()) + np.log(r_valid.max()))/2.0)
            i0 = np.argmin(np.abs(r_valid - r0))
            r0, v0 = r_valid[i0], v_valid[i0]
            xr = np.array([r_c.min(), r_c.max()])
            ax.plot(xr, v0*(xr/r0)**(-2.0), "k-",  lw=2.0, label="slope of -2")
            ax.plot(xr, v0*(xr/r0)**(-3.0), "k:",  lw=2.5, label="slope of -3")

        ax.set_xlabel(XLAB); ax.set_ylabel(YLAB); ax.set_title(TITLE)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(frameon=False, loc="upper right")
        out_png = os.path.join(IMG_DIR, "TIGER_decay_r_over_b_1_to_100.png")
        out_tif = os.path.splitext(out_png)[0] + ".tiff"
        fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
        fig.savefig(out_tif, dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print("Saved:", out_png, "and", out_tif)

    # Domain sanity note
    max_possible = np.max(np.sqrt(((np.array(v3d.shape)[1:] - 1)/2.0)**2)) / b_norm
    print(f"Max r/b reachable on this plane (approx): ~{max_possible:.2f}")
    if CLIP_TO_DOMAIN and (np.max(r_c[valid]) < 100.0):
        print("Note: bins beyond the reachable domain are left blank in the CSV.")
        print("      To fill them with power-law extrapolation, set EXTRAPOLATE_TAIL=True and CLIP_TO_DOMAIN=False.")

if __name__ == "__main__":
    main()
