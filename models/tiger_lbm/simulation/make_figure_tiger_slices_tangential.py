# make_figure_tiger_slices_tangential.py
# Slices of |u| (XY, XZ, YZ) + near-surface map on r≈Rb.
# Run:
#   python simulation/make_figure_tiger_slices_tangential.py

import os
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

# =============================== USER SETTINGS ===============================

# --- smoothing / de-grain ---
PHASE_LOCK      = True
PHASE_BINS      = 24
PHASE_TARGET    = 0.0
AVG_LAST_K      = 80

UPSAMPLE        = 2
INTERP          = "bicubic"   # "bicubic", "bilinear", "nearest"

# --- near-surface sampling ---
MAP_KIND        = "ut"        # "ut", "uphi", "utheta"
N_THETA         = 360
N_PHI           = 720

# --- cropping (centered boxes) ---
CROP_UNITS      = "rb"        # "rb" or "cells"
BOX_XY          = (1.0, 1.0)  # half-sizes (Δx,Δy)
BOX_XZ          = (1.0, 1.0)  # half-sizes (Δx,Δz)
BOX_YZ          = (1.0, 1.0)  # half-sizes (Δy,Δz)

ADD_GRID        = True
SHOW_BOX        = False

# --- typography / layout ---
FIG_DPI         = 300
FIG_SLICE_SIZE  = (6.6, 4.8)   # XY, XZ, YZ
FIG_SURF_SIZE   = (6.6, 4.8)   # match slice height → bottom-right == bottom-left
TITLE_PAD       = 14           # extra gap between title and plot

FONT_SIZES = {
    "base": 13,
    "title": 20,
    "label": 16,
    "tick": 12,
    "cbar_label": 16,
}

# ============================================================================

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

PROJ    = detect_project_root()
NUM_DIR = os.path.join(PROJ, "data", "numerical")
IMG_DIR = os.path.join(PROJ, "data", "image")
CFG_PATH = os.path.join(PROJ, "simulation", "tiger_simulation_config.yaml")
os.makedirs(IMG_DIR, exist_ok=True)

# Global typography + robust layout to avoid colorbar overlap
plt.rcParams.update({
    "font.size": FONT_SIZES["base"],
    "axes.titlesize": FONT_SIZES["title"],
    "axes.labelsize": FONT_SIZES["label"],
    "xtick.labelsize": FONT_SIZES["tick"],
    "ytick.labelsize": FONT_SIZES["tick"],
    "figure.constrained_layout.use": True,
    "axes.titlepad": TITLE_PAD,   # <— extra space above plots
})

def find_saved_outputs():
    d = os.path.join(PROJ, "data", "numerical")
    def p(n): return os.path.join(d, n)
    files = {
        "velocity_field":       p("velocity_field.npy"),
        "velocity_field_small": p("velocity_field_small.npy"),
        "bacterium_path":       p("bacterium_path.npy"),
        "rotation_meta":        p("rotation_meta.npy"),
    }
    for k,v in list(files.items()):
        if not os.path.exists(v): files[k] = None
    return files

def load_Rb_from_config(default=10.0):
    try:
        import yaml
        with open(CFG_PATH, "r") as f:
            cfg = yaml.safe_load(f)
        return float(cfg.get("bacterium", {}).get("Rb", default))
    except Exception:
        return float(default)

def load_field_and_center():
    files = find_saved_outputs()
    if files["velocity_field"]:
        U = np.load(files["velocity_field"], allow_pickle=True)
    elif files["velocity_field_small"]:
        U = np.load(files["velocity_field_small"], allow_pickle=True)
    else:
        nx, ny, nz = 64, 48, 48
        x = np.linspace(-1, 1, nx)
        y = np.linspace(-0.75, 0.75, ny)
        z = np.linspace(-0.75, 0.75, nz)
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        u = -Y * np.exp(-3*(X**2 + Y**2 + Z**2))
        v =  X * np.exp(-3*(X**2 + Y**2 + Z**2))
        w =  0.3*Z * np.exp(-3*(X**2 + Y**2 + Z**2))
        U = np.stack([u, v, w], axis=-1)[None, ...]
    # center
    files = find_saved_outputs()
    if files["bacterium_path"]:
        center = np.asarray(np.load(files["bacterium_path"])[-1], dtype=float)
        nx, ny, nz = U.shape[1:4]
        center = np.clip(center, [0,0,0], [nx-1, ny-1, nz-1])
    else:
        nx, ny, nz = U.shape[1:4]
        center = np.array([nx//2, ny//2, nz//2], dtype=float)
    return U, center

def magnitude(U): return np.sqrt((U**2).sum(axis=-1))

def slice_indices(center, shape):
    nx, ny, nz = shape
    ix = int(np.clip(round(center[0]), 0, nx-1))
    iy = int(np.clip(round(center[1]), 0, ny-1))
    iz = int(np.clip(round(center[2]), 0, nz-1))
    return ix, iy, iz

# ---------- phase-locked averaging ----------
def load_rotation_meta():
    files = find_saved_outputs()
    if not files["rotation_meta"]: return None
    try:
        omega_b, omega_body, dt = np.load(files["rotation_meta"])
        return float(omega_b), float(omega_body), float(dt)
    except Exception:
        return None

def phase_locked_average(U, omega, dt, nbins=24, target=0.0):
    T = U.shape[0]
    n = np.arange(T)
    phi = (omega * dt * n) % (2*np.pi)
    edges = np.linspace(0, 2*np.pi, nbins+1)
    centers = (edges[:-1] + edges[1:]) * 0.5
    j = int(np.argmin(np.abs(((centers - target + np.pi) % (2*np.pi)) - np.pi)))
    lo, hi = edges[j], edges[(j+1) % nbins]
    mask = (phi >= lo) & (phi < hi) if hi > lo else ((phi >= lo) | (phi < hi))
    idx = np.where(mask)[0]
    if len(idx) == 0:
        K = min(AVG_LAST_K, T)
        return U[-K:].mean(axis=0), idx
    return U[idx].mean(axis=0), idx

# ---------- near-surface sampling ----------
def trilinear_sample_vector(field, pts):
    nx, ny, nz, _ = field.shape
    x, y, z = pts[:,0], pts[:,1], pts[:,2]
    x = np.clip(x, 0, nx-1); y = np.clip(y, 0, ny-1); z = np.clip(z, 0, nz-1)
    x0 = np.floor(x).astype(int); y0 = np.floor(y).astype(int); z0 = np.floor(z).astype(int)
    x1 = np.clip(x0+1, 0, nx-1);  y1 = np.clip(y0+1, 0, ny-1);  z1 = np.clip(z0+1, 0, nz-1)
    xd = (x - x0)[:,None]; yd = (y - y0)[:,None]; zd = (z - z0)[:,None]
    c000 = field[x0, y0, z0]; c100 = field[x1, y0, z0]
    c010 = field[x0, y1, z0]; c110 = field[x1, y1, z0]
    c001 = field[x0, y0, z1]; c101 = field[x1, y0, z1]
    c011 = field[x0, y1, z1]; c111 = field[x1, y1, z1]
    c00 = c000*(1-xd) + c100*xd
    c10 = c010*(1-xd) + c110*xd
    c01 = c001*(1-xd) + c101*xd
    c11 = c011*(1-xd) + c111*xd
    c0  = c00*(1-yd) + c10*yd
    c1  = c01*(1-yd) + c11*yd
    return c0*(1-zd) + c1*zd

def tangential_map_trilinear(U_last, center, Rb, n_theta=N_THETA, n_phi=N_PHI, kind=MAP_KIND, smooth=True):
    theta = np.linspace(0.0, np.pi, n_theta, endpoint=False) + (np.pi/n_theta)/2
    phi   = np.linspace(-np.pi, np.pi, n_phi, endpoint=False) + (2*np.pi/n_phi)/2
    T, P  = np.meshgrid(theta, phi, indexing="ij")
    rx = np.cos(T);        ry = np.sin(T)*np.cos(P); rz = np.sin(T)*np.sin(P)
    tx = -np.sin(T);       ty =  np.cos(T)*np.cos(P); tz =  np.cos(T)*np.sin(P)
    px = 0.0*T;            py = -np.sin(P);          pz =  np.cos(P)
    sx = center[0] + Rb*rx; sy = center[1] + Rb*ry; sz = center[2] + Rb*rz
    pts = np.column_stack([sx.ravel(), sy.ravel(), sz.ravel()])
    U = trilinear_sample_vector(U_last, pts).reshape(n_theta, n_phi, 3)
    if kind == "ut":
        dot = U[...,0]*rx + U[...,1]*ry + U[...,2]*rz
        Ut  = np.sqrt((U[...,0]-dot*rx)**2 + (U[...,1]-dot*ry)**2 + (U[...,2]-dot*rz)**2)
    elif kind == "uphi":
        Ut = U[...,0]*px + U[...,1]*py + U[...,2]*pz
    elif kind == "utheta":
        Ut = U[...,0]*tx + U[...,1]*ty + U[...,2]*tz
    else:
        raise ValueError("MAP_KIND must be 'ut', 'uphi', or 'utheta'")
    if smooth:
        ker = np.array([[1,2,1],[2,4,2],[1,2,1]], float); ker /= ker.sum()
        Ut_pad = np.pad(Ut, ((1,1),(1,1)), mode='wrap')
        Ut = (ker[0,0]*Ut_pad[:-2,:-2] + ker[0,1]*Ut_pad[:-2,1:-1] + ker[0,2]*Ut_pad[:-2,2:] +
              ker[1,0]*Ut_pad[1:-1,:-2] + ker[1,1]*Ut_pad[1:-1,1:-1] + ker[1,2]*Ut_pad[1:-1,2:] +
              ker[2,0]*Ut_pad[2:,:-2] + ker[2,1]*Ut_pad[2:,1:-1] + ker[2,2]*Ut_pad[2:,2:])
    t_edges = np.linspace(0.0, np.pi, n_theta+1)
    p_edges = np.linspace(-np.pi, np.pi, n_phi+1)
    return Ut, t_edges, p_edges

# ---------- helpers ----------
def halfsizes_to_cells(half, units, Rb):
    ha, hb = half
    if units.lower() == "rb":
        return float(ha)*Rb, float(hb)*Rb
    return float(ha), float(hb)

def clamp_bounds(lo, hi, n):
    lo = max(0, int(np.floor(lo))); hi = min(n, int(np.ceil(hi)))
    if hi <= lo: hi = min(n, lo+1)
    return lo, hi

def upsample2d(A, k):
    if k <= 1: return A
    return np.repeat(np.repeat(A, k, axis=0), k, axis=1)

def panel_save(fig, outpath, dpi=FIG_DPI):
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def _colorbar(fig, im, ax):
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    label = "|u|" if ax.get_title().startswith("Flow speed") else ("|u_t|" if MAP_KIND=="ut" else MAP_KIND)
    cbar.set_label(label, fontsize=FONT_SIZES["cbar_label"])
    cbar.ax.tick_params(labelsize=FONT_SIZES["tick"])
    return cbar

# ---------- panel drawers ----------
def draw_slice_panels(U_last, center, Rb):
    spd = magnitude(U_last)
    nx, ny, nz = spd.shape
    ix, iy, iz = slice_indices(center, (nx, ny, nz))

    # XY
    S = spd[:, :, iz]
    hx, hy = halfsizes_to_cells(BOX_XY, CROP_UNITS, Rb)
    x0, x1 = clamp_bounds(center[0]-hx, center[0]+hx, nx)
    y0, y1 = clamp_bounds(center[1]-hy, center[1]+hy, ny)
    C = upsample2d(S[x0:x1, y0:y1], UPSAMPLE).T
    fig, ax = plt.subplots(figsize=FIG_SLICE_SIZE, constrained_layout=True)
    im = ax.imshow(C, origin="lower", aspect="equal",
                   extent=[x0, x1-1, y0, y1-1], interpolation=INTERP)
    ax.set_title(f"Flow speed |u| in XY slice (z = {iz:.0f})")
    ax.set_xlabel("x index"); ax.set_ylabel("y index")
    if ADD_GRID: ax.grid(True, alpha=0.20, linestyle='--', linewidth=0.5)
    if SHOW_BOX:
        ax.add_patch(plt.Rectangle((x0, y0), (x1-x0-1), (y1-y0-1),
                                   fill=False, edgecolor="red", linewidth=1.2))
    _colorbar(fig, im, ax)
    p_xy = os.path.join(IMG_DIR, "panel_A_XY_slice.png"); panel_save(fig, p_xy)

    # XZ
    S = spd[:, iy, :]
    hx, hz = halfsizes_to_cells(BOX_XZ, CROP_UNITS, Rb)
    x0, x1 = clamp_bounds(center[0]-hx, center[0]+hx, nx)
    z0, z1 = clamp_bounds(center[2]-hz, center[2]+hz, nz)
    C = upsample2d(S[x0:x1, z0:z1], UPSAMPLE).T
    fig, ax = plt.subplots(figsize=FIG_SLICE_SIZE, constrained_layout=True)
    im = ax.imshow(C, origin="lower", aspect="equal",
                   extent=[x0, x1-1, z0, z1-1], interpolation=INTERP)
    ax.set_title(f"Flow speed |u| in XZ slice (y = {iy:.0f})")
    ax.set_xlabel("x index"); ax.set_ylabel("z index")
    if ADD_GRID: ax.grid(True, alpha=0.20, linestyle='--', linewidth=0.5)
    if SHOW_BOX:
        ax.add_patch(plt.Rectangle((x0, z0), (x1-x0-1), (z1-z0-1),
                                   fill=False, edgecolor="red", linewidth=1.2))
    _colorbar(fig, im, ax)
    p_xz = os.path.join(IMG_DIR, "panel_B_XZ_slice.png"); panel_save(fig, p_xz)

    # YZ
    S = spd[ix, :, :]
    hy, hz = halfsizes_to_cells(BOX_YZ, CROP_UNITS, Rb)
    y0, y1 = clamp_bounds(center[1]-hy, center[1]+hy, ny)
    z0, z1 = clamp_bounds(center[2]-hz, center[2]+hz, nz)
    C = upsample2d(S[y0:y1, z0:z1], UPSAMPLE).T
    fig, ax = plt.subplots(figsize=FIG_SLICE_SIZE, constrained_layout=True)
    im = ax.imshow(C, origin="lower", aspect="equal",
                   extent=[y0, y1-1, z0, z1-1], interpolation=INTERP)
    ax.set_title(f"Flow speed |u| in YZ slice (x = {ix:.0f})")
    ax.set_xlabel("y index"); ax.set_ylabel("z index")
    if ADD_GRID: ax.grid(True, alpha=0.20, linestyle='--', linewidth=0.5)
    if SHOW_BOX:
        ax.add_patch(plt.Rectangle((y0, z0), (y1-y0-1), (z1-z0-1),
                                   fill=False, edgecolor="red", linewidth=1.2))
    _colorbar(fig, im, ax)
    p_yz = os.path.join(IMG_DIR, "panel_C_YZ_slice.png"); panel_save(fig, p_yz)

    return p_xy, p_xz, p_yz

def make_tangential_panel(U_last, center, Rb):
    tmap, t_edges, p_edges = tangential_map_trilinear(
        U_last, center, Rb, n_theta=N_THETA, n_phi=N_PHI, kind=MAP_KIND, smooth=True
    )
    fig, ax = plt.subplots(figsize=FIG_SURF_SIZE, constrained_layout=True)
    im = ax.imshow(tmap, origin="lower", aspect="auto",
                   extent=[p_edges[0], p_edges[-1], t_edges[0], t_edges[-1]],
                   interpolation=INTERP)
    title = {"ut":"Near-surface tangential speed |u_t| (r ≈ Rb)",
             "uphi":"Near-surface azimuthal component u_φ (r ≈ Rb)",
             "utheta":"Near-surface polar component u_θ (r ≈ Rb)"}[MAP_KIND]
    ax.set_title(title)
    ax.set_xlabel("azimuth φ (rad)")
    ax.set_ylabel("polar θ (rad)")
    _colorbar(fig, im, ax)
    p_tan = os.path.join(IMG_DIR, f"panel_D_surface_{MAP_KIND}.png")
    panel_save(fig, p_tan)
    return p_tan

# ---------- stitch with equal row heights ----------
def stitch_panels(panel_paths, labels=("A","B","C","D")):
    imgs = [Image.open(p).convert("RGB") for p in panel_paths]

    # Normalize width for all, then force equal heights within each row
    w_target = max(im.width for im in imgs)
    imgs = [im.resize((w_target, int(im.height * w_target / im.width))) for im in imgs]

    # Row heights (force both items in a row to the same height)
    row1_h = max(imgs[0].height, imgs[1].height)
    row2_h = max(imgs[2].height, imgs[3].height)

    imgs[0] = imgs[0].resize((w_target, row1_h))
    imgs[1] = imgs[1].resize((w_target, row1_h))
    imgs[2] = imgs[2].resize((w_target, row2_h))
    imgs[3] = imgs[3].resize((w_target, row2_h))

    W = w_target * 2
    H = row1_h + row2_h
    composite = Image.new("RGB", (W, H), (255, 255, 255))
    composite.paste(imgs[0], (0, 0))
    composite.paste(imgs[1], (w_target, 0))
    composite.paste(imgs[2], (0, row1_h))
    composite.paste(imgs[3], (w_target, row1_h))

    draw = ImageDraw.Draw(composite)
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", size=40)
    except Exception:
        font = ImageFont.load_default()
    positions = [(20, 15), (w_target+20, 15), (20, row1_h+15), (w_target+20, row1_h+15)]
    for label, pos in zip(labels, positions):
        draw.text(pos, label, fill=(0,0,0), font=font)

    out_png = os.path.join(IMG_DIR, "TIGER_flow_slices_tangential.png")
    composite.save(out_png, dpi=(FIG_DPI, FIG_DPI))
    out_tiff = os.path.splitext(out_png)[0] + ".tiff"
    composite.save(out_tiff, format="TIFF", compression="tiff_lzw", dpi=(FIG_DPI, FIG_DPI))
    return out_png, out_tiff

def main():
    U_all, center = load_field_and_center()
    Rb = load_Rb_from_config(default=10.0)

    if PHASE_LOCK:
        meta = load_rotation_meta()
        if meta is not None:
            omega_b, _, dt = meta
            U_last, used_idx = phase_locked_average(U_all, omega_b, dt, nbins=PHASE_BINS, target=PHASE_TARGET)
            print(f"Phase-locked average used {len(used_idx)} frames (bin of {PHASE_BINS}).")
        else:
            K = min(AVG_LAST_K, U_all.shape[0])
            print("No rotation_meta.npy found → averaging last", K, "frames.")
            U_last = U_all[-K:].mean(axis=0)
    else:
        K = min(AVG_LAST_K, U_all.shape[0]); U_last = U_all[-K:].mean(axis=0)

    p_xy, p_xz, p_yz = draw_slice_panels(U_last, center, Rb)
    p_tan = make_tangential_panel(U_last, center, Rb)
    out_png, out_tiff = stitch_panels([p_xy, p_xz, p_yz, p_tan])

    print("Wrote:")
    for p in [out_png, out_tiff, p_xy, p_xz, p_yz, p_tan]:
        print("  ", p)

if __name__ == "__main__":
    main()
