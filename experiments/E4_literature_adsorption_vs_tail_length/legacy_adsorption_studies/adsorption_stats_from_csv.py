# generate_boxplots_adsorption_scaled.py
# Same plots as before, but with MUCH better scaling controls:
# - Optional LOG scale (base 10) for y to handle wide dynamic range
# - Robust y-limits from percentiles (e.g., 1st–99th) so extreme fliers don’t crush the boxes
# - Configurable whiskers by percentiles (e.g., 5th–95th) to match the robust limits
# - Light grid + nicer tick formatting

import os, io, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter, FuncFormatter

# ======================
# CONFIG
# ======================
from pathlib import Path
CSV_PATH          = str(Path(__file__).resolve().parents[3] / "archive" / "data_legacy" / "adsorption_data_original.csv")

TITLE_LEFT        = "A. Bacterial Host vs adsorption"
TITLE_RIGHT       = "B. Phage Family vs adsorption"
SUPTITLE          = "Adsorption box plots"

# ---- Scaling knobs ----
YSCALE_MODE       = "log"   # "log" or "linear"
LOG_BASE          = 10
PCT_YLIMS         = (1, 99) # set y-limits to these percentiles of the panel’s data
WHIS_PCTS         = (5, 95) # boxplot whiskers as percentiles (robust to outliers)
SHOW_FLIERS       = True    # True keeps fliers; limits still come from PCT_YLIMS
MIN_SAMPLES       = 2

# ---- Aesthetics ----
ROT               = 25
FIGSIZE           = (14, 6)
DPI               = 220
GRID_ALPHA        = 0.25
TIGHT_LAYOUT      = True

Y_LABEL_LINEAR    = r"Adsorption rate  $d\ (mL/h)$"
Y_LABEL_LOG       = r"Adsorption rate  $d\ (mL/h)$ (log$_{10}$)"

# ======================
# Robust CSV loader + mapping
# ======================
def read_csv_smart(path):
    raw = open(path, "rb").read()
    head = raw[:4096].decode("cp1252", errors="ignore")
    first = next((ln for ln in head.splitlines() if ln.strip()), "")
    delim = max([",", "\t", ";", "|"], key=lambda d: first.count(d)) if first else ","
    for enc in ("utf-8-sig", "utf-8", "cp1252", "latin1"):
        try:
            return pd.read_csv(io.StringIO(raw.decode(enc)), sep=delim, engine="python")
        except Exception:
            pass
    return pd.read_csv(io.StringIO(raw.decode("cp1252", errors="replace")), sep=delim, engine="python")

def _norm(s):
    s = str(s).replace("\ufeff", "").strip().lower()
    return re.sub(r"[^a-z0-9]+", "", s)

def map_columns(df):
    cols = list(df.columns)
    norm_map = {c: _norm(c) for c in cols}
    host_col = next((c for c,n in norm_map.items() if n in {"bacterialhost","host"}), None)
    fam_col  = next((c for c,n in norm_map.items() if n in {"phagefamily","family"}), None)
    y_col    = None
    for c in cols:
        if re.search(r"d\s*\(\s*m[lL]\s*/\s*h\s*\)", str(c)):
            y_col = c; break
    if y_col is None:
        y_col = next((c for c,n in norm_map.items()
                      if any(tok in n for tok in ("adsorptionrate","dmlh","rate","k"))), None)
    missing = []
    if host_col is None: missing.append("Bacterial Host")
    if fam_col  is None: missing.append("Phage Family")
    if y_col    is None: missing.append("d(mL/h)")
    if missing:
        raise KeyError(f"Could not resolve columns. Headers: {list(df.columns)}. Missing: {missing}")
    return df.rename(columns={host_col:"Bacterial Host", fam_col:"Phage Family", y_col:"AdsorptionRate"})

_num_pat = re.compile(r"[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?")
def to_numeric_force(x):
    if pd.isna(x): return np.nan
    if isinstance(x, (int,float,np.number)): return float(x)
    m = _num_pat.search(str(x).replace(",",""))
    return float(m.group(0)) if m else np.nan

# ======================
# Plot helpers
# ======================
def sorted_categories_by_median(df, cat_col, val_col, min_n=1):
    stats = (df.groupby(cat_col)[val_col].agg(['median','count']).reset_index())
    stats = stats[stats['count'] >= min_n].sort_values('median', kind='mergesort')
    return stats[cat_col].tolist()

def panel_data(df, group_col, y_col, min_samples=1):
    order = sorted_categories_by_median(df, group_col, y_col, min_samples)
    data  = [df.loc[df[group_col]==g, y_col].dropna().values for g in order]
    positions = list(range(1, len(order)+1))
    return order, data, positions

def apply_y_scale(ax, all_values):
    vals = np.array([v for arr in all_values for v in arr if np.isfinite(v)])
    if vals.size == 0:
        return
    lo = np.percentile(vals, PCT_YLIMS[0])
    hi = np.percentile(vals, PCT_YLIMS[1])
    # avoid zero/negative on log scale
    if YSCALE_MODE == "log":
        lo = max(lo, np.nextafter(0, 1))  # tiny positive
        ax.set_yscale('log', base=LOG_BASE)
        ax.set_ylim(lo, hi)
        # tidy log ticks
        ax.yaxis.set_minor_locator(LogLocator(base=LOG_BASE, subs=range(2, LOG_BASE)))
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel(Y_LABEL_LOG)
    else:
        ax.set_ylim(lo, hi)
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.set_ylabel(Y_LABEL_LINEAR)
    ax.grid(True, which='major', axis='y', alpha=GRID_ALPHA)

def boxplot_panel(ax, df, group_col, y_col, title):
    order, data, positions = panel_data(df, group_col, y_col, MIN_SAMPLES)
    if not order:
        ax.set_visible(False); return
    bp = ax.boxplot(
        data,
        positions=positions,
        whis=WHIS_PCTS,
        showfliers=SHOW_FLIERS,
        patch_artist=False
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(order, rotation=ROT, ha='right')
    ax.set_title(title)
    apply_y_scale(ax, data)

# ======================
# Main
# ======================
def main():
    df = read_csv_smart(CSV_PATH)
    df.columns = [str(c).replace("\ufeff","").strip() for c in df.columns]
    df = map_columns(df)
    df["AdsorptionRate"] = df["AdsorptionRate"].apply(to_numeric_force)
    core = df[["Bacterial Host","Phage Family","AdsorptionRate"]].dropna().copy()

    fig, (axL, axR) = plt.subplots(1, 2, figsize=FIGSIZE)
    boxplot_panel(axL, core, "Bacterial Host", "AdsorptionRate", TITLE_LEFT)
    boxplot_panel(axR, core, "Phage Family",   "AdsorptionRate", TITLE_RIGHT)
    if SUPTITLE: fig.suptitle(SUPTITLE)
    if TIGHT_LAYOUT: fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig("adsorption_boxplots_1x2_scaled.png", dpi=DPI)
    print("[OK] wrote adsorption_boxplots_1x2_scaled.png")

if __name__ == "__main__":
    main()
