"""
This script reads a spreadsheet of lab measurements (adsorption_data.xlsx) and makes
four easy-to-compare violin plots. It groups the data in two ways:
by phage family (Podoviridae, Myoviridae, Siphoviridae) and by the bacteria they
infect (Enterobacteriaceae = "Entero.", Synechococcaceae = "Syn.", Vibrionaceae = "Vib.",
Pseudomonadaceae = "Psd."). It also shows simple summaries of phage tail and head sizes.

Because some numbers are very small, the charts use a log scale so everything fits
and can be compared fairly. The script produces one output figure with four
panels (A–D) as "figure_violin_panels.png".
"""

# Abbreviations:
#   Enterobacteriaceae - "Entero."
#   Synechococcaceae   - "Syn."
#   Vibrionaceae       - "Vib."
#   Pseudomonadaceae   - "Psd."

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from violinplot import plot_violins

# ------------------- load numeric matrix M -------------------
XLSX_FILE = "adsorption_data.xlsx" 

# Workbook encoding:
#   col 1: code k in {1..9} = families in triplets (1–3 sipho, 4–6 myo, 7–9 podo)
#   col 2: host code: 1=E. coli, 2=Synechococcus spp., 3=Ps. aeruginosa, 4=Vibrio spp.
#   col 3: adsorption rate (mL/h)
raw = pd.read_excel(XLSX_FILE, header=None)
M = raw.copy()
for j in M.columns:
    M[j] = pd.to_numeric(M[j], errors="coerce")
M = M.dropna(how="all", axis=1).dropna(how="all", axis=0).values
assert M.shape[1] >= 3, "Expected at least 3 numeric columns in adsorption_data.xlsx"

k_code    = M[:, 0].astype(int)
host_code = M[:, 1].astype(int)
rate      = M[:, 2].astype(float)

def select_rates(k_range, host_id):
    rows = np.isin(k_code, list(k_range)) & (host_code == host_id)
    vals = rate[rows]
    vals = vals[np.isfinite(vals) & (vals > 0)]
    return vals

# ------------------- (B) Phage Adsorption Rates by phage family -------------------
# sipho: k=1..3, myo: 4..6, podo: 7..9
sipho_ecoli = select_rates(range(1, 4), 1)
sipho_syn   = select_rates(range(1, 4), 2)
sipho_aeru  = select_rates(range(1, 4), 3)

myo_ecoli   = select_rates(range(4, 7), 1)
myo_syn     = select_rates(range(4, 7), 2)
myo_vibrio  = select_rates(range(4, 7), 4)

podo_ecoli  = select_rates(range(7, 10), 1)
podo_syn    = select_rates(range(7, 10), 2)
podo_aeru   = select_rates(range(7, 10), 3)
podo_vibrio = select_rates(range(7, 10), 4)

sipho_rates = np.concatenate([sipho_ecoli, sipho_syn, sipho_aeru])
myo_rates   = np.concatenate([myo_ecoli, myo_syn, myo_vibrio])
podo_rates  = np.concatenate([podo_ecoli, podo_syn, podo_aeru, podo_vibrio])

# ------------------- (A) Tail & Head arrays (nm) -------------------
sipho_tail = np.array([1.50e-07, 1.51e-07, 1.80e-07, 1.35e-07, 9.80e-08, 2.10e-07, 1.00e-07, 1.40e-07, 3.55e-07]) / 1e-9
myo_tail   = np.array([1.13e-07, 1.70e-07, 1.35e-07, 1.20e-07, 1.50e-07, 1.70e-07, 1.60e-07, 5.30e-08]) / 1e-9
podo_tail  = np.array([2.00e-08, 1.50e-08, 1.20e-08, 1.00e-08, 1.80e-08, 1.30e-08, 1.70e-08, 4.60e-08]) / 1e-9
sipho_head = 0.5 * np.array([60e-9, 60e-9, 80e-9, 60e-9, 55e-9, 53e-9, 60e-9, 81e-9]) / 1e-9
myo_head   = 0.5 * np.array([85e-9, 60e-9, 60e-9, 94e-9, 64e-9, 144.5e-9, 85.1e-9]) / 1e-9
podo_head  = 0.5 * np.array([60e-9, 70e-9, 63e-9, 70e-9, 61e-9, 60e-9, 67e-9, 60e-9]) / 1e-9

# ------------------- (D) Host–Phage aggregates -------------------
Ecoli         = np.concatenate([sipho_ecoli, myo_ecoli, podo_ecoli])
Synechococcus = np.concatenate([sipho_syn,   myo_syn,   podo_syn])
Vibrio        = np.concatenate([myo_vibrio,  podo_vibrio])
Aeruginosa    = np.concatenate([sipho_aeru,  podo_aeru])

# ------------------- build figure -------------------
plt.close("all")
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
(axA, axB), (axC, axD) = axes

# Common plotting kwargs for robust violins (log-scale + floor only)
KW = dict(y_floor=5e-9)

# ---- (A) Phage Tail & Head Summary ----
A_series = [
    podo_tail, podo_head, None,   # Podoviridae
    myo_tail,  myo_head,  None,   # Myoviridae
    sipho_tail, sipho_head        # Siphoviridae
]
A_labels = ["Tail", "Head", "", "Tail", "Head", "", "Tail", "Head"]
A_colors = ["#4F81BD", "#4F81BD", None, "#9BBB59", "#9BBB59", None, "#C0504D", "#C0504D"]
plot_violins(axA, A_series, A_labels, colors=A_colors, logy=True,
             title="(A) Phage Tail & Head Summary", xlabel="— Tail Length →", ylabel="Tail Length (nm)",
             **KW)
axA.set_ylim(5, 1e3)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(facecolor="#4F81BD", alpha=0.35, label="Podoviridae"),
                    Patch(facecolor="#9BBB59", alpha=0.35, label="Myoviridae"),
                    Patch(facecolor="#C0504D", alpha=0.35, label="Siphoviridae")],
           loc="upper left", frameon=False)

# ---- (B) Phage Adsorption Rates ----
B_series = [podo_rates, myo_rates, sipho_rates]
B_labels = ["Podoviridae", "Myoviridae", "Siphoviridae"]
B_colors = ["#4F81BD", "#9BBB59", "#C0504D"]
plot_violins(axB, B_series, B_labels, colors=B_colors, logy=True,
             title="(B) Phage Adsorption Rates", xlabel="— Tail Length →", ylabel="Adsorption Rate (mL/hr)",
             **KW)

# ---- (C) Phage–Host Adsorption by Bacterial Family (abbreviated) ----
C_series = [
    podo_ecoli, podo_aeru,  podo_vibrio, None,          # Podoviridae
    myo_ecoli,  myo_syn,    myo_vibrio,  None,          # Myoviridae
    sipho_ecoli, sipho_syn, sipho_aeru                   # Siphoviridae
]
C_labels = [
    "Entero.", "Psd.", "Vib.", "",                     # Podoviridae
    "Entero.", "Syn.", "Vib.", "",                     # Myoviridae
    "Entero.", "Syn.", "Psd."                          # Siphoviridae
]
C_colors = [
    "#4F81BD", "#4F81BD", "#4F81BD", None,
    "#9BBB59", "#9BBB59", "#9BBB59", None,
    "#C0504D", "#C0504D", "#C0504D"
]
plot_violins(
    axC, C_series, C_labels, colors=C_colors, logy=True,
    title="(C) Phage–Host Adsorption by Bacterial Family",
    xlabel="— Tail Length →", ylabel="Adsorption Rate (mL/hr)", **KW
)
lo, hi = axC.get_ylim()
axC.set_ylim(lo, hi * 1.5)
axC.margins(x=0.06)

# ---- (D) Host–Phage Adsorption by Bacterial Family ----
# NOTE: modeling features removed (no white median/IQR points, no red dashed line).
PODO = "#4F81BD"; MYO = "#9BBB59"; SIPHO = "#C0504D"

D_series = [Synechococcus, None, Ecoli, None, Vibrio, None, Aeruginosa]
D_labels = ["Syn.", "", "Entero.", "", "Vib.", "", "Psd."]
D_colors = ["#7F7F7F", None, "#7F7F7F", None, "#7F7F7F", None, "#7F7F7F"]

# Hide default black points here; we'll overplot colored ones only (raw data).
xmap = plot_violins(
    axD, D_series, D_labels, colors=D_colors, logy=True, show_points=False,
    title="(D) Host–Phage Adsorption by Bacterial Family",
    ylabel="Adsorption Rate (mL/hr)", **KW
)

lo, hi = axD.get_ylim()
axD.set_ylim(lo, hi * 1.5)
axD.margins(x=0.06)

# ---- Overplot colored raw points by phage family (kept) ----
rng = np.random.default_rng(4242)
jitter_w = 0.25
pt_kw = dict(s=22, edgecolors="k", linewidths=0.3, zorder=5)

def _plot_parts(x, parts):
    for vals, color in parts:
        v = np.asarray(vals, float)
        v = v[np.isfinite(v) & (v > 0)]
        if v.size == 0:
            continue
        jitter = (rng.random(v.size) - 0.5) * jitter_w
        axD.scatter(np.full(v.size, x) + jitter, v, c=color, **pt_kw)

# Syn. = Sipho + Myo (no Podo in your data)
_plot_parts(xmap["Syn."],   [(sipho_syn, SIPHO), (myo_syn, MYO)])
# Entero. = Sipho + Myo + Podo
_plot_parts(xmap["Entero."], [(sipho_ecoli, SIPHO), (myo_ecoli, MYO), (podo_ecoli, PODO)])
# Vib. = Myo + Podo
_plot_parts(xmap["Vib."],    [(myo_vibrio, MYO), (podo_vibrio, PODO)])
# Psd. = Sipho + Podo
_plot_parts(xmap["Psd."],    [(sipho_aeru, SIPHO), (podo_aeru, PODO)])

# ------------------- put 'Tail Length →' INSIDE the top of panels A, B, C -------------------
for ax in [axA, axB, axC]:
    ax.set_xlabel("")
    ax.text(0.5, 0.985, "— Tail Length →",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=10, fontstyle="italic", color="k")

# ------------------- cosmetics & save -------------------
for ax in axes.ravel():
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

fig.tight_layout()
fig.savefig("figure_violin_panels.png", dpi=300)
plt.show()

# ------------------- Statistical Verification Patch -------------------
from scipy.stats import mannwhitneyu

def describe_and_test(group_a, group_b, name_a, name_b):
    """Compute medians and Mann–Whitney U test between two groups."""
    a = np.asarray(group_a, float)
    b = np.asarray(group_b, float)
    a = a[np.isfinite(a) & (a > 0)]
    b = b[np.isfinite(b) & (b > 0)]
    if a.size == 0 or b.size == 0:
        return
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    print(f"\n{name_a} vs {name_b}")
    print(f"  n = {a.size}, {b.size}")
    print(f"  Median(a) = {np.median(a):.2e}, Median(b) = {np.median(b):.2e}")
    print(f"  Mann–Whitney U = {U:.2f}, p = {p:.3e}")

print("\n=== Statistical Verification: Adsorption Rates by Phage Family ===")
describe_and_test(sipho_rates, myo_rates, "Siphoviridae", "Myoviridae")
describe_and_test(sipho_rates, podo_rates, "Siphoviridae", "Podoviridae")
describe_and_test(myo_rates, podo_rates,  "Myoviridae",  "Podoviridae")

print("\n=== Statistical Verification: Adsorption Rates by Host ===")
describe_and_test(Ecoli, Synechococcus, "E. coli", "Synechococcus")
describe_and_test(Ecoli, Aeruginosa,    "E. coli", "P. aeruginosa")
describe_and_test(Ecoli, Vibrio,        "E. coli", "Vibrio")
describe_and_test(Synechococcus, Vibrio, "Synechococcus", "Vibrio")
describe_and_test(Synechococcus, Aeruginosa, "Synechococcus", "P. aeruginosa")
describe_and_test(Vibrio, Aeruginosa, "Vibrio", "P. aeruginosa")

print("\nVerification complete — results correspond to Mann–Whitney U-tests in the report.")
