"""E4 - four-panel violin figure of the literature adsorption data (manuscript Fig. 4).

Port of the 2025 ``data_violinplots.py`` (archive/violin_study_2025/) onto the cleaned
data in ../../data, so the figure and the statistics come from the same tables as the
rest of E4:

  (A) tail length and capsid size by phage family   <- phage_characteristics.xlsx
  (B) adsorption rate by phage family                <- adsorption_data.xlsx
  (C) adsorption rate by phage family x host family
  (D) adsorption rate by host family, points coloured by phage family

Families are ordered by tail length (Podo < Myo < Sipho). Mann-Whitney U tests between
families and between hosts are printed and written to results/violin_tests.csv.

    python violin_figure.py                   # capsid drawn as radius, as in the 2025 figure
    python violin_figure.py --capsid diameter
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu

from tiger.analysis import load_adsorption, load_phage_characteristics
from tiger.analysis.violin import FAMILY_COLORS, HOST_ABBREV, plot_violins

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
FAMILIES = ["Podoviridae", "Myoviridae", "Siphoviridae"]           # increasing tail length
HOSTS = ["Synechococcaceae", "Enterobacteriaceae", "Vibrionaceae", "Pseudomonadaceae"]

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--data", default=str(DATA))
ap.add_argument("--capsid", choices=["radius", "diameter"], default="radius")
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()
D, out = pathlib.Path(a.data), pathlib.Path(a.out)
out.mkdir(parents=True, exist_ok=True)

ads = load_adsorption(D / "adsorption_data.xlsx")
_, pools = load_phage_characteristics(D / "phage_characteristics.xlsx")
scale = 0.5 if a.capsid == "radius" else 1.0


def rates(family=None, host=None):
    m = np.ones(len(ads), bool)
    if family:
        m &= ads.family.eq(family).values
    if host:
        m &= ads.host.eq(host).values
    v = ads.rate.values[m]
    return v[np.isfinite(v) & (v > 0)]


KW = dict(y_floor=5e-9)
fig, axes = plt.subplots(2, 2, figsize=(10, 8))
(axA, axB), (axC, axD) = axes

# (A) tail and capsid sizes
sA, lA, cA = [], [], []
for i, f in enumerate(FAMILIES):
    sA += [pools[f]["tail"], np.asarray(pools[f]["capsid"]) * scale]
    lA += ["Tail", "Capsid"]
    cA += [FAMILY_COLORS[f]] * 2
    if i < len(FAMILIES) - 1:
        sA.append(None); lA.append(""); cA.append(None)
plot_violins(axA, sA, lA, colors=cA, title="(A) Phage Tail & Capsid Size",
             ylabel=f"Size (nm; capsid {a.capsid})", **KW)
axA.set_ylim(5, 2e3)
axA.legend(handles=[Patch(facecolor=FAMILY_COLORS[f], alpha=0.35, label=f) for f in FAMILIES],
           loc="upper left", frameon=False)

# (B) rate by family
plot_violins(axB, [rates(f) for f in FAMILIES], FAMILIES, colors=[FAMILY_COLORS[f] for f in FAMILIES],
             title="(B) Phage Adsorption Rates", ylabel="Adsorption Rate (mL/hr)", **KW)

# (C) family x host
sC, lC, cC = [], [], []
for i, f in enumerate(FAMILIES):
    for h in HOSTS:
        v = rates(f, h)
        if v.size:
            sC.append(v); lC.append(HOST_ABBREV[h]); cC.append(FAMILY_COLORS[f])
    if i < len(FAMILIES) - 1:
        sC.append(None); lC.append(""); cC.append(None)
xs = plot_violins(axC, sC, lC, colors=cC, title="(C) Phage–Host Adsorption by Bacterial Family",
                  ylabel="Adsorption Rate (mL/hr)", **KW)

# (D) host aggregates with points coloured by phage family
sD, lD = [], []
for i, h in enumerate(HOSTS):
    sD.append(rates(host=h)); lD.append(HOST_ABBREV[h])
    if i < len(HOSTS) - 1:
        sD.append(None); lD.append("")
xmap = plot_violins(axD, sD, lD, colors=["#7F7F7F" if s is not None else None for s in sD], show_points=False,
                    title="(D) Host–Phage Adsorption by Bacterial Family", ylabel="Adsorption Rate (mL/hr)", **KW)
rng = np.random.default_rng(4242)
for h in HOSTS:
    for f in FAMILIES:
        v = rates(f, h)
        if v.size:
            axD.scatter(xmap[HOST_ABBREV[h]] + (rng.random(v.size) - 0.5) * 0.25, v, c=FAMILY_COLORS[f],
                        s=22, edgecolors="k", linewidths=0.3, zorder=5)

for ax in (axC, axD):
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi * 1.5); ax.margins(x=0.06)
for ax in (axA, axB, axC):
    ax.text(0.5, 0.985, "— Tail Length →", transform=ax.transAxes, ha="center", va="top",
            fontsize=10, fontstyle="italic")
for ax in axes.ravel():
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.2)
fig.tight_layout()
fig.savefig(out / "figure_violin_panels.png", dpi=300)

# Mann-Whitney U tests
rows = []
for kind, groups in (("family", {f: rates(f) for f in FAMILIES}), ("host", {h: rates(host=h) for h in HOSTS})):
    for g1, g2 in combinations(groups, 2):
        x1, x2 = groups[g1], groups[g2]
        if x1.size and x2.size:
            U, p = mannwhitneyu(x1, x2, alternative="two-sided")
            rows.append(dict(kind=kind, a=g1, b=g2, n_a=x1.size, n_b=x2.size, median_a=np.median(x1),
                             median_b=np.median(x2), U=U, p=p))
tests = pd.DataFrame(rows)
tests.to_csv(out / "violin_tests.csv", index=False)
with pd.option_context("display.float_format", "{:.3g}".format, "display.width", 140):
    print(tests.to_string(index=False))
print(f"figure and tests in {out}")
