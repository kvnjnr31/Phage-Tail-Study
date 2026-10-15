"""E1 - Binding characteristics vs tail length, non-swimming cell.

Receptor model (tiger.model, port of encounter_rates_main.m) with the tail
length held fixed per run. For each tail length and replicate it records the
fraction of phage that attach to and bind on the cell, first-attachment /
first-binding times and their survival curves.

    python run.py                          # defaults below
    python run.py --tails 0.01 0.05 0.1 0.2 0.3 --reps 5 --nsteps 2400
    python run.py --rule slides            # binding rule from the update slides

Outputs (results/): summary.csv, survival.csv, binding_vs_tail.png, survival.png
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tiger.model import ReceptorSimParams
from tiger.sweeps import tail_length_sweep

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tails", type=float, nargs="+", default=[0.01, 0.05, 0.1, 0.2, 0.3], help="tail lengths (um)")
ap.add_argument("--reps", type=int, default=3)
ap.add_argument("--nsteps", type=int, default=1200)
ap.add_argument("--P", type=int, default=100)
ap.add_argument("--LX", type=float, default=15.0, help="box size (um); MATLAB default 15")
ap.add_argument("--rule", choices=["matlab", "slides"], default="matlab")
ap.add_argument("--tether", type=float, default=0.0)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()

out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
base = ReceptorSimParams(nsteps=a.nsteps, P=a.P, LX=a.LX, binding_rule=a.rule, tether_force_scale=a.tether)
summary, surv = tail_length_sweep(a.tails, a.reps, base, seed=a.seed)
summary.to_csv(out / "summary.csv", index=False)
surv.to_csv(out / "survival.csv", index=False)

g = summary.groupby("tail")
m, s = g.mean(numeric_only=True), g.std(numeric_only=True)
x = m.index.to_numpy() * 1000
fig, axs = plt.subplots(1, 3, figsize=(14, 4))
for ax, col, lab in ((axs[0], "frac_bound", "fraction bound"), (axs[0], "frac_attached", "fraction attached"),
                     (axs[1], "rate_bind", "binding rate (1/s per phage)"),
                     (axs[2], "median_t_bind", "median time to bind (s)")):
    ax.errorbar(x, m[col], yerr=s[col], marker="o", capsize=3, label=lab)
    ax.set_xlabel("tail length (nm)"); ax.legend()
fig.suptitle("E1: binding vs tail length (static cell)")
fig.tight_layout(); fig.savefig(out / "binding_vs_tail.png", dpi=150)

fig, ax = plt.subplots(figsize=(6, 4.5))
for L, d in surv.groupby("tail"):
    c = d.groupby("t").S_bind.mean()
    ax.plot(c.index, c.values, label=f"{L * 1000:.0f} nm")
ax.set(xlabel="time (s)", ylabel="fraction not yet bound", title="Binding survival")
ax.legend(title="tail"); fig.tight_layout(); fig.savefig(out / "survival.png", dpi=150)
print(m[["frac_attached", "frac_bound", "median_t_bind", "rate_bind"]].round(4).to_string())
print(f"results in {out}")
