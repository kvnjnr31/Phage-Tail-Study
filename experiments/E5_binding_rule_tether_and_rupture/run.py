"""E5 - Model choices that decide what "binding" means, and how long a bond lasts.

1. Binding rule: MATLAB (bind if rand > P_covered, i.e. prob 1 - P) vs the update
   slides (bind with prob = covered fraction of the tail cone), across tail lengths.
   The two give opposite tail-length trends, so this must be settled first.
2. Tether force on the surface walk: off (MATLAB, 0*F_proj) vs on.
3. Bond lifetime: gp37-LPS force curve (Brzozowska 2019, data/binding_force_curve_gp37_LPS.csv),
   sech-potential fit, and mean escape time over the truncated potential vs
   rupture force (first_passage.m).

    python run.py --tails 0.02 0.1 0.3 --reps 2

Outputs (results/): rule_comparison.csv, tether_comparison.csv, escape_times.csv, *.png
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from tiger.analysis import escape_times, fit_sech_potential, load_force_curve, plot_binding_analysis
from tiger.model import ReceptorSimParams
from tiger.sweeps import tail_length_sweep

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tails", type=float, nargs="+", default=[0.02, 0.1, 0.3])
ap.add_argument("--reps", type=int, default=2)
ap.add_argument("--nsteps", type=int, default=600)
ap.add_argument("--P", type=int, default=100)
ap.add_argument("--tether", type=float, default=1e-3, help="tether force scale for the 'on' case")
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)

frames = []
for rule in ("matlab", "slides"):
    s, _ = tail_length_sweep(a.tails, a.reps, ReceptorSimParams(nsteps=a.nsteps, P=a.P, binding_rule=rule), progress=False)
    s["rule"] = rule; frames.append(s)
rules = pd.concat(frames); rules.to_csv(out / "rule_comparison.csv", index=False)
frames = []
for scale in (0.0, a.tether):
    s, _ = tail_length_sweep(a.tails, a.reps, ReceptorSimParams(nsteps=a.nsteps, P=a.P, tether_force_scale=scale),
                             progress=False)
    s["tether_force_scale"] = scale; frames.append(s)
teth = pd.concat(frames); teth.to_csv(out / "tether_comparison.csv", index=False)

fig, axs = plt.subplots(1, 2, figsize=(11, 4))
for k, d in rules.groupby("rule"):
    m = d.groupby("tail").frac_bound.mean(); axs[0].plot(m.index * 1000, m.values, "o-", label=k)
axs[0].set(xlabel="tail length (nm)", ylabel="fraction bound", title="binding rule"); axs[0].legend()
for k, d in teth.groupby("tether_force_scale"):
    m = d.groupby("tail").median_search_time.mean(); axs[1].plot(m.index * 1000, m.values, "o-", label=f"tether x{k:g}")
axs[1].set(xlabel="tail length (nm)", ylabel="median surface search time (s)", title="tether force"); axs[1].legend()
fig.tight_layout(); fig.savefig(out / "rule_and_tether.png", dpi=150)

fd = load_force_curve(DATA / "binding_force_curve_gp37_LPS.csv")
fit = fit_sech_potential(fd)
f, tau = escape_times()
pd.DataFrame(dict(rupture_force_pN=f, mean_escape_s=tau)).to_csv(out / "escape_times.csv", index=False)
plot_binding_analysis(fd).savefig(out / "binding_energy.png", dpi=150)
print(rules.groupby(["rule", "tail"]).frac_bound.mean().round(3).to_string())
print(f"sech fit a = {fit['a']:.0f} kT, b = {fit['b']:.3f} /nm")
print(f"results in {out}")
