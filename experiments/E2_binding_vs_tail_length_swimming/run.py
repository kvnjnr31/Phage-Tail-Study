"""E2 - Binding vs tail length for phage around a swimming TIGER cell.

Same receptor model as E1, but free phage are carried by the flow around the
swimming bacterium:

  --flow stokes : analytic Stokes flow past a sphere of the capsule's volume
                  swimming along its long axis at --speed um/s (fast, default).
  --flow lbm    : the TIGER lattice-Boltzmann velocity field (models/tiger_lbm),
                  time-averaged, interpolated onto the phage positions. Point
                  --field at velocity_field_small.npy (424 MB, kept in iCloud:
                  phage_dynamics/Projects/agent_based_modeling/code/bacteria_ABM/data/numerical).
  --flow none   : static cell (= E1), for a paired comparison.

    python run.py --flow stokes --speed 30
    python run.py --flow lbm --field "C:/Users/kvnjn/iCloudDrive/phage_dynamics/Projects/agent_based_modeling/code/bacteria_ABM/data/numerical/velocity_field_small.npy"

Outputs (results/<flow>/): summary.csv, survival.csv, binding_vs_tail.png
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tiger.flow import StokesSwimmerFlow, TigerLBMFlow, equal_volume_radius
from tiger.model import ReceptorSimParams
from tiger.sweeps import tail_length_sweep

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--flow", choices=["stokes", "lbm", "none"], default="stokes")
ap.add_argument("--speed", type=float, default=30.0, help="swim speed (um/s) for --flow stokes")
ap.add_argument("--field", default=None, help="velocity_field(_small).npy for --flow lbm")
ap.add_argument("--lbm-target-speed", type=float, default=30.0,
                help="rescale the LBM field so its median speed 1.5-3 um from the cell is this (um/s); "
                     "0 = use the raw calibrated field (hundreds of um/s, see tiger/flow.py)")
ap.add_argument("--tails", type=float, nargs="+", default=[0.01, 0.05, 0.1, 0.2, 0.3])
ap.add_argument("--reps", type=int, default=3)
ap.add_argument("--nsteps", type=int, default=1200)
ap.add_argument("--P", type=int, default=100)
ap.add_argument("--LX", type=float, default=15.0)
ap.add_argument("--rule", choices=["matlab", "slides"], default="matlab")
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--out", default=None)
a = ap.parse_args()

base = ReceptorSimParams(nsteps=a.nsteps, P=a.P, LX=a.LX, binding_rule=a.rule)
if a.flow == "stokes":
    R = equal_volume_radius(base.Rb, base.bac_length)
    flow = StokesSwimmerFlow(U=(0.0, a.speed, 0.0), R=R)   # swims toward +y (head end)
    factory = lambda L: flow
elif a.flow == "lbm":
    if not a.field:
        sys.exit("--flow lbm needs --field path/to/velocity_field_small.npy")
    flow = TigerLBMFlow(a.field, target_speed=a.lbm_target_speed or None)
    print(flow.describe())
    factory = lambda L: flow
else:
    factory = None

out = pathlib.Path(a.out or pathlib.Path(__file__).parent / "results" / a.flow); out.mkdir(parents=True, exist_ok=True)
summary, surv = tail_length_sweep(a.tails, a.reps, base, flow_factory=factory, seed=a.seed)
summary.to_csv(out / "summary.csv", index=False)
surv.to_csv(out / "survival.csv", index=False)

g = summary.groupby("tail")
m, s = g.mean(numeric_only=True), g.std(numeric_only=True)
x = m.index.to_numpy() * 1000
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for col, lab in (("frac_bound", "bound"), ("frac_attached", "attached")):
    axs[0].errorbar(x, m[col], yerr=s[col], marker="o", capsize=3, label=lab)
axs[0].set(xlabel="tail length (nm)", ylabel="fraction of phage"); axs[0].legend()
axs[1].errorbar(x, m["rate_bind"], yerr=s["rate_bind"], marker="o", capsize=3)
axs[1].set(xlabel="tail length (nm)", ylabel="binding rate (1/s per phage)")
fig.suptitle(f"E2: binding vs tail length, flow = {a.flow}")
fig.tight_layout(); fig.savefig(out / "binding_vs_tail.png", dpi=150)
print(m[["frac_attached", "frac_bound", "median_t_bind", "rate_bind"]].round(4).to_string())
print(f"results in {out}")
