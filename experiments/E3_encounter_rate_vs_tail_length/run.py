"""E3 - Encounter (adsorption) rate constant vs tail length.

The tail extends the capture radius to Rc = Rb + Rp + L. This compares
  * diffusive limit   k_diff = 4 pi D0 Rc
  * ballistic limit   k_ball = pi Rc^2 v
  * simulation in the Stokes flow of a swimming cell (tiger.encounter.hydro)
for a range of tail lengths and swim speeds, in um^3/s and mL/hr.

    python run.py --tails 0 0.05 0.1 0.2 0.3 0.5 --speeds 0 15 30

Outputs (results/): encounter_vs_tail.csv, encounter_vs_tail.png
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from tiger.sweeps import encounter_vs_tail

UM3_S_TO_ML_HR = 1e-12 * 3600

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tails", type=float, nargs="+", default=[0.0, 0.05, 0.1, 0.2, 0.3, 0.5])
ap.add_argument("--speeds", type=float, nargs="+", default=[0.0, 15.0, 30.0])
ap.add_argument("--Rb", type=float, default=1.0)
ap.add_argument("--Rp", type=float, default=0.05)
ap.add_argument("--P", type=int, default=20000)
ap.add_argument("--nsteps", type=int, default=300)
ap.add_argument("--no-sim", action="store_true")
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)

frames = []
for v in a.speeds:
    d = encounter_vs_tail(a.tails, Rb=a.Rb, Rp=a.Rp, v=v, simulate=not a.no_sim, P=a.P, nsteps=a.nsteps)
    d["v"] = v
    frames.append(d)
df = pd.concat(frames, ignore_index=True)
for c in ("k_diff", "k_ball", "k_sim"):
    if c in df:
        df[c + "_mL_hr"] = df[c] * UM3_S_TO_ML_HR
df.to_csv(out / "encounter_vs_tail.csv", index=False)

fig, ax = plt.subplots(figsize=(7, 5))
d0 = df[df.v == df.v.iloc[0]]
ax.plot(d0["tail"] * 1000, d0.k_diff, "k--", label="diffusive 4 pi D Rc")
for v, d in df.groupby("v"):
    if v > 0:
        ax.plot(d["tail"] * 1000, d.k_ball, ":", label=f"ballistic, v = {v:g}")
    if "k_sim" in d:
        ax.plot(d["tail"] * 1000, d.k_sim, "o-", label=f"simulation, v = {v:g}")
ax.set(xlabel="tail length (nm)", ylabel="encounter kernel (um^3/s)", title="E3: encounter rate vs tail length")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out / "encounter_vs_tail.png", dpi=150)
print(df.round(3).to_string())
print(f"results in {out}")
