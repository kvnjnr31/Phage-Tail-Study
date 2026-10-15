"""E6 - Movie of phage searching and binding on the rotating TIGER cell.

Runs the receptor model (optionally in the swimming-cell flow) and renders it
with tiger.viz.render_tiger_movie (folded in from tiger_encounter_coupled.py):
rotating capsule + receptors + TIGER.m flagellum; free phage in orange,
attached-and-searching in red, bound in green.

    python run.py --tail 0.1 --nsteps 1200 --every 2
    python run.py --tail 0.1 --flow stokes --speed 30

Output: results/tiger_tail<nm>nm.mp4 (GIF if ffmpeg is not installed)
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
import argparse

from tiger.flow import StokesSwimmerFlow, equal_volume_radius
from tiger.model import ReceptorSimParams, run_receptor_simulation
from tiger.viz import render_tiger_movie

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--tail", type=float, default=0.05, help="tail length (um)")
ap.add_argument("--nsteps", type=int, default=1200)
ap.add_argument("--every", type=int, default=2, help="physics steps per frame")
ap.add_argument("--P", type=int, default=60)
ap.add_argument("--LX", type=float, default=5.0, help="box (um); small box keeps phage in view")
ap.add_argument("--flow", choices=["none", "stokes"], default="none")
ap.add_argument("--speed", type=float, default=30.0)
ap.add_argument("--fps", type=int, default=30)
ap.add_argument("--seed", type=int, default=7)
ap.add_argument("--out", default=str(pathlib.Path(__file__).parent / "results"))
a = ap.parse_args()

p = ReceptorSimParams(nsteps=a.nsteps, P=a.P, LX=a.LX, Ltail=a.tail, mutate_tails=False, seed=a.seed)
flow = StokesSwimmerFlow((0, a.speed, 0), equal_volume_radius(p.Rb, p.bac_length)) if a.flow == "stokes" else None
res = run_receptor_simulation(p, flow=flow, record_every=a.every, warn_overlap=False, progress=True)
out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
path = render_tiger_movie(res, out / f"tiger_tail{a.tail * 1000:.0f}nm.mp4", fps=a.fps)
print(f"wrote {path}")
