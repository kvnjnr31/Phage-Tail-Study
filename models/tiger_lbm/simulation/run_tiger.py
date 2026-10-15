"""Run the TIGER hybrid lattice Boltzmann model with low memory use.

``hybrid_simulation.py`` keeps every velocity frame in memory and saves them all
(2000 frames x 96x64x64x3 doubles is about 19 GB). This runner uses the same
physics (it imports the solver, forcing and body from ``hybrid_simulation``) but
keeps only:

  u_final.npy     velocity field at the last step            (nx, ny, nz, 3) float32
  u_mean.npy      field averaged over the last ``--avg`` steps (nx, ny, nz, 3) float32
  path.npy        body centre at every step                    (steps, 3)
  tail_final.npy  TIGER tail points at the last step (body-relative units)
  meta.npy        [dx (m/cell), dt (s/step), dt (tu/step), Rb (cells), n averaged]

Run from models/tiger_lbm:

    python simulation/run_tiger.py                                  # 2000 steps, average last 500
    python simulation/run_tiger.py --steps 400 --avg 100            # quick check
    python simulation/run_tiger.py --config simulation/my_config.yaml --out data/numerical/thrust

Then plot with ``simulation/plot_flow_slice.py``. The full 2000-step run takes about
35 min on one CPU core.
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--config", default=os.path.join(HERE, "tiger_simulation_config.yaml"),
                help="YAML config (default: simulation/tiger_simulation_config.yaml)")
ap.add_argument("--steps", type=int, default=None, help="number of steps (default: from config)")
ap.add_argument("--avg", type=int, default=None, help="average the last N steps (default: steps // 4)")
ap.add_argument("--out", default=os.path.join(ROOT, "data", "numerical", "run"))
ap.add_argument("--log-every", type=int, default=100)
a = ap.parse_args()

# hybrid_simulation reads its config at import time; point it at the chosen file.
os.environ["TIGER_CONFIG"] = os.path.abspath(a.config)
sys.path.insert(0, HERE)
import hybrid_simulation as H  # noqa: E402

steps = a.steps or H.steps
avg = min(a.avg if a.avg is not None else steps // 4, steps)
H.print_calibration_summary()
print(f"config: {a.config}\nsteps: {steps}, averaging the last {avg}\n")

f, rho, u = H.initialize_lbm(H.nx, H.ny, H.nz, H.tau)
c = H.config["bacterium"]
body = H.TIGERBacterium(center=np.array(c["center"], float), Rb=float(c["Rb"]), kE=float(c["kE"]),
                        k=float(c["k"]), w=float(c["w"]), h=float(c["h"]))
mean = np.zeros((H.nx, H.ny, H.nz, 3))
n_avg, path = 0, []
t0 = time.time()
for step in range(steps):
    rho, u = H.macroscopic(f)
    t = step * H.dt
    body.update_geometry(t)
    F = H.build_forces(body, u.shape, t)
    f = H.collision_with_guo(f, rho, u, H.tau, F, H.dt)
    f = H.streaming(f)
    body.advect(u, H.dt)
    path.append(body.center.copy())
    if step >= steps - avg:
        mean += u
        n_avg += 1
    if step % a.log_every == 0:
        print(f"[{step:05d}] |u|max={np.linalg.norm(u, axis=-1).max():.3e}  centre={body.center}  "
              f"{time.time() - t0:.0f} s", flush=True)

os.makedirs(a.out, exist_ok=True)
np.save(os.path.join(a.out, "u_final.npy"), u.astype(np.float32))
np.save(os.path.join(a.out, "u_mean.npy"), (mean / max(n_avg, 1)).astype(np.float32))
np.save(os.path.join(a.out, "path.npy"), np.array(path))
np.save(os.path.join(a.out, "tail_final.npy"), np.vstack([body.tail_x, body.tail_y, body.tail_z]).T)
np.save(os.path.join(a.out, "meta.npy"), np.array([H.DX_PHYS_M, H.DT_PHYS_S, H.dt, H.RB_LB, n_avg]))
print(f"saved {a.out} ({time.time() - t0:.0f} s)")
