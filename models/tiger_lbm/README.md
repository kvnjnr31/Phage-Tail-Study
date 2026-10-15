# TIGER lattice Boltzmann swimmer

**Use `simulation/tiger_capsule.py` for flow fields.** It has a solid capsule body (the same capsule as the phage
model) and a bundle that produces thrust, so the cell actually swims. The original `hybrid_simulation.py` (below)
has no body and, with its default config, no thrust.

## Capsule swimmer with thrust (`tiger_capsule.py`)

```bash
python simulation/tiger_capsule.py                 # ~35 min → data/numerical/capsule_thrust/
python simulation/plot_flow_slice.py --run data/numerical/capsule_thrust --frame lab   # or --frame body
```

- **Body:** no-slip capsule, radius 0.8 µm, total length 5.4 µm (cylinder half-length 1.9 µm), matching
  `tiger.model`. It counter-rotates about its axis (moving-wall bounce-back) at the rate set by torque balance.
- **Bundle:** a helix 5 µm long behind the rear pole that pushes fluid backward (thrust) and spins it (torque).
- **Frame:** the cell is fixed in a periodic box. In steady state the wall drag balances the thrust, so the cell
  and bundle together are force-free and the far fluid streams past at −U, the swim speed.
- **Units:** Re = 0.08, so the flow is Stokes-like and linear. Fields are given as u/U; multiply by the swim speed
  you want (e.g. 30 µm/s). Grid 100×48×48 at 0.16 µm (16 × 7.7 × 7.7 µm), τ = 1.
- Settings: `simulation/tiger_capsule_config.yaml`. `run.init_U` starts the fluid near the swim speed; the mean
  flow otherwise takes ~4000 steps to build up.

**Result (Oct 2026, 6000 steps, averaged over the last 1000):** U = 1.28×10⁻³ lattice units (converged to ~1%).
Body surface speed is 2.3 U. In the lab frame (`figures/tiger_capsule_flow_slice_lab.png`) the cell is a pusher:
the bundle drives a backward jet up to ~4.7 U, fluid is drawn in along the sides, and the swirl around the body
reaches ~2.3 U at the surface. Median |u| falls from 1.8 U at 0.24 µm from the surface to 0.3 U at 1.9 µm.
`figures/tiger_capsule_flow_slice_body.png` shows the same flow as the cell sees it.

**Limits.** The box is only 7.7 µm across, so periodic images confine the flow; the closed recirculation cells in
the lab-frame plot come partly from them. The capsule surface is voxelised at 0.16 µm (5 cells per radius). The
bundle is a prescribed force helix, not a resolved filament; its torque sets the body rotation, which is not yet
calibrated to measured E. coli rates.

## Original model (`hybrid_simulation.py`)

The **TIGER** model of a swimming *E. coli*: a D3Q19 BGK lattice Boltzmann fluid with Guo forcing, driven by a
torque ring on the cell body and a counter-rotating helical flagellar bundle (optional tail thrust and stresslet).
It produces the velocity fields that `tiger.flow.TigerLBMFlow` loads for the phage simulations (experiment E2).

All physical and forcing settings are in `simulation/tiger_simulation_config.yaml`: grid 96×64×64, τ = 0.9,
R_B = 10 cells = 1 µm (Δx = 0.1 µm), time scale set by a 100 Hz bundle rotation.

### Run and plot a flow slice

From this folder (`models/tiger_lbm`):

```bash
pip install -r requirements.txt
python simulation/run_tiger.py                    # 2000 steps, ~35 min on one core → data/numerical/run/
python simulation/plot_flow_slice.py              # → figures/tiger_flow_slice_mean.png
```

`run_tiger.py` keeps only the final field and the average over the last steps (`--avg`, default a quarter of
the run), so a run needs ~10 MB instead of the ~19 GB that `hybrid_simulation.py` writes. Use
`--config other.yaml` to run a variant (e.g. with thrust on) and `--out` to keep runs apart.
`plot_flow_slice.py` draws (A) the x–y plane through the body centre (swimming axis horizontal) and (B) the y–z
cross-section through the body, and prints |u| on shells around the cell, the mean axial flow and the body drift.
`.npy` outputs are not committed (`data/numerical/` is git-ignored).

### Result with the default config (Oct 2026)

`figures/tiger_flow_slice_rotation_only.png`: 2000 steps, averaged over the last 500.

- The default config has **no thrust** (`use_tail: false`, `use_stresslet: false`). Only the body torque ring
  and the counter-rotating bundle act on the fluid.
- The cell does not translate (drift < 10⁻⁴ cells) and the mean axial flow is zero. The field is two
  counter-rotating swirls (body, and bundle ~1.5 µm behind it), decaying roughly as r⁻³ (median |u| falls ~9×
  from 1.5 to 3 µm), as expected for a torque dipole.
- Peak |u| ≈ 1×10⁻⁵ lattice units. The older saved field used for E2 (median ~470 µm/s after calibration) must
  therefore come from a different configuration, probably with thrust on.

It also has two solver issues: the velocity omits the half-force correction of Guo forcing, and the forcing
is multiplied by `dt = 0.02`, which weakens it 50-fold. `tiger_capsule.py` uses standard Guo forcing.

## Scripts

| script | what it does |
|---|---|
| `simulation/tiger_capsule.py`, `tiger_capsule_config.yaml` | Capsule body + thrusting bundle, force- and torque-free swimmer (see above). |
| `simulation/hybrid_simulation.py` | The solver: lattice, Guo forcing, TIGER body, force deposition. Run directly it saves every frame (very large). Reads the config from `TIGER_CONFIG` if set. |
| `simulation/run_tiger.py` | Low-memory runner for the same model (final + time-averaged field). |
| `simulation/plot_flow_slice.py` | 2D slices (x–y and y–z) of a saved run, with diagnostics. Capsule runs are drawn with the capsule outline, in units of U, in the body or lab frame. |
| `simulation/make_figure_tiger_slices_tangential.py` | Phase-locked slices and the near-surface tangential map, from the full frame stack. |
| `simulation/make_figure_decay_verification.py` | ⟨|u|⟩(r) decay on the y–z plane (`data/image/decay_profile_r_over_b_1_to_100.csv`). |
| `simulation/plot_velocity_frame.py`, `plot_vorticity.py`, `run_vorticity_diagnostic.py`, `generate_*` , `full_3D_orbit_plot.py`, `make_3D_movie.py` | Frame, vorticity and movie rendering from the full frame stack. |
| `simulation/TIGER_simulation.py`, `lattice_boltzmann_3D.py` | Early stand-alone TIGER kinematics and a minimal LBM scaffold. |
| `extract_back_slice.py`, `recovery_script.py` | Rear-view slice and Berg-style recovery of rotation rate and drag. |

## Author

Kevin L. Joiner
