# TIGER lattice Boltzmann swimmer

The **TIGER** model of a swimming *E. coli*: a D3Q19 BGK lattice Boltzmann fluid with Guo forcing, driven by a
torque ring on the cell body and a counter-rotating helical flagellar bundle (optional tail thrust and stresslet).
It produces the velocity fields that `tiger.flow.TigerLBMFlow` loads for the phage simulations (experiment E2).

All physical and forcing settings are in `simulation/tiger_simulation_config.yaml`: grid 96×64×64, τ = 0.9,
R_B = 10 cells = 1 µm (Δx = 0.1 µm), time scale set by a 100 Hz bundle rotation.

## Run and plot a flow slice

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

## Result with the default config (Oct 2026)

`figures/tiger_flow_slice_rotation_only.png`: 2000 steps, averaged over the last 500.

- The default config has **no thrust** (`use_tail: false`, `use_stresslet: false`). Only the body torque ring
  and the counter-rotating bundle act on the fluid.
- The cell does not translate (drift < 10⁻⁴ cells) and the mean axial flow is zero. The field is two
  counter-rotating swirls (body, and bundle ~1.5 µm behind it), decaying roughly as r⁻³ (median |u| falls ~9×
  from 1.5 to 3 µm), as expected for a torque dipole.
- Peak |u| ≈ 1×10⁻⁵ lattice units. The older saved field used for E2 (median ~470 µm/s after calibration) must
  therefore come from a different configuration, probably with thrust on.

**Next:** rerun with thrust (`use_tail: true`, or `use_stresslet: true` with `f_stress > 0`) to get the flow of a
swimming (pusher) cell, then recalibrate `tiger.flow.TigerLBMFlow`.

## Scripts

| script | what it does |
|---|---|
| `simulation/hybrid_simulation.py` | The solver: lattice, Guo forcing, TIGER body, force deposition. Run directly it saves every frame (very large). Reads the config from `TIGER_CONFIG` if set. |
| `simulation/run_tiger.py` | Low-memory runner for the same model (final + time-averaged field). |
| `simulation/plot_flow_slice.py` | 2D slices (x–y and y–z) of a saved run, with shell diagnostics. |
| `simulation/make_figure_tiger_slices_tangential.py` | Phase-locked slices and the near-surface tangential map, from the full frame stack. |
| `simulation/make_figure_decay_verification.py` | ⟨|u|⟩(r) decay on the y–z plane (`data/image/decay_profile_r_over_b_1_to_100.csv`). |
| `simulation/plot_velocity_frame.py`, `plot_vorticity.py`, `run_vorticity_diagnostic.py`, `generate_*` , `full_3D_orbit_plot.py`, `make_3D_movie.py` | Frame, vorticity and movie rendering from the full frame stack. |
| `simulation/TIGER_simulation.py`, `lattice_boltzmann_3D.py` | Early stand-alone TIGER kinematics and a minimal LBM scaffold. |
| `extract_back_slice.py`, `recovery_script.py` | Rear-view slice and Berg-style recovery of rotation rate and drag. |

## Author

Kevin L. Joiner
