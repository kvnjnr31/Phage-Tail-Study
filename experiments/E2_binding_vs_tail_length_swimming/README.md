# E2 — Binding vs tail length around the swimming TIGER cell

**Question.** Does swimming change how binding depends on tail length? Phage are carried by the flow around the
cell instead of only diffusing.

**Model.** This is E1 plus advection of free phage by a `tiger.flow` field, expressed in the cell's frame.

| `--flow` | field |
|---|---|
| `stokes` (default) | Analytic Stokes flow past an equal-volume sphere (1.18 µm) swimming along the cell axis at `--speed` µm/s. |
| `lbm` | The TIGER lattice-Boltzmann field (`models/tiger_lbm`), time-averaged and interpolated. It is rescaled to `--lbm-target-speed` (default 30 µm/s) because the raw calibrated field is ~470 µm/s median (see the root README). |
| `none` | Static cell, the same as E1, for a paired comparison. |

**Run.**
```
python run.py --flow stokes --speed 30
python run.py --flow none
python run.py --flow lbm --field "C:/Users/kvnjn/iCloudDrive/phage_dynamics/Projects/agent_based_modeling/code/bacteria_ABM/data/numerical/velocity_field_small.npy"
```
The `velocity_field_small.npy` file is 424 MB and holds 45 frames on a 96×64×64 grid.

**Outputs.** `results/<flow>/summary.csv`, `survival.csv` and `binding_vs_tail.png`, with the same columns as E1.

**Caveats.** The flow moves only free phage; attached phage stay on the surface. The Stokes field is set to zero
between the capsule and the equivalent sphere. The LBM grid covers only about ±3 µm around the cell, and outside
it the flow is zero.
