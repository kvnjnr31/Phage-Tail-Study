# tiger — MATLAB → Python port notes (archived)

This package ports the phage–bacteria modeling MATLAB code in `archive/matlab_reference/` to Python
(numpy, scipy, pandas, matplotlib). The modules are reorganized by function rather
than by script, and the MATLAB `MatP`/`MatB` matrices are now named arrays. Each docstring
names the `.m` file it replaces.

```
pip install -r requirements.txt        # or: pip install -e .
python -m pytest tests -q              # 11 sanity tests (theory checks, constraints)
python experiments/E1_binding_vs_tail_length_static/run.py
python experiments/E3_encounter_rate_vs_tail_length/run.py
python experiments/E4_literature_adsorption_vs_tail_length/run.py
```

## Map: MATLAB → Python

| MATLAB (under `archive/matlab_reference/`) | Python |
|---|---|
| **Phage_Bacteria_2 / Simulation** | |
| `encounter_rates_main.m` | `tiger.model.simulation.run_receptor_simulation` (`ReceptorSimParams`) |
| `init_receptors.m`, `receptor_overlap.m` | `tiger.model.receptors` |
| `init_phage.m` | `tiger.model.phage.init_phage` (`PhageState` documents the MatP columns) |
| `Langevin_fun.m` (surface walk + tether) | `tiger.model.surface_walk.surface_step` |
| `noise_gen.m` (bulk CTRW + attachment) | `tiger.model.surface_walk.bulk_jumps` |
| `binding_fun.m` | `tiger.model.binding.binding_fun` |
| `boundarys_receptor.m`, `overlap.m` | `tiger.model.binding.boundarys_receptor`, `overlap` |
| `RATTLE.m` | `tiger.model.rattle.rattle_step` / `rattle` / `lj_force` |
| cylinder/sphere surface code (inline everywhere) | `tiger.model.geometry` (capsule level set, projection, mesh) |
| `simulation_plotter.m`, `sim_plotter.m` | `tiger.model.plotting` |
| **surface_walks/** and **encounter_rates_hydro/** (identical helpers) | |
| `encounter_rates_hydro.m`, `Langevin_fun.m`, `noise_gen.m`, `encounter_fun_hydro.m`, `boundarys_hydro.m` | `tiger.encounter.hydro.encounter_rates_hydro` (+ `hydro_advect`, `hydro_noise`, `hydro_encounter`, `hydro_boundaries`) |
| `plot_fun.m` (kernel theory) | `tiger.encounter.theory.diffusive_kernel`, `ballistic_kernel`, `plot_kernel` |
| **Phage_Bacteria_1** | |
| `Code for paper/encounter_rates_algorithm/*`, `Encounter Rates/encounter_rates_3D/*` | `tiger.encounter.simulators.encounter_rates_3d`, `boundarys_3d` |
| `Encounter Rates/encounter_rates_ballistic.m` | `encounter_rates_ballistic` |
| `Encounter Rates/encounter_rates_2D/*` | `encounter_rates_2d` |
| `Encounter Rates/encounter_rates_stratified/*` incl. `diffusion_vs_depth.m` | `encounter_rates_stratified`, `diffusion_vs_depth` |
| `Encounter Rates/stokes_vel.m` | `tiger.encounter.stokes` |
| `Encounter Rates/Diffusion Times/*` | `tiger.encounter.capture.diffusion_times`, `mean_time_to_capture` |
| `bacteria_phage_simulators/track_gen.m`, `phage_msd.m`, `phage_msd_ctrw.m`, `diffuser_msd.m` | `tiger.encounter.msd.ctrw_tracks`, `phage_msd_ctrw`, `phage_msd_sweep` |
| `bacteria_msd.m`, `propeller_msd.m` | `run_tumble_tracks`, `run_tumble_D` |
| `FIG2A_LEFT_PLOT.m`, `reflecting.m` | `fig2a_mucin_msd`, `reflect` |
| `bacteria_model/TIGER.m` (= `CODE/MATLAB/TIGER.m`) | `tiger.encounter.swimmer` (`animate_tiger` replaces movie2avi) |
| `adsorption_model.m`, `my_model.m`, `sherwood_numbers.m` | `tiger.encounter.theory.two_state_adsorption`, `sherwood`, `sherwood_kernel` |
| **Phage_Bacteria_2 top level + Data** | |
| `TheoreticalAdsorbtion.m` | `tiger.analysis.theoretical_adsorption` |
| `PhageTailExperiment.m` | `kdiff_size_map`, `kdiff_head_vs_tail`, `PHAGE_GEOMETRY_CM` |
| `DiffusionConstants.m` | `diffusion_constants` |
| `LinearModeling.m`, `LinearModelingMoldovan.m` | `linear_modeling_literature`, `linear_model_moldovan`, `t_even_rates` |
| `Data Modeling/phage_characteristics.m` | `family_diffusion` |
| `Data/first_passage.m`, `Data/data_analyzer.m` | `tiger.analysis.binding_energy` |
| `Data/data_converter.m` | `schwartz_saturation` |
| `Data/adsorption_model_bootstrap.m` | `bootstrap_tail_regression` (`fitlm` → `tiger.analysis.stats.ols`) |
| `Data/data_boxplots.m`, `AdsorptionRatesViolinPlots.m` | `family_ranksum`, `plot_family_boxplots`, `host_rates_log10` |
| `Data/model_fun.m`, `model_scratch.m` | `load_adsorption_data`, `kernel_predictors` |

## What changed relative to the MATLAB

**Programming bugs that were fixed.** These are unambiguous errors:
- `RATTLE.m`: both Newton loops were `while r_mag < eta` with `r_mag = 2` to start, so they
  never ran. There was also a distance typo (`xm(3)*MatB(:,3)`), the force had the wrong sign (+dU/dx),
  `Rb` was undefined, and the velocity stage used the old normal. The Python version is a working
  RATTLE, tested to keep the phage on the surface with tangent momentum. Note that `RATTLE.m` was never called by the driver.
- `Langevin_fun.m` (surface walk): the Heaviside `x/abs(x)` was a matrix division on vectors, and the
  rotation axis was the un-normalised gradient, which rescaled the kick length.
- `boundarys_receptor.m`: the code tested and changed row `kk` instead of `jj(kk)`, and "reactivate" set
  `MatP(kk,8) = 1`, which binds phage kk to receptor 1. Re-entering phage are now detached and unbound.
- `encounter_rates_hydro/noise_gen.m`: the code jumped rows `1..length(kk)` instead of `kk(p)`. This made no
  difference when every phage jumps every step (nu = 1).
- `adsorption_model_bootstrap.m`: heads were filtered with `tails > 0`.
- `model_scratch.m`: an infinite run time gave J_diff = Inf, because only NaN was caught.
- The stratified-mucus driver and the ballistic script called helpers with wrong argument lists and could not
  run. They have been reassembled into working functions.
- `fig2a_mucin_msd`: alpha > 1 makes Gamma(1-nu) negative, which gives complex times in MATLAB, so nu is capped at 0.999.

**Model switches that were off in the MATLAB.** These keep their MATLAB defaults and are now parameters:
- Tether force on the surface walk: `r1 = r0 + v + 0*F_proj` → `tether_force_scale=0.0`.
- Phage diffusion in the 3D "paper" encounter code: `0*sqrt(24*D0*dt)` → `diffusion_scale=0.0`.
- Bacterium tumbling in the hydro code: `0*Time >= run_time` → `tumble_angle_deg=None`.
- Hydro jump exclusion radius = Rb + Rp = Rcrit (`exclusion_radius=None`). This means diffusion can never reach
  the capture shell, so the kernel is about 0. Pass `exclusion_radius=Rb` to get contacts (about 27 µm³/s at v = 30).

**Questions about the model.** These kept the MATLAB behaviour:
- **Binding rule.** `binding_fun.m` binds when `rand > P`, where P = n·(Rr/Rcone)² is the fraction of the cone
  covered by receptors. So binding happens with probability 1 − P, which is almost certain once any receptor is in the cone.
  The update slides describe the opposite: bind with probability = covered fraction. `binding_rule="matlab"` is the
  default, and `"slides"` gives the slide version. With the default parameters, the slide rule gives essentially
  zero bindings (P ≈ 4e-7).
- **Which phage can bind.** Binding is attempted for every unbound phage, not only attached ones, so bound can exceed
  attached. Also, `Nact` in the kernel counts *attached* phage (`MatP(:,9) > 0`).
- **Kernel formula.** The `diffusive_kernel` theory term is kept as written (`Rcrit^2/sqrt(pi Dt t)`, `D = 2 D0 + ...`).
- **Spreadsheet columns.** The spreadsheets gained columns after the MATLAB was written, so columns are now read by header
  name. The xlsx column numbers in the MATLAB no longer line up.

## Not ported (judged outside TIGER, or plotting-only)
MAVERICK_GUI (55 kB GUI) and the population agent-based model (`Agent Based Model`, `Water Model_Old`),
`Alpha Omega/*`, mucin viscosity/diffusion scripts, `Encounter Rates/polymer_model`, `Galiean_Transform.m`,
`Miscellaneous/*` movie and plot scripts, `cluster_fun.m` (a percolation demo), Multigrid, third-party `grabit.m`,
`Violin.m`/`violinplot.m`/`boxplot2.m` (matplotlib covers these), and the `Documents\MATLAB` folder
(biofilm, fuel cells, perfect numbers, GNSS — no phage code). The sphere-based Python models already in
`../bacteriophage_ABM_surface_walks/` are separate and were left as they are.
