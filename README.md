# Phage-Tail-Study — phage around a swimming E. coli: binding vs tail length

Code, data and manuscript for the TIGER project (Joiner & Luque). The 2025 violin-plot deliverables that were
previously at the top of this repo are in `archive/violin_study_2025/`; the figure is now produced by
`experiments/E4_literature_adsorption_vs_tail_length/violin_figure.py` from the cleaned data.

**Project goal:** simulate T4-like phage around the TIGER bacterium model and find how their encounter and
binding depend on **tail length**.

```
pip install -r requirements.txt
python -m pytest tests -q
python experiments/E1_binding_vs_tail_length_static/run.py
```

## Layout

| folder | what it is |
|---|---|
| `tiger/` | Python package, ported from the MATLAB code. `tiger.model` is the receptor model: a capsule-shaped E. coli with a helical receptor strip, bulk CTRW, surface walk, tail-cone binding and RATTLE tether. `tiger.flow` is the flow around the swimmer, either Stokes or the TIGER LBM field. `tiger.sweeps` runs tail-length sweeps. `tiger.viz` renders the TIGER movie. `tiger.encounter` has the encounter-rate simulators and theory. `tiger.analysis` has adsorption theory and data analysis. |
| `experiments/` | One folder per question. Each has a `run.py`, a `README.md`, and writes its output to `results/`. |
| `data/` | Cleaned input data: `adsorption_data.xlsx` (66 rate constants, with `codes`, `data_issues` and `changes` sheets), `adsorption_references.xlsx` (source URLs), `phage_characteristics.xlsx`, `bacteria_characteristics.xlsx` (incl. Maeda swim/tumble data), the gp37–LPS force curve and the Schwartz receptor data. See `data/README.md`. |
| `manuscript/` | LaTeX draft (`main.tex`, `supplement.tex`, `references.bib`, `Figures/`). Overleaf-ready: upload the folder as a project. |
| `models/tiger_lbm/` | The **TIGER** swimmer itself: a lattice-Boltzmann fluid driven by the rotating body and flagellar bundle (`hybrid_simulation.py`). It produces the velocity fields that `tiger.flow.TigerLBMFlow` loads; `simulation/run_tiger.py` + `plot_flow_slice.py` run it and plot 2D slices (see its README). The large `.npy` outputs stay in iCloud: `phage_dynamics/Projects/agent_based_modeling/code/bacteria_ABM/data/numerical/`. |
| `models/sphere_surface_walk/` | A separate Python model from 2025: surface walks on a spherical cell with absorbing receptor patches, Mg²⁺-modulated binding and a beads-on-a-stick tail. It already supports tail-length and receptor-density sweeps, and is useful as a cross-check of E1. |
| `tests/` | Sanity tests: theory limits, constraint checks, receptor placement. |
| `archive/` | `matlab_reference/` holds the 95 MATLAB originals of the ported code. `PORT_NOTES.md` maps each MATLAB file to its Python replacement and lists the fixes and the switches that were off in the MATLAB. `data_legacy/` holds the raw spreadsheets the cleaned data came from (incl. the 2026 adsorption file and URL list). `violin_study_2025/` holds the 2025 violin scripts, figures, Word report and templates. Also here: `python_legacy/` and `c_hpc/`. |

## Experiments

| # | question | main output |
|---|---|---|
| **E1** | How do the attached fraction, the bound fraction and the time to bind change with tail length on a **non-swimming** cell? | `binding_vs_tail.png`, survival curves |
| **E2** | The same question for phage carried by the flow around the **swimming TIGER** cell (Stokes swimmer, or the TIGER LBM field). | binding vs tail length for each flow, compared with E1 |
| **E3** | Encounter rate constant vs tail length, with capture radius Rb + Rp + L: diffusive and ballistic theory against the flow simulation. | `encounter_vs_tail.png` |
| **E4** | What do the literature adsorption rates say about tail length and phage family? Covers the bootstrap regression, the 2026 dataset and k_diff theory. | regression histograms, family tests, `violin_figure.py` → four-panel violin figure + Mann–Whitney tests |
| **E5** | Model choices that set what "binding" means: the binding rule (MATLAB vs slides), the tether force, and bond lifetime vs rupture force. | `rule_and_tether.png`, escape times |
| **E6** | Movie of phage searching and binding on the rotating TIGER cell. | `.mp4` |

Suggested order: **E5 → E1 → E2 → E3**, then compare against **E4**. E5 comes first because the two binding rules
give opposite answers, so the rule should be settled before running long sweeps.

## Open issues to settle

1. **Binding rule.** The MATLAB binds with probability 1 − P (P = fraction of the tail cone covered by receptors). The
   update slides bind with probability P. With the default receptor size (Rr = 0.1 nm), the slides' rule almost never
   binds. Run E5, then choose the rule and a physical value for Rr.
2. **The tether force is off in the MATLAB** (`0*F_proj`). Turn it on with `--tether`. Its scale is not calibrated yet;
   the gp37–LPS stiffness from E5 is a starting point.
3. **The TIGER LBM field is not usable at its raw scale.** With the saved run, the calibrated speeds come out at
   ~470 µm/s median, which is far above E. coli's ~30 µm/s. E2 therefore rescales the field to a target swim speed
   (`--lbm-target-speed`, default 30 µm/s). The LBM grid is also small (about ±3 µm around the cell), and beyond it
   the flow is zero. A fresh run of the default config
   (Oct 2026) has no thrust: the cell does not move. `models/tiger_lbm/simulation/tiger_capsule.py` now gives a
   swimming capsule (same shape as the phage model) with thrust; its field, in units of the swim speed, is the one
   to load into `TigerLBMFlow` next.
4. **Data value to check.** `data/phage_characteristics.xlsx` (family_typical sheet, from `phage_characteristics.m`) lists a
   980 nm Siphoviridae tail. The 2025 violin script used 98 nm for the same entry. One of them is a units slip; 980 nm
   pulls the Siphoviridae mean tail to 281 nm (median 165.5 nm). The Abedon 2023 T4 label was fixed (now Myoviridae).
5. **Unrestricted binding.** Binding is attempted for every unbound phage near a receptor, not only phage attached
   to the surface. This is kept from the MATLAB. Decide whether binding should require attachment first.
