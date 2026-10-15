# E3 — Encounter rate constant vs tail length

**Question.** How much does a longer tail raise the encounter rate, purely by enlarging the capture radius
Rc = Rb + Rp + L? This is the upper bound that receptor binding (E1/E2) can only reduce.

**Method.** The script compares three estimates.
- Diffusive limit: 4πD₀Rc.
- Ballistic limit: πRc²v.
- Simulation: `tiger.encounter.encounter_rates_hydro`, phage in the Stokes flow of a swimming sphere, with capture at Rc.

It sweeps tail length and swim speed.

**Run.** `python run.py --tails 0 0.05 0.1 0.2 0.3 0.5 --speeds 0 15 30`

**Outputs.** `results/encounter_vs_tail.csv` (µm³/s and mL/hr) and `encounter_vs_tail.png`.

**Background.** `background_phage_tail_first_passage/` holds the 2025 project summary for this question, which covers
the ballistic and diffusive limits and the bootstrap. It also holds Inkscape figures and the phage characteristics and
bootstrap data.

**Note.** Short runs underestimate the steady-state kernel, because the kernel is still relaxing from its initial
transient. Use `--nsteps` of 1000 or more for numbers you report.
