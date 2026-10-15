# E5 — Binding rule, tether force, and bond lifetime

These are the model choices to settle **before** the long tail-length sweeps.

1. **Binding rule.** In the MATLAB, a phage binds when rand > P, where P = n·(Rr/Rcone)². So it binds with
   probability 1 − P, which is almost certain whenever any receptor lies in the cone. The update slides bind with
   probability P. With Rr = 0.1 nm, the slides' rule gives about zero binding. The smoke test gave a bound fraction
   of 0.35 under the MATLAB rule and 0.00 under the slides' rule.
2. **Tether force on the surface walk.** It is off in the MATLAB (`0*F_proj`). This experiment compares the surface
   search time with the tether off and on (`--tether` sets the scale).
3. **Bond lifetime.** This uses the gp37–LPS force curve (Brzozowska 2019, `data/binding_force_data.mat`).
   - The sech-potential fit gives a ≈ 2665 kT and b ≈ 0.12 per nm.
   - The mean escape time over the truncated linear potential is computed as a function of rupture force. It ranges
     from microseconds at −50 pN to over an hour at −120 pN.

**Run.** `python run.py --tails 0.02 0.1 0.3 --reps 3`

**Outputs.** `rule_comparison.csv`, `tether_comparison.csv`, `escape_times.csv`, `rule_and_tether.png` and
`binding_energy.png`.
