# E1 — Binding vs tail length, non-swimming cell

**Question.** On a static E. coli, how do attachment, binding probability and time to bind depend on phage tail length?

**Model.** `tiger.model.run_receptor_simulation` is the port of `encounter_rates_main.m`. The cell is a capsule
(Rb = 0.8 µm, half-length 1.9 µm) carrying 1000 receptors in a helical strip. Phage (Rp = 0.1 µm) make CTRW jumps
in the bulk, attach on contact, random-walk on a shell at Rb + (L + Rp)/2, and bind a receptor that lies inside
the tail cone (radius L + Rp). Tail mutation at the box walls is switched off, so each run keeps its own L.

**Run.**
```
python run.py --tails 0.01 0.05 0.1 0.2 0.3 --reps 5 --nsteps 2400
python run.py --rule slides          # binding rule from the update slides (see E5)
```

**Outputs (`results/`).**
- `summary.csv`: one row per tail length and replicate.
  - Fractions: attached, bound, bound given attached.
  - Times: median time to attach, median time to bind, median surface search time.
  - Event rates: attachment rate and binding rate per phage-second.
  - Mean encounter kernel.
- `survival.csv`: the fraction not yet attached or bound at each time.
- `binding_vs_tail.png` and `survival.png`.

**Smoke test (Oct 2026, 40 phage, 200 steps, 8 µm box).** Bound fraction rose from 0.63 to 0.76 and the binding
rate from 0.75 to 1.2 per second as L went from 20 to 300 nm. This used the MATLAB binding rule; use larger runs
before drawing conclusions.
