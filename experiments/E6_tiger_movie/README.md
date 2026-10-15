# E6 — TIGER movie

This renders the receptor model on a rotating cell with a TIGER.m-style flagellum, using `tiger.viz.render_tiger_movie`.
That renderer was folded in from `04_Phage_Adsorption_Tail_Length/tiger_encounter_coupled.py`. The physics now comes
from `tiger.model`, so the movie shows the same model that the experiments run.

Phage are coloured by state:
- orange: free
- gold triangle: bound while in the bulk
- red: attached and searching
- green: bound

The body rotation is for display only. It is an exact symmetry of the receptor geometry and doesn't affect binding.
To couple real swimming hydrodynamics, use `--flow stokes`.

**Run.**
```
python run.py --tail 0.1 --nsteps 1200 --every 2
python run.py --tail 0.1 --flow stokes --speed 30
```

**Output.** `results/tiger_tail<nm>nm.mp4`. If ffmpeg is missing, it writes a GIF instead.
