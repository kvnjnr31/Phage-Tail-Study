"""TIGER: Python port of Kevin Joiner's MATLAB phage-bacteria modeling code.

Subpackages
-----------
tiger.model      Capsule-shaped E. coli, receptor layout, phage surface walks,
                 tail-receptor binding, tethering (RATTLE) and the receptor
                 simulation driver (from Phage_Bacteria_2/Simulation).
tiger.encounter  Phage-bacterium encounter-rate simulators (2D, 3D, ballistic,
                 Stokes-flow "hydro", mucus-stratified), CTRW noise, MSD tools,
                 capture/first-passage times and the swimming TIGER bacterium
                 (from Phage_Bacteria_1 and surface_walks / encounter_rates_hydro).
tiger.analysis   Adsorption-rate theory and data analysis (diffusive/ballistic
                 kernels, temperature scaling, tail-length bootstrap regression,
                 rupture-force escape times) (from Phage_Bacteria_2 top level + Data).

Units follow the MATLAB originals: simulations use micrometres and seconds;
analysis scripts use the units noted in each function (often cm and hours).
"""

__version__ = "0.1.0"
