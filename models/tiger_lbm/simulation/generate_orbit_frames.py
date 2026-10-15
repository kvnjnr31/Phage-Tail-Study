import numpy as np
from full_3D_orbit_plot import full_3D_orbit_plot

velocity_data = np.load("data/numerical/velocity_field.npy", allow_pickle=True)
bacterium_path = np.load("data/numerical/bacterium_path.npy", allow_pickle=True)
tail_coords = np.load("data/numerical/tail_coords.npy", allow_pickle=True)

for t in range(0, 180, 2):  # ~180° of orbit
    full_3D_orbit_plot(
        t=t,
        velocity_data=velocity_data,
        bacterium_path=bacterium_path,
        tail_coords=tail_coords,
        output_dir="data/image/tiger_orbit_movie",
        dpi=200
    )
