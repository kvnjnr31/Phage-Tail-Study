# generate_panel_frames.py

import numpy as np
from plot_velocity_frame import plot_velocity_frame

# === Load saved simulation data ===
velocity_data = np.load("data/numerical/velocity_field.npy", allow_pickle=True)
bacterium_path = np.load("data/numerical/bacterium_path.npy", allow_pickle=True)
tail_coords = np.load("data/numerical/tail_coords.npy", allow_pickle=True)

# === Loop through frames ===
for t in range(0, len(velocity_data), 2):
    plot_velocity_frame(
        u=velocity_data[t],
        bacterium=bacterium_path[t],
        t=t,
        output_dir="data/image/velocity_frames",
        tail_points=tail_coords[t],
        dpi=200
    )
