
import numpy as np
from plot_vorticity import plot_vorticity_frame

u = np.load("data/velocity_fielddata/numerical/%s.npy")[164]
bacterium_center = np.load("data/bacterium_pathdata/numerical/%s.npy")[164]
plot_vorticity_frame(u, data/image/frame_idx=164, bacterium_center=bacterium_center, output_path="vorticity_data/image/frame_164data/image/%s.png")
