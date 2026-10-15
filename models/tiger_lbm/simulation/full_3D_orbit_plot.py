import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def full_3D_orbit_plot(t, velocity_data, bacterium_path, tail_coords,
                       output_dir="frames_orbit", dpi=300,
                       azim_base=-45, elev_base=20, azim_step=1.5):
    """
    Create a 3D plot of TIGER simulation at timestep `t`, rotating around the bacterium.
    """
    os.makedirs(output_dir, exist_ok=True)

    u = velocity_data[t]
    center = bacterium_path[t]
    tail_x = tail_coords[t][:, 0] + center[0]
    tail_y = tail_coords[t][:, 1] + center[1]
    tail_z = tail_coords[t][:, 2] + center[2]

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f"TIGER Simulation – Orbit Frame {t}")

    # Plot fluid vectors (sparse)
    step = 4
    X, Y, Z = np.meshgrid(np.arange(0, u.shape[0], step),
                          np.arange(0, u.shape[1], step),
                          np.arange(0, u.shape[2], step), indexing='ij')
    U = u[::step, ::step, ::step, 0]
    V = u[::step, ::step, ::step, 1]
    W = u[::step, ::step, ::step, 2]
    ax.quiver(X, Y, Z, U, V, W, length=1.5, normalize=True, color='purple', linewidth=0.5)

    # Plot head and tail
    ax.scatter(*center, c='blue', s=60, label='Head')
    ax.plot(tail_x, tail_y, tail_z, c='orange', label='Tail')
    ax.legend()

    # Axes limits
    ax.set_xlim(0, u.shape[0])
    ax.set_ylim(0, u.shape[1])
    ax.set_zlim(0, u.shape[2])

    # Set view angle for cinematic orbit
    azim = azim_base + azim_step * t
    ax.view_init(elev=elev_base, azim=azim)

    fig.tight_layout()
    save_path = os.path.join(output_dir, f"frame_{t:04d}.png")
    plt.savefig(save_path, dpi=dpi)
    plt.close()
